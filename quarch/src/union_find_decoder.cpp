#include "union_find_decoder.h"
#include "decoding_graph.h"
#include "logger.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <set>

namespace qrc {

    UFDecoder::UFDecoder(const stim::Circuit& circuit)
        : Decoder(circuit), cluster_id(0), current_instance(0) {
        
        n_detectors = circuit.count_detectors();
        
        log.debug("Building Union-Find decoder using DecodingGraph...");
        decoding_graph = to_decoding_graph(circuit);
        debug_print_raw_graph();
        
        debug_print_raw_edge_frames();
        
        log.debug("Decoding graph built with ", decoding_graph.count_detectors(), " detectors.");
        debug_print_raw_graph();
        log.debug("DECODER SETUP:");
        log.debug("  n_detectors = ", n_detectors);
        log.debug("  graph_vertices = ", decoding_graph.vertices().size());
        
        setup_boundaries_from_decoding_graph();
        build_edge_mapping();
        
        vertex_to_cluster.resize(decoding_graph.vertices().size(), nullptr);
    }

    void UFDecoder::debug_print_adjacency_matrix() {
        log.debug("=== Graph Adjacency Matrix ===");
        auto vertices = decoding_graph.vertices();
        
        // Print direct connections between all vertices
        for (uint v1 = 0; v1 < vertices.size(); v1++) {
            std::string connections = "Vertex " + std::to_string(v1) + " connects to: ";
            for (uint v2 = 0; v2 < vertices.size(); v2++) {
                // Check if edge exists
                bool connected = false;
                for (uint e = 0; e < edge_list.size(); e++) {
                    if ((edge_list[e].first == v1 && edge_list[e].second == v2) ||
                        (edge_list[e].first == v2 && edge_list[e].second == v1)) {
                        
                        // Add edge weight and temporal info
                        double weight = edge_weights[e];
                        bool is_temporal = (weight == 0.0);
                        
                        connections += std::to_string(v2) + "(" + 
                                    (is_temporal ? "T" : "S") + "), ";
                        connected = true;
                        break;
                    }
                }
            }
            log.debug(connections);
        }
    }


    void UFDecoder::setup_boundaries_from_decoding_graph() {
        log.debug("Setting up boundaries from DecodingGraph...");
        
        auto vertices = decoding_graph.vertices();
        vertex_boundaries.clear();
        
        // ONLY mark virtual boundary vertices
        for (uint i = 0; i < vertices.size(); i++) {
            const auto& vertex = vertices[i];
            
            if (vertex->detector == UINT_MAX) {
                vertex_boundaries[i] = 0; // Virtual boundary only
            }
        }
        
        log.debug("Boundary setup complete. Found ", vertex_boundaries.size(), " virtual boundary vertices");
    }

    void UFDecoder::build_edge_mapping() {
        log.debug("Building edge mapping...");
        auto vertices = decoding_graph.vertices();
        edge_list.clear();
        edge_observables_map.clear();
        edge_weights.clear();
        
        // Build detector-to-vertex mapping
        detector_to_vertex_id.clear();
        detector_to_vertex_id.resize(n_detectors, UINT_MAX);
        
        for (uint v = 0; v < vertices.size(); v++) {
            if (vertices[v]->detector < n_detectors) {
                detector_to_vertex_id[vertices[v]->detector] = v;
            }
        }
        
        // Process all edges
        for (uint v1 = 0; v1 < vertices.size(); v1++) {
            auto adj_list = decoding_graph.adjacency_list(vertices[v1]);
            for (auto v2_ptr : adj_list) {
                uint v2 = std::find(vertices.begin(), vertices.end(), v2_ptr) - vertices.begin();
                if (v1 < v2) {  // Each edge only once
                    edge_list.push_back({v1, v2});
                    
                    // Set edge weights: temporal = 0.0, spatial = 1.0
                    double edge_weight = 1.0;
                    auto edge = decoding_graph.get_edge(vertices[v1], vertices[v2]);
                    
                    if (edge != nullptr) {
                        // Check if this is a temporal edge
                        bool is_temporal = (v2_ptr == decoding_graph.get_next_round(vertices[v1]));
                        if (is_temporal) {
                            edge_weight = 0.0;
                            log.debug("TEMPORAL EDGE: ", v1, " -> ", v2,
                                    " between vertices ", vertices[v1]->detector, " and ", vertices[v2]->detector);
                        }
                        
                        // Store observables for correction application
                        // Always store with ordered vertex pairs
                        std::vector<uint> observables;
                        for (uint obs_id : edge->frames) {
                            observables.push_back(obs_id);
                        }
                        uint min_v = std::min(v1, v2);
                        uint max_v = std::max(v1, v2);
                        edge_observables_map[{min_v, max_v}] = observables;
                    } else {
                        uint min_v = std::min(v1, v2);
                        uint max_v = std::max(v1, v2);
                        edge_observables_map[{min_v, max_v}] = {};
                    }
                    
                    edge_weights.push_back(edge_weight);
                }
            }
        }
        
        edge_states.resize(edge_list.size(), EdgeSupport::NONE);
        
        // Sort edges by weight (temporal first, spatial second)
        sorted_edges.clear();
        for (uint e = 0; e < edge_list.size(); e++) {
            sorted_edges.push_back(e);
        }
        std::sort(sorted_edges.begin(), sorted_edges.end(),
                [this](uint a, uint b) { 
                    return edge_weights[a] < edge_weights[b]; 
                });
        
        log.debug("Built edge mapping: ", edge_list.size(), " edges");
        // Print adjacency matrix for debugging
        debug_print_adjacency_matrix();
        debug_print_edge_observables();
    }

    void UFDecoder::debug_print_edge_observables() {
        log.debug("=== EDGE OBSERVABLES MAP ===");
        log.debug("Total edges with observable effects: ", edge_observables_map.size());
        
        for (const auto& [edge_pair, obs_vec] : edge_observables_map) {
            uint v1 = edge_pair.first;
            uint v2 = edge_pair.second;
            
            if (!obs_vec.empty()) {
                std::string obs_str = "Observables: ";
                for (uint obs_id : obs_vec) {
                    obs_str += std::to_string(obs_id) + " ";
                }
                
                log.debug("Edge ", v1, " <-> ", v2, " affects ", obs_vec.size(), 
                        " observables: ", obs_str);
            }
        }
        log.debug("============================");
    }

    DecoderShotResult UFDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
        current_instance++;
        // Convert uint8_t vector to uint vector for logging
        std::vector<uint> syndrome_uint(syndrome.begin(), syndrome.end());
        std::string syndrome_str = vector_to_string(syndrome_uint);
        log.info("Input Syndrome: ", syndrome_str);
        // Clear previous state - delete ALL clusters from previous run
        for (auto cluster : clusters) {
            if (cluster != nullptr) {
                delete cluster;
            }
        }
        clusters.clear();
        cluster_id = 0;
        
        std::fill(vertex_to_cluster.begin(), vertex_to_cluster.end(), nullptr);
        vertex_has_syndrome.assign(decoding_graph.vertices().size(), false);
        edge_states.assign(edge_list.size(), EdgeSupport::NONE);
        
        extract_vertex_syndrome(syndrome);
        
        // PHASE TIMING START
        auto total_start = std::chrono::high_resolution_clock::now();
        
        // Phase 1: Initialize singleton clusters
        auto phase1_start = std::chrono::high_resolution_clock::now();
        initialize_clusters();
        auto phase1_end = std::chrono::high_resolution_clock::now();
        
        // Phase 2: Grow and merge clusters until frozen
        auto phase2_start = std::chrono::high_resolution_clock::now();
        grow_and_merge_clusters();
        auto phase2_end = std::chrono::high_resolution_clock::now();
        
        // Phase 3: Extract corrections from clusters
        auto phase3_start = std::chrono::high_resolution_clock::now();
        peel_clusters();
        auto phase3_end = std::chrono::high_resolution_clock::now();
        
        auto total_end = std::chrono::high_resolution_clock::now();
        
        // Calculate phase times in nanoseconds
        auto phase1_time = std::chrono::duration_cast<std::chrono::nanoseconds>(phase1_end - phase1_start).count();
        auto phase2_time = std::chrono::duration_cast<std::chrono::nanoseconds>(phase2_end - phase2_start).count();
        auto phase3_time = std::chrono::duration_cast<std::chrono::nanoseconds>(phase3_end - phase3_start).count();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(total_end - total_start).count();

        uint syndrome_weight = std::accumulate(syndrome.begin(), syndrome.end(), 0);
        
        if (total_time > 100000 || syndrome_weight > 10) {  // Log if >100μs or high weight
            log.debug("TIMING: HW=", syndrome_weight, 
                    " | Init=", phase1_time, "ns", 
                    " | Grow=", phase2_time, "ns", 
                    " | Peel=", phase3_time, "ns", 
                    " | Total=", total_time, "ns");
        }
        
        // Generate final correction
        DecoderShotResult result;
        result.correction.resize(circuit.count_observables(), 0);
        
        for (uint e = 0; e < edge_states.size(); e++) {
            if (edge_states[e] == EdgeSupport::MATCHING) {
                apply_edge_correction(e, result.correction);
            }
        }

        result.is_logical_error = is_logical_error(result.correction, syndrome, n_detectors, circuit.count_observables());
        log.info("Is logical error: ", result.is_logical_error? "TRUE" : "FALSE");
        return result;
    }

    void UFDecoder::extract_vertex_syndrome(const std::vector<uint8_t>& detector_syndrome) {
        for (uint det = 0; det < detector_syndrome.size(); det++) {
            if (detector_syndrome[det] == 1 && det < detector_to_vertex_id.size() && 
                detector_to_vertex_id[det] != UINT_MAX) {
                uint v = detector_to_vertex_id[det];
                vertex_has_syndrome[v] = true;
                log.debug("Detector ", det, " -> Vertex ", v, " syndrome=1");
            }
        }
    }

    void UFDecoder::initialize_clusters() {
        log.debug("Phase 1: Initializing singleton clusters");
        
        for (uint v = 0; v < vertex_has_syndrome.size(); v++) {
            if (vertex_has_syndrome[v]) {
                Cluster* cluster = new Cluster(cluster_id++, current_instance);
                cluster->parity = 0;
                cluster->is_frozen = false;
                cluster_add_vertex(cluster, v);  // This sets parity = 1
                clusters.push_back(cluster);
                vertex_to_cluster[v] = cluster;
                
                log.debug("Created cluster ", cluster->id, " for syndrome vertex ", v);
            }
        }
        
        log.debug("Initialized ", clusters.size(), " singleton clusters");
    }

    void UFDecoder::grow_and_merge_clusters() {
        log.debug("Phase 2: Growing and merging clusters");
        
        bool changes = true;
        int iteration = 0;
        
        // Track growth radius for all clusters
        int growth_radius = 1;
        
        // Track vertices by their distance from syndrome
        std::map<uint, int> vertex_distance;
        
        // Initialize: syndrome vertices at distance 0
        for (uint v = 0; v < vertex_has_syndrome.size(); v++) {
            if (vertex_has_syndrome[v]) {
                vertex_distance[v] = 0;
            }
        }
        
        while (changes) {
            changes = false;
            iteration++;
            log.debug("=== Growth iteration ", iteration, " (radius=", growth_radius, ") ===");
            
            // Check and freeze completed clusters
            freeze_completed_clusters();
            
            // Check if all clusters are frozen
            bool all_frozen = true;
            for (auto& cluster : clusters) {
                if (cluster && cluster->instance == current_instance) {
                    Cluster* root = cluster->find();
                    if (!root->is_frozen) {
                        all_frozen = false;
                        break;
                    }
                }
            }
            
            if (all_frozen) {
                log.debug("All clusters frozen - growth complete");
                break;
            }
            
            // STEP 1: Find vertices at current radius to grow from
            std::set<uint> vertices_at_current_radius;
            for (const auto& [vertex, distance] : vertex_distance) {
                if (distance == growth_radius - 1) {
                    vertices_at_current_radius.insert(vertex);
                }
            }
            
            log.debug("Found ", vertices_at_current_radius.size(), " vertices at radius ", (growth_radius - 1));
            
            // STEP 2: Process growth edges
            std::vector<std::pair<uint, uint>> merge_edges;
            
            for (uint sorted_idx = 0; sorted_idx < sorted_edges.size(); sorted_idx++) {
                uint e = sorted_edges[sorted_idx];
                
                if (edge_states[e] != EdgeSupport::NONE) continue;
                
                uint v1 = edge_list[e].first;
                uint v2 = edge_list[e].second;
                Cluster* c1 = vertex_to_cluster[v1];
                Cluster* c2 = vertex_to_cluster[v2];
                
                // Skip if either cluster is frozen
                if ((c1 && is_cluster_frozen(c1)) || (c2 && is_cluster_frozen(c2))) {
                    continue;
                }
                
                // Case 1: Growth - one vertex has cluster, other is free
                if ((c1 && !c2) || (!c1 && c2)) {
                    Cluster* active_cluster = c1 ? c1 : c2;
                    uint cluster_vertex = c1 ? v1 : v2;  // Vertex already in cluster
                    uint free_vertex = c1 ? v2 : v1;     // Vertex to add
                    
                    // STRICT RADIUS: Only grow from vertices at exactly (radius-1)
                    bool is_radius_vertex = vertices_at_current_radius.count(cluster_vertex) > 0;
                    
                    if (is_radius_vertex && !is_cluster_frozen(active_cluster)) {
                        cluster_add_vertex(active_cluster, free_vertex);
                        active_cluster->add_tree_edge(e);
                        edge_states[e] = EdgeSupport::GROWN;
                        vertex_to_cluster[free_vertex] = active_cluster;
                        vertex_distance[free_vertex] = growth_radius;  // Set distance of new vertex
                        changes = true;
                        
                        log.debug("Cluster ", active_cluster->id, " grew from vertex ", 
                                 cluster_vertex, " to vertex ", free_vertex, " (radius=", growth_radius, ")");
                    }
                }
                // Case 2: Potential merge - both vertices have different clusters
                else if (c1 && c2 && c1->find() != c2->find()) {
                    bool v1_at_radius = vertices_at_current_radius.count(v1) > 0;
                    bool v2_at_radius = vertices_at_current_radius.count(v2) > 0;
                    
                    if ((v1_at_radius || v2_at_radius) && 
                        !is_cluster_frozen(c1) && !is_cluster_frozen(c2)) {
                        // Store edge for merging AFTER growth
                        merge_edges.push_back({v1, v2});
                    }
                }
            }
            
            // STEP 3: Process merges AFTER growth
            for (auto& merge_pair : merge_edges) {
                uint v1 = merge_pair.first;
                uint v2 = merge_pair.second;
                
                Cluster* c1 = vertex_to_cluster[v1];
                Cluster* c2 = vertex_to_cluster[v2];
                
                if (c1 && c2 && c1->find() != c2->find() && 
                    !is_cluster_frozen(c1) && !is_cluster_frozen(c2)) {
                    
                    Cluster* root1 = c1->find();
                    Cluster* root2 = c2->find();
                    
                    log.debug("Merging clusters ", root1->id, " (parity=", root1->parity, 
                            ") and ", root2->id, " (parity=", root2->parity, ")");
                    
                    // Find merging edge
                    uint merge_edge = UINT_MAX;
                    for (uint e = 0; e < edge_list.size(); e++) {
                        if (edge_states[e] == EdgeSupport::NONE) {
                            uint ev1 = edge_list[e].first;
                            uint ev2 = edge_list[e].second;
                            if ((ev1 == v1 && ev2 == v2) || (ev1 == v2 && ev2 == v1)) {
                                merge_edge = e;
                                break;
                            }
                        }
                    }
                    
                    if (merge_edge != UINT_MAX) {
                        log.debug("MERGING via edge ", v1, " -> ", v2, 
                                " (temporal: ", (edge_weights[merge_edge] == 0.0 ? "YES" : "NO"), ")");
                        // Save all vertices from both clusters
                        std::vector<uint> all_vertices = root1->vertices;
                        all_vertices.insert(all_vertices.end(), root2->vertices.begin(), root2->vertices.end());
                        
                        // Save both clusters' tree edges first
                        std::vector<uint> combined_edges = root1->tree_edges;
                        combined_edges.insert(combined_edges.end(), root2->tree_edges.begin(), root2->tree_edges.end());
                        
                        // Do union
                        root1->union_with(root2);
                        Cluster* final_root = root1->find();
                        
                        final_root->vertices = all_vertices;
                        // Copy all edges to final root
                        final_root->tree_edges = combined_edges;
                        final_root->add_tree_edge(merge_edge);

                        // THEN update parity on the final root
                        final_root->parity = root1->parity ^ root2->parity;
                        
                        // Update all vertex mappings
                        for (uint v : all_vertices) {
                            vertex_to_cluster[v] = final_root;
                        }
                        
                        // Clean up old clusters
                        auto it1 = std::find(clusters.begin(), clusters.end(), root1);
                        if (it1 != clusters.end() && root1 != final_root) {
                            *it1 = nullptr;
                        }
                        auto it2 = std::find(clusters.begin(), clusters.end(), root2);
                        if (it2 != clusters.end() && root2 != final_root) {
                            *it2 = nullptr;
                        }
                        
                        edge_states[merge_edge] = EdgeSupport::GROWN;
                        changes = true;
                        
                        log.debug("Merged result: cluster ", final_root->id, " parity=", final_root->parity);
                    }
                }
            }
            
            // Increment radius for next iteration
            growth_radius++;
            
            if (iteration > 100) {
                log.error("Growth failed to converge - force freezing remaining clusters");
                for (auto& cluster : clusters) {
                    if (cluster && cluster->instance == current_instance) {
                        Cluster* root = cluster->find();
                        if (!root->is_frozen) {
                            root->is_frozen = true;
                        }
                    }
                }
                break;
            }
        }
        
        log.debug("Growth complete after ", iteration, " iterations");
    }

        uint UFDecoder::find_best_match_for_leaf(uint leaf, const std::map<uint, std::set<uint>>& adjacency,
            const std::set<uint>& remaining_vertices) {
        uint best_match = -1;
        
        // Priority 1: Regular neighbor with syndrome (non-boundary)
        for (uint neighbor : adjacency.at(leaf)) {
            if (remaining_vertices.count(neighbor) > 0 && 
                vertex_boundaries.count(neighbor) == 0 && 
                neighbor < vertex_has_syndrome.size() && 
                vertex_has_syndrome[neighbor]) {
                log.debug("  Found best match for leaf ", leaf, ": neighbor ", neighbor, " with syndrome");
                return neighbor;  // Return immediately if highest priority match found
            }
        }
        
        // Priority 2: Any boundary neighbor
        for (uint neighbor : adjacency.at(leaf)) {
                if (remaining_vertices.count(neighbor) > 0 && 
                    vertex_boundaries.count(neighbor) > 0) {
                    log.debug("  Found boundary match for leaf ", leaf, ": neighbor ", neighbor);
                    best_match = neighbor;
                    break;
                }
        }
        
        // Priority 3: Any regular neighbor (non-boundary) (only if no regular match found)
        if (best_match == -1) {
            for (uint neighbor : adjacency.at(leaf)) {
                if (remaining_vertices.count(neighbor) > 0 && 
                    vertex_boundaries.count(neighbor) == 0) {
                    log.debug("  Found regular match for leaf ", leaf, ": neighbor ", neighbor);
                    best_match = neighbor;
                    break;
                }
            }
        }
        
        return best_match;
    }

    void UFDecoder::peel_clusters() {
        log.debug("Phase 3: Peeling clusters to recover the error pattern");
        uint total_peel_iterations = 0;
        uint clusters_peeled = 0;
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;  // ← Skip nullified clusters
            Cluster* root = cluster->find();
            if (root->instance != current_instance) continue;
            
            log.debug("Peeling cluster ", root->id);
            clusters_peeled++;
            
            // Build edge lookup table from stored tree edges
            std::map<std::pair<uint, uint>, uint> vertex_pair_to_edge;
            std::map<uint, std::set<uint>> adjacency;
            
            // Initialize vertices
            for (uint vertex_id : root->vertices) {
                adjacency[vertex_id] = {};
            }
            
            // Build adjacency and edge lookup from stored edges (one pass)
            for (uint edge_id : root->tree_edges) {
                uint v1 = edge_list[edge_id].first;
                uint v2 = edge_list[edge_id].second;
                
                adjacency[v1].insert(v2);
                adjacency[v2].insert(v1);
                vertex_pair_to_edge[{std::min(v1,v2), std::max(v1,v2)}] = edge_id;
            }
            
            // Simple leaf peeling with direct edge access
            std::set<uint> remaining_vertices(root->vertices.begin(), root->vertices.end());
            uint cluster_peel_iterations = 0;
            log.debug("All remaining vertices: ");
            std::string vertices_str = "";
            for (uint v : remaining_vertices) {
                vertices_str += std::to_string(v) + " ";
            }
            log.debug(vertices_str);
            
            while (remaining_vertices.size() > 1) {
                cluster_peel_iterations++;
                // Find leaves - BUT EXCLUDE VIRTUAL BOUNDARIES
                print_spanning_tree(adjacency, remaining_vertices);
                std::vector<uint> leaves;
                for (uint vertex_id : remaining_vertices) {
                    uint degree = 0;
                    for (uint neighbor : adjacency[vertex_id]) {
                        if (remaining_vertices.count(neighbor) > 0) {
                            degree++;
                        }
                    }
                    // Only consider non-boundary vertices as leaves
                    if (degree == 1 && vertex_boundaries.count(vertex_id) == 0) {
                        leaves.push_back(vertex_id);
                    }
                }
                
                if (leaves.empty()) break;

                log.debug("  Peel iteration ", cluster_peel_iterations, ": removing ", leaves.size(), " leaves");
                // Process leaves
                for (uint leaf : leaves) {
                    // Check syndrome and mark corrections
                    bool leaf_has_syndrome = (leaf < vertex_has_syndrome.size() && 
                                            vertex_has_syndrome[leaf]);

                    if (leaf_has_syndrome) {
                        uint best_neighbor = find_best_match_for_leaf(leaf, adjacency, remaining_vertices);
                        
                        if (best_neighbor != -1) {
                            auto edge_it = vertex_pair_to_edge.find({std::min(leaf, best_neighbor), std::max(leaf, best_neighbor)});
                            if (edge_it != vertex_pair_to_edge.end()) {
                                edge_states[edge_it->second] = EdgeSupport::MATCHING;

                                if (vertex_boundaries.count(best_neighbor) > 0) {
                                    log.debug("  Matching leaf ", leaf, " with boundary ", best_neighbor);
                                } else {
                                    log.debug("  Matching leaf ", leaf, " with neighbor ", best_neighbor);
                                }
                            }
                            if (vertex_boundaries.count(best_neighbor) == 0) {
                                vertex_has_syndrome[best_neighbor] ^= 1;   // flip parent
                            }
                        }
                    }

                    // Remove leaf REGARDLESS of syndrome status
                    remaining_vertices.erase(leaf);  // ← Move this outside the if statement
                }
            }

            total_peel_iterations += cluster_peel_iterations;
            log.debug("Cluster ", root->id, " peeled in ", cluster_peel_iterations, " iterations");
        }
        log.debug("Peeling complete: ", clusters_peeled, " clusters, ", total_peel_iterations, " total iterations");
    }
    
    bool UFDecoder::should_freeze_cluster(Cluster* cluster) {
        Cluster* root = cluster->find();
        
        // Freeze condition 1: Even parity (all defects paired)
        if (root->parity % 2 == 0) {
            return true;
        }
        
        // Freeze condition 2: Odd parity + CONTAINS virtual boundary vertex
        if (root->parity % 2 == 1 && cluster_contains_virtual_boundary(root)) {
            return true;
        }
        
        // Allow odd parity clusters to keep growing if they don't contain virtual boundary
        return false;
    }

    bool UFDecoder::cluster_contains_virtual_boundary(Cluster* cluster) {
        // ONLY return true if cluster CONTAINS virtual boundary vertices
        for (uint vertex_id : cluster->vertices) {
            auto it = vertex_boundaries.find(vertex_id);
            if (it != vertex_boundaries.end() && it->second == 0) {
                return true; // Virtual boundary vertex is IN the cluster
            }
        }
        return false;
    }

    void UFDecoder::freeze_completed_clusters() {
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;  // ← Skip nullified clusters
            
            Cluster* root = cluster->find();
            if (root->instance != current_instance || root->is_frozen) continue;
            
            if (should_freeze_cluster(root)) {
                root->is_frozen = true;
                log.debug("Froze cluster ", root->id, " (parity=", root->parity, ")");
            }
        }
    }

    bool UFDecoder::is_cluster_frozen(Cluster* cluster) {
        if (cluster == nullptr) return false;
        return cluster->find()->is_frozen;
    }

    void UFDecoder::cluster_add_vertex(Cluster* cluster, uint vertex_id) {
        cluster->add_vertex(vertex_id);
        // Only increment parity if vertex has syndrome
        if (vertex_id < vertex_has_syndrome.size() && vertex_has_syndrome[vertex_id]) {
            cluster->parity ^= 1;
        }
        vertex_to_cluster[vertex_id] = cluster;
    }

    void UFDecoder::debug_observable_mapping() {
        log.debug("=== OBSERVABLE 0 MAPPING DEBUG ===");
        
        int edges_affecting_obs0 = 0;
        
        for (auto& [edge_pair, obs_vector] : edge_observables_map) {
            uint v1 = edge_pair.first, v2 = edge_pair.second;
            
            // Check if this edge affects Observable 0
            bool affects_obs0 = false;
            for (uint obs_id : obs_vector) {
                if (obs_id == 0) {
                    affects_obs0 = true;
                    break;
                }
            }
            
            if (affects_obs0) {
                edges_affecting_obs0++;
                log.debug("Edge (", v1, "→", v2, ") affects Observable 0");
                
                // Check if this edge is being used for correction
                auto edge_it = std::find(edge_list.begin(), edge_list.end(), std::make_pair(v1, v2));
                if (edge_it != edge_list.end()) {
                    uint edge_id = edge_it - edge_list.begin();
                    if (edge_states[edge_id] == EdgeSupport::MATCHING) {
                        log.debug("     This edge is ACTIVE in correction!");
                    }
                }
            }
        }
        
        log.debug("Total edges affecting Observable 0: ", edges_affecting_obs0);
        log.debug("Circuit observables: ", circuit.count_observables());
        
        if (circuit.count_observables() == 1) {
            log.debug("     WARNING: Observable 0 is likely the LOGICAL operator!");
            log.debug("     Corrections should NOT flip logical operators!");
        }
    }

    bool UFDecoder::check_logical_error() {
        log.debug("=== LOGICAL ERROR DETECTION (Observable-based) ===");
        
        debug_observable_mapping();
        
        // Generate correction from matching edges
        std::vector<uint8_t> correction(circuit.count_observables(), 0);
        
        for (uint e = 0; e < edge_states.size(); e++) {
            if (edge_states[e] == EdgeSupport::MATCHING) {
                apply_edge_correction(e, correction);
            }
        }
        
        // Check if any observable is affected (simplified logical error check)
        bool logical_error = false;
        for (uint obs = 0; obs < correction.size(); obs++) {
            if (correction[obs] == 1) {
                logical_error = true;
                log.debug("Observable ", obs, " affected by correction");
                break;
            }
        }
        
        // Additional check: isolated odd-parity clusters
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;
            Cluster* root = cluster->find();
            
            // Odd parity cluster that doesn't contain virtual boundary = logical error
            if (root->parity % 2 == 1 && !cluster_contains_virtual_boundary(root)) {
                logical_error = true;
                log.debug("Isolated odd-parity cluster ", root->id, " causes logical error");
            }
        }
        
        log.debug("Logical error detected: ", (logical_error ? "YES" : "NO"));
        log.debug("=================================");
        
        return logical_error;
    }

    void UFDecoder::apply_edge_correction(uint edge_id, std::vector<uint8_t>& correction) {
        if (edge_id >= edge_list.size()) return;
        
        // Get vertices and ensure proper ordering
        uint v1 = edge_list[edge_id].first;
        uint v2 = edge_list[edge_id].second;
        uint min_v = std::min(v1, v2);
        uint max_v = std::max(v1, v2);
        
        // Look up with consistent ordering
        auto it = edge_observables_map.find({min_v, max_v});
        if (it != edge_observables_map.end()) {
            // Add debug log here
            if (!it->second.empty()) {
                std::string obs_str = "";
                for (uint obs_id : it->second) {
                    obs_str += std::to_string(obs_id) + " ";
                }
                log.debug("Correction edge ", v1, " <-> ", v2, " flips observables: ", obs_str);
            }
            
            for (uint obs_id : it->second) {
                if (obs_id < correction.size()) {
                    correction[obs_id] ^= 1;
                }
            }
        }
    }

    void UFDecoder::debug_print_cluster_boundaries(Cluster* cluster) {
        log.debug("=== Cluster ", cluster->id, " Analysis ===");
        log.debug("Vertices: ", cluster->vertices.size());
        log.debug("Parity: ", cluster->parity, " (", (cluster->parity % 2 == 0 ? "EVEN" : "ODD"), ")");
        
        // Only check for virtual boundaries
        bool has_virtual = cluster_contains_virtual_boundary(cluster);
        log.debug("Contains virtual boundary: ", has_virtual);
        
        // Determine correctability
        if (cluster->parity % 2 == 0) {
            log.debug("Status: CORRECTABLE (even parity)");
        } else if (has_virtual) {
            log.debug("Status: CORRECTABLE (odd parity + virtual discharge)");
        } else {
            log.debug("Status: LOGICAL ERROR (isolated odd parity)");
        }
    }

    void UFDecoder::debug_print_surface_code_layout() {
        log.debug("=== Surface Code Layout Analysis ===");
        auto vertices = decoding_graph.vertices();
        
        int virtual_count = 0, physical_count = 0;
        for (uint i = 0; i < vertices.size(); i++) {
            if (vertices[i]->detector == UINT_MAX) {
                virtual_count++;
            } else {
                physical_count++;
            }
        }
        
        log.debug("Total vertices: ", vertices.size(), " (", physical_count, " physical, ", virtual_count, " virtual)");
        log.debug("Virtual boundary vertices: ", vertex_boundaries.size());
    }

    void UFDecoder::debug_print_ascii_grid() {
        log.debug("=== Syndrome Pattern ===");
        
        // Simple syndrome listing instead of complex grid
        std::vector<uint> syndrome_vertices;
        for (uint i = 0; i < vertex_has_syndrome.size(); i++) {
            if (vertex_has_syndrome[i]) {
                syndrome_vertices.push_back(i);
            }
        }
        
        log.debug("Syndrome vertices: ", vector_to_string(syndrome_vertices));
        log.debug("Virtual boundaries: ", vertex_boundaries.size());
    }

    std::string UFDecoder::vector_to_string(const std::vector<uint>& vec) {
        if (vec.empty()) return "[]";
        
        std::string result = "[";
        for (size_t i = 0; i < vec.size(); i++) {
            result += std::to_string(vec[i]);
            if (i < vec.size() - 1) result += ",";
        }
        result += "]";
        return result;
    }

    void UFDecoder::debug_print_raw_graph() {
        log.debug("=== RAW DECODING GRAPH STRUCTURE ===");
        auto vertices = decoding_graph.vertices();
        
        log.debug("Total vertices: ", vertices.size());
        
        // Print vertex information
        for (uint v_idx = 0; v_idx < vertices.size(); v_idx++) {
            auto v = vertices[v_idx];
            log.debug("Vertex ", v_idx, 
                    " [Det=", (v->detector == UINT_MAX ? "virtual" : std::to_string(v->detector)), 
                    ", Coords=", v->coord[0], ",", v->coord[1], ",", v->coord[2], "]");
            
            // Print raw adjacency from decoding graph
            auto adj_list = decoding_graph.adjacency_list(v);
            std::string adj_str = "  Raw adjacency: ";
            
            for (auto neighbor : adj_list) {
                // Find index of neighbor in vertices list
                uint n_idx = std::find(vertices.begin(), vertices.end(), neighbor) - vertices.begin();
                adj_str += std::to_string(n_idx) + "(" + 
                        (neighbor->detector == UINT_MAX ? "V" : std::to_string(neighbor->detector)) + "), ";
            }
            log.debug(adj_str);
        }
    }

    void UFDecoder::print_spanning_tree(const std::map<uint, std::set<uint>>& adjacency, 
                                   const std::set<uint>& remaining_vertices) {
    log.debug("  === SPANNING TREE STRUCTURE (remaining vertices: ", remaining_vertices.size(), ") ===");
    
    // Build a more visual representation with indentation
    std::map<uint, bool> visited;
    
    // Find a good root - preferably a non-leaf vertex
    uint root = *remaining_vertices.begin();
    for (uint v : remaining_vertices) {
        uint degree = 0;
        for (uint neighbor : adjacency.at(v)) {
            if (remaining_vertices.count(neighbor) > 0) {
                degree++;
            }
        }
        if (degree > 1) {
            root = v;
            break;
        }
    }
    
    // Print tree using DFS
    std::function<void(uint, int)> print_subtree = [&](uint vertex, int depth) {
        visited[vertex] = true;
        
        // Build indentation
        std::string indent(depth * 2, ' ');
        
        // Build node info
        std::string node_info = std::to_string(vertex);
        if (vertex_has_syndrome[vertex]) {
            node_info += " [S]"; // Syndrome
        }
        if (vertex_boundaries.count(vertex) > 0) {
            node_info += " [B]"; // Boundary
        }
        
        log.debug(indent, "└─ ", node_info);
        
        // Process children
        for (uint neighbor : adjacency.at(vertex)) {
            if (remaining_vertices.count(neighbor) > 0 && !visited[neighbor]) {
                print_subtree(neighbor, depth + 1);
            }
        }
    };
    
    // Start DFS from root
    print_subtree(root, 0);
    log.debug("  ===============================");
}

void UFDecoder::debug_print_raw_edge_frames() {
    log.debug("=== RAW EDGE FRAMES DEBUG ===");
    auto vertices = decoding_graph.vertices();
    
    int edges_with_frames = 0;
    
    for (uint v1_idx = 0; v1_idx < vertices.size(); v1_idx++) {
        auto v1 = vertices[v1_idx];
        auto adj_list = decoding_graph.adjacency_list(v1);
        
        for (auto v2 : adj_list) {
            uint v2_idx = std::find(vertices.begin(), vertices.end(), v2) - vertices.begin();
            
            // Only process each edge once
            if (v1_idx < v2_idx) {
                auto edge = decoding_graph.get_edge(v1, v2);
                
                if (edge != nullptr && !edge->frames.empty()) {
                    edges_with_frames++;
                    
                    std::string frames_str = "";
                    for (uint frame_id : edge->frames) {
                        frames_str += std::to_string(frame_id) + " ";
                    }
                    
                    log.debug("Edge ", v1_idx, " <-> ", v2_idx, 
                              " [Det: ", (v1->detector == UINT_MAX ? "virtual" : std::to_string(v1->detector)),
                              " <-> ", (v2->detector == UINT_MAX ? "virtual" : std::to_string(v2->detector)), "]",
                              " has raw frames: ", frames_str);
                }
            }
        }
    }
    
    log.debug("Total edges with non-empty frames: ", edges_with_frames);
    log.debug("===============================");
}
} // namespace qrc