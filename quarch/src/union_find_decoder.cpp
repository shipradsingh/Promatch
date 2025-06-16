#include "union_find_decoder.h"
#include "decoding_graph.h"
#include "logger.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <set>
#include <functional>
#include <map>
#include <vector>
#include <thread>

using namespace std;

namespace qrc {

    UFDecoder::UFDecoder(const stim::Circuit& circuit)
        : Decoder(circuit), cluster_id(0) {
        
        n_detectors = circuit.count_detectors();
        
        log.debug("Building Union-Find decoder using DecodingGraph...");
        decoding_graph = to_decoding_graph(circuit);

        // OPTIONAL: Uncomment to print debugging information
        //debug_print_raw_graph();
        //debug_print_raw_edge_frames();
        
        log.debug("Decoding graph built with ", decoding_graph.count_detectors(), " detectors.");
        debug_print_raw_graph();
        log.debug("DECODER SETUP:");
        log.debug("  n_detectors = ", n_detectors);
        log.debug("  graph_vertices = ", decoding_graph.vertices().size());
        
        setup_boundaries_from_decoding_graph();
        build_edge_mapping();
        
        vertex_to_cluster.resize(decoding_graph.vertices().size(), nullptr);
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
        detector_to_vertex_id.clear();
        detector_to_vertex_id.resize(n_detectors, UINT_MAX);
        
        for (uint v = 0; v < vertices.size(); v++) {
            if (vertices[v]->detector < n_detectors) {
                detector_to_vertex_id[vertices[v]->detector] = v;
            }
        }
        for (uint v1 = 0; v1 < vertices.size(); v1++) {
            auto adj_list = decoding_graph.adjacency_list(vertices[v1]);
            for (auto v2_ptr : adj_list) {
                uint v2 = find(vertices.begin(), vertices.end(), v2_ptr) - vertices.begin();
                if (v1 < v2) {
                    edge_list.push_back({v1, v2});
                    
                    // Set edge weights: temporal = 0.0, spatial = 1.0
                    double edge_weight = 1.0;
                    auto edge = decoding_graph.get_edge(vertices[v1], vertices[v2]);
                    
                    if (edge != nullptr) {
                        bool is_temporal = (v2_ptr == decoding_graph.get_next_round(vertices[v1]));
                        if (is_temporal) {
                            edge_weight = 0.0;
                            log.debug("TEMPORAL EDGE: ", v1, " -> ", v2,
                                    " between detectors ", vertices[v1]->detector, " and ", vertices[v2]->detector);
                        }
                        vector<uint> observables;
                        for (uint obs_id : edge->frames) {
                            observables.push_back(obs_id);
                        }
                        uint min_v = min(v1, v2);
                        uint max_v = max(v1, v2);
                        edge_observables_map[{min_v, max_v}] = observables;
                    } else {
                        uint min_v = min(v1, v2);
                        uint max_v = max(v1, v2);
                        edge_observables_map[{min_v, max_v}] = {};
                    }
                    
                    edge_weights.push_back(edge_weight);
                }
            }
        }
        
        edge_states.resize(edge_list.size(), EdgeSupport::NONE);
        
        log.debug("Built edge mapping: ", edge_list.size(), " edges");

        // Optional: Print adjacency matrix and edge observables for debugging
        //debug_print_adjacency_matrix();
        //debug_print_edge_observables();
    }

    DecoderShotResult UFDecoder::decode_error(const vector<uint8_t>& syndrome) {
        // Convert uint8_t vector to uint vector for logging
        vector<uint> syndrome_uint(syndrome.begin(), syndrome.end());
        string syndrome_str = vector_to_string(syndrome_uint);
        log.info("Input Syndrome: ", syndrome_str);
        // Clear previous state - delete ALL clusters from previous run
        for (auto cluster : clusters) {
            if (cluster != nullptr) {
                delete cluster;
            }
        }
        clusters.clear();
        cluster_id = 0;
        
        fill(vertex_to_cluster.begin(), vertex_to_cluster.end(), nullptr);
        vertex_has_syndrome.assign(decoding_graph.vertices().size(), false);
        edge_states.assign(edge_list.size(), EdgeSupport::NONE);
        
        extract_vertex_syndrome(syndrome);
        auto total_start = chrono::high_resolution_clock::now();
        
        // Phase 1: Initialize singleton clusters
        auto phase1_start = chrono::high_resolution_clock::now();
        initialize_clusters();
        auto phase1_end = chrono::high_resolution_clock::now();
        
        // Phase 2: Grow and merge clusters until frozen
        auto phase2_start = chrono::high_resolution_clock::now();
        grow_and_merge_clusters();
        auto phase2_end = chrono::high_resolution_clock::now();
        
        // Phase 3: Extract corrections from clusters
        auto phase3_start = chrono::high_resolution_clock::now();
        peel_clusters();
        auto phase3_end = chrono::high_resolution_clock::now();
        
        auto total_end = chrono::high_resolution_clock::now();
        
        // Calculate phase times in nanoseconds
        auto phase1_time = chrono::duration_cast<chrono::nanoseconds>(phase1_end - phase1_start).count();
        auto phase2_time = chrono::duration_cast<chrono::nanoseconds>(phase2_end - phase2_start).count();
        auto phase3_time = chrono::duration_cast<chrono::nanoseconds>(phase3_end - phase3_start).count();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(total_end - total_start).count();
        
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

    void UFDecoder::extract_vertex_syndrome(const vector<uint8_t>& detector_syndrome) {
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
        
        for (uint v = 1; v < vertex_has_syndrome.size(); v++) {
            if (vertex_has_syndrome[v]) {
                Cluster* cluster = new Cluster(cluster_id++);
                cluster->parity = 0; 
                cluster_add_vertex(cluster, v);
                
                // Initialize growth progress for the syndrome vertex
                cluster->growth_progress[v] = 0.0;
                
                clusters.push_back(cluster);
                vertex_to_cluster[v] = cluster;
                
                log.debug("Created syndrome cluster ", cluster->id, " with initial boundary vertex ", v);
            }
        }
        
        // Virtual boundary cluster
        if (vertex_boundaries.count(0) > 0) {
            Cluster* boundary_cluster = new Cluster(cluster_id++);
            boundary_cluster->parity = 0;
            boundary_cluster->is_virtual_boundary = true;
            cluster_add_vertex(boundary_cluster, 0);
            clusters.push_back(boundary_cluster);
            vertex_to_cluster[0] = boundary_cluster;
            log.debug("Created virtual boundary cluster ", boundary_cluster->id);
        }
    }

    void UFDecoder::grow_and_merge_clusters() {
        log.debug("Phase 2: Growing and merging clusters");
        volatile bool changes = true;
        int iteration = 1;
        double growth_unit = 0.5;
        
        // Initialize growth progress
        for (auto& cluster : clusters) {
            if (cluster) {
                for (uint v : cluster->vertices) {
                    cluster->growth_progress[v] = 0.0;
                }
            }
        }
        
        while (changes && iteration < 1000) {
            // Use index-based loop instead of range-based loop
            for (size_t i = 0; i < clusters.size(); ++i) {
                auto& cluster = clusters[i];
                
                if (!cluster || cluster->is_frozen || cluster->is_virtual_boundary) {
                    continue;
                }
                
                vector<uint> cluster_boundary = get_cluster_boundary(cluster);
                log.debug("Cluster ", cluster->id, " boundary vertices: ", vector_to_string(cluster_boundary));

                for (auto& v : cluster_boundary) {
                    if (vertex_to_cluster[v] == cluster) {
                        cluster->growth_progress[v] += growth_unit;
                        log.debug("Cluster ", cluster->id, " vertex ", v, " growth progress: ", 
                                cluster->growth_progress[v]);
                    }
                }
                
                process_cluster_boundary(cluster);
                freeze_completed_clusters();
                
                // After merge_clusters(), clusters vector might have changed size
                // The index-based loop handles this gracefully
            }
            
            // Check if any clusters remain unfrozen
            changes = false;
            iteration++;
            for (auto& cluster : clusters) {
                if (cluster && !cluster->is_frozen && !cluster->is_virtual_boundary) {
                    changes = true;
                    break;
                }
            }
        }

        if (iteration >= 1000) {
            log.error("Reached maximum iterations in grow_and_merge_clusters, potential infinite loop detected!");
        }
        log.debug("Phase 2 complete after ", iteration, " iterations");
    }

    vector<uint> UFDecoder::get_cluster_boundary(Cluster* cluster) {
        vector<uint> boundary_vertices;
        
        for (uint vertex_id : cluster->vertices) {
            // Check if this vertex has at least one neighbor with ungrown edge
            bool is_boundary = false;
            
            for (uint e = 0; e < edge_list.size(); e++) {
                if (edge_states[e] != EdgeSupport::NONE) continue;  // Skip grown edges
                
                uint v1 = edge_list[e].first;
                uint v2 = edge_list[e].second;
                
                // Check if this vertex is connected to an ungrown edge
                if (v1 == vertex_id || v2 == vertex_id) {
                    is_boundary = true;
                    break;
                }
            }
            
            if (is_boundary) {
                boundary_vertices.push_back(vertex_id);
            }
        }
        return boundary_vertices;
    }

    void UFDecoder::process_cluster_boundary(Cluster* cluster) {
        if (cluster->is_frozen) return;
        
        vector<uint> boundary_vertices = get_cluster_boundary(cluster);
        
        // Iterate over all boundary vertices
        for (uint boundary_vertex : boundary_vertices) {
            double growth_progress = cluster->growth_progress[boundary_vertex];
            
            // Check all edges connected to this boundary vertex
            for (uint e = 0; e < edge_list.size(); e++) {
                if (edge_states[e] != EdgeSupport::NONE) continue;  // Skip grown edges
                
                uint v1 = edge_list[e].first;
                uint v2 = edge_list[e].second;
                double edge_weight = edge_weights[e];
                
                uint neighbor_vertex = UINT_MAX;
                if (v1 == boundary_vertex) {
                    neighbor_vertex = v2;
                } else if (v2 == boundary_vertex) {
                    neighbor_vertex = v1;
                } else {
                    continue;  // This edge doesn't involve our boundary vertex
                }
                
                Cluster* neighbor_cluster = vertex_to_cluster[neighbor_vertex];
                
                // Case 1: Cluster-to-cluster merge (both need edge_weight/2.0)
                if (neighbor_cluster && neighbor_cluster->find() != cluster->find() && 
                    !neighbor_cluster->find()->is_virtual_boundary) {
                    
                    double required_growth = edge_weight / 2.0;  // Each cluster grows halfway
                    
                    if (growth_progress >= required_growth) {
                        // Check if neighbor cluster can also reach halfway
                        bool neighbor_can_reach = false;
                        if (neighbor_cluster->growth_progress.count(neighbor_vertex) > 0 &&
                            neighbor_cluster->growth_progress[neighbor_vertex] >= required_growth) {
                            neighbor_can_reach = true;
                        }
                        
                        if (neighbor_can_reach) {
                            log.debug("MERGE DETECTED: Clusters meet at midpoint - edge weight=", edge_weight,
                                    ", each needs=", required_growth);
                            merge_clusters(cluster, neighbor_cluster, e);
                            return;
                        }
                    }
                }
                
                // Case 2: Virtual boundary connection (cluster needs full edge_weight)
                else if (neighbor_cluster && neighbor_cluster->find()->is_virtual_boundary) {
                    if (growth_progress >= edge_weight) {  // Full weight needed for virtual boundary
                        log.debug("VIRTUAL BOUNDARY DETECTED: Cluster ", cluster->id, 
                                " reached virtual boundary (growth=", growth_progress, 
                                ", edge_weight=", edge_weight, ")");
                        merge_clusters(cluster, neighbor_cluster, e);
                        return;
                    }
                }
                
                // Case 3: Free vertex (cluster needs full edge_weight)
                else if (!neighbor_cluster) {
                    if (growth_progress >= edge_weight) {  // Full weight needed for free vertex
                        add_free_vertex_to_cluster(cluster, neighbor_vertex, e);
                        log.debug("FREE VERTEX: Cluster ", cluster->id, " absorbed vertex ", neighbor_vertex);
                    }
                }
            }
            
            if (cluster->is_frozen) break;
        }
    }

    void UFDecoder::add_free_vertex_to_cluster(Cluster* cluster, uint free_vertex, uint edge_id) {
        // Add vertex to cluster
        cluster_add_vertex(cluster, free_vertex);
        cluster->add_tree_edge(edge_id);
        edge_states[edge_id] = EdgeSupport::GROWN;
        vertex_to_cluster[free_vertex] = cluster;
        
        // Initialize growth progress for new vertex (it becomes a potential boundary)
        cluster->growth_progress[free_vertex] = 0.0;
        
        log.debug("Cluster ", cluster->id, " added free vertex ", free_vertex);
    }

    void UFDecoder::merge_clusters(Cluster* c1, Cluster* c2, uint edge_id) {
        Cluster* root1 = c1->find();
        Cluster* root2 = c2->find();
        
        log.debug("Merging clusters ", root1->id, " and ", root2->id);
        
        // Merge all content from root2 into root1
        for (uint v : root2->vertices) {
            root1->vertices.push_back(v);
            vertex_to_cluster[v] = root1;
        }
        for (uint e : root2->tree_edges) {
            root1->tree_edges.push_back(e);
        }
        root1->tree_edges.push_back(edge_id);
        
        // Update properties
        root1->parity ^= root2->parity;
        root1->is_virtual_boundary = root1->is_virtual_boundary || root2->is_virtual_boundary;
        root1->is_frozen = root1->is_frozen || root2->is_frozen;
        
        // Freeze if virtual boundary OR even parity
        if (root1->is_virtual_boundary || root1->parity % 2 == 0 || root1->is_frozen) {
            root1->is_frozen = true;
            root1->growth_progress.clear();
        }
        
        // Cleanup
        edge_states[edge_id] = EdgeSupport::GROWN;
        root2->parent = root1;
        clusters.erase(find(clusters.begin(), clusters.end(), root2));
        
        log.debug("Merge result: parity=", root1->parity, " frozen=", root1->is_frozen);
    }

    void UFDecoder::cluster_add_vertex(Cluster* cluster, uint vertex_id) {
        cluster->add_vertex(vertex_id);
        // Only increment parity if vertex has syndrome
        if (vertex_id < vertex_has_syndrome.size() && vertex_has_syndrome[vertex_id]) {
            cluster->parity ^= 1;
        }
        vertex_to_cluster[vertex_id] = cluster;
    }

    vector<SpanningTree> construct_spanning_forest(const map<uint, set<uint>>& adjacency,
                                                        const set<uint>& vertices) {
        vector<SpanningTree> forest;
        set<uint> visited;
        
        // DFS lambda to build a tree starting from 'start'.
        function<void(uint, SpanningTree&, uint)> dfs =
        [&](uint current, SpanningTree &tree, uint parent_node) {
            visited.insert(current);
            if (parent_node != UINT_MAX) {
                tree.parent[current] = parent_node;
            }
            // Traverse each neighbor in the adjacency (only within provided vertices)
            for (uint neighbor : adjacency.at(current)) {
                if (vertices.count(neighbor) && visited.find(neighbor) == visited.end()) {
                    dfs(neighbor, tree, current);
                }
            }
        };
        
        // For each vertex not visited, start a new spanning tree.
        for (uint v : vertices) {
            if (visited.find(v) == visited.end()) {
                SpanningTree tree;
                
                // Choose root: virtual vertex (0) if present, otherwise a leaf
                uint chosen_root = v;  // default
                
                // Priority 1: Look for virtual vertex (vertex 0)
                if (vertices.count(0) > 0) {
                    chosen_root = 0;
                } else {
                    // Priority 2: Look for a leaf (degree = 1)
                    for (uint candidate : vertices) {
                        if (visited.find(candidate) == visited.end()) {
                            uint degree = 0;
                            for (uint neighbor : adjacency.at(candidate)) {
                                if (vertices.count(neighbor)) {
                                    degree++;
                                }
                            }
                            if (degree == 1) {
                                chosen_root = candidate;
                                break;
                            }
                        }
                    }
                }
                
                tree.root = chosen_root;
                dfs(chosen_root, tree, UINT_MAX);
                forest.push_back(tree);
            }
        }
        return forest;
    }

    void UFDecoder::peel_clusters() {
        log.debug("Phase 3: Peeling clusters to recover the error pattern");
        uint clusters_peeled = 0;
        uint total_peel_iterations = 0;
        
        for (auto& cluster : clusters) {
            if (cluster == nullptr)
                continue;  // Skip nullified clusters
            Cluster* rootCluster = cluster->find();
            
            log.debug("Peeling cluster ", rootCluster->id);
            clusters_peeled++;
            
            // Build edge lookup table and adjacency from stored tree_edges.
            map<pair<uint, uint>, uint> vertex_pair_to_edge;
            map<uint, set<uint>> adjacency;
            set<uint> cluster_vertices(rootCluster->vertices.begin(), rootCluster->vertices.end());
            for (uint vertex_id : cluster_vertices) {
                adjacency[vertex_id] = {};
            }
            
            for (uint edge_id : rootCluster->tree_edges) {
                uint v1 = edge_list[edge_id].first;
                uint v2 = edge_list[edge_id].second;
                if (cluster_vertices.count(v1) && cluster_vertices.count(v2)) {
                    adjacency[v1].insert(v2);
                    adjacency[v2].insert(v1);
                    uint min_v = min(v1, v2);
                    uint max_v = max(v1, v2);
                    vertex_pair_to_edge[{min_v, max_v}] = edge_id;
                }
            }
            
            // Instead of one tree, create a spanning forest.
            vector<SpanningTree> forest = construct_spanning_forest(adjacency, cluster_vertices);
            print_spanning_forest(forest, adjacency, cluster_vertices);
            
            // Now use backward peeling over the spanning forest.
            peel_spanning_forest(forest, vertex_pair_to_edge, cluster_vertices);
        }
        log.debug("Peeling complete: ", clusters_peeled, " clusters processed");
    }

    void UFDecoder::peel_spanning_forest(const vector<SpanningTree>& forest,
                                        const map<pair<uint, uint>, uint>& vertex_pair_to_edge,
                                        set<uint>& cluster_vertices) {
        // Copy the current syndrome state.
        map<uint, bool> syndrome_status;
        for (uint v : cluster_vertices) {
            syndrome_status[v] = (v < vertex_has_syndrome.size() && vertex_has_syndrome[v]);
        }
        
        // Build parent and children mappings from the forest.
        map<uint, uint> parent_mapping;  // child -> parent
        map<uint, vector<uint>> children;
        for (const auto &tree : forest) {
            // tree.parent maps child -> parent.
            for (const auto &entry : tree.parent) {
                uint child = entry.first;
                uint parent = entry.second;
                parent_mapping[child] = parent;
                children[parent].push_back(child);
            }
        }
        
        // Build a children count map and initialize a queue with leaves.
        queue<uint> q;
        map<uint, int> childCount;
        for (uint v : cluster_vertices) {
            // If v is not in children map, children[v] returns an empty vector.
            childCount[v] = children[v].size();
            if (childCount[v] == 0) {
                q.push(v);
            }
        }
        
        // Process the queue.
        while (!q.empty()) {
            uint leaf = q.front();
            q.pop();
            
            // If the leaf has no syndrome, don't process its edge.
            if (!syndrome_status[leaf]) {
                if (parent_mapping.find(leaf) != parent_mapping.end()) {
                    uint par = parent_mapping[leaf];
                    childCount[par]--;
                    if (childCount[par] == 0) {
                        q.push(par);
                    }
                }
                continue;
            }
            
            // If the leaf has syndrome and has a parent, match the edge.
            if (parent_mapping.find(leaf) != parent_mapping.end()) {
                uint par = parent_mapping[leaf];
                uint min_v = min(par, leaf);
                uint max_v = max(par, leaf);
                auto it = vertex_pair_to_edge.find({min_v, max_v});
                if (it != vertex_pair_to_edge.end()) {
                    uint edge_id = it->second;
                    edge_states[edge_id] = EdgeSupport::MATCHING;
                    log.debug("Peeling BFS: Matching edge ", par, " <-> ", leaf, " (edge id ", edge_id, ")");
                }
                // Flip parent's syndrome and clear leaf's syndrome.
                syndrome_status[par] = !syndrome_status[par];
                syndrome_status[leaf] = false;
                
                // Process parent's children.
                childCount[par]--;
                if (childCount[par] == 0) {
                    q.push(par);
                }
            }
        }
    }

    bool UFDecoder::should_freeze_cluster(Cluster* cluster) {
        Cluster* root = cluster->find();        
        // Freeze if: virtual boundary OR even parity OR already frozen
        return (root->is_virtual_boundary || 
                root->parity % 2 == 0 || 
                root->is_frozen);
    }

    void UFDecoder::freeze_completed_clusters() {
        for (auto& cluster : clusters) {
            // Skip nullified clusters and virtual boundaries
            if (cluster == nullptr || cluster->is_virtual_boundary) continue;
            
            Cluster* root = cluster->find();
            if (root->is_frozen) continue;
            
            if (should_freeze_cluster(root)) {
                root->is_frozen = true;
                log.debug("Froze cluster ", root->id, " (parity=", root->parity, ")");
            }
        }
    }

    void UFDecoder::force_freeze_all_clusters() {
        for (auto& cluster : clusters) {
            if (cluster) {
                cluster->find()->is_frozen = true;
            }
        }
    }

    string UFDecoder::vector_to_string(const vector<uint>& vec) {
        if (vec.empty()) return "[]";
        
        string result = "[";
        for (size_t i = 0; i < vec.size(); i++) {
            result += to_string(vec[i]);
            if (i < vec.size() - 1) result += ",";
        }
        result += "]";
        return result;
    }

    // Add this helper function to UFDecoder class
    void UFDecoder::print_spanning_tree(const map<uint, set<uint>>& adjacency, 
                                const set<uint>& remaining_vertices) {
        log.debug("  === SPANNING TREE STRUCTURE (remaining vertices: ", remaining_vertices.size(), ") ===");
        
        // Build a more visual representation with indentation
        map<uint, bool> visited;
        
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
        function<void(uint, int)> print_subtree = [&](uint vertex, int depth) {
            visited[vertex] = true;
            
            // Build indentation
            string indent(depth * 2, ' ');
            
            // Build node info
            string node_info = to_string(vertex);
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

    void UFDecoder::print_spanning_forest(const vector<SpanningTree>& forest,
                                            const map<uint, set<uint>>& adjacency,
                                            const set<uint>& cluster_vertices) {
        log.debug("=== SPANNING FOREST ===");
        // Process each tree in the forest separately.
        for (const auto &tree : forest) {
            log.debug("Tree with root ", tree.root, ":");
            map<uint, bool> visited;
            // DFS lambda to print subtree structure starting at 'current'.
            function<void(uint, int)> print_subtree = [&](uint current, int depth) {
                visited[current] = true;
                string indent(depth * 2, ' ');
                string node_info = to_string(current);
                if (current < vertex_has_syndrome.size() && vertex_has_syndrome[current]) {
                    node_info += " [S]"; // Mark syndrome vertices.
                }
                if (vertex_boundaries.count(current) > 0) {
                    node_info += " [B]"; // Mark boundary vertices.
                }
                log.debug(indent, "└─ ", node_info);
                // Traverse neighbors that are still in the same connected component.
                for (uint neighbor : adjacency.at(current)) {
                    if (cluster_vertices.count(neighbor) && !visited[neighbor]) {
                        print_subtree(neighbor, depth + 1);
                    }
                }
            };
            print_subtree(tree.root, 0);
        }
        log.debug("=== END OF SPANNING FOREST ===");
    }

    void UFDecoder::apply_edge_correction(uint edge_id, vector<uint8_t>& correction) {
        uint v1 = edge_list[edge_id].first;
        uint v2 = edge_list[edge_id].second;
        
        // Use consistent ordering for lookup
        uint min_v = min(v1, v2);
        uint max_v = max(v1, v2);
        
        // Find the observables affected by this edge
        auto it = edge_observables_map.find({min_v, max_v});
        if (it != edge_observables_map.end()) {
            // Flip each observable affected by this edge
            for (uint obs_id : it->second) {
                if (obs_id < correction.size()) {
                    correction[obs_id] ^= 1;  // XOR flip
                    log.debug("Correction edge ", v1, " <-> ", v2, " flips observables: ", obs_id, " ");
                }
            }
        }
    }

    void UFDecoder::debug_print_raw_graph() {
        log.debug("=== RAW DECODING GRAPH STRUCTURE ===");
        auto vertices = decoding_graph.vertices();
        
        log.debug("Total vertices: ", vertices.size());
        
        // Print vertex information
        for (uint v_idx = 0; v_idx < vertices.size(); v_idx++) {
            auto v = vertices[v_idx];
            log.debug("Vertex ", v_idx, 
                    " [Det=", (v->detector == UINT_MAX ? "virtual" : to_string(v->detector)), 
                    ", Coords=", v->coord[0], ",", v->coord[1], ",", v->coord[2], "]");
            
            // Print raw adjacency from decoding graph
            auto adj_list = decoding_graph.adjacency_list(v);
            string adj_str = "  Raw adjacency: ";
            
            for (auto neighbor : adj_list) {
                // Find index of neighbor in vertices list
                uint n_idx = find(vertices.begin(), vertices.end(), neighbor) - vertices.begin();
                adj_str += to_string(n_idx) + "(" + 
                        (neighbor->detector == UINT_MAX ? "V" : to_string(neighbor->detector)) + "), ";
            }
            log.debug(adj_str);
        }
    }

    void UFDecoder::debug_print_adjacency_matrix() {
        log.debug("=== Graph Adjacency Matrix ===");
        auto vertices = decoding_graph.vertices();
        for (uint v1 = 0; v1 < vertices.size(); v1++) {
            string connections = "Vertex " + to_string(v1) + " connects to: ";
            for (uint v2 = 0; v2 < vertices.size(); v2++) {
                bool connected = false;
                for (uint e = 0; e < edge_list.size(); e++) {
                    if ((edge_list[e].first == v1 && edge_list[e].second == v2) ||
                        (edge_list[e].first == v2 && edge_list[e].second == v1)) {
                        double weight = edge_weights[e];
                        bool is_temporal = (weight == 0.0);
                        
                        connections += to_string(v2) + "(" + 
                                    (is_temporal ? "T" : "S") + "), ";
                        connected = true;
                        break;
                    }
                }
            }
            log.debug(connections);
        }
    }

    void UFDecoder::debug_print_raw_edge_frames() {
        log.debug("=== RAW EDGE FRAMES DEBUG ===");
        auto vertices = decoding_graph.vertices();
        
        int edges_with_frames = 0;
        
        for (uint v1_idx = 0; v1_idx < vertices.size(); v1_idx++) {
            auto v1 = vertices[v1_idx];
            auto adj_list = decoding_graph.adjacency_list(v1);
            
            for (auto v2 : adj_list) {
                uint v2_idx = find(vertices.begin(), vertices.end(), v2) - vertices.begin();
                
                // Only process each edge once
                if (v1_idx < v2_idx) {
                    auto edge = decoding_graph.get_edge(v1, v2);
                    
                    if (edge != nullptr && !edge->frames.empty()) {
                        edges_with_frames++;
                        
                        string frames_str = "";
                        for (uint frame_id : edge->frames) {
                            frames_str += to_string(frame_id) + " ";
                        }
                        
                        log.debug("Edge ", v1_idx, " <-> ", v2_idx, 
                                " [Det: ", (v1->detector == UINT_MAX ? "virtual" : to_string(v1->detector)),
                                " <-> ", (v2->detector == UINT_MAX ? "virtual" : to_string(v2->detector)), "]",
                                " has raw frames: ", frames_str);
                    }
                }
            }
        }
        
        log.debug("Total edges with non-empty frames: ", edges_with_frames);
        log.debug("===============================");
    }

    void UFDecoder::debug_print_edge_observables() {
        log.debug("=== EDGE OBSERVABLES MAP ===");
        log.debug("Total edges with observable effects: ", edge_observables_map.size());
        
        for (const auto& [edge_pair, obs_vec] : edge_observables_map) {
            uint v1 = edge_pair.first;
            uint v2 = edge_pair.second;
            
            if (!obs_vec.empty()) {
                string obs_str = "Observables: ";
                for (uint obs_id : obs_vec) {
                    obs_str += to_string(obs_id) + " ";
                }
                
                log.debug("Edge ", v1, " <-> ", v2, " affects ", obs_vec.size(), 
                        " observables: ", obs_str);
            }
        }
        log.debug("============================");
    }
}