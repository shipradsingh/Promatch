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
        
        // Find coordinate boundaries for logical error detection
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        
        // First pass: find min/max coordinates
        for (uint i = 0; i < vertices.size(); i++) {
            const auto& vertex = vertices[i];
            if (vertex->coord.size() >= 1) {
                double x = static_cast<double>(vertex->coord[0]);
                min_x = std::min(min_x, x);
                max_x = std::max(max_x, x);
            }
        }
        
        // Second pass: mark boundary vertices
        for (uint i = 0; i < vertices.size(); i++) {
            const auto& vertex = vertices[i];
            
            // Mark virtual boundary vertices
            if (vertex->detector == UINT_MAX) {
                vertex_boundaries[i] = 0; // Virtual boundary
                continue;
            }
            
            // Mark coordinate boundaries
            if (vertex->coord.size() >= 1) {
                double x = static_cast<double>(vertex->coord[0]);
                if (x == min_x) {
                    vertex_boundaries[i] = 1; // Left boundary
                } else if (x == max_x) {
                    vertex_boundaries[i] = 2; // Right boundary
                }
            }
        }
        
        log.debug("Boundary setup complete. Found ", vertex_boundaries.size(), " boundary vertices");
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
                        }
                        
                        // Store observables for correction application
                        std::vector<uint> observables;
                        for (uint obs_id : edge->frames) {
                            observables.push_back(obs_id);
                        }
                        edge_observables_map[{v1, v2}] = observables;
                    } else {
                        edge_observables_map[{v1, v2}] = {};
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
    }

    DecoderShotResult UFDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
        current_instance++;
        
        // Clear previous state
        for (auto cluster : clusters) {
            delete cluster;
        }
        clusters.clear();
        cluster_id = 0;
        
        std::fill(vertex_to_cluster.begin(), vertex_to_cluster.end(), nullptr);
        vertex_has_syndrome.assign(decoding_graph.vertices().size(), false);
        edge_states.assign(edge_list.size(), EdgeSupport::NONE);
        
        extract_vertex_syndrome(syndrome);
        
        // Phase 1: Initialize singleton clusters
        initialize_clusters();
        
        // Phase 2: Grow and merge clusters until frozen
        grow_and_merge_clusters();
        
        // Phase 3: Extract corrections from clusters
        peel_clusters();
        
        // Generate final correction
        DecoderShotResult result;
        result.correction.resize(circuit.count_observables(), 0);
        
        for (uint e = 0; e < edge_states.size(); e++) {
            if (edge_states[e] == EdgeSupport::MATCHING) {
                apply_edge_correction(e, result.correction);
            }
        }
        
        result.is_logical_error = check_logical_error();
        
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
        log.info("Phase 1: Initializing singleton clusters");
        
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
        
        log.info("Initialized ", clusters.size(), " singleton clusters");
    }

    void UFDecoder::grow_and_merge_clusters() {
        log.info("Phase 2: Growing and merging clusters");
        
        bool changes = true;
        int iteration = 0;
        
        while (changes) {
            changes = false;
            iteration++;
            log.debug("=== Growth iteration ", iteration, " ===");
            
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
            
            // Process edges in weight order
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
                
                // Case 1: One endpoint has cluster, other is free -> grow cluster
                if ((c1 && !c2) || (!c1 && c2)) {
                    Cluster* active_cluster = c1 ? c1 : c2;
                    uint free_vertex = c1 ? v2 : v1;
                    
                    if (!is_cluster_frozen(active_cluster)) {
                        cluster_add_vertex(active_cluster, free_vertex);
                        active_cluster->add_tree_edge(e);  // ← Track the edge
                        edge_states[e] = EdgeSupport::GROWN;
                        vertex_to_cluster[free_vertex] = active_cluster;
                        changes = true;
                        
                        log.debug("Cluster ", active_cluster->id, " grew to vertex ", free_vertex);
                    }
                }
                // Case 2: Both endpoints have different clusters -> merge
                else if (c1 && c2 && c1->find() != c2->find()) {
                    if (!is_cluster_frozen(c1) && !is_cluster_frozen(c2)) {
                        Cluster* root1 = c1->find();
                        Cluster* root2 = c2->find();
                        
                        log.debug("Merging clusters ", root1->id, " (parity=", root1->parity, 
                                ") and ", root2->id, " (parity=", root2->parity, ")");
                        
                        // IMPORTANT: Merge tree edges from root2 into root1
                        for (uint edge_id : root2->tree_edges) {
                            root1->tree_edges.push_back(edge_id);
                        }
                        
                        // Merge parity and union
                        root1->parity = root1->parity ^ root2->parity;
                        root1->union_with(root2);
                        root1->add_tree_edge(e);  // Add the merging edge
                        
                        // IMPORTANT: Remove root2 from clusters vector
                        auto it = std::find(clusters.begin(), clusters.end(), root2);
                        if (it != clusters.end()) {
                            clusters.erase(it);
                        }
                        
                        // IMPORTANT: Delete root2 - UNCOMMENT THIS LINE:
                        delete root2;  // ← This is the fix!
                        root2 = nullptr;
                        
                        // Update vertex mappings to point to final root
                        Cluster* final_root = root1->find();
                        for (uint v : final_root->vertices) {
                            vertex_to_cluster[v] = final_root;
                        }
                        
                        edge_states[e] = EdgeSupport::GROWN;
                        changes = true;
                        
                        log.debug("Merged result: cluster ", final_root->id, " parity=", final_root->parity);
                        log.debug("Removed cluster ", root2->id, " from clusters vector");
                    }
                }
            }
            
            if (iteration > 10000000) {
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
        
        log.info("Growth complete after ", iteration, " iterations");
    }

    void UFDecoder::peel_clusters() {
        log.info("Phase 3: Peeling clusters to recover the error pattern");
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;
            Cluster* root = cluster->find();
            if (root->instance != current_instance) continue;
            
            log.debug("Peeling cluster ", root->id);
            
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
            
            while (remaining_vertices.size() > 1) {
                // Find leaves
                std::vector<uint> leaves;
                for (uint vertex_id : remaining_vertices) {
                    uint degree = 0;
                    for (uint neighbor : adjacency[vertex_id]) {
                        if (remaining_vertices.count(neighbor) > 0) {
                            degree++;
                        }
                    }
                    if (degree == 1) {
                        leaves.push_back(vertex_id);
                    }
                }
                
                if (leaves.empty()) break;
                
                // Process leaves
                for (uint leaf : leaves) {
                    // Check syndrome and mark corrections
                    bool leaf_has_syndrome = (leaf < vertex_has_syndrome.size() && 
                                            vertex_has_syndrome[leaf]);
                    
                    if (leaf_has_syndrome) {
                        // Find neighbor and mark edge for correction
                        for (uint neighbor : adjacency[leaf]) {
                            if (remaining_vertices.count(neighbor) > 0) {
                                auto edge_it = vertex_pair_to_edge.find({std::min(leaf, neighbor), std::max(leaf, neighbor)});
                                if (edge_it != vertex_pair_to_edge.end()) {
                                    edge_states[edge_it->second] = EdgeSupport::MATCHING;
                                }
                                break; // Only one neighbor for leaf
                            }
                        }
                    }
                    
                    // Remove leaf REGARDLESS of syndrome status
                    remaining_vertices.erase(leaf);  // ← Move this outside the if statement
                }
            }
        }
    }
    
    bool UFDecoder::cluster_touches_virtual_boundary(Cluster* cluster) {
        auto vertices = decoding_graph.vertices();
        
        for (uint vertex_id : cluster->vertices) {
            if (vertex_id >= vertices.size()) continue;
            
            // Check if vertex is virtual boundary
            if (vertices[vertex_id]->detector == UINT_MAX) {
                return true;
            }
            
            // Check adjacency to virtual boundary
            auto adj_list = decoding_graph.adjacency_list(vertices[vertex_id]);
            for (auto adj : adj_list) {
                if (adj->detector == UINT_MAX) {
                    return true;
                }
            }
            
            // Check coordinate boundaries
            auto it = vertex_boundaries.find(vertex_id);
            if (it != vertex_boundaries.end()) {
                return true;  // Any boundary contact triggers freezing
            }
        }
        
        return false;
    }

    // Freezing logic as documented
    bool UFDecoder::should_freeze_cluster(Cluster* cluster) {
        Cluster* root = cluster->find();
        
        // Freeze condition 1: Even parity (all defects paired)
        if (root->parity % 2 == 0) {
            return true;
        }
        
        // Freeze condition 2: Touches virtual boundary
        if (cluster_touches_virtual_boundary(root)) {
            return true;
        }
        
        return false;
    }

    void UFDecoder::freeze_completed_clusters() {
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;
            
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

    bool UFDecoder::check_logical_error() {
        log.debug("=== LOGICAL ERROR DETECTION ===");
        
        // Print overall layout once
        debug_print_surface_code_layout();
        debug_print_ascii_grid();
        
        bool logical_error_found = false;
        
        for (auto& cluster : clusters) {
            if (cluster == nullptr) continue;
            Cluster* root = cluster->find();
            
            log.debug("\n--- Analyzing Cluster ", root->id, " ---");
            debug_print_cluster_boundaries(root);
            
            // Only odd-parity clusters can cause logical errors
            if (root->parity % 2 == 0) {
                log.debug("  SKIP: Even parity cluster cannot cause logical error");
                continue;
            }
            
            log.debug("  CHECKING: Odd parity cluster for logical error...");
            
            // Check if cluster spans left-right boundaries (X logical error)
            bool has_left = false, has_right = false;
            
            for (uint vertex_id : root->vertices) {
                auto it = vertex_boundaries.find(vertex_id);
                if (it != vertex_boundaries.end()) {
                    if (it->second == 1) has_left = true;   // Left boundary
                    if (it->second == 2) has_right = true;  // Right boundary
                }
            }
            
            if (has_left && has_right) {
                log.debug("LOGICAL ERROR: Cluster spans left-right boundaries (X logical error)");
                logical_error_found = true;
            }
            
            // Check virtual boundary contact
            if (cluster_touches_virtual_boundary(root)) {
                log.debug("LOGICAL ERROR: Cluster touches virtual boundary");
                logical_error_found = true;
            }
            
            if (!logical_error_found) {
                log.debug(" NO LOGICAL ERROR: Cluster is correctable");
            }
        }
        
        log.debug("\n=== FINAL RESULT ===");
        log.info("Logical error detected: ", (logical_error_found ? "YES" : "NO"));
        log.debug("====================");
        
        return logical_error_found;
    }

    void UFDecoder::apply_edge_correction(uint edge_id, std::vector<uint8_t>& correction) {
        if (edge_id >= edge_list.size()) return;
        
        uint v1 = edge_list[edge_id].first;
        uint v2 = edge_list[edge_id].second;
        
        auto it = edge_observables_map.find({v1, v2});
        if (it != edge_observables_map.end()) {
            for (uint obs_id : it->second) {
                if (obs_id < correction.size()) {
                    correction[obs_id] ^= 1;
                }
            }
        }
    }

    void UFDecoder::debug_print_cluster_boundaries(Cluster* cluster) {
        log.debug("=== Cluster ", cluster->id, " Boundary Analysis ===");
        log.debug("Vertices: [", cluster->vertices.size(), " total]");
        
        bool has_left = false, has_right = false, has_virtual = false;
        std::vector<uint> left_vertices, right_vertices, virtual_vertices, regular_vertices;
        
        for (uint vertex_id : cluster->vertices) {
            auto it = vertex_boundaries.find(vertex_id);
            if (it != vertex_boundaries.end()) {
                switch (it->second) {
                    case 0: 
                        virtual_vertices.push_back(vertex_id);
                        has_virtual = true;
                        break;
                    case 1: 
                        left_vertices.push_back(vertex_id);
                        has_left = true;
                        break;
                    case 2: 
                        right_vertices.push_back(vertex_id);
                        has_right = true;
                        break;
                }
            } else {
                regular_vertices.push_back(vertex_id);
            }
        }
        
        log.debug("  Left boundary vertices: ", left_vertices.size(), " -> ", vector_to_string(left_vertices));
        log.debug("  Right boundary vertices: ", right_vertices.size(), " -> ", vector_to_string(right_vertices));
        log.debug("  Virtual boundary vertices: ", virtual_vertices.size(), " -> ", vector_to_string(virtual_vertices));
        log.debug("  Regular vertices: ", regular_vertices.size(), " -> ", vector_to_string(regular_vertices));
        
        log.debug("  Boundary spans: Left=", has_left, ", Right=", has_right, ", Virtual=", has_virtual);
        log.debug("  Parity: ", cluster->parity, " (", (cluster->parity % 2 == 0 ? "EVEN" : "ODD"), ")");
        
        // Determine logical error status
        bool spans_boundaries = has_left && has_right;
        bool touches_virtual = has_virtual;
        bool should_be_logical = (cluster->parity % 2 == 1) && (spans_boundaries || touches_virtual);
        
        log.debug("  LOGICAL ERROR PREDICTION: ", (should_be_logical ? "YES" : "NO"));
        log.debug("    Reason: parity=", cluster->parity % 2, ", spans=", spans_boundaries, ", virtual=", touches_virtual);
    }

    void UFDecoder::debug_print_surface_code_layout() {
        log.debug("=== Surface Code Layout Analysis ===");
        auto vertices = decoding_graph.vertices();
        
        // Find physical coordinate ranges first
        double min_x = 1e6, max_x = -1e6, min_y = 1e6, max_y = -1e6;
        int physical_count = 0;
        
        for (uint i = 0; i < vertices.size(); i++) {
            if (vertices[i]->coord.size() >= 2) {
                double x = vertices[i]->coord[0];
                double y = vertices[i]->coord[1];
                
                // Skip virtual boundary coordinates (very large values)
                if (x > 1e6 || y > 1e6) {
                    continue;
                }
                
                min_x = std::min(min_x, x); max_x = std::max(max_x, x);
                min_y = std::min(min_y, y); max_y = std::max(max_y, y);
                physical_count++;
            }
        }
        
        log.debug("Physical coordinate ranges: X=[", min_x, ", ", max_x, "], Y=[", min_y, ", ", max_y, "]");
        log.debug("Total vertices: ", vertices.size(), " (", physical_count, " physical, ", 
                  vertices.size() - physical_count, " virtual)");
        log.debug("Boundary vertices: ", vertex_boundaries.size());
        
        // Print boundary summary with virtual positioning info
        int virtual_count = 0, left_count = 0, right_count = 0;
        for (auto& [vertex_id, boundary_type] : vertex_boundaries) {
            switch (boundary_type) {
                case 0: virtual_count++; break;
                case 1: left_count++; break;
                case 2: right_count++; break;
            }
        }
        
        log.debug("Boundary distribution: Virtual=", virtual_count, ", Left=", left_count, ", Right=", right_count);
    }

    void UFDecoder::debug_print_ascii_grid() {
        log.debug("=== ASCII Grid Visualization ===");
        auto vertices = decoding_graph.vertices();
        
        std::map<std::pair<int, int>, char> grid;
        std::map<std::pair<int, int>, uint> grid_vertex_id;
        
        // First pass: find physical coordinate bounds
        int min_x = 100, max_x = -100, min_y = 100, max_y = -100;
        
        for (uint i = 0; i < vertices.size(); i++) {
            if (vertices[i]->coord.size() >= 2) {
                double raw_x = vertices[i]->coord[0];
                double raw_y = vertices[i]->coord[1];
                
                // Skip virtual boundary coordinates in range calculation
                if (raw_x > 1e6 || raw_y > 1e6) {
                    continue;
                }
                
                int x = static_cast<int>(raw_x);
                int y = static_cast<int>(raw_y);
                
                min_x = std::min(min_x, x); max_x = std::max(max_x, x);
                min_y = std::min(min_y, y); max_y = std::max(max_y, y);
            }
        }
        
        // Calculate safe virtual boundary positions
        int virtual_left = min_x - 1;   // Place virtual boundaries 2 units outside
        int virtual_right = max_x + 1;
        int virtual_top = max_y + 1;
        int virtual_bottom = min_y - 1;
        
        // Second pass: place all vertices including virtual boundaries
        for (uint i = 0; i < vertices.size(); i++) {
            if (vertices[i]->coord.size() >= 2) {
                double raw_x = vertices[i]->coord[0];
                double raw_y = vertices[i]->coord[1];
                
                int x, y;
                bool is_virtual = false;
                
                // Handle virtual boundary coordinates
                if (raw_x > 1e6 || raw_y > 1e6) {
                    is_virtual = true;
                    
                    // Position virtual boundaries based on their connections or boundary type
                    auto boundary_it = vertex_boundaries.find(i);
                    if (boundary_it != vertex_boundaries.end() && boundary_it->second == 0) {
                        // This is a virtual boundary vertex
                        // Try to position it based on its connections to physical vertices
                        auto adj_list = decoding_graph.adjacency_list(vertices[i]);
                        
                        if (!adj_list.empty()) {
                            // Find the first physical neighbor to determine positioning
                            for (auto adj : adj_list) {
                                uint adj_idx = std::find(vertices.begin(), vertices.end(), adj) - vertices.begin();
                                if (adj_idx < vertices.size() && adj->coord.size() >= 2) {
                                    double adj_x = adj->coord[0];
                                    double adj_y = adj->coord[1];
                                    
                                    if (adj_x < 1e6 && adj_y < 1e6) {
                                        // Position virtual boundary relative to physical neighbor
                                        int phys_x = static_cast<int>(adj_x);
                                        int phys_y = static_cast<int>(adj_y);
                                        
                                        // Place virtual boundary outside the grid
                                        if (phys_x == min_x) {
                                            x = virtual_left; y = phys_y;
                                        } else if (phys_x == max_x) {
                                            x = virtual_right; y = phys_y;
                                        } else if (phys_y == min_y) {
                                            x = phys_x; y = virtual_bottom;
                                        } else if (phys_y == max_y) {
                                            x = phys_x; y = virtual_top;
                                        } else {
                                            // Default positioning
                                            x = virtual_left; y = phys_y;
                                        }
                                        break;
                                    }
                                }
                            }
                        } else {
                            // No connections found, use default position
                            x = virtual_left;
                            y = (min_y + max_y) / 2;
                        }
                    } else {
                        // Skip other types of large coordinates
                        continue;
                    }
                } else {
                    // Physical coordinates
                    x = static_cast<int>(raw_x);
                    y = static_cast<int>(raw_y);
                }
                
                // Determine symbol
                char symbol = '.';
                if (i < vertex_has_syndrome.size() && vertex_has_syndrome[i]) {
                    symbol = 'S';  // Syndrome
                } else if (vertex_boundaries.find(i) != vertex_boundaries.end()) {
                    switch (vertex_boundaries[i]) {
                        case 0: symbol = 'V'; break;  // Virtual
                        case 1: symbol = 'L'; break;  // Left
                        case 2: symbol = 'R'; break;  // Right
                    }
                }
                
                // Add virtual marker for repositioned virtual boundaries
                if (is_virtual && symbol != 'V') {
                    symbol = 'v';  // lowercase for repositioned virtual
                }
                
                grid[{x, y}] = symbol;
                grid_vertex_id[{x, y}] = i;
            }
        }
        
        // Update grid bounds to include virtual boundaries
        int grid_min_x = std::min(min_x, virtual_left);
        int grid_max_x = std::max(max_x, virtual_right);
        int grid_min_y = std::min(min_y, virtual_bottom);
        int grid_max_y = std::max(max_y, virtual_top);
        
        // Print grid
        log.debug("Legend: S=Syndrome, V=Virtual, L=Left, R=Right, .=Regular, v=repositioned virtual");
        log.debug("Physical bounds: X=[", min_x, ", ", max_x, "], Y=[", min_y, ", ", max_y, "]");
        log.debug("Display bounds: X=[", grid_min_x, ", ", grid_max_x, "], Y=[", grid_min_y, ", ", grid_max_y, "]");
        
        for (int y = grid_max_y; y >= grid_min_y; y--) {
            std::string row = "Y=" + std::to_string(y) + ": ";
            for (int x = grid_min_x; x <= grid_max_x; x++) {
                auto it = grid.find({x, y});
                if (it != grid.end()) {
                    row += it->second;
                } else {
                    row += ' ';
                }
            }
            log.debug(row);
        }
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

} // namespace qrc