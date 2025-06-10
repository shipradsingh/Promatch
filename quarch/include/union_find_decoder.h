#ifndef UNION_FIND_DECODER_H
#define UNION_FIND_DECODER_H

#include "decoder.h"
#include "decoding_graph.h"
#include <vector>
#include <map>
#include <set>
#include <limits>

namespace qrc {

    // Simple edge states for Union-Find algorithm
    enum class EdgeSupport {
        NONE = 0,           // Not processed
        GROWN = 1,          // Part of spanning forest
        MATCHING = 2        // Selected for correction
    };

    // Simplified cluster structure matching documentation
    struct Cluster {
        uint id;
        uint instance;
        uint parity;               // Number of defects in cluster (odd/even)
        std::vector<uint> vertices;
        std::vector<uint> tree_edges;  // Edges that form spanning tree
        bool is_frozen;            // Frozen when even parity OR touches boundary
        
        Cluster* parent;           // Union-Find parent pointer
        uint rank;                 // Union-Find rank
        
        Cluster(uint cluster_id, uint inst) 
            : id(cluster_id), instance(inst), parity(0), is_frozen(false), 
            parent(this), rank(0) {}
        
        void add_vertex(uint vertex_id) {
            vertices.push_back(vertex_id);
        }
        
        void add_tree_edge(uint edge_id) {
            tree_edges.push_back(edge_id);
        }
        
        Cluster* find() {
            if (parent != this) {
                parent = parent->find(); // Path compression
            }
            return parent;
        }
        
        void union_with(Cluster* other) {
            Cluster* root1 = this->find();
            Cluster* root2 = other->find();
            
            if (root1 == root2) return;
            
            // Preserve freezing state
            bool result_frozen = root1->is_frozen || root2->is_frozen;
            
            // Union by rank
            if (root1->rank < root2->rank) {
                root1->parent = root2;
                for (uint v : root1->vertices) {
                    root2->vertices.push_back(v);
                }
                root2->is_frozen = result_frozen;
            } else if (root1->rank > root2->rank) {
                root2->parent = root1;
                for (uint v : root2->vertices) {
                    root1->vertices.push_back(v);
                }
                root1->is_frozen = result_frozen;
            } else {
                root2->parent = root1;
                root1->rank++;
                for (uint v : root2->vertices) {
                    root1->vertices.push_back(v);
                }
                root1->is_frozen = result_frozen;
            }
        }
    };

    class UFDecoder : public Decoder {
    private:
        // Core algorithm data
        uint cluster_id;
        uint current_instance;
        uint n_detectors;
        DecodingGraph decoding_graph;
        
        // Clusters (as per documentation)
        std::vector<Cluster*> clusters;
        std::vector<Cluster*> vertex_to_cluster;
        std::vector<bool> vertex_has_syndrome;
        
        // Edge management (simplified)
        std::vector<EdgeSupport> edge_states;
        std::vector<std::pair<uint, uint>> edge_list;
        std::vector<double> edge_weights;
        std::vector<uint> sorted_edges;
        
        // Boundary detection
        std::map<uint, uint> vertex_boundaries;  // vertex_id -> boundary_type
        std::map<std::pair<uint, uint>, std::vector<uint>> edge_observables_map;
        std::vector<uint> detector_to_vertex_id;
        
        // Setup methods
        void setup_boundaries_from_decoding_graph();
        void build_edge_mapping();
        
        // Algorithm phases (matching documentation)
        void extract_vertex_syndrome(const std::vector<uint8_t>& detector_syndrome);
        void initialize_clusters();    // Phase 1: Create singleton clusters
        void grow_and_merge_clusters(); // Phase 2: Build spanning forest  
        void peel_clusters();          // Phase 3: Extract corrections
        
        // Freezing logic (as documented)
        bool cluster_touches_virtual_boundary(Cluster* cluster);
        bool should_freeze_cluster(Cluster* cluster);
        void freeze_completed_clusters();
        bool is_cluster_frozen(Cluster* cluster);
        
        // Helper methods
        void cluster_add_vertex(Cluster* cluster, uint vertex_id);
        bool check_logical_error();
        void apply_edge_correction(uint edge_id, std::vector<uint8_t>& correction);
        
        // Debug visualization functions
        void debug_print_cluster_boundaries(Cluster* cluster);
        void debug_print_surface_code_layout();
        void debug_print_ascii_grid();
        std::string vector_to_string(const std::vector<uint>& vec);
        
    public:
        UFDecoder(const stim::Circuit& circuit);
        ~UFDecoder() {
            // Clean up clusters
            for (auto cluster : clusters) {
                delete cluster;
            }
        }
        
        DecoderShotResult decode_error(const std::vector<uint8_t>& syndrome) override;
        
        std::string name() override { return "UnionFind"; }
        bool is_software() override { return true; }
        uint64_t sram_cost() override { 
            return clusters.size() * sizeof(Cluster) + 
                   edge_list.size() * sizeof(std::pair<uint, uint>);
        }
    };

} // namespace qrc

#endif