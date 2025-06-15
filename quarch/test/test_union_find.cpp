#include "union_find_decoder.h"
#include "logger.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <limits>

namespace qrc {

class UFLogicalErrorTester {
private:
    UFDecoder* decoder;
    uint distance;
    stim::Circuit circuit;
    
public:
    UFLogicalErrorTester(uint d) : distance(d) {
        stim::CircuitGenParameters params(1, distance, "rotated_memory_z");
        params.before_round_data_depolarization = 0.001;
        params.after_clifford_depolarization = 0.0; 
        params.after_reset_flip_probability = 0.001;
        params.before_measure_flip_probability = 0.001;

        std::cout << "DEBUG: Creating circuit with distance=" << distance << std::endl;
        
        circuit = stim::generate_surface_code_circuit(params).circuit;
        
        decoder = new UFDecoder(circuit);
        std::cout << "Created rotated memory Z surface code:" << std::endl;
        std::cout << "- Distance: " << distance << std::endl;
        std::cout << "- Detectors: " << circuit.count_detectors() << std::endl;
        std::cout << "- Observables: " << circuit.count_observables() << std::endl;
    }
    
    ~UFLogicalErrorTester() {
        delete decoder;
    }
    
    // Add this new test function to the UFLogicalErrorTester class
    void test_simulated_syndromes() {
        std::cout << "\n=== Test: Simulated Real-World Syndrome Patterns ===" << std::endl;
        
        // Define the additional syndrome patterns from logs
        std::vector<std::vector<uint8_t>> syndromes = {
            {0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1},
            {0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0},
            {1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,1},
            {0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0},
            {0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0},
            {0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,1},
            {0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,1},
            {0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0},
            {1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0},
            {1,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1},
            {1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1},
            {1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0},
            {0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1},
            {1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0},
            {0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0},
            {0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0},
            {0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
            {1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,1,0},
            {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1}
        };
        
        // Run each test pattern
        int logical_errors = 0;
        
        for (size_t i = 0; i < syndromes.size(); i++) {
            auto& syndrome = syndromes[i];
            
            // Ensure syndrome vector has right size
            if (syndrome.size() < circuit.count_detectors()) {
                syndrome.resize(circuit.count_detectors(), 0);
            }
            
            std::cout << "\n--- Pattern "<< i << "  ---" << std::endl;
            
            // Count syndrome weight
            int syndrome_weight = 0;
            for (auto bit : syndrome) syndrome_weight += bit;
            
            auto result =  decoder->decode_error(syndrome);
            if (result.is_logical_error) {
                logical_errors++;
            }
        }
        
        std::cout<<"\n Logical errors detected in simulated patterns: " << logical_errors << std::endl;
    }
    
    void run_logical_error_tests() {
        std::cout << "Running Union-Find Logical Error Tests" << std::endl;
        std::cout << "Rotated Memory Z Surface Code (Distance " << distance << ")" << std::endl;
        std::cout << "==============================================" << std::endl;
        test_simulated_syndromes();
    }
};

} // namespace qrc

int main() {
    try {
        qrc::log.setLevel(qrc::LogLevel::DEBUG);
        std::cout << "Union-Find Decoder: Rotated Memory Z Surface Code Tests" << std::endl;
        std::cout << "=======================================================" << std::endl;
        
        // Test different distances
        for (uint distance : {5}) {
            qrc::UFLogicalErrorTester tester(distance);
            tester.run_logical_error_tests();
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}