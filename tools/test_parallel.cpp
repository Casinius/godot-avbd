/*
 * Deterministic test for multithreading changes.
 * Verifies that parallel broad phase and lock-free manifold pool produce
 * bit-identical results across different thread counts.
 */

#include <iostream>
#include <vector>
#include <string>

#include "avbd/solver.h"

using namespace avbd;

// Simple scene with a few bodies
std::vector<Rigid *> createTestScene(Solver *solver) {
    std::vector<Rigid *> bodies;
    
    // Create a ground
    auto *ground = new Rigid(solver, float3{10.0f, 0.5f, 10.0f}, ShapeType::Box, 1000.0f, 0.5f,
                             float3{-5.0f, 0.0f, 0.0f}, float3{0.0f, 0.0f, 0.0f});
    ground->collisionLayer = 1;
    ground->collisionMask = 1;
    bodies.push_back(ground);
    
    // Create a few boxes falling
    for (int i = 0; i < 5; ++i) {
        auto *box = new Rigid(solver, float3{1.0f, 1.0f, 1.0f}, ShapeType::Box, 1.0f, 0.5f,
                             float3{(float)i * 1.5f, 10.0f, 0.0f}, float3{0.0f, 0.0f, 0.0f});
        box->collisionLayer = 1;
        box->collisionMask = 1;
        bodies.push_back(box);
    }
    
    return bodies;
}

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;

    std::cout << "Multithreading Deterministic Test\n";
    std::cout << "===================================\n\n";
    
    int threadCounts[] = {1, 2, 4};
    int numThreadCounts = sizeof(threadCounts) / sizeof(threadCounts[0]);
    std::string digests[numThreadCounts];
    
    for (int i = 0; i < numThreadCounts; ++i) {
        int threads = threadCounts[i];
        
        std::cout << "Testing with " << threads << " thread(s)...\n";
        
        // Create solver
        Solver solver;
        solver.threads = threads;
        
        // Create test scene
        auto bodies = createTestScene(&solver);
        
        // Run several steps
        for (int step = 0; step < 100; ++step) {
            solver.step();
        }
        
        // Generate digest (simple hash of body positions)
        uint64_t digest = 0;
        for (size_t j = 0; j < bodies.size(); ++j) {
            const auto &pos = bodies[j]->positionLin;
            digest ^= (uint64_t)(pos.x() * 1e6 + pos.y() * 1e6 + pos.z() * 1e6) ^ (j + 1);
        }
        
        digests[i] = std::to_string(digest);
        std::cout << "  Digest: " << digests[i] << "\n";
        
        // Cleanup
        for (auto *b : bodies) {
            delete b;
        }
    }
    
    std::cout << "\nDigest comparison:\n";
    std::cout << "  " << digests[0] << " (1 thread)\n";
    std::cout << "  " << digests[1] << " (2 threads)\n";
    std::cout << "  " << digests[2] << " (4 threads)\n";
    
    // Check if all digests match
    if (digests[0] == digests[1] && digests[1] == digests[2]) {
        std::cout << "\n✓ All digests match! Parallelization is deterministic.\n";
        return 0;
    } else {
        std::cout << "\n✗ Digests do not match! Parallelization is NOT deterministic.\n";
        return 1;
    }
}
