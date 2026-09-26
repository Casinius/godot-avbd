/*
 * Simple deterministic test for multithreading changes.
 * Verifies that parallel broad phase and lock-free manifold pool compile and run without errors.
 */

#include <iostream>

#include "avbd/solver.h"

using namespace avbd;

int main(int argc, char **argv) {
    (void)argc;
    (void)argv;

    std::cout << "Simple Multithreading Test\n";
    std::cout << "===========================\n\n";

    int threadCounts[] = {1, 2, 4};
    int numThreadCounts = sizeof(threadCounts) / sizeof(threadCounts[0]);
    int failures = 0;

    for (int i = 0; i < numThreadCounts; ++i) {
        int threads = threadCounts[i];
        
        std::cout << "Testing with " << threads << " thread(s)...\n";
        
        // Create solver
        Solver solver;
        solver.threads = threads;
        
        // Create test scene (minimal - just a box)
        auto *box = new Rigid(&solver, float3{1.0f, 1.0f, 1.0f}, ShapeType::Box, 1.0f, 0.5f,
                             float3{0.0f, 5.0f, 0.0f}, float3{0.0f, 0.0f, 0.0f});
        box->collisionLayer = 1;
        box->collisionMask = 1;
        
        // Run 10 steps
        for (int step = 0; step < 10; ++step) {
            solver.step();
        }
        
        // Check that the box is still on screen (not fallen off the world)
        if (box->positionLin.y() < -10.0f) {
            std::cout << "  ✗ Box fell off world (y=" << box->positionLin.y() << ")\n";
            failures++;
        } else {
            std::cout << "  ✓ Box is at y=" << box->positionLin.y() << "\n";
        }
        
        delete box;
    }
    
    std::cout << "\nResults: " << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures > 0 ? 1 : 0;
}
