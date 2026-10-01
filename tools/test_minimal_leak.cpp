/*
 * Simple test to verify LeakSanitizer works and solver is leak-free.
 */

#include <cstdio>
#include "avbd/solver.h"

using namespace avbd;

int main() {
    printf("=== Simple Leak Test ===\n");

    // Create a simple solver
    Solver s;

    // Add a ground and one box
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});

    // Run one step
    s.step();

    // Clean up
    s.clear();

    printf("Test completed successfully!\n");
    printf("Run with: xmake run test_simple_leak --leaksanitizer=y\n");

    return 0;
}
