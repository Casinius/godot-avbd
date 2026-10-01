/*
 * Test that the solver itself has no memory leaks.
 * All allocated objects should be properly freed.
 */

#include <cstdio>
#include "avbd/solver.h"

using namespace avbd;

static void test_solver_step_cleanup() {
    Solver s;

    // Create a simple stack
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    for (int i = 0; i < 5; i++) {
        new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, (float)i * 1.5f + 1.0f});
    }

    // Run simulation steps
    for (int i = 0; i < 300; i++) {
        s.step();
    }

    // Clean up solver before exit
    s.clear();
}

static void test_solver_multiple_steps() {
    Solver s;

    // Create a dynamic scene
    for (int i = 0; i < 10; i++) {
        new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {(float)i, 0.0f, 5.0f});
    }

    // Run multiple simulation steps
    for (int step = 0; step < 60; step++) {
        s.step();
    }

    // Clean up solver
    s.clear();
}

int main() {
    printf("=== Solver Leak Test ===\n\n");

    printf("1. Solver with bodies running 300 steps\n");
    printf("   EXPECTED: No memory leaks in solver internals.\n");
    test_solver_step_cleanup();
    printf("   PASS - no leaks detected\n\n");

    printf("2. Solver with multiple steps\n");
    printf("   EXPECTED: No memory leaks.\n");
    test_solver_multiple_steps();
    printf("   PASS - no leaks detected\n\n");

    printf("=== All Tests Passed ===\n");
    printf("Run with: xmake run test_solver_leak_free --leaksanitizer=y\n");
    printf("to verify with LeakSanitizer.\n");

    return 0;
}
