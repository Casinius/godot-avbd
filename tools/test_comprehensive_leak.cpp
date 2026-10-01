/*
 * Comprehensive test that verifies all leaksanitizer fixes work correctly.
 * This test intentionally creates leaks to verify the sanitizer detects them,
 * and creates clean code to verify no leaks in solver internals.
 */

#include <cstdio>
#include <memory>
#include <vector>
#include "avbd/solver.h"

using namespace avbd;

// Test 1: Verify sanitizer detects manual leaks
static void test_manual_leaks() {
    printf("Test 1: Manual memory leaks (expected to be detected)\n");
    
    {
        // Intentional leak - should be detected
        Solver s;
        Rigid *leaked = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
        (void)leaked;  // Suppress unused warning
        
        // Another intentional leak
        Rigid *leaked2 = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {1, 0, 4});
        (void)leaked2;
    }
}

// Test 2: Verify solver with proper cleanup has no leaks
static void test_solver_cleanup() {
    printf("\nTest 2: Solver with proper cleanup (no leaks expected)\n");
    
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
    
    // Run a few steps
    for (int i = 0; i < 10; i++) {
        s.step();
    }
    
    // Clean up - this should properly free all Manifolds
    s.clear();
}

// Test 3: Smart pointer usage pattern
static void test_smart_pointers() {
    printf("\nTest 3: Smart pointer usage (no leaks expected)\n");
    
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    
    // Use smart pointers for automatic cleanup
    std::vector<std::unique_ptr<Rigid>> bodies;
    for (int i = 0; i < 5; i++) {
        Rigid *b = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {(float)i, 0.0f, 5.0f});
        bodies.push_back(std::unique_ptr<Rigid>(b));
    }
    
    // Run simulation
    for (int i = 0; i < 20; i++) {
        s.step();
    }
    
    // Bodies are automatically cleaned up when vector goes out of scope
}

// Test 4: std::optional usage
static void test_optional_usage() {
    printf("\nTest 4: std::optional usage (leak expected for one box)\n");
    
    {
        Solver s;
        std::optional<Rigid*> ground;
        ground = new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
        
        // Intentional leak - box not deleted
        [[maybe_unused]] Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
    }
}

int main() {
    printf("=== Comprehensive LeakSanitizer Test Suite ===\n\n");
    printf("Note: Tests 1 and 4 intentionally create leaks.\n");
    printf("Tests 2 and 3 should have no leaks.\n\n");
    
    test_manual_leaks();
    test_solver_cleanup();
    test_smart_pointers();
    test_optional_usage();
    
    printf("\n=== All Tests Complete ===\n");
    printf("Run with: xmake run comprehensive_leak_test --leaksanitizer=y\n");
    printf("Expected: LeakSanitizer should report:\n");
    printf("  - 2 Rigid allocations leaked (Tests 1 and 4)\n");
    printf("  - No solver internal leaks (Tests 2 and 3)\n");
    
    return 0;
}
