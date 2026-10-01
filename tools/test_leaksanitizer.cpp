/*
 * Simple test to verify LeakSanitizer works with the AVBD solver.
 * Builds without leakanitizer flag: should have memory leaks.
 * Builds with --leaksanitizer=y: should show actual leaks.
 */

#include <cstdio>
#include <memory>
#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"

using namespace avbd;

// Test 1: Manually allocate and never delete (detectable by LeakSanitizer)
static void test_manual_leak() {
    Solver s;

    // Allocate without using smart pointers
    Rigid *leakedRigid = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    // Never delete leakedRigid - LEAK!
    (void)leakedRigid;
}

// Test 2: Smart pointer usage (should be leak-free)
static void test_smart_pointer_usage() {
    Solver s;

    // Use std::unique_ptr for automatic memory management
    // Rigid has 2 constructors; we use the one with explicit shape
    Rigid *body1 = new Rigid(&s, {1, 1, 1}, avbd::ShapeType::Box, 1.0f, 0.5f, {0, 0, 4});
    Rigid *body2 = new Rigid(&s, {1, 1, 1}, avbd::ShapeType::Box, 0.0f, 0.5f, {0, 0, 10});

    // Auto-delete using unique_ptr
    auto auto_body1 = std::unique_ptr<Rigid>(body1);
    auto auto_body2 = std::unique_ptr<Rigid>(body2);

    // No manual delete needed - smart pointers handle it
    // Simulate some steps
    for (int i = 0; i < 240; i++) {
        s.step();
    }
}

// Test 3: Vector of smart pointers (demonstrates modern C++ ownership)
static void test_vector_of_smart_pointers() {
    Solver s;

    std::vector<std::unique_ptr<Rigid>> bodies;

    // Add bodies to a vector
    for (int i = 0; i < 10; i++) {
        Rigid *b = new Rigid(&s, {1, 1, 1}, avbd::ShapeType::Box, 1.0f, 0.5f, {0, 0, (float)i * 1.5f + 1.0f});
        bodies.push_back(std::unique_ptr<Rigid>(b));
    }

    // Vector manages all bodies automatically
    for (int i = 0; i < 300; i++) {
        s.step();
    }
    // Bodies are automatically deleted when vector goes out of scope
}

// Test 4: Optional body reference (demonstrates use of <optional>)
static void test_optional_body() {
    Solver s;

    // Using optional for bodies that might not exist
    std::optional<Rigid*> ground;
    ground = new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    [[maybe_unused]] Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});

    // Use the optional body
    if (ground) {
        // Do something with ground
        (void)ground;
    }

    for (int i = 0; i < 240; i++) {
        s.step();
    }

    // Note: box is leaked in this test
    // In production, you would delete box or use smart pointers
}

int main() {
    printf("=== LeakSanitizer Test Suite ===\n\n");

    printf("1. Manual allocation leak test\n");
    printf("   Allocates memory without deletion.\n");
    printf("   EXPECTED: LeakSanitizer reports 1 allocation leaked.\n");
    test_manual_leak();
    printf("\n");

    printf("2. Smart pointer test\n");
    printf("   Uses std::unique_ptr for automatic memory management.\n");
    printf("   EXPECTED: No leaks detected.\n");
    test_smart_pointer_usage();
    printf("\n");

    printf("3. Vector of smart pointers test\n");
    printf("   Demonstrates modern C++ container ownership.\n");
    printf("   EXPECTED: No leaks detected.\n");
    test_vector_of_smart_pointers();
    printf("\n");

    printf("4. Optional body test\n");
    printf("   Uses std::optional for optional references.\n");
    printf("   EXPECTED: LeakSanitizer reports 1 allocation leaked (the box).\n");
    test_optional_body();
    printf("\n");

    printf("=== Test Suite Complete ===\n");
    printf("Run with: xmake run test_leaksanitizer --leaksanitizer=y\n");
    printf("to see detailed leak reports from LeakSanitizer.\n");

    return 0;
}
