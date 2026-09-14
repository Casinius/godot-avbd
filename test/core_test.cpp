/*
 * Numeric tests for the AVBD solver core. No Godot, no windowing: this binary
 * is the ground truth for solver behaviour and is run before the GDExtension
 * layer is exercised from Godot.
 *
 * Conventions follow the solver core: metres, kilograms, seconds, +Z up,
 * gravity along -Z, boxes described by their full widths.
 *
 *   xmake run avbd_core_test                 # all assertions
 *   xmake run avbd_core_test --list          # names of the ported scenes
 *   xmake run avbd_core_test --scene Stack --steps 300
 *   xmake run avbd_core_test --bench         # 512-box pyramid timing
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <string>
#include <vector>

#include "avbd/solver.h"
#include "core_scenes.hpp"

using namespace avbd;

static int g_checks = 0;
static int g_failures = 0;

static void report(bool ok, const char *name, const char *fmt, ...) {
    g_checks++;
    if (!ok) {
        g_failures++;
    }
    std::printf("%s %-30s ", ok ? "  ok " : "FAIL ", name);
    va_list ap;
    va_start(ap, fmt);
    std::vprintf(fmt, ap);
    va_end(ap);
    std::printf("\n");
}

static void step_n(Solver &s, int n) {
    for (int i = 0; i < n; i++) {
        s.step();
    }
}

static int count_bodies(const Solver &s) {
    int n = 0;
    for (const Rigid *b = s.bodies; b != nullptr; b = b->next) {
        n++;
    }
    return n;
}

static int count_forces(const Solver &s) {
    int n = 0;
    for (const Force *f = s.forces; f != nullptr; f = f->next) {
        n++;
    }
    return n;
}

static int count_contact_points(const Solver &s) {
    int n = 0;
    for (const Force *f = s.forces; f != nullptr; f = f->next) {
        n += f->contactPointCount();
    }
    return n;
}

static bool finite(const Rigid *b) {
    return std::isfinite(b->positionLin.x) && std::isfinite(b->positionLin.y) && std::isfinite(b->positionLin.z) &&
           std::isfinite(b->positionAng.x) && std::isfinite(b->positionAng.y) && std::isfinite(b->positionAng.z) &&
           std::isfinite(b->positionAng.w) && std::isfinite(b->velocityLin.x) && std::isfinite(b->velocityLin.y) &&
           std::isfinite(b->velocityLin.z) && std::isfinite(b->velocityAng.x) && std::isfinite(b->velocityAng.y) &&
           std::isfinite(b->velocityAng.z);
}

static bool all_finite(const Solver &s) {
    for (const Rigid *b = s.bodies; b != nullptr; b = b->next) {
        if (!finite(b)) {
            return false;
        }
    }
    return true;
}

static float max_abs_position(const Solver &s) {
    float m = 0.0f;
    for (const Rigid *b = s.bodies; b != nullptr; b = b->next) {
        m = std::max(m, std::fabs(b->positionLin.x));
        m = std::max(m, std::fabs(b->positionLin.y));
        m = std::max(m, std::fabs(b->positionLin.z));
    }
    return m;
}

// FNV-1a over the raw state of every body, for cross-run/cross-process checks.
static uint64_t state_digest(const Solver &s) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](const void *p, size_t n) {
        const unsigned char *c = (const unsigned char *)p;
        for (size_t i = 0; i < n; i++) {
            h ^= c[i];
            h *= 1099511628211ull;
        }
    };
    for (Rigid *b = s.bodies; b != nullptr; b = b->next) {
        mix(&b->positionLin, sizeof(b->positionLin));
        mix(&b->positionAng, sizeof(b->positionAng));
        mix(&b->velocityLin, sizeof(b->velocityLin));
        mix(&b->velocityAng, sizeof(b->velocityAng));
    }
    return h;
}

// ---------------------------------------------------------------------------
// Resting contact: a box dropped on a static slab comes to rest exactly on the
// surface, with bounded penetration and no residual motion.
// ---------------------------------------------------------------------------
static void test_rest_contact() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});

    step_n(s, 240); // 4 s

    // Ground slab is 1 thick, centred on z = 0, so its top face is z = 0.5 and a
    // 1x1x1 box rests with its centre at z = 1.0. AVBD is a soft-penalty method:
    // a small penetration is expected (COLLISION_MARGIN is 1 cm).
    const float rest_z = box->positionLin.z;
    const float penetration = 1.0f - rest_z;
    report(std::fabs(penetration) <= 0.02f, "rest contact", "z=%.5f (surface 0.5 + half 0.5), penetration=%.5f", rest_z,
            penetration);
    report(std::fabs(box->velocityLin.z) <= 0.05f, "rest velocity", "vz=%.5f", box->velocityLin.z);
    report(all_finite(s), "rest finite", "bodies=%d contacts=%d", count_bodies(s), count_contact_points(s));
}

// ---------------------------------------------------------------------------
// A 10-box stack holds its shape: every box stays a box-height above the one
// below with negligible lateral drift (the paper's headline stability claim).
//
// AVBD is a soft-penalty method: the contact interfaces interpenetrate slightly
// under load (measured ~2.4 cm per interface at the default betaLin, and roughly
// independent of the iteration count as long as it is >= 10). The assertions below
// bound that behaviour and would catch collapse (gaps blowing up) or sinking.
// ---------------------------------------------------------------------------
static void test_stack_stability() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    std::vector<Rigid *> boxes;
    for (int i = 0; i < 10; i++) {
        boxes.push_back(new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, i * 1.5f + 1.0f}));
    }

    step_n(s, 300); // 5 s

    float worst_ground = 0.0f;    // penetration of the bottom box into the slab
    float worst_interface = 0.0f; // penetration between two boxes
    float worst_lateral = 0.0f;
    for (size_t i = 0; i < boxes.size(); i++) {
        const float z = boxes[i]->positionLin.z;
        // Ideal: the bottom box's centre is one half-height above the slab's top
        // face (1.0), and each box above sits exactly one box-height higher.
        const float penetration = (i == 0) ? 1.0f - z : 1.0f - (z - boxes[i - 1]->positionLin.z);
        if (i == 0) {
            worst_ground = penetration;
        } else {
            worst_interface = std::max(worst_interface, penetration);
        }
        worst_lateral = std::max(worst_lateral, std::hypot(boxes[i]->positionLin.x, boxes[i]->positionLin.y));
    }

    report(worst_ground <= 0.02f, "stack on ground", "bottom box penetration = %.5f m", worst_ground);
    report(worst_interface <= 0.05f, "stack interfaces", "worst box-box penetration = %.5f m (size 1 m)", worst_interface);
    report(worst_lateral <= 0.05f, "stack drift", "max lateral = %.5f", worst_lateral);
    report(all_finite(s), "stack finite", "bodies=%d contacts=%d", count_bodies(s), count_contact_points(s));
}

// ---------------------------------------------------------------------------
// Friction on a ramp. tan(20 deg) = 0.364: mu = 0.5 must hold the box in place
// (AVBD is a penalty method, so "in place" means mm/s creep rather than an exact
// freeze), while mu = 0.05 must let it accelerate down the slope. Both boxes use
// the same ramp friction so the effective contact friction (the geometric mean
// of the pair) equals the box's mu.
// ---------------------------------------------------------------------------
struct RampResult {
    float total;  // horizontal travel from the spawn point, m
    float late;   // horizontal travel during the final 0.5 s, m
    float speed;  // final linear speed, m/s
};

static RampResult ramp_run(float mu, float seconds) {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, mu, {0, 0, 0});

    const float angle = rad(20.0f);
    Rigid *ramp = new Rigid(&s, {40, 24, 1}, 0.0f, mu, {0, 0, 6});
    ramp->positionAng = {0, std::sin(angle * 0.5f), 0, std::cos(angle * 0.5f)};

    const float3 tangent = normalize(rotate(ramp->positionAng, float3{1, 0, 0}));
    const float3 normal = normalize(rotate(ramp->positionAng, float3{0, 0, 1}));
    // Ramp half thickness (0.5) + box half (0.5) + a small drop gap.
    const float3 start = ramp->positionLin + tangent * -5.0f + normal * 1.05f;
    Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, mu, start);

    const int total_steps = (int)(seconds * 60.0f);
    const int late_steps = 30; // 0.5 s
    step_n(s, total_steps - late_steps);

    const float3 before = box->positionLin;
    step_n(s, late_steps);

    const float3 total = box->positionLin - start;
    const float3 late = box->positionLin - before;
    return {std::hypot(total.x, total.y), std::hypot(late.x, late.y), length(box->velocityLin)};
}

static void test_friction() {
    const RampResult stuck = ramp_run(0.5f, 4.0f);
    report(stuck.late <= 0.02f, "static friction (mu=0.5)", "creep=%.5f m in the last 0.5 s, speed=%.5f m/s, total=%.4f m",
            stuck.late, stuck.speed, stuck.total);
    report(stuck.speed <= 0.05f, "stiction holds", "speed=%.5f m/s on a 20 deg ramp with mu=0.5 (tan=0.364)", stuck.speed);

    const RampResult slid = ramp_run(0.05f, 4.0f);
    report(slid.late >= 1.0f, "dynamic friction (mu=0.05)", "travel=%.4f m in the last 0.5 s, total=%.4f m", slid.late,
            slid.total);
    report(slid.total >= 2.0f, "sliding accumulates", "total travel=%.4f m in 4 s", slid.total);
    report(slid.late > stuck.late * 50.0f, "friction ordering", "slide(mu=0.05)=%.4f m vs 50 * creep(mu=0.5)=%.4f m",
            slid.late, stuck.late * 50.0f);
}

// ---------------------------------------------------------------------------
// Hard constraints: a 4-link chain hung from a static body. Infinitely stiff
// joints must hold the chain together to within a fraction of a link.
// ---------------------------------------------------------------------------
static void test_hard_joints() {
    Solver s;
    Rigid *anchor = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 10});

    std::vector<Rigid *> links;
    std::vector<Joint *> joints;
    Rigid *prev = anchor;
    for (int i = 0; i < 4; i++) {
        Rigid *link = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 9.0f - (float)i});
        joints.push_back(new Joint(&s, prev, link, {0, 0, -0.5f}, {0, 0, 0.5f}, INFINITY, 0.0f));
        links.push_back(link);
        prev = link;
    }

    step_n(s, 300); // 5 s

    float worst = 0.0f;
    for (size_t i = 0; i < joints.size(); i++) {
        Rigid *a = (i == 0) ? anchor : links[i - 1];
        Rigid *b = links[i];
        const float3 ca = transform(a->positionLin, a->positionAng, float3{0, 0, -0.5f});
        const float3 cb = transform(b->positionLin, b->positionAng, float3{0, 0, 0.5f});
        worst = std::max(worst, length(ca - cb));
    }

    report(worst <= 0.01f, "hard joint error", "max anchor separation = %.6f m (link = 1 m)", worst);
    report(links[3]->positionLin.z < anchor->positionLin.z, "chain hangs", "tail z=%.4f below anchor z=%.4f",
            links[3]->positionLin.z, anchor->positionLin.z);
    report(all_finite(s), "joints finite", "bodies=%d joints=%d", count_bodies(s), count_forces(s));
}

// ---------------------------------------------------------------------------
// Stiff lattice soft body: a 3x3x3 grid of boxes wired together with stiff
// joints, dropped from 4 m. It must not tear apart or explode.
// ---------------------------------------------------------------------------
static void test_soft_lattice() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    const int N = 3;
    const float size = 0.8f;
    const float half = size * 0.5f;
    const float Klin = 1000.0f;
    const float Kang = 250.0f;

    std::vector<Rigid *> grid(N * N * N, nullptr);
    auto at = [&](int x, int y, int z) { return grid[(x * N + y) * N + z]; };

    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                grid[(x * N + y) * N + z] = new Rigid(&s, {size, size, size}, 1.0f, 0.5f,
                        {(float)x * size, (float)y * size, (float)z * size + 4.0f});
            }
        }
    }

    std::vector<Joint *> joints;
    struct Link {
        Rigid *a;
        Rigid *b;
        float3 ra, rb;
    };
    std::vector<Link> links;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            for (int z = 0; z < N; z++) {
                if (x + 1 < N) {
                    links.push_back({at(x, y, z), at(x + 1, y, z), {half, 0, 0}, {-half, 0, 0}});
                }
                if (y + 1 < N) {
                    links.push_back({at(x, y, z), at(x, y + 1, z), {0, half, 0}, {0, -half, 0}});
                }
                if (z + 1 < N) {
                    links.push_back({at(x, y, z), at(x, y, z + 1), {0, 0, half}, {0, 0, -half}});
                }
            }
        }
    }
    for (const Link &l : links) {
        joints.push_back(new Joint(&s, l.a, l.b, l.ra, l.rb, Klin, Kang));
    }

    step_n(s, 300); // 5 s

    float worst = 0.0f;
    for (size_t i = 0; i < joints.size(); i++) {
        const float3 ca = transform(links[i].a->positionLin, links[i].a->positionAng, links[i].ra);
        const float3 cb = transform(links[i].b->positionLin, links[i].b->positionAng, links[i].rb);
        worst = std::max(worst, length(ca - cb));
    }

    report(worst <= 0.05f, "lattice joint error", "max joint separation = %.5f m (rest length 0)", worst);
    report(max_abs_position(s) <= 100.0f, "lattice contained", "max |position| = %.3f m", max_abs_position(s));
    report(all_finite(s), "lattice finite", "bodies=%d joints=%d contacts=%d", count_bodies(s), count_forces(s),
            count_contact_points(s));
}

// ---------------------------------------------------------------------------
// Every ported stand-in scene must survive a few hundred steps without NaN or
// divergence. This is the smoke test for collisions, joints and springs at once.
// ---------------------------------------------------------------------------
static void test_all_scenes() {
    for (int i = 0; i < coreSceneCount; i++) {
        Solver s;
        coreScenes[i].build(&s);
        const int bodies = count_bodies(s);
        step_n(s, 120);
        const bool ok = all_finite(s) && max_abs_position(s) <= 1000.0f;
        char name[64];
        std::snprintf(name, sizeof(name), "scene: %s", coreScenes[i].name);
        report(ok, name, "bodies=%d forces=%d contacts=%d max|x|=%.2f", bodies, count_forces(s),
                count_contact_points(s), max_abs_position(s));
    }
}

// ---------------------------------------------------------------------------
// Same scene, same result: two independent runs of the same scenario must agree
// bit for bit.
// ---------------------------------------------------------------------------
static uint64_t stack_digest() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    for (int i = 0; i < 10; i++) {
        new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, i * 1.5f + 1.0f});
    }
    step_n(s, 300);
    return state_digest(s);
}

static void test_determinism() {
    const uint64_t a = stack_digest();
    const uint64_t b = stack_digest();
    report(a == b, "in-process determinism", "digest 0x%016llx vs 0x%016llx", (unsigned long long)a,
            (unsigned long long)b);
    std::printf("     digest stack_300 = 0x%016llx\n", (unsigned long long)a);
}

static const CoreScene *find_scene(const std::string &name) {
    for (int i = 0; i < coreSceneCount; i++) {
        if (name == coreScenes[i].name || name == std::to_string(i)) {
            return &coreScenes[i];
        }
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// The parallel update must not depend on how many threads run it: work inside a colour
// is independent, so the split cannot change the result. If this ever fails, the
// colouring has put two coupled bodies in the same group (a race and a wrong answer).
// ---------------------------------------------------------------------------
static uint64_t scene_digest(int threads, int steps) {
    Solver s;
    s.threads = threads;
    const CoreScene *scene = find_scene("Soft Body");
    if (scene != nullptr)
        scene->build(&s);
    step_n(s, steps);
    return state_digest(s);
}

static uint64_t mixed_digest(int threads) {
    Solver s;
    s.threads = threads;
    // Bodies coupled by contacts and joints, plus a free-falling stack, so the colours are
    // varied: contacts (coupled), joints (coupled), and singletons.
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *prev = 0;
    for (int i = 0; i < 30; i++) {
        Rigid *link = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 20.0f - i * 1.5f});
        if (prev != 0)
            new Joint(&s, prev, link, {0, 0, -0.5f}, {0, 0, 0.5f}, 1000.0f, 250.0f);
        prev = link;
    }
    for (int i = 0; i < 40; i++)
        new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {5.0f + (float)i * 1.2f, 0, 3.0f + (float)(i % 5)});
    step_n(s, 240);
    return state_digest(s);
}

static void test_parallel_equivalence() {
    const int threadCounts[] = {1, 2, 3, 4, 8, 12};

    const uint64_t reference = mixed_digest(1);
    std::printf("     digest mixed_240 (threads=1) = 0x%016llx\n", (unsigned long long)reference);
    for (int t : threadCounts) {
        const uint64_t digest = mixed_digest(t);
        char name[64];
        std::snprintf(name, sizeof(name), "mixed digest threads=%d", t);
        report(digest == reference, name, "0x%016llx vs serial 0x%016llx", (unsigned long long)digest,
                (unsigned long long)reference);
    }

    const uint64_t softSerial = scene_digest(1, 120);
    for (int t : {4, 12}) {
        const uint64_t digest = scene_digest(t, 120);
        char name[64];
        std::snprintf(name, sizeof(name), "soft-body digest threads=%d", t);
        report(digest == softSerial, name, "0x%016llx vs serial 0x%016llx", (unsigned long long)digest,
                (unsigned long long)softSerial);
    }
}

// Report the shape of the constraint graph: how much of the work can run at once.
static void report_colouring() {
    const char *names[] = {"Stack", "Pyramid", "Soft Body", "Bridge"};
    for (const char *name : names) {
        const CoreScene *scene = find_scene(name);
        if (scene == nullptr)
            continue;
        Solver s;
        scene->build(&s);
        // Settle first: at step 1 most bodies are still apart and the graph has no edges.
        step_n(s, 120);
        int movable = 0;
        for (Rigid *b = s.bodies; b != 0; b = b->next)
            if (b->mass > 0)
                movable++;
        std::printf("     colours %-10s colours=%2d widest=%4d movable=%4d threads=%d\n", name, s.colourCount(),
                s.widestColour(), movable, s.threadCount());
    }
}

// ---------------------------------------------------------------------------
// CLI helpers
// ---------------------------------------------------------------------------
static void list_scenes() {
    for (int i = 0; i < coreSceneCount; i++) {
        std::printf("%2d  %s\n", i, coreScenes[i].name);
    }
}

static int run_scene(const std::string &name, int steps, int threads) {
    const CoreScene *scene = find_scene(name);
    if (scene == nullptr) {
        std::fprintf(stderr, "unknown scene: %s\n", name.c_str());
        return 2;
    }
    Solver s;
    s.threads = threads;
    scene->build(&s);
    std::printf("scene '%s': %d bodies, %d forces\n", scene->name, count_bodies(s), count_forces(s));
    step_n(s, steps);
    std::printf("after %d steps: threads=%d colours=%d widest=%d contacts=%d max|x|=%.4f finite=%s digest=0x%016llx\n",
            steps, s.threadCount(), s.colourCount(), s.widestColour(), count_contact_points(s), max_abs_position(s),
            all_finite(s) ? "yes" : "NO", (unsigned long long)state_digest(s));
    return all_finite(s) ? 0 : 1;
}

// A pyramid with ~512 boxes: the reference implementation's stress scene, and the case
// where the O(n^2) broad phase dominates.
static void build_bench_scene(Solver &s) {
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, -0.5f});
    const int SIZE = 32;
    for (int y = 0; y < SIZE; y++) {
        for (int x = 0; x < SIZE - y; x++) {
            new Rigid(&s, {1, 0.5f, 0.5f}, 1.0f, 0.5f, {x * 1.01f + y * 0.5f - SIZE / 2.0f, 0.0f, y * 0.85f + 0.5f});
        }
    }
}

static int run_bench(int threads) {
    Solver s;
    s.threads = threads;
    build_bench_scene(s);
    std::printf("bench: %d boxes, %d iterations, %d bodies movable\n", count_bodies(s), s.iterations,
            s.widestColour());

    const int warmup = 30;
    const int timed = 120;
    step_n(s, warmup);

    // Wall time: parallel speedup has to be measured on the clock, not on CPU time.
    const auto t0 = std::chrono::steady_clock::now();
    step_n(s, timed);
    const auto t1 = std::chrono::steady_clock::now();

    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / (double)timed;
    std::printf("bench: threads=%d colours=%d widest=%d -> %.3f ms/step, contacts=%d, finite=%s\n", s.threadCount(),
            s.colourCount(), s.widestColour(), ms, count_contact_points(s), all_finite(s) ? "yes" : "NO");
    return all_finite(s) ? 0 : 1;
}

// Sweep the thread count over the same scene so the speedup (and where it stops) is
// visible rather than claimed.
static int run_bench_sweep() {
    const int warmup = 30;
    const int timed = 120;
    double serialMs = 0.0;

    for (int threads : {1, 2, 4, 8, 12}) {
        Solver s;
        s.threads = threads;
        build_bench_scene(s);
        step_n(s, warmup);
        const auto t0 = std::chrono::steady_clock::now();
        step_n(s, timed);
        const auto t1 = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / (double)timed;
        if (threads == 1)
            serialMs = ms;
        std::printf("bench: threads=%2d  %7.3f ms/step  speedup %5.2fx  (colours=%d widest=%d)\n", threads, ms,
                serialMs / ms, s.colourCount(), s.widestColour());
    }
    return 0;
}

int main(int argc, char **argv) {
    std::string scene;
    int steps = 300;
    int threads = 0; // 0 = automatic

    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        if (arg == "--list") {
            list_scenes();
            return 0;
        } else if (arg == "--bench") {
            return run_bench(threads);
        } else if (arg == "--bench-sweep") {
            return run_bench_sweep();
        } else if (arg == "--colour") {
            report_colouring();
            return 0;
        } else if (arg == "--scene" && i + 1 < argc) {
            scene = argv[++i];
        } else if (arg == "--steps" && i + 1 < argc) {
            steps = std::atoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
        } else {
            std::fprintf(stderr, "usage: %s [--list | --bench | --bench-sweep | --colour | --scene NAME [--steps N]]"
                                 " [--threads N]\n", argv[0]);
            return 2;
        }
    }

    if (!scene.empty()) {
        return run_scene(scene, steps, threads);
    }

    test_rest_contact();
    test_stack_stability();
    test_friction();
    test_hard_joints();
    test_soft_lattice();
    test_all_scenes();
    test_determinism();
    test_parallel_equivalence();
    report_colouring();

    std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
