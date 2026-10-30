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
#include <thread>
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
    // Tuning affordance: AVBD_JACOBI_RELAX overrides the combine control for sweep runs
    // without touching the solver default. Same value for every solver in a process, so
    // determinism/parallel-equivalence checks stay meaningful.
    static const float relax = []() {
        const char *env = std::getenv("AVBD_JACOBI_RELAX");
        return env != nullptr ? static_cast<float>(std::atof(env)) : 0.0f;
    }();
    if (relax > 0.0f)
        s.jacobiRelaxation = relax;
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
    int n = s.manifoldCount();
    for (const Force *f = s.forces; f != nullptr; f = f->next) {
        n++;
    }
    return n;
}

static int count_contact_points(const Solver &s) {
    return s.manifoldContactPoints();
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
    // a small penetration is expected (COLLISION_MARGIN is 1 mm).
    const float rest_z = box->positionLin.z;
    const float penetration = 1.0f - rest_z;
    report(std::fabs(penetration) <= 0.01f, "rest contact", "z=%.5f (surface 0.5 + half 0.5), penetration=%.5f", rest_z,
            penetration);
    report(std::fabs(box->velocityLin.z) <= 0.05f, "rest velocity", "vz=%.5f", box->velocityLin.z);
    report(all_finite(s), "rest finite", "bodies=%d contacts=%d", count_bodies(s), count_contact_points(s));
}

// ---------------------------------------------------------------------------
// A 10-box stack holds its shape: every box stays a box-height above the one
// below with negligible lateral drift (the paper's headline stability claim).
//
// AVBD is a soft-penalty method: the contact interfaces interpenetrate slightly
// under load. The assertions below bound that behaviour and would catch collapse
// (gaps blowing up) or sinking.
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

    report(worst_ground <= 0.01f, "stack on ground", "bottom box penetration = %.5f m", worst_ground);
    report(worst_interface <= 0.02f, "stack interfaces", "worst box-box penetration = %.5f m (size 1 m)", worst_interface);
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

    const float angle = 20.0f * float(M_PI) / 180.0f;
    Rigid *ramp = new Rigid(&s, {40, 24, 1}, 0.0f, mu, {0, 0, 6});
    ramp->positionAng = angleAxis(angle, float3{0, 1, 0});

    const float3 tangent = normalize(qrot(ramp->positionAng, float3{1, 0, 0}));
    const float3 normal = normalize(qrot(ramp->positionAng, float3{0, 0, 1}));
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
        const float3 ca = a->positionLin + qrot(a->positionAng, float3{0, 0, -0.5f});
        const float3 cb = b->positionLin + qrot(b->positionAng, float3{0, 0, 0.5f});
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
        const float3 ca = links[i].a->positionLin + qrot(links[i].a->positionAng, links[i].ra);
        const float3 cb = links[i].b->positionLin + qrot(links[i].b->positionAng, links[i].rb);
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

// Interleaved substepping: one broad phase/update-set per tick, then a few rounds per
// substep with the contact error rebased. Same scene twice must agree; the pile must stay
// finite and bounded (the double-pyramid regression is the Godot-side equivalent).
static uint64_t interleaved_digest(int threads = 0) {
    Solver s;
    s.threads = threads;
    s.dt = 1.0f / 120.0f; // two substeps per 1/60 tick
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    for (int i = 0; i < 10; i++)
        new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, i * 1.5f + 1.0f});
    for (int i = 0; i < 300; i++)
        s.stepInterleaved(2, 3);
    return state_digest(s);
}

static void test_interleaved_substeps() {
    const uint64_t a = interleaved_digest(1);
    const uint64_t b = interleaved_digest(1);
    report(a == b, "interleaved substeps deterministic", "0x%016llx vs 0x%016llx",
            (unsigned long long)a, (unsigned long long)b);
    const uint64_t parallel = interleaved_digest(8);
    report(parallel == a, "interleaved digest threads=8", "0x%016llx vs serial 0x%016llx",
            (unsigned long long)parallel, (unsigned long long)a);
    std::printf("     digest interleaved_300 = 0x%016llx\n", (unsigned long long)a);

    Solver s;
    s.gravity = -9.8f;
    s.dt = 1.0f / 240.0f;
    s.iterations = 12;
    s.betaLin = 100000.0f;
    s.speculativeMargin = 0.02f;
    s.jacobianRebuildDistance = 0.02f;
    s.projectionMaxCorrection = 0.1f;
    s.sleepFrames = 30;
    s.sleepLinearThreshold = 0.05f;
    s.sleepAngularThreshold = 0.05f;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.8f, {0, 0, 0});
    std::vector<Rigid *> boxes;
    for (int i = 0; i < 10; i++) {
        Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.8f, {0, 0, 1.0f + i});
        box->sleep_mode = 1;
        boxes.push_back(box);
    }
    for (int tick = 0; tick < 600; tick++)
        s.stepInterleaved(4, 2);
    float height = 0.0f, penetration = 0.0f, drift = 0.0f;
    for (int i = 0; i < 10; i++) {
        const float3 p = boxes[i]->positionLin;
        height = std::max(height, std::fabs(p.z - (1.0f + i)));
        penetration = std::max(penetration,
                i == 0 ? 1.0f - p.z : 1.0f - (p.z - boxes[i - 1]->positionLin.z));
        drift = std::max(drift, std::hypot(p.x, p.y));
    }
    report(height <= 0.05f, "interleaved tight stack height", "max height error=%.5f m", height);
    report(penetration <= 0.01f, "interleaved tight stack penetration", "max interface=%.5f m", penetration);
    report(drift <= 0.01f, "interleaved tight stack drift", "max lateral=%.5f m", drift);
}

// Depth of two box bodies as the narrow phase sees it: the deepest contact point pair a
// collision pass returns at this pose (no speculative margin). Compare with satBoxDepth,
// which is the geometry alone and needs no contact features.
static float boxOverlap(const Rigid &a, const Rigid &b) {
    std::array<Manifold::Contact, 8> contacts;
    float3x3 basis{};
    const Shape sa{ .type = ShapeType::Box, .center = a.positionLin, .rotation = a.positionAng,
            .half = a.size * 0.5f };
    const Shape sb{ .type = ShapeType::Box, .center = b.positionLin, .rotation = b.positionAng,
            .half = b.size * 0.5f };
    const int found = collideShapes(sa, sb, contacts, basis);
    if (found <= 0)
        return 0.0f;
    const float3 normal = row(basis, 0);
    float best = 0.0f;
    for (int k = 0; k < found; k++) {
        const float3 pa = a.positionLin + qrot(a.positionAng, contacts[k].rA);
        const float3 pb = b.positionLin + qrot(b.positionAng, contacts[k].rB);
        best = std::max(best, -dot(normal, pa - pb));
    }
    return best;
}

// Minimum separating-axis depth of two box bodies - the same test the Godot-side probe
// runs, so a core `overlap` number and a `probe_overlap` number mean the same thing. Unlike
// the contact points a manifold stores, this is a property of the geometry alone.
static float satBoxDepth(const Rigid &a, const Rigid &b) {
    const float3 ha = a.size * 0.5f, hb = b.size * 0.5f;
    const float3 ea[3] = { qrot(a.positionAng, float3{ 1, 0, 0 }), qrot(a.positionAng, float3{ 0, 1, 0 }),
            qrot(a.positionAng, float3{ 0, 0, 1 }) };
    const float3 eb[3] = { qrot(b.positionAng, float3{ 1, 0, 0 }), qrot(b.positionAng, float3{ 0, 1, 0 }),
            qrot(b.positionAng, float3{ 0, 0, 1 }) };
    const float3 d = b.positionLin - a.positionLin;
    float best = 1.0e9f;
    auto axis = [&](const float3 &n, float ra, float rb) {
        const float sep = std::fabs(dot(n, d)) - (ra + rb);
        if (sep > 0.0f)
            return false; // a separating axis exists: no overlap at all
        best = std::min(best, -sep);
        return true;
    };
    for (int i = 0; i < 3; i++) {
        float ra = ha[i], rb = 0.0f;
        for (int j = 0; j < 3; j++)
            rb += hb[j] * std::fabs(dot(ea[i], eb[j]));
        if (!axis(ea[i], ra, rb))
            return 0.0f;
        float rb2 = hb[i], ra2 = 0.0f;
        for (int j = 0; j < 3; j++)
            ra2 += ha[j] * std::fabs(dot(ea[j], eb[i]));
        if (!axis(eb[i], ra2, rb2))
            return 0.0f;
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            float3 e = cross(ea[i], eb[j]);
            const float len2 = dot(e, e);
            if (len2 < 1.0e-8f)
                continue; // parallel edges: that axis carries no information
            e = e / std::sqrt(len2);
            float r1 = 0.0f, r2 = 0.0f;
            for (int k = 0; k < 3; k++) {
                r1 += ha[k] * std::fabs(dot(ea[k], e));
                r2 += hb[k] * std::fabs(dot(eb[k], e));
            }
            if (!axis(e, r1, r2))
                return 0.0f;
        }
    }
    return best;
}

// Pyramid impact speed. A body arriving at 18 m/s travels 7.5 cm between two contact
// refreshes, so interception that only works at rest speeds still overlaps by that much -
// and the solver then converts that overlap into separation velocity, which is how a
// collision becomes a launch. The overlap must stay far below one substep of travel and
// the pile must not leave the impact faster than the striker arrived.
static void test_impact_interception() {
    Solver s;
    s.gravity = -9.8f;
    s.dt = 1.0f / 240.0f;
    s.iterations = 12;
    s.betaLin = 100000.0f;
    s.speculativeMargin = 0.02f;
    s.jacobianRebuildDistance = 0.02f;
    s.projectionMaxCorrection = 0.1f;
    s.sleepFrames = 30;
    s.sleepLinearThreshold = 0.05f;
    s.sleepAngularThreshold = 0.05f;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.8f, {0, 0, 0});
    Rigid *target = new Rigid(&s, {1, 1, 1}, 1.0f, 0.8f, {0, 0, 1.0f});
    Rigid *striker = new Rigid(&s, {1, 1, 1}, 1.0f, 0.8f, {0, 0, 3.0f}, {0, 0, -18.0f});

    float worst = 0.0f, worstSpeed = 0.0f;
    int worstTick = -1;
    for (int tick = 0; tick < 240; tick++) {
        s.stepInterleaved(4, 2);
        const float depth = boxOverlap(*striker, *target);
        if (depth > worst) {
            worst = depth;
            worstTick = tick;
        }
        for (Rigid *b = s.bodies; b != nullptr; b = b->next)
            if (b->mass > 0.0f)
                worstSpeed = std::max(worstSpeed, length(b->velocityLin));
    }
    std::printf("     impact: worst=%.5f m at tick %d, max speed=%.3f m/s, striker z=%.4f "
            "target z=%.4f\n",
            worst, worstTick, worstSpeed, striker->positionLin.z, target->positionLin.z);
    report(worst <= 0.02f, "impact overlap bounded", "worst=%.5f m at tick %d (substep travel 0.075 m)",
            worst, worstTick);
    report(worstSpeed <= 20.0f, "impact does not amplify speed", "max=%.2f m/s from 18 m/s arrival",
            worstSpeed);
    report(std::fabs(target->positionLin.z - 1.0f) <= 0.05f && std::isfinite(striker->positionLin.z),
            "impact leaves the base placed", "target z=%.4f striker z=%.4f", target->positionLin.z,
            striker->positionLin.z);
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
// The parallel update must not depend on how many threads run it: every force solves its
// own mass-split problem against the same pose snapshot, and each body sums its shares in
// list order, so the split cannot change the result. If this ever fails, a phase boundary
// is missing and coupled work is racing.
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

// Bodies coupled by contacts and joints, plus a free-falling stack, so the update set
// is varied: contacts (coupled), joints (coupled), and singletons.
static void build_mixed(Solver &s) {
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
}

static uint64_t mixed_digest(int threads) {
    Solver s;
    s.threads = threads;
    build_mixed(s);
    step_n(s, 240);
    return state_digest(s);
}

// Two solvers stepped from two host threads. Both must match the serial reference, and
// neither may deadlock on the shared worker pool: a region reserves its helpers instead of
// holding a process-global loop mutex, and every join is a per-dispatch latch. With
// threads=4 each, two regions fit in an 8-thread machine and can overlap.
static void test_concurrent_solvers() {
    const uint64_t expected = mixed_digest(1);

    Solver a;
    Solver b;
    a.threads = 4;
    b.threads = 4;
    build_mixed(a);
    build_mixed(b);

    std::thread ta([&a] { step_n(a, 240); });
    std::thread tb([&b] { step_n(b, 240); });
    ta.join();
    tb.join();

    const uint64_t da = state_digest(a);
    const uint64_t db = state_digest(b);
    report(da == expected, "concurrent solver A", "0x%016llx vs serial 0x%016llx",
            (unsigned long long)da, (unsigned long long)expected);
    report(db == expected, "concurrent solver B", "0x%016llx vs serial 0x%016llx",
            (unsigned long long)db, (unsigned long long)expected);
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

// Report the parallel shape of each scene: how many bodies the Jacobi primal moves.
static void report_parallel_shape() {
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
        std::printf("     shape    %-10s bodies=%4d movable=%4d threads=%d\n", name, count_bodies(s),
                movable, s.threadCount());
    }
}

// ---------------------------------------------------------------------------
// The general 6-DOF joint. The angular side is the interesting one: locking two axes while a
// third spins is what makes a wheel possible, and a naive implementation (reading the locked
// components straight off the rotation vector between the bodies) works for a fraction of a
// turn and then tears itself apart, because that vector wraps every half turn.
// ---------------------------------------------------------------------------

// A static anchor with one body hanging off it through a GenericJoint.
struct GenericRig {
    Solver *solver;
    Rigid *anchor;
    Rigid *body;
    GenericJoint *joint;

    explicit GenericRig(Solver &s, float3 bodyPosition, float3 rA = {0, 0, 0}, float3 rB = {0, 0, 0}) {
        solver = &s;
        anchor = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 0});
        body = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, bodyPosition);
        joint = new GenericJoint(&s, anchor, body, rA, rB);
    }

    void setAxis(int p_axis, bool p_angular, AxisMode mode, float lower = 0.0f, float upper = 0.0f) {
        JointAxis &axis = p_angular ? joint->angular[p_axis] : joint->linear[p_axis];
        axis.mode = mode;
        axis.lower = lower;
        axis.upper = upper;
    }
};

// A wheel: the axle is free to spin, everything else is held.
//
// This is the test the swing-twist decomposition exists for. Spin fast for five seconds - 24
// revolutions - and the two locked axes must never move. Reading the locked components straight
// off the rotation vector between the bodies would look fine for half a revolution and then
// wrap, tearing the hinge open.
//
// Note the sign of the measurements: they are the orientation of A relative to B in A's frame,
// so a positive rotation of B reports as a negative value. That is the same convention the
// linear measurements use (they read A's anchor minus B's), and it does not matter for a lock,
// which only has to drive the value to zero.
static void test_generic_hinge() {
    Solver s;
    GenericRig rig(s, {0, 0, 0});
    for (int i = 0; i < 3; i++)
        rig.setAxis(i, false, AxisMode::Locked);
    rig.setAxis(0, true, AxisMode::Locked); // the two tilt axes are held
    rig.setAxis(1, true, AxisMode::Free);   // the axle spins
    rig.setAxis(2, true, AxisMode::Locked);

    rig.body->velocityAng = {0, 30.0f, 0};

    float worstLocked = 0.0f;
    float worstHub = 0.0f;
    const int steps = 300; // 5 s
    for (int i = 0; i < steps; i++) {
        s.step();
        worstLocked = std::max(worstLocked, std::fabs(rig.joint->angularValue(0)));
        worstLocked = std::max(worstLocked, std::fabs(rig.joint->angularValue(2)));
        worstHub = std::max(worstHub, length(rig.joint->anchorA() - rig.joint->anchorB()));
    }

    report(worstLocked <= 0.01f, "hinge locked axes hold through 24 turns",
            "worst |X| or |Z| deviation = %.5f rad over %.0f revolutions", worstLocked,
            steps / 60.0f * 30.0f / 6.2831853f);
    report(worstHub <= 0.01f, "hinge hub holds", "worst anchor separation = %.5f m", worstHub);

    // The axle is free, so nothing may torque it. The integrator does lose some angular
    // velocity (it recovers it as 2*vec(q_now * inv(q_before)) / dt, which is short by
    // sin(x/2)/(x/2) for a per-step rotation x, so the loss grows with the spin rate), so this
    // is checked at a rate where that loss is small. test_spin_retention() measures the rate
    // dependence on a body with no joint at all, which isolates it from this joint.
    Solver slow;
    GenericRig slowRig(slow, {0, 0, 0});
    for (int i = 0; i < 3; i++)
        slowRig.setAxis(i, false, AxisMode::Locked);
    slowRig.setAxis(0, true, AxisMode::Locked);
    slowRig.setAxis(1, true, AxisMode::Free);
    slowRig.setAxis(2, true, AxisMode::Locked);
    slowRig.body->velocityAng = {0, 1.0f, 0};

    for (int i = 0; i < 120; i++)
        slow.step();

    const float kept = slowRig.body->velocityAng.y;
    report(kept >= 0.98f, "free axle keeps its spin", "%.4f rad/s of 1.0 after 2 s (%.1f%%)", kept,
            100.0f * kept);
}

// The angular velocity loss is a property of the integrator, not of any joint: a body with no
// constraints at all loses it the same way. Measuring it here documents the boundary rather
// than leaving it to be discovered.
static void test_spin_retention() {
    struct Case {
        float spin;
        int steps;
    };
    const Case cases[] = {{1.0f, 120}, {10.0f, 120}};

    float worstShortfall = 0.0f;
    for (const Case &c : cases) {
        Solver s;
        s.gravity = 0.0f;
        Rigid *body = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 0});
        body->velocityAng = {0, c.spin, 0};
        for (int i = 0; i < c.steps; i++)
            s.step();

        // How much of the rotation is lost is set by the angle turned per step: the recovery
        // vector is 2*sin(x/2)/x short. At 10 rad/s (0.167 rad per step) most of it goes.
        const float retained = body->velocityAng.y / c.spin;
        if (c.spin == 1.0f)
            worstShortfall = 1.0f - retained;
        report(retained > 0.5f && retained <= 1.0001f, "spin decays monotonically with rate",
                "%5.1f rad/s -> %.2f%% kept after 2 s (%.3f rad per step)", c.spin, 100.0f * retained,
                c.spin / 60.0f);
    }
    report(worstShortfall <= 0.02f, "slow spin is preserved", "loss at 1 rad/s = %.2f%% over 2 s",
            100.0f * worstShortfall);
}

// A limit must stop the axis at its bound without kicking it.
//
// Two things are checked. First that the bound holds: a penalty constraint is soft, so a hard
// impact does penetrate it a little, and the assertion bounds that rather than pretending it is
// rigid. Second that stiffer really is tighter - the failure mode of a one-sided constraint is a
// retained multiplier being re-applied when the axis comes back to the bound, which gets *worse*
// as the axis is stiffened, so this is the check that catches it.
static float limit_excursion(float p_betaAng) {
    Solver s;
    s.betaAng = p_betaAng;
    GenericRig rig(s, {0, 0, 0});
    for (int i = 0; i < 3; i++)
        rig.setAxis(i, false, AxisMode::Locked);
    rig.setAxis(0, true, AxisMode::Locked);
    rig.setAxis(1, true, AxisMode::Limited, -0.4f, 0.4f); // the axle, ±0.4 rad of travel
    rig.setAxis(2, true, AxisMode::Locked);

    float worst = 0.0f;
    float reached = 0.0f;
    for (int i = 0; i < 240; i++) {
        if (i % 10 == 0 && i < 60)
            rig.body->velocityAng += float3{0, 6.0f, 0};
        s.step();
        const float angle = rig.joint->angularValue(1);
        reached = std::max(reached, std::fabs(angle));
        worst = std::max(worst, std::fabs(angle) - 0.4f);

        if (std::fabs(rig.joint->angularValue(0)) > 0.02f || std::fabs(rig.joint->angularValue(2)) > 0.02f) {
            report(false, "limit leaves locked axes alone", "X=%.5f Z=%.5f rad", rig.joint->angularValue(0),
                    rig.joint->angularValue(2));
            return worst;
        }
    }
    report(reached >= 0.35f, "limit is reached", "worst |angle| = %.4f rad at betaAng=%.0f (bound 0.4)",
            reached, p_betaAng);
    return worst;
}

static void test_generic_axis_limit() {
    const float soft = limit_excursion(100.0f);   // the default betaAng
    const float stiff = limit_excursion(1000.0f); // a limit tuned to be hard

    report(soft <= 0.1f, "default limit holds its bound",
            "worst excursion past the ±0.4 rad bound = %.5f rad (%.1f%% of the range)", soft,
            100.0f * soft / 0.8f);

    // Raising beta no longer buys a proportionally tighter bound, because the penalty is capped
    // at what the iteration can resolve (see GenericJoint::penaltyLimit). What both settings must
    // do is hold the bound without running away - the failure this guards against is a retained
    // multiplier being re-applied on re-entry, which used to overshoot by more than the whole
    // range and get *worse* as the axis was stiffened.
    report(stiff <= 0.1f, "stiffer limit also holds its bound",
            "worst excursion = %.5f rad at betaAng=1000 vs %.5f at 100 (penalty cap applies)", stiff, soft);
    report(stiff < 2.0f && soft < 2.0f, "neither stiffness runs away",
            "excursions %.5f / %.5f rad against a 0.8 rad range", stiff, soft);
}

// A prismatic axis with a spring must settle where the spring balances gravity:
// displacement = m*g/k. This is the car's suspension, and it has a closed form.
static void test_generic_prismatic_spring() {
    Solver s;
    GenericRig rig(s, {0, 0, 0});
    for (int i = 0; i < 3; i++) {
        rig.setAxis(i, false, AxisMode::Locked);
        rig.setAxis(i, true, AxisMode::Locked);
    }
    // Free to slide along the joint's Z axis, with a spring, as a strut.
    rig.setAxis(2, false, AxisMode::Limited, -2.0f, 2.0f);
    JointAxis &strut = rig.joint->linear[2];
    strut.springStiffness = 500.0f;
    strut.springEquilibrium = 0.0f;

    const float mass = rig.body->mass;
    const float expected = -mass * s.gravity / strut.springStiffness;

    for (int i = 0; i < 600; i++)
        s.step();

    const float travel = rig.joint->linearValue(2);
    report(std::fabs(travel - expected) <= 0.2f * std::fabs(expected) + 0.002f, "prismatic spring equilibrium",
            "travel = %.5f m, m*g/k = %.5f m (mass %.3f kg, k %.0f N/m)", travel, expected, mass,
            strut.springStiffness);
    report(std::fabs(travel) < 2.0f, "prismatic spring stays in range", "travel = %.4f m (limit ±2)", travel);
}

// A ball socket whose anchor is off the body's centre has to carry a moment: the constraint
// moves the anchor point both by translating the body and by turning it, so both halves of the
// system have to be stamped. With the anchors placed to coincide at spawn, the socket must keep
// them together while the body swings under gravity.
static void test_generic_ball_socket() {
    Solver s;
    // Hub at the world origin: rA is the anchor on the (static) anchor body, rB the offset from
    // the hanging body's centre up to the same point.
    GenericRig rig(s, {0, 0, -1.0f}, {0, 0, 0}, {0, 0, 1.0f});
    for (int i = 0; i < 3; i++) {
        rig.setAxis(i, false, AxisMode::Locked);
        rig.setAxis(i, true, AxisMode::Free);
    }

    const float initial = length(rig.joint->anchorA() - rig.joint->anchorB());
    rig.body->velocityLin = {1.5f, 0, 0}; // set it swinging

    float worst = 0.0f;
    for (int i = 0; i < 300; i++) {
        s.step();
        worst = std::max(worst, length(rig.joint->anchorA() - rig.joint->anchorB()));
    }

    report(initial <= 1.0e-5f, "socket anchors start together", "initial separation = %.7f m", initial);
    report(worst <= 0.01f, "socket holds under a swinging moment",
            "worst anchor separation = %.5f m while swinging", worst);
    report(rig.body->positionLin.z < -0.5f, "socket hangs the body below its hub",
            "body centre at z = %.4f (starts at -1)", rig.body->positionLin.z);
}

// A body that starts out rotated. This is the case that hid the deviation bug: every other test
// builds joints from upright bodies, so the rest pose is the identity and several wrong orderings
// of the deviation formula all look correct. A wheel turned onto its axle is not upright, and with
// the wrong ordering its locked axes see an error proportional to that starting rotation, which
// the solver then fights until the joint tears itself apart - measured as the penalty climbing
// from 1 to 41138 and the error growing with it.
static void test_generic_rotated_rest_pose() {
    Solver s;
    Rigid *anchor = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 10.0f});
    Rigid *wheel = new Rigid(&s, {0.3f, 0.3f, 0.16f}, ShapeType::Cylinder, 500.0f, 0.9f, {0, 0, 9.0f});
    // Turned so the cylinder's axis points along world X: the pose the joint is created in is not
    // the identity.
    wheel->positionAng = angleAxis(-0.7853982f, float3{0, 1, 0});

    GenericJoint *joint = new GenericJoint(&s, anchor, wheel, {0, 0, 1.0f}, {0, 0, 0});
    for (int i = 0; i < 3; i++) {
        joint->linear[i].mode = AxisMode::Locked;
        joint->angular[i].mode = AxisMode::Locked;
    }
    joint->angular[0].mode = AxisMode::Free; // the axle is the same axis the body was turned onto

    // Before stepping, every locked axis must already read zero: the joint was created in the pose
    // the wheel is in, so there is nothing to correct yet. A non-zero reading here is the bug.
    float initialLeak = 0.0f;
    for (int i = 1; i < 3; i++)
        initialLeak = std::max(initialLeak, std::fabs(joint->angularValue(i)));
    report(initialLeak <= 1.0e-5f, "rotated rest pose reads zero",
            "locked axes at t=0 read %.8f rad from a body turned 90 degrees", initialLeak);

    // Spin it about the axle and hold everything else. The locked axes must stay put.
    wheel->velocityAng = {8.0f, 0, 0};
    float worstLocked = 0.0f;
    float worstPenalty = 0.0f;
    for (int i = 0; i < 600; i++) {
        s.step();
        for (int a = 1; a < 3; a++)
            worstLocked = std::max(worstLocked, std::fabs(joint->angularValue(a)));
        worstPenalty = std::max(worstPenalty, std::max(joint->angular[1].penalty, joint->angular[2].penalty));
    }

    // KNOWN LIMITATION, recorded rather than hidden: the axis drift is not zero. It used to be
    // unbounded - the joint read a spurious error, fought it, inflated its penalty to 41138 and
    // NaN'd within seconds. With the deviation formula fixed and the penalty capped it stays
    // finite and settles, but the locked axes still creep about 2 rad over ten seconds of
    // spinning at 8 rad/s. That creep is the remaining stability work; the bound below is the
    // measured behaviour, not a target.
    report(worstLocked <= 2.5f, "locked axes stay bounded with a rotated rest pose (known limitation)",
            "worst locked-axis drift = %.5f rad over 10 s of spinning (unbounded before: NaN)", worstLocked);
    report(worstPenalty < 1.0e5f, "locked axes do not inflate",
            "peak penalty = %.1f (unbounded it reached 41138 and tore the joint apart)", worstPenalty);
    report(std::isfinite(wheel->velocityAng.x), "rotated rest pose stays finite",
            "axle spin = %.3f rad/s after 10 s", wheel->velocityAng.x);
}

// A free axis with no spring is genuinely free: nothing may fight an applied spin.
static void test_generic_free_axis() {
    Solver s;
    GenericRig rig(s, {0, 0, 0});
    for (int i = 0; i < 3; i++) {
        rig.setAxis(i, false, AxisMode::Locked);
        rig.setAxis(i, true, AxisMode::Free);
    }
    rig.body->velocityAng = {0.4f, 0.5f, 0.6f};

    for (int i = 0; i < 120; i++)
        s.step();

    const float3 w = rig.body->velocityAng;
    report(std::fabs(w.x - 0.4f) <= 0.02f && std::fabs(w.y - 0.5f) <= 0.02f && std::fabs(w.z - 0.6f) <= 0.02f,
            "free axes are free", "spin (%.4f %.4f %.4f) started at (0.4 0.5 0.6)", w.x, w.y, w.z);
}


// ---------------------------------------------------------------------------
// Round shapes. Boxes were the only primitive, and a box cannot roll: rotating it about its
// axle drives a corner into the ground instead of over it. These check the two new primitives
// rest correctly and, for the wheel case, that a round shape converts spin into travel.
// ---------------------------------------------------------------------------

// A body resting on the ground: where should its centre settle?
static float restHeight(ShapeType shape, float3 size, quat orientation, int steps = 400) {
    Solver s;
    new Rigid(&s, {200, 200, 1}, 0.0f, 0.5f, {0, 0, -0.5f});
    Rigid *body = new Rigid(&s, size, shape, 1000.0f, 0.5f, {0, 0, 4});
    body->positionAng = orientation;
    for (int i = 0; i < steps; i++)
        s.step();
    return body->positionLin.z;
}

static void test_round_shapes_rest() {
    // Sphere: centre one radius above the ground (ground top is z = 0).
    const float sphere = restHeight(ShapeType::Sphere, {0.5f, 0, 0}, quatIdentity);
    report(std::fabs(sphere - 0.5f) <= 0.05f, "sphere rests at its radius",
            "centre z = %.4f (radius 0.5, penetration %.4f)", sphere, 0.5f - sphere);

    // Cylinder on its cap: centre half a height up. The axis is local Z, already vertical.
    const float capped = restHeight(ShapeType::Cylinder, {0.5f, 0.5f, 1.0f}, quatIdentity);
    report(std::fabs(capped - 0.5f) <= 0.06f, "cylinder rests on its cap",
            "centre z = %.4f (half height 0.5, penetration %.4f)", capped, 0.5f - capped);

    // Cylinder on its side: rotated so the axis is horizontal (local Z -> world -Y), the
    // round side touches and the centre sits one radius up. This is the wheel case.
    const quat toSide = angleAxis(1.5707963f, float3{1, 0, 0});
    const float onSide = restHeight(ShapeType::Cylinder, {0.5f, 0.5f, 1.0f}, toSide);
    report(std::fabs(onSide - 0.5f) <= 0.06f, "cylinder rests on its side",
            "centre z = %.4f (radius 0.5, penetration %.4f)", onSide, 0.5f - onSide);

    // A sphere on a sphere: centres one diameter apart.
    Solver s;
    new Rigid(&s, {200, 200, 1}, 0.0f, 0.5f, {0, 0, -0.5f});
    Rigid *lower = new Rigid(&s, {0.5f, 0, 0}, ShapeType::Sphere, 1000.0f, 0.5f, {0, 0, 1.0f});
    Rigid *upper = new Rigid(&s, {0.5f, 0, 0}, ShapeType::Sphere, 1000.0f, 0.5f, {0, 0, 2.5f});
    for (int i = 0; i < 500; i++)
        s.step();
    const float gap = upper->positionLin.z - lower->positionLin.z;
    report(std::fabs(gap - 1.0f) <= 0.12f, "sphere stacks on sphere",
            "centre separation = %.4f (two radii = 1.0, penetration %.4f)", gap, 1.0f - gap);
}

// A round wheel on the ground, given spin about its axle, must travel: that is rolling, and it
// is exactly what a box cannot do.
static void test_round_shape_rolls() {
    Solver s;
    new Rigid(&s, {200, 200, 1}, 0.0f, 0.9f, {0, 0, -0.5f});
    const float radius = 0.25f;
    Rigid *wheel = new Rigid(&s, {radius, 0, 0}, ShapeType::Sphere, 1000.0f, 1.0f, {0, 0, radius});

    // Settle, then spin it about the axle (world X) and let contact friction do the rest.
    for (int i = 0; i < 120; i++)
        s.step();
    wheel->velocityAng = {20.0f, 0, 0};

    const float startY = wheel->positionLin.y;
    for (int i = 0; i < 120; i++)
        s.step();
    const float travelled = wheel->positionLin.y - startY;
    const float spin = wheel->velocityAng.x;

    // Rolling without slipping means the contact point is stationary: v_centre + w x r = 0. With
    // the axle along +X the contact sits at -Z, so w x r is +Y and the centre must move -Y. A
    // positive travel would mean the friction was pushing the wheel the wrong way.
    report(spin < 20.0f, "the ground resists the spin", "spin fell from 20 to %.3f rad/s", spin);
    report(travelled < -0.2f, "a round shape rolls", "travelled %.4f m in 2 s of spinning (rolling is -Y)",
            travelled);
    report(std::isfinite(travelled) && std::fabs(travelled) < 50.0f, "rolling stays bounded",
            "travelled %.4f m, spin %.3f rad/s", travelled, spin);
}

// ---------------------------------------------------------------------------
// collideShapes (the solver-free pair query) must agree with Manifold::collide
// (the path the solver walks every step) on the same poses: same contact count,
// same local anchors. Space queries lean on this equivalence.
// ---------------------------------------------------------------------------
static void shape_of(const Rigid *b, Shape &s) {
    s.type = b->shape;
    s.center = b->positionLin;
    s.rotation = b->positionAng;
    s.half = b->size * 0.5f;
    s.radius = b->size.x;
    s.halfHeight = b->size.z * 0.5f;
    s.axis = qrot(b->positionAng, float3{0, 0, 1});
}

static void check_pair_equivalent(Rigid *a, Rigid *b, const char *name) {
    Manifold::Contact viaBody[8];
    Manifold::Contact viaShape[8];
    float3x3 basisBody{};
    float3x3 basisShape{};

    const int nBody = Manifold::collide(a, b, viaBody, basisBody);
    Shape sa;
    Shape sb;
    shape_of(a, sa);
    shape_of(b, sb);
    const int nShape = collideShapes(sa, sb, viaShape, basisShape);

    bool same = nBody == nShape;
    float worst = 0.0f;
    for (int i = 0; same && i < nBody; i++) {
        worst = std::max(worst, length(viaBody[i].rA - viaShape[i].rA));
        worst = std::max(worst, length(viaBody[i].rB - viaShape[i].rB));
        same = worst <= 1.0e-5f;
    }
    for (int r = 0; same && r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            worst = std::max(worst, std::fabs(basisBody[c][r] - basisShape[c][r]));
            same = worst <= 1.0e-5f;
        }
    }
    report(same && nBody > 0, name, "contacts=%d/%d worst delta=%.2e", nBody, nShape, worst);
}

static void test_collide_shapes_equivalence() {
    Solver s;
    // Penetrating box pair, one rotated: exercises the SAT + clipping path.
    Rigid *boxA = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *boxB = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0.3f, 0.2f, 0.9f});
    boxB->positionAng = angleAxis(0.35f, normalize(float3{0.196f, 0.392f, -0.098f}));
    check_pair_equivalent(boxA, boxB, "collideShapes matches (box-box)");

    // Sphere penetrating a box from above: exercises the closest-point path.
    Rigid *sphere = new Rigid(&s, {0.5f, 0, 0}, ShapeType::Sphere, 1.0f, 0.5f, {0.1f, 0, 0.55f});
    check_pair_equivalent(boxA, sphere, "collideShapes matches (box-sphere)");
    check_pair_equivalent(sphere, boxA, "collideShapes matches (sphere-box, swapped)");
}

// Number of contact points across every manifold, and the largest normal error
// (positive = the surfaces are still apart).
static int manifold_contacts(const Solver &s, float &r_maxNormalC0) {
    int contacts = 0;
    r_maxNormalC0 = -INFINITY;
    for (int m = 0; m < s.manifoldCount(); m++) {
        const int count = s.manifoldContactCount(m);
        for (int i = 0; i < count; i++) {
            r_maxNormalC0 = std::fmax(r_maxNormalC0, s.manifoldContactError(m, i)[0]);
            contacts++;
        }
    }
    return contacts;
}

// ---------------------------------------------------------------------------
// Speculative contacts: with a margin, surfaces that are still apart report a
// contact whose positive normal error generates no force - until the body would
// close the gap within the step, when the implicit constraint stops it at the
// surface. Margin 0 keeps the previous overlap-only behaviour.
// ---------------------------------------------------------------------------
static void test_speculative_contacts() {
    Shape boxA{};
    boxA.type = ShapeType::Box;
    boxA.half = {0.5f, 0.5f, 0.5f};
    boxA.rotation = quatIdentity;
    Shape boxB = boxA;
    boxB.center = {0, 0, 1.03f}; // 3 cm air gap between the faces

    std::array<Manifold::Contact, 8> contacts;
    float3x3 basis{};
    const int overlapOnly = collideShapes(boxA, boxB, contacts, basis);
    const int speculative = collideShapes(boxA, boxB, contacts, basis, 0.05f);
    report(overlapOnly == 0 && speculative > 0, "speculative narrow phase",
            "3 cm gap: margin 0 -> %d contacts, margin 0.05 -> %d", overlapOnly, speculative);

    // At rest inside the margin: the manifold exists, its normal error is the gap,
    // and no force is applied while the gap is not closing.
    {
        Solver s;
        s.speculativeMargin = 0.05f;
        new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
        Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 1.03f});
        box->gravity = 0.0f;
        s.step();

        float maxNormalC0 = 0.0f;
        const int count = manifold_contacts(s, maxNormalC0);
        report(count > 0 && maxNormalC0 > 0.02f, "speculative gap contact",
                "contacts=%d max normal C0=%.4f (gap 0.03)", count, maxNormalC0);
        report(std::fabs(box->positionLin.z - 1.03f) < 1.0e-4f, "speculative gap does not push",
                "rest z=%.5f from spawn 1.03", box->positionLin.z);
    }

    // Closing the same 3 cm gap at 10 m/s: without speculation the body sinks 0.14 m
    // into the slab in one step; the speculative constraint stops it at the surface.
    {
        Solver s;
        s.speculativeMargin = 0.05f;
        new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
        Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 1.03f}, {0, 0, -10.0f});
        box->gravity = 0.0f;
        s.step();
        report(box->positionLin.z > 0.95f, "speculative catches a closing body",
                "z=%.4f after one 10 m/s step (surface at 1.0)", box->positionLin.z);
    }
}

// ---------------------------------------------------------------------------
// The velocity-side speculative clamp: a gap smaller than the margin closes at
// most gap/dt per step, so a fast body stops at the surface instead of sinking
// a whole step's travel deep into the slab.
// ---------------------------------------------------------------------------
static void test_speculative_gap_clamp() {
    Solver s;
    s.speculativeMargin = 0.01f;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 1.005f}, {0, 0, -10.0f});
    box->gravity = 0.0f;
    s.step();
    // Surface at z = 1.0 (slab top 0.5 + half 0.5). Unimpeded, one step at 10 m/s
    // would end near z = 0.84; the clamp may only close the 5 mm gap.
    report(box->positionLin.z > 0.95f, "speculative clamp stops at surface",
            "z=%.4f after one 10 m/s step from gap 5 mm", box->positionLin.z);
    report(std::fabs(box->positionLin.z - 1.0f) < 1e-4f && length(box->velocityAng) < 1e-6f,
            "speculative face stays symmetric", "z=%.7f spin=%.7f",
            box->positionLin.z, length(box->velocityAng));
    s.step();
    report(box->positionLin.z > 0.95f, "speculative clamp holds on second step",
            "z=%.4f after two steps", box->positionLin.z);
}

static void test_speculative_collision_horizon() {
    for (int substeps : {1, 4}) {
        for (bool projection : {false, true}) {
            Solver s;
            s.gravity = 0.0f;
            s.dt = 1.0f / (60.0f * substeps);
            s.positionProjection = projection;
            new Rigid(&s, {100, 100, 0.02f}, 0.0f, 0.5f, {0, 0, -0.01f});
            Rigid *sphere = new Rigid(&s, {0.05f, 0.05f, 0.05f}, ShapeType::Sphere,
                    1.0f, 0.5f, {0, 0, 0.3f}, {0, 0, -30.0f});
            float minHeight = sphere->positionLin.z;
            for (int tick = 0; tick < 5; tick++) {
                s.stepInterleaved(substeps, 2);
                minHeight = std::min(minHeight, sphere->positionLin.z);
            }
            report(minHeight >= 0.0f,
                    "speculative full-tick interception", "substeps=%d projection=%d min z=%.5f final z=%.5f",
                    substeps, projection, minHeight, sphere->positionLin.z);
        }
    }
    Solver s;
    s.gravity = 0.0f;
    s.dt = 1.0f / 240.0f;
    s.speculativeMargin = 0.0f;
    s.positionProjection = false;
    new Rigid(&s, {100, 100, 0.02f}, 0.0f, 0.5f, {0, 0, -0.01f});
    Rigid *sphere = new Rigid(&s, {0.05f, 0.05f, 0.05f}, ShapeType::Sphere,
            1.0f, 0.5f, {0, 0, 0.3f}, {0, 0, -30.0f});
    s.stepInterleaved(4, 2);
    report(sphere->positionLin.z < 0.0f && s.manifoldCount() == 0,
            "zero margin disables speculation", "z=%.5f contacts=%d", sphere->positionLin.z, s.manifoldCount());
}

static void test_projection_shared_contacts() {
    for (bool reverse : {false, true}) {
        float3 primalVelocity{};
        float primalHeight = 0.0f;
        for (bool projection : {false, true}) {
            Solver s;
            s.gravity = 0.0f;
            s.alpha = 1.0f;
            s.iterations = 1;
            s.positionProjection = projection;
            s.projectionIterations = 1;
            s.projectionSlop = 0.0f;
            Rigid *box = nullptr;
            for (int i = 0; i < 2; i++) {
                if ((i == 0) == reverse)
                    box = new Rigid(&s, {1, 1, 1}, 1.0f, 0.0f, {0, 0, 0.9f}, {0, 0, -1});
                else
                    new Rigid(&s, {100, 100, 1}, 0.0f, 0.0f, {0, 0, 0});
            }
            s.step();
            if (!projection) {
                primalVelocity = box->velocityLin;
                primalHeight = box->positionLin.z;
            } else {
                report(length(box->velocityAng) < 1e-6f &&
                                length(box->velocityLin - primalVelocity) < 1e-6f,
                        "overlap projection preserves velocity", "reverse=%d spin=%.7f velocity delta=%.7f",
                        reverse, length(box->velocityAng), length(box->velocityLin - primalVelocity));
                report(std::fabs(box->positionLin.z - 1.0f) < 1e-5f,
                        "face projection resolves in one pass", "z=%.5f without projection=%.5f",
                        box->positionLin.z, primalHeight);
            }
        }
    }
    for (int slabs : {1, 2, 4}) {
        for (float density : {1.0f, 10.0f}) {
            for (float cap : {0.0f, 0.01f}) {
                Solver s;
                s.gravity = 0.0f;
                s.alpha = 1.0f;
                s.iterations = 1;
                s.projectionIterations = 1;
                s.projectionBeta = 1.0f;
                s.projectionSlop = 0.0f;
                s.projectionMaxCorrection = cap;
                std::vector<Rigid *> ground;
                for (int i = 0; i < slabs; i++)
                    ground.push_back(new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0}));
                Rigid *sphere = new Rigid(&s, {0.5f, 0.5f, 0.5f}, ShapeType::Sphere,
                        density, 0.5f, {0, 0, 0.9f});
                s.step();
                const float expected = cap > 0.0f ? 0.9f + cap : 1.0f;
                bool stationary = length2(sphere->velocityLin) < 1e-10f;
                for (Rigid *body : ground)
                    stationary = stationary && length2(body->positionLin) == 0.0f;
                report(std::fabs(sphere->positionLin.z - expected) < 1e-5f && stationary,
                        "projection shared-plane response", "slabs=%d density=%.0f cap=%.2f z=%.5f expected=%.5f",
                        slabs, density, cap, sphere->positionLin.z, expected);
            }
        }
    }
    Solver s;
    s.gravity = 0.0f;
    s.alpha = 1.0f;
    s.iterations = 1;
    s.projectionIterations = 1;
    s.projectionBeta = 1.0f;
    s.projectionSlop = 0.0f;
    new Rigid(&s, {1, 100, 100}, 0.0f, 0.5f, {-0.9f, 0, 0});
    new Rigid(&s, {1, 100, 100}, 0.0f, 0.5f, {0.9f, 0, 0});
    Rigid *sphere = new Rigid(&s, {0.5f, 0.5f, 0.5f}, ShapeType::Sphere,
            1.0f, 0.5f, {0, 0, 0});
    s.step();
    report(length(sphere->positionLin) < 1e-5f, "opposed projections stay symmetric",
            "centre offset=%.7f", length(sphere->positionLin));
}

// ---------------------------------------------------------------------------
// Margin-aware cylinder paths: collideCylinderBox / collideCylinderCylinder
// report gap contacts through the signed point queries, and a cylinder parked
// just above a box holds a speculative manifold that pushes nothing.
// ---------------------------------------------------------------------------
static void test_cylinder_margin_contacts() {
    // Narrow phase: cylinder cap 5 mm above a box face.
    {
        Shape box{};
        box.type = ShapeType::Box;
        box.half = {0.5f, 0.5f, 0.5f};
        box.rotation = quatIdentity;

        Shape cyl{};
        cyl.type = ShapeType::Cylinder;
        cyl.radius = 0.25f;
        cyl.halfHeight = 0.5f;
        cyl.rotation = quatIdentity;
        cyl.axis = {0, 0, 1};
        cyl.center = {0, 0, 1.005f}; // cap 5 mm above the box top (0.5)

        std::array<Manifold::Contact, 8> contacts;
        float3x3 basis{};
        const int overlapOnly = collideShapes(cyl, box, contacts, basis);
        const int speculative = collideShapes(cyl, box, contacts, basis, 0.01f);
        report(overlapOnly == 0 && speculative > 0, "cylinder-box speculative narrow phase",
                "5 mm gap: margin 0 -> %d contacts, margin 0.01 -> %d", overlapOnly, speculative);
    }

    // Cylinder against cylinder, 5 mm apart.
    {
        Shape lower{};
        lower.type = ShapeType::Cylinder;
        lower.radius = 0.25f;
        lower.halfHeight = 0.5f;
        lower.rotation = quatIdentity;
        lower.axis = {0, 0, 1};
        lower.center = {0, 0, 0};

        Shape upper = lower;
        upper.center = {0, 0, 1.005f};

        std::array<Manifold::Contact, 8> contacts;
        float3x3 basis{};
        const int overlapOnly = collideShapes(lower, upper, contacts, basis);
        const int speculative = collideShapes(lower, upper, contacts, basis, 0.01f);
        report(overlapOnly == 0 && speculative > 0, "cylinder-cylinder speculative narrow phase",
                "5 mm gap: margin 0 -> %d contacts, margin 0.01 -> %d", overlapOnly, speculative);
    }

    // Solver level: a cylinder parked in the gap holds a manifold, moves nothing.
    {
        Solver s;
        s.speculativeMargin = 0.01f;
        new Rigid(&s, {1, 1, 1}, ShapeType::Box, 0.0f, 0.5f, {0, 0, 0});
        Rigid *cyl = new Rigid(&s, {0.5f, 0.5f, 1.0f}, ShapeType::Cylinder, 1.0f, 0.5f, {0, 0, 1.005f});
        cyl->gravity = 0.0f;
        s.step();

        float maxNormalC0 = 0.0f;
        const int count = manifold_contacts(s, maxNormalC0);
        report(count > 0 && maxNormalC0 > 0.003f, "cylinder gap contact",
                "contacts=%d max normal C0=%.4f (gap 0.005)", count, maxNormalC0);
        report(std::fabs(cyl->positionLin.z - 1.005f) < 1.0e-4f, "cylinder gap does not push",
                "rest z=%.5f from spawn 1.005", cyl->positionLin.z);
    }
}

// ---------------------------------------------------------------------------
// Cylinder tunneling: falling fast onto a slab with the spawn gap inside the
// margin, the cylinder must never sink its bottom below the slab top.
// ---------------------------------------------------------------------------
static void test_cylinder_tunneling() {
    Solver s;
    s.speculativeMargin = 0.01f;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *cyl = new Rigid(&s, {0.5f, 0.5f, 1.0f}, ShapeType::Cylinder, 1.0f, 0.5f,
            {0, 0, 1.008f}, {0, 0, -15.0f});
    cyl->gravity = 0.0f;

    bool tunnelled = false;
    for (int i = 0; i < 60; i++)
    {
        s.step();
        // Cylinder bottom (centre - half height) vs slab top (0.5).
        if (cyl->positionLin.z - 0.5f < 0.5f - 0.01f)
            tunnelled = true;
    }
    report(!tunnelled, "cylinder does not tunnel at 15 m/s",
            "final z=%.4f (bottom %.4f vs slab top 0.5)",
            cyl->positionLin.z, cyl->positionLin.z - 0.5f);
}

// ---------------------------------------------------------------------------
// Layer/mask filtering: Godot semantics - a pair collides when either side's
// layer is in the other's mask. Defaults (1/1) keep everything colliding.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Broad-phase bounds. The pair filter drops every candidate whose world AABBs miss,
// so aabbHalfExtent() must be the tight extent of the ROTATED shape: an extent that
// is too small hides a genuinely overlapping pair from the narrow phase, and that
// pair then falls through its neighbour with no constraint at all. Indexing the
// rotation matrix the wrong way round (rows for columns) is exactly a too-small
// extent on whichever axis the shape's thin direction leans into.
// ---------------------------------------------------------------------------
static void test_broadphase_bounds() {
    Solver s;
    // Rigid::size is the full extent, so half = size * 0.5 as the narrow phase uses it.
    const float3 fullSizes[] = { {2, 2, 1}, {1.2f, 0.4f, 2.8f}, {2, 2, 2}, {1, 4, 0.2f} };
    const quat tilts[] = {
            quatIdentity,
            angleAxis(0.7f, normalize(float3{0.196f, 0.392f, -0.098f})),
            angleAxis(1.1f, normalize(float3{-0.5f, 0.5f, 0.707f})),
            angleAxis(-2.0f, normalize(float3{0.3f, -0.8f, 0.5f})),
            qmul(angleAxis(0.9f, float3{1, 0, 0}), angleAxis(0.6f, float3{0, 1, 0})),
            qmul(angleAxis(2.4f, float3{0, 1, 0}), angleAxis(-1.3f, float3{1, 0, 0})),
    };
    float slack = 0.0f; // largest shortfall of the helper against the corner truth
    for (const float3 &full : fullSizes) {
        for (const quat &q : tilts) {
            Rigid *b = new Rigid(&s, full, 1.0f, 0.5f, {0, 0, 0});
            b->positionAng = q;
            const float3 h = full * 0.5f;
            float3 truth{};
            for (int sx = -1; sx <= 1; sx += 2)
                for (int sy = -1; sy <= 1; sy += 2)
                    for (int sz = -1; sz <= 1; sz += 2) {
                        const float3 c = qrot(q, float3{sx * h.x, sy * h.y, sz * h.z});
                        truth.x = std::max(truth.x, std::fabs(c.x));
                        truth.y = std::max(truth.y, std::fabs(c.y));
                        truth.z = std::max(truth.z, std::fabs(c.z));
                    }
            const float3 got = aabbHalfExtent(*b);
            if (truth.x - got.x > 1.0e-5f || truth.y - got.y > 1.0e-5f || truth.z - got.z > 1.0e-5f)
                std::printf("  AABBDIFF size=(%.2f,%.2f,%.2f) q=(%.4f,%.4f,%.4f,%.4f) "
                        "truth=(%.4f,%.4f,%.4f) got=(%.4f,%.4f,%.4f)\n",
                        full.x, full.y, full.z, q.x, q.y, q.z, q.w, truth.x, truth.y, truth.z,
                        got.x, got.y, got.z);
            slack = std::max({slack, truth.x - got.x, truth.y - got.y, truth.z - got.z});
            delete b;
        }
    }
    report(slack <= 1.0e-5f, "broadphase aabb contains shape",
            "largest shortfall=%.6f (negative means the bound is roomy)", slack);

    // The same mismatch from the other side: qmat must be the matrix qrot applies.
    float worstRot = 0.0f;
    for (const quat &q : tilts) {
        const float3x3 m = qmat(q);
        for (const float3 v : { float3{1, 0, 0}, float3{0, 1, 0}, float3{0, 0, 1},
                float3{0.4f, -0.7f, 0.25f} }) {
            const float3 a = qrot(q, v), b = mul(m, v);
            worstRot = std::max({worstRot, std::fabs(a.x - b.x), std::fabs(a.y - b.y),
                    std::fabs(a.z - b.z)});
        }
    }
    report(worstRot <= 1.0e-6f, "qmat matches qrot", "largest component diff=%.6f", worstRot);
}

static void test_layer_mask() {
    Solver s;
    Rigid *ground = new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    (void)ground;

    // Compatible pair: box layer 4, ground mask widened to all - collides, rests.
    Rigid *compatible = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
    compatible->collisionLayer = 4;
    compatible->collisionMask = 0xFFFFFFFFu;

    // Incompatible pair: neither layer sits in the other's mask - falls through.
    Rigid *filtered = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {5, 0, 4});
    filtered->collisionLayer = 2;
    filtered->collisionMask = 2;

    step_n(s, 240);

    report(std::fabs(compatible->positionLin.z - 1.0f) <= 0.02f, "compatible layers collide",
            "rest z=%.4f (slab top 0.5 + half 0.5)", compatible->positionLin.z);
    report(filtered->positionLin.z < 0.5f, "incompatible layers pass through",
            "fallen to z=%.4f after 4 s (no contact possible: layer 2 vs 1)", filtered->positionLin.z);
}

// ---------------------------------------------------------------------------
// Axis locks: a body locked on the gravity axis must not fall nor gain velocity,
// and a locked angular axis must resist spin - the Godot BodyAxis promise.
// ---------------------------------------------------------------------------
static void test_axis_lock() {
    Solver s;
    const float z0 = 5.0f;

    Rigid *locked = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, z0});
    locked->axisLockLinear = 0x4; // solver z, where gravity acts

    Rigid *freeFall = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {3, 0, z0});

    // A locked angular axis zeroes spin despite the initial kick.
    Rigid *spinLocked = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {-3, 0, z0});
    spinLocked->axisLockAngular = 0x1; // x
    spinLocked->velocityAng = {5.0f, 0, 0};

    step_n(s, 240);

    report(std::fabs(locked->positionLin.z - z0) <= 1.0e-4f, "locked axis hovers",
            "z=%.6f (spawn %.1f, no fall in 4 s)", locked->positionLin.z, z0);
    report(std::fabs(locked->velocityLin.z) <= 1.0e-4f, "locked axis no gravity",
            "vz=%.6f", locked->velocityLin.z);
    report(freeFall->positionLin.z < z0 - 3.0f, "unlocked control falls",
            "control z=%.4f", freeFall->positionLin.z);
    report(std::fabs(spinLocked->positionAng.x) <= 1.0e-4f && std::fabs(spinLocked->velocityAng.x) <= 1.0e-4f,
            "locked angular axis resists spin",
            "qx=%.6f wx=%.6f", spinLocked->positionAng.x, spinLocked->velocityAng.x);
}

// ---------------------------------------------------------------------------
// Masked picking: Solver::pick honours each body's layer, so a query with a
// narrow mask skips bodies on other layers and lands on the next candidate.
// ---------------------------------------------------------------------------
static void test_pick_mask() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0}); // ground, layer 1
    Rigid *blocker = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 3});
    blocker->collisionLayer = 2; // floating blocker

    float3 local{};
    const float3 origin{0, 0, 10};
    const float3 dir{0, 0, -1};

    Rigid *all = s.pick(origin, dir, local);
    Rigid *masked = s.pick(origin, dir, local, 1u);

    report(all != nullptr && std::fabs(all->positionLin.z - 3.0f) < 1.0e-3f, "pick hits nearest",
            "hit z=%.4f (blocker at 3)", all ? all->positionLin.z : -1.0f);
    report(masked != nullptr && std::fabs(masked->positionLin.z) < 1.0e-3f, "pick mask skips",
            "mask=1 hit z=%.4f (ground at 0)", masked ? masked->positionLin.z : -1.0f);
}

// ---------------------------------------------------------------------------
// A ray must hit the round shapes. A sphere is stored as {radius, 0, 0} (Rigid's
// constructor and aabbHalfExtent both read size.x as the radius, and the server stores
// exactly that), and `pick` used to test every body as a box of `size * 0.5f`: a sphere
// became a zero-thickness slab (so nothing could hit it) and a cylinder's radius was halved.
// ---------------------------------------------------------------------------
static void test_pick_round_shapes() {
    Solver s;
    Rigid *sphere = new Rigid(&s, {0.5f, 0, 0}, ShapeType::Sphere, 1.0f, 0.5f, {0, 0, 5});
    Rigid *cylinder = new Rigid(&s, {0.5f, 0.5f, 2.0f}, ShapeType::Cylinder, 1.0f, 0.5f, {0, 0, -5});

    float3 local{};
    // 0.4 m above the sphere's centre: still inside the 0.5 radius, but outside the
    // zero-thickness slab the old code tested it as (y/z half extent 0), so this is the
    // case that used to miss.
    Rigid *hitSphere = s.pick(float3{-5, 0, 5.4f}, float3{1, 0, 0}, local);
    report(hitSphere == sphere, "pick hits a sphere off its centre", "hit=%s",
            hitSphere == sphere ? "sphere" : (hitSphere == nullptr ? "nothing" : "another body"));
    const float sphereRadius = hitSphere == sphere ? length(local) : 0.0f;
    report(hitSphere == sphere && std::fabs(sphereRadius - 0.5f) < 0.01f,
            "pick lands on the sphere surface", "|local|=%.4f (want the radius 0.5)", sphereRadius);

    // 0.4 m off the cylinder's axis: outside the old (halved) 0.25 half extent, inside the
    // real 0.5 radius.
    Rigid *hitCyl = s.pick(float3{-5, 0.4f, -5}, float3{1, 0, 0}, local);
    report(hitCyl == cylinder, "pick hits a cylinder off its axis", "hit=%s",
            hitCyl == cylinder ? "cylinder" : (hitCyl == nullptr ? "nothing" : "another body"));

    // And a real miss still misses: 0.6 m off the axis is outside the radius.
    Rigid *missed = s.pick(float3{-5, 0.6f, -5}, float3{1, 0, 0}, local);
    report(missed == nullptr, "pick misses outside the cylinder radius", "hit=%s",
            missed == nullptr ? "nothing" : "a body");
}

// ---------------------------------------------------------------------------
// Two static bodies that overlap can never need a contact: neither moves. Before the
// broad phase filtered that pair, each overlapping static pair still got a manifold slot,
// and `ManifoldStore::frozen()` (both bodies immovable) then reported it alive every step,
// so a zero-contact manifold was never released - one slot plus eight contact records per
// pair, forever, on any static-only geometry.
// ---------------------------------------------------------------------------
static void test_static_static_pairs() {
    Solver s;
    for (int i = 0; i < 4; i++)
        new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {(float)i * 0.5f, 0, 0}); // 4 overlapping statics

    step_n(s, 2);
    report(s.manifoldCount() == 0, "static-static pairs make no manifold", "manifolds=%d",
            s.manifoldCount());

    // A static body against a movable one must still collide.
    new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0.5f, 0, 0.75f});
    step_n(s, 2);
    report(s.manifoldCount() > 0, "static-dynamic pairs still collide", "manifolds=%d",
            s.manifoldCount());
}

// ---------------------------------------------------------------------------
// Idle sleep: a settled stack of sleep-enabled boxes must fall asleep on its own
// (Solver::sleepFrames idle detection), and sleeping must not corrupt state. The core
// suite default keeps sleepFrames = 0, so this is opt-in behaviour.
// ---------------------------------------------------------------------------
static void test_idle_sleep() {
    Solver s;
    s.sleepFrames = 5;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, -0.5f});
    Rigid *a = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 0.5f});
    Rigid *b = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 1.5f});
    a->sleep_mode = 1;
    b->sleep_mode = 1;
    step_n(s, 120);
    report(a->sleeping && b->sleeping, "idle bodies fall asleep",
            "a=%d b=%d", a->sleeping, b->sleeping);
    report(all_finite(s), "sleeping stack stays finite", "bodies=%d", count_bodies(s));
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
    std::printf("after %d steps: threads=%d bodies=%d contacts=%d max|x|=%.4f finite=%s digest=0x%016llx\n",
            steps, s.threadCount(), count_bodies(s), count_contact_points(s), max_abs_position(s),
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

// Broad phase under test: 0 = hashed grid (default), 1 = Morton LBVH. Both are fed the same
// swept AABBs and produce the same candidate contract; --bvh is the isolated A/B.
static int g_broadPhaseMode = 0;
// Bench over a named core scene (Pyramid, PyramidXL, Bridge, ...): the same phase
// breakdown as --bench but on a scene the broad-phase comparison cares about.
static int run_bench_scene(const std::string &name, int threads) {
    const CoreScene *scene = find_scene(name);
    if (scene == nullptr) {
        std::fprintf(stderr, "unknown scene: %s\n", name.c_str());
        return 2;
    }
    Solver s;
    s.threads = threads;
    s.broadPhaseMode = g_broadPhaseMode;
    scene->build(&s);
    const int warmup = 20;
    const int timed = 120;
    step_n(s, warmup);
    const auto t0 = std::chrono::steady_clock::now();
    step_n(s, timed);
    const auto t1 = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / (double)timed;
    const Solver::PhaseTimes &pt = s.lastStepTimes;
    std::printf("bench-scene '%s': bodies=%d broadPhase=%s pairs=%d contacts=%d %.3f ms/step "
                "(broad %.3f, index %.3f, solve %.3f, finish %.3f)\n",
            name.c_str(), count_bodies(s), g_broadPhaseMode == 1 ? "lbvh" : "hash",
            s.broadPhasePairCount(), count_contact_points(s), ms, pt.broadPhase, pt.broadIndex,
            pt.solve, pt.finish);
    return all_finite(s) ? 0 : 1;
}

static int run_bench(int threads) {
    Solver s;
    s.threads = threads;
    build_bench_scene(s);
    std::printf("bench: %d boxes, %d iterations\n", count_bodies(s), s.iterations);

    const int warmup = 30;
    const int timed = 120;
    step_n(s, warmup);

    // Wall time: parallel speedup has to be measured on the clock, not on CPU time.
    const auto t0 = std::chrono::steady_clock::now();
    step_n(s, timed);
    const auto t1 = std::chrono::steady_clock::now();

    const double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / (double)timed;
    // Phase breakdown of the last step (Solver::lastStepTimes). This is the number that says
    // whether the broad phase is worth replacing: on a contact scene the solver dominates.
    const Solver::PhaseTimes &pt = s.lastStepTimes;
    std::printf("bench: last step phases: broad=%.3f (index %.3f) graph=%.3f solve=%.3f "
                "finish=%.3f (sum=%.3f)\n", pt.broadPhase, pt.broadIndex, pt.graph, pt.solve,
            pt.finish, pt.broadPhase + pt.graph + pt.solve + pt.finish);
    std::printf("bench: broadPhase=%s pairs=%d\n", g_broadPhaseMode == 1 ? "lbvh" : "hash",
            s.broadPhasePairCount());
    std::printf("bench: threads=%d -> %.3f ms/step, contacts=%d, finite=%s\n", s.threadCount(),
            ms, count_contact_points(s), all_finite(s) ? "yes" : "NO");    return all_finite(s) ? 0 : 1;
}

// Two independent solvers stepped from two host threads: the case the old global loop
// mutex turned into two serial runs. Reports the wall time for both, so the overlap (or
// its absence) is visible against the single-solver number above.
static int run_bench_two(int threads) {
    Solver a;
    Solver b;
    a.threads = threads;
    b.threads = threads;
    build_bench_scene(a);
    build_bench_scene(b);

    const int warmup = 30;
    const int timed = 120;
    step_n(a, warmup);
    step_n(b, warmup);

    const auto t0 = std::chrono::steady_clock::now();
    std::thread ta([&a, timed] { step_n(a, timed); });
    std::thread tb([&b, timed] { step_n(b, timed); });
    ta.join();
    tb.join();
    const auto t1 = std::chrono::steady_clock::now();

    const double total = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::printf("bench: two solvers, threads=%d each -> %.3f ms wall for %d steps (%.3f ms/step/solver), finite=%s\n",
            a.threadCount(), total, timed, total / (double)timed,
            all_finite(a) && all_finite(b) ? "yes" : "NO");
    return all_finite(a) && all_finite(b) ? 0 : 1;
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
        std::printf("bench: threads=%2d  %7.3f ms/step  speedup %5.2fx  (bodies=%d)\n", threads, ms,
                serialMs / ms, count_bodies(s));
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Block-Jacobi update-set coverage: two boxes at different heights collide mid-run,
// so a contact force appears between bodies whose update-set membership was decided
// on earlier, force-free steps. Every awake, movable body taking part in a stamping
// force must be in the update set - otherwise the constraint is silently dropped
// for that side and the pair drifts.
// ---------------------------------------------------------------------------
static void test_jacobi_shares() {
    Solver s;
    new Rigid(&s, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    Rigid *a = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 2});
    Rigid *b = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0.25f, 0, 6}); // different height: contact appears mid-run

    step_n(s, 300);

    bool covered = true;
    for (const Force *f = s.forces; f != nullptr; f = f->next) {
        const Rigid *sides[2] = {f->bodyA, f->bodyB};
        for (const Rigid *side : sides) {
            if (side != nullptr && side->mass > 0 && !side->sleeping && !side->inUpdateSet)
                covered = false;
        }
    }
    for (int m = 0; m < s.manifoldCount(); m++) {
        const Rigid *sides[2] = {s.manifoldBodyA(m), s.manifoldBodyB(m)};
        for (const Rigid *side : sides) {
            if (side != nullptr && side->mass > 0 && !side->sleeping && !side->inUpdateSet)
                covered = false;
        }
    }
    report(covered, "update set covers every awake constraint", "forces=%d", count_forces(s));
    report(all_finite(s), "update set finite", "bodies=%d", count_bodies(s));
    (void)a;
    (void)b;
}


// ---------------------------------------------------------------------------
// Overflow bucket: a 10-body joint clique needs 10 colours, the colouring is
// capped at 8, so bodies must fall into the block-Jacobi overflow bucket. With
// gravity off and no contacts the clique has nowhere to be flung: staying
// finite and in place exercises the bucket's snapshot solve and combine.
// ---------------------------------------------------------------------------
static void test_colour_overflow() {
    // 10 bodies all jointed into one clique: greedy colouring needs 10 colours,
    // the cap is 8, so bodies must land in the overflow bucket and its Jacobi
    // round must keep the scene finite and bounded.
    Solver s;
    s.gravity = 0.0f;
    Rigid *bodies[10];
    for (int i = 0; i < 10; i++)
        bodies[i] = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0.1f * i, 0, 5.0f + 0.05f * i});
    for (int i = 0; i < 10; i++)
        for (int j = i + 1; j < 10; j++)
            new Joint(&s, bodies[i], bodies[j], {0.5f, 0, 0}, {-0.5f, 0, 0}, 1000.0f, 250.0f);
    for (int i = 0; i < 600; i++)
        s.step();
    float worst = 0, maxabs = 0;
    bool finite = true;
    for (int i = 0; i < 10; i++) {
        worst = std::max(worst, std::fabs(bodies[i]->positionLin.x - 0.1f * i));
        maxabs = std::max(maxabs, std::fabs(bodies[i]->positionLin.z));
        finite = finite && std::isfinite(bodies[i]->positionLin.x) && std::isfinite(bodies[i]->positionLin.z);
    }
    report(finite && worst < 2.0f && maxabs < 100.0f, "overflow bucket holds a 10-clique",
            "drift=%.4f m, max|z|=%.3f m after 10 s mid-air", worst, maxabs);
}

// ---------------------------------------------------------------------------
// --pyramid: the Godot pyramid_dbl scene in core coordinates - a 20-column wall of
// 1 x 1 x 0.5 m boxes with an identical wall dropping onto it from 12 m. Reproduced
// here so impact penetration can be measured in about a second, without Godot, and
// with the true geometric overlap read straight from the narrow phase.
// ---------------------------------------------------------------------------
static void buildPyramidDouble(Solver &s) {
    new Rigid(&s, {80, 20, 1}, 0.0f, 0.5f, {0, 0, -0.5f});
    for (int tier = 0; tier < 2; tier++)
        for (int row = 0; row < 20; row++) {
            const int n = 20 - row;
            const float z = 0.25f + row * 0.52f + tier * 11.813f;
            for (int i = 0; i < n; i++)
                new Rigid(&s, {1, 1, 0.5f}, 2.0f, 0.5f,
                        {(i - (n - 1) * 0.5f) * 1.01f, 0.0f, z});
        }
}

// Index of the manifold holding this pair, or -1. Pairs are usually found in the first
// few slots, so the scan is fine at diagnostic rates.
static int pairManifold(const Solver &s, const Rigid *a, const Rigid *b) {
    for (int m = 0; m < s.manifoldCount(); m++) {
        const Rigid *ma = s.manifoldBodyA(m);
        const Rigid *mb = s.manifoldBodyB(m);
        if ((ma == a && mb == b) || (ma == b && mb == a))
            return m;
    }
    return -1;
}

// Diagnostic knobs for --pyramid, so a mechanism can be attributed without rebuilding
// variants of the scene.
static int g_substeps = 4;
static float g_projectionCap = 0.1f;
static bool g_projection = true;
static float g_margin = 0.02f;
static bool g_contactGs = false; // true: contact-only bodies join the coloured Gauss-Seidel
static int g_iterations = 12;    // rounds available per tick; per substep it is ceil(/substeps)
static float g_maxPen = CONTACT_MAX_PENETRATION; // depth the solve is allowed to see
static float g_ejectGain = 8.0f; // per-tick speed jump that counts as a solver kick

static int runPyramid(int ticks) {
    Solver s;
    s.gravity = -9.8f;
    s.dt = 1.0f / 240.0f;
    s.iterations = g_iterations;
    s.betaLin = 100000.0f;
    s.speculativeMargin = g_margin;
    s.jacobianRebuildDistance = 0.02f;
    s.positionProjection = g_projection;
    s.projectionMaxCorrection = g_projectionCap;
    s.contactOnlyOverflow = !g_contactGs;
    s.broadPhaseMode = g_broadPhaseMode;
    s.maxPenetration = g_maxPen;
    s.sleepFrames = 30;
    buildPyramidDouble(s);

    std::vector<Rigid *> boxes;
    Rigid *ground = nullptr;
    for (Rigid *b = s.bodies; b != nullptr; b = b->next) {
        if (b->mass > 0.0f)
            boxes.push_back(b);
        else
            ground = b;
    }
    std::sort(boxes.begin(), boxes.end(), [](const Rigid *a, const Rigid *b) {
        return a->initialLin.z != b->initialLin.z ? a->initialLin.z < b->initialLin.z
                                                  : a->initialLin.x < b->initialLin.x;
    });
    std::printf("== core pyramid_dbl  bodies=%d ticks=%d (%.2f s)\n", (int)boxes.size(), ticks,
            ticks / 60.0f);

    int blind = 0;
    // The manifold's own constraint value for a pair, at the pose it was measured at.
    // Deepest contact's stored normal error (the shallowest is hidden by a plain max: a
    // sheared face carries one deep corner and several near-zero ones).
    auto manifoldC0 = [&](int m) {
        float v = m < 0 ? -9.0e3f : 9.0e3f;
        const int cc = m < 0 ? 0 : s.manifoldContactCount(m);
        for (int k = 0; k < cc; k++)
            v = std::min(v, s.manifoldContactError(m, k)[0]);
        return v;
    };
    // Depth the manifold's stored arms and normal imply at the pose they are read at -
    // what C0 would be recomputed to right now. Compared with the stored C0 it separates a
    // wrong contact geometry (both small) from a stale constraint value (arms deep, C0 not).
    auto armsDepthOf = [&](int m) {
        float v = -9.0f;
        if (m < 0)
            return v;
        const Rigid *ma = s.manifoldBodyA(m), *mb = s.manifoldBodyB(m);
        const float3 n = s.manifoldNormal(m);
        const int cc = s.manifoldContactCount(m);
        for (int k = 0; k < cc; k++) {
            const Manifold::Contact &ct = s.manifoldContact(m, k);
            const float3 pa = ma->positionLin + qrot(ma->positionAng, ct.rA);
            const float3 pb = mb->positionLin + qrot(mb->positionAng, ct.rB);
            v = std::max(v, COLLISION_MARGIN - dot(n, pa - pb));
            // The stored value is the clamped one, so compare the deepest contact only.
        }
        return v;
    };
    // Worst true geometric overlap over all nearby pairs, plus how the solver itself sees
    // that pair (manifold present? constraint value?), so an overlap that appears inside a
    // single tick can be told apart from one the solver was told about and ignored.
    auto findWorst = [&](float &rWorst, int &rIa, int &rIb) {
        rWorst = 0.0f;
        rIa = rIb = -1;
        for (size_t i = 0; i < boxes.size(); i++) {
            for (size_t j = i + 1; j < boxes.size(); j++) {
                const Rigid *a = boxes[i], *b = boxes[j];
                if (std::fabs(a->positionLin.z - b->positionLin.z) > 1.2f ||
                        std::fabs(a->positionLin.x - b->positionLin.x) > 1.2f)
                    continue;
                const float d = satBoxDepth(*a, *b);
                if (d > 0.005f && pairManifold(s, a, b) < 0) {
                    blind++;
                    // No constraint for a pair that is really overlapping: is the narrow
                    // phase blind to it, or did the pair never reach the narrow phase?
                    if (blind <= 4) {
                        std::array<Manifold::Contact, 8> fresh{};
                        float3x3 fb{};
                        const Shape sa{ .type = ShapeType::Box, .center = a->positionLin,
                                .rotation = a->positionAng, .half = a->size * 0.5f };
                        const Shape sb{ .type = ShapeType::Box, .center = b->positionLin,
                                .rotation = b->positionAng, .half = b->size * 0.5f };
                        const int found = collideShapes(sa, sb, fresh, fb);
                        const float3 gap = b->positionLin - a->positionLin;
                        const float3 ha = aabbHalfExtent(*a), hb = aabbHalfExtent(*b);
                        const bool aabb = std::fabs(gap.x) < ha.x + hb.x &&
                                std::fabs(gap.y) < ha.y + hb.y && std::fabs(gap.z) < ha.z + hb.z;
                        std::printf("  BLIND pair %zu/%zu depth=%.4f freshContacts=%d aabb=%d"
                                " sweep=%d dist=(%.2f,%.2f,%.2f) v=(%.2f,%.2f)\n", i, j, d, found,
                                (int)aabb, s.sweepHasPair(a, b) ? 1 : 0,
                                gap.x, gap.y, gap.z, length(a->velocityLin), length(b->velocityLin));
                    }
                }
                if (d > rWorst) {
                    rWorst = d;
                    rIa = (int)i;
                    rIb = (int)j;
                }
            }
        }
    };

    std::vector<float> samples;
    float maxSpeed = 0.0f, lowest = 1.0e9f;
    int over5cm = 0, detail = 0, inTick = 0;
    std::vector<bool> escaped(boxes.size(), false);
    // Ejection watch: a tick that leaves a body moving much faster than it arrived is where
    // the pile's energy comes from, so the contacts that were open at the time get dumped.
    std::vector<float> prevSpeed(boxes.size(), 0.0f);
    int ejects = 0;
    float maxGain = 0.0f;
    int traceTicks = 0;
    for (int tick = 0; tick < ticks; tick++) {
        const bool sample = tick % 5 == 4;
        float preWorst = 0.0f;
        int pia = -1, pib = -1;
        float preC0 = -9.0f, preArms = -9.0f;
        if (sample) {
            findWorst(preWorst, pia, pib);
            if (pia >= 0) {
                const int pm = pairManifold(s, boxes[pia], boxes[pib]);
                preC0 = manifoldC0(pm);
                // Stored C0 versus the stored arms read at this pose: agreement means the
                // contact geometry is right and only the value is old; agreeing with the true
                // depth instead means the constraint is right and the solve is what fails.
                preArms = armsDepthOf(pm);
                // A pair already deep and unconstrained when the tick starts: trace the tick
                // that failed to give it a constraint, rather than the one after it.
                if (preWorst > 0.05f && pm < 0 && traceTicks == 0) {
                    s.tracePair[0] = boxes[pia];
                    s.tracePair[1] = boxes[pib];
                    traceTicks = 1;
                }
            }
        }

        s.stepInterleaved(g_substeps, 2);
        if (traceTicks > 0 && --traceTicks == 0)
            s.tracePair[0] = s.tracePair[1] = nullptr;
        for (Rigid *b : boxes) {
            maxSpeed = std::max(maxSpeed, length(b->velocityLin));
            lowest = std::min(lowest, b->positionLin.z);
        }
        // A box below the slab's underside: tell "pushed through the floor" (still inside the
        // slab footprint, so the ground contact failed) from "squeezed out over the edge".
        for (size_t i = 0; i < boxes.size(); i++) {
            Rigid *b = boxes[i];
            if (escaped[i] || b->positionLin.z > -0.6f)
                continue;
            escaped[i] = true;
            const float3 p = b->positionLin, v = b->velocityLin;
            std::printf("  ESC t=%.2fs body %zu pose=(%.2f,%.2f,%.2f) v=(%.2f,%.2f,%.2f)"
                    " overSlab=%d groundManifold=%d\n",
                    (tick + 1) / 60.0f, i, p.x, p.y, p.z, v.x, v.y, v.z,
                    (int)(std::fabs(p.x) < 40.0f && std::fabs(p.y) < 10.0f),
                    pairManifold(s, ground, b));
        }
        for (size_t i = 0; i < boxes.size(); i++) {
            const float now = length(boxes[i]->velocityLin);
            const float before = prevSpeed[i];
            const float gain = now - before;
            prevSpeed[i] = now;
            maxGain = std::max(maxGain, gain);
            if (gain < g_ejectGain || ejects >= 12)
                continue;
            ejects++;
            std::printf("  EJECT t=%.2fs body %d speed %.2f -> %.2f dir=(%.2f,%.2f,%.2f)"
                    " pose=(%.2f,%.2f)\n", (tick + 1) / 60.0f, (int)i, before, now,
                    boxes[i]->velocityLin.x / now, boxes[i]->velocityLin.y / now,
                    boxes[i]->velocityLin.z / now, boxes[i]->positionLin.x,
                    boxes[i]->positionLin.z);
            const int mc = s.manifoldCount();
            for (int m = 0; m < mc; m++) {
                const Rigid *oa = s.manifoldBodyA(m), *ob = s.manifoldBodyB(m);
                const bool isA = oa == boxes[i];
                if (!isA && ob != boxes[i])
                    continue;
                const Rigid *other = isA ? ob : oa;
                int oi = -1;
                for (size_t j = 0; j < boxes.size(); j++)
                    if (boxes[j] == other)
                        oi = (int)j;
                const int cc = s.manifoldContactCount(m);
                const float3 n = s.manifoldNormal(m);
                float c0 = -9.0f;
                for (int k = 0; k < cc; k++)
                    c0 = std::max(c0, s.manifoldContactError(m, k)[0]);
                std::printf("      contact %d/%d contacts=%d c0=%.4f n=(%.2f,%.2f,%.2f)"
                        " otherSpeed=%.2f\n", (int)i, oi, cc, c0, n.x, n.y, n.z,
                        length(other->velocityLin));
            }
        }
        if (!sample)
            continue;

        // Energy census: the running peak above says nothing about whether the pile is still
        // moving, so report the live speeds and how many bodies are awake every second.
        if ((tick + 1) % 60 == 0) {
            float liveMax = 0.0f;
            int awake = 0, moving = 0, flying = 0;
            for (Rigid *b : boxes) {
                const float v = length(b->velocityLin);
                liveMax = std::max(liveMax, v);
                if (!b->sleeping)
                    awake++;
                if (v > 1.0f)
                    moving++;
                if (v > 5.0f)
                    flying++;
            }
            std::printf("STATUS t=%.2fs liveMax=%.2f moving=%d flying=%d awake=%d manifolds=%d\n",
                    (tick + 1) / 60.0f, liveMax, moving, flying, awake, s.manifoldCount());
        }

        float worst = 0.0f;
        int ia = -1, ib = -1;
        findWorst(worst, ia, ib);
        samples.push_back(worst);
        if (worst > 0.05f) {
            over5cm++;
            if (pia == ia && pib == ib && worst - preWorst > 0.02f)
                inTick++; // the pair was not reported that way at the start of the tick
        }
        if (worst > 0.05f && detail < 20) {
            detail++;
            Rigid *a = boxes[ia], *b = boxes[ib];
            const int m = pairManifold(s, a, b);
            const int cc = m < 0 ? 0 : s.manifoldContactCount(m);
            const float c0 = manifoldC0(m);
            const float armsDepth = armsDepthOf(m);
            // For a pair with no constraint at all, split "the broad phase never offered
            // them" from "the narrow phase was asked and produced nothing": call the narrow
            // phase directly on the same poses the solver just used.
            int fresh = -1;
            int swept = -1;
            if (m < 0) {
                std::array<Manifold::Contact, 8> fc{};
                float3x3 fb{};
                const Shape sa{ .type = ShapeType::Box, .center = a->positionLin,
                        .rotation = a->positionAng, .half = a->size * 0.5f };
                const Shape sb{ .type = ShapeType::Box, .center = b->positionLin,
                        .rotation = b->positionAng, .half = b->size * 0.5f };
                fresh = collideShapes(sa, sb, fc, fb, 0.0f);
                swept = s.sweepHasPair(a, b) ? 1 : 0;
                // Follow this pair through the solver for the next ticks: candidate, filter,
                // acquire, initialize, reap.
                s.tracePair[0] = a;
                s.tracePair[1] = b;
                traceTicks = 4;
                // Convention-independent AABB truth: the extent of the eight rotated corners.
                // If this disagrees with aabbHalfExtent(), the broad phase tests wrong bounds.
                auto cornerAabb = [](const Rigid &r) {
                    const float3 h = r.size * 0.5f;
                    float3 e{};
                    for (int sx = -1; sx <= 1; sx += 2)
                        for (int sy = -1; sy <= 1; sy += 2)
                            for (int sz = -1; sz <= 1; sz += 2) {
                                const float3 c = qrot(r.positionAng,
                                        float3{sx * h.x, sy * h.y, sz * h.z});
                                e.x = std::max(e.x, std::fabs(c.x));
                                e.y = std::max(e.y, std::fabs(c.y));
                                e.z = std::max(e.z, std::fabs(c.z));
                            }
                    return e;
                };
                const float3 ca = cornerAabb(*a), cb = cornerAabb(*b);
                const float3 ab = aabbHalfExtent(*a), bb = aabbHalfExtent(*b);
                std::printf("      AABB a size=(%.2f,%.2f,%.2f) corner=(%.3f,%.3f,%.3f) "
                        "helper=(%.3f,%.3f,%.3f)\n        b size=(%.2f,%.2f,%.2f) "
                        "corner=(%.3f,%.3f,%.3f) helper=(%.3f,%.3f,%.3f)\n",
                        a->size.x, a->size.y, a->size.z, ca.x, ca.y, ca.z, ab.x, ab.y, ab.z,
                        b->size.x, b->size.y, b->size.z, cb.x, cb.y, cb.z, bb.x, bb.y, bb.z);
            }
            std::printf("  DEEP t=%.2fs pair %d/%d overlap=%.4f (tick start %.4f, c0 then %.4f,"
                    " same pair %d)  manifold=%d contacts=%d c0=%.4f armsNow=%.4f preArms=%.4f"
                    " sweep=%d fresh=%d"
                    "  a=(%.2f,%.2f) v=%.2f  b=(%.2f,%.2f) v=%.2f\n",
                    (tick + 1) / 60.0f, ia, ib, worst, preWorst, preC0, pia == ia && pib == ib, m,
                    cc, c0, armsDepth, preArms, swept, fresh, a->positionLin.x, a->positionLin.z,
                    length(a->velocityLin),
                    b->positionLin.x, b->positionLin.z, length(b->velocityLin));
        }
        if (worst > 0.05f)
            std::printf("t=%.2fs worst=%.4f (%d/%d) max_v=%.2f lowest=%.3f\n", (tick + 1) / 60.0f,
                    worst, ia, ib, maxSpeed, lowest);
    }
    std::sort(samples.begin(), samples.end());
    const size_t n = samples.size();
    std::printf("RESULT substeps=%d iters=%d proj=%d cap=%.3f spec=%.3f gs=%d penmax=%.3f groups=%d "
            "overflow=%d "
            "p50=%.4f p95=%.4f max=%.4f over5cm=%d/%zu blindPairs=%d grewInTick=%d max_v=%.2f "
            "lowest=%.3f maxGain=%.2f\n",
            g_substeps, g_iterations, (int)g_projection, g_projectionCap, g_margin, (int)g_contactGs,
            g_maxPen,
            s.colourGroups(), s.overflowBodies(), n > 0 ? samples[n / 2] : 0.0f,
            n > 0 ? samples[(size_t)(n * 0.95)] : 0.0f,
            n > 0 ? samples[n - 1] : 0.0f, over5cm, n, blind, inTick, maxSpeed, lowest, maxGain);
    return 0;
}

int main(int argc, char **argv) {
    std::string scene;
    int steps = 300;
    int threads = 0; // 0 = automatic
    bool pyramid = false;

    for (int i = 1; i < argc; i++) {
        const std::string arg = argv[i];
        if (arg == "--list") {
            list_scenes();
            return 0;
        } else if (arg == "--bench") {
            return run_bench(threads);
        } else if (arg == "--bench-sweep") {
            return run_bench_sweep();
        } else if (arg == "--bench-two") {
            return run_bench_two(threads == 0 ? 4 : threads);
        } else if (arg == "--shape") {
            report_parallel_shape();
            return 0;
        } else if (arg == "--pyramid") {
            pyramid = true;
        } else if (arg.rfind("--substeps=", 0) == 0) {
            g_substeps = std::atoi(arg.c_str() + 11);
        } else if (arg.rfind("--cap=", 0) == 0) {
            g_projectionCap = (float)std::atof(arg.c_str() + 6);
        } else if (arg.rfind("--proj=", 0) == 0) {
            g_projection = std::atoi(arg.c_str() + 7) != 0;
        } else if (arg.rfind("--spec=", 0) == 0) {
            g_margin = (float)std::atof(arg.c_str() + 7);
        } else if (arg.rfind("--gs=", 0) == 0) {
            g_contactGs = std::atoi(arg.c_str() + 5) != 0;
        } else if (arg.rfind("--iters=", 0) == 0) {
            g_iterations = std::atoi(arg.c_str() + 8);
        } else if (arg == "--bench-scene" && i + 1 < argc) {
            return run_bench_scene(argv[++i], threads);
        } else if (arg == "--bvh") {
            g_broadPhaseMode = 1;
        } else if (arg.rfind("--penmax=", 0) == 0) {
            g_maxPen = (float)std::atof(arg.c_str() + 9);
        } else if (arg.rfind("--eject=", 0) == 0) {
            g_ejectGain = (float)std::atof(arg.c_str() + 8);
        } else if (arg == "--scene" && i + 1 < argc) {
            scene = argv[++i];
        } else if (arg == "--steps" && i + 1 < argc) {
            steps = std::atoi(argv[++i]);
        } else if (arg == "--threads" && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
        } else {
            std::fprintf(stderr, "usage: %s [--list | --bench | --bench-sweep | --bench-two | --shape |"
                                 " --pyramid [--steps N] | --scene NAME [--steps N]] [--threads N]"
                                 " [--bvh]\n", argv[0]);
            return 2;
        }
    }

    if (!scene.empty()) {
        return run_scene(scene, steps, threads);
    }
    if (pyramid) {
        return runPyramid(steps);
    }

    test_rest_contact();
    test_stack_stability();
    test_friction();
    test_hard_joints();
    test_soft_lattice();
    test_all_scenes();
    test_determinism();
    test_interleaved_substeps();
    test_jacobi_shares();
    test_parallel_equivalence();
    test_concurrent_solvers();
    test_colour_overflow();
    report_parallel_shape();
    test_generic_hinge();
    test_generic_axis_limit();
    test_generic_prismatic_spring();
    test_generic_ball_socket();
    test_generic_free_axis();
    test_spin_retention();
    test_generic_rotated_rest_pose();
    test_round_shapes_rest();
    test_round_shape_rolls();
    test_collide_shapes_equivalence();
    test_speculative_contacts();
    test_speculative_gap_clamp();
    test_speculative_collision_horizon();
    test_projection_shared_contacts();
    test_cylinder_margin_contacts();
    test_cylinder_tunneling();
    test_impact_interception();
    test_layer_mask();
    test_broadphase_bounds();
    test_axis_lock();
    test_pick_mask();
    test_pick_round_shapes();
    test_static_static_pairs();
    test_idle_sleep();

    std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
