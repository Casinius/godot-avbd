/*
 * Solver-only benchmark: the grid scene from demo/tests/test_pulse_bench.gd, stepped
 * directly through avbd::Solver with no Godot scene tree, no node sync, and no frame
 * pacing. What it measures is the solver, so it is the fast loop for optimisation work;
 * use tools/pulse_bench.py when the question involves the engine or a Jolt comparison.
 *
 *   xmake f -m release && xmake avbd_bench
 *   ./build/linux/x86_64/release/avbd_bench --grid=23 --ticks=480 --steps=1200
 *   gprof ./build/linux/x86_64/release/avbd_bench gmon.out
 */
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "avbd/solver.h"
#include "core_scenes.hpp"

using namespace avbd;

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

double percentile(std::vector<double> values, double p) {
    if (values.empty())
        return 0.0;
    std::sort(values.begin(), values.end());
    const size_t index = static_cast<size_t>(p * static_cast<double>(values.size() - 1));
    return values[index];
}

bool matches(const char *arg, const char *name) {
    return std::strncmp(arg, name, std::strlen(name)) == 0;
}

} // namespace

int main(int argc, char **argv) {
    int grid = 23;
    int steps = 1200;
    int warmup = 120;
    int iterations = 4;
    int threads = 0;
    float ticks = 240.0f;
    float gap = 1.02f;
    bool adaptive = false;

    for (int i = 1; i < argc; ++i) {
        const char *a = argv[i];
        if (matches(a, "--grid="))
            grid = std::atoi(a + 7);
        else if (matches(a, "--steps="))
            steps = std::atoi(a + 8);
        else if (matches(a, "--warmup="))
            warmup = std::atoi(a + 9);
        else if (matches(a, "--iterations="))
            iterations = std::atoi(a + 13);
        else if (matches(a, "--threads="))
            threads = std::atoi(a + 10);
        else if (matches(a, "--ticks="))
            ticks = static_cast<float>(std::atof(a + 8));
        else if (matches(a, "--gap="))
            gap = static_cast<float>(std::atof(a + 6));
        else if (std::strcmp(a, "--adaptive") == 0)
            adaptive = true;
        else if (std::strcmp(a, "--help") == 0) {
            std::printf("usage: avbd_bench [--grid=N] [--ticks=HZ] [--steps=N] [--warmup=N]\n"
                        "                  [--iterations=N] [--threads=N] [--gap=F] [--adaptive]\n");
            return 0;
        } else {
            std::fprintf(stderr, "unknown argument: %s\n", a);
            return 2;
        }
    }

    if (grid < 2 || steps < 1 || ticks <= 0.0f) {
        std::fprintf(stderr, "grid >= 2, steps >= 1 and ticks > 0 are required\n");
        return 2;
    }

    Solver solver;
    solver.threads = threads;
    solver.iterations = iterations;
    solver.adaptiveIterations = adaptive;
    solver.dt = 1.0f / ticks;
    solver.profile = true;
    sceneBoxGrid(&solver, grid, gap);

    for (int i = 0; i < warmup; ++i)
        solver.step();

    std::vector<double> samples;
    samples.reserve(steps);
    Solver::Profile total;
    for (int i = 0; i < steps; ++i) {
        const Clock::time_point start = Clock::now();
        solver.step();
        samples.push_back(elapsed_ms(start));
        const Solver::Profile &p = solver.profileLast;
        total.broadPhase += p.broadPhase;
        total.warmstartForces += p.warmstartForces;
        total.colourGraph += p.colourGraph;
        total.warmstartBodies += p.warmstartBodies;
        total.solveIterations += p.solveIterations;
        total.finishVelocities += p.finishVelocities;
        total.rounds += p.rounds;
    }

    const double inv = 1.0 / static_cast<double>(steps);
    double mean = 0.0;
    for (double v : samples)
        mean += v;
    mean *= inv;
    const std::vector<double> ordered = samples;

    std::printf("avbd_bench: grid=%d bodies=%d ticks=%.0f steps=%d iterations=%d threads=%d adaptive=%s\n",
            grid, grid * grid, ticks, steps, iterations, threads, adaptive ? "yes" : "no");
    std::printf("avbd_bench: step mean=%.4f median=%.4f min=%.4f p95=%.4f ms\n",
            mean, percentile(ordered, 0.50), percentile(ordered, 0.0), percentile(ordered, 0.95));
    std::printf("avbd_bench: rounds mean=%.2f bodies=%d forces=%d colours=%d pairs=%d\n",
            static_cast<double>(total.rounds) * inv, solver.profileLast.bodies,
            solver.profileLast.forces, solver.profileLast.colours,
            solver.profileLast.pairCandidates);
    std::printf("avbd_bench: phases broad=%.4f warm_forces=%.4f colour=%.4f warm_bodies=%.4f "
                "solve=%.4f finish=%.4f ms\n",
            total.broadPhase * inv, total.warmstartForces * inv, total.colourGraph * inv,
            total.warmstartBodies * inv, total.solveIterations * inv,
            total.finishVelocities * inv);
    return 0;
}
