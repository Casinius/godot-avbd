/*
 * Throwaway repro: a chain of fully-locked GenericJoint welds hanging from a static
 * anchor, run at the demo's substep schedule. Isolates the solver core from the Godot
 * server so the runaway can be attributed to the joint code.
 *
 *   g++ -O2 -std=c++23 -Isrc tools/repro_joint.cpp src/avbd/*.cpp -o /tmp/repro_joint -pthread
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "avbd/maths.h"
#include "avbd/solver.h"

using namespace avbd;

// One welded link pair: all six axes Locked.
static void lockAll(GenericJoint *j)
{
    for (int i = 0; i < 3; i++)
    {
        j->linear[i].mode = AxisMode::Locked;
        j->angular[i].mode = AxisMode::Locked;
    }
}

int main(int argc, char **argv)
{
    // Arguments: [substeps] [roundsPerSubstep] [links] [lockedAngular]
    const int substeps = argc > 1 ? std::atoi(argv[1]) : 4;
    const int rounds = argc > 2 ? std::atoi(argv[2]) : 3;
    const int linkCount = argc > 3 ? std::atoi(argv[3]) : 4;
    const int lockAngular = argc > 4 ? std::atoi(argv[4]) : 1;
    const float dtDiv = argc > 5 ? (float)std::atof(argv[5]) : (float)(60 * substeps);
    const float betaLin = argc > 6 ? (float)std::atof(argv[6]) : 100000.0f;

    Solver s;
    s.dt = 1.0f / dtDiv; // demo default: 60 Hz tick split across the substeps
    s.iterations = 12;
    s.betaLin = betaLin;  // demo project.godot ships 1e5
    s.betaAng = 100.0f;
    s.speculativeMargin = 0.02f;

    Rigid *anchor = new Rigid(&s, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 10});
    Rigid *previous = anchor;
    for (int i = 0; i < linkCount; i++)
    {
        Rigid *link = new Rigid(&s, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 9.0f - (float)i});
        GenericJoint *joint = new GenericJoint(&s, previous, link, {0, 0, -0.5f}, {0, 0, 0.5f});
        for (int a = 0; a < 3; a++)
        {
            joint->linear[a].mode = AxisMode::Locked;
            joint->angular[a].mode = lockAngular ? AxisMode::Locked : AxisMode::Free;
        }
        (void)lockAll;
        previous = link;
    }

    printf("substeps=%d rounds=%d links=%d lockAngular=%d dt=%.6f betaLin=%.0f\n", substeps, rounds,
            linkCount, lockAngular, (double)s.dt, (double)s.betaLin);

    // First joint: the one from the static anchor, so it carries the whole chain's weight.
    GenericJoint *j0 = nullptr;
    for (Force *f = s.forces; f != nullptr; f = f->next)
        if (auto *gj = dynamic_cast<GenericJoint *>(f)) { j0 = gj; break; }

    int tick = 0;
    for (; tick < 300; tick++)
    {
        if (substeps > 1)
            s.stepInterleaved(substeps, rounds);
        else
            s.step();
        if (tick < 5 || tick % 20 == 0)
        {
            printf("t=%4d", tick);
            for (Rigid *b = s.bodies; b != nullptr; b = b->next)
                printf("  z=%9.3f vz=%9.2f", (double)b->positionLin.z, (double)b->velocityLin.z);
            if (j0 != nullptr)
                printf("  | kZ=%.0f lamZ=%.3f C=%9.5f", (double)j0->linear[2].penalty,
                        (double)j0->linear[2].lambda,
                        (double)(j0->anchorA().z - j0->anchorB().z));
            printf("\n");
        }
    }
    return 0;
}
