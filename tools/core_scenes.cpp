/*
 * Standard AVBD test scenes, ported from the authors' reference implementation
 * (github.com/savant117/avbd-demo3d).
 *
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 *
 * Units: metres, kilograms, seconds. Gravity acts along -Z; +Z is up.
 */

#include "core_scenes.hpp"

#include "avbd/maths.h"
#include "avbd/solver.h"

namespace avbd {

void sceneEmpty(Solver *solver)
{
    solver->clear();
}

void sceneGround(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0}, {0, 0, 0});
    new Rigid(solver, {1, 1, 1}, 1.0f, 0.5f, {0, 0, 4});
}

void sceneDynamicFriction(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0}, {0, 0, 0});
    for (int x = 0; x <= 10; x++)
        new Rigid(solver, {1, 1, 0.5f}, 1.0f, 5.0f - (x / 10.0f * 5.0f), {0, -30.0f + x * 2.0f, 0.75f}, {10.0f, 0, 0});
}

void sceneStaticFriction(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    const float angle = rad(30.0f);
    Rigid *ramp = new Rigid(solver, {40, 24, 1}, 0.0f, 1.0f, {0, 0, 3});
    ramp->positionAng = {0, sinf(angle * 0.5f), 0, cosf(angle * 0.5f)};

    float3 rampTangent = normalize(rotate(ramp->positionAng, float3{1, 0, 0}));
    float3 rampNormal = normalize(rotate(ramp->positionAng, float3{0, 0, 1}));

    for (int i = 0; i <= 10; i++)
    {
        float friction = i / 10.0f * 0.25f + 0.25f;
        float y = -10.0f + i * 2.0f;
        float3 pos = ramp->positionLin + rampTangent * -12.0f + float3{0, y, 0} + rampNormal * 1.05f;
        new Rigid(solver, {1, 1, 1}, 1.0f, friction, pos);
    }
}

void scenePyramid(Solver *solver)
{
    const int SIZE = 16;
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0.0f, 0.0f, -0.5f});

    for (int y = 0; y < SIZE; y++)
        for (int x = 0; x < SIZE - y; x++)
            new Rigid(solver, {1, 0.5f, 0.5f}, 1.0f, 0.5f, {x * 1.01f + y * 0.5f - SIZE / 2.0f, 0.0f, y * 0.85f + 0.5f});
}

void sceneRope(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, -20});

    Rigid *prev = 0;
    for (int i = 0; i < 20; i++)
    {
        Rigid *curr = new Rigid(solver, {1, 0.5f, 0.5f}, i == 0 ? 0.0f : 1.0f, 0.5f, {(float)i, 0.0f, 10.0f});
        if (prev)
            new Joint(solver, prev, curr, {0.5f, 0, 0}, {-0.5f, 0, 0});
        prev = curr;
    }
}

void sceneHeavyRope(Solver *solver)
{
    const int N = 20;
    const float SIZE = 5;
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, -20});

    Rigid *prev = 0;
    for (int i = 0; i < N; i++)
    {
        Rigid *curr = new Rigid(solver, i == N - 1 ? float3{SIZE, SIZE, SIZE} : float3{1, 0.5f, 0.5f},
                                i == 0 ? 0.0f : 1.0f, 0.5f, {(float)i + (i == N - 1 ? SIZE / 2 : 0), 0.0f, 10.0f});
        if (prev)
            new Joint(solver, prev, curr, {0.5f, 0, 0}, i == N - 1 ? float3{-SIZE / 2, 0, 0} : float3{-0.5f, 0, 0});
        prev = curr;
    }
}

void sceneSpring(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    Rigid *anchor = new Rigid(solver, {1, 1, 1}, 0.0f, 0.5f, {0, 0, 14.0f});
    Rigid *block = new Rigid(solver, {2, 2, 2}, 1.0f, 0.5f, {0, 0, 8.0f});
    new Spring(solver, anchor, block, {0, 0, 0}, {0, 0, 0}, 100.0f, 4.0f);
}

void sceneSpringsRatio(Solver *solver)
{
    const int N = 8;
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, -10});

    Rigid *prev = 0;
    for (int i = 0; i < N; i++)
    {
        float x = (i - (N - 1) * 0.5f) * 3.0f;
        Rigid *curr = new Rigid(solver, {1, 0.75f, 0.75f}, i == 0 || i == N - 1 ? 0.0f : 1.0f, 0.5f, {x, 0.0f, 12.0f});
        if (prev)
            new Spring(solver, prev, curr, {0.5f, 0, 0}, {-0.5f, 0, 0}, i % 2 == 0 ? 10.0f : 10000.0f, 3.0f);
        prev = curr;
    }
}

void sceneStack(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});
    for (int i = 0; i < 10; i++)
        new Rigid(solver, {1, 1, 1}, 1.0f, 0.5f, {0, 0, i * 1.5f + 1.0f});
}

void sceneStackRatio(Solver *solver)
{
    solver->clear();
    const float groundThickness = 1.0f;
    new Rigid(solver, {100, 100, groundThickness}, 0.0f, 0.5f, {0, 0, 0});

    float topZ = groundThickness * 0.5f;
    float s = 1.0f;
    for (int i = 0; i < 4; i++)
    {
        float half = s * 0.5f;
        float centerZ = topZ + half;
        new Rigid(solver, {s, s, s}, 1.0f, 0.5f, {0, 0, centerZ});
        topZ = centerZ + half;
        s *= 2.0f;
    }
}

void sceneSoftBody(Solver *solver)
{
    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    const float Klin = 1000.0f;
    const float Kang = 250.0f;
    const int W = 4;
    const int D = 4;
    const int H = 4;
    const int N = 3;
    const float size = 0.8f;
    const float half = size * 0.5f;
    const float baseZ = 8.0f;
    const float stackGap = 2.0f;

    for (int i = 0; i < N; i++)
    {
        Rigid *grid[W][D][H];
        float stackZ = i * (H * size + stackGap);

        for (int x = 0; x < W; x++)
        {
            for (int y = 0; y < D; y++)
            {
                for (int z = 0; z < H; z++)
                {
                    float px = (x - (W - 1) * 0.5f) * size;
                    float py = (y - (D - 1) * 0.5f) * size;
                    float pz = baseZ + stackZ + z * size;
                    grid[x][y][z] = new Rigid(solver, {size, size, size}, 1.0f, 0.5f, {px, py, pz});
                }
            }
        }

        for (int x = 1; x < W; x++)
        {
            for (int y = 0; y < D; y++)
            {
                for (int z = 0; z < H; z++)
                {
                    new Joint(solver, grid[x - 1][y][z], grid[x][y][z], {half, 0, 0}, {-half, 0, 0}, Klin, Kang);
                }
            }
        }

        for (int x = 0; x < W; x++)
        {
            for (int y = 1; y < D; y++)
            {
                for (int z = 0; z < H; z++)
                {
                    new Joint(solver, grid[x][y - 1][z], grid[x][y][z], {0, half, 0}, {0, -half, 0}, Klin, Kang);
                }
            }
        }

        for (int x = 0; x < W; x++)
        {
            for (int y = 0; y < D; y++)
            {
                for (int z = 1; z < H; z++)
                {
                    new Joint(solver, grid[x][y][z - 1], grid[x][y][z], {0, 0, half}, {0, 0, -half}, Klin, Kang);
                }
            }
        }

        for (int x = 1; x < W; x++)
        {
            for (int y = 0; y < D; y++)
            {
                for (int z = 1; z < H; z++)
                {
                    new IgnoreCollision(solver, grid[x - 1][y][z - 1], grid[x][y][z]);
                    new IgnoreCollision(solver, grid[x][y][z - 1], grid[x - 1][y][z]);
                }
            }
        }

        for (int x = 0; x < W; x++)
        {
            for (int y = 1; y < D; y++)
            {
                for (int z = 1; z < H; z++)
                {
                    new IgnoreCollision(solver, grid[x][y - 1][z - 1], grid[x][y][z]);
                    new IgnoreCollision(solver, grid[x][y][z - 1], grid[x][y - 1][z]);
                }
            }
        }

        for (int x = 1; x < W; x++)
        {
            for (int y = 1; y < D; y++)
            {
                for (int z = 0; z < H; z++)
                {
                    new IgnoreCollision(solver, grid[x - 1][y - 1][z], grid[x][y][z]);
                    new IgnoreCollision(solver, grid[x][y - 1][z], grid[x - 1][y][z]);
                }
            }
        }
    }
}

void sceneBridge(Solver *solver)
{
    const int N = 40;
    const float plankLength = 1.0f;
    const float plankWidth = 4.0f;
    const float plankHeight = 0.5f;
    const float halfLength = plankLength * 0.5f;
    const float halfWidth = plankWidth * 0.5f;

    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    Rigid *prev = 0;
    for (int i = 0; i < N; i++)
    {
        Rigid *curr = new Rigid(solver, {plankLength, plankWidth, plankHeight}, i == 0 || i == N - 1 ? 0.0f : 1.0f, 0.5f, {(float)i - N / 2.0f, 0.0f, 10.0f});
        if (prev)
        {
            new Joint(solver, prev, curr, {halfLength, halfWidth, 0}, {-halfLength, halfWidth, 0}, INFINITY, 0.0f);
            new Joint(solver, prev, curr, {halfLength, -halfWidth, 0}, {-halfLength, -halfWidth, 0}, INFINITY, 0.0f);
        }
        prev = curr;
    }

    for (int x = 0; x < N / 4; x++)
    {
        for (int y = 0; y < N / 8; y++)
        {
            new Rigid(solver, {1, 1, 1}, 1.0f, 0.5f, {(float)x - N / 8.0f, 0.0f, (float)y + 12.0f});
        }
    }
}

void sceneBreakable(Solver *solver)
{
    const int N = 10;
    const int M = 5;
    const float breakForce = 90.0f;

    solver->clear();
    new Rigid(solver, {100, 100, 1}, 0.0f, 0.5f, {0, 0, 0});

    Rigid *prev = 0;
    for (int i = 0; i <= N; i++)
    {
        Rigid *curr = new Rigid(solver, {1, 1, 0.5f}, 1.0f, 0.5f, {(float)i - N / 2.0f, 0.0f, 6.0f});
        if (prev)
            new Joint(solver, prev, curr, {0.5f, 0, 0}, {-0.5f, 0, 0}, INFINITY, INFINITY, breakForce);
        prev = curr;
    }

    new Rigid(solver, {1, 1, 5}, 0.0f, 0.5f, {-N / 2.0f, 0, 2.5f});
    new Rigid(solver, {1, 1, 5}, 0.0f, 0.5f, {N / 2.0f, 0, 2.5f});

    for (int i = 0; i < M; i++)
        new Rigid(solver, {2, 1, 1}, 1.0f, 0.5f, {0, 0, i * 2.0f + 8.0f});
}

const CoreScene coreScenes[] = {
    {"Empty", sceneEmpty},
    {"Ground", sceneGround},
    {"Dynamic Friction", sceneDynamicFriction},
    {"Static Friction", sceneStaticFriction},
    {"Pyramid", scenePyramid},
    {"Rope", sceneRope},
    {"Heavy Rope", sceneHeavyRope},
    {"Spring", sceneSpring},
    {"Spring Ratio", sceneSpringsRatio},
    {"Stack", sceneStack},
    {"Stack Ratio", sceneStackRatio},
    {"Soft Body", sceneSoftBody},
    {"Bridge", sceneBridge},
    {"Breakable", sceneBreakable},
};

const int coreSceneCount = 14;

} // namespace avbd
