/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#include <climits>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <span>
#include <vector>

#include "avbd/solver.h"

#include "avbd/job_pool.hpp"
#include "avbd/list_range.hpp"
#include "avbd/fixed_stack.hpp"
#include "radixsort.hpp"

namespace avbd {

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local)
{
    return pick(origin, dir, local, 0xFFFFFFFFu);
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local, uint32_t p_mask)
{
    const float epsilon = 1.0e-6f;
    float bestT = std::numeric_limits<float>::infinity();
    Rigid *bestBody = 0;
    float3 bestLocal = {0, 0, 0};

    // Exact per-body test against that body's own shape; picking scans bodies, not the LBVH.
    // Halve extents per shape: for a sphere size.x is already the radius, so size*0.5 made it unhittable.
    const auto localHalfExtent = [](const Rigid &body) noexcept {
        switch (body.shape)
        {
            case ShapeType::Sphere:
            {
                const float r = body.size.x;
                return float3{r, r, r};
            }
            case ShapeType::Cylinder:
                // Local +Z is the axis: radius spans X/Y, half the height spans Z.
                return float3{body.size.x, body.size.x, body.size.z * 0.5f};
            case ShapeType::Box:
            default:
                return body.size * 0.5f;
        }
    };

    const auto testBody = [&](Rigid *body) {
        if ((body->collisionLayer & p_mask) == 0)
            return;
        quat invRot = qconj(body->positionAng);
        float3 o = qrot(invRot, origin - body->positionLin);
        float3 d = qrot(invRot, dir);

        if (body->shape == ShapeType::Sphere)
        {
            // Exact ray/sphere. The box below would only be a proxy with corners that are
            // not on the surface.
            const float r = body->size.x;
            const float b = dot(o, d);
            const float c = length2(o) - r * r;
            const float disc = b * b - c;
            if (disc < 0.0f)
                return;
            const float sq = std::sqrt(disc);
            float t = -b - sq;
            if (t < 0.0f)
                t = -b + sq;
            if (t < 0.0f || t >= bestT)
                return;
            bestT = t;
            bestBody = body;
            bestLocal = o + d * t;
            return;
        }

        // Box, and the box that bounds a cylinder. `dir` is not normalised by the callers,
        // so `t` is in units of `dir` for every shape path.
        const float3 half = localHalfExtent(*body);

        float tEnter = 0.0f;
        float tExit = std::numeric_limits<float>::infinity();
        bool hit = true;

        for (int i = 0; i < 3; ++i)
        {
            if (std::fabs(d[i]) < epsilon)
            {
                if (o[i] < -half[i] || o[i] > half[i])
                {
                    hit = false;
                    break;
                }
                continue;
            }

            float invD = 1.0f / d[i];
            float t0 = (-half[i] - o[i]) * invD;
            float t1 = (half[i] - o[i]) * invD;
            if (t0 > t1)
            {
                float tmp = t0;
                t0 = t1;
                t1 = tmp;
            }

            tEnter = std::max(tEnter, t0);
            tExit = std::min(tExit, t1);
            if (tEnter > tExit)
            {
                hit = false;
                break;
            }
        }

        if (!hit)
            return;

        float tHit = tEnter >= 0.0f ? tEnter : tExit;
        if (tHit < 0.0f)
            return;

        if (tHit < bestT)
        {
            bestT = tHit;
            bestBody = body;
            bestLocal = o + d * tHit;
        }
    };

    // Unlike the upstream demo (which picks bodies for dragging) this also tests
    // static bodies, so rays hit the ground and ramps.
    for (Rigid *body = bodies; body != 0; body = body->next)
        testBody(body);

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
}

void Solver::collectAABB(float3 min, float3 max, uint32_t p_mask, std::vector<Rigid *> &out) const
{
    // Shape/point queries walk the body list: each body's bound is centre plus its tight
    // world half-extent (the same formula the broad phase derives its leaves from).
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if ((body->collisionLayer & p_mask) == 0)
            continue;
        const float3 wx = aabbHalfExtent(*body);
        const float3 bmin = body->positionLin - wx;
        const float3 bmax = body->positionLin + wx;
        if (!(bmin.x > max.x || bmin.y > max.y || bmin.z > max.z ||
                min.x > bmax.x || min.y > bmax.y || min.z > bmax.z))
            out.push_back(body);
    }
}


} // namespace avbd
