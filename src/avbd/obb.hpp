/*
 * Conservative OBB overlap test for the broad phase.
 *
 * The spatial hash still buckets bodies by their tight world AABB (a cell must fit the
 * AABB at its level), but the exact candidate filter uses an oriented box instead of the
 * world AABB: a rotated box's AABB can be several times its true footprint, and every
 * false positive becomes a Manifold that collide() has to reject and warmstartForces()
 * has to delete.
 *
 * Spheres are represented by an axis-aligned cube of half extent = radius (a sphere is
 * rotation invariant, so this is its exact AABB). Cylinders use their local enclosing box.
 * Both enclosures contain the real shape, so the test can only keep too many pairs, never
 * drop a real overlap.
 */

#pragma once

#include "avbd/maths.h"

namespace avbd {

struct Obb
{
    float3 center;
    float3 axis[3]; // world-space unit axes (the rotation matrix's columns)
    float3 half;
    bool isSphere = false;
    float radius = 0.0f;
};

// 15-axis SAT (6 face axes + 9 edge axes). `axes` are unit and orthogonal, so the
// projections reduce to table reads; no normalisation is needed. A small epsilon biases
// the test towards overlap (a false positive is cheap, a false negative is a lost contact).
[[nodiscard]] inline bool obbOverlaps(const Obb &a, const Obb &b) noexcept
{
    constexpr float kSlop = 1.0e-5f;

    if (a.isSphere && b.isSphere)
    {
        const float r = a.radius + b.radius + kSlop;
        return length2(a.center - b.center) <= r * r;
    }

    // A sphere only needs the closest point on the other box; this is exact when the other
    // side is a box and conservative when it is a cylinder's enclosing box.
    if (a.isSphere || b.isSphere)
    {
        const Obb &sphere = a.isSphere ? a : b;
        const Obb &box = a.isSphere ? b : a;
        const float3 d = sphere.center - box.center;
        float3 u{dot(d, box.axis[0]), dot(d, box.axis[1]), dot(d, box.axis[2])};
        u.x = std::clamp(u.x, -box.half.x, box.half.x);
        u.y = std::clamp(u.y, -box.half.y, box.half.y);
        u.z = std::clamp(u.z, -box.half.z, box.half.z);
        const float3 closest = box.center + box.axis[0] * u.x + box.axis[1] * u.y + box.axis[2] * u.z;
        const float r = sphere.radius + kSlop;
        return length2(sphere.center - closest) <= r * r;
    }

    // R[i][j] = axisA_i . axisB_j, and the translation in A's frame.
    float r[3][3];
    float ar[3][3];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            r[i][j] = dot(a.axis[i], b.axis[j]);
            ar[i][j] = std::fabs(r[i][j]) + kSlop;
        }
    }

    const float3 d = b.center - a.center;
    const float t[3]{dot(d, a.axis[0]), dot(d, a.axis[1]), dot(d, a.axis[2])};

    // A's face axes.
    for (int i = 0; i < 3; i++)
    {
        const float projected = b.half.x * ar[i][0] + b.half.y * ar[i][1] + b.half.z * ar[i][2];
        if (std::fabs(t[i]) > a.half[i] + projected + kSlop)
            return false;
    }

    // B's face axes.
    for (int j = 0; j < 3; j++)
    {
        const float projected = a.half.x * ar[0][j] + a.half.y * ar[1][j] + a.half.z * ar[2][j];
        const float tb = t[0] * r[0][j] + t[1] * r[1][j] + t[2] * r[2][j];
        if (std::fabs(tb) > b.half[j] + projected + kSlop)
            return false;
    }

    // Edge axes A_i x B_j.
    for (int i = 0; i < 3; i++)
    {
        const int i1 = (i + 1) % 3;
        const int i2 = (i + 2) % 3;
        for (int j = 0; j < 3; j++)
        {
            const int j1 = (j + 1) % 3;
            const int j2 = (j + 2) % 3;
            const float ra = a.half[i1] * ar[i2][j] + a.half[i2] * ar[i1][j];
            const float rb = b.half[j1] * ar[i][j2] + b.half[j2] * ar[i][j1];
            const float projected = t[i2] * r[i1][j] - t[i1] * r[i2][j];
            if (std::fabs(projected) > ra + rb + kSlop)
                return false;
        }
    }
    return true;
}

} // namespace avbd
