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

#include "avbd/solver.h"
#include <cmath>

namespace avbd {

Rigid::Rigid(Solver *p_solver, float3 p_size, float p_density, float p_friction, float3 p_position, float3 p_velocity)
    : Rigid(p_solver, p_size, ShapeType::Box, p_density, p_friction, p_position, p_velocity)
{
}

Rigid::Rigid(Solver *p_solver, float3 p_size, ShapeType p_shape, float p_density, float p_friction,
        float3 p_position, float3 p_velocity)
    : solver(p_solver), forces(0), next(0), positionLin(p_position), positionAng(quatIdentity),
    velocityLin(p_velocity), velocityAng({ 0, 0, 0 }), prevVelocityLin(p_velocity), size(p_size),
    shape(p_shape), friction(p_friction), gravity(p_solver->gravity)
{
    // Add to linked list
    next = p_solver->bodies;
    p_solver->bodies = this;

    // Mass properties and bounding radius, per shape. The caller sets `shape` before this
    // runs, so switch on it here rather than on a separate "set shape" call afterwards.
    switch (shape)
    {
        case ShapeType::Sphere:
        {
            const float r = p_size.x;
            mass = (4.0f / 3.0f) * 3.14159265358979f * r * r * r * p_density;
            const float i = 0.4f * mass * r * r; // 2/5 m r^2
            moment = float3{i, i, i};
            radius = r;
            break;
        }

        case ShapeType::Cylinder:
        {
            const float r = p_size.x;
            const float h = p_size.z;
            mass = 3.14159265358979f * r * r * h * p_density;
            const float axial = 0.5f * mass * r * r;                    // about the local Z axis
            const float radial = mass * (3.0f * r * r + h * h) / 12.0f; // about local X and Y
            moment = float3{radial, radial, axial};
            // Bounding radius of the whole cylinder: half the diagonal of its box.
            radius = std::sqrt(dot(float3(r, h * 0.5f, r), float3(r, h * 0.5f, r)));
            break;
        }

        case ShapeType::Box:
        default:
        {
            mass = p_size.x * p_size.y * p_size.z * p_density;
            moment = float3 {
                (p_size.y * p_size.y + p_size.z * p_size.z) / 12.0f * mass,
                (p_size.x * p_size.x + p_size.z * p_size.z) / 12.0f * mass,
                (p_size.x * p_size.x + p_size.y * p_size.y) / 12.0f * mass
            };
            radius = std::sqrt(dot(p_size * 0.5f, p_size * 0.5f));
            break;
        }
    }
}

Rigid::~Rigid()
{
    // Remove from linked list
    Rigid** p = &solver->bodies;
    while (*p != this)
        p = &(*p)->next;
    *p = next;
}

float3 aabbHalfExtent(const Rigid &body) noexcept
{
    float3 half;
    switch (body.shape)
    {
        case ShapeType::Sphere:
            half = float3{body.size.x, body.size.x, body.size.x};
            break;
        case ShapeType::Cylinder:
            // Local +Z is the axis: radius spans X/Y, height spans Z.
            half = float3{body.size.x, body.size.x, body.size.z * 0.5f};
            break;
        case ShapeType::Box:
        default:
            half = body.size * 0.5f;
            break;
    }

    const float3x3 rot = qmat(body.positionAng);
    return float3{
            std::abs(rot[0][0]) * half.x + std::abs(rot[1][0]) * half.y + std::abs(rot[2][0]) * half.z,
            std::abs(rot[0][1]) * half.x + std::abs(rot[1][1]) * half.y + std::abs(rot[2][1]) * half.z,
            std::abs(rot[0][2]) * half.x + std::abs(rot[1][2]) * half.y + std::abs(rot[2][2]) * half.z};
}

} // namespace avbd
