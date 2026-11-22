/* Godot Y-up <-> solver Z-up: sim=(x,-z,y)=+90 deg about X (Godot +Y->sim +Z, +Z->sim -Y);
 * solver gravity is -Z and boxes are centres with full widths. Core: savant117/avbd-demo3d (MIT). */

#ifndef AVBD_GODOT_CONVERT_HPP
#define AVBD_GODOT_CONVERT_HPP

#include <cmath>
#include <limits>

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/variant/basis.hpp>
#include <godot_cpp/variant/quaternion.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/maths.h"

namespace avbd_godot {

// --- vectors -----------------------------------------------------------------

inline avbd::float3 to_sim(const godot::Vector3 &v) {
    return avbd::float3{v.x, -v.z, v.y};
}

inline godot::Vector3 to_godot(const avbd::float3 &v) {
    return godot::Vector3(v.x, v.z, -v.y);
}

// Axis indices (per-axis joint modes, limits, springs) must go through here or they drive
// the wrong axis: Godot Y is solver Z, Godot Z is solver -Y.
inline avbd::float3 axis_to_sim(int p_godot_axis) {
    return to_sim(p_godot_axis == 0 ? godot::Vector3(1, 0, 0)
                                    : (p_godot_axis == 1 ? godot::Vector3(0, 1, 0)
                                                         : godot::Vector3(0, 0, 1)));
}

// Which solver axis a Godot axis becomes, and whether its positive direction is preserved.
// Godot X -> solver X (+), Godot Y -> solver Z (+), Godot Z -> solver -Y (-).
struct AxisMapping {
    int solver;
    float sign;
};

inline constexpr AxisMapping kAxisMapping[3] = {{0, 1.0f}, {2, 1.0f}, {1, -1.0f}};

inline avbd::float3 to_sim_extents(const godot::Vector3 &v) {
    return avbd::float3{std::fabs(v.x), std::fabs(v.z), std::fabs(v.y)};
}

// --- rotations ---------------------------------------------------------------

inline avbd::quat to_sim(const godot::Quaternion &q) {
    // linalg quaternions are {x, y, z, w}; the component mapping is
    // solver (x, y, z) = Godot (x, -z, y).
    return avbd::quat{q.x, -q.z, q.y, q.w};
}

inline godot::Quaternion to_godot(const avbd::quat &q) {
    return godot::Quaternion(q.x, q.z, -q.y, q.w);
}

// A basis rotated into solver space, normalised (AVBD quaternions integrate by
// addition and are only renormalised, so they can drift by a few ulps).
inline godot::Basis to_godot_basis(const avbd::quat &q) {
    return godot::Basis(to_godot(q).normalized());
}

// --- transforms --------------------------------------------------------------

// Transform of `node` expressed in the simulation space of `world`.
inline godot::Transform3D node_in_sim(const godot::Node3D *world, const godot::Node3D *node) {
    return world->get_global_transform().affine_inverse() * node->get_global_transform();
}

inline godot::Transform3D sim_transform(const avbd::float3 &position, const avbd::quat &rotation) {
    return godot::Transform3D(to_godot_basis(rotation), to_godot(position));
}

// World direction (Godot axes, no translation) in solver space.
inline avbd::float3 direction_to_sim(const godot::Transform3D &world_global, const godot::Vector3 &v) {
    return to_sim(world_global.affine_inverse().basis.xform(v));
}

inline godot::Vector3 direction_to_godot(const godot::Transform3D &world_global, const avbd::float3 &v) {
    return world_global.basis.xform(to_godot(v));
}

// A point in Godot world space expressed in solver space.
inline avbd::float3 point_to_sim(const godot::Transform3D &world_global, const godot::Vector3 &p) {
    return to_sim(world_global.affine_inverse().xform(p));
}

// Stiffness: strictly negative is the "infinitely stiff" sentinel (hard constraint), since
// Godot floats/scenes cannot carry infinity; zero disables the constraint.
inline float stiffness_from_property(double p_value) {
    return p_value < 0.0 ? std::numeric_limits<float>::infinity() : (float)p_value;
}

} // namespace avbd_godot

#endif // AVBD_GODOT_CONVERT_HPP
