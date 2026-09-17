/*
 * Conversions between Godot's Y-up scene space and the AVBD solver's Z-up space.
 *
 * The solver core (src/avbd, ported from github.com/savant117/avbd-demo3d) works in
 * Z-up coordinates with gravity along -Z, and treats boxes as centres with full
 * widths. Godot is Y-up with gravity along -Y. The mapping between the two is a
 * +90 degree rotation about X:
 *
 *     sim = (x, -z, y)          godot = (x, z, -y)
 *
 * which maps Godot +Y onto sim +Z and Godot +Z onto sim -Y, preserving handedness.
 *
 * "Sim space" below means the coordinate system defined by the soft-body world node's
 * transform (AVBDSoftWorld3D's local space), or the server's world for the physics-server
 * path - the two agree because a server scene has no intermediate world transform.
 */

#ifndef AVBD_GODOT_CONVERT_HPP
#define AVBD_GODOT_CONVERT_HPP

#include <cmath>

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

// A direction expressed on a Godot axis, in solver space. Anything that names an *axis index*
// rather than passing a vector - per-axis joint modes, limits and springs - has to go through
// this, or it will drive the wrong axis: Godot Y is solver Z, and Godot Z is solver -Y.
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

// Property value to solver stiffness: strictly negative means "infinitely stiff"
// (a hard constraint). Godot floats and scene files cannot carry INFINITY, so the
// sentinel keeps the property serialisable. Zero disables the constraint.
inline float stiffness_from_property(double p_value) {
    return p_value < 0.0 ? INFINITY : (float)p_value;
}

} // namespace avbd_godot

#endif // AVBD_GODOT_CONVERT_HPP
