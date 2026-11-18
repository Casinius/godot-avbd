/*
 * modern_convert.hpp - Godot <-> avbd-mo type conversions for the server path.
 *
 * Both sides are Y-up, right-handed, with the same quaternion convention, so this is a
 * pure component copy. That is the whole point of the modern server path: unlike the
 * legacy Z-up core used by the soft-body nodes, no axis permutation belongs here.
 */

#ifndef AVBD_MODERN_CONVERT_HPP
#define AVBD_MODERN_CONVERT_HPP

#include <godot_cpp/variant/quaternion.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "phys_core/state.hpp"

namespace godot {

inline Vec3 to_modern(const Vector3 &p_v) noexcept {
    return Vec3{static_cast<Scalar_t>(p_v.x), static_cast<Scalar_t>(p_v.y), static_cast<Scalar_t>(p_v.z)};
}

inline Vector3 to_godot(const Vec3 &p_v) noexcept {
    return Vector3(static_cast<double>(p_v.x), static_cast<double>(p_v.y), static_cast<double>(p_v.z));
}

inline Quat to_modern(const Quaternion &p_q) noexcept {
    return Quat{static_cast<Scalar_t>(p_q.w), static_cast<Scalar_t>(p_q.x),
                static_cast<Scalar_t>(p_q.y), static_cast<Scalar_t>(p_q.z)};
}

inline Quaternion to_godot(const Quat &p_q) noexcept {
    return Quaternion(static_cast<double>(p_q.x), static_cast<double>(p_q.y),
                      static_cast<double>(p_q.z), static_cast<double>(p_q.w));
}

} // namespace godot

#endif // AVBD_MODERN_CONVERT_HPP
