/*
 * AVBDSpring3D - linear spring between two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * Both endpoints must be bodies: the solver core reads both bodies when it builds
 * the spring (to derive its rest length). Anchor one end to a static body to hang
 * a spring from the world.
 */

#ifndef AVBD_SPRING3D_HPP
#define AVBD_SPRING3D_HPP

#include <godot_cpp/variant/vector3.hpp>

#include "nodes/avbd_constraint3d.hpp"

namespace godot {

class AVBDSpring3D : public AVBDConstraint3D {
    GDCLASS(AVBDSpring3D, AVBDConstraint3D)

    Vector3 anchor_a;
    Vector3 anchor_b;
    double stiffness = 1000.0;
    double rest_length = -1.0;

protected:
    static void _bind_methods();

public:
    AVBDSpring3D() = default;
    ~AVBDSpring3D() override = default;

    void set_anchor_a(const Vector3 &p_anchor);
    Vector3 get_anchor_a() const;
    void set_anchor_b(const Vector3 &p_anchor);
    Vector3 get_anchor_b() const;

    // Spring constant in N/m.
    void set_stiffness(double p_stiffness);
    double get_stiffness() const;

    // Rest length in metres. Negative means "derive it from the bodies' pose when
    // the spring is created".
    void set_rest_length(double p_length);
    double get_rest_length() const;

    avbd::Force *create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies,
            const Transform3D &p_world_global) const override;
};

} // namespace godot

#endif // AVBD_SPRING3D_HPP
