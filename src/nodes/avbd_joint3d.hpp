/*
 * AVBDJoint3D - ball-and-socket constraint (optionally angular and breakable)
 * between two bodies, or between one body and a fixed point in space.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */

#ifndef AVBD_JOINT3D_HPP
#define AVBD_JOINT3D_HPP

#include <godot_cpp/variant/vector3.hpp>

#include "nodes/avbd_constraint3d.hpp"

namespace godot {

class AVBDJoint3D : public AVBDConstraint3D {
    GDCLASS(AVBDJoint3D, AVBDConstraint3D)

    Vector3 anchor_a;
    Vector3 anchor_b;
    double linear_stiffness = -1.0;
    double angular_stiffness = 0.0;
    double fracture_force = -1.0;
    bool broken = false;

protected:
    static void _bind_methods();

public:
    AVBDJoint3D() = default;
    ~AVBDJoint3D() override = default;

    // Anchor on body A. A body-local offset, or a world-space position when
    // node_a is empty (then the joint pins body B to that point in space).
    void set_anchor_a(const Vector3 &p_anchor);
    Vector3 get_anchor_a() const;

    // Anchor on body B, as a body-local offset.
    void set_anchor_b(const Vector3 &p_anchor);
    Vector3 get_anchor_b() const;

    // Linear stiffness in N/m. Negative means an infinitely stiff (hard)
    // constraint, 0 disables the linear constraint, positive values give a soft
    // joint whose penalty is capped at this stiffness.
    void set_linear_stiffness(double p_stiffness);
    double get_linear_stiffness() const;

    // Angular stiffness. Negative is rigid (the joint also locks rotation), 0
    // leaves rotation free, positive values give a soft angular spring.
    void set_angular_stiffness(double p_stiffness);
    double get_angular_stiffness() const;

    // Angular force above which the joint breaks and stops constraining the
    // bodies. Negative means unbreakable.
    void set_fracture_force(double p_force);
    double get_fracture_force() const;

    // Whether the joint has broken. Breaking is permanent for this node: the solver
    // deletes a broken joint, so the flag is re-applied whenever the world rebuilds
    // the simulation. Setting it back to false mends the joint (the next rebuild
    // creates it unbroken).
    bool is_broken() const;
    void set_broken(bool p_broken);
    void break_joint();

    avbd::Force *create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies,
            const Transform3D &p_world_global) const override;

    // A broken joint gets no solver force: the solver has already removed it, and
    // leaving it out lets the two bodies collide again.
    bool _is_inactive() const override { return broken; }

protected:
    void _on_force_removed() override;
};

} // namespace godot

#endif // AVBD_JOINT3D_HPP
