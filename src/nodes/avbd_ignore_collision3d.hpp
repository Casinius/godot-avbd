/*
 * AVBDIgnoreCollision3D - suppresses contacts between two bodies.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * Bodies joined by any constraint already skip each other when the solver looks
 * for contacts, so this exists for pairs that would otherwise collide spuriously:
 * diagonal neighbours in a lattice, segments of a rope that pass close by, wheels
 * versus the chassis they are mounted on.
 */

#ifndef AVBD_IGNORE_COLLISION3D_HPP
#define AVBD_IGNORE_COLLISION3D_HPP

#include "nodes/avbd_constraint3d.hpp"

namespace godot {

class AVBDIgnoreCollision3D : public AVBDConstraint3D {
    GDCLASS(AVBDIgnoreCollision3D, AVBDConstraint3D)

protected:
    static void _bind_methods();

public:
    AVBDIgnoreCollision3D() = default;
    ~AVBDIgnoreCollision3D() override = default;

    avbd::Force *create_force(avbd::Solver &p_solver, const AVBDBodyMap &p_bodies,
            const Transform3D &p_world_global) const override;
};

} // namespace godot

#endif // AVBD_IGNORE_COLLISION3D_HPP
