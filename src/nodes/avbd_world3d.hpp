/*
 * AVBDWorld3D - root of an AVBD simulation.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * Owns the solver and drives it once per physics tick. Every AVBDRigidBody3D,
 * AVBDJoint3D, AVBDSpring3D, AVBDIgnoreCollision3D and AVBDSoftBody3D in its
 * subtree is collected automatically; adding or removing one rebuilds the solver
 * while keeping the current pose of the bodies that stay. This node's transform
 * defines the simulation's coordinate system.
 *
 * The simulation is a fixed-step, single-threaded, deterministic Gauss-Seidel
 * solver: identical node graphs stepped with the same tick rate produce identical
 * results.
 */

#ifndef AVBD_WORLD3D_HPP
#define AVBD_WORLD3D_HPP

#include <vector>

#include <godot_cpp/classes/node3d.hpp>
#include <godot_cpp/variant/dictionary.hpp>
#include <godot_cpp/variant/transform3d.hpp>
#include <godot_cpp/variant/vector3.hpp>

#include "avbd/solver.h"

namespace godot {

class AVBDRigidBody3D;
class AVBDConstraint3D;
class AVBDSoftBody3D;

class AVBDWorld3D : public Node3D {
    GDCLASS(AVBDWorld3D, Node3D)

    double gravity = 9.8;
    int iterations = 10;
    double alpha = 0.99;
    double beta_linear = 10000.0;
    double beta_angular = 100.0;
    double gamma = 0.999;
    int substeps = 1;
    bool paused = false;
    int threads = 0;

    struct BodyEntry {
        AVBDRigidBody3D *node = nullptr;
        avbd::Rigid *rigid = nullptr;
        // Last known solver state, so rebuilding keeps the current pose.
        avbd::float3 positionLin;
        avbd::quat positionAng;
        avbd::float3 velocityLin;
        avbd::float3 velocityAng;
        // Node transform last written back, used to detect user edits of static bodies.
        Transform3D pushed_transform;
        bool pushed = false;
    };

    struct SoftBodyEntry {
        AVBDSoftBody3D *node = nullptr;
    };

    avbd::Solver solver;
    std::vector<BodyEntry> bodies;
    std::vector<AVBDConstraint3D *> constraints;
    std::vector<SoftBodyEntry> soft_bodies;
    uint64_t last_step_usec = 0;
    uint64_t step_count = 0;

    // Scratch buffers for the per-tick scene scan (kept to avoid per-frame allocation).
    std::vector<AVBDRigidBody3D *> scan_bodies;
    std::vector<AVBDConstraint3D *> scan_constraints;
    std::vector<AVBDSoftBody3D *> scan_soft_bodies;

protected:
    static void _bind_methods();
    void _notification(int p_what);

public:
    AVBDWorld3D();
    ~AVBDWorld3D() override;

    // Gravity in m/s^2, applied along -Y (Godot's down direction).
    void set_gravity(double p_gravity);
    double get_gravity() const;

    // Solver iterations per step. More iterations mean stiffer contacts and
    // joints at a linear cost; 10 is the paper's default.
    void set_iterations(int p_iterations);
    int get_iterations() const;

    // Constraint stabilisation in (0, 1]: higher is smoother and less energetic.
    void set_alpha(double p_alpha);
    double get_alpha() const;

    // Penalty ramping parameters for linear (N/m) and angular (N.m) constraints.
    // Larger values make contacts and joints stiffer; the useful range depends on
    // the mass and length scales of the scene.
    void set_beta_linear(double p_beta);
    double get_beta_linear() const;
    void set_beta_angular(double p_beta);
    double get_beta_angular() const;

    // Warm-start decay for the dual variables, < 1.
    void set_gamma(double p_gamma);
    double get_gamma() const;

    // Solver steps per physics tick (each with dt / substeps).
    void set_substeps(int p_substeps);
    int get_substeps() const;

    void set_paused(bool p_paused);
    bool is_paused() const;

    // Worker threads for the solver's per-body phases. 0 (the default) uses one per hardware
    // thread, 1 runs everything on the calling thread. The result does not depend on this
    // value: bodies are split into groups that share no constraint, so splitting them across
    // more threads cannot change any number. Reported by get_thread_count().
    void set_threads(int p_threads);
    int get_threads() const;
    // Worker threads actually in use.
    int get_thread_count() const;

    // --- introspection ------------------------------------------------------
    int get_body_count() const;
    int get_force_count() const;
    int get_contact_count() const; // contact manifolds: broad-phase pairs
    int get_contact_point_count() const;
    uint64_t get_step_time_usec() const;
    // Number of solver steps performed since the node entered the tree. Useful for
    // offline/batch simulation and for tests that need an exact tick count.
    uint64_t get_step_count() const;

    // Re-create the solver from the scene graph, keeping the pose of bodies that
    // are already simulated. Called automatically when the subtree changes.
    void rebuild();

    // The result of a ray query, in Godot's world space. `node` is null for a miss; for a
    // lattice soft body `cell` names the cell that was hit, and is -1 otherwise.
    struct Hit {
        Node3D *node = nullptr;
        Vector3 position;
        double distance = 0.0;
        int cell = -1;
    };

    // Trace a ray against the simulated bodies. `direction` need not be normalised, and
    // static bodies are hit as well.
    [[nodiscard]] Hit trace(const Vector3 &p_origin, const Vector3 &p_direction, double p_max_distance = 10000.0);

    // Same query, as the Dictionary a script sees. Keys: body, position, distance, and cell
    // for lattice soft bodies. Empty on a miss.
    Dictionary raycast(const Vector3 &p_origin, const Vector3 &p_direction, double p_max_distance = 10000.0);

    // --- simulation interface, used by the nodes above -----------------------
    avbd::Rigid *rigid_for(const AVBDRigidBody3D *p_body) const;

    void _physics_process(double p_delta) override;
    void _ready() override;

private:
    void _scan();
    bool _scan_changed();
    void _collect(Node *p_node);
    void _rebuild();
    void _release();
    void _apply_params();
    void _step(double p_delta);
    void _resolve_removed_forces();
    [[nodiscard]] Hit identify(const avbd::Rigid *p_rigid) const;
    void _sync_in_static_bodies();
    void _sync_out();
};

} // namespace godot

#endif // AVBD_WORLD3D_HPP
