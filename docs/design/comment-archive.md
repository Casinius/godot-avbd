# Comment archive

Full text of every comment block longer than two lines, as it stood before the compression pass. Kept for the rationale (measured results, known limitations, why an expression must not be reordered). Source files now carry at most two lines each.

### `src/avbd/broadphase_pair.hpp:1`

```cpp
/*
 * Shared broad-phase pair encoding.
 *
 * The Morton LBVH emits this packed candidate key, and the sort/dedupe/filter downstream
 * consumes it. The key is (smaller body index << 32) | larger body index: lexicographic
 * order over std::pair is plain integer order, which is what lets the pair list use a radix
 * sort.
 */
```

### `src/avbd/bvh_broadphase.cpp:1`

```cpp
/*
 * Morton LBVH broad phase. See bvh_broadphase.hpp for what the port changed and why.
 */
```

### `src/avbd/bvh_broadphase.cpp:117`

```cpp
    // Leaves first (in sorted order), so the merge below is a simple range walk. A range
    // of leaves [begin, end) is replaced by the internal node that spans it; the stack of
    // open ranges is rebuilt without recursion (a 4k-body scene would otherwise recurse
    // log2(n) deep per level with a heap allocation at every step).
```

### `src/avbd/bvh_broadphase.cpp:142`

```cpp
    // Every node was created in a build order that puts children before parents inside a
    // pair but not globally, so the walk is a post-order over the stored child indices.
    // Iterative with an explicit stack: same reason as the build.
```

### `src/avbd/bvh_broadphase.hpp:1`

```cpp
/*
 * Morton LBVH broad phase, ported from `main:src/avbd/bvh/`. It is the solver's only
 * broad phase; the hierarchical hashed grid it used to share the pipeline with has been
 * removed.
 *
 * The port keeps the algorithm and changes three things, all of which were needed for the
 * comparison that selected it:
 *
 *   1. Eigen -> linalg (the `main` tree is the pre-port branch).
 *   2. Leaves are sized from the solver's **swept world AABBs** (`half`/`centre`) instead
 *      of `body->radius` bounding spheres. A sphere-vs-AABB difference is a bigger effect
 *      than the tree-vs-grid difference, so sharing the input keeps the broad phase honest
 *      about the bodies it is given.
 *   3. `NodeStorage` is folded into one class. The original spread the arrays over a
 *      separate translation unit with accessors per node; that costs nothing measurable but
 *      makes the file count four instead of two.
 *
 * Output contract: packed `packPair(i, j)` candidates, unsorted and with each overlapping
 * pair emitted once. The caller sorts, dedupes and applies the exact AABB/layer/constraint
 * filter.
 *
 * Determinism: the build sorts by (morton code, body index), the traversal pushes children
 * in a fixed order, and the pair list is integer-sorted by the caller. The build and walk
 * are serial - the `pool` parameter is reserved for future parallel work.
 */
```

### `src/avbd/bvh_broadphase.hpp:56`

```cpp
    // Refits allowed before the topology is rebuilt from fresh Morton codes. A refit is
    // cheap (one post-order walk) but the topology ages: bodies that move far from where
    // their code said they were leave the tree unbalanced and the overlap tests worse.
```

### `src/avbd/collide.cpp:1`

```cpp
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
```

### `src/avbd/collide.cpp:235`

```cpp
// SAT with the frame table. The caller precomputes m[i][j] = axisA_i·axisB_j and
// dA/dB = delta in each box's frame once per pair; every axis projection then reduces
// to table reads:
//   face axis A_i: |delta·A_i| = |dA[i]| and the projection radius of A is exactly
//     halfA[i] (its own face normal) - no sqrt, no per-axis absDot storm;
//   edge axis A_i×B_j: the cross product is perpendicular to both parents, so each
//     box's radius reads the two off-axis entries of m (perpendicularity identity).
// The generic per-axis form (normalise, absDot both boxes' three axes) did ~30
// multiplies plus a sqrt per axis, 15 axes per pair.
```

### `src/avbd/collide.cpp:297`

```cpp
    // (A_i x B_j) is perpendicular to A_i and B_j, so each box's support along it is
    // the support along the two adjacent axes, read off the frame table:
    //   |A_k . (A_i x B_j)| = |m[(i+2)%3][j]| for k=(i+1)%3, |m[(i+1)%3][j]| for k=(i+2)%3
```

### `src/avbd/collide.cpp:461`

```cpp
        // The clip window caught nothing: the faces only meet at a corner. The support
        // points are the deepest points along the normal, so their normal separation is
        // the contact error - but their *tangential* offset is arbitrary (it is the
        // distance between whichever corners the support direction happens to pick, up to
        // a full box width). That offset becomes a tangential constraint error, ramps the
        // friction penalty against a phantom, and after the next rebase rotates into the
        // normal: box-sized penetration and an explosive impulse. Keep the normal
        // separation, put both anchors at A's support point.
```

### `src/avbd/collide.cpp:507`

```cpp
// -----------------------------------------------------------------------------
// Shapes
//
// Box-box is the reference implementation's SAT + clipping, above. The other pairs are built on
// closest points, which is exact for a sphere against anything and for the curved side of a
// cylinder; the flat caps of a cylinder are handled by testing the rim and cap samples against
// the other shape. Sampling a curved surface is the usual way to get a manifold (rather than a
// single point) out of it, and the samples are placed deterministically, so contacts are stable
// from frame to frame and warm-starting works.
// -----------------------------------------------------------------------------
```

### `src/avbd/collide.cpp:569`

```cpp
// Speculative variant of pointInBox: a point outside but within `margin` of the box reports a
// gap - r_depth negative, r_push pointing out of the box toward the point, so the existing
// anchor math `other = sample + push * depth` lands on the box surface for both signs. The
// closest point is the componentwise clamp of the local coordinates, which is exact for box
// edges. Returns false when the point is inside (use pointInBox) or farther than margin out.
```

### `src/avbd/collide.cpp:677`

```cpp
// Speculative variant of pointInCylinder: a point outside but within `margin` of the cylinder
// reports a gap - r_depth negative, r_push pointing out of the cylinder toward the point, so
// `other = sample + push * depth` lands on the cylinder surface for both signs. Returns false
// when the point is inside (use pointInCylinder) or farther than margin out.
```

### `src/avbd/collide.cpp:699`

```cpp
// The validity test both cylinder paths run on a sample: inside the shape, or within the
// speculative margin of it. The order matters: the exact inside query runs first and only its
// failure lets the near query run, exactly as the open-coded `!in && !(margin > 0 && near)` did.
```

### `src/avbd/collide.cpp:722`

```cpp
    // Stable identity of the sample that produced this contact: a path tag, which shape
    // owns the sample, and its local index. Warm-start matching keys on this, so it must
    // not depend on the depth ranking (a rank-keyed contact resets its state whenever the
    // ranking moves, which for sampled cylinders is almost every tick).
```

### `src/avbd/collide.cpp:729`

```cpp
// Emit the deepest MAX_CONTACTS candidates. Ranking by depth rather than taking them in sampling
// order is what keeps a contact set sensible when there are more samples than the budget: the
// shallow ones are the ones to drop.
//
// The ranking is an insertion sort into a stack array: stable (equal depths keep sampling
// order) and allocation-free, unlike the vector + stable_sort copy this replaces.
```

### `src/avbd/collide.cpp:765`

```cpp
// Samples a cylinder's surface: the two cap centres, the two rim circles, and rings partway
// along the side (which is what gives a lying cylinder a line of contacts rather than just its
// two ends).
//
// The angles are visited in an interleaved order rather than 0, 1, 2, ... because the contact
// budget is smaller than the sample count: adding contacts stops at MAX_CONTACTS, and taking the
// first N of an ordered ring would cover one arc and leave the other side unsupported - which
// tips the body over. With a stride coprime to the sample count, any prefix is spread around the
// circle.
// The rim angles are the same on every call, so the sine and cosine table is built once.
// The entries are produced by the same `std::cos` / `std::sin` of the same float angle the
// per-sample code used, which keeps the sampled positions bit-identical.
```

### `src/avbd/collide.cpp:930`

```cpp
// Cylinder against box, in either order. Both directions are tested - the cylinder's surface
// samples against the box, and the box's corners against the cylinder - so neither shape can slip
// through the other however they are arranged.
```

### `src/avbd/collide.cpp:1070`

```cpp
        // The sampled cylinder paths accept gap contacts through the same signed point
        // queries: a sample within the margin of the other shape reports a negative depth
        // (a gap), and the solver's speculative machinery turns that into a contact that
        // only pushes once the pair would close the gap within the step.
```

### `src/avbd/fixed_stack.hpp:1`

```cpp
/*
 * FixedStack - a bounded LIFO for hot traversals that must not allocate.
 *
 * The BVH descents (broad-phase self-query, pick, AABB queries) walk trees of known
 * maximum depth, where a raw `int stack[64]` plus a hand-kept index works until the
 * day the depth bound changes and the overflow silently corrupts the caller's frame.
 * This is that array with the two failure modes turned into assertions (compiled out
 * under NDEBUG, so the release path is the same instructions as the raw array) and
 * the index bookkeeping in one place.
 */
```

### `src/avbd/force.cpp:1`

```cpp
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
```

### `src/avbd/forces.hpp:1`

```cpp
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
```

### `src/avbd/forces.hpp:76`

```cpp
    // Break the joint immediately, discarding its dual state. Used when a host
    // application wants the fracture decision to survive being re-created (the
    // solver has no way to remember a broken joint once it is deleted).
```

### `src/avbd/forces.hpp:143`

```cpp
// A general 6-DOF joint: three linear and three angular degrees of freedom, each Free, Locked
// or Limited, each with an optional spring. The degrees of freedom are the axes of body A's
// local frame (world axes when bodyA is null).
//
// Angular degrees of freedom are measured from the pose the joint was created in. Before the
// locked axes are evaluated, the twist about every non-locked axis is removed: without that, a
// wheel spinning about its axle would drag the locked axes with it, because a rotation vector
// wraps every half turn and so cannot be read directly off the deviation. Removing one twist is
// exact; with several non-locked axes it is applied axis by axis, which is the usual approximation.
```

### `src/avbd/forces.hpp:175`

```cpp
    // Current displacement of one degree of freedom, measured in the joint frame: metres for a
    // linear axis, radians for an angular one, both relative to the pose the joint was created
    // in. Reporting only - the solver never reads these back.
```

### `src/avbd/forces.hpp:194`

```cpp
    // Stamp one degree of freedom into the body's system: its constraint, and its spring if it
    // has one.
    //
    // A constraint is stamped through *both* halves of the system: `p_jLin` is how it responds
    // to the body moving and `p_jAng` to the body turning. They are the two Jacobian rows of one
    // constraint, so a linear axis anchored off the centre of mass needs both - it moves the
    // anchor point as much by rotating the body as by translating it.
```

### `src/avbd/generic_joint.cpp:1`

```cpp
/*
 * General 6-DOF joint.
 *
 * Part of godot-avbd; the augmented-Lagrangian structure follows the reference
 * implementation (github.com/savant117/avbd-demo3d, MIT), as the rest of src/avbd does.
 *
 * Each of the six degrees of freedom is Free, Locked or Limited, and may carry a spring. A
 * degree of freedom is measured in the joint frame, which is body A's local frame anchored at
 * rA - or the world when bodyA is null, which is how the joint pins a body to a fixed pose.
 * Measurements are in physical units: metres for a linear axis, radians for an angular one,
 * and radians for the limits, so they mean the same thing whatever size the bodies are.
 *
 * Two things here are not a direct copy of anything in the reference implementation.
 *
 * The angular side has to survive a spinning axle. A locked axis cannot read its component of
 * the rotation vector between the two bodies: a wheel turning about its axle sweeps those
 * components through the full circle every half turn, so a hinge built that way tears itself
 * apart after a few revolutions. Instead the deviation from the joint's rest orientation is
 * split into the twist about each non-locked axis and the remaining swing, and only the swing
 * is handed to the locked axes.
 *
 * And a constraint anchored away from a body's centre acts through a moment arm, so it is
 * stamped through both halves of the system at once - see stampAxis.
 */
```

### `src/avbd/generic_joint.cpp:71`

```cpp
// Remove the twist about axis `p_axis` from `p_dev`, writing its angle to `r_angle`.
//
// A quaternion's twist about n is (n * (v.n), w) normalised, so the angle is 2*atan2(v.n, w).
// That stays continuous through a full turn, which is what a spinning axle needs. Near a half
// turn the twist quaternion degenerates (both the projected vector part and w vanish), and
// there is then no twist to remove.
```

### `src/avbd/generic_joint.cpp:98`

```cpp
    // Deviation from the rest configuration, in the joint frame (body A's).
    //
    // Both quantities below are B relative to A as seen *by A*: the current one, and the one
    // captured when the joint was created. The deviation is then the rest pose taken away from
    // the current one, which is identity exactly when the bodies are in their rest arrangement.
    //
    // The order matters and is easy to get wrong in a way nothing notices: putting the rest term
    // anywhere else happens to be harmless while the rest pose is the identity (every joint built
    // from upright bodies), and silently wrong as soon as a body starts out rotated - a wheel
    // turned onto its axle, say. Then the locked axes see an error proportional to the rest
    // rotation, which grows as the body spins and tears the joint apart.
    // B relative to A now, and the same quantity at rest. Taking the rest arrangement away from
    // the current one gives identity exactly when the bodies are back where they started, and its
    // axis-angle in A's frame is what each degree of freedom measures.
```

### `src/avbd/generic_joint.cpp:116`

```cpp
    // The twist about every non-locked axis belongs to that axis - a spring or a limit reads
    // its angle - and must not be visible to the locked ones. Removing one twist is exact; with
    // several non-locked axes it is applied one axis at a time, the usual approximation.
```

### `src/avbd/generic_joint.cpp:176`

```cpp
// The force one degree of freedom applies this pass, and the constraint value behind it.
//
// Both passes must agree on this value: the dual pass stores the force it produces as the next
// multiplier and grows the penalty from the same error. Computing it in one place is what keeps
// the augmented Lagrangian consistent - a lock stabilised in only one of the two would have the
// primal fighting the full error while the dual believed it was correcting a fraction of it,
// which converges to the wrong answer rather than merely converging slowly.
```

### `src/avbd/generic_joint.cpp:189`

```cpp
    // A lock is stabilised against where the step started, so its error is worked out over
    // several steps instead of being fought in one. A limit is not: biasing it would let it
    // settle short of its bound.
```

### `src/avbd/generic_joint.cpp:208`

```cpp
    // Constraint term of the augmented Lagrangian, stamped through both Jacobian rows. The
    // cross block is what couples the body's translation and rotation, and it is not optional:
    // an anchor away from the centre of mass moves as much by turning the body as by moving it,
    // so dropping that term constrains the wrong point.
```

### `src/avbd/generic_joint.cpp:237`

```cpp
// How stiff a degree of freedom may become before the solver can no longer resolve it.
//
// The augmented Lagrangian grows a constraint's penalty while it is violated, which is what makes
// it converge to zero error. Left unbounded that is a trap: the constraint's own frequency is
// sqrt(penalty / effective inertia), and once it approaches the step rate the iteration stops
// converging, so the error grows, so the penalty grows further. A degree of freedom that is
// *already* satisfied has nothing to gain from a large penalty and everything to lose, and a
// body spinning fast is exactly where it shows up first.
//
// The cap is the penalty at which the constraint's frequency is still a comfortable fraction of
// the step rate, using the smaller of the two bodies' effective inertia about the axis - the
// lighter side is the one that limits the timestep.
```

### `src/avbd/generic_joint.cpp:360`

```cpp
                // Complementary slackness: a limit that is not being pushed against applies no
                // force, so its multiplier is zero. Keeping the previous one would re-apply it
                // the moment the axis came back to the bound, which reads as a violent kick -
                // and the stiffer the axis, the harder the kick, so a stiffer limit would
                // behave *worse* than a soft one.
```

### `src/avbd/job_pool.hpp:1`

```cpp
/*
 * JobPool - the parallel machinery the solver runs its phases on.
 *
 * Thin adapter over BS::thread_pool v5.1.0 (https://github.com/bshoshany/thread-pool,
 * MIT, header-only, C++17). The library owns the worker queue and lifetime; what is left
 * to us is the part no library can know: how many independent pieces of work there are,
 * which of them may overlap, and where the phase boundaries are.
 *
 * TASK CLASSES (the solver sequences them directly):
 * - Independent work: per-item functions with no internal synchronization. Any worker may
 *   run any item, and several independent tasks may be in flight at once. Joined by a
 *   per-dispatch latch, never by a global queue drain.
 * - Region work: a body that runs on every participant and fences its own phases with
 *   LoopSync barriers (the Gauss-Seidel primal). Regions hold a budget of participant
 *   slots for their lifetime, so two solvers can run their regions concurrently without
 *   either silently falling back to serial. Regions must not nest.
 * - Serial work: runs on the calling thread; used for mutations (manifold linking,
 *   deletion, colouring) that must not race with anything else.
 *
 * WHY THE GLOBAL MUTEX IS GONE:
 * The previous design serialized all solvers with a try-locked mutex around the persistent
 * loop and joined with `sharedPool().wait()`. That `wait()` drains the whole queue, so it
 * cannot coexist with region participants blocked on a barrier: it would wait for them
 * forever. Every join here is a per-dispatch `std::latch`, which also lets independent
 * tasks from another solver run while a region is parked at a barrier.
 *
 * BUDGET:
 * - One process-wide worker pool, lazily created on first parallel use; serial solvers
 *   (threads=1) never create it.
 * - Each JobPool's `threadCount()` is its requested budget, clamped to the machine.
 * - A region needs `threadCount()-1` workers (the calling thread is participant 0) and
 *   reserves them from a counting semaphore. Two regions at full width therefore run one
 *   after the other (by budget, not by silent serialization); setting `threads=4` on two
 *   solvers on an 8-thread machine lets both run at once.
 *
 * MEMORY MODEL:
 * - `std::latch` release/acquire covers the dispatch join; `std::barrier` covers every
 *   region phase boundary. Work items are disjoint by contract, so no floating-point
 *   reduction happens across threads and the result is independent of the split.
 */
```

### `src/avbd/job_pool.hpp:63`

```cpp
// The phase-boundary synchronisers handed to a region body.
//
// The loop body is written once against `auto &sync` and must treat every
// `arriveAndWait()` as a hard phase fence: all participants arrive, then all proceed. Both
// types provide `forItems`, which spreads one phase's independent items over the
// participants (or runs them serially) - the barrier call itself stays in the solver's
// loop so the Gauss-Seidel structure reads plainly.
```

### `src/avbd/job_pool.hpp:73`

```cpp
    // Rearm the work cursor for the next phase. Runs exactly once per completed phase,
    // after every participant has arrived (so every pull of the finished phase is done)
    // and before any participant is released (so every pull of the next phase sees 0).
```

### `src/avbd/job_pool.hpp:95`

```cpp
    // Run `fn(i)` for every i in [0, count), self-balanced: workers grab chunks from the
    // shared cursor until this phase's range is drained. Items must be independent - the
    // split varies with timing by design. The cursor is rearmed by the barrier's
    // completion, so consecutive phases must be separated by `arriveAndWait()`, which is
    // the loop contract anyway.
```

### `src/avbd/job_pool.hpp:144`

```cpp
// All-or-nothing worker reservation for regions. A region asks for its whole helper
// width at once, before dispatching anything: holding a partial reservation while waiting
// for the rest is exactly the deadlock two concurrent regions can create, and a plain
// counting semaphore does not offer a bulk acquire. `notify_all` with a size predicate is
// not FIFO, which is fine here: regions are short and at most a handful exist. If a large
// region could be starved by a stream of small ones, add a ticket queue; the solver does
// not produce that pattern.
```

### `src/avbd/job_pool.hpp:156`

```cpp
    // The common case (slots are free) is a single atomic CAS, no lock. Only a region that
    // has to wait touches the mutex/condition variable, and the predicate is re-checked
    // under that mutex so a release in the check-to-wait window cannot be lost.
```

### `src/avbd/job_pool.hpp:197`

```cpp
// The process-wide pool and the region-width budget. Materialised lazily so a serial
// solver never spawns threads; static lifetime is safe because no Solver touches it
// during teardown.
```

### `src/avbd/job_pool.hpp:212`

```cpp
    // Reserve `count` region workers. Blocking: a region that cannot get its full width
    // waits for another region to finish instead of silently running serially. A region
    // never holds slots while waiting for another region, so this cannot deadlock.
```

### `src/avbd/job_pool.hpp:259`

```cpp
    // Run `fn(i)` for every i in [0, count). Independent items: the split only affects which
    // worker computes an item, never the result. The calling thread runs one chunk and then
    // joins the others on a latch.
    //
    // `minItems` is the smallest batch worth dispatching: one dispatch costs tens of
    // microseconds (waking workers and synchronising), which would dominate a scene
    // with a handful of bodies.
```

### `src/avbd/job_pool.hpp:277`

```cpp
        // Target ~4 items per chunk: enough chop points for an uneven per-body cost, few
        // enough that the latch and task bookkeeping stay in the noise. Clamp to
        // [workers, 4*workers] so a batch near the break-even size does not turn into one
        // task per item while a large batch still splits evenly.
```

### `src/avbd/job_pool.hpp:290`

```cpp
    // Run `fn(chunk, first, last)` over contiguous chunks of [0, count). Unlike `forCount`
    // this hands the chunk index to the callback, which is what lets a phase give each
    // participant its own output buffer (the broad phase's pair emission). Chunk boundaries
    // depend on the thread count, so the per-chunk order is only fixed for a given
    // `threads` value: only use this for work whose final order is imposed afterwards
    // (the pair list is integer-sorted), never for reduction order.
```

### `src/avbd/job_pool.hpp:307`

```cpp
    // Run a region body on `threadCount()` participants and return when every participant
    // has left the region. `body(worker, sync)` must fence its phases with
    // `sync.arriveAndWait()`. Regions reserve workers for their whole lifetime, so two
    // solvers can run concurrently (up to the machine's width) without the old global
    // mutex, and never silently degrade to serial.
```

### `src/avbd/job_pool.hpp:338`

```cpp
    // The small-scene escape hatch: run the region body on the calling thread with the
    // serial synchroniser. Kept separate from `runRegion` so the flow can decide per step
    // whether the work justifies waking workers at all.
```

### `src/avbd/joint.cpp:1`

```cpp
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
```

### `src/avbd/joint.cpp:149`

```cpp
        // Choose jacobian depending on input body. The reference implementation's
        // body-B half is diag(0, -1, -1) (the x row is carried by A's identity half);
        // changing it to diag(1, -1, -1) doubles the x correction and tears stiff
        // joint lattices apart, so keep the reference matrix.
        //
        // Known limitation: that split assumes A exists. With bodyA == nullptr (a joint to
        // the world) nothing carries the x row, so the x component of the angular constraint
        // produces no force while its multiplier keeps ramping. Use GenericJoint for
        // world-anchored angular constraints - it reads the world frame as body A's frame
        // and constrains all three axes. This class is kept for the core tests' joints.
```

### `src/avbd/joint.cpp:203`

```cpp
        // The reference implementation also ramps penaltyAng here, with *this* block's C -
        // a length in metres fed to the angular beta, and then the angular block below ramps
        // it again with the real angular error, so a joint with both halves enabled grew its
        // angular penalty at up to twice the rate and used the wrong units for the first
        // half. Only the angular block may touch penaltyAng.
```

### `src/avbd/list_range.hpp:1`

```cpp
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
```

### `src/avbd/list_range.hpp:19`

```cpp
// A forward range over an intrusive singly-linked list: Node must expose a public
// `Node *next`. Empty iteration = null head; the end iterator is a null-node iterator,
// so begin()/end() compare by pointer and the range is trivially multi-pass.
```

### `src/avbd/manifold.cpp:1`

```cpp
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
```

### `src/avbd/manifold.cpp:175`

```cpp
// S = B^T K B in the six independent entries of a symmetric 3x3. The basis rows are the
// linear jacobian rows (up to sign) and K is diagonal, so S is a weighted sum of the
// three basis-row outer products. The expression order matches the old per-stamp code
// exactly, so caching it does not change the result.
```

### `src/avbd/manifold.cpp:198`

```cpp
// The constraint value the primal and dual work against: the step-start error damped by
// (1-alpha), plus the pose bracket in the contact basis. Penetration is damped so a deep
// impact does not overshoot; a speculative contact that is still geometrically apart is
// NOT penetrating yet, and damping its gap to 1% makes any millimetre of jitter look like
// a violation - the pair starts pushing from up to `speculativeMargin` away, an invisible
// inflation layer that measurably deepens every jam in a dense pile. Use the full gap for
// those, so the pair only acts once the bracket predicts a real overlap. With the margin
// disabled the geometry never reports a separated pair, so the branch is unreachable and
// the expression is the original one bit for bit.
```

### `src/avbd/manifold.cpp:214`

```cpp
// Cap the stored penetration (see CONTACT_MAX_PENETRATION) on the normal row only. The two
// tangential rows are the friction error, and their size is a property of where the contact
// anchors sit on the two faces, not of how deep the pair is: in a sheared pile the anchors of
// one contact are routinely decimetres apart along the friction axes. Scaling the whole
// vector by its norm then divides the normal depth down by that tangential distance -
// measured on the double pyramid, a 9 cm penetration arrived at the solve as 0.65 cm, which
// is under the projection slop and buys almost no penalty force, so the pair sank as far as
// it could and stayed there. The tangential rows need no cap here: the friction cone limits
// the force they can produce.
```

### `src/avbd/manifold.cpp:237`

```cpp
// Normal force is one-sided and the tangential pair must stay inside the friction cone.
// The primal force cache, the dual force and the dual's look-ahead force all need the same
// clamp; the expressions and their order are kept verbatim so the solver stays bit-exact.
```

### `src/avbd/manifold.cpp:322`

```cpp
                    // Stick heal: a sticky pair's anchors are glued to each other, so their
                    // separation (what STICK_THRESH checks every dual update) stays ~zero no
                    // matter how far the pair has slid. Compare the frozen anchors against the
                    // fresh contact geometry at the current pose instead: past Solver::stickDriftHeal (default STICK_DRIFT_HEAL)
                    // the "static" contact is really sliding, the old anchors no longer sit on
                    // the touching surfaces, and C0 measures anchor separation - phantom
                    // penetration the dual then answers with real force. Re-anchor from the
                    // current pass and drop the force history (the penalty ramp survives, so
                    // support rebuilds in-round).
```

### `src/avbd/manifold.cpp:347`

```cpp
                        // Keep the normal support impulse: it is a function of how hard the
                        // pair is pressed together, not of where the anchors sit. Dropping it
                        // made every sliding contact lose its support each tick - the pile
                        // never settled, it just kept re-sinking. Only friction is tied to the
                        // anchors that just moved.
```

### `src/avbd/manifold.cpp:366`

```cpp
    // Fresh-contact seed: a contact that matched no feature key starts its penalty at a
    // fraction of betaLin instead of the floor, so an impact-time contact does not spend its
    // first ticks ramping from k=1 while the load runs through it. Gated on the grace
    // setting (default 0 keeps the historical floor start).
```

### `src/avbd/manifold.cpp:389`

```cpp
    // Grace dormancy: the pair separated and no contacts came back. Parking the manifold
    // (instead of reporting it dead for warmstartForces to reap) keeps the hot state of
    // `liveCount` contacts untouched in the arrays above, so a pair that returns within
    // Solver::contactGraceTicks re-matches by feature key instead of rebuilding from the
    // penalty floor. numContacts stays 0, which makes every solve loop over this manifold
    // a no-op.
```

### `src/avbd/manifold.cpp:399`

```cpp
            // Park the stiffness, not the force. What a revived pair needs is the penalty
            // ramp it earned; the lambda from before the separation is support for a pose
            // that no longer exists. Kept alive, it is applied in full by the first primal
            // round after revival - a kick of lambda*dt that a settled box feels as a
            // several-m/s shove out of the pile. The release rule in the live path exists
            // for the same reason: a separated pair may carry no force.
```

### `src/avbd/manifold.cpp:426`

```cpp
        // Cap the depth the solve sees. The BDF1 finish turns the step's displacement into
        // velocity (v = dp/dt), so an unclamped 10 cm impact overlap is expelled at ~24 m/s,
        // that body crosses several neighbours before their next collision pass, and each of
        // those hits injects another ejection - the pile never settles (measured: a double
        // pyramid at 0.6 m overlaps, boxes leaving at 18 m/s). The residual is re-measured and
        // driven out over the following steps.
```

### `src/avbd/manifold.cpp:434`

```cpp
        // Adaptive beta (the augmented-Lagrangian kappa update, per contact): a sticking
        // contact whose penetration did not shrink against last step doubles its personal
        // ramp multiplier; a converging one halves it back toward 1. Fresh contacts start
        // neutral. The residual is penetration depth (the margin is not an error to
        // converge), and a residual at or below the floor counts as solved so float noise
        // on a resting contact cannot ratchet the multiplier to the cap.
```

### `src/avbd/manifold.cpp:460`

```cpp
        // Contact release: the anchors say the pair is geometrically apart, so no contact
        // force can act - whatever impulse the dual still holds came from an earlier pose
        // and must not be applied across the gap. The warm-start decay alone would keep a
        // big impact impulse alive for hundreds of ticks, and the penalty term keeps
        // stamping it into the solve every round; the pile then pops a box out at 10+ m/s
        // long after it was still. The penalty (stiffness) survives, so the support
        // rebuilds in one round when the pair closes again.
```

### `src/avbd/manifold.cpp:473`

```cpp
        // Initial penalty (Eq. 19): the original computed base * 10^(2*log10(pen)) —
        // algebraically pen^2 — and clamped to [PENALTY_MIN, PENALTY_MAX]. pen is capped
        // at 0.02 above, so the product is <= 4e-4 < PENALTY_MIN: the clamp always
        // saturates at the floor and the log10/pow chain was dead computation
        // (~300 us/step across a full manifold set). Penalty starts at PENALTY_MIN and is
        // ramped per iteration in updateDual (betaLin * |C|, capped by PENALTY_MAX).
        // The warmstart decay below is kept bit-for-bit.
```

### `src/avbd/manifold.cpp:481`

```cpp
        // Warmstart the penalty too (Eq. 19). Resetting it to the floor every frame
        // threw away the ramp the previous step had already built, so a resting body
        // had to re-grow its support from PENALTY_MIN each frame and sank until the
        // dual happened to catch up. The reference implementation decays and clamps
        // the carried value instead; a fresh contact starts from the floor because
        // its slots were value-initialised.
```

### `src/avbd/manifold.cpp:506`

```cpp
// Recompute C0 at the current pose, keeping the contact set (rA/rB/basis) from the tick's
// collision pass. The prepared force cache belongs to the old C0, so invalidate it; the
// next primal touch will precompute against the rebased error.
```

### `src/avbd/manifold.cpp:524`

```cpp
// The primal runs in two phases so the two bodies of a manifold share one assembly:
// preCompute (parallel over forces, once per round) evaluates everything both stamps
// need - world arms, angular jacobians, and the clamped contact forces - from the
// round's pose snapshot; stampCached (parallel over bodies) stamps each body's side
// from that cache. Before the split, each of the two per-round updatePrimal calls
// recomputed the whole shared half and used only its own side of the stamps. Same
// inputs in the same order, so the cached values are bit-identical to what the stamps
// used to recompute.
```

### `src/avbd/manifold.cpp:548`

```cpp
    // Constraint value in the basis frame. The linear jacobian rows are the basis rows
    // (±), and each angular row is rWorld x that basis row, so the four matrix-vector
    // products the old assembly built collapse into one bracket the basis then rotates:
    //   C = C0(1-alpha) + B * [ (dqA_lin - dqB_lin) + (dqA_ang x rA) - (dqB_ang x rB) ]
    // No angular jacobian matrices are needed: the basis is the only frame C lives in.
```

### `src/avbd/manifold.cpp:580`

```cpp
    // Both bodies share the linear jacobian up to sign, Jlin = ±basis, so the linear
    // Hessian Jlin^T K Jlin = B^T K B is the same for A and B. K is diagonal, so that
    // "basis sandwich" S is just a weighted sum of the three basis-row outer products.
    // With R = skew(rWorld) the angular jacobian is Jang = ±B R^T, which gives the whole
    // stamp in closed form:
    //   lhsLin   = S
    //   lhsAng   = Jang^T K Jang = R S R^T
    //   lhsCross = Jang^T K Jlin = R S
    //   rhsLin   = ±B^T F
    //   rhsAng   = ±R (B^T F)
    // This replaces the three generic 3x3 products (and their transposes) the old stamp
    // built per contact, per body. The symmetric results are written mirrored so the
    // pivot fallback in Block::solveDelta still sees a full matrix.
```

### `src/avbd/manifold.cpp:597`

```cpp
    // S is refreshed by initialize()/updateDual() whenever the penalty changes, so the
    // per-body stamp does not recompute the same sandwich for the second body.
    // Accumulate in registers, then write back once.
```

### `src/avbd/manifold.cpp:606`

```cpp
    // A zero cached force does not make the contact inactive: the normal force is clamped
    // at zero while the contact can still be compressed within the margin, and its penalty
    // still belongs in the Hessian. Skipping it removes the support stiffness, which lets
    // loaded stacks sink through each other and, once deep enough, launch.
```

### `src/avbd/manifold.cpp:659`

```cpp
    // The first body to touch this manifold in a round consumes the force/jacobian cache
    // the previous round's dual left behind: no pose in this computation has moved since
    // then (the other body has not been solved yet either), so recomputing would produce
    // the same floats. The second body's colour runs after poses have moved, so it
    // recomputes - that is the round's one fresh evaluation.
```

### `src/avbd/manifold.cpp:683`

```cpp
    // Both contact rotations share the two body quats, which never change inside this
    // loop: one rotation-matrix conversion (~20 ns) replaces two quat multiplies per
    // contact (~16 ns each at 8 contacts = ~256 ns/manifold).
```

### `src/avbd/manifold.cpp:689`

```cpp
    // Drift-gated arm reuse (the jacobianRebuildDistance semantics, now living in the
    // dual): when neither body's rotation moved further than the threshold since the last
    // capture, the cached world moment arms are reused for C and for the next round's
    // stamp. Only the arms are cached now - the angular jacobian matrices are gone.
```

### `src/avbd/manifold.cpp:746`

```cpp
        // Leave the primal's cache ready for the next Newton round. The pose is frozen
        // between here and the next stampCached (the solve runs after the next
        // preCompute would), so the next preCompute would recompute exactly this bracket
        // - and, with the ramped penalty and the just-stored lambda, exactly this force.
        // Doing it here removes the standalone preCompute from every round after the first.
```

### `src/avbd/manifold.cpp:767`

```cpp
// Accumulate this manifold's penetration error for the convergence probe. With
// speculative contacts on, a positive normal error is a gap the bodies have not closed
// yet, not penetration: counting it would skew the residual probe.
```

### `src/avbd/manifold.hpp:1`

```cpp
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
```

### `src/avbd/manifold.hpp:27`

```cpp
// Collision manifold types. Manifolds themselves live in the solver's ManifoldStore
// (structure of arrays, stable slots) rather than in per-pair objects, so this struct is
// a namespace for the shared records plus the solver-free pair query.
```

### `src/avbd/manifold.hpp:46`

```cpp
    // Contact point geometry. This is what collision detection writes; the hot solver
    // state lives in ManifoldStore's parallel arrays, so the per-contact stride the primal
    // and dual walk stays small.
```

### `src/avbd/manifold.hpp:56`

```cpp
    // The symmetric linear Hessian B^T K B, six independent entries: K is diagonal and
    // the basis is orthonormal, so every pair of bodies of this manifold shares it. It
    // changes only when the penalty changes, so initialize()/updateDual() refresh it once
    // per round and both body stamps read it.
```

### `src/avbd/manifold.hpp:70`

```cpp
// Every collision manifold of one solver, as parallel arrays over stable slots.
//
// A slot is allocated from `freeList` (or appended) and stays put until the manifold is
// released; `live` lists the slots that currently own a manifold, in creation order, so
// every per-manifold phase is a compact walk over `live`. Contact state is indexed
// `slot * kMaxManifoldContacts + k`; the fixed stride wastes a partly unused tail on
// manifolds with few contacts, but keeps warm-start matching, parallel per-slot writes
// and the per-slot passes trivial.
```

### `src/avbd/manifold.hpp:127`

```cpp
    // Refresh the contact set for the step and warm-start the hot state. Returns false
    // when the pair separated (and is not parked in the contact-grace window); the slot
    // stays allocated either way and the caller decides whether to release it.
```

### `src/avbd/manifold.hpp:134`

```cpp
    // Two consumers share one cache. updatePrimal (the coloured Gauss-Seidel path) uses
    // the force/moment-arm cache the previous round's dual left behind: the first body to
    // touch the manifold this round consumes it for free, and the second colour - by then
    // poses have moved - recomputes once (preCompute) and stamps its side. The overflow
    // bucket's Jacobi path instead runs preCompute over slots, then stampCached over
    // bodies, so one assembly serves both sides against the same snapshot.
```

### `src/avbd/manifold.hpp:170`

```cpp
// Per-body constraint adjacency in CSR form, keyed by Rigid::broadIndex. Rebuilt every
// step (before the broad phase for the pair filter, after warmstart for the solve): the
// colouring, the primal walk, the sleep wake-up and the island grouping read it instead
// of following the intrusive body lists. Manifold references carry a store slot.
```

### `src/avbd/manifold.hpp:195`

```cpp
// Discrete shape-pair collision query: up to 8 contact points between two `Shape`s,
// independent of any solver. `contacts` must hold at least 8 slots (`Manifold`'s
// capacity) - exported because space queries reuse the exact path `Manifold` walks.
// `basis` receives an orthonormal frame whose first row is the contact normal (A -> B).
// `speculativeMargin` lets separated shapes report contacts up to that distance apart
// (the constraint error is then positive); 0 keeps the query strictly overlap-only.
// Returns the contact count (0 = separated).
```

### `src/avbd/maths.h:1`

```cpp
/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 *
 * Math layer over linalg.h (sgorsten, Unlicense). The types are linalg aliases
 * with the storage layout the rest of the engine always assumed:
 *   float3   linalg::vec<float,3>          - {x, y, z}
 *   float3x3 linalg::mat<float,3,3>        - column-major, m[col][row]
 *   quat     linalg::vec<float,4>          - {x, y, z, w}, Hamilton product
 *
 * Differences from Eigen that the port had to respect:
 *   - linalg's `operator*` is component-wise between two operands of the same
 *     shape (cmul); algebraic products are `mul(...)`, `dot(...)`, `qmul(...)`,
 *     `qrot(...)` and `outerprod(...)`.
 *   - Eigen's `m(row, col)` maps to linalg's `m[col][row]`.
 *   - Eigen's Quaternion ctor was (w, x, y, z); linalg stores {x, y, z, w}.
 */
```

### `src/avbd/maths.h:40`

```cpp
// The linalg free functions the engine calls unqualified, re-exported so the
// call sites read `dot(a, b)` / `normalize(v)` without a linalg:: qualifier.
// Scalar overloads of std::min/std::max/std::clamp are always named with the
// std:: prefix at the call sites, so nothing here collides.
```

### `src/avbd/maths.h:58`

```cpp
// Eigen reduces a 4-vector's squared norm as (x^2+z^2)+(y^2+w^2) (its SSE lane
// pairing), while linalg folds left-to-right; the orders round differently, so the
// quaternion helpers below keep the reference rounding. Non-template overloads win
// over the linalg templates for vec<float,4> arguments; 2- and 3-component vectors
// fall through to linalg's, which already match Eigen's scalar order.
```

### `src/avbd/maths.h:78`

```cpp
// Quaternion rotate and rotation matrix, written with Eigen's operation order
// (t = 2*(z x v); v + w*t + z x t, and the tx/ty/... factorisation). linalg's own
// qrot/qmat expand different-but-equivalent expressions whose rounding differs by
// ulps; on stiff, chaotic scenes (joint lattices, pyramids) that drift compounds
// into visibly different behaviour, so the port keeps the reference rounding.
```

### `src/avbd/maths.h:97`

```cpp
    // Column-major: each braced float3 is a COLUMN, so column j reads (R0j, R1j, R2j).
    // The R12/R21 pair is (tyz - twx, tyz + twx); exchanging them makes the matrix neither
    // a rotation nor a transpose, and everything built on it silently uses the wrong
    // orientation - the broad-phase bounds then exclude pairs that really overlap.
```

### `src/avbd/maths.h:126`

```cpp
// Right-handed orthonormal frame with the given unit vector as its first ROW (the contact
// normal). The manifold code reads the frame axes with row()/mul(), which use the matrix's
// rows, so the axes must end up there: linalg matrices are column-major, hence transpose.
```

### `src/avbd/maths.h:165`

```cpp
// First-order quaternion integration, as the original custom maths did: half of omega (as a
// pure quaternion) left-multiplied, then renormalised. The added rotation vector is the same
// small-angle quantity `operator-(quat, quat)` measures.
```

### `src/avbd/maths.h:176`

```cpp
// Solve the symmetric 6x6 system [aLin  aCross^T; aCross  aAng] x = [bLin;
// bAng].
//
// Two kernels, one bounded-update policy:
//   1. Fixed-size LDL^T - stack-only, ~200 ns kernel. The workhorse.
//   2. Fixed-size pivoted LDL^T - diagonal pivoting, different rounding. Retry
//      path when (1) is not sane.
//   3. Zero update - last resort when neither kernel is sane.
// "Sane" means finite and < 10 m per component: an update is velocity*dt plus a position
// correction, so 10 m is orders past anything physical. Degenerate blocks (thin cylinders
// standing on end produce near-singular stiffness ratios ~1e6) used to let one wild pivot
// compound into e15-scale positions; the policy caps every solve at the physical bound.
// The fixed kernel's rounding differs from the old dynamic-Eigen implementation, which
// re-anchors every state digest - accepted on the real-time branch.
// Fixed-size 6x6 LDL^T without pivoting, over an already-assembled system. The
// system is dominated by the M/dt^2 diagonal (a 1 kg body at 240 Hz carries a
// 57,600 kg/s^2 diagonal against ~1e4 penalty terms), so the pivot growth symmetric
// pivoting exists for does not occur here, and dropping it lands the kernel an order
// of magnitude cheaper. Plain sequential float ops keep the result bit-identical for
// every thread count. Returns false on a degenerate pivot (non-positive D, which a
// near-singular block can produce) or an unsane solution - callers retry with the
// pivoted path in solve() below.
```

### `src/avbd/maths.h:242`

```cpp
// Pivoted LDL^T, a faithful port of Eigen's LDLT (ldlt_inplace<Lower>::unblocked
// plus LDLT::_solve_impl_transposed) over the fixed 6x6 arrays. This is the
// rounding fallback the Eigen dynamic-LDLT retry used to provide, and it must
// reproduce Eigen's behaviour bit for bit, because the retry fires on indefinite
// blocks (negative geometric stiffness can push a diagonal below zero) where the
// two factorisations otherwise pick different pivots and the solutions diverge:
//   - pivot = largest |diagonal| of the trailing block, first maximum wins;
//   - the pivot threshold is exactly zero (a negative pivot is a valid pivot);
//   - a zero pivot fails only if the remaining column is non-zero, or when a
//     later step is non-zero after a zero one was accepted;
//   - the solve divides by |d| > numeric_limits<float>::min() and zeroes the
//     component otherwise (Eigen's pseudo-inverse of D);
//   - the permutation is applied as running swaps, P forward then P^T backward.
```

### `src/avbd/obb.hpp:1`

```cpp
/*
 * Conservative OBB overlap test for the broad phase.
 *
 * The broad phase buckets bodies by their tight world AABB, but the exact candidate
 * filter uses an oriented box instead of the world AABB: a rotated box's AABB can be
 * several times its true footprint, and every false positive becomes a Manifold that
 * collide() has to reject and warmstartForces() has to delete.
 *
 * Spheres are represented by an axis-aligned cube of half extent = radius (a sphere is
 * rotation invariant, so this is its exact AABB). Cylinders use their local enclosing box.
 * Both enclosures contain the real shape, so the test can only keep too many pairs, never
 * drop a real overlap.
 */
```

### `src/avbd/obb.hpp:30`

```cpp
// 15-axis SAT (6 face axes + 9 edge axes). `axes` are unit and orthogonal, so the
// projections reduce to table reads; no normalisation is needed. A small epsilon biases
// the test towards overlap (a false positive is cheap, a false negative is a lost contact).
```

### `src/avbd/query/query_shape.cpp:451`

```cpp
    // The AABB interval is conservative. Walk it once to find the first contact, then
    // bisect that bracket for a stable safe/unsafe pair. The interval already bounds the
    // sweep, so a fast motion cannot skip past the target.
```

### `src/avbd/query/query_shape.hpp:1`

```cpp
/*
 * Space-query geometry for the AVBD solver core.
 *
 * This is an extension layer, not part of the solver. It only reads the public
 * avbd::Shape description and calls the existing collideShapes() narrow phase, so the
 * teaching-oriented core (solver.h/solver.cpp/collide.cpp) stays untouched.
 *
 * Shapes use exactly the proxy the solver simulates: a capsule query is a cylinder, a
 * convex-hull query is the hull's local AABB box. Keeping the query and the contact
 * model on the same geometry is deliberate - a raycast may not report a surface the
 * solver does not collide against.
 */
```

### `src/avbd/query/query_shape.hpp:50`

```cpp
// Penetration along a contact normal: support(a)+support(b)-|centre separation|. This is
// the SAT overlap on that axis and gives rest queries a real depth instead of a distance
// between two arbitrary anchor points.
```

### `src/avbd/query/query_shape.hpp:64`

```cpp
// Translation sweep of `base` along `motion` against a static `target`.
//
// `safe` is the last tested fraction without a contact, `unsafe` the first with one;
// `margin` has the same meaning as collideShapes()' speculative margin. The sweep is
// bounded by the swept AABB overlap interval, so a fast motion cannot tunnel through a
// target just because its final pose is past it.
```

### `src/avbd/rigid.cpp:1`

```cpp
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
```

### `src/avbd/solver.cpp:1`

```cpp
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
```

### `src/avbd/solver.h:1`

```cpp
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
```

### `src/avbd/solver.h:38`

```cpp
// Work per item differs by orders of magnitude between the phases below, so each one gets
// its own break-even batch size. Measured on a 12-core machine: one dispatch into the pool
// costs ~45 us, a body's primal update ~1 us, one manifold's contact solve ~85 ns, and a
// dual update ~200 ns. Break-even is therefore roughly 45 bodies, 500 contact solves and
// 250 dual updates - i.e. a scene has to be a few hundred bodies before splitting a phase
// across threads pays for the dispatch itself.
```

### `src/avbd/solver.h:46`

```cpp
// The broad phase's per-body preparation (AABB bound, level, cell key) costs tens of
// nanoseconds per body, so a dispatch only pays off at a much larger batch than a body's
// primal update does (~1 us). Below this the loop runs inline and saves the round trip.
```

### `src/avbd/solver.h:67`

```cpp
    // Primal rounds the current step actually runs (set by step()/stepInterleaved before
    // any constraint is touched). A constraint penalty is only as stiff as the iteration
    // budget can resolve, so the joint penalty cap reads this: the augmented-Lagrangian
    // dual gains `penalty * C` every round, and a penalty that a long iteration set
    // converges with makes a short one diverge into a launcher. See
    // GenericJoint::penaltyLimit.
```

### `src/avbd/solver.h:75`

```cpp
    // Round split: the first ceil(newtonRatio * rounds) rounds run the full primal (the
    // dense 6x6 block solve); the rest skip the primal and run dual-only relaxation on the
    // existing forces, with contact stiffness decayed by `stiffnessDecay` per round so the
    // cheap rounds finish what the Newton rounds started instead of fighting a stiff ramp.
    // 1.0 / 1.0 = every round is a Newton round (historical behaviour).
```

### `src/avbd/solver.h:83`

```cpp
    // Adaptive beta: a per-contact multiplier on betaLin, carried across steps with the
    // warm start. A contact whose penetration did not shrink against last step doubles its
    // multiplier (capped at ADAPTIVE_BETA_MAX); a converging one halves it back toward 1.
    // This is the augmented-Lagrangian kappa update applied per constraint: a loaded stack
    // contact stiffens until it actually carries its share, instead of every contact
    // sharing one global ramp rate.
```

### `src/avbd/solver.h:90`

```cpp
    // Residual penetration (metres) at or below which a contact counts as solved: the
    // multiplier decays instead of doubling, so float noise on a resting contact cannot
    // ratchet it to the cap.
```

### `src/avbd/solver.h:108`

```cpp
    // Early-exit threshold on the mean penetration residual: when a Newton round
    // improves the average residual by less than this (or worsens it), the iteration set
    // ends early - remaining rounds would polish float noise at full price. 0 disables
    // the check. Realistic value ~1e-5.
```

### `src/avbd/solver.h:114`

```cpp
    // Combine-pass relaxation for the overflow bucket of the coloured primal - the
    // bodies the greedy colouring could not fit inside kMaxColours, which run as one
    // block-Jacobi round against a pose snapshot. Such a snapshot solve can overshoot
    // where the coloured Gauss-Seidel updates do not; scaling the applied delta by
    // omega < 1 trades per-iteration progress for stability. The
    // bucket is empty in ordinary scenes, so these are tuning knobs for pathological
    // constraint cliques. Suite sweep: 0.5-0.7 hold every rigid bucket scene; 1.0
    // drifts a fully-cliqued stack.
```

### `src/avbd/solver.h:124`

```cpp
    // Route every body whose forces are all contact manifolds (Rigid::contactOnly) into the
    // block-Jacobi overflow bucket instead of the coloured Gauss-Seidel. True in a pile of
    // boxes - which makes the bucket the whole primal for a contact-only scene, not an
    // exception. False runs those bodies through the colouring.
    //
    // A/B on the 420-box double pyramid (`test_penetration.gd --pen-seconds=30`, identical
    // wall time, 626 ms / 600-frame pyramid run either way):
    //                        p50 pen   p99 pen   worst *settled*   max speed
    //   true  (Jacobi)        8.0 mm    11.3 mm        15.3 mm        17.9 m/s
    //   false (colour GS)     7.4 mm     8.8 mm        36.3 mm        18.4 m/s
    // Gauss-Seidel improves the median (it is the better solver) but leaves a systematic
    // tail: a pair that the colour order revisits last settles 3.6 cm deep, because the
    // update order is a property of the scene, not of the load. Jacobi averages those
    // order effects away and has the better worst case, so it stays the default until the
    // colouring gets an ordering pass (lowest body first, or a random sweep per round).
    // Exposed as `physics/avbd/contact_gs` on the Godot side.
```

### `src/avbd/solver.h:142`

```cpp
    // Anchor drift (metres) at which a sticky contact's frozen arms are re-adopted from the
    // current collision pass. STICK_DRIFT_HEAL is the default; lowering it re-anchors more
    // eagerly, which trades warm-start continuity for a C0 that always measures the surfaces
    // that are actually touching.
```

### `src/avbd/solver.h:148`

```cpp
    // Depth (metres) the constraint solve is allowed to see for one contact. A position-level
    // step turns the correction it applies into velocity (v = dp/dt), so an unbounded deep
    // overlap is expelled as a launch. This is the historical cap; see Solver docs.
```

### `src/avbd/solver.h:153`

```cpp
    // Incremental arm-rebuild distance (metres): later rounds reuse a contact's cached
    // world moment arms while the body rotations have drifted less than this from the
    // captured pose. 0 disables caching (fresh arms every round; default).
    // Realistic value 1e-4 (0.1 mm): in-iteration arm drift is millimetre-scale.
```

### `src/avbd/solver.h:164`

```cpp
    // Jacobi passes measure all contacts before moving bodies, avoiding manifold-order pose bias.
    // `projectionIterations` is an upper bound: a pass whose worst remaining depth is below
    // 1% of projectionSlop stops the loop, and sleeping/parked pairs are not walked at all.
```

### `src/avbd/solver.h:174`

```cpp
    // Ticks a manifold whose pair slipped out of the swept AABBs stays dormant instead of
    // being deleted. A reaped manifold loses its whole hot state - penalty ramp, duals,
    // stick flags - so the impact-time feature churn (contacts flickering in and out as
    // faces toggle) restarts every returning contact from the penalty floor, and a loaded
    // stack pays that restart hundreds of times per tick. A dormant manifold keeps its
    // contacts and hot state (it contributes nothing while separated: every solve loop
    // walks contacts) and re-matches by feature key when the pair returns. Freshly formed
    // contacts inside a live manifold also seed their penalty from 5% of betaLin instead
    // of the floor while this is enabled. 0 = delete on first separation, as before.
```

### `src/avbd/solver.h:185`

```cpp
    // Worker threads for the per-body phases. 0 = one per hardware thread, 1 = run
    // everything inline on the calling thread. The result does not depend on this value:
    // every body solves its own block system against the same pose snapshot, and applies
    // its own delta.
```

### `src/avbd/solver.h:191`

```cpp
    // Idle detection: a movable body whose BDF1 velocities stay under both thresholds for
    // this many consecutive steps falls asleep (velocity zeroed, Godot semantics).
    // 0 disables automatic sleeping; the core tests rely on the default for bit-identical
    // digests. The horizon is a TIME (Godot's 0.5 s) converted to steps by the host, not a
    // fixed 60 frames. The server layer overrides all three from physics/avbd/sleep_*
    // settings (defaults 0.05 m/s, 0.05 rad/s, 1.0x horizon); the core defaults here stay
    // as they were so pinned digests are unaffected.
```

### `src/avbd/solver.h:209`

```cpp
    // Bodies and forces are linked into the solver and into each other by raw pointers, and
    // a step deletes forces (see warmstartForces), so neither can be copied or moved: a copy
    // would share the list heads and double-delete them.
```

### `src/avbd/solver.h:219`

```cpp
    // Masked variant with Godot's pair semantics: a body is a candidate when
    // (body.layer & p_mask) != 0. The two-argument form above forwards with an
    // accept-everything mask. Both walk the last built broad-phase tree (ray/AABB
    // pruned); before any step has built one, they fall back to a body-list walk.
```

### `src/avbd/solver.h:224`

```cpp
    // Bodies whose last-built broad-phase bounds overlap the query AABB, filtered by
    // (layer & p_mask) != 0. Backs the direct space state's shape and point queries;
    // the same no-tree fallback applies.
```

### `src/avbd/solver.h:235`

```cpp
    // Diagnostics only: did the last broad phase emit this pair as a narrow-phase candidate?
    // A pair that genuinely overlaps and answers no here is a broad-phase miss, not a solve
    // failure. Linear scan - the callers run at diagnostic rates.
```

### `src/avbd/solver.h:246`

```cpp
    // Diagnostics for a probe run: when a pair matches these two bodies (either order), the
    // broad phase, the manifold pass and the reaper print every lifecycle event for it. A
    // pair with a real overlap and no constraint needs exactly this to be attributed.
```

### `src/avbd/solver.h:298`

```cpp
    // Decide which bodies the primal phases may move this step (wake-up, island sleep and
    // the sleep-mode gate), size the per-body delta buffers, and colour the constraint
    // graph for the Gauss-Seidel primal.
```

### `src/avbd/solver.h:302`

```cpp
    // Interleaved substeps: one broad phase / force warmstart / update-set selection per
    // tick, then `substeps` rounds of (body warmstart, primal+dual, velocity finish) with
    // the contact error rebased between them. `roundsPerSubstep` 0 splits the tick's
    // `iterations` evenly. `dt` must already be the substep dt.
```

### `src/avbd/solver.h:314`

```cpp
    // Per-phase wall time of the last step, milliseconds. Diagnostics only - nothing in the
    // solver reads it - and it is what the host's get_*_time() accessors report, so a profile
    // can attribute a step without a second copy of the schedule. A host that reads it from
    // another thread while a step runs gets the previous step's numbers.
```

### `src/avbd/solver.h:321`

```cpp
        // The algorithm under the broad phase alone: the LBVH build/refit plus candidate
        // emission, excluding the bounds/AABB pass, the pair sort, the exact filter and
        // manifold slot creation. This separates the broad-phase index cost from the
        // shared pipeline around it.
```

### `src/avbd/solver.h:341`

```cpp
    // Every body in list order. The inertial/warmstart and BDF1 phases run over this,
    // because they also have to keep the cached state of static bodies current: contact
    // and joint code reads `initialLin`/`initialAng` of *both* bodies of a constraint, and
    // the ground is usually one of them.
```

### `src/avbd/solver.h:347`

```cpp
    // Broad-phase: the body list flattened to index-addressable order (rebuilt every
    // step), the Morton LBVH over those bodies, and the candidate pair list (sorted
    // unique (min index, other), the manifold creation order).
```

### `src/avbd/solver.h:351`

```cpp
    // Swept AABBs per bodiesInOrder: the tight bounds grown by this step's linear
    // displacement, the rotation sweep and the speculative margin. center is the
    // current position plus half the displacement; half is the matching extent.
```

### `src/avbd/solver.h:363`

```cpp
    // Bodies that take part in the primal update (movable, awake ones), in list order.
    // The coloured Gauss-Seidel primal solves one block system per body of each
    // colour and applies it immediately; bodies the greedy colouring cannot fit
    // inside kMaxColours land in an overflow bucket that runs one block-Jacobi
    // round (solve every delta against the same snapshot, then combine), using
    // the delta buffers below.
```

### `src/avbd/solver.h:379`

```cpp
    // Constraint-graph colouring, rebuilt by colourGraph() every step. Bodies of
    // one colour share no force, so a colour's whole set can solve and apply in
    // parallel (Gauss-Seidel between colours, one barrier each). bodyColour is
    // indexed by updateOrder position and holds either a colour < kMaxColours or
    // kMaxColours itself for the block-Jacobi overflow bucket.
```

### `src/avbd/solver.h:398`

```cpp
    // The non-manifold force list rebuilt each step so the phases can address joints,
    // springs and ignore pairs by index, plus the per-force result of initialize() so the
    // deletion can be done in one thread after the parallel initialisation.
```

### `src/avbd/solver.h:426`

```cpp
    // The whole iteration set as one persistent-worker region: `targetIterations` rounds of
    // the block-Jacobi primal (one full block solve per body against the iteration-start
    // pose, then a combine pass) plus dual passes, fenced by LoopSync barriers. Both sync
    // types run the identical body; NullSync is the inline/reference path.
```

### `src/avbd/solver.h:448`

```cpp
    // Island-wide sleep bookkeeping (runs when sleeping is enabled): movable bodies are
    // grouped into islands through force connectivity, and an island whose every body
    // cleared idle detection this step puts all its members to sleep together. Prevents
    // the partially-asleep contact chain from thrashing between states.
```

### `src/avbd/solver_broadphase.cpp:1`

```cpp
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
```

### `src/avbd/solver_broadphase.cpp:46`

```cpp
    // Two static bodies can never need a contact: neither can move, and the solver has no
    // update for either. Without this the pair still got a manifold slot in broadPhase,
    // and `ManifoldStore::frozen()` then reported it alive every step (both sides are
    // immovable), so `warmstartForces` kept a zero-contact manifold forever - one store
    // slot and eight contact records per overlapping pair of static bodies, on a terrain
    // made of static boxes for instance. A pair that both *sleep* is a different thing: its
    // manifold carries the state that has to survive until it wakes, so only mass 0 counts.
```

### `src/avbd/solver_broadphase.cpp:72`

```cpp
// A sleeping body whose swept AABB touches a moving movable body wakes up. AABBs are
// already grown by the speculative margin, so a touch here means "these two meet within the
// tick"; no real contact is needed. Serial by construction (called from broadPhase's pair
// loop): the fields are plain bools/ints and two pairs sharing a body used to write them
// from two threads.
```

### `src/avbd/solver_broadphase.cpp:91`

```cpp
// Contact detection: the Morton LBVH broad phase, rebuilt or refit every step.
//
// It indexes the swept world AABBs and emits packed pair candidates. The surrounding
// pipeline - pair sort, layer/constraint filter, manifold creation - is shared.
//
// Determinism: pair collection is integer-sorted and de-duplicated, and the final
// (min-index, other) sort fixes the manifold creation order.
```

### `src/avbd/solver_broadphase.cpp:107`

```cpp
    // Pair existence and AABB wake-up read the constraint graph, so it has to be current
    // before the candidate filter runs. Built again after warm-starting, when the live
    // manifold set has changed.
```

### `src/avbd/solver_broadphase.cpp:115`

```cpp
    // Swept bounds: the union of the body's AABB at the step start and at the step end.
    // For linear motion that union is exact - centred at p + disp/2 with half extent
    // h + |disp|/2 - and the rotation term conservatively grows it by the arc a point at
    // the bounding radius sweeps. The speculative margin is the last growth, which is
    // what lets a pair whose surfaces are about to meet enter the narrow phase.
```

### `src/avbd/solver_broadphase.cpp:172`

```cpp
        // 睡眠 + awake 的 AABB 重叠就唤醒。不需要真实接触——AABB 已经含 margin，
        // 重叠意味着"下几个 substep 内会碰"。
        //
        // 这一步必须在串行阶段：`sleeping`/`stillFrames` 是普通字段，两个候选 pair
        // 在并行 filter 里同时写同一个 body 是数据竞争（C++ 里就是 UB，TSAN 会报）。
        // 放在这里还有一个好处：顺序由排好序的 pair 表决定，唤醒结果与线程数无关。
```

### `src/avbd/solver_broadphase.cpp:181`

```cpp
        // A dormant manifold (contact grace) let this pair back into the candidate set
        // (pairConstrained skips dormant ones so overlap can wake sleepers); it will
        // revive in its own initialize(), so creating a second manifold here would double
        // the pair's constraints.
```

### `src/avbd/solver_internal.hpp:1`

```cpp
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
```

### `src/avbd/solver_internal.hpp:18`

```cpp
// Godot BodyAxis bits: one helper for every place a locked component must read back as
// motionless (the positional projection, the warm start and the final velocity). Keeping the
// bit tests here means the three sites cannot silently disagree about which axis is locked.
```

### `src/avbd/solver_internal.hpp:69`

```cpp
// True when a force already joins the pair. Mirrors Rigid::constrainedTo: any joint,
// spring or ignore pair counts, and a live manifold counts, but a dormant one (contact
// grace) does not - the pair must still wake on AABB overlap and re-enter the narrow
// phase, or a sleeping body behind a dormant manifold is invisible to an incoming body.
```

### `src/avbd/solver_projection.cpp:1`

```cpp
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
```

### `src/avbd/solver_projection.cpp:56`

```cpp
// A locked axis does not move, so the projection drops that component. The effective mass
// was computed as if it did, which makes the correction slightly soft on locked bodies -
// the same documented approximation the rest of the solver makes for axis locks.
```

### `src/avbd/solver_projection.cpp:79`

```cpp
    // Only manifolds whose pair can still move are worth walking. A parked or sleeping
    // pair cannot change velocity or pose this step, so every contact test against it
    // would be a no-op; compacting the live set once keeps the brake and all projection
    // iterations over the same small list instead of re-filtering every slot.
```

### `src/avbd/solver_projection.cpp:97`

```cpp
    // Interception. A pair may only close the gap it has at this pose within the time left
    // before the contact set is refreshed again; otherwise a fast body travels a whole
    // substep into its neighbour and the solver spends the next round undoing that overlap,
    // which is what converts a collision into a launch.
```

### `src/avbd/solver_projection.cpp:181`

```cpp
    // Nothing below this residual is worth another full pass: the depth is already
    // measured past the slop, so a fraction of the slop is a correction that cannot
    // show up above float noise in the next pose. Resting and sleeping piles therefore
    // stop after the first pass that has nothing left to move.
```

### `src/avbd/solver_query.cpp:1`

```cpp
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
```

### `src/avbd/solver_query.cpp:41`

```cpp
    // Exact per-body test against the body's own shape. Picking scans the body list rather
    // than walking the broad phase; the exact test decides.
    //
    // The local half extents are per shape, and this used to be `size * 0.5f` for all of
    // them. For a sphere `size.x` is the *radius* (Rigid's constructor and aabbHalfExtent
    // both read it that way, and the server stores {radius, 0, 0}), so a sphere was tested
    // as a zero-thickness box of half-extent r/2 in y and z - a ray could not hit one at
    // all. A cylinder had its radius halved the same way.
```

### `src/avbd/solver_solve.cpp:1`

```cpp
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
```

### `src/avbd/solver_solve.cpp:39`

```cpp
    // One pass both reaps and compacts: `forceOrder` is list order, and deleting a force
    // only unlinks it from that list, so the survivors keep the exact order a second
    // collectForces() would produce.
```

### `src/avbd/solver_solve.cpp:67`

```cpp
    // Dormant manifolds (contact grace) keep their slot and hot state for a few ticks
    // instead of being reaped. initialize only enters dormancy while the grace setting is
    // positive, so with the default 0 this branch never fires.
```

### `src/avbd/solver_solve.cpp:109`

```cpp
// Compute one body's inertial position and warm-start it (Eq. 2 and the adaptive
// warm-start of the original VBD paper).
//
// Per-body and independent, so `runSubstep` dispatches it as an independent batch. It
// covers static bodies as well: their poses do not move, but the constraint code reads
// their cached `initialLin`/`initialAng`.
```

### `src/avbd/solver_solve.cpp:119`

```cpp
    // Axis locks (Godot BodyAxis): a locked axis cannot carry velocity into the step,
    // and gravity does not act along a locked linear axis. Constraint coupling can
    // still displace a locked axis within a step; documented approximation.
```

### `src/avbd/solver_solve.cpp:167`

```cpp
    // Idle detection: a movable, non-sleeping body whose velocities stayed under both
    // thresholds for sleepFrames consecutive steps falls asleep (velocity zeroed,
    // Godot semantics). Only touches this body's fields, so it is safe inside the
    // parallel dispatch and deterministic regardless of thread count.
    // Automatic sleep only for bodies whose host allows it (sleep_mode == SLEEP).
    // A body with can_sleep=false used to enter idle sleep here anyway; the next
    // selection then simulated it again (sleep_mode NEVER), but the finish stage
    // skipped its velocity update because `sleeping` was still set - the body
    // moved with a permanently zero velocity, and contacts could pump it upward
    // indefinitely. That was the pyramid_dbl runaway.
```

### `src/avbd/solver_solve.cpp:192`

```cpp
// Island-wide sleep bookkeeping. Movable, sleep-enabled bodies are grouped into islands
// through force connectivity (union-find over the flattened body list; static bodies and
// sleep_mode == NEVER bodies join nothing and never propagate). An island where every
// member cleared idle detection this step (stillFrames >= sleepFrames or already
// sleeping) is put to sleep as a unit: Godot semantics, and it removes the
// partially-assembled-chain edge case where the bottom box qualifies a step before the
// boxes above it. Deterministic: the union-find walks the force list in list order and
// the minimum-body-index island representative is independent of visit order.
```

### `src/avbd/solver_solve.cpp:202`

```cpp
    // Only bodies that can sleep join an island. With the default sleepFrames = 0 this
    // whole pass is skipped by the caller, so the core suite (which never sets
    // sleep_mode) is unaffected.
```

### `src/avbd/solver_solve.cpp:275`

```cpp
// The primal half of one iteration for a single body (Eqs. 4-6): assemble the 6x6 system
// from every force acting on the body, then solve it and move the body.
//
// Solve the assembled system and move the body by the result (Eq. 4).
```

### `src/avbd/solver_solve.cpp:293`

```cpp
    // Wake-up pass: a sleeping body (sleep_mode enabled) sharing a force with a moving
    // movable body is woken before selection, so it re-enters this step's update set.
    // A neighbour that is itself below the idle thresholds is not a disturbance: stacks
    // must be allowed to fall asleep together, one body at a time. Serial and read-only
    // over the force adjacency apart from the examined body itself, so the outcome is
    // fixed regardless of thread count.
```

### `src/avbd/solver_solve.cpp:334`

```cpp
        // Sleep mode support (Godot 4.7 API)
        // SLEEP_MODE_NEVER: always update
        // SLEEP_MODE_SLEEP: update unless idle detection has put it to sleep
        // SLEEP_MODE_START_IN_SLEEP: start in sleep
```

### `src/avbd/solver_solve.cpp:344`

```cpp
                // Idle detection (finishBody) is the only way into sleep here;
                // a sleeping body with contacts stays out of the update set until the
                // wake-up pass above re-admits it.
```

### `src/avbd/solver_solve.cpp:382`

```cpp
    // Overflow-bucket eligibility: a body carrying any non-manifold force (joint, spring,
    // world anchor) is routed into the overflow bucket by colourGraph - the contact map
    // is clamped, hard joints do not tolerate the snapshot solve. Every force in the
    // list is non-manifold; contacts do not count.
```

### `src/avbd/solver_step.cpp:1`

```cpp
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
```

### `src/avbd/solver_step.cpp:42`

```cpp
// One tick. The broad phase / force warmstart / update-set selection run once on the
// calling thread (they parallelise their own stages), then `runSubstep` runs the per-body
// warmstart, the colour-parallel iteration region and the velocity finish. Interleaved
// substepping calls `runSubstep` more than once per tick.
```

### `src/avbd/solver_step.cpp:79`

```cpp
// One substep: body warmstart (independent), the colour-parallel Gauss-Seidel region,
// then the BDF1 velocity finish (independent). Idle detection only runs on the tick's
// final substep, so substeps cannot sleep bodies faster than ticks.
```

### `src/avbd/solver_step.cpp:110`

```cpp
// Interleaved substepping: collision/forces/update-set once per tick, then a small number
// of iterations per substep with the contact error rebased at each substep start. This is
// the cheaper alternative to running the full iteration set once per substep; it trades
// some constraint freshness (the contact set is not refreshed inside the tick) for a
// much lower round count. `dt` is already the substep dt.
```

### `src/avbd/solver_step.cpp:138`

```cpp
    // Rounds per substep. An explicit `roundsPerSubstep` is an override: the caller asked
    // for exactly that many primal rounds in each substep, and honouring it is what makes
    // the setting mean something (it used to be a floor, so asking for 2 while `iterations`
    // implied 3 silently ran 3 and the knob could not be lowered below the split). 0 asks
    // for the automatic split of the tick's `iterations` across the substeps.
    //
    // This number is a *stability* parameter, not just a cost one: at a fixed 12 rounds per
    // tick a four-link welded chain settles with 1 or 2 substeps and diverges with 4 or more
    // (measured, tools/repro_joint.cpp: |vz| = 0.07 / 0.06 / 0.31 / 140 / 413 / 1310 m/s for
    // 1x12 / 2x6 / 3x4 / 4x3 / 6x2 / 12x1). Fewer substeps means more rounds inside one
    // linearisation, and the augmented-Lagrangian dual converges inside that linearisation.
```

### `src/avbd/solver_step.cpp:161`

```cpp
        // The contact set is only refreshed once per tick, so every substep re-brakes
        // approach speeds against the gap measured at its own start pose, over the time
        // still left in the tick.
```

### `src/avbd/solver_step.cpp:171`

```cpp
    // The tick's stored error is the one the solver saw when its contacts were initialized;
    // the bodies have moved since. Rebase it at the final pose so what the host reads back
    // (the penetration diagnostic) is the geometry the pair actually ended the tick with,
    // not a one-tick-old pose. The next tick's initialize() replaces it from its own
    // collision pass either way; this only feeds the diagnostics and the host.
```

### `src/avbd/solver_step.cpp:193`

```cpp
    // Pass the value straight through: 0 means "one worker per hardware thread" and 1 means
    // a single worker, which both the pool and `threads` document the same way. (Clamping 0
    // up to 1 here would silently turn the default into a serial run.)
```

### `src/avbd/solver_step.cpp:213`

```cpp
// One region body for every iteration: colour-by-colour Gauss-Seidel primal (a colour is an
// independent set, so the barrier carries the ordering), then the dual pass. Bodies that do not
// fit in kMaxColours run one block-Jacobi round at the end. NullSync runs the identical body
// inline in list order - the threads=1 reference the digest tests compare against.
```

### `src/avbd/solver_step.cpp:221`

```cpp
    // Convergence early exit: measure the residual at the end of a round; when it
    // stops improving (delta below the threshold or rising), the remaining rounds would
    // burn cycles on float noise. 0 disables the check.
```

### `src/avbd/solver_step.cpp:229`

```cpp
        // Colour-parallel Gauss-Seidel. One pass per colour: assemble
        // (manifolds stamp self-contained against the live poses, other
        // forces stamp directly), solve, apply. Bodies of one colour
        // share no force, so the pass is race-free without a snapshot.
```

### `src/avbd/solver_step.cpp:256`

```cpp
        // Overflow bucket (> kMaxColours cliques): block Jacobi against
        // the post-colour snapshot, in three phases. Phase A precomputes
        // each manifold's shared per-contact assembly (world arms,
        // angular jacobians, clamped forces) once - work the two body
        // stamps used to duplicate. Nobody writes a pose until the
        // combine pass, so every phase reads the same frozen snapshot.
```

### `src/avbd/solver_step.cpp:277`

```cpp
            // Phase B: each bucket body assembles its full block system
            // (own inertia plus every incident constraint - manifolds
            // stamped from the phase-A cache, other forces stamping
            // directly) and solves for its delta.
```

### `src/avbd/solver_step.cpp:297`

```cpp
            // Combine: each bucket body applies its own delta, scaled by
            // jacobiRelaxation. A snapshot solve can overshoot where the
            // coloured Gauss-Seidel updates do not; omega < 1 trades
            // per-iteration progress for stability.
```

### `src/avbd/solver_step.cpp:310`

```cpp
        // Dual update: one item per force, and each force only reads the (frozen)
        // body poses and writes its own dual state, so the whole pass is independent.
        // Frozen constraints have nothing to update: every input they read is parked.
```

### `src/avbd/solver_step.cpp:333`

```cpp
        // Convergence probe after a Newton round: a full primal pass is expensive, so
        // skip further rounds when the residual stops shrinking. Only the calling
        // thread (worker 0) evaluates the reduction, and the verdict is stored in a
        // shared atomic that EVERY participant reads BEFORE the barrier - each worker
        // that sees `r_converged` skips its own remaining work but still participates
        // in every barrier (draining instead of returning), so the barrier
        // participation count stays equal and no worker deadlocks.
```

### `src/avbd/spring.cpp:1`

```cpp
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
```

### `src/avbd/types.hpp:1`

```cpp
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
```

### `src/avbd/types.hpp:29`

```cpp
// Base collision margin: a 1 mm skin keeps contacts from flickering on/off at the
// exact surface. It is deliberately small - the margin is added to the contact
// error, so a large value is a resting penetration floor.
```

### `src/avbd/types.hpp:37`

```cpp
// A sticky contact's anchors are material points frozen when static friction engaged, and
// because the two anchors stay glued to each other while the pair slides, their separation
// (what STICK_THRESH measures) never grows - stick never clears on its own, even as the
// frozen anchors slide metres off the touching surfaces and turn C0 into phantom depth.
// The merge heals instead: when a sticky pair's anchor world positions drift this far from
// where the current collision pass puts the contact, the anchors are re-adopted and the
// force history resets. The threshold sits at the collision-margin scale: below it the
// "drift" is numerical noise on a genuinely static contact and re-anchoring would churn
// resting manifolds for nothing.
```

### `src/avbd/types.hpp:47`

```cpp
// Geometric separation at which a contact's warm-started dual is dropped. The dual is a
// force, and the alpha*gamma warm-start decay (~0.99/tick) bleeds an impact impulse off
// over hundreds of ticks; while any of it is left the contact keeps applying it even
// though the pair has visibly parted - measured as a box popping out of an otherwise
// still pile at 10+ m/s seconds after the impact that loaded it. Anything above the
// collision-margin skin plus this tolerance is a real gap, and a real gap carries no
// contact force. Below it, the skin's own noise must not churn resting contacts.
```

### `src/avbd/types.hpp:56`

```cpp
// Cap on the penetration the constraint solve sees in one step (metres). A position-level
// step turns displacement into velocity at finishBody (v = dp/dt), so removing a deep
// impact overlap in a single step launches the pair at tens of m/s; that body then
// crosses other boxes before their next collision pass and the pile never settles. The
// cap bounds the per-step expulsion to 2 cm, i.e. <= 2 cm/dt, and the residual depth is
// resolved over the following steps.
```

### `src/avbd/types.hpp:69`

```cpp
// The collision and inertia shape of a body.
//
// `size` means different things per shape, and this documents it once:
//   Box      - full extents in x, y and z
//   Sphere   - size.x is the radius; y and z are ignored
//   Cylinder - size.x and size.y are the radius, size.z is the full height, and the axis is the
//              body's local +Z. (Godot's CylinderShape3D is Y-up; the server rotates it into
//              sim space, see godot_convert.hpp - the solver itself is always +Z.)
```

### `src/avbd/types.hpp:84`

```cpp
// A shape collision query: the collision and inertia shape of one body, decoupled from
// the body itself so shape pairs can be tested without a solver (space queries do this).
//
// `half` / `radius` / `halfHeight` mean different things per type, documented once:
//   Box      - half is half the full extents in x, y and z
//   Sphere   - radius; half and halfHeight unused
//   Cylinder - radius and halfHeight; the axis is the shape's local +Z, unit length,
//              expressed in the same space as `center` (the solver's cylinder is Z-axis)
```

### `src/avbd/types.hpp:124`

```cpp
    // Per-body gravity: the signed magnitude applied along -Z. Defaults to the
    // solver's global gravity; a host may stamp a different value per body (this
    // is how per-body gravity modes are expressed).
```

### `src/avbd/types.hpp:129`

```cpp
    // Consecutive steps with both velocities under the solver's sleep
    // thresholds; reaching sleepFrames flips `sleeping` on. Only meaningful when the body
    // is movable and sleep_mode enables sleeping.
```

### `src/avbd/types.hpp:134`

```cpp
    // Godot-style collision filtering: a pair collides when either side's layer is in the
    // other's mask. Defaults keep every pair colliding, which is what the solver did
    // before these existed.
```

### `src/avbd/types.hpp:140`

```cpp
    // Per-axis motion locks, bit 0..2 = solver x, y, z. A locked axis gets no gravity and
    // cannot carry velocity into a step, which is what keeps a locked body hanging in
    // place. Constraint coupling can still move it within a step; that residual is a
    // documented approximation.
```

### `src/avbd/types.hpp:147`

```cpp
    // Mass and inertia follow from the shape, so it is a constructor argument rather than
    // something set afterwards: a body that exists with the wrong inertia is a bug waiting to
    // happen. The five-argument form is a box.
```

### `src/avbd/types.hpp:162`

```cpp
    // Block-Jacobi update-set bookkeeping, rewritten by Solver::selectUpdateBodies() each
    // step: `inUpdateSet` marks the movable, awake bodies whose pose the primal phases may
    // change. Bodies outside the set (static, sleeping) are never solved or moved.
```

### `src/avbd/types.hpp:169`

```cpp
    // Position of this body in Solver::bodiesInOrder, rewritten by the broad phase every
    // step. The constraint adjacency and the pair-existence check address bodies through
    // it, so a rebuild never has to hash pointers.
```

### `src/avbd/types.hpp:173`

```cpp
    // True when every force on the body is a contact manifold. Only such bodies are
    // routed into the overflow bucket's block-Jacobi round by colourGraph: the contact
    // map is clamped, hard joints do not tolerate the snapshot solve.
```

### `src/avbd/types.hpp:184`

```cpp
// The 6x6 Newton system a single body solves each iteration, split into its two 3x3
// diagonal blocks plus the coupling block, and the matching right-hand side.
//
// This is the output of the primal phase: every constraint the body takes part in stamps
// its gradient and Hessian into one of these, and then it is solved for the update. Passing
// it as one object replaced six separate reference parameters on `updatePrimal`.
```

### `src/avbd/types.hpp:198`

```cpp
    // Mass term of the system: M / dt^2, with the body's current error as the right-hand
    // side. `dtSq` is dt*dt: the division is kept as a division on purpose, since
    // multiplying by a reciprocal would not be bit-identical.
    //
    // The angular half is the *world* inertia R I R^T, not the body-frame diagonal I. The
    // rest of the system is world-frame: the constraint angle Jacobians are `r x n` with r
    // in world space (ManifoldStore::stampCached builds R S R^T, Joint stamps `r x jLin`),
    // and the update `positionAng + dAng` is a world-frame rotation vector. Pairing that with
    // a body-frame diagonal inertia makes the effective angular mass wrong by whatever
    // rotation the body happens to have: only an isotropic inertia (a sphere, a cube) is
    // unaffected, and the demo's own 1 x 0.5 x 1 box is not (0.104 / 0.167 / 0.104).
    // projectContacts' angularCompliance/angularResponse already did this correctly.
```

### `src/avbd/types.hpp:227`

```cpp
    // Solve for the update without applying it (block-Jacobi path: each body stores its
    // delta and a later combine pass applies it). Assembles the 6x6 directly from the
    // block's members - no by-value matrix copies on the hot path.
```

### `src/nodes/avbd_soft_body3d.cpp:1`

```cpp
/*
 * AVBDSoftBody3D - a lattice of jointed boxes that behaves like a soft body.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 */
```

### `src/nodes/avbd_soft_body3d.cpp:28`

```cpp
// Offset of `p_distance` along a Godot axis, expressed in solver space. The
// conversion is the same one positions use, because the mapping flips one sign:
// Godot +Z points along solver -Y, so the neighbour offsets are not symmetric in
// the two frames.
```

### `src/nodes/avbd_soft_body3d.cpp:208`

```cpp
    // Rigid links along every axis, plus collision suppression for the diagonal
    // neighbours that the links do not already cover (constraint-connected pairs
    // never collide, see Rigid::constrainedTo).
```

### `src/nodes/avbd_soft_body3d.hpp:1`

```cpp
/*
 * AVBDSoftBody3D - a lattice of jointed boxes that behaves like a soft body.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * This is the formulation the AVBD authors use for soft bodies: a W x D x H grid of
 * rigid boxes wired together with stiff joints, with diagonal neighbours whose
 * collisions are suppressed. The boxes are real solver bodies, so the lattice
 * collides with the rest of the scene (and with itself along its outer shell). The
 * lattice is drawn through a single MultiMesh, so a 4x4x4 blob costs one draw call
 * and one node.
 *
 * The lattice is centred on this node's origin: the node's transform places it.
 */
```

### `src/nodes/avbd_soft_body3d.hpp:43`

```cpp
    // Latest simulated pose of each cell. A rebuild recreates the solver, so the poses are
    // kept here to put the lattice back where it was rather than resetting it to its spawn
    // layout. Absent until the lattice has been built once.
```

### `src/nodes/avbd_soft_world3d.cpp:1`

```cpp
/*
 * AVBDSoftWorld3D implementation. See the header.
 *
 * The rebuild / step / sync pattern follows the soft-body branches of the retired
 * AVBDWorld3D: snapshot each lattice's cell poses before the solver drops them so a
 * rebuild does not reset the scene, rebuild only when the node set changes, and push
 * cell poses back through _sync_visuals every simulated frame.
 */
```

### `src/nodes/avbd_soft_world3d.cpp:85`

```cpp
	// Scan children directly: shape-owner registration is timing-sensitive during
	// scene instantiation, and the ground contract is "a CollisionShape3D child
	// carrying one BoxShape3D".
```

### `src/nodes/avbd_soft_world3d.cpp:153`

```cpp
    // Solver gravity is the signed magnitude along the simulation's -Z; the exposed
    // property is a positive magnitude. Set BEFORE the rebuild: new rigids snapshot
    // the solver's gravity at construction.
```

### `src/nodes/avbd_soft_world3d.hpp:1`

```cpp
/*
 * AVBDSoftWorld3D - drives AVBDSoftBody3D nodes through a solver of their own.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * This is the minimal surviving piece of the old AVBDWorld3D node layer: soft bodies
 * have no PhysicsServer3D representation this round, so a node has to own their solver.
 * Rigid bodies and joints do NOT go here - a scene with those uses the AVBD physics
 * server (project setting physics/3d/physics_engine = "AVBD") and standard Godot nodes.
 *
 * Responsibilities, deliberately narrow:
 *   - collect the AVBDSoftBody3D nodes in this subtree (re-scan on tree changes)
 *   - collect child StaticBody3D "ground" nodes (one box shape each) as static solver
 *     bodies, so a lattice has something to rest on: soft cells live in this node's
 *     solver and cannot collide with bodies the physics server owns
 *   - step their solver in _physics_process (following the game: a paused SceneTree
 *     pauses this node through the engine's own pause handling)
 *   - write results back to the soft body visuals
 * Properties: gravity and threads only. Everything else is the soft body's own.
 */
```

### `src/nodes/godot_convert.hpp:1`

```cpp
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
```

### `src/nodes/godot_convert.hpp:44`

```cpp
// A direction expressed on a Godot axis, in solver space. Anything that names an *axis index*
// rather than passing a vector - per-axis joint modes, limits and springs - has to go through
// this, or it will drive the wrong axis: Godot Y is solver Z, and Godot Z is solver -Y.
```

### `src/nodes/godot_convert.hpp:109`

```cpp
// Property value to solver stiffness: strictly negative means "infinitely stiff"
// (a hard constraint). Godot floats and scene files cannot carry infinity, so the
// sentinel keeps the property serialisable. Zero disables the constraint.
```

### `src/nodes/register_types.cpp:1`

```cpp
/*
 * AVBD (Augmented Vertex Block Descent) for Godot 4 - GDExtension entry point.
 *
 * Part of godot-avbd. The solver core in src/avbd is ported from the authors'
 * reference implementation (github.com/savant117/avbd-demo3d, MIT).
 *
 * Rigid bodies, joints, areas and space queries run through the AVBD physics server
 * (standard Godot nodes, physics/3d/physics_engine = "AVBD"). The only node classes
 * left here are the soft-body pair: AVBDSoftBody3D (no server representation yet)
 * and AVBDSoftWorld3D, which owns their solver.
 */
```

### `src/nodes/register_types.cpp:29`

```cpp
// Offer AVBD to the engine as a selectable 3D physics engine. This is what lets a scene built
// from Godot's own nodes - RigidBody3D, CollisionShape3D and so on - be simulated by AVBD, without
// the project adopting any AVBD node classes. Selecting it is a project setting:
//
//   [physics]
//   3d/physics_engine="AVBD"
// The physics server factory the manager's creation Callable points at; freed at
// deinitialization because nothing else owns it.
```

### `src/nodes/register_types.cpp:45`

```cpp
    // The manager wants a callable that produces the server, so a small factory object holds one.
    // The manager only keeps the Callable, which does not own the object - the factory is
    // freed at deinitialization (see uninitialize_avbd_module).
```

### `src/nodes/register_types.cpp:53`

```cpp
// Register every physics/avbd/* solver setting with ProjectSettings, so the editor's
// Project Settings dialog shows an "AVBD" group instead of the settings existing only
// after a script writes them. Defaults must match read_project_params in the server.
```

### `src/nodes/register_types.cpp:106`

```cpp
// Everything registers at the SERVERS level (as godot-jolt does): the physics server
// object must keep its ClassDB entry until PhysicsServer3D::finish() has called
// _finish on it at exit. Registered at SCENE level the classes deregister with the
// main loop first, and teardown fails with "Required virtual method ::_finish must be
// overridden" plus a leaked instance.
```

### `src/nodes/register_types.cpp:122`

```cpp
    // SERVERS: the engine-facing classes (the physics server must exist before the
    // engine picks its physics engine, and its ClassDB entry must outlive
    // PhysicsServer3D::finish()).
```

### `src/server/avbd_direct_body_state3d.cpp:1`

```cpp
/*
 * AVBDDirectBodyState3D implementation. See the header for the contract.
 *
 * Frame conversion: the solver is Z-up, Godot is Y-up; to_sim/to_godot in
 * src/nodes/godot_convert.hpp do the mapping for vectors, quaternions and
 * transforms. Everything the engine reads comes straight off the bound
 * avbd::Rigid; writes go straight back onto it.
 *
 * Contact data is filled by the server after each step (frame_contacts on the
 * BodyData); until the server binds a body, this object reports neutral values
 * rather than dereferencing a null Rigid.
 */
```

### `src/server/avbd_direct_body_state3d.hpp:1`

```cpp
/*
 * AVBDDirectBodyState3D - the PhysicsDirectBodyState3D the server hands the engine.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * Godot's body state sync callback contract (verified against Godot 4.7,
 * scene/3d/physics/rigid_body_3d.cpp `_body_state_changed` and
 * modules/godot_physics_3d/godot_body_3d.cpp `call_queries`) is a *single argument*:
 * a PhysicsDirectBodyState3D. The same instance is also the first argument of the
 * force integration callback, followed by the stored userdata.
 *
 * One instance per space is cached by the server and re-pointed at each body before
 * the callbacks run, so the engine never holds a stale view. Every getter reads the
 * avbd::Rigid live; setters write straight through to it, mirroring what the engine's
 * own `_integrate_forces` users expect.
 *
 * The converter lives in src/nodes/godot_convert.hpp and maps Godot's Y-up scene
 * space onto the solver's Z-up space.
 */
```

### `src/server/avbd_direct_space_state3d.cpp:1`

```cpp
/*
 * AVBDDirectSpaceState3D implementation. See the header.
 *
 * All coordinates arriving here are Godot world space; the queries convert to the
 * solver's Z-up frame, walk the space's bodies and answer from collideShapes() - the
 * same pair test the solver's broad phase uses, so query results and contact behaviour
 * cannot disagree.
 *
 * The engine reaches this through World3D.direct_space_state (RayCast3D, intersect_*,
 * ShapeCast3D, move_and_collide's rest info, ...). Queries are answered from the solver's
 * current state, i.e. the poses as of the last step, exactly like the built-in server.
 */
```

### `src/server/avbd_direct_space_state3d.cpp:52`

```cpp
// // World-space (solver frame) position of body A's contact anchor.
// avbd::float3 contact_world_a(const avbd::Rigid &a, const avbd::Manifold::Contact &c) {
//     return a.positionLin + qrot(a.positionAng, c.rA);
// }
```

### `src/server/avbd_direct_space_state3d.cpp:59`

```cpp
    // Row 0 of the basis is the contact normal pointing from B to A; Godot's convention
    // (PhysicsServer3DExtensionRayResult.normal) is the surface normal facing the query,
    // which is from B towards A here as well - it is what the built-in server reports.
```

### `src/server/avbd_direct_space_state3d.hpp:1`

```cpp
/*
 * AVBDDirectSpaceState3D - the PhysicsDirectSpaceState3D behind space queries.
 *
 * Part of godot-avbd. Solver core ported from github.com/savant117/avbd-demo3d (MIT).
 *
 * The engine reaches this through World3D.direct_space_state (RayCast3D, intersect_*,
 * ShapeCast3D, ...). One instance per space, created lazily by the server.
 *
 * Queries walk the space's solver bodies and areas, coarse-filter by the query AABB and
 * answer with the same collideShapes() the solver uses, so query results and contact
 * behaviour cannot disagree. Exclusion lists (`params.exclude`) are honoured through
 * PhysicsDirectSpaceState3DExtension::is_body_excluded_from_query(); `ray_pickable` is
 * read for body and area hits. `cast_motion` sweeps a shape instead of testing only its
 * final pose; `_body_test_motion` uses the same sweep for CharacterBody3D/test_move.
 */
```

### `src/server/avbd_physics_server3d.cpp:1`

```cpp
/*
 * AVBD as a Godot physics server. See the header for the design and the unsupported list.
 */
```

### `src/server/avbd_physics_server3d.cpp:39`

```cpp
// --- step profiling (AVBD_PROFILE=1) -----------------------------------------
//
// Env-gated aggregate timers for headless runs that need to attribute per-tick
// cost to a phase: rebuild, project-settings reads, force integration, the
// solver itself, contact collection, sync dispatch. Off by default; the on
// path is a handful of steady_clock reads per step.
```

### `src/server/avbd_physics_server3d.cpp:64`

```cpp
    // Spike attribution: per-phase maxima within the reporting window. The average
    // hides exactly what a player feels - the one tick in a hundred that stalls the
    // frame (a deep-impact tick, a broad-phase pair storm, a rebuild).
```

### `src/server/avbd_physics_server3d.cpp:92`

```cpp
// Godot 4.3 added per-axis springs to the generic 6-DOF joint and grew the server enums; the
// godot-cpp headers this extension builds against predate that, so the new slots are named
// here with the engine's values (servers/physics_server_3d.h, G6DOFJointAxisParam/Flag).
```

### `src/server/avbd_physics_server3d.cpp:106`

```cpp
// -----------------------------------------------------------------------------
// Handles
//
// The extension interface offers no way to allocate a RID, so the server mints its own: an
// integer id written into the RID's opaque bytes through the pointer godot-cpp does expose. The
// engine only ever hands these values back to us, so the encoding is ours to choose - but it must
// round-trip, which is what id_of() checks by reading them back.
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:147`

```cpp
// -----------------------------------------------------------------------------
// Lifecycle
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:213`

```cpp
    // The engine flushes before stepping, so this is where the last step's poses reach
    // the scene: the per-body state-sync callable, then scripted integration. Both read
    // and write the solver's live state, so forces a script applies in _integrate_forces
    // land in the step about to run (the GodotPhysics contract, and godot-jolt's shape).
    // The first flush of a scene has no step behind it yet: nothing to deliver.
```

### `src/server/avbd_physics_server3d.cpp:331`

```cpp
// -----------------------------------------------------------------------------
// Spaces
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:347`

```cpp
// Solver parameters come from project settings the same way gravity does, so a project can
// tune the solver without a node layer. Keys absent from a project fall back to the solver's
// defaults.
```

### `src/server/avbd_physics_server3d.cpp:408`

```cpp
    // The space parameters are contact tolerances and sleep thresholds, none of which map onto
    // AVBD's model. They are accepted and ignored rather than rejected, so a project's settings do
    // not error out on load. Gravity does not come through here at all: Godot keeps it in the
    // project settings, and a server is expected to read it (see read_project_gravity).
```

### `src/server/avbd_physics_server3d.cpp:440`

```cpp
// One state object per space, created on first ask. PhysicsDirectSpaceState3D is
// reference-counted by the engine, so the server holds a raw pointer to the object it
// owns for the space's lifetime.
```

### `src/server/avbd_physics_server3d.cpp:500`

```cpp
// -----------------------------------------------------------------------------
// Shapes
//
// The engine sends a shape's parameters as a Variant from the matching *Shape3D resource, which
// is also what it can read back, so the getters rebuild the same object.
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:599`

```cpp
// -----------------------------------------------------------------------------
// Soft-body server API: documented unsupported. Godot only drives these through
// SoftBody3D nodes; every method answers neutrally and nothing is simulated.
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:772`

```cpp
    // Godot 4.x payload per shape type (scene/resources/3d/*_shape_3d.cpp::_update_shape):
    //   box      Vector3 - HALF extents (the resource sends `size / 2`)
    //   sphere   float radius
    //   cylinder Dictionary { radius: float, height: float } - full height
```

### `src/server/avbd_physics_server3d.cpp:930`

```cpp
// Godot BodyAxis bits (linear X/Y/Z = 1/2/4, angular X/Y/Z = 8/16/32) to the solver's
// sim-axis bits (0..2 = sim x/y/z). Godot Y is sim Z, Godot Z is sim -Y; the lock is a
// directionless bit so the sign does not matter, only the permutation.
```

### `src/server/avbd_physics_server3d.cpp:956`

```cpp
// -----------------------------------------------------------------------------
// The solver
//
// A space's solver is rebuilt from its bodies whenever any of them changes: AVBD holds one shape
// per body and a flat body list, so there is nothing to incrementally patch. Bodies are few
// enough per space that rebuilding is cheap, and it keeps the solver a pure function of the scene
// rather than of the order the engine happened to call things in.
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:965`

```cpp
// Density the solver should see for a body: its mass divided by its first shape's volume. The
// solver derives mass and inertia from shape and density, so the density has to carry exactly
// the engine's mass - a body that never set one is 1 kg (BODY_PARAM_MASS's documented default),
// whatever its size. Mass 0 means a body with no volume to divide into, or one the engine
// explicitly made massless; both stay massless in the solver, which reads that as immovable.
```

### `src/server/avbd_physics_server3d.cpp:1018`

```cpp
    // Deterministic order: ids are minted in creation order, so sorting by id makes the
    // solver's body list - and with it the Gauss-Seidel update order - a function of the
    // scene, not of unordered_map bucket layout. Two identical scenes rebuild identically.
```

### `src/server/avbd_physics_server3d.cpp:1060`

```cpp
                // Godot: x and z are the radius, y the height, axis local Y. The solver reads
                // the radius from size.x and the height from size.z (its cylinder axis is its
                // local Z, which is where Godot's Y lands in sim space), so order {r, r, h}.
```

### `src/server/avbd_physics_server3d.cpp:1092`

```cpp
            // No neutral proxy exists (a trimesh ground's AABB would be a giant solid
            // box): the body stays out of the solver entirely, per the unsupported-list
            // contract in the header. Swap the scene to a box shape instead.
```

### `src/server/avbd_physics_server3d.cpp:1113`

```cpp
        //body.rigid->sleeping = false;
        // Solver sleep modes: 0 = never sleeps, 1 = may sleep. Idle detection lives in the
        // solver (Solver::sleepFrames); can_sleep only opens the gate.
```

### `src/server/avbd_physics_server3d.cpp:1139`

```cpp
    // Solver-body -> server id index for the query paths: the BVH hands candidates back
    // as Rigid*, and this is the O(1) bridge to RIDs. Rebuilt here because every Rigid
    // was just re-created.
```

### `src/server/avbd_physics_server3d.cpp:1149`

```cpp
    // Joints second, once every body exists. One solver GenericJoint per engine joint
    // whose both bodies live in this space (or whose A side is the world). Every Godot
    // joint type maps onto the solver's GenericJoint; the specialised types fix which
    // axes are locked (conventions verified against Godot 4.7 sources):
    //   Pin       - all three linear axes locked, all angular free
    //   Hinge     - linear locked; angular free about the hinge axis, locked elsewhere
    //               (the hinge axis is the Z column of the hinge frame, see
    //               modules/godot_physics_3d/joints/godot_hinge_joint_3d.cpp:160)
    //   Slider    - free (or limited) along the slide axis (the X column of the joint
    //               frame, godot_slider_joint_3d.cpp:125), locked on the orthogonal
    //               linear axes; angular free about the slide axis, locked elsewhere
    //   ConeTwist - linear locked; angular free about the twist axis, swing limited
    //               on the two orthogonal axes (godot_cone_twist_joint_3d.cpp:117-139)
    // A joint axis that does not line up with a solver axis snaps to the nearest one:
    // the GenericJoint's degrees of freedom are principal axes only.
```

### `src/server/avbd_physics_server3d.cpp:1246`

```cpp
    // Anchors: Godot's local_ref is a full transform whose basis carries the joint frame;
    // AVBD wants the offset in each body's frame and measures angles from the relative
    // orientation at creation, which is what the two bases encode. With no A body the A
    // anchor is world-space already (the server's sim space is the world).
```

### `src/server/avbd_physics_server3d.cpp:1256`

```cpp
    // Six axes, Godot order X..Z linear then angular. Each Godot axis maps to a solver axis
    // with a sign (Godot X -> solver X +, Y -> Z +, Z -> Y -); the same kAxisMapping the
    // node layer uses, so both frontends mean the same joint.
```

### `src/server/avbd_physics_server3d.cpp:1328`

```cpp
    // The hinge axis: for _joint_make_hinge_simple it is the vector the engine stored in
    // ref_a's basis Z column (see _joint_make_hinge_simple below); for the full form it is
    // the Z column of the hinge frame. Expressed in body A's local frame.
```

### `src/server/avbd_physics_server3d.cpp:1455`

```cpp
// Run every active space for one step, then tell each body what happened.
//
// Stepping follows the game: while the engine has the server inactive (SceneTree paused)
// nothing steps. Solver parameters are stamped from project settings every step, and each
// space runs `substeps` solver steps per engine step with dt split across them.
```

### `src/server/avbd_physics_server3d.cpp:1494`

```cpp
        // A body joined/left this space since the last step; rebuild now that the
        // engine has finished sending its setup (shapes, mass, material). Runs before
        // the empty check: a space may be pending precisely because its first body is
        // only now arriving.
```

### `src/server/avbd_physics_server3d.cpp:1535`

```cpp
        // Godot's physics tick is the sampling rate. If an internal rate is configured,
        // accumulate engine time and advance the solver on its own fixed tick; the engine
        // tick then only samples the state. 0 keeps the legacy one-solver-tick-per-engine-
        // tick behaviour. The catch-up cap drops a backlog instead of spiralling.
```

### `src/server/avbd_physics_server3d.cpp:1554`

```cpp
        // Idle-sleep horizon: Godot's 0.5 s default scaled by sleep_time, expressed in
        // internal ticks. A fixed step count would make a high-tick-rate project fall
        // asleep sooner than intended - the horizon is a TIME, not 60 frames.
```

### `src/server/avbd_physics_server3d.cpp:1570`

```cpp
            // Forces accumulated since the last solver tick (persistent constant_* plus
            // one-shot frame ones) become velocity changes now, so the solver integrates
            // them with everything else. With an internal tick, a one-shot force waits
            // until the next internal tick rather than being consumed by a sample.
```

### `src/server/avbd_physics_server3d.cpp:1616`

```cpp
        // Snapshot the deepest contact on the stepping thread, where the force list is
        // stable. The getters run on the main thread; walking solver.forces there races
        // with manifold creation, revival and reaping, and a manifold read after it was
        // unlinked (or reused from the pool) reports a depth for a pair that is no longer
        // touching - the "phantom depth" that does not match any geometry.
        //
        // The stored error alone is not the answer either: it belongs to the contact set
        // the tick's collision pass found, and the bodies have moved since (a fast pair can
        // have separated within the tick, and a contact that rotated makes its stored
        // normal measure a separation that no longer exists). Report the pair the solver
        // believes is deepest, but measure it with a real SAT pass at the current pose, so
        // the number always describes the geometry the host can see. One extra narrow
        // phase per step, on the stepping thread.
```

### `src/server/avbd_physics_server3d.cpp:1675`

```cpp
        // Collect each body's contacts for the step (used by the DirectBodyState contact
        // getters, i.e. RigidBody3D contact monitoring). Manifolds are the forces with
        // contact points; a contact's world position comes from the anchor of the body
        // being filled. Capped by max_contacts_reported, as Godot does: without the cap
        // request no contacts are recorded. The manifold walk (and its two hash maps)
        // only runs when some body in this space actually reports contacts.
```

### `src/server/avbd_physics_server3d.cpp:1746`

```cpp
            // Snapshot the solver's output into the host record; _flush_queries delivers
            // it to the scene before the next step, through the callback the engine
            // registered for this body. The server does not touch the scene itself.
```

### `src/server/avbd_physics_server3d.cpp:1758`

```cpp
            // Keep the host copy in step with the solver's sleep state: the pre-step sync
            // writes body.sleeping back into body.rigid->sleeping, so a solver-side sleep
            // decision would otherwise be overwritten on the next step.
```

### `src/server/avbd_physics_server3d.cpp:1831`

```cpp
// -----------------------------------------------------------------------------
// Joints
//
// The engine hands over a joint RID, then a type (joint_make_*), then axis parameters and flags
// one call at a time. Everything is recorded in JointData; the solver object is built only in
// rebuild_space, so a joint whose bodies are still being set up costs nothing.
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:1883`

```cpp
    // A solver joint keeps its bodies from colliding no matter what (they share a force), so
    // Godot's default exclusion comes free; asking for contact between the pair cannot be
    // honoured and the record only reflects the request.
```

### `src/server/avbd_physics_server3d.cpp:2217`

```cpp
// -----------------------------------------------------------------------------
// Bodies
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:2238`

```cpp
    // The engine calls this the moment a node enters the tree - before its shapes and
    // physics material arrive. Rebuilding now would freeze an empty body into the
    // solver; deferring to the next step lets the whole setup land first.
```

### `src/server/avbd_physics_server3d.cpp:2385`

```cpp
// --- DirectBodyState accessors -------------------------------------------------
//
// The state object is the engine's view of one body during (or after) a step; these
// helpers do the BodyData lookups it cannot do itself. All of them tolerate an unknown
// RID: the engine may ask about a body this server never saw.
```

### `src/server/avbd_physics_server3d.cpp:2509`

```cpp
// The deepest normal overlap across every live manifold, as snapshotted at the end of the
// last step. A contact's C0 is stamped at the start of the step (Manifold::initialize), so
// this is the geometric penetration the last solve was working against, not a
// solver-internal residual.
```

### `src/server/avbd_physics_server3d.cpp:2527`

```cpp
// Sleep/wake state of the last step, aggregated over every space. Mirrors the solver's own
// predicates: a body is asleep when `sleeping` is set, a manifold is frozen when neither side
// is a movable awake body (ManifoldStore::frozen), and the update set is what the solve
// actually moved.
```

### `src/server/avbd_physics_server3d.cpp:2755`

```cpp
// -----------------------------------------------------------------------------
// Body parameters and forces
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:2928`

```cpp
// Godot BodyAxis bits (linear X/Y/Z = 1/2/4, angular X/Y/Z = 8/16/32) to the solver's
// sim-axis bits (0..2 = sim x/y/z). Godot Y is sim Z, Godot Z is sim -Y; the lock is a
// directionless bit so the sign does not matter, only the permutation.
```

### `src/server/avbd_physics_server3d.cpp:3148`

```cpp
    // Godot's contract for this call is "restore the default inertia and centre of mass to
    // cancel custom values". Both are always derived from the shape here, so the rebuild is
    // the whole of the reset - it must not touch the engine's mass.
```

### `src/server/avbd_physics_server3d.cpp:3211`

```cpp
    // The engine keeps the shape RIDs; this server only records the parameters it was given, so
    // there is no handle to give back. Returning an invalid RID is honest - the caller is asking
    // for something this implementation does not track.
```

### `src/server/avbd_physics_server3d.cpp:3250`

```cpp
// -----------------------------------------------------------------------------
// Areas
//
// Every storage virtual keeps the engine's description in an AreaData record, and
// _monitor_areas diffs each area against the bodies of its space every step, calling
// the monitor callbacks with Godot's four-argument contract (verified against 4.7
// scene/3d/physics/area_3d.cpp): (status, body_rid, instance_id, body_shape, area_shape)
// - status 0 = entered (AREA_BODY_ADDED), 1 = exited (AREA_BODY_REMOVED).
// -----------------------------------------------------------------------------
```

### `src/server/avbd_physics_server3d.cpp:3448`

```cpp
// One step of area monitoring: diff every area's inside-set against the previous one and
// fire the callbacks on the transitions. Godot's body-shape/area-shape arguments are the
// first overlapping pair (this server keeps one solver shape per body, so 0 on the body
// side and the area's first enabled shape on its side).
```

### `src/server/avbd_physics_server3d.hpp:1`

```cpp
/*
 * AVBD as a Godot physics server.
 *
 * Part of godot-avbd. This implements PhysicsServer3DExtension, so a scene built from Godot's own
 * nodes - RigidBody3D, StaticBody3D, CollisionShape3D, Generic6DOFJoint3D - is simulated by AVBD
 * instead of Godot's built-in engine. Project settings > Physics > 3D > Physics Engine selects it
 * by name once the extension is loaded.
 *
 * This is the main path of the extension: standard Godot nodes, no AVBD node classes. The only
 * node classes left in src/nodes are the soft-body pair (AVBDSoftBody3D + AVBDSoftWorld3D), which
 * have no server representation this round.
 *
 * Design
 * ------
 * - RIDs are minted here. The extension interface has no RID factory, so each handle is an
 *   integer id written into the RID's opaque bytes (`RID::_native_ptr()`); the engine only ever
 *   hands them back to us, so their internal representation is ours to choose. Every id maps to
 *   one of the records below, and the map is the single place an object's lifetime is decided.
 * - One Solver per space. Bodies are `avbd::Rigid`; the shapes attached to a body become that
 *   body's size and shape type, since AVBD gives each body exactly one collision shape. A body
 *   with several shapes uses the first, and says so rather than pretending.
 * - Transforms flow one way per step: the engine writes a body's transform through
 *   `_body_set_state`, `_step` runs the solver, and the result is pushed back by calling the sync
 *   callback the engine registered for that body with a single argument - a
 *   PhysicsDirectBodyState3D (the Godot 4.x contract, verified against 4.7 sources).
 * - Stepping follows the game: `_set_active(false)` (what SceneTree pause triggers) stops the
 *   simulation until the game resumes.
 *
 * Unsupported this round (README has the full list and the reasons): exact capsule/convex/
 * concave geometry (queries use the same cylinder/AABB proxies the solver simulates), the
 * soft-body server API, CCD, damping, automatic sleeping, area gravity overrides. Motion
 * queries (`body_test_motion`, `move_and_collide`, `test_move`) and space queries are
 * implemented through the shared `src/avbd/query/` extension.
 */
```

### `src/server/avbd_physics_server3d.hpp:71`

```cpp
    // --- records -----------------------------------------------------------------
    //
    // One shape, as the engine described it. AVBD bodies have a single shape, so a body's shape
    // list is collapsed to one of these when the body is simulated. Public because the space
    // queries read query shape descriptors through it.
```

### `src/server/avbd_physics_server3d.hpp:87`

```cpp
    // --- query shape conversion -------------------------------------------------
    // One place decides how a Godot shape becomes a solver query shape. Bodies and
    // areas use the same proxy the solver simulates: capsule -> cylinder, convex ->
    // the hull's local AABB box, unsupported -> false. Keeping this shared is what
    // keeps ray/shape-cast results consistent with contact behaviour.
```

### `src/server/avbd_physics_server3d.hpp:109`

```cpp
        // Body parameters the engine sets through _body_set_param. Mass decides the density the
        // solver is built with; friction is passed through; bounce has no AVBD counterpart (the
        // contact model is friction-only) and is remembered for the getter only. The default is
        // BODY_PARAM_MASS's documented 1 kg: a RigidBody3D whose scene never mentions `mass`
        // never calls the setter, and the solver must still simulate the engine's 1 kg.
```

### `src/server/avbd_physics_server3d.hpp:176`

```cpp
    // One joint, as the engine described it. Parameters arrive axis by axis after the joint is
    // made, so the record keeps everything and the solver is rebuilt from it. Godot axis order
    // (X, Y, Z) is preserved here; the solver permutation happens at rebuild time, once.
```

### `src/server/avbd_physics_server3d.hpp:193`

```cpp
        // [axis][param/flag], axis 0..2 = Godot X..Z, first three rows linear, last three
        // angular. Param slot count covers Godot 4.3+'s G6DOFJointAxisParam (spring entries at
        // 7..9 and 19..21); flag slots cover its G6DOFJointAxisFlag (spring toggles at 2 and 3).
```

### `src/server/avbd_physics_server3d.hpp:212`

```cpp
        // Re-read the project-setting solver parameters on the next step. Set at creation
        // and by the settings-changed signal, so the common tick does not hash nine
        // setting names for nothing.
```

### `src/server/avbd_physics_server3d.hpp:231`

```cpp
        // Rotation drift (unit-column distance, ~2 sin(angle/2)) above which a manifold's
        // cached angular jacobians are rebuilt. 0 disables the cache and reassembles them
        // for every contact on every iteration (original behaviour).
```

### `src/server/avbd_physics_server3d.hpp:239`

```cpp
        // Godot's physics tick is a sampling rate. When internal_ticks_per_second is
        // non-zero the solver advances on its own fixed tick and the engine tick only
        // samples the state; 0 keeps the legacy one-solver-tick-per-engine-tick.
```

### `src/server/avbd_physics_server3d.hpp:245`

```cpp
        // Positional projection of overlapping contacts, run after the velocity finish
        // (see Solver::positionProjection). The penalty solve alone leaves deep overlaps
        // in a large pile, and BDF1 turns any later correction into velocity; this bounds
        // the overlap in position instead, so under-solving cannot become kinetic energy.
```

### `src/server/avbd_physics_server3d.hpp:255`

```cpp
        // Idle-sleep thresholds. The horizon (how long a body must stay under both
        // thresholds) is a TIME - Godot's 0.5 s - converted to internal ticks per step
        // (`sleepFrames = 0.5 s / dt`), NOT a fixed 60 frames: at 120 internal Hz that is
        // 60 internal ticks, at 240 Hz it is 120. Only the velocity thresholds and a
        // multiplier on the horizon are project settings; the core Solver defaults stay
        // as they were so pinned digests are unaffected.
```

### `src/server/avbd_physics_server3d.hpp:268`

```cpp
        // Route contact-only bodies through the coloured Gauss-Seidel instead of the
        // block-Jacobi bucket (Solver::contactOnlyOverflow). Off by default: see the A/B in
        // solver.h.
```

### `src/server/avbd_physics_server3d.hpp:273`

```cpp
        // Contact generation distance, in metres (see Solver::speculativeMargin). A contact
        // is created while the shapes are still separated by up to this much, and the
        // velocity pass clamps the approach so the gap closes within one step at most -
        // a fast body is stopped at the surface instead of first penetrating several
        // centimetres and paying for it with a stiff-penalty impulse. 0 = only real overlaps.
```

### `src/server/avbd_physics_server3d.hpp:290`

```cpp
        // Solver-body -> server id, rebuilt with the space. The direct space state's
        // queries get their candidates back as Rigid* from the BVH; this is the O(1)
        // bridge back to RIDs (body_id_of used to scan every body per result).
```

### `src/server/avbd_physics_server3d.hpp:320`

```cpp
    // Object kinds are kept in separate maps so an id that arrives for the wrong kind is a miss
    // rather than a reinterpretation. `ankerl::unordered_dense` stores values contiguously and
    // relocates them on rehash, which is fine for the plain-data records. SpaceData owns a
    // Solver - a non-movable RAII holder of intrusive lists whose nodes point back at it - so
    // spaces keeps it behind a stable unique_ptr: the flat map only ever moves the handle,
    // never the solver, and the space's address survives any rehash.
```

### `src/server/avbd_physics_server3d.hpp:342`

```cpp
    // Deepest contact of the last completed step, snapshotted on the stepping thread.
    // The getters run on the main thread, where walking solver.forces races with the
    // solver's manifold creation/revival/reaping and can report depths for pooled
    // manifolds that no longer describe a touching pair.
```

### `src/server/avbd_physics_server3d.hpp:370`

```cpp
    // Rebuild the solver for one space from the bodies and joints assigned to it. Needed
    // whenever a body, its shapes, its mode or a joint changes, because those decide what the
    // solver holds.
    // The source location is the rebuild reason: every call site is a distinct engine
    // entry point, and knowing which one fired (and when) is what tells a mid-run
    // rebuild - which destroys every manifold - apart from scene setup.
```

### `src/server/avbd_physics_server3d.hpp:384`

```cpp
    // Why the rebuilds happened, and when the last one landed relative to _step() calls. A
    // rebuild destroys every manifold, so one arriving mid-run silently discards all
    // warm-started contact state; these make that visible from a script.
```

### `src/server/avbd_physics_server3d.hpp:424`

```cpp
    // The same contract, narrowed by the last built broad-phase tree: only bodies whose
    // bound overlaps the query AABB (sim space). Before any step has built a tree this
    // answers from the body list with per-body bounds.
```

### `src/server/avbd_physics_server3d.hpp:465`

```cpp
    // Penetration diagnostics for the headless stress suites. A manifold's contact error
    // C0 carries the geometric overlap captured at the start of the step: its normal
    // component is `COLLISION_MARGIN - depth`, so the depth in metres is margin - C0.x.
    // Returns the deepest contact across every space (0 when nothing penetrates).
```

### `src/server/avbd_physics_server3d.hpp:473`

```cpp
    // Sleep/wake diagnostics for the stress suites. A sleeping body is skipped by the
    // velocity update in Solver::finishBody, and a manifold with neither side awake is
    // frozen (its contact state is deliberately left stale), so both counts say whether a
    // penetration measurement was taken against a live solve or a parked one. `rebuilds` is
    // cumulative over the server's lifetime: a rebuild destroys every manifold, so a rising
    // count during a run means warm-started contact state is being thrown away.
```

### `src/server/avbd_physics_server3d.hpp:560`

```cpp
    // --- soft-body server API: documented unsupported this round -----------------
    // The engine calls these only for SoftBody3D nodes, which the demo scenes do not
    // use (the AVBDSoftBody3D node class runs its own solver). Neutral answers keep a
    // stray SoftBody3D from erroring; nothing is simulated.
```

### `src/server/avbd_physics_server3d.hpp:660`

```cpp
    // --- joints -------------------------------------------------------------------
    // Only the generic 6-DOF joint gets a solver counterpart: AVBD's GenericJoint covers every
    // configuration the car and the demos need, and the other joint types are accepted and
    // tracked (so scenes load) without acting. Their setters record and their getters replay.
```

### `src/server/avbd_physics_server3d.hpp:714`

```cpp
    // --- areas --------------------------------------------------------------------
    // Created and tracked; every storage virtual is implemented and monitoring runs each
    // step (see _monitor_areas).
```

### `src/server/physics_server_factory.cpp:1`

```cpp
/*
 * Factory for the AVBD physics server. See the header.
 */
```

### `src/server/physics_server_factory.hpp:1`

```cpp
/*
 * The physics server manager expects a Callable that produces a server instance, so this is the
 * object that holds that method. It exists only for registration.
 */
```

