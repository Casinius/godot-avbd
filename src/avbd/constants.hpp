/*
 * avbd::constants — named physical and numerical constants used across the solver.
 *
 * Replaces raw literals (3.14159…, 9.8, 1e-4, etc.) with named, typed constexpr values
 * so every call site is self-documenting and there is a single place to adjust precision
 * or default tuning.
 *
 * Copyright (c) 2026 Chris Giles
 * SPDX-License-Identifier: MIT
 */

#pragma once

#ifndef AVBD_CONSTANTS_HPP
#define AVBD_CONSTANTS_HPP

#include <cmath>

namespace avbd {
namespace constants {

// --- Mathematical constants (float precision, physics-grade) ---

/// π as float.  Use this instead of typing out digits.
constexpr float pi_f = 3.14159265358979323846f;

/// π as double, for the rare double-precision volume calculations in the server layer.
constexpr double pi_d = 3.14159265358979323846;

/// (4/3)π  — sphere volume coefficient
constexpr float four_thirds_pi_f = (4.0f / 3.0f) * pi_f;
constexpr double four_thirds_pi_d = (4.0 / 3.0) * pi_d;

// --- Gravity defaults ---
// The canonical default is 9.8 m/s² (Earth surface).  The solver stores the *signed*
// value along −Z; consumers that need the sign apply it themselves.

/// Default gravitational acceleration magnitude (m/s²).
constexpr float default_gravity = 9.8f;

// --- Solver tolerances and limits ---

/// Maximum penetration error tolerance (metres).  Default 0.1 mm.
constexpr float max_penetration_error = 1e-4f;

/// BVH Morton-code quantisation: 2²¹−1 (21-bit unsigned maximum).
constexpr float morton_scale_divisor = 2097151.0f; // 2^21 - 1

/// Small epsilon to avoid division by zero in BVH extent normalization.
constexpr float bvh_extent_epsilon = 1e-6f;

/// Upper bound on a sane Newton-solve update (metres).  Updates exceeding
/// this are clamped to zero (degenerate block freeze).
constexpr float solve_sanity_bound = 1.0e1f;

} // namespace constants
} // namespace avbd

#endif // AVBD_CONSTANTS_HPP
