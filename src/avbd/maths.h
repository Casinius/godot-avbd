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
 * Port note: an earlier revision of this file contained `using namespace std;`,
 * which leaked every standard-library name into any translation unit that
 * included it. The handful of places that relied on it now name `std::`
 * explicitly.
 */

#pragma once

#include "Eigen/Core"
#include "Eigen/Geometry"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <array>
#include <cmath>
#include <cstddef>

namespace avbd {
using float3x3 = Eigen::Matrix<float, 3, 3>;
using float3 = Eigen::Vector3<float>;
using float2 = Eigen::Vector2<float>;
using quat = Eigen::Quaternion<float>;
[[nodiscard]] inline float3x3 diagonal(float m00, float m11, float m22) noexcept
{
    float3x3 m;
    m << m00, 0.0f, 0.0f,
         0.0f, m11, 0.0f,
         0.0f, 0.0f, m22;
    return m;
}

// Right-handed orthonormal frame with the given unit vector as its first ROW (the contact
// normal). Row-major on purpose: the manifold code treats basis.row(0) as the normal.
[[nodiscard]] inline float3x3 orthonormal(float3 normal) noexcept
{
    float3 t1 = std::fabs(normal.x()) > std::fabs(normal.z())
            ? float3{-normal.y(), normal.x(), 0}
            : float3{0, -normal.z(), normal.y()};
    t1 = t1.normalized();
    const float3 t2 = normal.cross(t1);
    float3x3 m;
    m.row(0) = normal;
    m.row(1) = t1;
    m.row(2) = t2;
    return m;
}
[[nodiscard]] inline float3 operator-(quat a, quat b) noexcept {
  const quat d = a * b.inverse();
  return float3{d.x(), d.y(), d.z()} * 2.0f;
}

// First-order quaternion integration, as the original custom maths did: half of omega (as a
// pure quaternion) left-multiplied, then renormalised. The added rotation vector is the same
// small-angle quantity `operator-(quat, quat)` measures.
[[nodiscard]] inline quat operator+(quat a, float3 b) noexcept {
  // Eigen's four-scalar Quaternion constructor is (w, x, y, z): a pure quaternion has w = 0.
  const quat omega(0.0f, b.x(), b.y(), b.z());
  quat result;
  result.coeffs() = a.coeffs() + (omega * a).coeffs() * 0.5f;
  return result.normalized();
}

// Solve the symmetric 6x6 system [aLin  aCross^T; aCross  aAng] x = [bLin;
// bAng].
//
// Two kernels, one bounded-update policy:
//   1. Eigen fixed-size LDLT - stack-only, ~200 ns kernel. The workhorse.
//   2. Eigen dynamic LDLT - thread_local scratch (no steady-state allocation), different
//      pivot rounding. Retry path when (1) is not sane.
//   3. Zero update - last resort when neither kernel is sane.
// "Sane" means finite and < 10 m per component: an update is velocity*dt plus a position
// correction, so 10 m is orders past anything physical. Degenerate blocks (thin cylinders
// standing on end produce near-singular stiffness ratios ~1e6) used to let one wild pivot
// compound into e15-scale positions; the policy caps every solve at the physical bound.
// The fixed kernel's rounding differs from the old dynamic-only implementation, which
// re-anchors every state digest - accepted on the real-time branch.
// Fixed-size 6x6 LDL^T without pivoting, over an already-assembled system. The
// system is dominated by the M/dt^2 diagonal (a 1 kg body at 240 Hz carries a
// 57,600 kg/s^2 diagonal against ~1e4 penalty terms), so the pivot growth Eigen's
// symmetric pivoting exists for does not occur here, and dropping it - plus Eigen's
// expression and Transpositions plumbing - lands the kernel an order of magnitude
// cheaper. Plain sequential float ops keep the result bit-identical for every thread
// count. Returns false on a degenerate pivot (non-positive D, which a near-singular
// block can produce) or an unsane solution - callers retry with the pivoted dynamic
// path in solve() below.
inline bool ldlt6_solve(const std::array<std::array<float, 6>, 6> &a, const std::array<float, 6> &b,
                        std::array<float, 6> &x) noexcept {
  std::array<std::array<float, 6>, 6> l{};
  std::array<float, 6> d{};
  for (int j = 0; j < 6; j++) {
    float dj = a[j][j];
    for (int s = 0; s < j; s++)
      dj -= l[j][s] * l[j][s] * d[s];
    if (!(dj > 1.0e-12f))
      return false;
    d[j] = dj;
    l[j][j] = 1.0f;
    for (int i = j + 1; i < 6; i++) {
      float lij = a[i][j];
      for (int s = 0; s < j; s++)
        lij -= l[i][s] * l[j][s] * d[s];
      l[i][j] = lij / d[j];
    }
  }

  std::array<float, 6> z{};
  // Forward solve L z = b (unit lower triangle).
  for (int i = 0; i < 6; i++) {
    float v = b[i];
    for (int s = 0; s < i; s++)
      v -= l[i][s] * z[s];
    z[i] = v;
  }
  // Diagonal scaling y = D^-1 z, then back solve L^T x = y, all in one array.
  for (int i = 0; i < 6; i++)
    x[i] = z[i] / d[i];
  for (int i = 5; i >= 0; i--) {
    float v = x[i];
    for (int s = i + 1; s < 6; s++)
      v -= l[s][i] * x[s];
    x[i] = v;
  }

  float worst = 0.0f;
  for (int i = 0; i < 6; i++)
    worst = std::max(worst, std::abs(x[i]));
  return std::isfinite(worst) && worst < 1.0e1f;
}

inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin,
                  float3 bAng, float3 &xLin, float3 &xAng) noexcept {
  {
    std::array<std::array<float, 6>, 6> a{};
    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
        a[r][c] = aLin(r, c);
        a[r + 3][c + 3] = aAng(r, c);
        a[r + 3][c] = aCross(r, c);     // [aLin aCross^T; aCross aAng]
        a[r][c + 3] = aCross(c, r);
      }
    }

    const std::array<float, 6> b{bLin[0], bLin[1], bLin[2], bAng[0], bAng[1], bAng[2]};
    std::array<float, 6> x{};
    if (ldlt6_solve(a, b, x)) {
      xLin = float3{x[0], x[1], x[2]};
      xAng = float3{x[3], x[4], x[5]};
      return;
    }
  }

  thread_local Eigen::MatrixXf A(6, 6);
  thread_local Eigen::MatrixXf b(6, 1);
  thread_local Eigen::MatrixXf x(6, 1);
  thread_local Eigen::LDLT<Eigen::MatrixXf> solver;

  A.setZero();
  A.block<3, 3>(0, 0) = aLin;
  A.block<3, 3>(3, 3) = aAng;
  A.block<3, 3>(0, 3) = aCross.transpose();
  A.block<3, 3>(3, 0) = aCross;
  b << bLin[0], bLin[1], bLin[2], bAng[0], bAng[1], bAng[2];

  solver.compute(A);
  x = solver.solve(b);
  if (x.allFinite() && x.cwiseAbs().maxCoeff() < 1.0e1f)
  {
    xLin[0] = x(0, 0);
    xLin[1] = x(1, 0);
    xLin[2] = x(2, 0);
    xAng[0] = x(3, 0);
    xAng[1] = x(4, 0);
    xAng[2] = x(5, 0);
    return;
  }

  // Degenerate block: freeze this DOF for one iteration instead of teleporting.
  xLin = float3{0, 0, 0};
  xAng = float3{0, 0, 0};
}

} // namespace avbd
