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
inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin,
                  float3 bAng, float3 &xLin, float3 &xAng) noexcept {
  {
    Eigen::Matrix<float, 6, 6> A;
    A << aLin, aCross.transpose(), aCross, aAng; // fills row-major by 3x3 blocks: [aLin aCross^T; aCross aAng]
    Eigen::Matrix<float, 6, 1> b;
    b << bLin, bAng;
    const Eigen::LDLT<Eigen::Matrix<float, 6, 6>> solver(A);
    const Eigen::Matrix<float, 6, 1> x = solver.solve(b);
    if (x.allFinite() && x.cwiseAbs().maxCoeff() < 1.0e1f)
    {
      xLin = x.head<3>();
      xAng = x.tail<3>();
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
