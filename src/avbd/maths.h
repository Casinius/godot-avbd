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
//   1. Hand-unrolled scalar LDL^T over the lower triangle - no matrix object, no
//      temporary, no lazy products: the ~70 flops a 6x6 needs, straight into registers.
//      Measured ~2.5x cheaper than the generic factorisation it replaced (the whole
//      step got 8-18% cheaper depending on body count), and it does not move a single
//      division, so its rounding is the same one this file used before Eigen.
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
    // (1) Hand-unrolled scalar LDL^T over the lower triangle: no matrix object, no
    // temporary, no lazy-product machinery - just the ~70 flops the 6x6 actually needs,
    // straight into registers. This kernel runs once per movable body per iteration, so
    // the generic Eigen factorisation below it costs more than the solve itself.
    float A11 = aLin(0, 0);
    float A21 = aLin(1, 0), A22 = aLin(1, 1);
    float A31 = aLin(2, 0), A32 = aLin(2, 1), A33 = aLin(2, 2);
    float A41 = aCross(0, 0), A42 = aCross(0, 1), A43 = aCross(0, 2), A44 = aAng(0, 0);
    float A51 = aCross(1, 0), A52 = aCross(1, 1), A53 = aCross(1, 2), A54 = aAng(1, 0), A55 = aAng(1, 1);
    float A61 = aCross(2, 0), A62 = aCross(2, 1), A63 = aCross(2, 2), A64 = aAng(2, 0), A65 = aAng(2, 1), A66 = aAng(2, 2);

    float L21 = A21 / A11;
    float L31 = A31 / A11;
    float L41 = A41 / A11;
    float L51 = A51 / A11;
    float L61 = A61 / A11;

    float D1 = A11;

    float D2 = A22 - L21 * L21 * D1;

    float L32 = (A32 - L21 * L31 * D1) / D2;
    float L42 = (A42 - L21 * L41 * D1) / D2;
    float L52 = (A52 - L21 * L51 * D1) / D2;
    float L62 = (A62 - L21 * L61 * D1) / D2;

    float D3 = A33 - (L31 * L31 * D1 + L32 * L32 * D2);

    float L43 = (A43 - L31 * L41 * D1 - L32 * L42 * D2) / D3;
    float L53 = (A53 - L31 * L51 * D1 - L32 * L52 * D2) / D3;
    float L63 = (A63 - L31 * L61 * D1 - L32 * L62 * D2) / D3;

    float D4 = A44 - (L41 * L41 * D1 + L42 * L42 * D2 + L43 * L43 * D3);

    float L54 = (A54 - L41 * L51 * D1 - L42 * L52 * D2 - L43 * L53 * D3) / D4;
    float L64 = (A64 - L41 * L61 * D1 - L42 * L62 * D2 - L43 * L63 * D3) / D4;

    float D5 = A55 - (L51 * L51 * D1 + L52 * L52 * D2 + L53 * L53 * D3 + L54 * L54 * D4);

    float L65 = (A65 - L51 * L61 * D1 - L52 * L62 * D2 - L53 * L63 * D3 - L54 * L64 * D4) / D5;

    float D6 = A66 - (L61 * L61 * D1 + L62 * L62 * D2 + L63 * L63 * D3 + L64 * L64 * D4 + L65 * L65 * D5);

    // Forward substitution Ly = b, then Dz = y, then back substitution L^T x = z, fused:
    // the scalars stay live, so the compiler never spills them to the stack frame.
    float y1 = bLin[0];
    float y2 = bLin[1] - L21 * y1;
    float y3 = bLin[2] - L31 * y1 - L32 * y2;
    float y4 = bAng[0] - L41 * y1 - L42 * y2 - L43 * y3;
    float y5 = bAng[1] - L51 * y1 - L52 * y2 - L53 * y3 - L54 * y4;
    float y6 = bAng[2] - L61 * y1 - L62 * y2 - L63 * y3 - L64 * y4 - L65 * y5;

    float z6 = y6 / D6;
    float xa2 = z6;
    float xa1 = y5 / D5 - L65 * xa2;
    float xa0 = y4 / D4 - L54 * xa1 - L64 * xa2;
    float xl2 = y3 / D3 - L43 * xa0 - L53 * xa1 - L63 * xa2;
    float xl1 = y2 / D2 - L32 * xl2 - L42 * xa0 - L52 * xa1 - L62 * xa2;
    float xl0 = y1 / D1 - L21 * xl1 - L31 * xl2 - L41 * xa0 - L51 * xa1 - L61 * xa2;

    const bool sane = std::isfinite(xl0) && std::isfinite(xl1) && std::isfinite(xl2) &&
            std::isfinite(xa0) && std::isfinite(xa1) && std::isfinite(xa2) &&
            std::fabs(xl0) < 1.0e1f && std::fabs(xl1) < 1.0e1f && std::fabs(xl2) < 1.0e1f &&
            std::fabs(xa0) < 1.0e1f && std::fabs(xa1) < 1.0e1f && std::fabs(xa2) < 1.0e1f;
    if (sane)
    {
      xLin[0] = xl0;
      xLin[1] = xl1;
      xLin[2] = xl2;
      xAng[0] = xa0;
      xAng[1] = xa1;
      xAng[2] = xa2;
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
