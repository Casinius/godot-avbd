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
// Three kernels, one bounded-update policy:
//   1. Schur complement on the linear half, 3x3 inverses - ~23 ns. The workhorse.
//   2. Hand-rolled 6x6 LDL^T without pivoting - ~115 ns. Retry when (1) is not sane.
//   3. Eigen dynamic LDLT - thread_local scratch (no steady-state allocation), different
//      pivot rounding. Retry path when (1) and (2) are not sane.
//   4. Zero update - last resort when none of them is.
// "Sane" means finite and < 10 m per component: an update is velocity*dt plus a position
// correction, so 10 m is orders past anything physical. Degenerate blocks (thin cylinders
// standing on end produce near-singular stiffness ratios ~1e6) used to let one wild pivot
// compound into e15-scale positions; the policy caps every solve at the physical bound.
// The kernel's rounding re-anchors every state digest whenever it changes - accepted on the
// real-time branch, which is why the cheapest kernel is the one that runs first.
inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin,
                  float3 bAng, float3 &xLin, float3 &xAng) noexcept {
  {
    // Block elimination on the linear half. That block is M/dt^2 plus linear stiffness gathered
    // from PSD stamps, so it is SPD and is by far the better conditioned half of the system -
    // eliminating it first costs two 3x3 inverses and four 3x3 products, against the 6x6
    // factorization's 56 pivots and two triangular sweeps. Measured on a representative block (a
    // 1 kg body at 240 Hz with four contacts): 23 ns against 115, agreeing with the LDL^T to 1e-7
    // relative across the stiffness range the scenes use.
    //
    // The guards are the same bounded-update policy the other kernels carry. A block whose
    // determinant collapsed - a body with a zero principal moment, or a stiffness ratio past what
    // float resolves - falls through to the factorizations instead of producing a wild update, and
    // the physical bound below is the final say either way.
    {
      const float3x3 aLinInv = aLin.inverse();
      const float3x3 aCrossT = aCross.transpose();
      const float3x3 schur = aAng - aCross * aLinInv * aCrossT;

      const float scaleLin = aLin.cwiseAbs().maxCoeff();
      const float scaleSchur = schur.cwiseAbs().maxCoeff();
      const float detLin = aLin.determinant();
      const float detSchur = schur.determinant();
      const bool invertible = std::isfinite(detLin) && std::isfinite(detSchur) &&
              std::fabs(detLin) > scaleLin * scaleLin * scaleLin * 1.0e-12f &&
              std::fabs(detSchur) > scaleSchur * scaleSchur * scaleSchur * 1.0e-12f;

      if (invertible) {
        const float3 candAng = schur.inverse() * (bAng - aCross * (aLinInv * bLin));
        const float3 candLin = aLinInv * (bLin - aCrossT * candAng);
        const float worst = std::max(candLin.cwiseAbs().maxCoeff(), candAng.cwiseAbs().maxCoeff());
        if (std::isfinite(worst) && worst < 1.0e1f) {
          xLin = candLin;
          xAng = candAng;
          return;
        }
      }
    }

    // Hand-rolled 6x6 LDL^T without pivoting. The system is dominated by the M/dt^2
    // Hand-rolled 6x6 LDL^T without pivoting. The system is dominated by the M/dt^2
    // diagonal (a 1 kg body at 240 Hz carries a 57,600 kg/s^2 diagonal against ~1e4
    // penalty terms), so the pivot growth Eigen's symmetric pivoting exists for does
    // not occur here, and dropping it - plus Eigen's expression and Transpositions
    // plumbing - lands the kernel an order of magnitude cheaper. Plain sequential
    // float ops keep the result bit-identical for every thread count. A degenerate
    // pivot (non-positive D, which a near-singular block can produce) bails to the
    // pivoted dynamic retry below.
    std::array<std::array<float, 6>, 6> a{};
    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++) {
        a[r][c] = aLin(r, c);
        a[r + 3][c + 3] = aAng(r, c);
        a[r + 3][c] = aCross(r, c);     // [aLin aCross^T; aCross aAng]
        a[r][c + 3] = aCross(c, r);
      }
    }

    std::array<std::array<float, 6>, 6> l{};
    std::array<float, 6> d{};
    bool factored = true;
    for (int j = 0; j < 6; j++) {
      float dj = a[j][j];
      for (int s = 0; s < j; s++)
        dj -= l[j][s] * l[j][s] * d[s];
      if (!(dj > 1.0e-12f)) {
        factored = false;
        break;
      }
      d[j] = dj;
      l[j][j] = 1.0f;
      for (int i = j + 1; i < 6; i++) {
        float lij = a[i][j];
        for (int s = 0; s < j; s++)
          lij -= l[i][s] * l[j][s] * d[s];
        l[i][j] = lij / d[j];
      }
    }

    if (factored) {
      std::array<float, 6> z{};
      // Forward solve L z = b (unit lower triangle).
      for (int i = 0; i < 6; i++) {
        float v = i < 3 ? bLin[i] : bAng[i - 3];
        for (int s = 0; s < i; s++)
          v -= l[i][s] * z[s];
        z[i] = v;
      }
      // Diagonal scaling y = D^-1 z, then back solve L^T x = y, all in one array.
      std::array<float, 6> x{};
      for (int i = 0; i < 6; i++)
        x[i] = z[i] / d[i];
      for (int i = 5; i >= 0; i--) {
        float v = x[i];
        for (int s = i + 1; s < 6; s++)
          v -= l[s][i] * x[s];
        x[i] = v;
      }

      const float3 candLin{x[0], x[1], x[2]};
      const float3 candAng{x[3], x[4], x[5]};
      const float worst = std::max(candLin.cwiseAbs().maxCoeff(), candAng.cwiseAbs().maxCoeff());
      if (std::isfinite(worst) && worst < 1.0e1f) {
        xLin = candLin;
        xAng = candAng;
        return;
      }
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
