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
 */
// Math over linalg.h: engine layout (float3x3 m[col][row], quat {x,y,z,w}); operator* is
// component-wise, use mul/dot/qmul/qrot/outerprod; transpose for rows (Eigen term order kept where noted).

#pragma once

#include "linalg.h"
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <limits>

namespace avbd {
using float2 = linalg::aliases::float2;
using float3 = linalg::aliases::float3;
using float3x3 = linalg::aliases::float3x3;
using quat = linalg::aliases::float4;

// Re-export linalg free functions so call sites read `dot(a,b)` unqualified; scalar std::min/max/
// clamp are always std::-qualified at call sites, so nothing collides.
using linalg::abs;
using linalg::cross;
using linalg::dot;
using linalg::length;
using linalg::length2;
using linalg::max;
using linalg::min;
using linalg::mul;
using linalg::normalize;
using linalg::outerprod;
using linalg::qconj;
using linalg::rotation_quat;
using linalg::transpose;

// Keep Eigen's (x^2+z^2)+(y^2+w^2) reduction order (linalg folds left-to-right; different
// rounding). Non-template overloads win for vec<float,4>; smaller vectors use linalg's matching order.
[[nodiscard]] inline float qlength2(quat q) noexcept
{
    return (q.x * q.x + q.z * q.z) + (q.y * q.y + q.w * q.w);
}
[[nodiscard]] inline float length(quat q) noexcept { return std::sqrt(qlength2(q)); }
[[nodiscard]] inline quat normalize(quat q) noexcept
{
    const float z = qlength2(q);
    return z > 0.0f ? q / std::sqrt(z) : q;
}
[[nodiscard]] inline quat qinv(quat q) noexcept
{
    return quat{-q.x, -q.y, -q.z, q.w} / qlength2(q);
}

// Use Eigen's operation order (t=2*(z x v); v+w*t+z x t, tx/ty/... factorisation): linalg's
// equivalent expressions round by differing ulps, and on stiff scenes (lattices, pyramids) drift compounds.
[[nodiscard]] inline float3 qrot(quat q, float3 v) noexcept
{
    const float3 z{q.x, q.y, q.z};
    float3 uv = cross(z, v);
    uv += uv; // x2 as a doubling, like Eigen's `uv += uv`
    return v + uv * q.w + cross(z, uv);
}

[[nodiscard]] inline float3x3 qmat(quat q) noexcept
{
    const float tx = 2.0f * q.x, ty = 2.0f * q.y, tz = 2.0f * q.z;
    const float twx = tx * q.w, twy = ty * q.w, twz = tz * q.w;
    const float txx = tx * q.x, txy = ty * q.x, txz = tz * q.x;
    const float tyy = ty * q.y, tyz = tz * q.y, tzz = tz * q.z;
    // Each braced float3 is a COLUMN. R12/R21 = (tyz - twx, tyz + twx): do not swap - that is
    // neither a rotation nor a transpose, and broad-phase bounds would exclude real overlaps.
    return float3x3{
            {1.0f - (tyy + tzz), txy + twz, txz - twy},
            {txy - twz, 1.0f - (txx + tzz), tyz + twx},
            {txz + twy, tyz - twx, 1.0f - (txx + tyy)}};
}

// Hamilton product with Eigen's exact term order per component (linalg's qmul adds
// the same four terms in a different order, which rounds differently).
[[nodiscard]] inline quat qmul(quat a, quat b) noexcept
{
    return quat{
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
            a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
            a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

// Component-wise helpers standing in for the Eigen cwise* family.
using linalg::cmul;

[[nodiscard]] inline float3x3 diagonal(float m00, float m11, float m22) noexcept
{
    return float3x3{float3{m00, 0.0f, 0.0f}, float3{0.0f, m11, 0.0f}, float3{0.0f, 0.0f, m22}};
}

// Right-handed frame with the normal as its first ROW. The manifold reads rows via row()/mul(),
// so linalg's column-major storage is transposed to put the axes in rows.
[[nodiscard]] inline float3x3 orthonormal(float3 normal) noexcept
{
    float3 t1 = std::fabs(normal.x) > std::fabs(normal.z)
            ? float3{-normal.y, normal.x, 0}
            : float3{0, -normal.z, normal.y};
    t1 = normalize(t1);
    const float3 t2 = cross(normal, t1);
    return transpose(float3x3{normal, t1, t2});
}

// Row read/write. linalg matrices index by column (`m[col][row]`); the contact code
// is written in terms of rows (the normal is row 0), so these keep that vocabulary.
[[nodiscard]] inline float3 row(const float3x3 &m, int i) noexcept { return m.row(i); }
inline void setRow(float3x3 &m, int i, float3 v) noexcept
{
    m.x[i] = v.x;
    m.y[i] = v.y;
    m.z[i] = v.z;
}

// The multiplicative identities, in the spots Eigen spelled ::Identity().
inline constexpr float3x3 identity3x3{float3{1, 0, 0}, float3{0, 1, 0}, float3{0, 0, 1}};
inline constexpr quat quatIdentity{0.0f, 0.0f, 0.0f, 1.0f};

// Eigen::AngleAxisf(angle, axis) with the arguments the engine always passed.
[[nodiscard]] inline quat angleAxis(float angle, float3 axis) noexcept
{
    return rotation_quat(axis, angle);
}

// The rotation vector from a to b: twice the vector part of a * conj(b).
[[nodiscard]] inline float3 operator-(quat a, quat b) noexcept {
  const quat d = qmul(a, qconj(b));
  return float3{d.x, d.y, d.z} * 2.0f;
}

// First-order integration: half omega (pure quaternion) left-multiplied, then renormalised; the
// added vector is the same small-angle quantity `operator-(quat,quat)` measures.
[[nodiscard]] inline quat operator+(quat a, float3 b) noexcept {
  // A pure quaternion has w = 0; the product's coefficients are added to a's, w included,
  // exactly as Eigen's coeffs()+ chain did, then the result is renormalised.
  const quat omega{b.x, b.y, b.z, 0.0f};
  const quat p = qmul(omega, a);
  return normalize(quat{a.x + p.x * 0.5f, a.y + p.y * 0.5f, a.z + p.z * 0.5f, a.w + p.w * 0.5f});
}

// Solve 6x6: LDL^T (~200 ns) -> pivoted LDL^T -> zero; sane = finite <10 m (~1e6 stiffness gave e15).
// M/dt^2 dominates (1 kg@240 Hz: 57,600 vs ~1e4); sequential float ops bit-identical, false -> retry pivoted.
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

// Pivoted LDL^T: bit-for-bit port of Eigen's LDLT (do not change); for indefinite blocks (pivots diverge).
// largest |diag| first-max, zero threshold; zero pivot fails on non-zero remainder; |d|>min else 0; swaps P then P^T.
inline bool ldlt6_pivoted_solve(const std::array<std::array<float, 6>, 6> &a_in, const std::array<float, 6> &b,
                                std::array<float, 6> &x) noexcept {
  std::array<std::array<float, 6>, 6> a = a_in; // working copy; lower triangle becomes L and D
  std::array<float, 6> d{};
  std::array<int, 6> t{0, 1, 2, 3, 4, 5};
  bool foundZeroPivot = false;

  for (int k = 0; k < 6; k++) {
    // Largest |diagonal| of the trailing block, first maximum (Eigen maxCoeff).
    int j = k;
    float best = std::fabs(a[k][k]);
    for (int i = k + 1; i < 6; i++) {
      const float candidate = std::fabs(a[i][i]);
      if (candidate > best) {
        best = candidate;
        j = i;
      }
    }
    t[k] = j;
    if (j != k) {
      // Symmetric row/column transposition (Eigen swaps the stored triangle; on
      // the mirrored full matrix this is the same assignment).
      for (int c = 0; c < 6; c++)
        std::swap(a[k][c], a[j][c]);
      for (int r = 0; r < 6; r++)
        std::swap(a[r][k], a[r][j]);
    }

    if (k > 0) {
      // temp = D(0..k-1) * A10; d(k,k) -= A10 . temp; A21 -= A20 * temp.
      std::array<float, 6> temp{};
      for (int i = 0; i < k; i++)
        temp[i] = d[i] * a[k][i];
      // A10 . temp with Eigen's redux pairing: full 4-lane packets reduce as
      // (t0+t2)+(t1+t3), scalars accumulate left to right after them.
      std::array<float, 6> prod{};
      for (int i = 0; i < k; i++)
        prod[i] = a[k][i] * temp[i];
      float dot;
      switch (k) {
        case 1: dot = prod[0]; break;
        case 2: dot = prod[0] + prod[1]; break;
        case 3: dot = (prod[0] + prod[1]) + prod[2]; break;
        case 4: dot = (prod[0] + prod[2]) + (prod[1] + prod[3]); break;
        default: dot = ((prod[0] + prod[2]) + (prod[1] + prod[3])) + prod[4]; break;
      }
      a[k][k] -= dot;
      // A21 -= A20 * temp, evaluated as Eigen's column-major GEMV does: one
      // saxpy per column, not a per-row dot product.
      for (int i = 0; i < k; i++) {
        const float ti = temp[i];
        for (int r = k + 1; r < 6; r++)
          a[r][k] -= a[r][i] * ti;
      }
    }

    const float pivot = a[k][k];
    const bool valid = std::fabs(pivot) > 0.0f;
    if (!valid) {
      for (int r = k + 1; r < 6; r++)
        if (a[r][k] != 0.0f)
          return false;
      foundZeroPivot = true;
    } else if (foundZeroPivot) {
      return false;
    }
    if (valid) {
      for (int r = k + 1; r < 6; r++)
        a[r][k] /= pivot;
    }
    d[k] = pivot;
  }

  // x = P b (running swaps in factorisation order), then L^-1, D^+, L^-T, then P^T.
  x = b;
  for (int k = 0; k < 6; k++)
    std::swap(x[k], x[t[k]]);

  // L^-1 y, column-oriented forward substitution (Eigen triangular solveInPlace).
  for (int jj = 0; jj < 5; jj++) {
    const float v = x[jj];
    if (v != 0.0f)
      for (int i = jj + 1; i < 6; i++)
        x[i] -= a[i][jj] * v;
  }

  const float tolerance = std::numeric_limits<float>::min();
  for (int i = 0; i < 6; i++)
    x[i] = std::fabs(d[i]) > tolerance ? x[i] / d[i] : 0.0f;

  // L^-T z, column-oriented back substitution over the transposed unit triangle.
  for (int jj = 5; jj >= 1; jj--) {
    const float v = x[jj];
    if (v != 0.0f)
      for (int i = jj - 1; i >= 0; i--)
        x[i] -= a[jj][i] * v;
  }

  for (int k = 5; k >= 0; k--)
    std::swap(x[k], x[t[k]]);

  float worst = 0.0f;
  for (int i = 0; i < 6; i++)
    worst = std::max(worst, std::abs(x[i]));
  return std::isfinite(worst) && worst < 1.0e1f;
}

#ifdef AVBD_DIAG
// Root-cause instrumentation: both 6x6 kernels failed sanity and the whole DOF update
// was zeroed. Reset by the driver between ticks.
inline std::atomic<long> diagUnsane{0};
#endif

inline void solve(float3x3 aLin, float3x3 aAng, float3x3 aCross, float3 bLin,
                  float3 bAng, float3 &xLin, float3 &xAng) noexcept {
  std::array<std::array<float, 6>, 6> a{};
  for (int r = 0; r < 3; r++) {
    for (int c = 0; c < 3; c++) {
      a[r][c] = aLin[c][r];
      a[r + 3][c + 3] = aAng[c][r];
      a[r + 3][c] = aCross[c][r];  // [aLin aCross^T; aCross aAng]
      a[r][c + 3] = aCross[r][c];
    }
  }

  const std::array<float, 6> b{bLin[0], bLin[1], bLin[2], bAng[0], bAng[1], bAng[2]};
  std::array<float, 6> x{};

  
  if (ldlt6_solve(a, b, x) || ldlt6_pivoted_solve(a, b, x)) {
    xLin = float3{x[0], x[1], x[2]};
    xAng = float3{x[3], x[4], x[5]};
    return;
  }

  // Degenerate block: freeze this DOF for one iteration instead of teleporting.
#ifdef AVBD_DIAG
  diagUnsane.fetch_add(1, std::memory_order_relaxed);
#endif
  xLin = float3{0, 0, 0};
  xAng = float3{0, 0, 0};
}

} // namespace avbd
