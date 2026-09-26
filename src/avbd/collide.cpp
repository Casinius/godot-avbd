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

#include "avbd/solver.h"
#include "avbd/bvh/node_storage.hpp"
#include "avbd/maths.h"
// cppcheck-suppress missingIncludeSystem
#include <algorithm>
// cppcheck-suppress missingIncludeSystem
#include <cfloat>
// cppcheck-suppress missingIncludeSystem
#include <cmath>
// cppcheck-suppress missingIncludeSystem
#include <cstddef>
// cppcheck-suppress missingIncludeSystem
#include <ranges>
// cppcheck-suppress missingIncludeSystem
#include <span>
// cppcheck-suppress missingIncludeSystem
#include <array>
// cppcheck-suppress missingIncludeSystem
#include <vector>

namespace avbd {

namespace
{
constexpr int MAX_CONTACTS = 8;
constexpr int MAX_POLY_VERTS = 16;
constexpr float SAT_AXIS_EPSILON = 1.0e-6f;
constexpr float PLANE_EPSILON = 1.0e-5f;
constexpr float CONTACT_MERGE_DIST_SQ = 1.0e-6f;

enum AxisType
{
    AXIS_FACE_A,
    AXIS_FACE_B,
    AXIS_EDGE
};

struct OBB
{
    float3 center;
    quat rotation;
    float3 half;
    float3 axis[3];
};

struct SatAxis
{
    AxisType type{};
    int indexA{};
    int indexB{};
    float separation{};
    float3 normalAB{};
    bool valid{};
};

struct FaceFrame
{
    int axisIndex{};
    float3 normal{};
    float3 center{};
    float3 u{};
    float3 v{};
    float extentU{};
    float extentV{};
};

inline OBB makeOBB(const Shape& shape)
{
    OBB box{};
    box.center = shape.center;
    box.rotation = shape.rotation;
    box.half = shape.half;
    box.axis[0] = shape.rotation * float3{1.0f, 0.0f, 0.0f};
    box.axis[1] = shape.rotation * float3{0.0f, 1.0f, 0.0f};
    box.axis[2] = shape.rotation * float3{0.0f, 0.0f, 1.0f};
    return box;
}

inline float absDot(float3 a, float3 b)
{
    return std::fabs(a.dot(b));
}

inline float3 supportPoint(const OBB& box, const float3& dir)
{
    float sx = dir.dot(box.axis[0]) >= 0.0f ? 1.0f : -1.0f;
    float sy = dir.dot(box.axis[1]) >= 0.0f ? 1.0f : -1.0f;
    float sz = dir.dot(box.axis[2]) >= 0.0f ? 1.0f : -1.0f;

    return box.center
        + box.axis[0] * (box.half.x() * sx)
        + box.axis[1] * (box.half.y() * sy)
        + box.axis[2] * (box.half.z() * sz);
}

inline void getFaceAxes(const OBB& box, int axisIndex, float3& u, float3& v, float& extentU, float& extentV)
{
    if (axisIndex == 0)
    {
        u = box.axis[1];
        v = box.axis[2];
        extentU = box.half.y();
        extentV = box.half.z();
    }
    else if (axisIndex == 1)
    {
        u = box.axis[0];
        v = box.axis[2];
        extentU = box.half.x();
        extentV = box.half.z();
    }
    else
    {
        u = box.axis[0];
        v = box.axis[1];
        extentU = box.half.x();
        extentV = box.half.y();
    }
}

inline void buildFaceFrame(const OBB& box, int axisIndex, const float3& outwardNormal, FaceFrame& frame)
{
    float sign = outwardNormal.dot(box.axis[axisIndex]) >= 0.0f ? 1.0f : -1.0f;
    frame.axisIndex = axisIndex;
    frame.normal = box.axis[axisIndex] * sign;
    frame.center = box.center + frame.normal * box.half[axisIndex];
    getFaceAxes(box, axisIndex, frame.u, frame.v, frame.extentU, frame.extentV);
}

inline int chooseIncidentFaceAxis(const OBB& box, const float3& referenceNormal)
{
    int axis = 0;
    float best = -FLT_MAX;

    for (size_t i = 0; i < 3; ++i)
    {
        float d = absDot(box.axis[i], referenceNormal);
        if (d > best)
        {
            best = d;
            axis = i;
        }
    }

    return axis;
}

inline void buildIncidentFace(const OBB& box, int axisIndex, const float3& referenceNormal, float3 outVerts[4])
{
    float sign = box.axis[axisIndex].dot(referenceNormal) > 0.0f ? -1.0f : 1.0f;
    float3 faceNormal = box.axis[axisIndex] * sign;
    float3 faceCenter = box.center + faceNormal * box.half[axisIndex];

    float3 u;
    float3 v;
    float extentU;
    float extentV;
    getFaceAxes(box, axisIndex, u, v, extentU, extentV);

    outVerts[0] = faceCenter + u * extentU + v * extentV;
    outVerts[1] = faceCenter - u * extentU + v * extentV;
    outVerts[2] = faceCenter - u * extentU - v * extentV;
    outVerts[3] = faceCenter + u * extentU - v * extentV;
}

inline int clipPolygonAgainstPlane(const float3* inVerts, int inCount, const float3& planeNormal, float planeOffset, float3* outVerts)
{
    if (inCount <= 0)
        return 0;

    int outCount = 0;
    float3 a = inVerts[inCount - 1];
    float da = planeNormal.dot(a) - planeOffset;

    for (size_t i = 0; i < static_cast<size_t>(inCount); ++i)
    {
        float3 b = inVerts[i];
        float db = planeNormal.dot(b) - planeOffset;

        bool aInside = da <= PLANE_EPSILON;
        bool bInside = db <= PLANE_EPSILON;

        if (aInside != bInside)
        {
            float t = 0.0f;
            float denom = da - db;
            if (std::fabs(denom) > SAT_AXIS_EPSILON)
                t = std::clamp(da / denom, 0.0f, 1.0f);

            if (outCount < MAX_POLY_VERTS)
                outVerts[outCount++] = a + (b - a) * t;
        }

        if (bInside && outCount < MAX_POLY_VERTS)
            outVerts[outCount++] = b;

        a = b;
        da = db;
    }

    return outCount;
}

inline bool addContact(const Shape& shapeA, const Shape& shapeB, std::span<Manifold::Contact> contacts, int& contactCount, float3* contactMidpoints, float3 xA, float3 xB, int featureKey)
{
    float3 midpoint = (xA + xB) * 0.5f;

    for (size_t i = 0; i < static_cast<size_t>(contactCount); ++i)
    {
        float3 d = midpoint - contactMidpoints[i];
        if (d.squaredNorm() < CONTACT_MERGE_DIST_SQ)
            return false;
    }

    if (contactCount >= MAX_CONTACTS)
        return false;

    Manifold::FeaturePair feature{};
    feature.key = featureKey;

    Manifold::Contact& c = contacts[contactCount];
    c.feature = feature;
    c.rA = shapeA.rotation.conjugate() * (xA - shapeA.center);
    c.rB = shapeB.rotation.conjugate() * (xB - shapeB.center);
    contactMidpoints[contactCount] = midpoint;
    ++contactCount;

    return true;
}

inline bool testAxis(const OBB& boxA, const OBB& boxB, const float3& delta, const float3& axis, AxisType type, int indexA, int indexB, SatAxis& best)
{
    float lenSq = axis.squaredNorm();
    if (lenSq < SAT_AXIS_EPSILON)
        return true;

    float invLen = 1.0f / std::sqrt(lenSq);
    float3 n = axis * invLen;
    if (n.dot(delta) < 0.0f)
        n = -n;

    float distance = std::abs(delta.dot(n));

    float rA =
        boxA.half.x() * absDot(n, boxA.axis[0]) +
        boxA.half.y() * absDot(n, boxA.axis[1]) +
        boxA.half.z() * absDot(n, boxA.axis[2]);

    float rB =
        boxB.half.x() * absDot(n, boxB.axis[0]) +
        boxB.half.y() * absDot(n, boxB.axis[1]) +
        boxB.half.z() * absDot(n, boxB.axis[2]);

    float separation = distance - (rA + rB);
    if (separation > 0.0f)
        return false;

    if (!best.valid || separation > best.separation)
    {
        best.valid = true;
        best.type = type;
        best.indexA = indexA;
        best.indexB = indexB;
        best.separation = separation;
        best.normalAB = n;
    }

    return true;
}

inline void supportEdge(const OBB& box, int axisIndex, const float3& dir, float3& edgeA, float3& edgeB)
{
    int axis1 = (axisIndex + 1) % 3;
    int axis2 = (axisIndex + 2) % 3;

    float sign1 = dir.dot(box.axis[axis1]) >= 0.0f ? 1.0f : -1.0f;
    float sign2 = dir.dot(box.axis[axis2]) >= 0.0f ? 1.0f : -1.0f;

    float3 edgeCenter = box.center
        + box.axis[axis1] * (box.half[axis1] * sign1)
        + box.axis[axis2] * (box.half[axis2] * sign2);

    edgeA = edgeCenter - box.axis[axisIndex] * box.half[axisIndex];
    edgeB = edgeCenter + box.axis[axisIndex] * box.half[axisIndex];
}

inline void closestPointsOnSegments(const float3& p0, const float3& p1, const float3& q0, const float3& q1, float3& c0, float3& c1)
{
    float3 d1 = p1 - p0;
    float3 d2 = q1 - q0;
    float3 r = p0 - q0;
    float a = d1.dot(d1);
    float e = d2.dot(d2);
    float f = d2.dot(r);

    float s = 0.0f;
    float t = 0.0f;

    if (a <= SAT_AXIS_EPSILON && e <= SAT_AXIS_EPSILON)
    {
        c0 = p0;
        c1 = q0;
        return;
    }

    if (a <= SAT_AXIS_EPSILON)
    {
        t = std::clamp(f / e, 0.0f, 1.0f);
    }
    else
    {
        float c = d1.dot(r);
        if (e <= SAT_AXIS_EPSILON)
        {
            s = std::clamp(-c / a, 0.0f, 1.0f);
        }
        else
        {
            float b = d1.dot(d2);
            float denom = a * e - b * b;

            if (std::fabs(denom) > SAT_AXIS_EPSILON)
                s = std::clamp((b * f - c * e) / denom, 0.0f, 1.0f);

            t = (b * s + f) / e;

            if (t < 0.0f)
            {
                t = 0.0f;
                s = std::clamp(-c / a, 0.0f, 1.0f);
            }
            else if (t > 1.0f)
            {
                t = 1.0f;
                s = std::clamp((b - c) / a, 0.0f, 1.0f);
            }
        }
    }

    c0 = p0 + d1 * s;
    c1 = q0 + d2 * t;
}

inline int buildFaceManifold(const Shape& shapeA, const Shape& shapeB, const OBB& boxA, const OBB& boxB, bool referenceIsA, int referenceAxis, const float3& normalAB, std::span<Manifold::Contact> contacts)
{
    const OBB& referenceBox = referenceIsA ? boxA : boxB;
    const OBB& incidentBox = referenceIsA ? boxB : boxA;
    float3 referenceOutward = referenceIsA ? normalAB : -normalAB;

    FaceFrame referenceFace{};
    buildFaceFrame(referenceBox, referenceAxis, referenceOutward, referenceFace);

    int incidentAxis = chooseIncidentFaceAxis(incidentBox, referenceFace.normal);

    std::array<float3, MAX_POLY_VERTS> clip0;
    std::array<float3, MAX_POLY_VERTS> clip1;
    buildIncidentFace(incidentBox, incidentAxis, referenceFace.normal, clip0.data());
    int count = 4;

    float3 n0 = referenceFace.u;
    float o0 = n0.dot(referenceFace.center) + referenceFace.extentU;
    count = clipPolygonAgainstPlane(clip0.data(), count, n0, o0, clip1.data());
    if (!count)
        return 0;

    float3 n1 = -referenceFace.u;
    float o1 = n1.dot(referenceFace.center) + referenceFace.extentU;
    count = clipPolygonAgainstPlane(clip1.data(), count, n1, o1, clip0.data());
    if (!count)
        return 0;

    float3 n2 = referenceFace.v;
    float o2 = n2.dot(referenceFace.center) + referenceFace.extentV;
    count = clipPolygonAgainstPlane(clip0.data(), count, n2, o2, clip1.data());
    if (!count)
        return 0;

    float3 n3 = -referenceFace.v;
    float o3 = n3.dot(referenceFace.center) + referenceFace.extentV;
    count = clipPolygonAgainstPlane(clip1.data(), count, n3, o3, clip0.data());
    if (!count)
        return 0;

    int contactCount = 0;
    std::array<float3, MAX_CONTACTS> contactMidpoints;
    int featurePrefix = (referenceIsA ? AXIS_FACE_A : AXIS_FACE_B) << 24;
    featurePrefix |= (referenceAxis & 0xFF) << 16;
    featurePrefix |= (incidentAxis & 0xFF) << 8;

    for (size_t i = 0; i < static_cast<size_t>(count) && contactCount < MAX_CONTACTS; ++i)
    {
        float3 pIncident = clip0[i];
        float distance = (pIncident - referenceFace.center).dot(referenceFace.normal);
        if (distance > PLANE_EPSILON)
            continue;

        float3 pReference = pIncident - referenceFace.normal * distance;
        float3 xA = referenceIsA ? pReference : pIncident;
        float3 xB = referenceIsA ? pIncident : pReference;

        addContact(shapeA, shapeB, contacts, contactCount, contactMidpoints.data(), xA, xB, featurePrefix | (i & 0xFF));
    }

    if (!contactCount)
    {
        float3 xA = supportPoint(boxA, normalAB);
        float3 xB = supportPoint(boxB, -normalAB);
        addContact(shapeA, shapeB, contacts, contactCount, contactMidpoints.data(), xA, xB, featurePrefix);
    }

    return contactCount;
}

inline int buildEdgeContact(const Shape& shapeA, const Shape& shapeB, const OBB& boxA, const OBB& boxB, int axisA, int axisB, const float3& normalAB, std::span<Manifold::Contact> contacts)
{
    float3 a0;
    float3 a1;
    float3 b0;
    float3 b1;
    supportEdge(boxA, axisA, normalAB, a0, a1);
    supportEdge(boxB, axisB, -normalAB, b0, b1);

    float3 xA;
    float3 xB;
    closestPointsOnSegments(a0, a1, b0, b1, xA, xB);

    int contactCount = 0;
    std::array<float3, MAX_CONTACTS> contactMidpoints;
    int featureKey = (AXIS_EDGE << 24) | ((axisA & 0xFF) << 8) | (axisB & 0xFF);
    addContact(shapeA, shapeB, contacts, contactCount, contactMidpoints.data(), xA, xB, featureKey);

    if (!contactCount)
    {
        xA = supportPoint(boxA, normalAB);
        xB = supportPoint(boxB, -normalAB);
        addContact(shapeA, shapeB, contacts, contactCount, contactMidpoints.data(), xA, xB, featureKey);
    }

    return contactCount;
}


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

constexpr int RIM_SAMPLES = 12;
// Interior rings along a cylinder's side, between the two caps.
constexpr int AXIAL_SAMPLES = 3;
// Upper bound on the candidate contacts one shape pair can produce before they are ranked.
constexpr int MAX_CANDIDATES = 96;

inline Shape makeShape(const Rigid* body)
{
    Shape shape;
    shape.type = body->shape;
    shape.center = body->positionLin;
    shape.rotation = body->positionAng;
    shape.half = body->size * 0.5f;
    shape.radius = body->size.x();
    shape.halfHeight = body->size.z() * 0.5f;
    shape.axis = body->positionAng * float3{0, 0, 1};
    return shape;
}

// Whether a point is inside a box, and if so the shortest way out: the axis of least
// penetration. Returns false when the point is outside.
inline bool pointInBox(const Shape& box, const float3& p, float3& r_push, float& r_depth)
{
    const float3 local = box.rotation.conjugate() * (p - box.center);

    float bestDepth = FLT_MAX;
    int bestAxis = 0;
    float bestSign = 1.0f;
    for (size_t i = 0; i < 3; ++i)
    {
        const float depth = box.half[i] - std::fabs(local[i]);
        if (depth <= 0.0f)
            return false;
        if (depth < bestDepth)
        {
            bestDepth = depth;
            bestAxis = i;
            bestSign = local[i] >= 0.0f ? 1.0f : -1.0f;
        }
    }

    float3 localPush{0, 0, 0};
    localPush[bestAxis] = bestSign;
    r_push = box.rotation * localPush;
    r_depth = bestDepth;
    return true;
}

// Closest point on a solid capped cylinder to a point, plus the outward direction from the
// cylinder's surface to that point and the distance.
inline float3 closestPointOnCylinder(const Shape& cyl, const float3& p, float3& r_outward, float& r_distance)
{
    const float3 local = cyl.rotation.conjugate() * (p - cyl.center);

    // Clamp along the axis, then radially.
    const float t = std::clamp(local.y(), -cyl.halfHeight, cyl.halfHeight);
    float2 radial(local.x(), local.z());
    const float lenSq = radial.squaredNorm();

    float3 closest{0, t, 0};
    if (lenSq > cyl.radius * cyl.radius)
    {
        const float scale = cyl.radius / std::sqrt(lenSq);
        closest.x() = radial.x() * scale;
        closest.y() = radial.y() * scale;
    }
    else if (lenSq > 1.0e-12f)
    {
        closest.x() = radial.x();
        closest.y() = radial.y();
    }
    else
    {
        // On the axis: the nearest surface point is on the radius, pick a fixed direction.
        closest.x() = cyl.radius;
    }

    const float3 closestLocal = cyl.rotation * closest + cyl.center;
    const float3 d = p - closestLocal;
    r_distance = d.norm();
    r_outward = r_distance > 1.0e-9f ? d / r_distance : cyl.rotation * float3{1, 0, 0};
    return closestLocal;
}

// The shortest way from a point inside a solid cylinder out to its surface: along the cap, or
// radially, whichever is closer. Returns false when the point is outside.
inline bool pointInCylinder(const Shape& cyl, const float3& p, float3& r_push, float& r_surface)
{
    const float3 local = cyl.rotation.conjugate() * (p - cyl.center);
    const float radialSq = local.x() * local.x() + local.y() * local.y();
    if (radialSq >= cyl.radius * cyl.radius || std::fabs(local.z()) >= cyl.halfHeight)
        return false;

    const float radialDepth = cyl.radius - std::sqrt(radialSq);
    const float capDepth = cyl.halfHeight - std::fabs(local.z());

    float3 pushLocal{0, 0, 0};
    if (capDepth < radialDepth)
    {
        pushLocal.z() = local.z() >= 0.0f ? 1.0f : -1.0f;
        r_surface = capDepth;
    }
    else
    {
        const float len = std::sqrt(radialSq);
        if (len > 1.0e-6f)
        {
            pushLocal.x() = local.x() / len;
            pushLocal.y() = local.y() / len;
        }
        else
        {
            pushLocal.x() = 1.0f; // on the axis: any radial direction will do
        }
        r_surface = radialDepth;
    }

    r_push = cyl.rotation * pushLocal;
    return true;
}

// A candidate contact before it is handed to the manifold, with the penetration depth that
// decides whether it is worth keeping.
struct Candidate
{
    float3 xA{};
    float3 xB{};
    float3 normal{}; // points from the shape the sample belongs to, outwards
    float depth{};
};

// Emit the deepest MAX_CONTACTS candidates. Ranking by depth rather than taking them in sampling
// order is what keeps a contact set sensible when there are more samples than the budget: the
// shallow ones are the ones to drop.
inline void emitRanked(const Shape& shapeA, const Shape& shapeB,
                       std::span<const Candidate> candidates,
                       std::span<Manifold::Contact> contacts, int& contactCount)
{
    // Make a copy since we need to sort the candidates
    std::vector<Candidate> sortedCandidates(candidates.begin(), candidates.end());

    std::ranges::stable_sort(sortedCandidates, std::greater{},
            [](const Candidate& c) { return c.depth; });

    std::array<float3, MAX_CONTACTS> midpoints;
    contactCount = 0;
    for (const Candidate& c : std::views::take(sortedCandidates, MAX_CONTACTS))
    {
        if (!addContact(shapeA, shapeB, contacts, contactCount, midpoints.data(), c.xA, c.xB, contactCount + 1))
            break;
    }
}

// Samples a cylinder's surface: the two cap centres, the two rim circles, and rings partway
// along the side (which is what gives a lying cylinder a line of contacts rather than just its
// two ends).
//
// The angles are visited in an interleaved order rather than 0, 1, 2, ... because the contact
// budget is smaller than the sample count: adding contacts stops at MAX_CONTACTS, and taking the
// first N of an ordered ring would cover one arc and leave the other side unsupported - which
// tips the body over. With a stride coprime to the sample count, any prefix is spread around the
// circle.
inline void cylinderSamples(const Shape& cyl, float3* out, int& count)
{
    count = 0;
    const float3 x = cyl.rotation * float3{1, 0, 0};
    const float3 y = cyl.rotation * float3{0, 1, 0};

    // Both caps, centres first so a flat contact has a central point.
    for (int end = -1; end <= 1; end += 2)
    {
        const float3 capCentre = cyl.center + cyl.axis * (cyl.halfHeight * (float)end);
        out[count++] = capCentre;
        for (size_t i = 0; i < RIM_SAMPLES; ++i)
        {
            const int step = (i * 5) % RIM_SAMPLES;
            const float a = 6.28318530718f * (float)step / (float)RIM_SAMPLES;
            out[count++] = capCentre + (x * std::cos(a) + y * std::sin(a)) * cyl.radius;
        }
    }

    // Interior rings along the side.
    for (int ring = 1; ring < AXIAL_SAMPLES; ++ring)
    {
        const float t = -cyl.halfHeight + 2.0f * cyl.halfHeight * (float)ring / (float)AXIAL_SAMPLES;
        const float3 ringCentre = cyl.center + cyl.axis * t;
        for (size_t i = 0; i < RIM_SAMPLES; ++i)
        {
            const int step = (i * 5) % RIM_SAMPLES;
            const float a = 6.28318530718f * (float)step / (float)RIM_SAMPLES;
            out[count++] = ringCentre + (x * std::cos(a) + y * std::sin(a)) * cyl.radius;
        }
    }
}

inline int collideSphereSphere(const Shape& a, const Shape& b,
        std::span<Manifold::Contact> contacts, float3x3& basisOut)
{
    const float3 d = b.center - a.center;
    const float distSq = d.squaredNorm();
    const float sum = a.radius + b.radius;
    if (distSq >= sum * sum)
        return 0;

    const float dist = std::sqrt(distSq);
    const float3 normalAB = dist > 1.0e-6f ? d / dist : float3{0, 0, 1};
    basisOut = orthonormal(-normalAB);

    int count = 0;
    std::array<float3, MAX_CONTACTS> midpoints;
    addContact(a, b, contacts, count, midpoints.data(), a.center + normalAB * a.radius,
            b.center - normalAB * b.radius, 1);
    return count;
}

// Sphere against box, in either order. Exact.
inline int collideSphereBox(const Shape& a, const Shape& b, bool sphereIsA,
        std::span<Manifold::Contact> contacts, float3x3& basisOut)
{
    const Shape& sphere = sphereIsA ? a : b;
    const Shape& box = sphereIsA ? b : a;

    const float3 local = box.rotation.conjugate() * (sphere.center - box.center);
    float3 clamped{std::clamp(local.x(), -box.half.x(), box.half.x()), std::clamp(local.y(), -box.half.y(), box.half.y()),
            std::clamp(local.z(), -box.half.z(), box.half.z())};
    const float3 closestLocal = box.rotation * clamped + box.center;
    const float3 d = sphere.center - closestLocal;
    const float distSq = d.squaredNorm();

    if (distSq >= sphere.radius * sphere.radius)
        return 0;

    const float dist = std::sqrt(distSq);
    float3 normalSphereToBox;
    if (dist > 1.0e-6f)
    {
        normalSphereToBox = d / dist; // points from the box surface towards the sphere centre
    }
    else
    {
        // The centre is inside the box: push out along the axis of least penetration. This is
        // the case a sphere spawned inside geometry hits, and it needs a defined direction.
        float3 push;
        float depth;
        if (!pointInBox(box, sphere.center, push, depth))
            return 0;
        normalSphereToBox = push;
    }

    // `normalSphereToBox` leaves the box towards the sphere. The contact normal runs from A to
    // B, so it is the opposite of that when the sphere is A.
    const float3 normalAB = sphereIsA ? -normalSphereToBox : normalSphereToBox;
    basisOut = orthonormal(-normalAB);

    // The point of the sphere nearest the box, not the far side: the pair of points then spans
    // the penetration, which is what the solver's constraint value measures.
    const float3 onSphere = sphere.center - normalSphereToBox * sphere.radius;
    const float3 onBox = closestLocal;
    const float3 xA = sphereIsA ? onSphere : onBox;
    const float3 xB = sphereIsA ? onBox : onSphere;

    int count = 0;
    std::array<float3, MAX_CONTACTS> midpoints;
    addContact(a, b, contacts, count, midpoints.data(), xA, xB, 1);
    return count;
}

// Sphere against the closest point on a solid cylinder. Exact for the curved side and both caps.
inline int collideSphereCylinder(const Shape& a, const Shape& b, bool sphereIsA,
        std::span<Manifold::Contact> contacts, float3x3& basisOut)
{
    const Shape& sphere = sphereIsA ? a : b;
    const Shape& cyl = sphereIsA ? b : a;

    float3 outward;
    float distance;
    const float3 onCylinder = closestPointOnCylinder(cyl, sphere.center, outward, distance);
    if (distance >= sphere.radius)
        return 0;

// `outward` points from the cylinder surface towards the sphere centre, so it is the A-to-B
// normal only when the cylinder is A.
    const float3 normalAB = sphereIsA ? -outward : outward;
    basisOut = orthonormal(-normalAB);

    const float3 onSphere = sphere.center - outward * sphere.radius;
    const float3 xA = sphereIsA ? onSphere : onCylinder;
    const float3 xB = sphereIsA ? onCylinder : onSphere;

    int count = 0;
    std::array<float3, MAX_CONTACTS> midpoints;
    addContact(a, b, contacts, count, midpoints.data(), xA, xB, 1);
    return count;
}

// Cylinder against box, in either order. Both directions are tested - the cylinder's surface
// samples against the box, and the box's corners against the cylinder - so neither shape can slip
// through the other however they are arranged.
inline int collideCylinderBox(const Shape& a, const Shape& b, bool cylinderIsA,
        std::span<Manifold::Contact> contacts, float3x3& basisOut)
{
    const Shape& cyl = cylinderIsA ? a : b;
    const Shape& box = cylinderIsA ? b : a;

    std::array<Candidate, MAX_CANDIDATES> candidates;
    int candidateCount = 0;
    float3 normalAB{0, 0, 0};

    // Cylinder surface points that are inside the box.
    float3 samples[2 * (RIM_SAMPLES + 1) + RIM_SAMPLES * (AXIAL_SAMPLES - 1)];
    int sampleCount = 0;
    cylinderSamples(cyl, samples, sampleCount);

    for (size_t i = 0; i < static_cast<size_t>(sampleCount) && candidateCount < MAX_CANDIDATES; ++i)
    {
        float3 pushOut;
        float depth;
        if (!pointInBox(box, samples[i], pushOut, depth))
            continue;

        Candidate& c = candidates[candidateCount++];
        c.xA = cylinderIsA ? samples[i] : samples[i] + pushOut * depth;
        c.xB = cylinderIsA ? samples[i] + pushOut * depth : samples[i];
        c.normal = cylinderIsA ? -pushOut : pushOut;
        c.depth = depth;
        normalAB += c.normal;
    }

    // Box corners inside the cylinder.
    for (int corner = 0; corner < 8 && candidateCount < MAX_CANDIDATES; ++corner)
    {
        const float3 cornerLocal{float((corner & 1) ? box.half.x() : -box.half.x()),
                float((corner & 2) ? box.half.y() : -box.half.y()),
                float((corner & 4) ? box.half.z() : -box.half.z())};
        const float3 cornerWorld = box.rotation * cornerLocal + box.center;

        float3 pushOut;
        float surfaceDist;
        if (!pointInCylinder(cyl, cornerWorld, pushOut, surfaceDist))
            continue;

        Candidate& c = candidates[candidateCount++];
        c.xA = cylinderIsA ? cornerWorld + pushOut * surfaceDist : cornerWorld;
        c.xB = cylinderIsA ? cornerWorld : cornerWorld + pushOut * surfaceDist;
        c.normal = cylinderIsA ? pushOut : -pushOut;
        c.depth = surfaceDist;
        normalAB += c.normal;
    }

    if (!candidateCount)
        return 0;

    basisOut = orthonormal(-(normalAB.squaredNorm() > 1.0e-12f ? normalAB.normalized() : cyl.axis));

    int count = 0;
    emitRanked(a, b, candidates, contacts, count);
    return count;
}

// Cylinder against cylinder: each one's surface samples tested against the other.
inline int collideCylinderCylinder(const Shape& a, const Shape& b,
        std::span<Manifold::Contact> contacts, float3x3& basisOut)
{
    std::array<Candidate, MAX_CANDIDATES> candidates;
    int candidateCount = 0;
    float3 normalAB{0, 0, 0};

    float3 samples[2 * (RIM_SAMPLES + 1) + RIM_SAMPLES * (AXIAL_SAMPLES - 1)];
    const Shape* pair[2] = {&a, &b};

    for (int pass = 0; pass < 2; ++pass)
    {
        const Shape& from = *pair[pass];
        const Shape& into = *pair[pass == 0 ? 1 : 0];
        const bool fromIsA = pass == 0;

        int sampleCount = 0;
        cylinderSamples(from, samples, sampleCount);

        for (size_t i = 0; i < static_cast<size_t>(sampleCount) && candidateCount < MAX_CANDIDATES; ++i)
        {
            float3 pushOut;
            float surfaceDist;
            if (!pointInCylinder(into, samples[i], pushOut, surfaceDist))
                continue;

            Candidate& c = candidates[candidateCount++];
            c.xA = fromIsA ? samples[i] : samples[i] + pushOut * surfaceDist;
            c.xB = fromIsA ? samples[i] + pushOut * surfaceDist : samples[i];
            c.normal = fromIsA ? -pushOut : pushOut;
            c.depth = surfaceDist;
            normalAB += c.normal;
        }
    }

    if (!candidateCount)
        return 0;

    basisOut = orthonormal(-(normalAB.squaredNorm() > 1.0e-12f ? normalAB.normalized() : a.axis));

    int count = 0;
    emitRanked(a, b, candidates, contacts, count);
    return count;
}
} // namespace

int collideShapes(const Shape& a, const Shape& b,
                  std::span<Manifold::Contact> contacts,
                  float3x3& basisOut)
{
    // Shapes other than boxes take their own paths; box-box keeps the SAT implementation below
    // untouched, so its behaviour is bit-for-bit what it was.
    if (a.type != ShapeType::Box || b.type != ShapeType::Box)
    {
        if (a.type == ShapeType::Sphere && b.type == ShapeType::Sphere)
            return collideSphereSphere(a, b, contacts, basisOut);

        if (a.type == ShapeType::Sphere || b.type == ShapeType::Sphere)
        {
            if (a.type == ShapeType::Cylinder || b.type == ShapeType::Cylinder)
            {
                const bool sphereIsA = a.type == ShapeType::Sphere;
                return collideSphereCylinder(a, b, sphereIsA, contacts, basisOut);
            }
            const bool sphereIsA = a.type == ShapeType::Sphere;
            return collideSphereBox(a, b, sphereIsA, contacts, basisOut);
        }

        if (a.type == ShapeType::Cylinder && b.type == ShapeType::Cylinder)
            return collideCylinderCylinder(a, b, contacts, basisOut);

        // One cylinder, one box.
        const bool cylinderIsA = a.type == ShapeType::Cylinder;
        return collideCylinderBox(a, b, cylinderIsA, contacts, basisOut);
    }

    OBB boxA = makeOBB(a);
    OBB boxB = makeOBB(b);
    float3 delta = boxB.center - boxA.center;

    SatAxis bestFace{};
    bestFace.separation = -FLT_MAX;
    bestFace.valid = false;

    SatAxis bestEdge{};
    bestEdge.separation = -FLT_MAX;
    bestEdge.valid = false;

    for (size_t i = 0; i < 3; ++i)
    {
        if (!testAxis(boxA, boxB, delta, boxA.axis[i], AXIS_FACE_A, i, -1, bestFace))
            return 0;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        if (!testAxis(boxA, boxB, delta, boxB.axis[i], AXIS_FACE_B, -1, i, bestFace))
            return 0;
    }

    for (size_t i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            float3 axis = boxA.axis[i].cross(boxB.axis[j]);
            if (!testAxis(boxA, boxB, delta, axis, AXIS_EDGE, i, j, bestEdge))
                return 0;
        }
    }

    if (!bestFace.valid)
        return 0;

    SatAxis best = bestFace;
    if (bestEdge.valid)
    {
        const float edgeRelTol = 0.95f;
        const float edgeAbsTol = 0.01f;
        if (edgeRelTol * bestEdge.separation > bestFace.separation + edgeAbsTol)
            best = bestEdge;
    }

    basisOut = orthonormal(-best.normalAB);

    if (best.type == AXIS_EDGE)
        return buildEdgeContact(a, b, boxA, boxB, best.indexA, best.indexB, best.normalAB, contacts);

    if (best.type == AXIS_FACE_A)
        return buildFaceManifold(a, b, boxA, boxB, true, best.indexA, best.normalAB, contacts);

    return buildFaceManifold(a, b, boxA, boxB, false, best.indexB, best.normalAB, contacts);
}

int Manifold::collide(Rigid* bodyA, Rigid* bodyB,
                      std::span<Manifold::Contact> contacts,
                      float3x3& basisOut)
{
    // Thin wrapper: flatten both bodies into shape queries and reuse the solver-free path.
    return collideShapes(makeShape(bodyA), makeShape(bodyB), contacts, basisOut);
}

} // namespace avbd
