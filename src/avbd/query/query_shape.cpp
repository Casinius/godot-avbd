#include "avbd/query/query_shape.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace avbd::query {
namespace {

constexpr float kEpsilon = 1.0e-6f;

[[nodiscard]] float3 localPoint(const Shape &shape, const float3 &world) noexcept
{
    return qrot(qconj(shape.rotation), world - shape.center);
}

[[nodiscard]] float3 worldVector(const Shape &shape, const float3 &local) noexcept
{
    return qrot(shape.rotation, local);
}

[[nodiscard]] float3 boxHalf(const Shape &shape) noexcept
{
    return shape.half;
}

[[nodiscard]] bool insideLocal(const Shape &shape, const float3 &p, float margin) noexcept
{
    switch (shape.type)
    {
        case ShapeType::Sphere:
            return length2(p) <= (shape.radius + margin) * (shape.radius + margin);
        case ShapeType::Cylinder:
            return (p.x * p.x + p.y * p.y) <= (shape.radius + margin) * (shape.radius + margin)
                    && std::fabs(p.z) <= shape.halfHeight + margin;
        case ShapeType::Box:
        default: {
            const float3 h = boxHalf(shape);
            return std::fabs(p.x) <= h.x + margin && std::fabs(p.y) <= h.y + margin
                    && std::fabs(p.z) <= h.z + margin;
        }
    }
}

[[nodiscard]] float3 boxNormal(const float3 &p, const float3 &half) noexcept
{
    const float3 q{std::fabs(p.x) / (half.x > 0.0f ? half.x : 1.0f),
            std::fabs(p.y) / (half.y > 0.0f ? half.y : 1.0f),
            std::fabs(p.z) / (half.z > 0.0f ? half.z : 1.0f)};
    if (q.x >= q.y && q.x >= q.z)
        return float3{p.x < 0.0f ? -1.0f : 1.0f, 0, 0};
    if (q.y >= q.z)
        return float3{0, p.y < 0.0f ? -1.0f : 1.0f, 0};
    return float3{0, 0, p.z < 0.0f ? -1.0f : 1.0f};
}

[[nodiscard]] float3 cylinderNormal(const float3 &p, float halfHeight) noexcept
{
    // A point on a cap is at |z| == halfHeight; every other surface point is on the side.
    if (std::fabs(std::fabs(p.z) - halfHeight) <= 1.0e-4f)
        return float3{0, 0, p.z < 0.0f ? -1.0f : 1.0f};
    const float radial = std::sqrt(p.x * p.x + p.y * p.y);
    if (radial <= kEpsilon)
        return float3{1, 0, 0};
    return float3{p.x / radial, p.y / radial, 0};
}

[[nodiscard]] bool raySphere(const Shape &shape, const float3 &o, const float3 &d, float maxT,
        RayHit &r_hit) noexcept
{
    const float a = dot(d, d);
    if (a <= kEpsilon)
        return false;
    const float b = dot(o, d);
    const float c = dot(o, o) - shape.radius * shape.radius;
    const float disc = b * b - a * c;
    if (disc < 0.0f)
        return false;
    const float root = std::sqrt(disc);
    float t0 = (-b - root) / a;
    float t1 = (-b + root) / a;
    if (t0 > t1)
        std::swap(t0, t1);
    const bool inside = c < 0.0f;
    float t = inside ? t1 : t0;
    if (t < 0.0f)
        t = t1;
    if (t < 0.0f || t > maxT)
        return false;
    r_hit.t = t;
    r_hit.point = shape.center + qrot(shape.rotation, o + d * t);
    const float3 n = o + d * t;
    r_hit.normal = worldVector(shape, normalize(n));
    r_hit.inside = inside;
    return true;
}

[[nodiscard]] bool rayBox(const Shape &shape, const float3 &o, const float3 &d, float maxT,
        RayHit &r_hit) noexcept
{
    const float3 half = boxHalf(shape);
    float tEnter = 0.0f;
    float tExit = std::numeric_limits<float>::infinity();
    for (int i = 0; i < 3; i++)
    {
        if (std::fabs(d[i]) <= kEpsilon)
        {
            if (o[i] < -half[i] || o[i] > half[i])
                return false;
            continue;
        }
        const float inv = 1.0f / d[i];
        float t0 = (-half[i] - o[i]) * inv;
        float t1 = (half[i] - o[i]) * inv;
        if (t0 > t1)
            std::swap(t0, t1);
        tEnter = std::max(tEnter, t0);
        tExit = std::min(tExit, t1);
        if (tEnter > tExit)
            return false;
    }
    const bool inside = insideLocal(shape, o, 0.0f);
    float t = tEnter > 0.0f ? tEnter : tExit;
    if (t < 0.0f || t > maxT)
        return false;
    r_hit.t = t;
    const float3 p = o + d * t;
    r_hit.point = shape.center + qrot(shape.rotation, p);
    r_hit.normal = worldVector(shape, boxNormal(p, half));
    r_hit.inside = inside;
    return true;
}

[[nodiscard]] bool rayCylinder(const Shape &shape, const float3 &o, const float3 &d, float maxT,
        RayHit &r_hit) noexcept
{
    const float r = shape.radius;
    const float hh = shape.halfHeight;
    std::array<float, 4> hits{};
    int count = 0;

    const auto add = [&](float t) {
        if (t >= -kEpsilon && t <= maxT + kEpsilon && count < static_cast<int>(hits.size()))
            hits[static_cast<size_t>(count++)] = std::max(0.0f, t);
    };

    // Infinite cylinder side, clipped to the caps.
    const float a = d.x * d.x + d.y * d.y;
    if (a > kEpsilon)
    {
        const float b = o.x * d.x + o.y * d.y;
        const float c = o.x * o.x + o.y * o.y - r * r;
        const float disc = b * b - a * c;
        if (disc >= 0.0f)
        {
            const float root = std::sqrt(disc);
            const float t0 = (-b - root) / a;
            const float t1 = (-b + root) / a;
            for (const float t : {t0, t1})
            {
                const float z = o.z + d.z * t;
                if (std::fabs(z) <= hh + kEpsilon)
                    add(t);
            }
        }
    }
    else if (o.x * o.x + o.y * o.y > r * r + kEpsilon)
    {
        return false; // parallel to the axis and outside the radius
    }

    // Caps.
    if (std::fabs(d.z) > kEpsilon)
    {
        for (const float z : {hh, -hh})
        {
            const float t = (z - o.z) / d.z;
            const float px = o.x + d.x * t;
            const float py = o.y + d.y * t;
            if (px * px + py * py <= r * r + kEpsilon)
                add(t);
        }
    }

    if (count == 0)
        return false;
    for (int i = 1; i < count; i++)
    {
        const float value = hits[static_cast<size_t>(i)];
        int j = i - 1;
        while (j >= 0 && hits[static_cast<size_t>(j)] > value)
        {
            hits[static_cast<size_t>(j + 1)] = hits[static_cast<size_t>(j)];
            j--;
        }
        hits[static_cast<size_t>(j + 1)] = value;
    }

    const bool inside = insideLocal(shape, o, 0.0f);
    const float t = inside ? hits[static_cast<size_t>(count - 1)] : hits[0];
    if (t < 0.0f || t > maxT)
        return false;
    r_hit.t = t;
    const float3 p = o + d * t;
    r_hit.point = shape.center + qrot(shape.rotation, p);
    r_hit.normal = worldVector(shape, cylinderNormal(p, hh));
    r_hit.inside = inside;
    return true;
}

// Per-axis interval of t in [0, 1] where a moving AABB overlaps a static AABB.
[[nodiscard]] bool sweptAabbInterval(const float3 &aMin, const float3 &aMax,
        const float3 &bMin, const float3 &bMax, const float3 &motion, float &r_enter, float &r_exit) noexcept
{
    float enter = 0.0f;
    float exit = 1.0f;
    for (int i = 0; i < 3; i++)
    {
        const float m = motion[i];
        if (std::fabs(m) <= kEpsilon)
        {
            if (aMax[i] < bMin[i] || aMin[i] > bMax[i])
                return false;
            continue;
        }
        float t0 = (bMin[i] - aMax[i]) / m;
        float t1 = (bMax[i] - aMin[i]) / m;
        if (t0 > t1)
            std::swap(t0, t1);
        enter = std::max(enter, t0);
        exit = std::min(exit, t1);
        if (enter > exit)
            return false;
    }
    if (exit < 0.0f || enter > 1.0f)
        return false;
    r_enter = std::max(0.0f, enter);
    r_exit = std::min(1.0f, exit);
    return r_enter <= r_exit;
}

[[nodiscard]] bool contactAt(const Shape &base, const float3 &motion, float t, const Shape &target,
        float margin, float3x3 &r_basis, Manifold::Contact &r_contact) noexcept
{
    Shape posed = base;
    posed.center = base.center + motion * t;
    Manifold::Contact contacts[8];
    const int count = collideShapes(posed, target, std::span<Manifold::Contact>(contacts, 8), r_basis, margin);
    if (count <= 0)
        return false;
    const auto span = std::span<const Manifold::Contact>(contacts, static_cast<size_t>(count));
    const auto best = std::max_element(span.begin(), span.end(), [](const auto &a, const auto &b) {
        return length2(a.rA - a.rB) < length2(b.rA - b.rB);
    });
    r_contact = *best;
    return true;
}

} // namespace

bool intersectRay(const Shape &shape, const Ray &ray, RayHit &r_hit) noexcept
{
    const float maxT = ray.maxT > 0.0f ? ray.maxT : std::numeric_limits<float>::infinity();
    const float3 o = qrot(qconj(shape.rotation), ray.origin - shape.center);
    const float3 d = qrot(qconj(shape.rotation), ray.direction);

    switch (shape.type)
    {
        case ShapeType::Sphere:
            return raySphere(shape, o, d, maxT, r_hit);
        case ShapeType::Cylinder:
            return rayCylinder(shape, o, d, maxT, r_hit);
        case ShapeType::Box:
        default:
            return rayBox(shape, o, d, maxT, r_hit);
    }
}

float3 closestPoint(const Shape &shape, const float3 &point) noexcept
{
    const float3 p = localPoint(shape, point);
    float3 local = p;

    switch (shape.type)
    {
        case ShapeType::Sphere: {
            const float len = length(p);
            local = len > kEpsilon ? p * (shape.radius / len) : float3{shape.radius, 0, 0};
            break;
        }
        case ShapeType::Cylinder: {
            const float hh = shape.halfHeight;
            const float radial = std::sqrt(p.x * p.x + p.y * p.y);
            const float clampedZ = std::clamp(p.z, -hh, hh);
            if (radial > shape.radius)
            {
                const float scale = shape.radius / radial;
                local = float3{p.x * scale, p.y * scale, clampedZ};
            }
            else if (radial <= kEpsilon)
            {
                // On the axis: choose the nearer cap or the side at the same z.
                const float capDistance = hh - std::fabs(p.z);
                if (capDistance < shape.radius)
                    local = float3{0, 0, p.z < 0.0f ? -hh : hh};
                else
                    local = float3{shape.radius, 0, p.z};
            }
            else
            {
                const float sideDistance = shape.radius - radial;
                const float capDistance = hh - std::fabs(p.z);
                if (capDistance < sideDistance)
                    local = float3{p.x, p.y, p.z < 0.0f ? -hh : hh};
                else
                {
                    const float scale = shape.radius / radial;
                    local = float3{p.x * scale, p.y * scale, p.z};
                }
            }
            break;
        }
        case ShapeType::Box:
        default: {
            const float3 h = boxHalf(shape);
            local = float3{std::clamp(p.x, -h.x, h.x), std::clamp(p.y, -h.y, h.y),
                    std::clamp(p.z, -h.z, h.z)};
            break;
        }
    }
    return shape.center + qrot(shape.rotation, local);
}

bool containsPoint(const Shape &shape, const float3 &point, float margin) noexcept
{
    return insideLocal(shape, localPoint(shape, point), margin);
}

float supportRadius(const Shape &shape, const float3 &axis) noexcept
{
    const float3 n = normalize(axis);
    switch (shape.type)
    {
        case ShapeType::Sphere:
            return shape.radius;
        case ShapeType::Cylinder: {
            const float axial = std::fabs(dot(shape.axis, n));
            const float radial = std::sqrt(std::max(0.0f, 1.0f - axial * axial));
            return shape.radius * radial + shape.halfHeight * axial;
        }
        case ShapeType::Box:
        default: {
            const float3x3 rot = qmat(shape.rotation);
            return shape.half.x * std::fabs(dot(row(rot, 0), n))
                    + shape.half.y * std::fabs(dot(row(rot, 1), n))
                    + shape.half.z * std::fabs(dot(row(rot, 2), n));
        }
    }
}

float penetrationDepth(const Shape &a, const Shape &b, const float3 &normal) noexcept
{
    const float len = length(normal);
    if (len <= kEpsilon)
        return 0.0f;
    const float3 n = normal * (1.0f / len);
    return supportRadius(a, n) + supportRadius(b, n) - std::fabs(dot(b.center - a.center, n));
}

void shapeAabb(const Shape &shape, float3 &r_min, float3 &r_max) noexcept
{
    float3 half;
    switch (shape.type)
    {
        case ShapeType::Sphere:
            half = float3{shape.radius, shape.radius, shape.radius};
            break;
        case ShapeType::Cylinder:
            half = float3{shape.radius, shape.radius, shape.halfHeight};
            break;
        case ShapeType::Box:
        default:
            half = shape.half;
            break;
    }
    const float3x3 rot = qmat(shape.rotation);
    const float3 wx{
            std::fabs(rot[0][0]) * half.x + std::fabs(rot[1][0]) * half.y + std::fabs(rot[2][0]) * half.z,
            std::fabs(rot[0][1]) * half.x + std::fabs(rot[1][1]) * half.y + std::fabs(rot[2][1]) * half.z,
            std::fabs(rot[0][2]) * half.x + std::fabs(rot[1][2]) * half.y + std::fabs(rot[2][2]) * half.z};
    r_min = shape.center - wx;
    r_max = shape.center + wx;
}

int shapeContacts(const Shape &a, const Shape &b, std::span<Manifold::Contact> contacts,
        float3x3 &r_basis, float margin)
{
    return collideShapes(a, b, contacts, r_basis, margin);
}

SweepHit sweepShape(const Shape &base, const float3 &motion, const Shape &target, float margin)
{
    SweepHit out;
    const float3 zero{0, 0, 0};
    if (length2(motion) <= kEpsilon * kEpsilon)
    {
        float3x3 basis{};
        Manifold::Contact contact{};
        if (!contactAt(base, zero, 0.0f, target, margin, basis, contact))
            return out;
        out.hit = true;
        out.safe = 0.0f;
        out.unsafe = 0.0f;
        out.depth = penetrationDepth(base, target, row(basis, 0));
        out.point = base.center + qrot(base.rotation, contact.rA);
        out.normal = row(basis, 0);
        return out;
    }

    float3 aMin, aMax, bMin, bMax;
    shapeAabb(base, aMin, aMax);
    shapeAabb(target, bMin, bMax);
    if (margin > 0.0f)
    {
        const float3 m{margin, margin, margin};
        aMin -= m;
        aMax += m;
    }

    float enter = 0.0f;
    float exit = 1.0f;
    if (!sweptAabbInterval(aMin, aMax, bMin, bMax, motion, enter, exit))
        return out;

    float3x3 basis{};
    Manifold::Contact contact{};
    if (contactAt(base, motion, enter, target, margin, basis, contact))
    {
        out.hit = true;
        out.safe = enter;
        out.unsafe = enter;
        Shape posed = base;
        posed.center = base.center + motion * enter;
        out.depth = penetrationDepth(posed, target, row(basis, 0));
        out.point = posed.center + qrot(posed.rotation, contact.rA);
        out.normal = row(basis, 0);
        return out;
    }

    // The AABB interval is conservative. Walk it once to find the first contact, then
    // bisect that bracket for a stable safe/unsafe pair. The interval already bounds the
    // sweep, so a fast motion cannot skip past the target.
    constexpr int kSamples = 32;
    float low = enter;
    float high = exit;
    bool found = false;
    float previous = enter;
    for (int i = 1; i <= kSamples; i++)
    {
        const float t = enter + (exit - enter) * (static_cast<float>(i) / static_cast<float>(kSamples));
        if (contactAt(base, motion, t, target, margin, basis, contact))
        {
            low = previous;
            high = t;
            found = true;
            break;
        }
        previous = t;
    }
    if (!found && contactAt(base, motion, exit, target, margin, basis, contact))
    {
        low = previous;
        high = exit;
        found = true;
    }
    if (!found)
        return out;

    for (int i = 0; i < 12; i++)
    {
        const float mid = (low + high) * 0.5f;
        if (contactAt(base, motion, mid, target, margin, basis, contact))
            high = mid;
        else
            low = mid;
    }

    // Re-read the contact at the first colliding pose for the reported point and normal.
    if (!contactAt(base, motion, high, target, margin, basis, contact))
        return out;
    out.hit = true;
    out.safe = low;
    out.unsafe = high;
    Shape posed = base;
    posed.center = base.center + motion * high;
    out.depth = penetrationDepth(posed, target, row(basis, 0));
    out.point = posed.center + qrot(posed.rotation, contact.rA);
    out.normal = row(basis, 0);
    return out;
}

} // namespace avbd::query
