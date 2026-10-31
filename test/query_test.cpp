/*
 * Numeric tests for the query extension (src/avbd/query).
 *
 * No Godot, no solver step: exact ray/shape, closest point, AABB and sweep behaviour.
 *
 *   xmake run avbd_query_test
 */

#include <cmath>
#include <cstdarg>
#include <cstdio>

#include "avbd/query/query_shape.hpp"

using namespace avbd;

static int g_checks = 0;
static int g_failures = 0;

static void report(bool ok, const char *name, const char *fmt, ...) {
    g_checks++;
    if (!ok) {
        g_failures++;
    }
    std::printf("%s %-34s ", ok ? "  ok " : "FAIL ", name);
    va_list ap;
    va_start(ap, fmt);
    std::vprintf(fmt, ap);
    va_end(ap);
    std::printf("\n");
}

static bool near(float a, float b, float eps) {
    return std::fabs(a - b) <= eps;
}

static Shape box(const float3 &center, const float3 &half) {
    Shape s;
    s.type = ShapeType::Box;
    s.center = center;
    s.rotation = quatIdentity;
    s.half = half;
    s.axis = float3{0, 0, 1};
    return s;
}

static Shape sphere(const float3 &center, float radius) {
    Shape s;
    s.type = ShapeType::Sphere;
    s.center = center;
    s.rotation = quatIdentity;
    s.radius = radius;
    s.axis = float3{0, 0, 1};
    return s;
}

static Shape cylinder(const float3 &center, float radius, float halfHeight, quat rotation = quatIdentity) {
    Shape s;
    s.type = ShapeType::Cylinder;
    s.center = center;
    s.rotation = rotation;
    s.radius = radius;
    s.halfHeight = halfHeight;
    s.axis = qrot(rotation, float3{0, 0, 1});
    return s;
}

static void test_rays() {
    query::RayHit hit{};

    query::Ray ray;
    ray.origin = {-5, 0, 0};
    ray.direction = {1, 0, 0};
    ray.maxT = 10.0f;

    const Shape ball = sphere({0, 0, 0}, 1.0f);
    bool ok = query::intersectRay(ball, ray, hit);
    report(ok && near(hit.t, 4.0f, 1e-5f) && near(hit.point.x, -1.0f, 1e-5f)
                    && near(hit.normal.x, -1.0f, 1e-5f) && !hit.inside,
            "ray/sphere entry", "t=%.5f point=(%.4f,%.4f,%.4f)", hit.t, hit.point.x, hit.point.y,
            hit.point.z);

    query::Ray short_ray = ray;
    short_ray.maxT = 3.0f;
    report(!query::intersectRay(ball, short_ray, hit), "ray/sphere maxT rejects",
            "hit t=%.4f beyond maxT=3", hit.t);

    query::Ray inside_ray;
    inside_ray.origin = {0, 0, 0};
    inside_ray.direction = {1, 0, 0};
    inside_ray.maxT = 10.0f;
    ok = query::intersectRay(ball, inside_ray, hit);
    report(ok && hit.inside && near(hit.t, 1.0f, 1e-5f), "ray/sphere from inside",
            "t=%.5f inside=%d", hit.t, hit.inside);

    const Shape cube = box({0, 0, 0}, {1, 1, 1});
    ok = query::intersectRay(cube, ray, hit);
    report(ok && near(hit.t, 4.0f, 1e-5f) && near(hit.normal.x, -1.0f, 1e-5f), "ray/box entry",
            "t=%.5f normal=(%.3f,%.3f,%.3f)", hit.t, hit.normal.x, hit.normal.y, hit.normal.z);

    const Shape cyl = cylinder({0, 0, 0}, 1.0f, 2.0f);
    ok = query::intersectRay(cyl, ray, hit);
    report(ok && near(hit.t, 4.0f, 1e-5f) && near(hit.normal.x, -1.0f, 1e-5f),
            "ray/cylinder side", "t=%.5f normal=(%.3f,%.3f,%.3f)", hit.t, hit.normal.x,
            hit.normal.y, hit.normal.z);

    query::Ray cap_ray;
    cap_ray.origin = {0, 0, 5};
    cap_ray.direction = {0, 0, -1};
    cap_ray.maxT = 10.0f;
    ok = query::intersectRay(cyl, cap_ray, hit);
    report(ok && near(hit.t, 3.0f, 1e-5f) && near(hit.normal.z, 1.0f, 1e-5f),
            "ray/cylinder cap", "t=%.5f normal=(%.3f,%.3f,%.3f)", hit.t, hit.normal.x,
            hit.normal.y, hit.normal.z);

    query::Ray miss_ray;
    miss_ray.origin = {0, 2, 0};
    miss_ray.direction = {0, 0, 1};
    miss_ray.maxT = 10.0f;
    report(!query::intersectRay(ball, miss_ray, hit), "ray/sphere miss", "off-axis ray");

    // Rotating the cylinder by 90 degrees about X makes its axis +Y: a ray down +Y now
    // hits the cap at y=+2, not the side.
    const Shape rotated = cylinder({0, 0, 0}, 1.0f, 2.0f, angleAxis(1.57079632679f, float3{1, 0, 0}));
    query::Ray axis_ray;
    axis_ray.origin = {0, 5, 0};
    axis_ray.direction = {0, -1, 0};
    axis_ray.maxT = 10.0f;
    ok = query::intersectRay(rotated, axis_ray, hit);
    report(ok && near(hit.t, 3.0f, 1e-4f), "ray/cylinder rotated cap",
            "t=%.5f (want 3)", hit.t);
}

static void test_closest_and_contains() {
    const Shape cube = box({0, 0, 0}, {1, 2, 3});
    const float3 p = query::closestPoint(cube, float3{2, 1, 0});
    report(near(p.x, 1.0f, 1e-5f) && near(p.y, 1.0f, 1e-5f) && near(p.z, 0.0f, 1e-5f),
            "closest point on box", "p=(%.3f,%.3f,%.3f)", p.x, p.y, p.z);

    const Shape ball = sphere({0, 0, 0}, 2.0f);
    const float3 b = query::closestPoint(ball, float3{10, 0, 0});
    report(near(b.x, 2.0f, 1e-5f), "closest point on sphere", "p=(%.3f,%.3f,%.3f)", b.x, b.y, b.z);

    const Shape cyl = cylinder({0, 0, 0}, 2.0f, 3.0f);
    const float3 c = query::closestPoint(cyl, float3{10, 0, 0});
    report(near(c.x, 2.0f, 1e-5f) && near(c.z, 0.0f, 1e-5f), "closest point on cylinder side",
            "p=(%.3f,%.3f,%.3f)", c.x, c.y, c.z);
    const float3 cap = query::closestPoint(cyl, float3{0, 0, 10});
    report(near(cap.z, 3.0f, 1e-5f), "closest point on cylinder cap", "p=(%.3f,%.3f,%.3f)", cap.x,
            cap.y, cap.z);
    const float3 near_cap = query::closestPoint(cyl, float3{0, 0, 2.9f});
    report(near(near_cap.z, 3.0f, 1e-4f) && near(near_cap.x, 0.0f, 1e-4f),
            "closest point near cylinder cap", "p=(%.3f,%.3f,%.3f)", near_cap.x, near_cap.y,
            near_cap.z);

    report(query::containsPoint(cube, float3{0.5f, 1.5f, 2.5f}), "contains box inside", "inside");
    report(!query::containsPoint(cube, float3{2, 0, 0}), "contains box outside", "outside");
    report(query::containsPoint(ball, float3{1, 0, 0}), "contains sphere inside", "inside");
    report(query::containsPoint(cyl, float3{1, 0, 2}), "contains cylinder inside", "inside");
    report(!query::containsPoint(cyl, float3{1, 0, 4}), "contains cylinder outside cap", "outside");
}

static void test_aabb() {
    const Shape cube = box({1, 2, 3}, {1, 2, 3});
    float3 mn, mx;
    query::shapeAabb(cube, mn, mx);
    report(near(mn.x, 0, 1e-5f) && near(mn.y, 0, 1e-5f) && near(mn.z, 0, 1e-5f)
                    && near(mx.x, 2, 1e-5f) && near(mx.y, 4, 1e-5f) && near(mx.z, 6, 1e-5f),
            "aabb box", "min=(%.2f,%.2f,%.2f) max=(%.2f,%.2f,%.2f)", mn.x, mn.y, mn.z, mx.x, mx.y,
            mx.z);

    // 90 degrees about Z swaps the x/y half extents.
    Shape rotated = box({0, 0, 0}, {1, 2, 0.5f});
    rotated.rotation = angleAxis(1.57079632679f, float3{0, 0, 1});
    query::shapeAabb(rotated, mn, mx);
    report(near(mx.x, 2.0f, 1e-4f) && near(mx.y, 1.0f, 1e-4f) && near(mx.z, 0.5f, 1e-4f),
            "aabb rotated box", "max=(%.3f,%.3f,%.3f)", mx.x, mx.y, mx.z);
}

static void test_sweep() {
    const Shape ball = sphere({-5, 0, 0}, 0.5f);
    const Shape cube = box({0, 0, 0}, {0.5f, 0.5f, 0.5f});

    query::SweepHit hit = query::sweepShape(ball, float3{10, 0, 0}, cube, 0.0f);
    report(hit.hit && hit.safe > 0.35f && hit.safe < 0.42f && hit.unsafe >= hit.safe,
            "sweep sphere into box", "safe=%.4f unsafe=%.4f (want board 0.40)", hit.safe, hit.unsafe);

    query::SweepHit miss = query::sweepShape(ball, float3{0, 10, 0}, cube, 0.0f);
    report(!miss.hit, "sweep miss", "motion perpendicular");

    query::SweepHit margin = query::sweepShape(ball, float3{10, 0, 0}, cube, 0.5f);
    report(margin.hit && margin.safe > 0.30f && margin.safe < 0.37f,
            "sweep margin stops early", "safe=%.4f (want board 0.35)", margin.safe);

    // Fast motion whose endpoint is past the target: the AABB interval must still find it.
    const Shape fast = sphere({0, 20, 0}, 0.2f);
    const Shape target = box({0, 3, 0}, {0.5f, 0.5f, 0.5f});
    query::SweepHit tunnel = query::sweepShape(fast, float3{0, -30, 0}, target, 0.0f);
    report(tunnel.hit && tunnel.safe > 0.50f && tunnel.safe < 0.56f,
            "sweep cannot tunnel", "safe=%.4f (board 0.533)", tunnel.safe);

    // Starting already overlapping: safe fraction is zero.
    const Shape overlapping = sphere({0.5f, 0, 0}, 0.5f);
    query::SweepHit start = query::sweepShape(overlapping, float3{10, 0, 0}, cube, 0.0f);
    report(start.hit && near(start.safe, 0.0f, 1e-6f), "sweep starts in contact", "safe=%.5f",
            start.safe);

    report(std::isfinite(hit.point.x) && std::isfinite(hit.normal.x), "sweep hit data finite",
            "point=(%.3f,%.3f,%.3f) normal=(%.3f,%.3f,%.3f)", hit.point.x, hit.point.y, hit.point.z,
            hit.normal.x, hit.normal.y, hit.normal.z);
}

int main() {
    std::printf("== query geometry\n");
    test_rays();
    test_closest_and_contains();
    test_aabb();
    test_sweep();
    std::printf("\n%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
