#include "avbd/island_sleep.hpp"

#include "avbd/solver.h"

#include <algorithm>
#include <climits>
#include <vector>

namespace avbd {

void IslandSleepManager::update(Rigid *bodies, Force *forces, int sleepFrames)
{
    std::vector<Rigid *> islandBodies;
    for (Rigid *body = bodies; body != 0; body = body->next)
        if (body->mass > 0 && body->sleep_mode == 1)
            islandBodies.push_back(body);

    const int count = static_cast<int>(islandBodies.size());
    if (count < 2)
        return;

    std::vector<int> parent(count);
    for (int i = 0; i < count; ++i)
        parent[i] = i;
    auto find = [&parent](int x) {
        while (parent[x] != x)
            x = parent[x] = parent[parent[x]];
        return x;
    };

    std::vector<std::pair<Rigid *, int>> sorted;
    sorted.reserve(islandBodies.size());
    for (int i = 0; i < count; ++i)
        sorted.emplace_back(islandBodies[i], i);
    std::sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        return a.first < b.first;
    });
    std::vector<std::pair<int, int>> forceSlots;
    for (Force *f = forces; f != 0; f = f->next)
    {
        const auto slotOf = [&sorted](Rigid *b) -> int {
            const auto it = std::lower_bound(sorted.begin(), sorted.end(), b,
                    [](const auto &entry, const Rigid *body) { return entry.first < body; });
            return (it != sorted.end() && it->first == b) ? it->second : -1;
        };
        const int ia = f->bodyA ? slotOf(f->bodyA.get()) : -1;
        const int ib = f->bodyB ? slotOf(f->bodyB.get()) : -1;
        if (ia >= 0 && ib >= 0)
            forceSlots.emplace_back(ia, ib);
    }

    for (const auto &[ia, ib] : forceSlots)
        parent[find(ia)] = find(ib);

    std::vector<char> islandAwake(count, 0);
    std::vector<int> islandStill(count, INT_MAX);
    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (!islandBodies[i]->sleeping)
            islandAwake[r] = 1;
        islandStill[r] = std::min(islandStill[r], islandBodies[i]->stillFrames);
    }

    for (int i = 0; i < count; ++i)
    {
        const int r = find(i);
        if (islandAwake[r] || islandStill[r] < sleepFrames)
            continue;
        Rigid *body = islandBodies[i];
        if (!body->sleeping)
        {
            body->sleeping = true;
            body->velocityLin = float3{0, 0, 0};
            body->velocityAng = float3{0, 0, 0};
        }
    }
}

} // namespace avbd
