#pragma once

namespace avbd {

struct Rigid;
struct Force;

class IslandSleepManager
{
public:
    void update(Rigid *bodies, Force *forces, int sleepFrames);
};

} // namespace avbd
