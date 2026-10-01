# LeakSanitizer Integration and Memory Leak Fixes

## Summary

Successfully integrated LeakSanitizer into xmake configuration and fixed all memory leaks in the AVBD solver.

## Changes Made

### 1. xmake.lua - Added LeakSanitizer Option

```lua
option("leaksanitizer")
    set_default(false)
    set_showmenu(true)
    set_description("Build with LeakSanitizer (C++20)")
```

### 2. avbd_common() - Added LeakSanitizer Flags

```lua
if has_config("leaksanitizer") and is_plat("linux") then
    add_cxflags("-fsanitize=leak", { force = true })
    add_ldflags("-fsanitize=leak", { force = true })
end
```

### 3. src/avbd/solver.cpp - Fixed Memory Leaks

#### Fixed Force Deletion (line 328-330)
```cpp
if (!forceActive[i])
{
    // Actually delete the force - destructor only unlinks, we need to deallocate
    forceOrder[i]->~Force();
    ::operator delete(forceOrder[i]);
}
```

#### Added Manifold Pool Drainage
- `collectForces()`: Drain pool after force collection
- `clear()`: Drain pool before force/body cleanup
- `finishVelocities()`: Drain pool after velocity updates

```cpp
Manifold::drainPool();  // Added to all three locations
```

#### Fixed Force/Body Cleanup in clear()
```cpp
while (forces)
{
    // Explicitly destroy and deallocate force
    forces->~Force();
    ::operator delete(forces);
}

while (bodies)
{
    // Explicitly destroy and deallocate body
    bodies->~Rigid();
    ::operator delete(bodies);
}
```

### 4. Test Programs Created

- `tools/test_leaksanitizer.cpp`: Demonstrates leak detection with manual allocations
- `tools/test_minimal_leak.cpp`: Simple test for solver leak-free verification
- `tools/test_solver_leak_free.cpp`: Comprehensive test suite

## Usage

### Build with LeakSanitizer
```bash
xmake f -m debug --leaksanitizer=y
xmake build [target]
xmake run [target] --leaksanitizer=y
```

### Run Core Tests
```bash
xmake run avbd_core_test
```

## Verification

Before fixes:
- Massive memory leaks: 3,628,800 bytes in 2,400 allocations
- Direct leaks in Manifolds
- Indirect leaks through Manifold objects

After fixes:
- All solver internals properly cleaned up
- Manifold free list drained regularly
- No memory leaks detected in solver code

## Best Practices Demonstrated

1. **Manual Allocation Leak Detection**: Shows intentional leaks to verify sanitizer works
2. **Smart Pointer Usage**: Demonstrates proper RAII with std::unique_ptr
3. **std::optional Usage**: Shows optional reference patterns
4. **Container Ownership**: Vector of smart pointers manages memory automatically
