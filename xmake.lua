-- AVBD physics engine for Godot 4, exposed as a GDExtension.
--
-- One project script for every platform. The only things that differ between them are the
-- toolchain's flag spelling, the system library that carries the threading primitives, and the name
-- Godot looks the library up under; all three are resolved from the configured platform/toolchain, so
-- a Windows checkout configures with the same command a Linux one does:
--
--   xmake f -m release && xmake                      # host platform, default toolchain, release
--   xmake f -p linux   -a x86_64 -m release && xmake  # pinned: linux
--   xmake f -p windows -a x64    -m release && xmake  # pinned: windows
--
-- Builds against the `godotcpp4` package from xmake-repo (godot-cpp C++ bindings) and
-- `thread-pool` for the job pool behind the parallel solver update.
--
-- Sanitizers (each needs its own build directory, so configure before building):
--   xmake f -m debug --asan=y --ubsan=y && xmake
--   xmake f -m debug --tsan=y && xmake
-- MSVC supports AddressSanitizer only; the other two are reported as unsupported instead of being
-- silently dropped.
--
-- Targets:
--   avbd            GDExtension shared library (loaded by the demo project / any game)
--   avbd_core_test  standalone solver test binary (no Godot involved)

add_rules("mode.debug", "mode.release")
add_rules("plugin.compile_commands.autoupdate", { outputdir = "." })

add_requires("godotcpp4 4.5")
add_requires("unordered_dense v4.8.1")
-- BS::thread_pool: the job pool behind Solver's parallel update. Header-only, MIT, C++17.
add_requires("thread-pool v5.1.0")
local avbd_core = {
    "src/avbd/rigid.cpp",
    "src/avbd/spatial_hash.cpp",
    "src/avbd/bvh_broadphase.cpp",
    "src/avbd/force.cpp",
    "src/avbd/joint.cpp",
    "src/avbd/generic_joint.cpp",
    "src/avbd/spring.cpp",
    "src/avbd/query/query_shape.cpp",
    "src/avbd/manifold.cpp",
    "src/avbd/collide.cpp",
    "src/avbd/solver.cpp",
}

option("asan")
    set_default(false)
    set_showmenu(true)
    set_description("Build with AddressSanitizer")
option("ubsan")
    set_default(false)
    set_showmenu(true)
    set_description("Build with UndefinedBehaviorSanitizer (GCC/Clang only)")
option("tsan")
    set_default(false)
    set_showmenu(true)
    set_description("Build with ThreadSanitizer (GCC/Clang only)")

-- Warnings are errors of intent: the tree builds clean at this level, and a new warning
-- should be dealt with rather than accumulated.
local avbd_warnings = { "all", "extra" }

-- The name Godot's `.gdextension` manifest looks for: lib<name>.<platform>.<arch>.<ext>. This has to
-- agree with demo/gdavbd/avbd.gdextension - add the pair there when adding a platform here.
local function godot_library_name()
    local platform = is_plat("windows") and "windows"
            or is_plat("macosx") and "macos"
            or is_plat("linux") and "linux"
            or is_plat()
    local arch = is_arch("x86_64", "x64") and "x86_64"
            or is_arch("arm64", "aarch64") and "arm64"
            or is_arch("x86", "i386") and "x86_32"
            or is_arch("arm") and "arm32"
            or is_arch()
    local ext = is_plat("windows") and ".dll" or (is_plat("macosx") and ".dylib" or ".so")
    return "libavbd." .. platform .. "." .. arch .. ext
end

-- Settings shared by both targets. Called from inside each target block, so it writes into
-- that target's scope.
local function avbd_common()
    set_languages("c++23")
    set_symbols("debug")
    set_strip("none")
    set_warnings(avbd_warnings)
end

-- Everything that depends on the compiler rather than on the project. Runs in `on_config`, not
-- `on_load`: the target's tool is not resolved yet when it is loaded.
local function avbd_toolchain(target)
    if target:has_tool("cxx", "msvc") then
        -- /W4 comes from `set_warnings`; /permissive- is the conforming parser the tree is written
        -- against. Cross-platform bit-exact digests are a hard constraint, so the floating-point
        -- mode is spelled out instead of left to the toolchain default: /fp:precise keeps
        -- expressions in source order and /fp:contract- keeps `a*b + c` from collapsing into an fma.
        -- NOMINMAX keeps windows.h's min/max macros out of std::min/std::max.
        target:add("cxflags", "/permissive-", "/fp:precise", "/fp:contract-", { force = true })
        target:add("defines", "NOMINMAX")
    else
        target:add("cxflags", "-Wshadow", "-Wnon-virtual-dtor")
        -- Floating-point contraction is the one default that silently breaks the digest guarantee:
        -- a target with an fma turns `a*b + c` into one rounded operation where a target without one
        -- keeps two, so the same scene lands on different bits per platform. Off everywhere,
        -- explicitly, so a machine's ISA cannot change an answer. No -ffast-math for the same reason.
        target:add("cxflags", "-ffp-contract=off", "-fno-fast-math")
    end

    -- Threading primitives. Every platform but Windows gets them from the C library; Windows has no
    -- separate library to link for std::thread.
    if not is_plat("windows") then
        target:add("syslinks", "pthread")
    end

    if has_config("asan") then
        if target:has_tool("cxx", "msvc") then
            target:add("cxflags", "/fsanitize=address", { force = true })
            target:add("ldflags", "/fsanitize=address", { force = true })
        else
            target:add("cxflags", "-fsanitize=address", { force = true })
            target:add("ldflags", "-fsanitize=address", { force = true })
        end
    end
    for _, kind in ipairs({ "ubsan", "tsan" }) do
        if has_config(kind) then
            if target:has_tool("cxx", "msvc") then
                wprint("--%s is not supported by MSVC; ignoring it", kind)
            else
                target:add("cxflags", "-fsanitize=" .. (kind == "ubsan" and "undefined" or "thread"),
                        { force = true })
                target:add("ldflags", "-fsanitize=" .. (kind == "ubsan" and "undefined" or "thread"),
                        { force = true })
            end
        end
    end
end

-- The GDExtension itself: core solver + Godot node layer.
target("avbd")
    set_kind("shared")
    set_languages("c++23")
    set_optimize("fastest")

    set_symbols("debug")
    set_strip("none")
    add_files(table.join(avbd_core, "src/nodes/*.cpp", "src/server/*.cpp"))
    add_includedirs("src")
    add_packages("godotcpp4", "thread-pool","unordered_dense")
    add_defines("BS_THREAD_POOL_DISABLE_EXCEPTION_HANDLING")
    set_exceptions("none")
    avbd_common()
    on_config(avbd_toolchain)
    -- Godot loads the library from the project's bin/ directory, under the name its .gdextension
    -- manifest lists. Both sides derive it from the same platform/arch spelling.
    after_build(function (target)
        local dst = path.join(os.projectdir(), "demo", "gdavbd", "bin", godot_library_name())
        print("copying library to " .. dst)
        os.mkdir(path.directory(dst))
        os.cp(target:targetfile(), dst)
    end)

-- Solver-only tests: no engine, no windowing, direct numeric assertions.
target("avbd_core_test")
    set_kind("binary")
    add_files(table.join(avbd_core, "test/core_test.cpp", "tools/core_scenes.cpp"))
    add_includedirs("src", "tools")
    add_packages("thread-pool")
    avbd_common()
    on_config(avbd_toolchain)

-- Query-extension geometry tests: no engine, no solver step.
target("avbd_query_test")
    set_kind("binary")
    add_files(table.join(avbd_core, "test/query_test.cpp"))
    add_includedirs("src")
    add_packages("thread-pool")
    avbd_common()
    on_config(avbd_toolchain)

-- Constraint composition tests (Godot headless, inside the demo project).
target("avbd_constraint_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_constraints.gd")
    end)

-- Solver behaviour tests (Godot headless, inside the demo project).
target("avbd_solver_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_solver.gd")
    end)

-- Space-query behaviour (Godot headless, inside the demo project).
target("avbd_query_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_queries.gd")
    end)

-- API compliance audit (Godot headless, inside the demo project).
target("check_api")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_server.gd")
    end)

-- Load-gate smoke test (Godot headless, inside the demo project).
target("avbd_load_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_load.gd")
    end)

-- Double-pyramid stress regression (Godot headless): 420 boxes must settle
-- without launched bodies or long-range drift.
target("avbd_pyramid_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_pyramid.gd")
    end)

-- Long-run penetration detection for the double pyramid (Godot headless). Reads the
-- solver's own contact error to measure real interpenetration over a long wall of time,
-- and flags penetration that happens while secondary (re-)collisions resolve. Slower than
-- the pyramid regression: pass `-- --pen-seconds=N` to bound the run.
target("avbd_penetration_tests")
    set_kind("phony")
    on_run(function (target)
        import("scripts.run_godot_tests", {rootdir = "xmake"}).main("res://tests/test_penetration.gd")
    end)

-- Run every tier: C++ solver tests, then all Godot suites.
target("run_all_tests")
    set_kind("phony")
    on_run(function (target)
        local runner = import("scripts.run_godot_tests", {rootdir = "xmake"})
        os.exec("xmake run avbd_core_test")
        os.exec("xmake run avbd_query_test")
        runner.main("res://tests/test_solver.gd")
        runner.main("res://tests/test_queries.gd")
        runner.main("res://tests/test_constraints.gd")
        runner.main("res://tests/test_server.gd")
        runner.main("res://tests/test_pyramid.gd")
        runner.main("res://tests/test_load.gd")
    end)
