-- AVBD physics engine for Godot 4, exposed as a GDExtension.
--
-- Built against the `godotcpp4` package from xmake-repo (godot-cpp C++ bindings) and
-- `thread-pool` for the job pool behind the parallel solver update.
--
--   xmake f -m release && xmake
--   xmake run avbd_core_test
--
-- Sanitizers (each needs its own build directory, so configure before building):
--   xmake f -m debug --asan=y --ubsan=y && xmake
--   xmake f -m debug --tsan=y && xmake
--
-- Targets:
--   avbd            GDExtension shared library (loaded by the demo project / any game)
--   avbd_core_test  standalone solver test binary (no Godot involved)

add_rules("mode.debug", "mode.release")
add_rules("plugin.compile_commands.autoupdate", { outputdir = "." })

add_requires("godotcpp4 4.1")
-- BS::thread_pool: the job pool behind Solver's parallel update. Header-only, MIT, C++17.
add_requires("thread-pool v5.1.0", {configs = {headeronly = true}})
add_requires("eigen 5.0.1")
local avbd_core = {
    "src/avbd/rigid.cpp",
    "src/avbd/bvh/node_storage.cpp",
    "src/avbd/bvh/bvh.cpp",
    "src/avbd/force.cpp",
    "src/avbd/joint.cpp",
    "src/avbd/generic_joint.cpp",
    "src/avbd/spring.cpp",
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
    set_description("Build with UndefinedBehaviorSanitizer")
option("tsan")
    set_default(false)
    set_showmenu(true)
    set_description("Build with ThreadSanitizer")

-- Warnings are errors of intent: the tree builds clean at this level, and a new warning
-- should be dealt with rather than accumulated.
-- xmake's own warning presets stop at "extra"; the two clang/gcc flags the tree is also
-- clean at are added explicitly.
local avbd_warnings = { "all", "extra" }

-- Settings shared by both targets. Called from inside each target block, so it writes into
-- that target's scope.
local function avbd_common()
    set_languages("c++23")
    set_warnings(avbd_warnings)
    add_cxflags("-Wshadow", "-Wnon-virtual-dtor")
    add_syslinks("pthread")
    if has_config("asan") and is_plat("linux") then
        add_cxflags("-fsanitize=address", { force = true })
        add_ldflags("-fsanitize=address", { force = true })
    end
    if has_config("ubsan") and is_plat("linux") then
        add_cxflags("-fsanitize=undefined", { force = true })
        add_ldflags("-fsanitize=undefined", { force = true })
    end
    if has_config("tsan") and is_plat("linux") then
        add_cxflags("-fsanitize=thread", { force = true })
        add_ldflags("-fsanitize=thread", { force = true })
    end
end

-- The GDExtension itself: core solver + Godot node layer.
target("avbd")
    avbd_common()
    set_kind("shared")
    set_languages("c++23")
    set_symbols("debug")
    set_strip("none")
    add_files(table.join(avbd_core, "src/nodes/*.cpp", "src/server/*.cpp"))
    add_includedirs("src")
    add_packages("godotcpp4", "thread-pool","eigen")
    add_defines("BS_THREAD_POOL_DISABLE_EXCEPTION_HANDLING")
    set_exceptions("none")
    -- Godot loads the library from the project's bin/ directory.
    after_build(function (target)
        local ext = is_plat("windows") and ".dll" or (is_plat("macosx") and ".dylib" or ".so")
        local dst = path.join(os.projectdir(), "demo","gdavbd", "bin", "libavbd.linux.x86_64" .. ext)
        os.mkdir(path.directory(dst))
        os.cp(target:targetfile(), dst)
    end)

-- Solver-only tests: no engine, no windowing, direct numeric assertions.
target("avbd_core_test")
    avbd_common()
    set_kind("binary")
    add_files(table.join(avbd_core, "test/core_test.cpp", "tools/core_scenes.cpp"))
    add_includedirs("src", "tools")
    add_packages("thread-pool","eigen")

-- Parallelization deterministic test: verifies bit-identical results across thread counts.
target("test_parallel")
    avbd_common()
    set_kind("binary")
    add_files("tools/test_parallel.cpp")
    add_includedirs("src")
    add_packages("thread-pool","eigen")
    add_deps("avbd")

-- Simple parallelization test: verifies parallel code compiles and runs without errors.
target("test_parallel_simple")
    avbd_common()
    set_kind("binary")
    add_files("tools/test_parallel_simple.cpp")
    add_includedirs("src")
    add_packages("thread-pool","eigen")
    add_deps("avbd")

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

-- Run every tier: C++ solver tests, then all Godot suites.
target("run_all_tests")
    set_kind("phony")
    on_run(function (target)
        local runner = import("scripts.run_godot_tests", {rootdir = "xmake"})
        os.exec("xmake run avbd_core_test")
        runner.main("res://tests/test_solver.gd")
        runner.main("res://tests/test_constraints.gd")
        runner.main("res://tests/test_server.gd")
        runner.main("res://tests/test_load.gd")
    end)
