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
-- There is no second project script and no per-platform wrapper: xmake reads exactly one project
-- script (XMAKE_PROJECT_FILE can rename it, but every target would have to be defined twice), and the
-- only thing that differs per platform is the configure command above - which is xmake's own
-- interface, not a script layered on top of it.
--
-- Builds against the `godotcpp4` package from xmake-repo (godot-cpp C++ bindings) and `thread-pool`
-- for the job pool behind the parallel solver update.
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
-- BS::thread_pool: the job pool behind Solver's parallel update. Header-only, MIT, C++17.
add_requires("thread-pool v5.1.0")
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
    set_description("Build with UndefinedBehaviorSanitizer (GCC/Clang only)")
option("tsan")
    set_default(false)
    set_showmenu(true)
    set_description("Build with ThreadSanitizer (GCC/Clang only)")
-- Instrumentation. `--profile=y` fills Solver::StepProfile with per-phase wall-clock totals
-- (worker 0: the thread that waits when a phase is imbalanced); `--diag=y` fills Solver::Diag with
-- the projection counters. Both are compile-time: the counters sit in the per-contact hot path and an
-- always-on atomic there would change the thing being measured.
-- Note the name: xmake has a global `--profile` option of its own (saved configuration profiles), and
-- a project option of the same name shadows it. Nothing in this project uses those, but a script that
-- wants one can spell it as `-c --profile=...` for xmake's and read this one from the project.
option("profile")
    set_default(false)
    set_showmenu(true)
    set_description("Per-phase timing (AVBD_PROFILE)")
option("diag")
    set_default(false)
    set_showmenu(true)
    set_description("Projection counters (AVBD_DIAG)")

-- Warnings are errors of intent: the tree builds clean at this level, and a new warning should be
-- dealt with rather than accumulated.
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

-- Settings shared by both targets. Called from inside each target block, so it writes into that
-- target's scope.
local function avbd_common()
    set_languages("c++23")
    set_symbols("debug")
    set_strip("none")
    set_warnings(avbd_warnings)
    if has_config("profile") then
        add_defines("AVBD_PROFILE")
    end
    if has_config("diag") then
        add_defines("AVBD_DIAG")
    end
end

-- Everything that depends on the compiler rather than on the project. Runs in `on_config`, not
-- `on_load`: the target's tool is not resolved yet when it is loaded (xmake refuses to answer
-- `has_tool` there), and this is the hook it points at.
local function avbd_toolchain(target)
    if target:has_tool("cxx", "msvc") then
        -- /W4 comes from `set_warnings`; /permissive- is the conforming parser the tree is written
        -- against. Cross-platform bit-exact digests are a hard constraint (see the GCC branch), so the
        -- floating-point mode is spelled out instead of left to the toolchain default: /fp:precise is
        -- the mode that keeps expressions in source order and /fp:contract- keeps `a*b + c` from
        -- collapsing into an fma. NOMINMAX keeps windows.h's min/max macros out of
        -- std::min/std::max. The math constants come from <numbers>/<limits>, so the MSVC
        -- M_PI workaround (_USE_MATH_DEFINES) is not needed.
        target:add("cxflags", "/permissive-", "/fp:precise", "/fp:contract-", { force = true })
        target:add("defines", "NOMINMAX")
    else
        target:add("cxflags", "-Wshadow", "-Wnon-virtual-dtor")
        -- Floating-point contraction is the one default that silently breaks the digest guarantee: a
        -- target with an fma turns `a*b + c` into one rounded operation where a target without one keeps
        -- two, so the same scene lands on different bits per platform. Off everywhere, explicitly, so a
        -- machine's ISA cannot change an answer. No -ffast-math for the same reason: no reassociation,
        -- no reciprocal swaps.
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
    add_files(table.join(avbd_core, "src/nodes/*.cpp", "src/server/*.cpp"))
    add_includedirs("src")
    add_packages("godotcpp4", "thread-pool", "eigen")
    add_defines("BS_THREAD_POOL_DISABLE_EXCEPTION_HANDLING")
    set_exceptions("none")
    avbd_common()
    on_config(avbd_toolchain)
    -- Godot loads the library from the project's bin/ directory, under the name its .gdextension
    -- manifest lists. Both sides derive it from the same platform/arch spelling.
    after_build(function (target)
        local dst = path.join(os.projectdir(), "demo", "gdavbd", "bin", godot_library_name())
        os.mkdir(path.directory(dst))
        os.cp(target:targetfile(), dst)
    end)

-- Solver-only tests: no engine, no windowing, direct numeric assertions.
target("avbd_core_test")
    set_kind("binary")
    add_files(table.join(avbd_core, "test/core_test.cpp", "tools/core_scenes.cpp"))
    add_includedirs("src", "tools")
    add_packages("thread-pool", "eigen")
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

-- The projection accounting matrix, as one command: both ledgers (multiplier charge vs stabilisation
-- re-anchor) on the same scenes, so which one a constraint kind takes is a number instead of an
-- argument. Configure with --profile=y --diag=y first - the interesting part of the output is the
-- `stage:`/`diag:`/`broad:` line each run prints. AVBD_MATRIX_STEPS overrides the tick count for a
-- quick pass.
target("projection_matrix")
    set_kind("phony")
    on_run(function (target)
        local project = import("core.project.project")
        local bin = project.target("avbd_core_test"):targetfile()
        if not os.isfile(bin) then
            raise("%s is not built; run `xmake` first", bin)
        end
        local steps = os.getenv("AVBD_MATRIX_STEPS") or "1000"
        local threads = os.getenv("AVBD_MATRIX_THREADS") or "8"
        local rounds = os.getenv("AVBD_MATRIX_ROUNDS") or "8"

        local function row(label, scene, envs)
            print("\n== %s", label)
            -- try = true: one scene that reports non-finite state should not abort the matrix, the
            -- exit code is printed instead.
            local code = os.execv(bin, { "--scene", scene, "--steps", steps, "--threads", threads },
                    { envs = envs, try = true })
            if code ~= 0 then
                print("   (exit %s)", tostring(code))
            end
        end

        row("jointed: the ledger decides whether the chain is an oscillator - baseline",
                "Bridge", { AVBD_ITERATIONS = rounds })
        row("jointed: joints charge the multiplier (old behaviour), Bridge", "Bridge",
                { AVBD_ITERATIONS = rounds, AVBD_PROJECTION_ROUNDS = "2", AVBD_PROJECTION_REBASE = "0" })
        row("jointed: joints re-anchor the reference (default), Bridge", "Bridge",
                { AVBD_ITERATIONS = rounds, AVBD_PROJECTION_ROUNDS = "2" })
        row("jointed: pure chain - baseline, Rope", "Rope", { AVBD_ITERATIONS = rounds })
        row("jointed: pure chain - rebase, Rope", "Rope",
                { AVBD_ITERATIONS = rounds, AVBD_PROJECTION_ROUNDS = "2" })

        row("damping: contact-only bodies keep 0.7 (no projection), Rope", "Rope",
                { AVBD_ITERATIONS = rounds, AVBD_JOINT_RELAX = "0.7" })
        row("damping: jointed bodies 0.5 (default), Rope", "Rope", { AVBD_ITERATIONS = rounds })

        row("contacts: the charge is what holds the lift - baseline, Stack", "Stack",
                { AVBD_ITERATIONS = rounds })
        row("contacts: charge, Stack", "Stack", { AVBD_ITERATIONS = rounds, AVBD_PROJECTION_ROUNDS = "2" })

        row("round budget: 6 rounds + projection, Pyramid XL", "Pyramid XL",
                { AVBD_ITERATIONS = "6", AVBD_PROJECTION_ROUNDS = "2" })
        row("round budget: 8 rounds + projection, Pyramid XL", "Pyramid XL",
                { AVBD_ITERATIONS = "8", AVBD_PROJECTION_ROUNDS = "2" })
        row("round budget: 12 rounds, no projection, Pyramid XL", "Pyramid XL", { AVBD_ITERATIONS = "12" })

        row("broad phase: incremental pairs off, Pyramid XL", "Pyramid XL",
                { AVBD_ITERATIONS = rounds, AVBD_INC_EPS = "0" })
        row("broad phase: incremental pairs, 1 mm drift, Pyramid XL", "Pyramid XL",
                { AVBD_ITERATIONS = rounds })
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
