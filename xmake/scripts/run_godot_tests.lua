-- Run a Godot test script headlessly inside the demo project (which loads the AVBD GDExtension from
-- demo/bin/). `...` forwards extra command-line arguments.
--
-- Exits nonzero (xmake aborts) when Godot exits nonzero.
--
-- The executable is looked up in this order, so a machine with several Godot builds (or one whose
-- binary carries the release name) does not need a PATH edit:
--   1. $AVBD_GODOT / $GODOT - an explicit path, which is what CI and the wrappers set;
--   2. `godot`, `godot4`, `Godot` on PATH - the .exe suffix is added by find_program on Windows;
--   3. the release-style names Godot's own downloads use, e.g. Godot_v4.7-stable_win64.exe.
import("lib.detect.find_program")

local function find_godot()
    for _, var in ipairs({ "AVBD_GODOT", "GODOT" }) do
        local explicit = os.getenv(var)
        if explicit and explicit ~= "" then
            local found = find_program(explicit)
            if found then
                return found
            end
            -- A path that exists but is not on PATH: hand it to find_program as an absolute name.
            if os.isfile(explicit) then
                return explicit
            end
        end
    end

    for _, name in ipairs({ "godot", "godot4", "Godot" }) do
        local found = find_program(name)
        if found then
            return found
        end
    end

    -- Release downloads dropped into the project root: Godot_v4.7-stable_win64.exe on Windows,
    -- Godot_v4.7-stable_linux.x86_64 (or .app bundle) elsewhere.
    for _, pattern in ipairs({ "Godot_v*.exe", "Godot_v*" }) do
        for _, file in ipairs(os.files(path.join(os.projectdir(), pattern)) or {}) do
            return file
        end
    end
    return nil
end

function main(script, ...)
    local godot = find_godot()
    if not godot then
        raise("godot not found; install Godot 4.x (https://godotengine.org/download) or set "
                .. "AVBD_GODOT=/path/to/godot")
    end

    local args = {godot, "--headless", "--path", path.join(os.projectdir(), "demo"), "--script", script}
    for _, a in ipairs({...}) do
        table.insert(args, a)
    end

    print("[godot] %s", script)
    os.execv(godot, args)
end
