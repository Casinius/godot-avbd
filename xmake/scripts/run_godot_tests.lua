-- Run a Godot test script headlessly inside the demo project (which loads the
-- AVBD GDExtension from demo/bin/). `...` forwards extra command-line arguments.
--
-- Exits nonzero (xmake aborts) when Godot exits nonzero.
import("lib.detect.find_program")

function main(script, ...)
    local godot = find_program("godot")
    if not godot then
        raise("godot not found on PATH; install Godot 4.x (https://godotengine.org/download)")
    end

    local args = {godot, "--headless", "--path", path.join(os.projectdir(), "demo"), "--script", script}
    for _, a in ipairs({...}) do
        table.insert(args, a)
    end

    print("[godot] %s", script)
    os.execv(godot, args)
end
