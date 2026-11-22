#!/usr/bin/env bash
#
# Refactor safety gate. Runs the solver-free test binaries and compares the
# state digests they print with tools/baseline/digests.txt. Bit-identical
# digests are the contract: a refactor must move code, not change float
# expression order or loop order.
#
#   tools/verify_refactor.sh            check against the recorded baseline
#   tools/verify_refactor.sh --update   re-record the baseline (behaviour change)
#
# The Godot-side suites (avbd_query_tests, avbd_solver_tests, ...) need a Godot
# binary and are not covered here.
set -euo pipefail

cd "$(dirname "$0")/.."

update=0
if [[ "${1:-}" == "--update" ]]; then
    update=1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

fail=0

run_suite() {
    local target=$1 log=$2
    if ! xmake run "$target" >"$log" 2>&1; then
        echo "FAIL $target exited non-zero (see $log)"
        tail -20 "$log"
        fail=1
        return
    fi
    local summary
    summary=$(grep -E '^[0-9]+ checks, [0-9]+ failures$' "$log" | tail -1 || true)
    if [[ -z "$summary" ]]; then
        echo "FAIL $target printed no check summary (see $log)"
        fail=1
        return
    fi
    echo "     $target: $summary"
    if [[ "$summary" != *", 0 failures" ]]; then
        echo "FAIL $target reported failures:"
        grep -E '^FAIL ' "$log" || true
        fail=1
    fi
}

run_suite avbd_core_test "$tmp/core.txt"
run_suite avbd_query_test "$tmp/query.txt"

# The core test prints one digest per scenario; cross-thread equality is
# asserted in-process, so the serial line is the full record.
grep -E '^ +digest ' "$tmp/core.txt" | sed 's/^ *//' >"$tmp/digests.txt" || true

if [[ ! -s "$tmp/digests.txt" ]]; then
    echo "FAIL no digest lines in avbd_core_test output"
    fail=1
elif [[ "$update" == 1 ]]; then
    cp "$tmp/digests.txt" tools/baseline/digests.txt
    echo "     baseline updated:"
    sed 's/^/       /' tools/baseline/digests.txt
elif ! diff -u tools/baseline/digests.txt "$tmp/digests.txt"; then
    echo "FAIL digests changed: the refactor altered simulation state"
    fail=1
else
    echo "     digests match baseline"
fi

if [[ "$fail" != 0 ]]; then
    exit 1
fi

echo "OK refactor gate passed"
