#!/bin/bash
# The projection accounting matrix, as one command.
#
# Why it exists: the layer moves a pose, and there are exactly two places that can record the
# move - the constraint's multiplier, or its stabilisation reference. Which of the two a
# constraint kind takes is not a preference, it is a property of the kind: a contact is unilateral
# and clamped (the charge holds the lift and cannot run away), a joint is bilateral and unbounded
# (a charge has no clamp to hold it, and becomes an oscillator). This script runs both ledgers on
# the same scene so the difference is a number instead of an argument.
#
#   ./tools/projection_matrix.sh                # Bridge + Rope + stack, defaults
#   STEPS=300 ./tools/projection_matrix.sh      # shorter
#
# Needs a build with the counters: xmake f -m release --profile=y --diag=y && xmake build avbd_core_test
set -u
cd "$(dirname "$0")/.." || exit 1
BIN=./build/linux/x86_64/release/avbd_core_test
STEPS=${STEPS:-1000}
THREADS=${THREADS:-8}
IT=${IT:-12}

[ -x "$BIN" ] || { echo "no test binary; build it first (see header)"; exit 1; }

# Prints the one line that matters for a scenario: pose extent, worst penetration, cost, and the
# projection counters. `diag` needs --diag=y; without it the second line is simply absent.
row() {
    local label=$1 scene=$2
    shift 2
    printf '%-34s ' "$label"
    env AVBD_ITERATIONS="$IT" "$@" \
        "$BIN" --scene "$scene" --steps "$STEPS" --threads "$THREADS" 2>&1 |
        grep -E 'after|diag' |
        sed -E 's/.*contacts=([0-9]+) max\|x\|=([0-9.]+).*worst=([0-9.]+) +([0-9.]+) ms.*/c=\1 |x|=\2 pen=\3 t=\4ms/; s/^  diag: maxJoint\|lambda\|=([0-9.e+]+).*charge=([0-9.]+\/[0-9.]+) asked=([0-9.]+\/[0-9.]+).*/joineff lam=\1 charge=\2 asked=\3/' |
        paste -sd' '
}

echo "== adaptive iterations off, $IT rounds/substep, threads=$THREADS, $STEPS steps =="
echo
echo "-- jointed: the ledger decides whether the chain is an oscillator --"
row "Bridge  baseline (no projection)"      Bridge
row "Bridge  joints: charge (rebase=0)"     Bridge AVBD_PROJECTION_ROUNDS=2 AVBD_PROJECTION_REBASE=0
row "Bridge  joints: rebase (default)"      Bridge AVBD_PROJECTION_ROUNDS=2
row "Rope    baseline (no projection)"      Rope
row "Rope    joints: charge (rebase=0)"     Rope   AVBD_PROJECTION_ROUNDS=2 AVBD_PROJECTION_REBASE=0
row "Rope    joints: rebase (default)"      Rope   AVBD_PROJECTION_ROUNDS=2
echo
echo "-- with joint_relaxation defaulting to 0.5 the chain holds without any projection --"
row "Rope    baseline (no projection)"      Rope
for jr in 0.7 0.6 0.5 0.4; do
    row "Rope    joint_relaxation=$jr" Rope AVBD_JOINT_RELAX="$jr"
done
echo
echo "-- contacts: the charge is what holds the lift (re-anchoring them costs 4x penetration) --"
row "Stack   baseline (no projection)"      Stack
row "Stack   contacts: charge (default)"    Stack AVBD_PROJECTION_ROUNDS=2
echo
echo "-- the same on the loaded bridge: the chain's damping is a per-kind constant --"
for jr in 0.7 0.5 0.4; do
    row "Bridge  joint_relaxation=$jr" Bridge AVBD_PROJECTION_ROUNDS=2 AVBD_JOINT_RELAX="$jr"
done
echo
echo "-- the dual-only tail has no fixed point: newton_ratio < 1 winds the multiplier up --"
row "Pyramid nr=1.0 decay=1.0 (default)" Pyramid AVBD_NEWTON_RATIO=1.0 AVBD_STIFFNESS_DECAY=1.0
row "Pyramid nr=0.5 decay=1.0"           Pyramid AVBD_NEWTON_RATIO=0.5 AVBD_STIFFNESS_DECAY=1.0
row "Pyramid nr=0.5 decay=0.2 (old demo)" Pyramid AVBD_NEWTON_RATIO=0.5 AVBD_STIFFNESS_DECAY=0.2
echo
echo "-- round budget is the host's number now; projection buys rounds back --"
for it in 8 12 16; do
    row "Stack   iterations=$it, projection off" Stack AVBD_ITERATIONS="$it"
    row "Stack   iterations=$it, projection"     Stack AVBD_ITERATIONS="$it" AVBD_PROJECTION_ROUNDS=2
done
