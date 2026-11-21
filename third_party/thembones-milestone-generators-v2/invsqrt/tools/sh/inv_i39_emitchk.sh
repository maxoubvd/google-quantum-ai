#!/usr/bin/env bash
# I39 (third holder), LOCAL WSL: emit a hold-route variant (the env of tools/sh/inv_i39_s17hold.sh) and check it on
# ALL 3,133,452 inputs against a base emit (GHOST's r63_inv_ghost_chk.py check: outcomes 0 / 1 / two seeds; x
# restored, no dirty wire, phase equal to the base's on every input).
#   bash tools/sh/inv_i39_emitchk.sh NAME SRC CFG CORES [BASE.kmx]     (default base: Egq6's raw emit)
# Writes ~/qcc/ws/inv_i39/hold/files18/NAME/.  Safe to run: other seats' tools are only run.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/inv_i39/hold/files18
NAME=$1; S=$2; CFG=$3; CORES=$4; BASE=${5:-~/qcc/ws/inv_i39/hold/files17/Egq6/invsqrt.kmx}
O=$F/$NAME; mkdir -p $O; rm -f $O/invsqrt.kmx
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=7,3,0;9,3,1" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/i39nb_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== emit $(date -u +%H:%M:%SZ)"
(cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S $O $O/angles.json $O/cfg.json 2>&1 \
  | grep -E "emit |Traceback|Error|assert" | cut -c1-260)
[ -f $O/invsqrt.kmx ] || { echo "EMIT FAILED"; exit 1; }
echo "== all-input check against $BASE $(date -u +%H:%M:%SZ)"
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $BASE $O/invsqrt.kmx 0 2>&1 | tail -6 | cut -c1-220)
echo "== done $(date -u +%H:%M:%SZ)"
