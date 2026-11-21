#!/usr/bin/env bash
# FC eighteenth holder, LOCAL WSL: EMIT ONLY (count-only probe) of a hold-route variant with R64's xcrun (per-tag caps
# from env I39_CAPS), r59b0's router and angles.  Prints the emit line, per-tag caps and any error.
#   [I39_CAPS=H1=38,...] bash tools/sh/invsqrt_fc_s18emit.sh NAME SRC CFG CORES
# Writes ~/qcc/ws/invsqrt_fc/s18/files/NAME/ only.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/invsqrt_fc/s18/files
NAME=$1; S=$2; CFG=$3; CORES=$4
O=$F/$NAME
mkdir -p $O; rm -f $O/invsqrt.kmb $O/invsqrt.kmx
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=7,3,0;9,3,1" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/fc18nb_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== emit $NAME caps=${I39_CAPS:-none} $(date -u +%H:%M:%SZ)"
(cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S $O $O/angles.json $O/cfg.json 2>&1 \
  | grep -E "R55 net|emit |per-tag|r63:|Traceback|Error|assert|File \"" | cut -c1-300 | tail -14)
[ -f $O/invsqrt.kmx ] && echo "EMIT OK $(date -u +%H:%M:%SZ)" || echo "EMIT FAILED (no .kmx) $(date -u +%H:%M:%SZ)"
