#!/usr/bin/env bash
# Workstream invsqrt_pp (P1), local WSL: the exact AFFINE pass on one invsqrt circuit, one command.
#   bash tools/sh/invsqrt_pp_affpass.sh <circuit invsqrt.kmb> <fresh out dir> [cores=0,1]
# 1 audit: every Toffoli-class op against the exact affine span of the wires just before it, all 3,133,452 inputs
#   (tools/invsqrt_pp_affine.py -> <out>/affine.json)
# 2 rewrite: every top-level site replaced by CNOT / X gates (or HMR + conditioned Clifford repair when the target is
#   in the relation), then the child checked on all inputs (state, Z-phase on 4 transcripts, every Z_POW value)
#   (tools/invsqrt_pp_affwrite.py -> <out>/invsqrt.kmb, <out>/affwrite.log)
# Then verify: bash tools/sh/box_verify.sh <out> invsqrt <name>
set -u
IN=$(readlink -f "${1:?circuit.kmb}"); OUT=${2:?out dir}; CORES=${3:-0,1}
T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
[ -e "$OUT/invsqrt.kmb" ] && { echo "$OUT/invsqrt.kmb exists: give a fresh folder"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
P="nice -n 10 taskset -c $CORES $PY"
$P $T/invsqrt_pp_affine.py "$IN" affine.json > affine.log 2>&1
head -n 1 affine.log
$P $T/invsqrt_pp_affwrite.py "$IN" affine.json invsqrt.kmb.tmp > affwrite.log 2>&1
grep -E "^sites|expected Toffoli|CHECK" affwrite.log
if grep -q "CHECK PASS" affwrite.log; then mv invsqrt.kmb.tmp invsqrt.kmb; sha256sum invsqrt.kmb | cut -c1-16
else echo "rewrite FAILED its check: invsqrt.kmb not written"; rm -f invsqrt.kmb.tmp; exit 1; fi
