#!/bin/bash
# invsqrt_gen: run the lookup workstream's program generator (tools/lut_synth_gen.py, its own venv, both
# read-only for us) on the block specs of the last local build and write the programs into OUR folder
# ~/qcc/ws/invsqrt_gen/lut/ext/. Nothing is written into the lookup workstream's folders.
#   invsqrt_gen_ext.sh <budget> <block> [<block> ...]          method A (any block)
#   METHOD=designA invsqrt_gen_ext.sh <budget> LUT1            consumed-algebra design (address inputs only)
# usage from PowerShell:
#   wsl -d Ubuntu -e bash -lc "tr -d '\r' < <HOME>/Quantum/tools/sh/invsqrt_gen_ext.sh | bash -s -- 14 LUT3 LUT1b"
set -u
T=<HOME>/Quantum/tools
W=~/qcc/ws/invsqrt_gen
PY=~/qcc/ws/lut_synth/venv/bin/python
export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 NUMBA_NUM_THREADS=3
budget="$1"; shift
# one folder per (method, budget, seed), so that a later run never overwrites a better program
O=lut/ext/${METHOD:-A}_b${budget}_s${SEED:-0}
mkdir -p $W/$O
cd $W || exit 1
for b in "$@"; do
  nice -n 5 $PY $T/lut_synth_gen.py lut/$b.spec.json --budget $budget --seed ${SEED:-0} --method ${METHOD:-A} ${EXTRA:-} --out $O < /dev/null 2>&1 | grep -E "RESULT|Error|error|Traceback" | cut -c1-220
done
