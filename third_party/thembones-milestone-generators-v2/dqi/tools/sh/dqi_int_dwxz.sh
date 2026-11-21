#!/bin/bash
# dqi_int_dwxz.sh -- B5 (19th holder, Nov 17; runs ON a box): P1's dwx (tools/sh/dqi_pp_dwx.sh, copied, not edited) with
# R55's EXACT block-aware zero intervals (tools/r55_dqi_gpeak_zero.py, read only: "a drop-in for every census / pass",
# 0 soundness violations on 4,096 samples) given to P1's links pass as its extra zero file (dqi_pp_exact.sh's 5th
# argument; dqi_pp_ret.py merges overlapping intervals per wire). On t54ch's parked file R55's analysis frees 3,589.4M
# wire-ops against P1's 3,268.6M (+9.8%): un-synthesis, un-extraction, scan and unconvert's lower levels.
#   bash dqi_int_dwxz.sh <raw dqi.kmb> <fresh out dir> [cores=4-11]
set -u
IN=$(readlink -f "$1"); OUT=$2; C=${3:-4-11}; T=~/Quantum/tools; PY=~/qcc/.venv_work/bin/python
[ -e "$OUT" ] && { echo "$OUT exists"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
bash $T/sh/dqi_pp_deadw.sh $IN dw1 $C 200 || exit 1
nice -n 10 taskset -c $C $PY -B $T/r55_dqi_gpeak_zero.py $(readlink -f dw1/dqi.kmb) dw1_r55zero.npz > dw1_r55zero.log 2>&1 || exit 1
grep -a -E "exact idle-zero|wire-ops free" dw1_r55zero.log
tr -d '\r' < $T/sh/dqi_pp_exact.sh > exact_frozen.sh
bash exact_frozen.sh $(readlink -f dw1/dqi.kmb) links $C 200 $(readlink -f dw1_r55zero.npz) | tail -n 5
bash $T/sh/dqi_pp_deadw.sh $(readlink -f links/dqi.kmb) final $C 200
echo "dwxz done $(date -u +%H:%M)Z"
