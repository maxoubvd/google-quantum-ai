#!/bin/bash
# r89_chem_cascade_go.sh (seat R89-CHEM-CASCADE): one k-loader build with the forced group-4 relation space
# (r89_chem_cascade_build.py = C4's chem_c4_build.py + my hook) + C4's cubic block, groups 0-1 injected from a C4 pair,
# optional group-2 plans (C4's annealed 42s). Output: OUT/angle_k_K.kmx, OUT/table_dead.npz, OUT/build.log, OUT/dumpg.json
#   bash r89_chem_cascade_go.sh K TABLE PAIRDIR ARGSFILE OUT CORE [PL42]      (env R89_G, R89_T2, R89_TRIALS pass through)
K=$1; TB=$2; P=$3; AF=$4; OUT=$5; CORE=$6; PL=${7:-}
mkdir -p $OUT
Q=<HOME>/Quantum; [ -d $Q ] || Q=$HOME/Quantum         # this PC, or a box
PY=${PY:-$HOME/qcc/.venv/bin/python}; [ -x $PY ] || PY=$HOME/qcc/.venv_work/bin/python
cd $Q/tools || exit 1
export CHEM_BROT=7 CHEM_HB=77 C4_DROP=1 C4_DEAD=6 C4_CUBIC=3 C4_CUBIC_TRIES=${C4_CUBIC_TRIES:-1000} C4_DEAD_OUT=EMIT
INJ2=""; [ -n "$PL" ] && INJ2="--inj 2=$PL"
nice -n 10 taskset -c $CORE $PY -I r89_chem_cascade_build.py $TB --k $K --wirebasis \
  --plan 6,6,7:1:0:0,4:1:0:1 $(cat $AF) --inj 0=$P/g0_best.json --inj 1=$P/g1_best.json $INJ2 --modes=-,-,-,- \
  --tries 3000 --emit $OUT --dumpg $OUT/dumpg.json > $OUT/build.log 2>&1
echo "$(date -u +%H:%M:%SZ) done: $(grep -v '^  attempt' $OUT/build.log | grep -E 'CLOSED|FAILED' | cut -c1-200)" >> $OUT/build.log
