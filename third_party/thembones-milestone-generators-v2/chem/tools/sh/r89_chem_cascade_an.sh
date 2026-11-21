#!/bin/bash
# r89_chem_cascade_an.sh (seat R89-CHEM-CASCADE): C4's group-3 step (builder C2's annealer chem_ld2_e19grp.py, read-only)
# on one of my builds, then the rebuild with the annealed 42-plans injected (my form in group 4), as chem_c4_loader2.sh does.
#   bash r89_chem_cascade_an.sh K TABLE PAIRDIR ARGSFILE B1DIR "cores" [secs]   -> B1DIR/b2_<r>/angle_k_K.kmx, an.log
K=$1; TB=$2; P=$3; AF=$4; B=$5; CORES=($6); SECS=${7:-300}; R0=${RSEED0:-0}   # RSEED0: first annealer seed
Q=<HOME>/Quantum; [ -d $Q ] || Q=$HOME/Quantum         # this PC, or a box
PY=${PY:-$HOME/qcc/.venv/bin/python}; [ -x $PY ] || PY=$HOME/qcc/.venv_work/bin/python
cd $Q/tools || exit 1
export CHEM_BROT=7 CHEM_HB=77
log() { echo "$(date -u +%H:%M:%SZ) $*" >> $B/an.log; }
log "anneal start on ${CORES[*]} ($SECS s)"
i=$R0
for c in ${CORES[@]}; do
  nice -n 10 taskset -c $c $PY -I chem_ld2_e19grp.py anneal $B/dumpg.json --g 2 --P 42 --secs $SECS --every 20000 \
     --price 40 --rseed $i --save $B/an_$i.json > $B/an_$i.log 2>&1 &
  i=$((i + 1))
done
wait
for r in $(seq $R0 $((i - 1))); do
  [ -s $B/an_$r.json ] || { log "rseed $r: no space"; continue; }
  $PY -I chem_ld2_e19grp.py plan $B/dumpg.json $B/an_$r.json --g 2 --tries 20000 --out $B/pl42_$r.json > $B/pl42_$r.log 2>&1
  log "rseed $r: $(tail -n 1 $B/pl42_$r.log | cut -c1-120)"
  [ -s $B/pl42_$r.json ] || continue
  bash sh/r89_chem_cascade_go.sh $K $TB $P $AF $B/b2_$r ${CORES[0]} $B/pl42_$r.json
  log "b2_$r: $(grep -v '^  attempt' $B/b2_$r/build.log | grep -E 'CLOSED|FAILED' | head -n 1 | cut -c1-160)"
  grep -q "CLOSED 142 products" $B/b2_$r/build.log && { log "SUCCESS 142 at b2_$r"; exit 0; }
done
log "end without 142"
