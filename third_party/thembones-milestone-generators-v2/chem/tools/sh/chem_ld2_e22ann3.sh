#!/bin/bash
# chem_ld2_e22ann3.sh RUNDIR TABLE K SEED "cores" [SECS]: builder C2 (22nd holder). ON A BOX. A dead-word OWNER WITHOUT
# TIES, three groups re-planned from annealed relation spaces (tools/chem_ld2_e22dci.py hooks over my pinned copy of
# R63's tool; tools/chem_ld2_e19grp.py anneal / plan), plan 6,6,6:1:3:1,7:2:0:2, modes -,-,-,-, the closure on:
#  1. dump at SEED; groups 1 and 2 (index 0, 1) annealed for one product less (39 -> 38 each), half the cores each;
#  2. rebuilt with both injected and dumped again (the later groups' pools are the earlier products);
#  3. the last group (index 3) annealed for one product less (40 -> 39); 4. rebuilt with the three injected:
#     RUNDIR/b<K>x/angle_k_<K>.kmx + table_k<K>.npz (38 / 38 / 36 / 39 = 151 = 225 when all three land).
# Lines in RUNDIR/b<K>x.log; RUNDIR/b<K>x.done at the end. Every stage that finds nothing keeps the engine's plan.
cd "$1" || exit 1
TB=$2; K=$3; S=$4; CS=($5); SECS=${6:-600}
PY=${PY:-~/qcc/.venv_work/bin/python}; T=${CT:-~/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
D=${DX:-b${K}x}; L=$D.log; mkdir -p $D; rm -f $D.done   # DX: another folder (plans_g0 / plans_g1 put there are kept)
log() { echo "$(date -u +%H:%M:%SZ) $*" >> $L; }
ARGS="$TB --k $K --own cur,prev --worder auto --plan 6,6,6:1:3:1,7:2:0:2 --modes=-,-,-,- --avoid 0,1 --seeds $S-$S --pad"
ann() {   # ann DUMP G P "cores" -> D/plans_g<G>.json (or nothing)
  local dump=$1 g=$2 p=$3; local cs=($4); local j=0
  for c in "${cs[@]}"; do
    ( nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $dump --g $g --P $p --secs $SECS \
        --save $D/sp_g${g}_$j.json --multi --rseed $((17 + j)) > $D/ann_g${g}_$j.log 2>&1 ) &
    j=$((j + 1))
  done
  wait
  local n=$(python3 -c "
import glob, json
out = []
for f in sorted(glob.glob('$D/sp_g${g}_*.json')):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    out += d if isinstance(d, list) else [d]
json.dump(out, open('$D/sp_g$g.json', 'w'))
print(len(out))")
  log "group $g for $p: ${n:-0} spaces ($(for f in $D/ann_g${g}_*.log; do tail -n 1 $f | grep -o -E 'best planes [0-9]+'; done | sort -t' ' -k3 -n | tail -n 1))"
  if [ "${n:-0}" -gt 0 ]; then
    nice -n 10 taskset -c ${cs[0]} $PY -I $T/chem_ld2_e19grp.py plan $dump $D/sp_g$g.json --g $g --out $D/plans_g$g.json \
      --tries 3000 > $D/plan_g$g.log 2>&1
    log "group $g plans: $(tail -n 1 $D/plan_g$g.log)"
  fi
}
made() { python3 -c "import json, sys; print(json.load(open(sys.argv[1]))['groups'][int(sys.argv[2])]['made'])" $1 $2; }
# 1. dump; groups 1 and 2 at once (half the cores each)
DCI_DUMP=$D/g1.json nice -n 10 taskset -c ${CS[0]} $PY -I $T/chem_ld2_e22dci.py $ARGS > $D/dump1.log 2>&1
log "dump 1: $(grep -h '^seed' $D/dump1.log | cut -c1-150)"
[ -s $D/g1.json ] || { log "no dump"; touch $D.done; exit 1; }
H=$(( ${#CS[@]} / 2 ))
# 16:36Z: a group whose plans are already in the folder (e.g. from extra anneals) is not annealed again
[ -s $D/plans_g0.json ] && log "group 0: plans given" || ( ann $D/g1.json 0 $(( $(made $D/g1.json 0) - 1 )) "${CS[*]:0:$H}" ) &
[ -s $D/plans_g1.json ] && log "group 1: plans given" || ( ann $D/g1.json 1 $(( $(made $D/g1.json 1) - 1 )) "${CS[*]:$H}" ) &
wait
INJ=""; for g in 0 1; do [ -s $D/plans_g$g.json ] && INJ="$INJ;$g=$D/plans_g$g.json"; done; INJ=${INJ#;}
# 2. rebuilt with groups 1-2 injected, dumped again
DCI_INJECT="$INJ" DCI_DUMP=$D/g2.json nice -n 10 taskset -c ${CS[0]} $PY -I $T/chem_ld2_e22dci.py $ARGS > $D/dump2.log 2>&1
log "dump 2 (injected: ${INJ:-none}): $(grep -h '^seed' $D/dump2.log | cut -c1-150)"
[ -s $D/g2.json ] || { log "no dump 2"; touch $D.done; exit 1; }
# 3. the last group
ann $D/g2.json 3 $(( $(made $D/g2.json 3) - 1 )) "${CS[*]}"
[ -s $D/plans_g3.json ] && INJ="${INJ:+$INJ;}3=$D/plans_g3.json"
# 4. the final build
DCI_INJECT="$INJ" nice -n 10 taskset -c ${CS[0]} $PY -I $T/chem_ld2_e22dci.py $ARGS --emit $D > $D/build.log 2>&1
log "build (injected: ${INJ:-none}): $(grep -h -E '^seed|emitted|differs' $D/build.log | tr '\n' ' ' | cut -c1-300)"
touch $D.done
