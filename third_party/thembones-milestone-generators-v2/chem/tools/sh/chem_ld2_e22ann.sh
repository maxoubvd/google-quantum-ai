#!/bin/bash
# chem_ld2_e22ann.sh RUNDIR TABLE K SEED CORE [PLAN] [MODES] [G] [P]: builder C2 (22nd holder). ON A BOX (or WSL with PY
# set). R63-CHEM-WIRED's next step (a) for ONE dead-word owner: the engine run at SEED dumped (tools/chem_ld2_e22dci.py,
# DCI_DUMP), the relation space of generic group G (default 3 = the last 7-group) annealed for P products (default
# the group's count - 1; tools/chem_ld2_e19grp.py anneal, ASECS seconds, two restarts), explicit plans (plan), then the
# SAME seed rebuilt with them injected (DCI_INJECT) and the closure: RUNDIR/a<K>/angle_k_<K>.kmx + table_k<K>.npz.
# Lines in RUNDIR/a<K>.log; RUNDIR/a<K>.done at the end.
cd "$1" || exit 1
TB=$2; K=$3; S=$4; C=$5; PLAN=${6:-6,6,6:1:3:1,7:2:0:2}; MODES=${7:-reg,tie3,-,-}; G=${8:-3}; P=${9:-}
PY=${PY:-~/qcc/.venv_work/bin/python}; T=${CT:-~/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
L=a$K.log; D=a$K; mkdir -p $D
log() { echo "$(date -u +%H:%M:%SZ) $*" >> $L; }
run() { nice -n 10 taskset -c $C $PY -I "$@"; }
ARGS="$TB --k $K --own cur,prev --worder auto --plan $PLAN --modes $MODES --avoid 0,1 --seeds $S-$S --pad"
DCI_DUMP=$D/g.json run $T/chem_ld2_e22dci.py $ARGS > $D/dump.log 2>&1
log "dump: $(grep -h '^seed' $D/dump.log | cut -c1-160)"
[ -s $D/g.json ] || { log "no dump"; touch a$K.done; exit 1; }
M=$(python3 -c "import json, sys; print(json.load(open(sys.argv[1]))['groups'][int(sys.argv[2])]['made'])" $D/g.json $G)
P=${P:-$((M - 1))}
# 15:12Z: one spaces file per anneal (--save overwrites its file), merged before the plans; PAR="c1 c2 .." = that many
# anneals at once on those cores (rseeds RSEEDS in order), else one after another on CORE until NMIN spaces
merge() { python3 -c "
import glob, json, sys
out = []
for f in sorted(glob.glob(sys.argv[1] + '/sp_*.json')):
    try:
        d = json.load(open(f))
    except Exception:
        continue
    out += d if isinstance(d, list) else [d]
json.dump(out, open(sys.argv[1] + '/sp.json', 'w'))
print(len(out))" $D; }
RS=(${RSEEDS:-1 2 3 4}); PC=(${PAR:-})
if [ ${#PC[@]} -gt 0 ]; then
  j=0
  for r in "${RS[@]}"; do
    c=${PC[$((j % ${#PC[@]}))]}; j=$((j + 1))
    ( nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $D/g.json --g $G --P $P --secs ${ASECS:-240} --save $D/sp_$r.json --multi --rseed $r > $D/ann$r.log 2>&1 ) &
    [ $((j % ${#PC[@]})) = 0 ] && wait
  done
  wait
  NS=$(merge); log "anneals ${RS[*]} on cores ${PC[*]}: group $G made $M, target $P: ${NS:-0} spaces"
else
  for r in "${RS[@]}"; do
    run $T/chem_ld2_e19grp.py anneal $D/g.json --g $G --P $P --secs ${ASECS:-240} --save $D/sp_$r.json --multi --rseed $r > $D/ann$r.log 2>&1
    NS=$(merge); log "anneal rseed $r: group $G made $M, target $P: ${NS:-0} spaces"
    [ "${NS:-0}" -ge ${NMIN:-6} ] && break
  done
fi
[ "${NS:-0}" -gt 0 ] || { log "no space at $P"; touch a$K.done; exit 0; }
run $T/chem_ld2_e19grp.py plan $D/g.json $D/sp.json --g $G --out $D/plans.json --tries 3000 > $D/plan.log 2>&1
log "plans: $(tail -n 1 $D/plan.log)"
DCI_INJECT=$G=$D/plans.json run $T/chem_ld2_e22dci.py $ARGS --emit $D > $D/build.log 2>&1
log "build: $(grep -h -E '^seed|emitted|differs' $D/build.log | tr '\n' ' ' | cut -c1-300)"
touch a$K.done
