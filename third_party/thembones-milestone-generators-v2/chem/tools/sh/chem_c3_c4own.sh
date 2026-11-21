#!/bin/bash
# Builder C3 (sixth holder, Nov 18 09:55Z): ON A BOX. Builder C4's FIRST-GROUP SCREEN on one R65 owner + my LAST-GROUP FIX:
#  1. R65's placement of the dumped seed (DD/dump.log: 'outputs placed: w:b,...': 25 wires in R65's group order);
#  2. tools/chem_c4_ownscan.py (C4, read-only): every 6 of R65's group-1 + group-4 wires as a pool-free first group at 38,
#     non-silent on coordinate 14 (one part a core);
#  3. for the first NTRY sets found: ORDER = the set + R65's groups 2 and 3 + the rest of the pool; R65's build with
#     --oorder (dump), group 0 planned at 38 (tools/chem_c4_g0plans.py), the build with --inj 0 (C4's route: k10 of
#     t10_g10k1 at 224 [M C4]) and its dump; then the last group annealed one lower (builder C2's chem_ld2_e19grp.py) and
#     the build with both injections (aimed at 38|36|36|39 = 223 on such a loader).
#  Writes OUT/c<i>/s0 and OUT/c<i>/s3 (angle_k_K.kmx, table_kK.npz, build.log) and OUT/RESULT.
#   TB=w.npz K=10 PRF=prices.npz DD=<refine dir> PLAN=6,6:1:3:0,6:1:3:1,7:2:0:2 SEED=13 TAU=2e-3 CORES="0 1 2 3" OUT=... \
#     [NTRY=2] [MAXP=0.6] bash chem_c3_c4own.sh
set -u
C=(${CORES:-0 1 2 3}); N=${#C[@]}; PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
VR="$PY -I $T/chem_c3_vrun.py"
mkdir -p $OUT; rm -f $OUT/RESULT
log() { echo "$(date -u +%H:%M:%SZ) k$K $*" >> $OUT/c4own.log; }
run1() { local c=$1 lf=$2; shift 2; (cd $T/.. && nice -n 10 taskset -c $c "$@" > $lf 2>&1); }
PLACE=$(sed -n 's/.*outputs placed: \([0-9:,]*\).*/\1/p' $DD/dump.log | head -n 1)
[ -n "$PLACE" ] || { log "no placement in $DD/dump.log"; echo "none: no placement" > $OUT/RESULT; exit 0; }
log "start: placement $PLACE"
# 2. the screen, one part a core
for i in $(seq 0 $((N - 1))); do
  run1 ${C[$i]} $OUT/scan_$i.log $PY -I $T/chem_c4_ownscan.py $TB --k $K --place "$PLACE" --t 14 --part $i/$N --out $OUT/scan_$i.json &
done
wait
SETS=$($PY - $OUT "$PLACE" ${NTRY:-2} << 'PYEOF'
import glob, json, sys
d, place, ntry = sys.argv[1], sys.argv[2].split(","), int(sys.argv[3])
found = []
for f in sorted(glob.glob(d + "/scan_*.json")):
    try:
        found += json.load(open(f))
    except Exception:
        pass
g = [place[0:6], place[6:12], place[12:18], place[18:25]]
# 10:05Z: C4's rule (measured on k10: taking two of R65's last-group wires cost that group a kill, 40 -> 43; mine too:
# 38 | 36 | 36 | 43 on the set with the most fan-outs [M]): sets that keep the most of R65's first-group wires first
found.sort(key=lambda x: (-sum(1 for w in x["g1"] if w in g[0]), -x["fans"]))
found = [x for x in found if sum(1 for w in x["g1"] if w in g[0]) >= 5] or found
for x in found[:ntry]:
    g1 = x["g1"]
    rest = [w for w in g[0] + g[3] if w not in g1]
    print(",".join(g1 + g[1] + g[2] + rest))
PYEOF
)
log "screen: $(cat $OUT/scan_*.json 2>/dev/null | grep -o '"g1"' | wc -l) first groups at 38 found; trying $(echo "$SETS" | grep -c ,)"
[ -n "$SETS" ] || { echo "none: no first group at 38" > $OUT/RESULT; exit 0; }
i=0
for ORDER in $SETS; do
  W=$OUT/c$i; mkdir -p $W/s0 $W/s3
  A="--k $K --plan $PLAN --oorder $ORDER --avoid 0 --seeds $SEED-$SEED --tries 400 --pad --prices $PRF --tau $TAU --maxbit 3 --maxprice ${MAXP:-0.6}"
  run1 ${C[0]} $W/dump.log $VR r65_chem_wired221_build $TB $A --dump $W/d.json
  [ -s $W/d.json ] || { log "c$i: no dump"; i=$((i + 1)); continue; }
  run1 ${C[0]} $W/g0.log $PY -I $T/chem_c4_g0plans.py $W/d.json --g 0 --t 14 --out $W/p0.json
  [ -s $W/p0.json ] || { log "c$i: no group-0 plans ($(tail -n 1 $W/g0.log | cut -c1-100))"; i=$((i + 1)); continue; }
  run1 ${C[0]} $W/s0/build.log $VR r65_chem_wired221_build $TB $A --inj "0=$W/p0.json" --dump $W/d2.json --emit $W/s0
  t0=$(grep -h -o -E 'toffoli [0-9]+' $W/s0/build.log | tail -n 1 | awk '{print $2}')
  gw=$(sed -n 's/^DUMP.*groups \[\([0-9, ]*\)\].*/\1/p' $W/s0/build.log | head -n 1 | tr -d ' ')
  log "c$i: first group at 38 injected: toffoli ${t0:-none}, groups ${gw:-?}"
  [ -s $W/d2.json ] && [ -n "$gw" ] || { i=$((i + 1)); continue; }
  gl=${gw##*,}; wl=$(echo $gw | awk -F, -v OFS=, '{$NF = $NF - 1; print}')
  r=1
  for c in "${C[@]}"; do
    (cd $T/.. && nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $W/d2.json --g 3 --P $((gl - 1)) --price $((gl - 3)) \
       --secs ${G3SECS:-300} --multi --rseed $r --save $W/an_$r.json > $W/an_$r.log 2>&1) &
    r=$((r + 1))
  done
  wait
  $PY - $W << 'PYEOF'
import glob, json, sys
d = sys.argv[1]
out = []
for f in sorted(glob.glob(d + "/an_*.json")):
    try:
        x = json.load(open(f))
    except Exception:
        continue
    out += x if isinstance(x, list) else [x]
json.dump(out[::max(1, len(out) // 30)][:30], open(d + "/spaces3.json", "w"))
PYEOF
  n3=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" $W/spaces3.json 2>/dev/null || echo 0)
  log "c$i: last-group spaces at $((gl - 1)): $n3"
  if [ "${n3:-0}" -gt 0 ]; then
    run1 ${C[0]} $W/plan3.log $PY -I $T/chem_ld2_e19grp.py plan $W/d2.json $W/spaces3.json --g 3 --tries 3000 --out $W/p3.json
    [ -s $W/p3.json ] && run1 ${C[0]} $W/s3/build.log $VR r65_chem_wired221_build $TB $A --want $wl --inj "0=$W/p0.json;3=$W/p3.json" \
        --emit $W/s3
    log "c$i: both injections: toffoli $(grep -h -o -E 'toffoli [0-9]+' $W/s3/build.log 2>/dev/null | tail -n 1 | awk '{print $2}')"
  fi
  i=$((i + 1))
done
B=$(for f in $OUT/c*/s*/build.log; do t=$(grep -h -o -E 'toffoli [0-9]+' $f | tail -n 1 | awk '{print $2}'); [ -n "$t" ] && echo "$t $(dirname $f)"; done | sort -n | head -n 1)
echo "best ${B:-none}" > $OUT/RESULT; log "RESULT best ${B:-none}"
