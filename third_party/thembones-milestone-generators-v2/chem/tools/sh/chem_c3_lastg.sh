#!/bin/bash
# Builder C3 (fifth holder, Nov 18 08:35Z): ON A BOX. R65's LAST-GROUP REFINEMENT WITHOUT THE GROUP-0 STEP: tools/sh/chem_c3_r65own.sh
# runs the last-group anneal only after a group-0 anneal at 38 succeeded; on t10_g10k1 group 0 never does (0 spaces on k8 / k10),
# so the last group (40 products on k8 / k10 / k12) was never annealed. Here: builder C2's annealer on group 3 of the given dump at
# P3 (chem_ld2_e19grp.py anneal / plan, read-only), then the build with that injection alone (r65_chem_wired221_build through
# tools/chem_c3_vrun.py, owner mode, the dump's own plan / order / seed / tau). Writes OUT/{an_*.json, spaces.json, plans_g3.json,
# s3/angle_k_K.kmx + table_kK.npz, RESULT, cand.txt ('toffoli price dir', for tools/sh/chem_c3_r65merge.sh)}.
#   TB=w1.npz K=10 PRF=prices.npz DUMP=d.json PLAN=6,6:1:3:0,6:1:3:1,7:2:0:2 OAUTO=0,3,1,2 SEED=18 TAU=2e-3 WANT=39,36,36,39 \
#   P3=39 PR3=37 OWNERS="4 8 10 12 14" CORES="0 1 2 3" OUT=... bash chem_c3_lastg.sh
set -u
C=(${CORES:-0 1 2 3}); PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OWNERS
VR="$PY -I $T/chem_c3_vrun.py"
mkdir -p $OUT; cd $OUT || exit 1; rm -f RESULT
log() { echo "$(date -u +%H:%M:%SZ) k$K $*" >> $OUT/lastg.log; }
ARGS="--k $K --plan $PLAN --oauto $OAUTO --avoid 0 --seeds $SEED-$SEED --tries 100 --pad --prices $PRF --tau $TAU --maxbit 3"
log "start: dump $DUMP, plan $PLAN, want $WANT, last group at $P3"
r=1
for c in "${C[@]}"; do
  (cd $T/.. && nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $DUMP --g 3 --P $P3 --price $PR3 --secs ${G3SECS:-300} \
     --multi --rseed $r --save $OUT/an_$r.json > $OUT/an_$r.log 2>&1) &
  r=$((r + 1))
done
wait
$PY - $OUT << 'PYEOF'
import glob, json, sys
d = sys.argv[1]
out = []
for f in sorted(glob.glob(d + "/an_*.json")):
    try:
        x = json.load(open(f))
    except Exception:
        continue
    out += x if isinstance(x, list) else [x]
out = out[::max(1, len(out) // 30)][:30]
json.dump(out, open(d + "/spaces.json", "w"))
PYEOF
n3=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" $OUT/spaces.json 2>/dev/null || echo 0)
log "last-group spaces at $P3: $n3"
[ "${n3:-0}" -gt 0 ] || { echo "none: no last-group space at $P3" > RESULT; exit 0; }
(cd $T/.. && nice -n 10 taskset -c ${C[0]} $PY -I $T/chem_ld2_e19grp.py plan $DUMP $OUT/spaces.json --g 3 --tries 3000 --out $OUT/plans_g3.json \
   > $OUT/plan3.log 2>&1)
[ -s $OUT/plans_g3.json ] || { echo "none: no group-3 plan" > RESULT; exit 0; }
mkdir -p $OUT/s3
(cd $T/.. && nice -n 10 taskset -c ${C[0]} $VR r65_chem_wired221_build $TB $ARGS --want $WANT --maxprice ${MAXP:-0.6} --inj "3=$OUT/plans_g3.json" \
   --emit $OUT/s3 > $OUT/s3/build.log 2>&1)
t=$(grep -h -o -E 'toffoli [0-9]+' $OUT/s3/build.log | tail -n 1 | awk '{print $2}')
p=$(grep -h -o -E 'formula price[a-z ]* [-+][0-9.e]+' $OUT/s3/build.log | tail -n 1 | awk '{print $NF}')
log "built: toffoli ${t:-none}, price ${p:-?}"
echo "built ${t:-none} price ${p:-?}" > RESULT
[ -n "$t" ] && [ -s $OUT/s3/angle_k_$K.kmx ] && echo "$t ${p:-0} $OUT/s3" > $OUT/cand.txt
