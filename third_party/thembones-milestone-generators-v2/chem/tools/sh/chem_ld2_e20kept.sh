#!/bin/bash
# chem_ld2_e20kept.sh <staircase in-folder> <out dir> <core> K [K ...]: builder C2 (twentieth holder, Nov 17).
# R58-CHEM-KEPT's patch (tools/r58_chem_kept_*.py, read-only) on 30-output loaders: my engine on four blocks + the
# point stage with class-mask helpers. Per loader and (rep, D) choice: FIRST the 68-AND cache with 8 helpers
# (c2a -> c4chk --V <the loader's block set from proj.json> -> emit68: 74 + 182 = 256 [V on GCL k1-k3, M on G8F k2 by
# C3]); if that does not emit, the 69-AND cache with 7 helpers (c2a -> mask -> emit: 75 + 183 = 258 [V on G8F k1, k5]).
# Keeps the cheapest circuit as <out>/k<K>/angle_k_<K>.kmx; one line per try in <out>/RESULT.txt.
# FIRST=1 (default): stop a loader at its first 256 (or, after every choice, keep its best).
IN=$1; OUT=$2; C=$3; shift 3
T=${T:-<HOME>/Quantum/tools}
PY=${PY:-$HOME/qcc/ws/chem_own/venv/bin/python}
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $OUT
# twenty-first holder (Nov 17 11:05Z; LEAD_NOTES 10:46Z; research/panel/R58-CHEM-KEPT-reprice.md section 5): THE
# DON'T-CARE RELIEF: the 68-AND point stage is R58's re-pricer's r58_chem_kept_rp_mask.py --dc --eqw (block 76's
# never-run rows free at the points of D, equal helper classes; 100 partitions) + r58_chem_kept_rp_emit68.py: 180
# products = 254 a loader [V on G8F and GCL by the re-pricer]. NREAL = the rows the verifier runs (N + R B of the
# table, 1,219 on every (9, 127, C) table); D stays inside low points 3-15. KTARGET (default 254) ends a loader.
NREAL=${NREAL:-$($PY -c "import numpy as np, sys; z = np.load(sys.argv[1]); print(int(z['U'].shape[-1] + z['W'].shape[0] * z['W'].shape[1]) if 'W' in z.files and 'U' in z.files else 1219)" $IN/pad.npz 2>/dev/null)}
NREAL=${NREAL:-1219}
KT=${KTARGET:-254}
# twenty-first holder (11:12Z): by default the relief kept step runs ON A BOX, one core a loader, all loaders at once
# (tools/sh/chem_ld2_e21keptbox.sh; cores claimed from KPOOL by lock folders: box 5's idle cores 0-9, 52-57 = the
# LEAD_NOTES 10:33Z allowance of up to 16 idle cores); KBOX=0 = the old local loop below
if [ "${KBOX:-5}" != 0 ]; then
  exec bash "$(dirname "$0")/chem_ld2_e21keptbox2.sh" "$IN" "$OUT" "${KBOX:-5}" "${KPOOL:-0 1 2 3 4 5 6 7 8 9 52 53 54 55 56 57}" "$@"
fi
for K in "$@"; do
  D=$OUT/k$K; mkdir -p $D
  V=$($PY -c "import json, sys; print(json.load(open(sys.argv[1])).get('V', {}).get(sys.argv[2], ''))" $IN/proj.json $K 2>/dev/null)
  best=""
  for rd in "0:13,14,15" "2:13,14,15" "1:13,14,15" "0:4,10,14" "2:5,10,14" "1:4,10,14" "0:5,10,14" "1:5,10,14" "2:4,10,14"; do
    R=${rd%%:*}; DD=${rd#*:}; tag=r${R}_$(echo $DD | tr ',' '_')
    nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_c2a.py $IN/pad.npz --k $K --pfile $IN/k$K.json --rep $R --D $DD \
      --out $D/planA_$tag.json > $D/c2a_$tag.log 2>&1 || { echo "k$K $tag: c2a stopped: $(tail -n 1 $D/c2a_$tag.log | cut -c1-110)" >> $OUT/RESULT.txt; continue; }
    tof=""
    if [ -n "$V" ]; then          # 68-AND cache, 8 helpers; helper distributions given one by one (R58's tool stops
      for dist in ${DISTS:-6,1,1 5,2,1 4,2,2 4,3,1 3,3,2}; do   # on a distribution with no admissible block)
        dt=$(echo $dist | tr ',' '_')
        if nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_rp_mask.py $IN/pad.npz --planA $D/planA_$tag.json --K 8 --onlyK \
             --seed 0 --pseeds ${PSEEDS:-100} --dist $dist --V $V --dc --eqw --nreal $NREAL \
             --plan $D/plan68_${tag}_$dt.json > $D/c4chk_${tag}_$dt.log 2>&1 && \
           nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_rp_emit68.py $IN/pad.npz $D/plan68_${tag}_$dt.json --V $V \
             --nreal $NREAL --out $D/k${K}_${tag}_${dt}_68.kmx > $D/emit68_${tag}_$dt.log 2>&1; then
          tof=$(grep -h -o -E "[0-9]+ Toffolis" $D/emit68_${tag}_$dt.log | tail -n 1 | grep -o -E "^[0-9]+")
          echo "k$K $tag dist $dist: 68-AND emitted (relief), toffoli ${tof:-?}" >> $OUT/RESULT.txt
          if [ -n "$tof" ] && { [ -z "$best" ] || [ "$tof" -lt "$best" ]; }; then best=$tof; cp $D/k${K}_${tag}_${dt}_68.kmx $D/angle_k_$K.kmx; fi
          [ -n "$best" ] && [ "$best" -le $KT ] && break
        else
          echo "k$K $tag dist $dist: 68-AND not emitted: $(cat $D/emit68_${tag}_$dt.log $D/c4chk_${tag}_$dt.log 2>/dev/null | grep -v '^$' | tail -n 1 | cut -c1-90)" >> $OUT/RESULT.txt
        fi
      done
    fi
    if [ -z "$best" ] || [ "$best" -gt 258 ]; then     # 69-AND cache, 7 helpers
      nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_mask.py $IN/pad.npz --planA $D/planA_$tag.json --K 7 --onlyK \
        --pseeds 2 --plan $D/plan_$tag.json > $D/mask_$tag.log 2>&1 && \
      nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_emit.py $IN/pad.npz $D/plan_$tag.json --out $D/k${K}_$tag.kmx \
        > $D/emit_$tag.log 2>&1 && {
        t9=$(grep -h -o -E "[0-9]+ Toffolis" $D/emit_$tag.log | tail -n 1 | grep -o -E "^[0-9]+")
        echo "k$K $tag: 69-AND emitted, toffoli ${t9:-?}" >> $OUT/RESULT.txt
        if [ -n "$t9" ] && { [ -z "$best" ] || [ "$t9" -lt "$best" ]; }; then best=$t9; cp $D/k${K}_$tag.kmx $D/angle_k_$K.kmx; fi
      }
    fi
    [ -n "$best" ] && [ "$best" -le $KT ] && [ "${FIRST:-1}" = 1 ] && break
  done
  echo "k$K: BEST ${best:-none}" >> $OUT/RESULT.txt
done
