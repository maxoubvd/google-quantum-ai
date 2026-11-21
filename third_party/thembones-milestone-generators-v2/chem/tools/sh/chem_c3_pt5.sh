#!/usr/bin/env bash
# Builder C3 (second holder, Nov 17 13:10Z): DESIGN SEAT R62-CHEM-BILIN's PATCH 2 (S5) on a GIVEN stage-A plan: the
# point products read the helpers of the first point (FIVE helpers there, 7 in all), point 13 last:
# 252 = 74 + 149 + 7 + 22 on a 149-product stage A, 127 qubits [M R62 on G8F]. Its two tools, read-only:
# r62_chem_bilin_mask2.py --K 7 --dist 5,1,1 --s5 and r62_chem_bilin_emit2.py. Partition seeds SEED0.. are tried one
# call each until one emits (a seed can fail in phase 1: 'no common output set').
# Appends the watcher's line 'rc=0 <tag>: EMITTED ...' to OUT/bsrun_k<K>.txt.
#   bash chem_c3_pt5.sh PAD PROJ PLANA K OUT TAG [SEED0] [NSEEDS]
set -u
PAD=$1; PROJ=$2; PLANA=$3; K=$4; OUT=$5; TAG=$6; S0=${7:-0}; NS=${8:-4}
PY=${PY:-$HOME/qcc/ws/chem_own/venv/bin/python}
T=${CT:-<HOME>/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=77
mkdir -p "$OUT"
V=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" "$PROJ" "$K")
# third holder (15:50Z): every seed is tried until one emits 22 point products (the S5 count on a 149 stage A); the
# cheapest emitted circuit is kept (on a rewritten table the first emitting seed gave 23 on k5, seeds 2 and 6 gave 22 [M])
best=999
for s in $(seq "$S0" $((S0 + NS - 1))); do
  $PY -I $T/r62_chem_bilin_mask2.py "$PAD" --planA "$PLANA" --K 7 --onlyK --seed "$s" --pseeds 1 --dist 5,1,1 --V "$V" \
      --dc --eqw --s5 --plan "$OUT/plan5_${TAG}_t.json" 2>&1 | grep -E "^SUMMARY|cells emptied|S5|C-4 cond|skipped|rror" \
      | cut -c1-300 > "$OUT/${TAG}_t.log"
  [ -s "$OUT/plan5_${TAG}_t.json" ] || continue
  $PY -I $T/r62_chem_bilin_emit2.py "$PAD" "$OUT/plan5_${TAG}_t.json" --V "$V" --out "$OUT/${TAG}_t.kmx" >> "$OUT/${TAG}_t.log" 2>&1
  n=$(grep -h -o "EMITTED k [0-9]* on the 68-AND cache: [0-9]*" "$OUT/${TAG}_t.log" | tail -n 1 | awk '{print $NF}')
  if [ -n "$n" ] && [ "$n" -lt "$best" ]; then
    best=$n; mv "$OUT/${TAG}_t.kmx" "$OUT/$TAG.kmx"; mv "$OUT/plan5_${TAG}_t.json" "$OUT/plan5_$TAG.json"
    mv "$OUT/${TAG}_t.log" "$OUT/$TAG.log"
  fi
  rm -f "$OUT/plan5_${TAG}_t.json" "$OUT/${TAG}_t.kmx"
  grep -q " 22 point" "$OUT/$TAG.log" 2>/dev/null && break
done
if [ "$best" -lt 999 ]; then
  echo "$(date -u +%H:%M:%SZ) rc=0 $TAG: $(grep EMITTED "$OUT/$TAG.log" | tail -1)" >> "$OUT/bsrun_k$K.txt"
  exit 0
fi
echo "$(date -u +%H:%M:%SZ) END k$K S5: no plan in seeds $S0..$((S0 + NS - 1))" >> "$OUT/bsrun_k$K.txt"
