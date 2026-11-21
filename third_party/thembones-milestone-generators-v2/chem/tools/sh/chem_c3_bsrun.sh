#!/usr/bin/env bash
# Builder C3: the kept-product loader with a SEARCHED STAGE-A BASIS (tools/chem_c3_bsearch.py: blocks of rank 74 =
# codim 16) for ONE loader: NS basis searches, stage A (C2's engine via r58_chem_kept_c2a.py) under a few engine
# seeds / efforts, then the point stage (chem_c3_c4chk.py, 8 helpers, C-4 conditions) and the 68-AND emitter on the
# smallest stage A. Lines in OUT/bsrun_k<K>.txt.  Run pinned: taskset -c <core> nice -n 10 bash chem_c3_bsrun.sh ...
#   bash chem_c3_bsrun.sh IN K OUT [NS] [D] [TRIES]
set -u
IN=$1; K=$2; OUT=$3; NS=${4:-3}; D=${5:-13,14,15}; TRIES=${6:-300000}
PY=${PY:-$HOME/qcc/ws/chem_own/venv/bin/python}
T=${CT:-<HOME>/Quantum/tools}      # on a box: CT=~/Quantum/tools PY=~/qcc/.venv_work/bin/python
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77}
B=$OUT/in_k$K; mkdir -p "$B"
ln -sf "$IN/pad.npz" "$B/pad.npz"; ln -sf "$IN/proj.json" "$B/proj.json"
R=$OUT/bsrun_k$K.txt
best=999; bestpf=""; bestA=""
S0=${SEED0:-1}
for s in $(seq "$S0" $((S0 + NS - 1))); do
  PF=k${K}bs${s}_D${D//,/-}.json
  $PY -I $T/chem_c3_bsearch.py "$B/pad.npz" --k "$K" --D "$D" --tries "$TRIES" --seed "$s" --out "$B/$PF" 2>&1 | tail -2 >> "$R"
  [ -s "$B/$PF" ] || continue
  IFS='|' read -r -a CMB <<< "${COMBOS:-0 1|1 1|2 1|0 2|1 2}"
  for se in "${CMB[@]}"; do
    set -- $se
    PA=$OUT/planA_${PF%.json}_s$1e$2.json
    $PY -I $T/r58_chem_kept_c2a.py "$B/pad.npz" --k "$K" --pfile "$B/$PF" --rep 4 --D "$D" --seed "$1" --effort "$2" \
        --out "$PA" > "$PA.log" 2>&1
    n=$(sed -n 's/.*four blocks) \([0-9]*\) products.*/\1/p' "$PA.log")
    echo "$(date -u +%H:%M:%SZ) k$K $PF seed $1 effort $2: stage A ${n:-fail} $(grep -o 'codim \[[0-9, ]*\]' "$PA.log")" >> "$R"
    if [ -n "$n" ] && [ "$n" -lt "$best" ]; then best=$n; bestpf=$PF; bestA=$PA; fi
  done
done
[ -n "$bestA" ] || { echo "$(date -u +%H:%M:%SZ) END k$K: no stage A" >> "$R"; exit 4; }
V=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" "$IN/proj.json" "$K")
TAG=$(basename "${bestA%.json}" | sed 's/^planA_//')
if [ "${RELIEF:-1}" = "1" ]; then
  # since 11:10Z: the re-pricer's DON'T-CARE RELIEF (R58-CHEM-KEPT-reprice.md section 5): block 76 free at D,
  # equal helper classes, split 6,1,1, many partition seeds; its emitter checks the rows the verifier runs
  # since 12:55Z (second holder): design seat R62-CHEM-BILIN's point stage (class-by-class helper classes, capacity at
  # 75 equations a cell): 30 + 8 + 8 cells on a 149-product stage A = 253 [M R62 on G8F's six]; --cls 200, 12 seeds
  TAG=${TAG}_bl
  $PY -I $T/r62_chem_bilin_mask.py "$B/pad.npz" --planA "$bestA" --K 8 --onlyK --seed 0 --pseeds "${BLSEEDS:-12}" \
      --dist 6,1,1 --V "$V" --dc --eqw --cls "${BLCLS:-200}" --plan "$OUT/plan68_$TAG.json" 2>&1 \
      | grep -E "^SUMMARY|cells emptied|C-4 cond|skipped|^R62" | cut -c1-300 > "$OUT/$TAG.log"
  $PY -I $T/r58_chem_kept_rp_emit68.py "$B/pad.npz" "$OUT/plan68_$TAG.json" --V "$V" --out "$OUT/$TAG.kmx" >> "$OUT/$TAG.log" 2>&1
  echo "$(date -u +%H:%M:%SZ) rc=0 $TAG: $(tail -1 "$OUT/$TAG.log")" >> "$R"
  # 13:10Z: R62-CHEM-BILIN's PATCH 2 (S5: 7 helpers, the point products read the first point's helpers): 252 on a
  # 149-product stage A [M R62]; up to four partition seeds; its line joins the same file (the watcher takes the min)
  [ "${S5:-1}" = "1" ] && CT=$T PY=$PY bash $T/sh/chem_c3_pt5.sh "$B/pad.npz" "$IN/proj.json" "$bestA" "$K" "$OUT" \
      "$(basename "${bestA%.json}" | sed 's/^planA_//')_s5" 0 4
  echo "$(date -u +%H:%M:%SZ) DONE k$K: stage A $best" >> "$R"
  exit 0
else
  $PY -I $T/chem_c3_c4chk.py "$B/pad.npz" --planA "$bestA" --K 8 --onlyK --pseeds "${PSEEDS:-4}" --V "$V" \
      --plan "$OUT/plan68_$TAG.json" 2>&1 | grep -v "^   point" | cut -c1-300 > "$OUT/$TAG.log"
  $PY -I $T/chem_c3_emit68.py "$B/pad.npz" "$OUT/plan68_$TAG.json" --V "$V" --out "$OUT/$TAG.kmx" >> "$OUT/$TAG.log" 2>&1
fi
echo "$(date -u +%H:%M:%SZ) rc=0 $TAG: $(tail -1 "$OUT/$TAG.log")" >> "$R"
echo "$(date -u +%H:%M:%SZ) DONE k$K: stage A $best" >> "$R"
