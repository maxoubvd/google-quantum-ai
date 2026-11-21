#!/usr/bin/env bash
# Builder C3 (third holder, Nov 17 16:40Z): ON A BOX, design seat R66-CHEM-PLAN244's S5T (S5 + the two-point PAIR
# helper on the 8th wire: tools/r66_chem_plan244_mask3.py --tp + tools/r66_chem_plan244_emit3.py, read-only) for ONE
# kept loader on a GIVEN stage A and padding, with the seat's settings (F 10, c12 31, c3 20), NP processes of the same
# seed and different lseeds sharing the lucky classes (as the seat's r66_chem_plan244_job.sh), into MY folder.
# Lines in OUT/s5t_k<K>.txt ('EMITTED ...'); circuits OUT/kmx_<TAG>_l<LSEED>.kmx.
#   bash chem_c3_s5t.sh PAD PROJ PLANA K OUT TAG "CORES" [SEED] [F] [C12] [C3]
set -u
PAD=$1; PROJ=$2; PLANA=$3; K=$4; OUT=$5; TAG=$6; CO=($7); SEED=${8:-0}; FN=${9:-10}; C12=${10:-31}; C3=${11:-20}
PY=${PY:-$HOME/qcc/.venv_work/bin/python}; T=${CT:-$HOME/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p "$OUT/share_$TAG"
V=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" "$PROJ" "$K")
R=$OUT/s5t_k$K.txt
echo "$(date -u +%H:%M:%SZ) k$K V $V planA $(basename $PLANA) seed $SEED: ${#CO[@]} processes" >> "$R"
l=${LSEED0:-0}
for c in "${CO[@]}"; do
  ( O=$OUT/plan_${TAG}_l$l.json
    # 19:38Z (fourth holder): a time limit (S5T_SECS, default 900 s since 20:37Z): a stalled class search ends and the caller keeps S5
    # Nov 18 12:15Z (seventh holder): 1500 s (ch2 on t10_g10k1: k1 / k7 stopped at 900 s with 2 of 4 classes; the
    # successes found their 4th class at 313-645 s [M]: a miss is +1 Toffoli a pass)
    nice -n 10 taskset -c $c timeout ${S5T_SECS:-1500} $PY -I $T/r66_chem_plan244_mask3.py "$PAD" --planA "$PLANA" --K 8 --onlyK --seed $SEED \
        --pseeds 1 --dist 5,1,1 --V $V --dc --eqw --s5 --tp --F $FN --c12 $C12 --c3 $C3 --s5bases 4000000 \
        --lshare "$OUT/share_$TAG" --lseed $l --plan $O > "$OUT/mask_${TAG}_l$l.log" 2>&1
    if [ -s $O ]; then
      nice -n 10 taskset -c $c $PY -I $T/r66_chem_plan244_emit3.py "$PAD" $O --V $V --out "$OUT/kmx_${TAG}_l$l.kmx" \
          > "$OUT/emit_${TAG}_l$l.log" 2>&1
      echo "$(date -u +%H:%M:%SZ) l$l: $(grep -h -E 'EMITTED|rror' "$OUT/emit_${TAG}_l$l.log" | tail -n 1 | cut -c1-200)" >> "$R"
    else
      echo "$(date -u +%H:%M:%SZ) l$l: no plan ($(grep -h -E '^SUMMARY|NO PLAN|rror' "$OUT/mask_${TAG}_l$l.log" | tail -n 1 | cut -c1-150))" >> "$R"
    fi ) &
  l=$((l + 1))
done
wait
echo "$(date -u +%H:%M:%SZ) DONE k$K" >> "$R"
