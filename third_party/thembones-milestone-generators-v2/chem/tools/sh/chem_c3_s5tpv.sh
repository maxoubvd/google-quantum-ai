#!/bin/bash
# Fifth holder (Nov 18 05:00Z) VARIANT chem_c3_s5tpv.sh: DEADUSE=own|prev|either (R66 mask4dead --deaduse); with prev the rebuilt
# 07:50Z: P1COMBOS (phase-1 'slotpoint seed' pairs separated by |, e.g. "1 0|2 0|1 2|2 2") TPTRIES (default 60) and QSAMPLES (default 1500000) as parameters.
# 07:48Z: VPROJ=<proj.json> = the block set V for loader K taken from that file (r58_chem_kept_vpad.py --proj): another valid V (ledger O29).
# neighbour is K-1 (batch K-1 changed), else K+1. Original header follows.
# Builder C3 (fourth holder, Nov 17 19:55Z): ON A BOX. Design seat R66-CHEM-PLAN244's "S5T+" (S5T + the dead words of
# the loader's own batches + the staircase stage A with its slot: 250 = 74 + 176 [V on G9s2's k1, c3_g9s2tp65c]) for
# loader K (normally 1: batches 0, 1 are never an owner's) of a table, the seat's tools read-only and every file in MY
# folder (RUNDIR/s5tp); then the neighbour K+1 rebuilt on the changed words (R59's g42 driver, a plain wired loader).
# Steps as the seat's 'For C3' block: vpad + my basis search (chem_c3_bsrun.sh, its S5 = the fallback), staircase
# stage A (four last blocks, the first at 149), phase 1 (slot point 14 / 15 x seeds), the full plan on N processes
# sharing the lucky classes, r62_chem_cache_rp_apply.py + r66_chem_plan244_emit4dead.py, words copied into the table.
# Writes RUNDIR/s5tp/{angle_k_K.kmx, words.npz, n<K+1>/, DONE ('ok <toffoli>' or 'fail <step>')}.
#   CORES="17 18 19" bash chem_c3_s5tp.sh RUNDIR TABLE K WIRED      (WIRED: the table's wired loaders, e.g. 2,4,6,8,9)
set -u
R=$1/s5tp; TB=$2; K=$3; WIRED=$4; CO=(${CORES:-17 18 19}); N=${#CO[@]}
PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools; SECS=${S5TP_SECS:-1500}
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $R; cd $R || exit 1; rm -f DONE
log() { echo "$(date -u +%H:%M:%SZ) k$K $*" >> $R/s5tp.log; }
fail() { log "FAIL $1"; echo "fail $1" > $R/DONE; exit 0; }
log "start: table $TB, wired $WIRED, cores ${CO[*]}"
mkdir -p vin
$PY -I $T/r58_chem_kept_vpad.py $TB --ks $K ${VPROJ:+--proj $VPROJ} --out vin/pad.npz > vpad.log 2>&1
python3 -c "import json, sys; json.dump({'V': json.load(open(sys.argv[1]))['V']}, open(sys.argv[2], 'w'))" vin/pad.npz.json vin/proj.json \
  || fail vpad
V=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" vin/proj.json $K)
# the searched basis only (tools/chem_c3_bsearch.py, as chem_c3_bsrun.sh's first step; seconds)
mkdir -p b$K/in_k$K; ln -sf $R/vin/pad.npz b$K/in_k$K/pad.npz
nice -n 10 taskset -c ${CO[0]} $PY -I $T/chem_c3_bsearch.py b$K/in_k$K/pad.npz --k $K --D 13,14,15 --tries 300000 --seed 1 \
  --out $R/b$K/in_k$K/k${K}bs1_D13-14-15.json > bsearch.log 2>&1
PF=$(ls $R/b$K/in_k$K/k${K}bs1_D13-14-15.json 2>/dev/null)
[ -n "$PF" ] || fail basis
log "basis $PF"
# staircase stage A on the four last blocks (10 to 25 s each), the first at 149 exact
i=0
for st in "0,2,3|1" "0,1,2|3" "0,1,3|2" "1,2,3|0"; do
  t=sa_$(echo "$st" | tr -d ',|')
  ( nice -n 10 taskset -c ${CO[$((i % N))]} timeout 300 $PY -I $T/r66_chem_plan244_stairA.py vin/pad.npz --k $K --pfile $PF --rep 4 \
      --stages "$st" --kills 2:2 --seed 0 --out $t.json > $t.log 2>&1 ) &
  i=$((i + 1)); [ $((i % N)) = 0 ] && wait
done
wait
SA=""
for t in sa_0231 sa_0123 sa_0132 sa_1230; do
  grep -q "STAGE A 149" $t.log 2>/dev/null && grep -q "outside D: True" $t.log && [ -s $t.json ] && { SA=$t.json; break; }
done
[ -n "$SA" ] || fail "stairA ($(grep -h -o 'STAGE A [0-9]*' sa_*.log 2>/dev/null | tr '\n' ' '))"
log "staircase stage A: $SA"
PH="--planA $SA --K 8 --onlyK --pseeds 1 --dist 5,1,1 --V $V --dc --eqw --s5 --tp --dead --deaduse ${DEADUSE:-own} --wired $WIRED --F 11 --c12 29 --c3 20"
# phase 1: slot point 1 (= 14) / 2 (= 15) x seeds 0, 1: the first CHECKED EXACT
i=0
IFS="|" read -r -a P1L <<< "${P1COMBOS:-1 0|2 0|1 1|2 1}"
for sp in "${P1L[@]}"; do
  set -- $sp
  ( nice -n 10 taskset -c ${CO[$((i % N))]} timeout 600 $PY -I $T/r66_chem_plan244_mask4dead.py vin/pad.npz $PH --seed $2 --slots --slotpt $1 \
      --tptries ${TPTRIES:-60} --qsamples ${QSAMPLES:-1500000} --phase1only > p1_$1_$2.log 2>&1 ) &
  i=$((i + 1)); [ $((i % N)) = 0 ] && wait
done
wait
SPT=""; SEED=""
IFS="|" read -r -a P1L <<< "${P1COMBOS:-1 0|2 0|1 1|2 1}"
for sp in "${P1L[@]}"; do
  set -- $sp
  grep -q "CHECKED EXACT" p1_$1_$2.log 2>/dev/null && { SPT=$1; SEED=$2; break; }
done
[ -n "$SPT" ] || fail "phase 1 ($(grep -h -o 'rank [0-9]*' p1_*.log 2>/dev/null | sort -k2 -n | tail -n 1))"
log "phase 1 exact: slot point $SPT, seed $SEED"
# the full plan: N processes, one seed, lseeds 0..N-1 sharing the lucky classes
mkdir -p share
l=0
for c in "${CO[@]}"; do
  ( nice -n 10 taskset -c $c timeout $SECS $PY -I $T/r66_chem_plan244_mask4dead.py vin/pad.npz $PH --seed $SEED --slots --slotpt $SPT \
      --tptries ${TPTRIES:-60} --qsamples ${QSAMPLES:-1500000} --s5bases 4000000 --lshare $R/share --lseed $l --plan plan_l$l.json --chg chg_l$l.json \
      > full_l$l.log 2>&1 ) &
  l=$((l + 1))
done
wait
best=999; BL=""
for l in $(seq 0 $((N - 1))); do
  [ -s plan_l$l.json ] && [ -s chg_l$l.json ] || continue
  $PY -I $T/r62_chem_cache_rp_apply.py vin/pad.npz chg_l$l.json --out pad_mod_l$l.npz --noverify --noneighbours > apply_l$l.log 2>&1 || continue
  $PY -I $T/r66_chem_plan244_emit4dead.py vin/pad.npz plan_l$l.json --padmod pad_mod_l$l.npz --V $V --out k${K}_l$l.kmx > emit_l$l.log 2>&1
  n=$(grep -h -o 'EMITTED k [0-9]* on the 68-AND cache: [0-9]* Toffolis' emit_l$l.log | grep -o '[0-9]* Toffolis' | awk '{print $1}')
  log "lseed $l: $(grep -h -E 'EMITTED|rror' emit_l$l.log | tail -n 1 | cut -c1-160)"
  [ -n "$n" ] && [ "$n" -lt "$best" ] && { best=$n; BL=$l; }
done
[ -n "$BL" ] || fail "no plan emitted ($(grep -h '^SUMMARY' full_l*.log 2>/dev/null | tail -n 1 | cut -c1-120))"
cp k${K}_l$BL.kmx angle_k_$K.kmx
# the words: the changed words of the padded table copied into the table (dead words of batches K-1, K only)
$PY - $TB pad_mod_l$BL.npz words.npz $K << 'PYEOF' > words.log 2>&1 || fail words
import os, sys
import numpy as np
sys.path.insert(0, os.path.expanduser("~/Quantum/tools"))
from r62_chem_cache_dead import dead_mask
tb, pm, out, k = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
z = np.load(tb)
arrs = {a: z[a] for a in z.files}
w = np.asarray(arrs["angle_words"], np.int64)
p = np.asarray(np.load(pm)["angle_words"], np.int64)
n = w.shape[0]
assert p.shape[0] >= n and p.shape[1] == w.shape[1], (p.shape, w.shape)
d = p[:n] != w
dead = dead_mask(w, np.asarray(arrs["tree_edges"], np.int64), 7)
assert not (d & ~dead).any(), "S5T+ changed a word the decoder reads"
cols = [c for c in range(w.shape[1]) if d[:, c].any()]
assert all(c // 5 in (k - 1, k) for c in cols), f"changed columns outside batches {k - 1}, {k}: {cols}"
w2 = w.copy()
w2[d] = p[:n][d]
arrs["angle_words"] = w2.astype(z["angle_words"].dtype)
np.savez(out, **arrs)
print("dead words changed per batch", [int(d[:, 5 * b:5 * b + 5].sum()) for b in range(15)])
PYEOF
log "words: $(tail -n 1 words.log)"
K2=$((K + 1)); [ "${DEADUSE:-own}" = prev ] && K2=$((K - 1))
( cd $T/.. && nice -n 10 $PY -I tools/r59_chem_groups_g42.py --table $R/words.npz --k $K2 --seeds 0-63 --cores ${CO[0]},${CO[$((1 % N))]} \
    --rounds 6 --secs 600 --out $R/n$K2 > $R/n$K2.log 2>&1 )
log "neighbour k$K2: $(tail -n 1 n$K2/RESULT.txt 2>/dev/null | cut -c1-120)"
{ ls n$K2/s*/angle_k_$K2.kmx > /dev/null 2>&1 || ls n$K2/angle_k_$K2.kmx > /dev/null 2>&1; } || fail "neighbour k$K2"
echo "ok $best" > DONE
log "ok $best"
