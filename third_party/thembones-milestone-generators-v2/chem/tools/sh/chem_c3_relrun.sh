#!/usr/bin/env bash
# Builder C3 (second holder, Nov 17): ONE RELATION-BEARING WIRED LOADER (a wired loader whose table carries bit-1
# relations: 24 or 23 generic outputs, WR7's P51 / P52 layouts) through R59-CHEM-GROUPS' tie engine with the closure,
# on the plan of its generic count (C2's tie step needs a 7-group at group 2 and does not take these plans).
# Writes OUT/relrun_k<K>.txt with a line 'rc=0 <tag>: EMITTED k K ...: N Toffolis ..., qubits Q' (the watcher's
# format), OUT/<tag>.kmx, then OUT/fin_k<K>.  Run pinned: taskset -c <core> nice -n 10 bash chem_c3_relrun.sh ...
#   bash chem_c3_relrun.sh TABLE(unpadded sub.npz) K PLAN MODES STOP OUT [SEEDS]
set -u
TB=$1; K=$2; PL=$3; MO=$4; ST=$5; OUT=$6; DRV=${7:-tie}; WANTG=${8:--}; CORE=${9:-0}; SEEDS=${SEEDS:-0-23}
PY=${PY:-$HOME/qcc/ws/chem_own/venv/bin/python}
T=${CT:-<HOME>/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1
mkdir -p "$OUT/rel_k$K"
R=$OUT/relrun_k$K.txt
TAG=relk${K}
if [ "$DRV" = "g42" ]; then
  # R59's driver: ties + C2's annealed 7-group at group 2 (the plan must have one there), one core
  $PY -I $T/r59_chem_groups_g42.py --table "$TB" --k "$K" --plan "$PL" --modes "$MO" --wantg "$WANTG" --P 42 \
      --seeds "$SEEDS" --cores "$CORE" --rounds 4 --secs 150 --out "$OUT/rel_k$K" > "$OUT/rel_k$K.log" 2>&1
  KMX=$OUT/rel_k$K/s1/angle_k_$K.kmx
  [ -s "$KMX" ] || KMX=$OUT/rel_k$K/angle_k_$K.kmx
  if [ ! -s "$KMX" ]; then
    # 12:46Z: the 23-output plan with the 7-group does not keep its kills on every layout (RK7 run 24's k13: 39 | 38 |
    # 46 | 28 | 7): fall back to the tie engine on C2's plan of 23 outputs (one-output groups; closed 226 there)
    FBP=${FBPLAN:-6,6,1:1:0:0,1:1:0:1,1:1:0:0,1:1:0:1,7:1:0:0}; FBM=${FBMODES:-reg,tie3,-,-,-,-,-}
    $PY -I $T/r59_chem_groups_tie.py "$TB" --k "$K" --plan "$FBP" --modes "$FBM" --seeds "$SEEDS" --pad --stop 148 \
        --emit "$OUT/rel_k${K}fb" 2>&1 | grep -E "^k |BEST|emitted|no plan|rror" | cut -c1-200 >> "$OUT/rel_k$K.log"
    KMX=$OUT/rel_k${K}fb/angle_k_$K.kmx
  fi
else
  $PY -I $T/r59_chem_groups_tie.py "$TB" --k "$K" --plan "$PL" --modes "$MO" --seeds "$SEEDS" --pad --stop "$ST" \
      --emit "$OUT/rel_k$K" 2>&1 | grep -E "^k |BEST|emitted|no plan|rror" | cut -c1-200 > "$OUT/rel_k$K.log"
  KMX=$OUT/rel_k$K/angle_k_$K.kmx
fi
if [ -s "$KMX" ]; then
  cp "$KMX" "$OUT/$TAG.kmx"
  # Toffolis and qubits of the emitted circuit, read by the simulator's own loader
  LINE=$($PY -I - "$OUT/$TAG.kmx" <<'EOF'
import os, sys
for p in (os.path.expanduser("~/qcc/data/kickmix"),):
    sys.path.insert(0, p)
import kickmix as km
c = km.Circuit.from_file(sys.argv[1])
print(int(c.operation_statistics().num_toffoli_like), int(c.num_qubits))
EOF
)
  set -- $LINE
  echo "$(date -u +%H:%M:%SZ) rc=0 $TAG: EMITTED k $K relation-wired (R59 tie engine, plan $PL): ${1:-0} Toffolis, qubits ${2:-999}" >> "$R"
else
  echo "$(date -u +%H:%M:%SZ) END k$K: no closed circuit ($(tail -n 1 "$OUT/rel_k$K.log"))" >> "$R"
fi
echo FIN > "$OUT/fin_k$K"
