#!/usr/bin/env bash
# Builder C3: ONE kept-product loader on the 68-AND cache = design seat R58-CHEM-KEPT's three steps (its tools,
# imported read-only): stage A by C2's engine (r58_chem_kept_c2a.py), the point stage with 8 helpers under the C-4
# conditions (r58_chem_kept_c4chk.py --V), the 68-AND emitter (r58_chem_kept_emit68.py; it checks every wire on
# every address). Local WSL.
#   bash chem_c3_one.sh IN K PFILE REP D OUTDIR [PSEEDS]
# IN = C2's in-folder (pad.npz, proj.json, kK*.json; read-only). Writes OUTDIR/<tag>.log, planA_<tag>.json,
# plan68_<tag>.json, <tag>.kmx; prints the last line of the log. Exit 3: singular basis; 4: no 256 emitted.
set -u
IN=$1; K=$2; PF=$3; REP=$4; D=$5; OUT=$6; PS=${7:-2}
PY=${PY:-$HOME/qcc/ws/chem_own/venv/bin/python}
T=${CT:-<HOME>/Quantum/tools}
export CHEM_BROT=7 CHEM_HB=77
mkdir -p "$OUT"
ASEED=${ASEED:-0}; AEFF=${AEFF:-1}
TAG=k${K}_${PF%.json}_r${REP}_D${D//,/-}
if [ "$ASEED" != "0" ] || [ "$AEFF" != "1" ]; then TAG=${TAG}_s${ASEED}e${AEFF}; fi
LOG=$OUT/$TAG.log
V=$(python3 -c "import json,sys; v=json.load(open(sys.argv[1]))['V'].get(sys.argv[2]); print(v if v else '')" "$IN/proj.json" "$K")
if [ -z "$V" ]; then echo "$TAG: loader $K has no block set V in proj.json"; exit 5; fi
{
  echo "V of loader $K: $V"
  $PY -I $T/r58_chem_kept_c2a.py "$IN/pad.npz" --k "$K" --pfile "$IN/$PF" --rep "$REP" --D "$D" --seed "$ASEED" --effort "$AEFF" --out "$OUT/planA_$TAG.json"
} > "$LOG" 2>&1
if [ ! -s "$OUT/planA_$TAG.json" ]; then echo "$TAG: $(tail -1 "$LOG")"; exit 3; fi
$PY -I $T/chem_c3_c4chk.py "$IN/pad.npz" --planA "$OUT/planA_$TAG.json" --K 8 --onlyK --pseeds "$PS" --V "$V" \
    --plan "$OUT/plan68_$TAG.json" 2>&1 | grep -v "^   point" | cut -c1-300 >> "$LOG"
$PY -I $T/chem_c3_emit68.py "$IN/pad.npz" "$OUT/plan68_$TAG.json" --V "$V" --out "$OUT/$TAG.kmx" ${EMITX:-} >> "$LOG" 2>&1
L=$(tail -1 "$LOG")
echo "$TAG: $L"
case "$L" in *"EMITTED"*) exit 0;; esac
exit 4
