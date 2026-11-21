#!/usr/bin/env bash
# R66-CHEM-PLAN244 (box): S5 + the two-point helper on the 8th wire (tools/r66_chem_plan244_mask3.py --tp) on a GIVEN
# stage A and padding, then the emitter (tools/r66_chem_plan244_emit3.py). One process; several processes with the
# same SEED and different LSEED share the lucky classes of phase 2 through out/share_<TAG>/.
# Output: ~/qcc/ws/r66_chem_plan244/out/{plan,kmx,log}_<TAG>_l<LSEED>.*
#   bash r66_chem_plan244_job.sh PAD PROJ_JSON PLANA LOADER SEED LSEED TAG [F] [C12] [C3] [EXTRA]
set -u
PAD=$1; PROJ=$2; PLANA=$3; K=$4; SEED=$5; LSEED=$6; TAG=$7; FN=${8:-10}; C12=${9:-31}; C3=${10:-20}; EXTRA=${11:-}
R=~/qcc/ws/r66_chem_plan244; PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=77
mkdir -p $R/out/share_$TAG
V=$(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" $PROJ $K)
O=$R/out/plan_${TAG}_l$LSEED.json
echo "V of loader $K: $V   planA $(basename $PLANA)   seed $SEED lseed $LSEED   TP F $FN c12 $C12 c3 $C3   $(date -u +%H:%M:%SZ)"
$PY -I $T/r66_chem_plan244_mask3.py $PAD --planA $PLANA --K 8 --onlyK --seed $SEED --pseeds 1 --dist 5,1,1 --V $V \
    --dc --eqw --s5 --tp --F $FN --c12 $C12 --c3 $C3 --s5bases 4000000 --lshare $R/out/share_$TAG --lseed $LSEED \
    $EXTRA --plan $O 2>&1 | grep --line-buffered -E "^SUMMARY|cells emptied per point|^RP|^R62|S5|TP|C-4|^k |skipped|NO PLAN|rror|ssert|Trace"
[ -s $O ] && $PY -I $T/r66_chem_plan244_emit3.py $PAD $O --V $V --out $R/out/kmx_${TAG}_l$LSEED.kmx
echo "END $(date -u +%H:%M:%SZ)"
