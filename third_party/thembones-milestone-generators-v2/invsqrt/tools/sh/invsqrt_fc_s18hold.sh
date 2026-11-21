#!/usr/bin/env bash
# FC eighteenth holder (copy of invsqrt_fc_s17hold.sh, s18 folders), LOCAL WSL: a hold-route 39-qubit variant end to end (copy of invsqrt_fc16_holdB.sh with a
# choice of width step and the all-input check before P1):
#   emit (R64's xcrun, r59b0 router + angles) -> WIDTH STEP -> ALL-input check against the emit (r63_inv_ghost_chk.py,
#   four outcome modes) -> .kmb -> P1's affine pass -> FC's post-chain (-> final.kmb).
#   WIDTH STEP: MODE=relabel: GHOST's early tool q18-q21 -> I39's X-pair peephole -> I39's multi-outcome relabel;
#               MODE=compact: FILE38's compaction pass (r65_inv_file38_compact.py: by-value liveness, moves, packing).
#   bash tools/sh/invsqrt_fc_s18hold.sh NAME SRC CFG CORES MODE
# Writes ~/qcc/ws/invsqrt_fc/s17/files/NAME{,_w}/ and NAME_w.aff, NAME_w.pc.  Other seats' tools are only run.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/invsqrt_fc/s18/files
NAME=$1; S=$2; CFG=$3; CORES=$4; MODE=${5:-relabel}
O=$F/$NAME; WD=${O}_w
mkdir -p $O $WD; rm -f $O/invsqrt.kmb $WD/invsqrt.kmx
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=7,3,0;9,3,1" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/fc18nb_h_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== emit $(date -u +%H:%M:%SZ)"
(cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S $O $O/angles.json $O/cfg.json 2>&1 \
  | grep -E "R55 net|emit |per-tag|r63:|Traceback|Error|assert" | cut -c1-260)
[ -f $O/invsqrt.kmx ] || { echo "EMIT FAILED (no .kmx)"; exit 1; }
echo "== width step ($MODE) $(date -u +%H:%M:%SZ)"
if [ "$MODE" = compact ]; then
  (cd $Q/tools && nice -n 10 taskset -c $CORES $P r65_inv_file38_compact.py $O/invsqrt.kmx $WD/invsqrt.kmx 2>&1 | tail -4 | cut -c1-260)
else
  mkdir -p ${O}_e ${O}_x
  (cd $Q/tools && $P r63_inv_ghost_early.py $O/invsqrt.kmx ${O}_e/invsqrt.kmx q18,q19,q20,q21 2>&1 | tail -1)
  (cd $Q/tools && $P inv_i39_xpair.py ${O}_e/invsqrt.kmx ${O}_x/invsqrt.kmx 2>&1 | tail -1 | cut -c1-80)
  (cd $Q/tools && nice -n 10 taskset -c $CORES $P inv_i39_relabel.py ${O}_x/invsqrt.kmx $WD/invsqrt.kmx 2>&1 | tail -1)
fi
[ -f $WD/invsqrt.kmx ] || { echo "WIDTH STEP FAILED"; exit 1; }
echo "== all-input check $(date -u +%H:%M:%SZ)"
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $O/invsqrt.kmx $WD/invsqrt.kmx 0 2>&1 | tail -2 | cut -c1-200)
(cd $Q/tools && $P -c "import invsqrt_fw_rw as R; R.to_kmb('$WD/invsqrt.kmx', '$WD/invsqrt.kmb')" 2>&1 | tail -2)
cp $O/cfg.json $O/angles.json $WD/
echo "== P1 pass + post-chain $(date -u +%H:%M:%SZ)"
cd $Q && bash tools/sh/invsqrt_pp_affpass.sh $WD/invsqrt.kmb ${WD}.aff $CORES 2>&1 | grep -E "expected Toffoli-class|sites [0-9]+, rewritten|CHECK" | tr '\n' ' '; echo
A=${WD}.aff/invsqrt.kmb; [ -f $A ] || A=$WD/invsqrt.kmb
bash tools/sh/invsqrt_fc_s7pubchain.sh ${WD}.pc $CORES $A 2>&1 | tail -2
FIN=${WD}.pc/$(basename $(dirname $A))/final.kmb
echo "final: $FIN  sha $(sha256sum $FIN | cut -c1-16)  $(date -u +%H:%M:%SZ)"
