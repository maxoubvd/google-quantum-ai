#!/usr/bin/env bash
# FC nineteenth holder (copy of invsqrt_fc_s18hold.sh, MODE relabel only), LOCAL WSL: the hold-route generator end to end, plus
# the R91 SPLIT-DEBT PASS (cfg key "r91_splitdebt", tools/invsqrt_fc_s19r91pass.py) after FC's post-chain:
#   emit (R64's xcrun, r59b0 router + angles) -> GHOST early q18-q21 -> I39 X-pair -> I39 relabel -> ALL-input check against
#   the emit -> .kmb -> P1's affine pass -> FC's post-chain (final.kmb) -> R91 pass (NAME_w.r91/invsqrt.kmb) -> ALL-input check
#   of the R91 output against the post-chain's text (r63_inv_ghost_chk.py, four outcome modes; must PASS).
#   bash tools/sh/invsqrt_fc_s19hold.sh NAME SRC CFG CORES
# Writes ~/qcc/ws/invsqrt_fc/s18/files/NAME{,_e,_x,_w}/, NAME_w.aff, NAME_w.pc, NAME_w.r91. Other seats' tools are only run.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/invsqrt_fc/s18/files
NAME=$1; S=$2; CFG=$3; CORES=$4
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
echo "== width step (relabel) $(date -u +%H:%M:%SZ)"
mkdir -p ${O}_e ${O}_x
(cd $Q/tools && $P r63_inv_ghost_early.py $O/invsqrt.kmx ${O}_e/invsqrt.kmx q18,q19,q20,q21 2>&1 | tail -1)
(cd $Q/tools && $P inv_i39_xpair.py ${O}_e/invsqrt.kmx ${O}_x/invsqrt.kmx 2>&1 | tail -1 | cut -c1-80)
(cd $Q/tools && nice -n 10 taskset -c $CORES $P inv_i39_relabel.py ${O}_x/invsqrt.kmx $WD/invsqrt.kmx 2>&1 | tail -1)
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
echo "post-chain final: $FIN  sha $(sha256sum $FIN | cut -c1-16)  $(date -u +%H:%M:%SZ)"
echo "== R91 split-debt pass $(date -u +%H:%M:%SZ)"
rm -rf ${WD}.r91
(cd $Q/tools && nice -n 10 taskset -c $CORES $P invsqrt_fc_s19r91pass.py $FIN $O/cfg.json ${WD}.r91 2>&1 | tail -6 | cut -c1-240)
[ -f ${WD}.r91/invsqrt.kmb ] || { echo "R91 PASS FAILED"; exit 1; }
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check ${WD}.r91/x.kmx ${WD}.r91/invsqrt.kmx 0 2>&1 | tail -2 | cut -c1-200) | tee ${WD}.r91/chk.log
grep -q "CHECK PASS" ${WD}.r91/chk.log || { echo "R91 ALL-INPUT CHECK FAILED"; exit 1; }
echo "final: ${WD}.r91/invsqrt.kmb  sha $(sha256sum ${WD}.r91/invsqrt.kmb | cut -c1-16)  $(date -u +%H:%M:%SZ)"
