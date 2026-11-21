#!/usr/bin/env bash
# S5-INV-H1END: FC's hold-route pipeline (a copy of tools/sh/invsqrt_fc_s19hold.sh, unchanged steps) on the PATCHED generator
# copy srcB (tools/s5_inv_h1end_genpatch.py: H1's word un-add lends q11 and runs at one more clean wire), env S5LEND=1.
# Output ONLY under ~/qcc/ws/s5/inv_h1end/files/NAME{,_e,_x,_w,_w.aff,_w.pc,_w.r91}.
#   bash tools/s5_inv_h1end_build.sh NAME [CORES=11]
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/s5/inv_h1end/files
NAME=$1; S=~/qcc/ws/s5/inv_h1end/srcB; CFG=~/qcc/ws/invsqrt_fc/s18/cfg/cfg_b312g.json; CORES=${2:-11}
grep -q S5LEND $S/invsqrt_fc2_patch_h1c.py || { echo "srcB not patched"; exit 1; }
O=$F/$NAME; WD=${O}_w
mkdir -p $O $WD; rm -f $O/invsqrt.kmb $WD/invsqrt.kmx
export S5LEND=${S5LEND:-1}
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=${FCGEN:-7,3,0;9,3,1}" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/s5h1nb_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== emit $(date -u +%H:%M:%SZ) S5LEND=$S5LEND"
(cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S $O $O/angles.json $O/cfg.json 2>&1 \
  | grep -E "S5LEND|R55 net|emit |per-tag|r63:|Traceback|Error|assert" | cut -c1-260)
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
