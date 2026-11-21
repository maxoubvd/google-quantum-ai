#!/usr/bin/env bash
# FC sixteenth holder, LOCAL WSL: the hold-route 39-qubit FILE by I39's FILE-A chain (I39 16:24Z), run in FC's folder:
# emit (R64's xcrun, I39's env: R41_XC=DIV1, no R56_LEV, r59b0's router and angles) -> GHOST's early tool (q18-q21 on
# the .kmx) -> I39's X-pair peephole -> I39's multi-outcome relabel (qubits = peak live width) -> .kmx -> .kmb ->
# P1's affine pass -> FC's post-chain.  Tools of other seats are only READ and run (I39, R63-GHOST, R64, P1, F1).
#   bash tools/sh/invsqrt_fc16_holdB.sh NAME SRC CFG CORES
# Writes ~/qcc/ws/invsqrt_fc/s16/files/NAME{,_e,_x,_r}/ and NAME_r.aff, NAME_r.pc; prints the final .kmb.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
R59=~/qcc/ws/invsqrt_fc/best/r59b0; F=~/qcc/ws/invsqrt_fc/s16/files
NAME=$1; S=$2; CFG=$3; CORES=$4
O=$F/$NAME
mkdir -p $O ${O}_e ${O}_x ${O}_r; rm -f $O/invsqrt.kmb
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=7,3,0;9,3,1" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/fc16nb_hb_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== emit $(date -u +%H:%M:%SZ)"
(cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S $O $O/angles.json $O/cfg.json 2>&1 \
  | grep -E "R55 net|emit |per-tag|r63:|Traceback|Error|assert" | cut -c1-260)
[ -f $O/invsqrt.kmx ] || { echo "EMIT FAILED (no .kmx)"; exit 1; }
echo "== early / xpair / relabel $(date -u +%H:%M:%SZ)"
(cd $Q/tools && $P r63_inv_ghost_early.py $O/invsqrt.kmx ${O}_e/invsqrt.kmx q18,q19,q20,q21 2>&1 | tail -4)
(cd $Q/tools && $P inv_i39_xpair.py ${O}_e/invsqrt.kmx ${O}_x/invsqrt.kmx 2>&1 | tail -3)
(cd $Q/tools && nice -n 10 taskset -c $CORES $P inv_i39_relabel.py ${O}_x/invsqrt.kmx ${O}_r/invsqrt.kmx 2>&1 | tail -4)
[ -f ${O}_r/invsqrt.kmx ] || { echo "RELABEL FAILED"; exit 1; }
(cd $Q/tools && $P -c "import invsqrt_fw_rw as R; R.to_kmb('${O}_r/invsqrt.kmx', '${O}_r/invsqrt.kmb')" 2>&1 | tail -2)
cp $O/cfg.json $O/angles.json ${O}_r/
echo "== P1 pass + post-chain $(date -u +%H:%M:%SZ)"
cd $Q && bash tools/sh/invsqrt_pp_affpass.sh ${O}_r/invsqrt.kmb ${O}_r.aff $CORES 2>&1 | grep -E "expected Toffoli-class|sites [0-9]+, rewritten|CHECK" | tr '\n' ' '; echo
A=${O}_r.aff/invsqrt.kmb; [ -f $A ] || A=${O}_r/invsqrt.kmb
bash tools/sh/invsqrt_fc_s7pubchain.sh ${O}_r.pc $CORES $A 2>&1 | tail -2
FIN=${O}_r.pc/$(basename $(dirname $A))/final.kmb
echo "final: $FIN  sha $(sha256sum $FIN | cut -c1-16)  $(date -u +%H:%M:%SZ)"
