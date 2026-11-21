#!/usr/bin/env bash
# FC seventeenth holder (local WSL): pack a chained hold-route variant into ~/qcc/ws/invsqrt_fc/best/fc17_hold<NAME>
# (invsqrt.kmb = the post-chain's final, cfg, angles, SOURCE.txt = the given note + the parent's SOURCE.txt) and queue
# the official verifier on both salts (box 4 lane; log ~/qcc/ws/invsqrt_fc/s17/logs/verify_<NAME>.log).
#   bash tools/sh/invsqrt_fc_s17pack.sh NAME PARENT_BEST_DIR "one-line note on what changed" "reproduce line"
set -eu
Q=<HOME>/Quantum
NAME=$1; PAR=$2; NOTE=$3; REPRO=$4
F=~/qcc/ws/invsqrt_fc/s17/files; B=~/qcc/ws/invsqrt_fc/best/fc17_hold$NAME
FIN=$F/${NAME}_w.pc/${NAME}_w.aff/final.kmb
[ -f $FIN ] || { echo "no final: $FIN"; exit 1; }
mkdir -p $B; cp $FIN $B/invsqrt.kmb; cp $F/$NAME/cfg.json $F/$NAME/angles.json $B/
{ echo "invsqrt: FILE $NAME (FC seventeenth holder, $(date -u +%Y-%m-%d\ %H:%MZ)): $NOTE"
  echo "CHAIN: emit (r64_inv_stack_xcrun.py, r59b0 router and angles) -> GHOST early q18-q21 -> I39 X-pair -> I39 relabel (39 qubits;"
  echo "  r63_inv_ghost_chk.py against the emit on ALL 3,133,452 inputs, outcomes zero / one / random 1 / random 2: CHECK PASS) -> P1 affine"
  echo "  pass (CHECK PASS) -> FC post-chain (CHECK PASS; T fold)."
  echo "REPRODUCE: $REPRO"
  echo "  then: bash tools/sh/box_verify.sh $B invsqrt fc17_hold$NAME"
  echo; echo "--- parent: $PAR ---"; cat $PAR/SOURCE.txt; } > $B/SOURCE.txt
echo "packed $B sha $(sha256sum $B/invsqrt.kmb | cut -c1-16)"
cd $Q && setsid nohup bash tools/sh/box_verify.sh $B invsqrt fc17_hold$NAME > ~/qcc/ws/invsqrt_fc/s17/logs/verify_$NAME.log 2>&1 < /dev/null &
echo "verify queued: ~/qcc/ws/invsqrt_fc/s17/logs/verify_$NAME.log"
