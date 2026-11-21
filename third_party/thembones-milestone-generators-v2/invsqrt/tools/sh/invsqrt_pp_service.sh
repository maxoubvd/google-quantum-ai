#!/usr/bin/env bash
# Workstream invsqrt_pp (P1): THE SERVICE CHAIN for one invsqrt circuit (FC's family C and any other), exact passes only.
#   bash tools/sh/invsqrt_pp_service.sh <circuit invsqrt.kmb> <fresh out dir> [core]
# 1 affine audit + rewrite (invsqrt_pp_affine.py, invsqrt_pp_affwrite.py)       -> s1.kmb (or the input if no site)
# 2 squeeze: drop idle values (affine rebuilds, then one-CCX rebuilds if SQ_P1=1), re-allocate wires; repeated
#   while it finds a child (invsqrt_pp_squeeze.py)                               -> s2_<k>.kmb
# 3 exact links + hopping holders (invsqrt_pp_links.py, invsqrt_pp_hop.py)        -> final.kmb
# Each step writes its child only after its own all-input check (state, Z-phase on 4 transcripts, every Z_POW).
# Works on this PC (WSL) and on the boxes (tools under ~/Quantum/tools there). Then: tools/sh/box_verify.sh.
set -u
IN=$(readlink -f "${1:?circuit.kmb}"); OUT=${2:?out dir}; CORE=${3:-1}
if [ -d <HOME>/Quantum/tools ]; then T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
else T=~/Quantum/tools; PY=~/qcc/.venv_work/bin/python; fi
mkdir -p "$OUT" && cd "$OUT" || exit 1
P="nice -n 10 taskset -c $CORE $PY"
cur=$IN
$P $T/invsqrt_pp_affine.py "$cur" affine.json > affine.log 2>&1
$P $T/invsqrt_pp_affwrite.py "$cur" affine.json s1.kmb.tmp > affwrite.log 2>&1
if grep -q "CHECK PASS" affwrite.log; then mv s1.kmb.tmp s1.kmb; cur=$PWD/s1.kmb; grep -E "^expected" affwrite.log; else rm -f s1.kmb.tmp; echo "affine: no child"; fi
k=0
while :; do
  k=$((k+1))
  SQ_P1=${SQ_P1:-0} $P $T/invsqrt_pp_squeeze.py "$cur" s2_$k.kmb > squeeze_$k.log 2>&1
  grep -E "^declared|^chosen" squeeze_$k.log
  if [ -f s2_$k.kmb ]; then cur=$PWD/s2_$k.kmb; else break; fi
done
$P $T/invsqrt_pp_links.py "$cur" links.json 0 > links.log 2>&1
$P $T/invsqrt_pp_hop.py "$cur" links.json final.kmb > hop.log 2>&1
grep -E "^admitted|CHECK" hop.log
[ -f final.kmb ] || cp "$cur" final.kmb
echo "final $(sha256sum final.kmb | cut -c1-16)"
