#!/bin/bash
# P1 (sixth holder): full exact post-chain on invsqrt files: public fixed point -> service chain (affine, squeeze,
# links + hops) -> public fixed point.  Every step all-input checked.  Prints the fair-coin Toffoli-class delta.
#   bash tools/sh/invsqrt_pp_pubchain.sh <out root> <cores> <invsqrt.kmb> [more .kmb ...]   -> <out root>/<parent dir name>/final.kmb
set -u
ROOT=$1; C=$2; shift 2
T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
for f in "$@"; do
  f=$(readlink -f "$f"); n=$(basename $(dirname "$f")); O=$ROOT/$n; rm -rf "$O"; mkdir -p "$O"
  bash $T/sh/invsqrt_pp_pubfull.sh "$f" "$O/p1" "$C" > "$O/p1.log" 2>&1
  a=$f; [ -f "$O/p1/invsqrt.kmb" ] && a=$O/p1/invsqrt.kmb
  SQ_P1=1 bash $T/sh/invsqrt_pp_service.sh "$a" "$O/svc" "$C" > "$O/svc.log" 2>&1
  b=$a; [ -f "$O/svc/final.kmb" ] && b=$O/svc/final.kmb
  bash $T/sh/invsqrt_pp_pubfull.sh "$b" "$O/p2" "$C" > "$O/p2.log" 2>&1
  z=$b; [ -f "$O/p2/invsqrt.kmb" ] && z=$O/p2/invsqrt.kmb
  cp "$z" "$O/final.kmb"
  nice -n 10 taskset -c $C $PY -c "
import sys; sys.path.insert(0, '$T')
from invsqrt_pp_pubpass import check
ok, a, b = check('$f', '$O/final.kmb'); print('$n: fair-coin Toffoli-class %.3f -> %.3f (%.3f); CHECK against the input %s' % (a, b, b - a, 'PASS' if ok else 'FAIL'))" | tee "$O/chain.txt"
done
