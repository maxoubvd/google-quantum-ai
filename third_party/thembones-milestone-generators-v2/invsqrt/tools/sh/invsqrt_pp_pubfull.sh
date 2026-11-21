#!/bin/bash
# P1 (sixth holder): the public exact rewrites to a fixed point on one invsqrt circuit:
#   repeat { invsqrt_pp_pubpass.py (P2, c3, c2 rounds) ; invsqrt_pp_pubc3q.py (dead writes paid by a quadratic form) }
#   until neither changes the file.  Every step checked on all 3,133,452 inputs; the final file checked against the input.
#   bash tools/sh/invsqrt_pp_pubfull.sh <invsqrt.kmb> <fresh out dir> [cores]   -> <out>/invsqrt.kmb, full.log
set -u
IN=$(readlink -f "$1"); OUT=$2; C=${3:-1}
if [ -d <HOME>/Quantum/tools ]; then T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
else T=~/Quantum/tools; PY=~/qcc/.venv_work/bin/python; fi
[ -e "$OUT/invsqrt.kmb" ] && { echo "$OUT/invsqrt.kmb exists"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
P="nice -n 10 taskset -c $C $PY"
cp "$IN" cur.kmb
for it in 1 2 3 4 5 6; do
  ch=0
  $P $T/invsqrt_pp_pubpass.py cur.kmb $PWD/a.kmb 6 > pass_$it.log 2>&1
  if grep -q "FINAL.*CHECK against the input PASS" pass_$it.log; then mv a.kmb cur.kmb; ch=1; fi
  rm -f a.kmb a.kmb.r*
  $P $T/invsqrt_pp_pubc3q.py cur.kmb $PWD/b.kmb > c3q_$it.log 2>&1
  if grep -q "CHECK PASS" c3q_$it.log; then mv b.kmb cur.kmb; ch=1; fi
  echo "iteration $it: $(grep -h -E 'FINAL|no change' pass_$it.log | sed 's/.*Toffoli-class/pass/') | c3q: $(tail -1 c3q_$it.log)"
  [ $ch = 1 ] || break
done
$P -c "
import sys; sys.path.insert(0, '$T')
from invsqrt_pp_pubpass import check
ok, a, b = check('$IN', 'cur.kmb'); print('FULL: expected Toffoli-class %.3f -> %.3f; CHECK against the input %s' % (a, b, 'PASS' if ok else 'FAIL'))" | tee full.log
grep -q "CHECK against the input PASS" full.log && mv cur.kmb invsqrt.kmb && sha256sum invsqrt.kmb | cut -c1-16
