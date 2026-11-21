#!/bin/bash
# chem_c4_r86go25.sh: builder C4 (holder 6, Nov 19 05:35Z). ON A BOX OR LOCAL, one core. Design seat R86's group-4 form on a
# 25-OUTPUT wired loader (a wired loader with no bit-1 relation: G10P57x's k3, every wired loader of G11w11). Copy of
# design seat R89's strict tools/sh/r89_chem_cascade_r86go.sh (itself a copy of my chem_c4_r86go.sh) with two changes:
#  (1) the plan is the 25-output one, 6,6,7:1:0:0,6:1:0:1 (groups 1-2 the pool-free 38s, group 3 the 7-output group,
#      group 4 the 6-output last group with the column kill and NO row kills, R86's form);
#  (2) no group-3 row drop (the 7-group is not row-bound); its annealed 42-plan is injected when PL42 is given.
# Every R86 / R89 tool is used read-only. Log OUT/r86go.log; OUT/{angle_k_K.kmx, table_dead.npz} on success.
#   [PL42=<pairdir>/pl42_<r>.json] [R86_BATCHES=own|both|prev] [CENSUS_ONLY=1] bash chem_c4_r86go25.sh K TABLE PAIR OUT [CORE]
K=$1; TB=$2; PAIR=$3; OUT=$4; CORE=${5:-10}
mkdir -p $OUT
if [ -d <HOME>/Quantum/tools ]; then cd <HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
else cd ~/Quantum/tools || exit 1; PY=~/qcc/.venv_work/bin/python; fi
# R86's dump / hook attach through chem_c4_dead.install_row, which the engine calls only when C4_DEAD_ROW is set: the spec
# below installs that interface with its own drop aimed at a group that does not exist (9), so no group gets my drop and
# R86's wrapper answers group 4 (R86_GI=3). R89's strict wrapper is NOT used here (it would count every group-3 return
# as a miss).
export C4_DEAD_BATCHES=${R86_BATCHES:-both} CHEM_BROT=7 CHEM_HB=77 C4_DROP=1 C4_DEAD_ROW=3:9:72 C4_DEAD_ROW_TRIALS=100
INJ2=""; [ -n "${PL42:-}" ] && [ -s "$PL42" ] && INJ2="--inj 2=$PL42"
ARGS="--k $K --wirebasis --plan 6,6,7:1:0:0,6:1:0:1 $(cat $PAIR/args.txt) --inj 0=$PAIR/g0_best.json --inj 1=$PAIR/g1_best.json $INJ2 --modes=-,-,-,- --tries 3000"
say() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/r86go.log; }
say "start k$K (25 outputs) table $TB pair $PAIR ${INJ2:+with $INJ2}"
R86_GI=3 R86_DUMP=$OUT/k${K}_g4.npz nice -n 10 taskset -c $CORE $PY -I r86_chem_g35_dump.py $TB $ARGS > $OUT/dump.log 2>&1
[ -s $OUT/k${K}_g4.npz ] || { say "no dump: $(tail -n 2 $OUT/dump.log | tr '\n' ' ' | cut -c1-200)"; echo "fail dump" > $OUT/DONE; exit 2; }
nice -n 10 taskset -c $CORE $PY -I r86_chem_g35_census.py $OUT/k${K}_g4.npz --psis 16 --payer-min 0 > $OUT/census.log 2>&1
LAY=$($PY -c "
import json,sys
d = json.loads(open(sys.argv[1]).read().strip().splitlines()[-1])
f = d.get('first_found') or []
print(len(f), d.get('consistent_layouts'), d.get('best'))
for x in f:
    print(x['payer'], ','.join('+'.join(str(o) for o in c) for c in x['rels']), x['psi'])" $OUT/census.log 2>/dev/null)
say "census: $(echo "$LAY" | head -n 1 | cut -c1-200)"
[ "${CENSUS_ONLY:-0}" = 1 ] && { echo "census $(echo "$LAY" | head -n 1 | cut -c1-60)" > $OUT/DONE; exit 0; }
n=0
echo "$LAY" | tail -n +2 > $OUT/layouts.txt
while read -r P R S; do
  [ -n "$P" ] || continue
  n=$((n + 1))
  nice -n 10 taskset -c $CORE $PY -I r86_chem_g35_solve.py $OUT/k${K}_g4.npz --payer $P --rel $R --psi $S --nplans 400 \
     --out $OUT/k${K}_plan35_$n.json > $OUT/solve_$n.log 2>&1
  say "solve $n (payer $P rel $R psi $S): $(grep SUMMARY $OUT/solve_$n.log | cut -c1-200)"
  [ -s $OUT/k${K}_plan35_${n}_plans.json ] || continue
  R86_SOL=$OUT/k${K}_plan35_$n.json R86_DUMPMETA=$OUT/k${K}_g4.json C4_DEAD_OUT=EMIT nice -n 10 taskset -c $CORE \
     $PY -I r86_chem_g35_hook.py $TB $ARGS --inj 3=$OUT/k${K}_plan35_${n}_plans.json --emit $OUT/e$n > $OUT/hook_$n.log 2>&1
  say "loader $n: $(grep -h -o 'CLOSED [0-9]* products, emit_c4 toffoli [0-9]* ([a-z]*)' $OUT/hook_$n.log | tail -n 1)"
  if [ -f $OUT/e$n/angle_k_$K.kmx ] && grep -q "emit_c4 toffoli" $OUT/hook_$n.log; then
    cp $OUT/e$n/angle_k_$K.kmx $OUT/angle_k_$K.kmx; [ -f $OUT/e$n/table_dead.npz ] && cp $OUT/e$n/table_dead.npz $OUT/table_dead.npz
    say "DONE $(grep -h -o 'emit_c4 toffoli [0-9]*' $OUT/hook_$n.log | tail -n 1)"; echo ok > $OUT/DONE; exit 0
  fi
done < $OUT/layouts.txt
say "no loader"; echo fail > $OUT/DONE
