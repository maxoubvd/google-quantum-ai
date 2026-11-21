#!/bin/bash
# r89_chem_cascade_r86go.sh (seat R89-CHEM-CASCADE): COPY of builder C4 tools/sh/chem_c4_r86go.sh with ONE change: the dump
# and the hook run through tools/r89_chem_cascade_strict.py (group 3 row drop STRICT: a miss retries the group), so the
# dump and the loader both reach group 4 with the drop paid. Original header:
# chem_c4_r86go.sh: builder C4 (holder 5, Nov 18 21:52Z). LOCAL WSL, one core. Design seat R86-CHEM-G35's group-4-at-35
# form for ONE 24-output wired loader K on TABLE with a pool-free pair PAIR (args.txt, g0_best.json, g1_best.json), the
# seat's four steps as written in its 'For C4' block (its tools read-only): dump -> census (first deliverable layout) ->
# solve (400 plans) -> the loader through C4's engine (hook) -> OUT/{angle_k_K.kmx, table_dead.npz}. Log OUT/r86go.log.
#   bash chem_c4_r86go.sh K TABLE PAIR OUT [CORE]
K=$1; TB=$2; PAIR=$3; OUT=$4; CORE=${5:-10}
mkdir -p $OUT
if [ -d <HOME>/Quantum/tools ]; then cd <HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
else cd ~/Quantum/tools || exit 1; PY=~/qcc/.venv_work/bin/python; fi
export C4_DEAD_BATCHES=${R86_BATCHES:-both} CHEM_BROT=7 CHEM_HB=77 C4_DROP=1 C4_DEAD_ROW=3:2:72 C4_DEAD_ROW_TRIALS=100
ARGS="--k $K --wirebasis --plan 6,6,6:1:0:0,6:1:0:1 $(cat $PAIR/args.txt) --inj 0=$PAIR/g0_best.json --inj 1=$PAIR/g1_best.json --modes=-,-,-,- --tries 3000"
say() { echo "$(date -u +%H:%M:%S) $*" >> $OUT/r86go.log; }
say "start k$K table $TB pair $PAIR"
R86_GI=3 R86_DUMP=$OUT/k${K}_g4.npz nice -n 10 taskset -c $CORE $PY -I r89_chem_cascade_strict.py r86_chem_g35_dump.py $TB $ARGS > $OUT/dump.log 2>&1
[ -s $OUT/k${K}_g4.npz ] || { say "no dump: $(tail -n 2 $OUT/dump.log | tr '\n' ' ' | cut -c1-200)"; exit 2; }
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
     $PY -I r89_chem_cascade_strict.py r86_chem_g35_hook.py $TB $ARGS --inj 3=$OUT/k${K}_plan35_${n}_plans.json --emit $OUT/e$n > $OUT/hook_$n.log 2>&1
  say "loader $n: $(grep -h -o 'CLOSED [0-9]* products, emit_c4 toffoli [0-9]* ([a-z]*)' $OUT/hook_$n.log | tail -n 1)"
  if [ -f $OUT/e$n/angle_k_$K.kmx ] && grep -q "emit_c4 toffoli" $OUT/hook_$n.log; then
    cp $OUT/e$n/angle_k_$K.kmx $OUT/angle_k_$K.kmx; cp $OUT/e$n/table_dead.npz $OUT/table_dead.npz
    say "DONE $(grep -h -o 'emit_c4 toffoli [0-9]*' $OUT/hook_$n.log | tail -n 1)"; echo ok > $OUT/DONE; exit 0
  fi
done < $OUT/layouts.txt
say "no loader"; echo fail > $OUT/DONE
