#!/bin/bash
# COPY by builder C3 of builder C4's tools/sh/chem_c4_loader.sh (read-only original): the same steps, the work folder
# under ~/qcc/ws/chem_c3 (WD=<folder>, default ~/qcc/ws/chem_c3/c4L) and the table argument required.
# chem_c4_loader.sh: builder C4. "TWO POOL-FREE 38s" for one wired loader K, on cores C0 .. C0+N-1 of this box:
#   groups 1 and 2 = two disjoint 6-output subsets at 38 by exact hyperplane plans, every product non-silent on group
#   3's kill coordinate t (wire basis: chem_c4_wb.py), groups 3 and 4 by the production engine, then builder C2's
#   annealer on group 3 (42). Target 154 products = 228 Toffolis (today's plain wired loaders: 229).
#   bash chem_c4_loader.sh <K> <C0> <N> [table]        (PY defaults to python3; folder ~/qcc/ws/chem_c4/L<K>)
K=$1; C0=$2; N=$3; TABLE=$4
# CL="a b c ..." (sixth holder): an explicit core list (any cores); then C0 / N are ignored
if [ -n "${CL:-}" ]; then CA=($CL); else CA=($(seq $C0 $((C0 + N - 1)))); fi; N=${#CA[@]}; C0=0
core() { echo ${CA[$(( $1 % N ))]}; }
PY=${PY:-python3}; W=${WD:-$HOME/qcc/ws/chem_c3/c4L}/L$K; mkdir -p $W/best; cd ~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77}
say() { echo "$(date -u +%H:%M:%S) $*" >> $W/run.log; }
say "start K=$K cores $C0+$N"
[ -f $W/d.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_dump.py $TABLE --k $K --seeds 0-15 --out $W/d.json > $W/dump.log 2>&1
[ -f $W/d_wb.json ] || $PY -I chem_c4_wb.py $TABLE --k $K --dump $W/d.json --out $W/d_wb.json >> $W/run.log 2>&1
say "dump: $(tail -n 1 $W/dump.log)"
for i in $(seq 1 $((N-1))); do
  nice -n 10 taskset -c $(core $i) $PY -I chem_c4_ns.py $W/d_wb.json --pool 0-24 --samples ${SAMPLES:-90} \
     --seed $((1000*K+i)) --top 6 --out $W/ns$i.json > $W/ns$i.log 2>&1 &
done
wait
say "scans: $(cat $W/ns*.log | grep -c FOUND) finds"
$PY -I chem_c4_pair.py $W/ns*.json --outdir $W/pairs --max ${MAXP:-24} --t2 14,13 > $W/pair.log 2>&1
say "pairs: $(tail -n 1 $W/pair.log | cut -c1-120)"
# pre-check every pair (parallel), then build the ones with a closing, full-rank plan pair
ls -d $W/pairs/p* | while read p; do echo $p; done > $W/plist
i=0
while read p; do
  ( nice -n 10 taskset -c $(core $i) $PY -I chem_c4_pre.py $TABLE --k $K --wirebasis --pair $p > $p/pre.log 2>&1 ) &
  i=$((i+1)); if [ $((i % N)) -eq 0 ]; then wait; fi
done < $W/plist
wait
say "pre-checks: $(grep -l WROTE $W/pairs/p*/pre.log | wc -l) of $(wc -l < $W/plist) pairs ready"
ok=""
i=0
for p in $(grep -l WROTE $W/pairs/p*/pre.log | xargs -n1 dirname); do
  ( nice -n 10 taskset -c $(core $i) $PY -I chem_c4_build.py $TABLE --k $K --wirebasis $(cat $p/args.txt) \
      --inj 0=$p/g0_best.json --inj 1=$p/g1_best.json --modes=-,-,-,- --tries 3000 --emit $p/b155 \
      --dumpg $p/dump155.json > $p/b155.log 2>&1 ) &
  i=$((i+1)); if [ $((i % N)) -eq 0 ]; then wait; fi
done
wait
for p in $(grep -l "CLOSED 155 products" $W/pairs/p*/b155.log 2>/dev/null | xargs -n1 dirname); do ok="$ok $p"; done
say "builds at 155 products: $(echo $ok | wc -w)"
for p in $ok; do
  # builder C2's annealer on group 3 (index 2) for 42, N runs of ${ASECS:-300} s, then explicit plans and the build
  for r in $(seq 0 $((N-1))); do
    nice -n 10 taskset -c $(core $r) $PY -I chem_ld2_e19grp.py anneal $p/dump155.json --g 2 --P 42 --secs ${ASECS:-300} \
       --every 20000 --price 40 --rseed $r --save $p/an_$r.json > $p/an_$r.log 2>&1 &
  done
  wait
  for r in $(seq 0 $((N-1))); do
    [ -s $p/an_$r.json ] || continue
    $PY -I chem_ld2_e19grp.py plan $p/dump155.json $p/an_$r.json --g 2 --tries 20000 --out $p/pl42_$r.json > $p/pl42_$r.log 2>&1
    [ -s $p/pl42_$r.json ] || continue
    nice -n 10 taskset -c $(core 0) $PY -I chem_c4_build.py $TABLE --k $K --wirebasis $(cat $p/args.txt) \
       --inj 0=$p/g0_best.json --inj 1=$p/g1_best.json --inj 2=$p/pl42_$r.json --modes=-,-,-,- --tries 3000 \
       --emit $p/b154 > $p/b154.log 2>&1
    if grep -q "CLOSED 154 products" $p/b154.log; then
      cp $p/b154/angle_k_$K.kmx $W/best/angle_k_$K.kmx
      say "SUCCESS $p: $(tail -n 2 $p/b154.log | head -n 1 | cut -c1-140)"
      exit 0
    fi
  done
  say "no 154 from $p"
done
# fall back: keep a 155-product (229) circuit as the tie
for p in $ok; do cp $p/b155/angle_k_$K.kmx $W/best/angle_k_${K}_229.kmx; say "kept 229 from $p"; break; done
say "end without 154"
