#!/bin/bash
# chem_c4_loader2.sh: builder C4 (holder 2). "TWO POOL-FREE 38s" for one wired loader K of ANY generic count:
#   25 outputs (plain wired): plan 6 | 6 | 7:1 | 6:1:3 -> 155, then C2's 42-annealer on group 3 -> 154 = 228;
#   24 outputs (one bit-1 relation): plan 6 | 6 | 6:1 | 6:1:2 (C2's PLANS_BY_N[24], the tie engine reads 152 = 226):
#      two 38s + 38 + 37 -> 151 = 225;
#   23 outputs (two relations): plan 6 | 6 | 7:1 | 4:1:2 (C3's relrun job: 39|38|43|28, 147 with the annealed 7-group):
#      two 38s + 43 + 28 -> 147, then the 42-annealer on group 3 -> 146 = 220.
# Same steps and tools as chem_c4_loader.sh (that file stays as it is for C3's copy).
#   bash chem_c4_loader2.sh <K> <C0> <N> <table> [workdir]       (CL="a b c" for an explicit core list)
K=$1; C0=$2; N=$3; TABLE=$4; WD=${5:-$HOME/qcc/ws/chem_c4/L2}
if [ -n "${CL:-}" ]; then CA=($CL); else CA=($(seq $C0 $((C0 + N - 1)))); fi; N=${#CA[@]}
core() { echo ${CA[$(( $1 % N ))]}; }
PY=${PY:-python3}; W=$WD/L$K; mkdir -p $W/best; cd ~/Quantum/tools || exit 1
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77}
say() { echo "$(date -u +%H:%M:%S) $*" >> $W/run.log; }
NG=$($PY -I - "$TABLE" $K "$PWD" << 'PYEOF'
import sys
sys.path.insert(0, sys.argv[3])
import numpy as np
import chem_loaders_core as core
import chem_ld2_e19grp as eg
core.J, core.BROT = 5, 7
w = np.asarray(np.load(sys.argv[1])["angle_words"], np.int64)
_, wp = eg.block_set(w, int(sys.argv[2]))
print(int(core.analyze(wp, int(sys.argv[2]), 4).X0.shape[1]))
PYEOF
)
case "$NG" in
  25) PLAN="6,6,7:1:0:0,6:1:3:1"; STOP=156; SIZES="6,6,7,6"; T1=155; T2=154; AN=1;;
  24) PLAN="6,6,6:1:0:0,6:1:2:1"; STOP=153; SIZES="6,6,6,6"; T1=151; T2=151; AN=0;;
  23) PLAN="6,6,7:1:0:0,4:1:2:1"; STOP=149; SIZES="6,6,7,4"; T1=147; T2=146; AN=1;;
  *) say "k$K: $NG generic outputs: no plan"; exit 0;;
esac
# holder 3 (12:45Z): with C4_DROP=1 the 23-output group 4 (4 outputs, 56 cells: bound by its cells) keeps its pool's free
# bits for two COLUMN relations (tools/chem_c4_drop.py) instead of two row kills that cannot lower its rank: 28 -> 27
if [ "${C4_DROP:-0}" = 1 ] && [ "$NG" = 23 ]; then PLAN="6,6,7:1:0:0,4:1:0:1"; T1=146; T2=145; fi
# holder 4: targets from the environment (dead toggles / cubic blocks reach below the plain targets)
[ -n "${C4_T1:-}" ] && T1=$C4_T1
[ -n "${C4_T2:-}" ] && T2=$C4_T2
say "start K=$K ($NG generic outputs, plan $PLAN, targets $T1 / $T2) cores ${CA[*]}"
[ -f $W/d.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_dump.py $TABLE --k $K --plan $PLAN --stop $STOP \
   --seeds 0-15 --out $W/d.json > $W/dump.log 2>&1
[ -f $W/d_wb.json ] || $PY -I chem_c4_wb.py $TABLE --k $K --dump $W/d.json --out $W/d_wb.json >> $W/run.log 2>&1
say "dump: $(tail -n 1 $W/dump.log | cut -c1-150)"
[ -s $W/d_wb.json ] || { say "no dump: stop"; exit 0; }
for i in $(seq 1 $((N - 1))); do
  nice -n 10 taskset -c $(core $i) $PY -I chem_c4_ns.py $W/d_wb.json --pool 0-$((NG - 1)) --samples ${SAMPLES:-90} \
     --seed $((1000 * K + i)) --top 6 --out $W/ns$i.json > $W/ns$i.log 2>&1 &
done
wait
say "scans: $(cat $W/ns*.log | grep -c FOUND) finds"
$PY -I chem_c4_pair.py $W/ns*.json --nb $NG --outdir $W/pairs --max ${MAXP:-24} --t2 14,13 > $W/pair.log 2>&1
say "pairs: $(tail -n 1 $W/pair.log | cut -c1-120)"
# holder 4 (19:18Z): C4_REST_FIRST="o,o,...": after the two pool-free groups, these outputs come first (group 3), e.g. the
# outputs whose dead cells pay group 3's row drop (the pair step orders the rest by index)
if [ -n "${C4_REST_FIRST:-}" ]; then
  for f in $W/pairs/p*/args.txt; do
    $PY -I - "$f" "$C4_REST_FIRST" << 'PYEOF'
import sys
f, first = sys.argv[1], [int(x) for x in sys.argv[2].split(",")]
a = open(f).read().split()
i = a.index("--order")
o = [int(x) for x in a[i + 1].split(",")]
S, rest = o[:12], o[12:]
rest = [x for x in first if x in rest] + [x for x in rest if x not in first]
a[i + 1] = ",".join(str(x) for x in S + rest)
open(f, "w").write(" ".join(a) + "\n")
PYEOF
  done
  say "rest reordered: $C4_REST_FIRST first"
fi
ls -d $W/pairs/p* 2>/dev/null > $W/plist
i=0
while read p; do
  ( nice -n 10 taskset -c $(core $i) $PY -I chem_c4_pre.py $TABLE --k $K --wirebasis --sizes $SIZES --pair $p \
      > $p/pre.log 2>&1 ) &
  i=$((i + 1)); if [ $((i % N)) -eq 0 ]; then wait; fi
done < $W/plist
wait
say "pre-checks: $(grep -l WROTE $W/pairs/p*/pre.log 2>/dev/null | wc -l) of $(wc -l < $W/plist) pairs ready"
i=0
for p in $(grep -l WROTE $W/pairs/p*/pre.log 2>/dev/null | xargs -r -n1 dirname); do
  ( nice -n 10 taskset -c $(core $i) $PY -I chem_c4_build.py $TABLE --k $K --wirebasis --plan $PLAN $(cat $p/args.txt) \
      --inj 0=$p/g0_best.json --inj 1=$p/g1_best.json --modes=-,-,-,- --tries 3000 --emit $p/b1 \
      --dumpg $p/dumpb1.json > $p/b1.log 2>&1 ) &
  i=$((i + 1)); if [ $((i % N)) -eq 0 ]; then wait; fi
done
wait
# every closed build, fewest products first (11:20Z: k9 of the P57x projection closed one pair at 150 under the target)
ok=$(for f in $W/pairs/p*/b1.log; do n=$(grep -h -o -E 'CLOSED [0-9]+ products' $f | awk '{print $2}'); \
     [ -n "$n" ] && [ "$n" -le $T1 ] && echo "$n $(dirname $f)"; done 2>/dev/null | sort -n | awk '{print $2}')
say "builds at $T1 products or fewer: $(echo $ok | wc -w); best seen: $(grep -h -o -E 'CLOSED [0-9]+ products' $W/pairs/p*/b1.log 2>/dev/null | sort -k2 -n | head -n 1)"
for p in $ok; do
  n=$(grep -h -o -E 'CLOSED [0-9]+ products' $p/b1.log | awk '{print $2}')
  if [ $AN = 0 ] || [ "$n" -le $T2 ]; then
    cp $p/b1/angle_k_$K.kmx $W/best/angle_k_$K.kmx; [ -f $p/b1/table_dead.npz ] && cp $p/b1/table_dead.npz $W/best/table_dead.npz
    say "SUCCESS $p: $(grep CLOSED $p/b1.log | cut -c1-150)"; exit 0
  fi
done
[ $AN = 0 ] && { say "end without $T1"; exit 0; }
for p in $ok; do
  for r in $(seq 0 $((N - 1))); do
    nice -n 10 taskset -c $(core $r) $PY -I chem_ld2_e19grp.py anneal $p/dumpb1.json --g 2 --P 42 --secs ${ASECS:-300} \
       --every 20000 --price 40 --rseed $r --save $p/an_$r.json > $p/an_$r.log 2>&1 &
  done
  wait
  for r in $(seq 0 $((N - 1))); do
    [ -s $p/an_$r.json ] || continue
    $PY -I chem_ld2_e19grp.py plan $p/dumpb1.json $p/an_$r.json --g 2 --tries 20000 --out $p/pl42_$r.json > $p/pl42_$r.log 2>&1
    [ -s $p/pl42_$r.json ] || continue
    nice -n 10 taskset -c $(core 0) $PY -I chem_c4_build.py $TABLE --k $K --wirebasis --plan $PLAN $(cat $p/args.txt) \
       --inj 0=$p/g0_best.json --inj 1=$p/g1_best.json --inj 2=$p/pl42_$r.json --modes=-,-,-,- --tries 3000 \
       --emit $p/b2 > $p/b2.log 2>&1
    if grep -q "CLOSED $T2 products" $p/b2.log; then
      cp $p/b2/angle_k_$K.kmx $W/best/angle_k_$K.kmx; [ -f $p/b2/table_dead.npz ] && cp $p/b2/table_dead.npz $W/best/table_dead.npz
      say "SUCCESS $p: $(grep CLOSED $p/b2.log | cut -c1-150)"
      exit 0
    fi
  done
  say "no $T2 from $p"
done
for p in $ok; do cp $p/b1/angle_k_$K.kmx $W/best/angle_k_${K}_t1.kmx; say "kept $T1 from $p"; break; done
say "end without $T2"
