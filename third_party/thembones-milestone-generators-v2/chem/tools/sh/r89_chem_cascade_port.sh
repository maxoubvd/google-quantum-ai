#!/bin/bash
# r89_chem_cascade_port.sh (seat R89-CHEM-CASCADE): ONE COMMAND for a 23-output loader K (k14 first: batch 14 has no other
# reader) on any table: C4's pair pipeline (chem_c4_dump / wb / ns / pre, read-only), my SUPPLY-AWARE pair choice
# (r89_chem_cascade_pair.py: group 4 = the four dead-richest outputs of the rest), a scan of group 4's killed column
# for a lossless pool compensation with 6 relations consistent (C4's chem_c4_dead.py census), then my form (full second
# column kill + 2 random relations + C4's cubic block) and C2's annealed group 3. Target 142 products (216 Toffolis).
#   bash r89_chem_cascade_port.sh K TABLE OUTDIR "core core ..." [npairs] [kill columns to scan]
# Result: OUTDIR/best/{angle_k_K.kmx,table_dead.npz} and a line 'SUCCESS ...' in OUTDIR/port.log.
K=$1; TB=$2; O=$3; CA=($4); NP=${5:-4}; XS=${6:-"13 12 11 10 14 9 8 7 6 5 4 2 1 0"}
N=${#CA[@]}; core() { echo ${CA[$(( $1 % N ))]}; }
Q=<HOME>/Quantum; [ -d $Q ] || Q=$HOME/Quantum         # this PC, or a box
cd $Q/tools || exit 1
PY=${PY:-$HOME/qcc/.venv/bin/python}; [ -x $PY ] || PY=$HOME/qcc/.venv_work/bin/python
export PY CHEM_BROT=7 CHEM_HB=77 C4_DROP=1
mkdir -p $O/best; log() { echo "$(date -u +%H:%M:%SZ) $*" >> $O/port.log; }
log "start K=$K table $TB cores ${CA[*]}"
# PAIRSRC="dir dir ...": reuse pair folders (args.txt with the order already supply-aware, g0.json, g1.json) instead
# of new scans; C4's pre-check below re-validates every one on TABLE. ONLY for the same target words: measured 04:46Z,
# k14's pairs from other batch-13 words fail the pre-check ('row outside the plan's span') on 6 of 6
if [ -n "${PAIRSRC:-}" ]; then
  mkdir -p $O/pairs; n=0
  for d in $PAIRSRC; do
    q=$O/pairs/p$n; mkdir -p $q; cp $d/g0.json $d/g1.json $q/; a=$d/args.txt; [ -f $d/args_r89.txt ] && a=$d/args_r89.txt
    cp $a $q/args.txt; n=$((n + 1))
  done
  log "pairs reused: $n from PAIRSRC"
else
[ -s $O/d.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_dump.py $TB --k $K --plan 6,6,7:1:0:0,4:1:0:1 --stop 149 \
   --seeds 0-15 --out $O/d.json > $O/dump.log 2>&1
[ -s $O/d_wb.json ] || $PY -I chem_c4_wb.py $TB --k $K --dump $O/d.json --out $O/d_wb.json > $O/wb.log 2>&1
[ -s $O/d_wb.json ] || { log "no dump: stop"; exit 0; }
# SCANROUNDS (default 1): rounds of one scan a core (2 cores x 1 round read 98 finds on k14 and no pair with all four
# slot-3 wires; 6 cores x 1 round read 358 and four such pairs)
for r in $(seq 0 $((${SCANROUNDS:-1} - 1))); do
  for i in $(seq 0 $((N - 1))); do
    s=$((r * N + i))
    [ -s $O/ns$s.json ] || nice -n 10 taskset -c $(core $i) $PY -I chem_c4_ns.py $O/d_wb.json --pool 0-22 --samples ${SAMPLES:-150} \
       --seed $((1000 * K + s + 1)) --top 6 --out $O/ns$s.json > $O/ns$s.log 2>&1 &
  done
  wait
done
log "scans: $(cat $O/ns*.log | grep -c FOUND) finds"
$PY -I r89_chem_cascade_cells.py $TB --k $K --json $O/cells.json > $O/cells.log 2>&1
$PY -I r89_chem_cascade_pair.py $O/ns*.json --cells $O/cells.json --nb 23 --max $NP --outdir $O/pairs > $O/pair.log 2>&1
log "pairs: $(head -n $NP $O/pair.log | cut -c1-90 | tr '\n' ';')"
fi
for p in $(ls -d $O/pairs/p* 2>/dev/null); do
  [ -s $p/g0_best.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_pre.py $TB --k $K --wirebasis --sizes 6,6,7,4 \
     --pair $p > $p/pre.log 2>&1
done
for p in $(grep -l WROTE $O/pairs/p*/pre.log 2>/dev/null | xargs -r -n1 dirname); do
  read -r ORD G3 <<< "$($PY -c "a=open('$p/args.txt').read().split(); print(a[a.index('--order')+1], a[a.index('--kills')+1].split(',')[0])")"
  X=""
  for x in $XS; do
    [ "$x" = "$G3" ] && continue
    echo "--order $ORD --kills $G3,$x" > $p/args_$x.txt; rm -f $p/g4_$x.npz
    ( export C4_DROPDUMP=$p/g4_$x.npz
      exec nice -n 10 taskset -c $(core 0) $PY -I chem_c4_build.py $TB --k $K --wirebasis --plan 6,6,7:1:0:0,4:1:0:1 \
        $(cat $p/args_$x.txt) --inj 0=$p/g0_best.json --inj 1=$p/g1_best.json --modes=-,-,-,- --tries 3000 > $p/kd_$x.log 2>&1 ) &
    bp=$!
    for i in $(seq 1 90); do [ -s $p/g4_$x.npz ] && break; kill -0 $bp 2>/dev/null || break; sleep 2; done
    kill $bp 2>/dev/null; wait $bp 2>/dev/null
    [ -s $p/g4_$x.npz ] || continue
    C4_DEAD_MS=6 nice -n 10 taskset -c $(core 0) $PY -I chem_c4_dead.py $TB --k $K --dump $p/g4_$x.npz --draws 6 > $p/kc_$x.log 2>&1
    ok=$(sed -n 's/.*m = 6 relations: consistent on \([0-9]*\) of.*/\1/p' $p/kc_$x.log)
    log "$(basename $p) X=$x: $(grep -o 'unknowns by output \[[0-9, ]*\]' $p/kc_$x.log | tail -n 1), m=6 consistent ${ok:-?}/6"
    [ -n "$ok" ] && [ "$ok" -ge 4 ] && { X=$x; break; }
  done
  [ -n "$X" ] || { log "$(basename $p): no kill column with 6 relations"; continue; }
  B=$O/$(basename $p)_x$X
  R89_G=2 bash $Q/tools/sh/r89_chem_cascade_go.sh $K $TB $p $p/args_$X.txt $B $(core 0)
  log "$(basename $p) b1: $(grep -v '^  attempt' $B/build.log | grep -E 'CLOSED|FAILED' | head -n 1 | cut -c1-60) | $(grep -o 'b=4 a=14 kills=(1,0) dimW=[0-9]* amb=56 bonus=12 -> [0-9]*' $B/build.log | tail -n 1)"
  grep -q "dimW=44 amb=56 bonus=12 -> 22" $B/build.log || continue
  bash $Q/tools/sh/r89_chem_cascade_an.sh $K $TB $p $p/args_$X.txt $B "${CA[*]}" ${ASECS:-300}
  b2=$(grep -l "CLOSED 142 products" $B/b2_*/build.log 2>/dev/null | head -n 1)
  if [ -n "$b2" ]; then
    d=$(dirname $b2); cp $d/angle_k_$K.kmx $d/table_dead.npz $O/best/
    log "SUCCESS $d: $(grep -v '^  attempt' $b2 | grep CLOSED | head -n 1 | cut -c1-120)"; exit 0
  fi
  log "$(basename $p): no 142 after the anneal"
done
log "end without 142"
# FALLBACK=1: builder C4's own dead23 (tools/sh/chem_c4_loader2.sh, read-only) in OUTDIR/fb, copied into OUTDIR/best,
# so that OUTDIR = <chain WD>/L<K> is a drop-in for chem_c4_chain.sh's 'K:dead23' step either way
if [ "${FALLBACK:-0}" = 1 ]; then
  ( export PY MAXP=${MAXP:-48} SAMPLES=${SAMPLES:-150} CL="${CA[*]}" C4_T1=144 C4_T2=143 C4_DEAD=4 C4_CUBIC=3 \
      C4_CUBIC_TRIES=1000 C4_DEAD_OUT=EMIT; bash $Q/tools/sh/chem_c4_loader2.sh $K 0 0 $TB $O/fb ) > /dev/null 2>&1
  if [ -s $O/fb/L$K/best/angle_k_$K.kmx ]; then
    cp $O/fb/L$K/best/angle_k_$K.kmx $O/best/; [ -f $O/fb/L$K/best/table_dead.npz ] && cp $O/fb/L$K/best/table_dead.npz $O/best/
    log "FALLBACK dead23: $(grep -h -E 'SUCCESS|kept' $O/fb/L$K/run.log | tail -n 1 | cut -c1-120)"
  else
    log "FALLBACK dead23: no circuit"
  fi
fi
