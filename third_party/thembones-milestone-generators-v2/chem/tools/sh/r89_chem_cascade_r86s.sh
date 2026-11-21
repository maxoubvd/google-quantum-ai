#!/bin/bash
# r89_chem_cascade_r86s.sh (seat R89-CHEM-CASCADE): design seat R86's group-4 form (+ C4's group-3 row drop) on a
# 24-output wired loader with SUPPLY-AWARE pairs: builder C4's pair pipeline (chem_c4_dump / wb / ns / pre, read-only),
# my selector (r89_chem_cascade_pair.py: the richest output of the rest kept in group 3 as the drop's payer, the next six
# richest = group 4), then my strict copy of C4 s chem_c4_r86go.sh (census on the OWN batch, then the loader).
#   bash r89_chem_cascade_r86s.sh K TABLE OUTDIR "core core" [npairs]      (LOCAL WSL; r86go.sh runs here)
K=$1; TB=$2; O=$3; CA=($4); NP=${5:-8}; N=${#CA[@]}; core() { echo ${CA[$(( $1 % N ))]}; }
Q=<HOME>/Quantum; [ -d $Q ] || Q=$HOME/Quantum
PY=${PY:-$HOME/qcc/.venv/bin/python}; [ -x $PY ] || PY=$HOME/qcc/.venv_work/bin/python
cd $Q/tools || exit 1; export CHEM_BROT=7 CHEM_HB=77 C4_DROP=1
mkdir -p $O; log() { echo "$(date -u +%H:%M:%SZ) $*" >> $O/r86s.log; }
log "start K=$K table $TB"
[ -s $O/d.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_dump.py $TB --k $K --plan 6,6,6:1:0:0,6:1:2:1 --stop 153 \
   --seeds 0-15 --out $O/d.json > $O/dump.log 2>&1
[ -s $O/d_wb.json ] || $PY -I chem_c4_wb.py $TB --k $K --dump $O/d.json --out $O/d_wb.json > $O/wb.log 2>&1
for r in 0 1 2; do
  for i in $(seq 0 $((N - 1))); do
    s=$((r * N + i)); [ -s $O/ns$s.json ] || nice -n 10 taskset -c $(core $i) $PY -I chem_c4_ns.py $O/d_wb.json --pool 0-23 \
       --samples ${SAMPLES:-150} --seed $((1000 * K + s + 1)) --top 6 --out $O/ns$s.json > $O/ns$s.log 2>&1 &
  done
  wait
done
log "scans: $(cat $O/ns*.log | grep -c FOUND) finds"
[ -s $O/cells.json ] || $PY -I r89_chem_cascade_cells.py $TB --k $K --json $O/cells.json > $O/cells.log 2>&1
$PY -I r89_chem_cascade_pair.py $O/ns*.json --cells $O/cells.json --nb 24 --g4 6 --g3first 1 --x 14 --max $NP \
   --outdir $O/pairs > $O/pair.log 2>&1
log "pairs: $(tail -n 1 $O/pair.log)"
for p in $(ls -d $O/pairs/p* 2>/dev/null); do
  [ -s $p/g0_best.json ] || nice -n 10 taskset -c $(core 0) $PY -I chem_c4_pre.py $TB --k $K --wirebasis --sizes 6,6,6,6 \
     --pair $p > $p/pre.log 2>&1
done
i=0
for p in $(grep -l WROTE $O/pairs/p*/pre.log 2>/dev/null | xargs -r -n1 dirname); do
  ( R86_BATCHES=own CENSUS_ONLY=1 bash $Q/tools/sh/r89_chem_cascade_r86go.sh $K $TB $p $O/s/$(basename $p) $(core $i) > /dev/null 2>&1 ) &
  i=$((i + 1)); [ $((i % N)) = 0 ] && wait
done
wait
for d in $O/s/p*; do log "$(basename $d) $(head -n 1 $p/../$(basename $d)/r89.json 2>/dev/null) $(grep census $d/r86go.log | cut -c1-90)"; done
best=$(for d in $O/s/p*; do n=$(sed -n 's/.*census: \([0-9]*\) \([0-9]*\) .*/\2/p' $d/r86go.log 2>/dev/null | head -n 1); \
  [ -n "$n" ] && [ "$n" -gt 0 ] && echo "$n $(basename $d)"; done | sort -n -r | head -n $N | awk '{print $2}')
log "deliverable pairs (own batch): ${best:-none}"
j=0
for p in $best; do
  ( R86_BATCHES=own bash $Q/tools/sh/r89_chem_cascade_r86go.sh $K $TB $O/pairs/$p $O/b/$p $(core $j) > /dev/null 2>&1 ) &
  j=$((j + 1))
done
wait
for p in $best; do
  log "$p: $(tail -n 1 $O/b/$p/r86go.log 2>/dev/null | cut -c1-120)"
  grep -q "DONE emit_c4 toffoli 221" $O/b/$p/r86go.log 2>/dev/null && { log "SUCCESS $O/b/$p (221)"; exit 0; }
done
log "end without 221"
