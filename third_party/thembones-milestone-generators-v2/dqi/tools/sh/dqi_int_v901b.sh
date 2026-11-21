#!/usr/bin/env bash
# dqi_int_v901b.sh -- B5 (17th holder, Nov 16; runs ON box 6, detached): a 901 variant from its EXISTING marker build
# r46/<name>p: census + flag checks, then the spec port; while runs stay uncovered, the cover flags whose re-test is not
# exact are excluded (B5_IDLE_EXCLUDE) and the census is made again (up to 4 rounds); then dqi_int_b901.sh with INV17
# before the re-assignment, P1's dwx at 901, quick.
#   bash dqi_int_v901b.sh <name> <c1,c2> <c1-c2>          lines -> ~/qcc/ws/dqi_int/v901b_<name>.txt
set -u
N=${1:?name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PX=~/qcc/.venv_work/bin/python; PW=~/qcc/.venv/bin/python
G=$D/v901b_$N.txt; MB=${N}p; O=$D/p901/$MB; K=$D/r46/$MB/dqi.kmb
echo "v901b $N start $(date -u +%H:%M:%SZ) marker $(sha256sum $K | cut -c1-16)" >> $G
EX=""
for it in 1 2 3 4; do
  B5_IDLE_EXCLUDE=$EX bash $T/sh/dqi_int_c901.sh $K $CC 902
  cd $T
  L=$(grep -a '^unconvert' $O/idle902.txt | grep -a ' - \[conv > rank_update' | head -1)
  R=$(echo "$L" | awk '{print $2}'); RA=${R%%..*}; RB=${R##*..}; HM=""
  if [ -n "$R" ]; then
    [ -s $O/linpark_X.txt ] || nice -n 10 taskset -c $CC $PX -B dqi_int_linpark.py $K $O/zero.npz $RA:$RB > $O/linpark_X.txt 2>&1
    W=$(grep -a 'support 1:' $O/linpark_X.txt | head -1 | awk '{print $3}')
    [ -n "$W" ] && HM="--hmr $($PX -B dqi_int_lasttouch.py $K $W $((RA-20000)) $RA | tail -1 | awk '{print $2}') $W"
  fi
  $PX -B dqi_int_spec901.py $D/p901/idle_l9r4dp_902.txt $D/p901/spec_k1.txt $O/idle902.txt $O/flagchk902.txt $O/spec.txt $HM > $O/port.log 2>&1 || { echo "port FAILED" >> $G; exit 1; }
  echo "round $it exclude '$EX': $(grep -a 'coverage' $O/port.log) | $(grep -c '^park' $O/spec.txt) parks" >> $G
  grep -a -q UNCOVERED $O/port.log || break
  BAD=$(grep -a '^f ' $O/flagchk902.txt | grep -v 'RE-TEST AT THE CLOSE EXACT ON SAMPLES' | awk '{print $2}' | sort -u | tr '\n' ',' | sed 's/,$//')
  [ -z "$BAD" ] && { echo "uncovered runs and no non-exact cover flag to exclude" >> $G; exit 1; }
  EX=${EX:+$EX,}$BAD
done
grep -a -q UNCOVERED $O/port.log && { echo "still uncovered after 4 rounds" >> $G; exit 1; }
cp $O/spec.txt $D/p901/spec_${N}.txt
rm -rf $D/p901/w$N
B5_INV17=1 bash $T/sh/dqi_int_b901.sh $MB $D/p901/spec_${N}.txt $CC w$N 901
tail -3 $D/p901/w$N/log.txt >> $G
R9=$D/p901/w$N/r901/dqi.kmb
[ -s $R9 ] || { echo "no 901 file" >> $G; exit 1; }
cd ~/Quantum && nice -n 10 bash tools/sh/dqi_pp_dwx.sh $R9 $D/p901/w$N/dwx $PC > $D/p901/w$N/dwx.log 2>&1
F=$D/p901/w$N/dwx/final
echo "dwx: $(grep -a -o 'links admitted [0-9]*' $D/p901/w$N/dwx/links/report.txt | tail -1)" >> $G
cd $T && nice -n 10 taskset -c $CC $PW qcc_verify.py $F --tasks dqi --quick --workers 2 --json $F/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/dwx quick (INV17 in): /" >> $G
rm -f $D/p901/w$N/k.kmb $D/p901/w$N/ks.kmb $D/p901/w$N/ksi.kmb; rm -rf $D/p901/w$N/dwx/dw1 $D/p901/w$N/dwx/links/*.npz
echo "DONE $(date -u +%H:%M:%SZ) sha $(sha256sum $F/dqi.kmb | cut -c1-16)" >> $G
