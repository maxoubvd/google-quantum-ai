#!/usr/bin/env bash
# dqi_int_r901repro.sh -- B5 (17th holder, Nov 16; runs ON box 6): REPRODUCE a 901 file end to end into FRESH folders:
# marker build <fresh>p (dqi_int_r50c.sh, the tail module in) -> 902 census + flag checks (dqi_int_c901.sh) -> X's HMR
# (linpark + last touch) -> park spec ported from l9r4dp's p901/spec_k1.txt (dqi_int_spec901.py, coverage-checked; while
# runs stay uncovered, the cover flags whose re-test is not exact are excluded: B5_IDLE_EXCLUDE, as dqi_int_v901b.sh) ->
# dqi_int_b901.sh with the park pass dqi_int_park901c.py (B5_PARK_MOVEZ=1), INV17 before the re-assignment (B5_INV17=1),
# re-assignment pinning only the syndrome wires' ends (B5_REALLOC_PIN=ends) -> dqi_int_f901.sh B5_NOINV=1 (P1's dwx at
# 901, quick, official both salts, paired vs best_q24k3a) -> optional re-roll tag (B5_REROLL="<label> <expected sha16>").
# Refuses to start if any of its folders exists.
#   bash dqi_int_r901repro.sh <fresh> <c1,c2> <c1-c2> <variant> [VAR=V ...]     lines -> ~/qcc/ws/dqi_int/repro901_<fresh>.txt
set -u
N=${1:?fresh}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; VAR=${4:?variant}; shift 4
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PX=~/qcc/.venv_work/bin/python; PW=~/qcc/.venv/bin/python; G=$D/repro901_$N.txt; MB=${N}p
for x in $D/pinvar/$MB $D/r46/$MB $D/p901/$MB $D/p901/$N $D/f901_$N.txt $G; do [ -e $x ] && { echo "$x exists: give a fresh name"; exit 1; }; done
FX="B5_ROOT4_FORM=d B5_ROOT4_EVICT=1 R49_M_DIRS=both R49_M_RHOLD=1 R49_M_HFWD=1 R49_M_RSHARE=1"
echo "repro901 $N variant $VAR extra $* start $(date -u +%H:%M:%SZ)" >> $G
B5_FILES=dqi_int_tail901.py${B5_FILES_ADD:+,$B5_FILES_ADD} B5_PRE_EXT=dqi_int_tail901 nice -n 10 bash $T/sh/dqi_int_r50c.sh $MB $CC $PC ${VAR}p $FX "$@" >> $G 2>&1
until [ -s $D/r46/$MB/dqi.kmb ] || grep -a -q '^end_' $D/pinvar/$MB/RESULT.txt 2>/dev/null; do sleep 15; done
BQ=$D/pinvar/$MB/home/qcc/ws/dqi_int/r12/b_q8u
if [ ! -s $D/r46/$MB/dqi.kmb ] && [ -s $BQ/quick.log ]; then mkdir -p $D/r46/$MB; ln -f $BQ/dqi.kmb $D/r46/$MB/dqi.kmb; cp $BQ/segments.json $D/r46/$MB/; fi
K=$D/r46/$MB/dqi.kmb; O=$D/p901/$MB
echo "marker $(sha256sum $K | cut -c1-16) quick: $(grep -a cost $BQ/quick.log | tail -1)" >> $G
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
  echo "round $it exclude '$EX': $(grep -a 'coverage' $O/port.log)" >> $G
  grep -a -q UNCOVERED $O/port.log || break
  BAD=$(grep -a '^f ' $O/flagchk902.txt | grep -v 'RE-TEST AT THE CLOSE EXACT ON SAMPLES' | awk '{print $2}' | sort -u | tr '\n' ',' | sed 's/,$//')
  [ -z "$BAD" ] && { echo "uncovered runs and no non-exact cover flag to exclude" >> $G; exit 1; }
  EX=${EX:+$EX,}$BAD
done
echo "spec $(grep -v '^#' $O/spec.txt | sha256sum | cut -c1-16) ($(grep -c . $O/spec.txt) lines)" >> $G
B5_PARK_TOOL=dqi_int_park901c.py B5_PARK_MOVEZ=1 B5_REALLOC_PIN=ends B5_INV17=1 bash $T/sh/dqi_int_b901.sh $MB $O/spec.txt $CC $N 901
grep -a 'park:\|realloc:\|quick r901\|sha r901' $D/p901/$N/log.txt | cut -c1-200 >> $G
B5_NOINV=1 bash $T/sh/dqi_int_f901.sh $N $MB $CC 901
grep -a 'dwx:\|file A\|A quick\|official\|paired' $D/f901_$N.txt >> $G
if [ -n "${B5_REROLL:-}" ]; then
  set -- $B5_REROLL
  mkdir -p $D/p901/$N/rr && (cd $T && $PW dqi_int_pptag.py $D/p901/$N/A/dqi.kmb $D/p901/$N/rr/dqi.kmb $1 > /dev/null 2>&1)
  S=$(sha256sum $D/p901/$N/rr/dqi.kmb | cut -c1-16)
  echo "re-roll $1: $S expected $2: $([ "$S" = "$2" ] && echo MATCH || echo DIFFER)" >> $G
fi
echo "DONE $(date -u +%H:%M:%SZ)" >> $G
