#!/usr/bin/env bash
# dqi_int_f901.sh -- B5 (17th holder, Nov 16; runs ON box 6, detached): finish a 901 file made by dqi_int_b901.sh:
# P1's dwx on the re-assigned file (its links pass caps at the file's own Q), quick; INV17's list A through a site map
# made from the marker build (as dqi_int_r48fin.sh 'make'); build2 -> file A; quick, official on both salts, paired
# against best_q24k3a on both official salts.
#   bash dqi_int_f901.sh <b901 tag> <marker build name> <cores> [Q=901]
# Lines -> ~/qcc/ws/dqi_int/f901_<tag>.txt; file A -> ~/qcc/ws/dqi_int/p901/<tag>/A/dqi.kmb
set -u
TG=$1; MBN=$2; C=$3; Q=${4:-901}
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PW=~/qcc/.venv/bin/python; PX=~/qcc/.venv_work/bin/python
O=$D/p901/$TG; G=$D/f901_$TG.txt; LINE=$D/r46/r46k3A/dqi.kmb; X=~/qcc/ws/r46_dqi_ext/prof1
R=$O/r$Q/dqi.kmb; F=$O/dwx/final
avail() { awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo; }
pickp() { a=$(avail); if [ $a -gt 40 ]; then echo 2; else echo 1; fi; }
pickv() { a=$(avail); if [ $a -gt 40 ]; then echo 4; elif [ $a -gt 25 ]; then echo 2; else echo 1; fi; }
echo "f901 $TG: input $(sha256sum $R | cut -c1-16) start $(date -u +%H:%M:%SZ)" >> $G
cd ~/Quantum && nice -n 10 bash tools/sh/${B5_DWX:-dqi_pp_dwx.sh} $R $O/dwx $C > $O/dwx.log 2>&1   # B5_DWX=dqi_int_dwxz.sh: R55's exact zeros to P1
echo "dwx: $(grep -a -o 'links admitted [0-9]*' $O/dwx/links/report.txt 2>/dev/null | tail -1) sha $(sha256sum $F/dqi.kmb | cut -c1-16) $(date -u +%H:%M:%SZ)" >> $G
cd $T
nice -n 10 taskset -c $C $PW qcc_verify.py $F --tasks dqi --quick --workers $(pickv) --json $F/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/dwx quick: /" >> $G
A=$O/A; mkdir -p $A
if [ "${B5_NOINV:-0}" = 1 ]; then      # INV17 already applied before the re-assignment (dqi_int_b901.sh B5_INV17=1)
  cp $F/dqi.kmb $A/dqi.kmb; echo "file A = dwx final (INV17 applied before re-assignment)" >> $G
else
MB=$D/r46/$MBN; M=$D/prof_$MBN; mkdir -p $M $D/maps
[ -s $M/segments_real.json ] || nice -n 10 taskset -c $C $PX -B r46_dqi_ext_segreal.py $MB/dqi.kmb $MB/segments.json $M/segments_real.json > $M/segreal.log 2>&1
S=""; for f in stat_s5 stat_s6 stat_s7 stat_s8 s3u_0 s3u_1 s3u_2 s3u_3; do [ -f $X/$f.npz ] && S="$S $X/$f.npz"; done
nice -n 10 taskset -c $C $PX -B dqi_int_inv17map2.py $X/dqi.kmb $X/segments_real.json $MB/dqi.kmb $M/segments_real.json $F/dqi.kmb $D/maps/f901_$TG $S 2>&1 | cut -c1-220 >> $G
nice -n 10 taskset -c $C $PX -B r44_xfer_dqidrop.py build2 $F/dqi.kmb $D/maps/f901_${TG}_A.npz $A/dqi.kmb > $A/build2.log 2>&1
echo "build2 A: $(tail -1 $A/build2.log | cut -c1-160)" >> $G
fi
echo "file A $(sha256sum $A/dqi.kmb | cut -c1-16) $(date -u +%H:%M:%SZ)" >> $G
nice -n 10 taskset -c $C $PW qcc_verify.py $A --tasks dqi --quick --workers $(pickv) --json $A/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/A quick: /" >> $G
for salt in "" alt-salt-1; do
  tag=${salt:-salt0}
  nice -n 10 taskset -c $C $PW qcc_verify.py $A --tasks dqi --salt "$salt" --workers $(pickv) --json $A/verify_$tag.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/A official $tag: /" >> $G
done
nice -n 10 taskset -c $C $PW -B dqi_int_rate.py $A/dqi.kmb $LINE salt0,alt-salt-1 $(pickp) $A/paired.tsv > $A/paired.log 2>&1
cut -f1-8,10 $A/paired.tsv | sed "s/^/paired A vs best_q24k3a: /" >> $G
echo "DONE $(date -u +%H:%M:%SZ)" >> $G
