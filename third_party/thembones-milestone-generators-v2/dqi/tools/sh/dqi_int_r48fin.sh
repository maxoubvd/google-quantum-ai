#!/usr/bin/env bash
# dqi_int_r48fin.sh -- B5 (15th holder, Nov 16; runs ON box 6, detached): finish a dqi_int_r48c.sh build <name>:
# wait for its dwx file, then INV17's list A (the 28 proved sites) through a site map, build2, quick + official on both
# salts, and the paired test against best_q24k3a (r46k3A) on both official salts.  RET2's chain7.sh steps, my folders.
#   bash dqi_int_r48fin.sh <name> <map> <cores> [<marker build name>]
#     <map> = a ready map file (e.g. RET2's ~/qcc/ws/r48_dqi_ret2/maps/lf6_A.npz, valid only on the dwx sha it was
#             made for), or "make": map made here from the marker build <marker build name> (dqi_int_r48c.sh ...p)
# Lines -> ~/qcc/ws/dqi_int/fin48_<name>.txt; file A -> ~/qcc/ws/dqi_int/r46/<name>A/dqi.kmb
set -u
N=$1; MAP=$2; C=$3; MBN=${4:-}
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PW=~/qcc/.venv/bin/python; PX=~/qcc/.venv_work/bin/python
G=$D/fin48_$N.txt; LINE=$D/r46/r46k3A/dqi.kmb; X=~/qcc/ws/r46_dqi_ext/prof1
avail() { awk '/MemAvailable/ {print int($2/1048576)}' /proc/meminfo; }
pickp() { a=$(avail); if [ $a -gt 40 ]; then echo 2; else echo 1; fi; }
pickv() { a=$(avail); if [ $a -gt 40 ]; then echo 4; elif [ $a -gt 25 ]; then echo 2; else echo 1; fi; }
until grep -a -q "final" $D/rebuild_$N.txt 2>/dev/null; do sleep 20; done
echo "rebuild line: $(cat $D/rebuild_$N.txt)" >> $G
F=$D/dwx/$N/final
echo "raw $(sha256sum $D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u/dqi.kmb | cut -c1-16) dwx $(sha256sum $F/dqi.kmb | cut -c1-16) $(date -u +%H:%MZ)" >> $G
cd $T
nice -n 10 taskset -c $C $PW qcc_verify.py $F --tasks dqi --quick --workers $(pickv) --json $F/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/dwx quick: /" >> $G
if [ "$MAP" = make ]; then
  MB=$D/pinvar/$MBN/home/qcc/ws/dqi_int/r12/b_q8u; O=$D/prof_$MBN
  until [ -s $MB/quick.log ] || grep -a -q '^end_' $D/pinvar/$MBN/RESULT.txt 2>/dev/null; do sleep 10; done
  mkdir -p $O $D/maps
  nice -n 10 taskset -c $C $PX -B r46_dqi_ext_segreal.py $MB/dqi.kmb $MB/segments.json $O/segments_real.json > $O/segreal.log 2>&1
  S=""; for f in stat_s5 stat_s6 stat_s7 stat_s8 s3u_0 s3u_1 s3u_2 s3u_3; do [ -f $X/$f.npz ] && S="$S $X/$f.npz"; done
  nice -n 10 taskset -c $C $PX -B dqi_int_inv17map2.py $X/dqi.kmb $X/segments_real.json $MB/dqi.kmb $O/segments_real.json $F/dqi.kmb $D/maps/$N $S 2>&1 | cut -c1-220 >> $G
  MAP=$D/maps/${N}_A.npz
fi
A=$D/r46/${N}A; mkdir -p $A
nice -n 10 taskset -c $C $PX -B r44_xfer_dqidrop.py build2 $F/dqi.kmb $MAP $A/dqi.kmb > $A/build2.log 2>&1
echo "build2 A ($MAP): $(tail -1 $A/build2.log | cut -c1-140)" >> $G
echo "file A $(sha256sum $A/dqi.kmb | cut -c1-16) $(date -u +%H:%MZ)" >> $G
nice -n 10 taskset -c $C $PW qcc_verify.py $A --tasks dqi --quick --workers $(pickv) --json $A/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/A quick: /" >> $G
for salt in "" alt-salt-1; do
  tag=${salt:-salt0}
  nice -n 10 taskset -c $C $PW qcc_verify.py $A --tasks dqi --salt "$salt" --workers $(pickv) --json $A/verify_$tag.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/A official $tag: /" >> $G
done
nice -n 10 taskset -c $C $PW -B dqi_int_rate.py $A/dqi.kmb $LINE salt0,alt-salt-1 $(pickp) $A/paired.tsv > $A/paired.log 2>&1
cut -f1-8,10 $A/paired.tsv | sed "s/^/paired A vs best_q24k3a: /" >> $G
echo "DONE $(date -u +%H:%MZ)" >> $G
