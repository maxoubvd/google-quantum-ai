#!/usr/bin/env bash
# dqi_int_b901.sh -- B5 (17th holder, Nov 16; runs ON box 6): the 901 file from one MARKER build and a park spec:
# converter parks (dqi_int_park901b.py), markers stripped (dqi_int_strip.py), exact zero intervals (P1's
# dqi_pp_zero.py), wire re-assignment to Q (dqi_int_realloc.py), quick verify.
# B5_PARK_TOOL=dqi_int_park901c.py (with B5_PARK_MOVEZ=1: diagonal trans ops moved after the remake);
# B5_REALLOC_PIN=ends (re-assignment pins only a syndrome wire's input and output segments).
# B5_INV17_MAP=dqi_int_inv17map3.py (23rd holder, Nov 18): the site map that matches stages by frames (default map2).
# B5_INV17=1: INV17's list A applied on the stripped parked file (its wire ids are the marker build's, so the site map
# aligns through the converter's parks), before the zero intervals and the re-assignment; then skip INV17 after P1.
#   bash dqi_int_b901.sh <marker build name> <spec file> <cores> <out tag> [Q=901]
# Output: ~/qcc/ws/dqi_int/p901/<out tag>/{k.kmb,ks.kmb,ks_zero.npz,r<Q>/dqi.kmb,log.txt}
set -u
N=${1:?marker build}; SPEC=$(readlink -f ${2:?spec}); C=${3:?cores}; TAG=${4:?out tag}; Q=${5:-901}
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PX=~/qcc/.venv_work/bin/python; PW=~/qcc/.venv/bin/python
K=$D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u/dqi.kmb
O=$D/p901/$TAG; mkdir -p $O; L=$O/log.txt; cd $T
P="nice -n 10 taskset -c $C"
echo "b901 $N spec $SPEC ($(sha256sum $SPEC | cut -c1-8)) Q $Q start $(date -u +%H:%M:%SZ)" >> $L
$P $PX -B ${B5_PARK_TOOL:-dqi_int_park901b.py} $K $SPEC $O/k.kmb > $O/park.log 2>&1 || { echo "park FAILED: $(tail -3 $O/park.log)" >> $L; exit 1; }
echo "park: $(tail -2 $O/park.log | tr '\n' ' ' | cut -c1-300)" >> $L
$P $PX -B dqi_int_strip.py $O/k.kmb $O/ks.kmb > $O/strip.log 2>&1 || { echo "strip FAILED" >> $L; exit 1; }
echo "strip: $(tail -2 $O/strip.log | tr '\n' ' ' | cut -c1-200)" >> $L
KS=$O/ks.kmb
if [ "${B5_INV17:-0}" = 1 ]; then
  X=~/qcc/ws/r46_dqi_ext/prof1; MS=$D/p901/$N/segments_real.json
  [ -s $MS ] || { mkdir -p $D/p901/$N; $P $PX -B r46_dqi_ext_segreal.py $D/r46/$N/dqi.kmb $D/r46/$N/segments.json $MS > $D/p901/$N/segreal.log 2>&1; }
  S=""; for f in stat_s5 stat_s6 stat_s7 stat_s8 s3u_0 s3u_1 s3u_2 s3u_3; do [ -f $X/$f.npz ] && S="$S $X/$f.npz"; done
  mkdir -p $D/maps
  $P $PX -B ${B5_INV17_MAP:-dqi_int_inv17map2.py} $X/dqi.kmb $X/segments_real.json $D/r46/$N/dqi.kmb $MS $O/ks.kmb $D/maps/b901_$TAG $S > $O/inv17.log 2>&1
  echo "inv17 map: $(grep -a -e aligned -e '_A.npz' $O/inv17.log | tr '\n' ' ' | cut -c1-250)" >> $L
  $P $PX -B r44_xfer_dqidrop.py build2 $O/ks.kmb $D/maps/b901_${TAG}_A.npz $O/ksi.kmb > $O/build2.log 2>&1 || { echo "build2 FAILED" >> $L; exit 1; }
  echo "build2: $(tail -1 $O/build2.log | cut -c1-200)" >> $L
  KS=$O/ksi.kmb
fi
$P $PX -B dqi_pp_zero.py $KS $O/ks_zero.npz > $O/zero.log 2>&1 || { echo "zero FAILED" >> $L; exit 1; }
mkdir -p $O/r$Q
$P $PX -B dqi_int_realloc.py $KS $O/ks_zero.npz $O/r$Q/dqi.kmb $Q > $O/realloc.log 2>&1 || { echo "realloc FAILED: $(tail -3 $O/realloc.log | tr '\n' ' ')" >> $L; exit 1; }
echo "realloc: $(tail -2 $O/realloc.log | tr '\n' ' ' | cut -c1-300)" >> $L
$P $PW qcc_verify.py $O/r$Q --tasks dqi --quick --workers 2 --json $O/r$Q/quick.json 2>&1 | grep -a "dqi  *cost" | tail -1 | sed "s/^/quick r$Q: /" >> $L
echo "sha r$Q $(sha256sum $O/r$Q/dqi.kmb | cut -c1-16) done $(date -u +%H:%M:%SZ)" >> $L
