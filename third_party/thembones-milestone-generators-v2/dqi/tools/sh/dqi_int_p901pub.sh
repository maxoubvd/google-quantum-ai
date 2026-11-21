#!/usr/bin/env bash
# dqi_int_p901pub.sh -- B5 (17th holder, Nov 16; runs ON box 6, detached): publication steps for a 901 file A made by
# dqi_int_f901.sh (~/qcc/ws/dqi_int/p901/<tag>/A/dqi.kmb): waits for f901_<tag>.txt DONE, hard-links file A as
# ~/qcc/ws/dqi_int/r46/<tag>A/dqi.kmb, then re-rolls (dqi_int_r46rr.sh, two clean on both official salts wanted) and 10
# further salts (alt-salt-2..11) on file A.
#   bash dqi_int_p901pub.sh <tag> <cores> <workers>          lines -> ~/qcc/ws/dqi_int/pub901_<tag>.txt
set -u
N=$1; C=$2; W=$3
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PW=~/qcc/.venv/bin/python; G=$D/pub901_$N.txt; A=$D/r46/${N}A
until grep -a -q '^DONE' $D/f901_$N.txt 2>/dev/null; do sleep 30; done
mkdir -p $A; ln -f $D/p901/$N/A/dqi.kmb $A/dqi.kmb
echo "file A $(sha256sum $A/dqi.kmb | cut -c1-16) $(date -u +%H:%MZ)" >> $G
bash $T/sh/dqi_int_r46rr.sh ${N}A $C $W 2 a b c d e f g h
echo "re-rolls: $(tail -1 $D/rr_${N}A.txt)" >> $G
cd $T
for i in 2 3 4 5 6 7 8 9 10 11; do
  nice -n 10 taskset -c $C $PW qcc_verify.py $A --tasks dqi --salt alt-salt-$i --workers $W --json $A/verify_alt$i.json 2>&1 \
    | grep -a "dqi  *cost" | tail -1 | sed "s/^/alt-salt-$i: /" >> $G
done
echo "DONE $(date -u +%H:%MZ)" >> $G
