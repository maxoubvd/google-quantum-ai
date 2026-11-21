#!/usr/bin/env bash
# dqi_int_mhrpair.sh -- B5 (21st holder, Nov 17; runs ON box 6, detached): waits for a 901 chain's f901 DONE, then the
# paired failure test of its file A against the published line's file A (kcgABC, ba4862cd = best_q29kcgABC before its
# tag) on twelve salts (tools/dqi_int_rate.py: the child's own official draws, every failing shot named; CHILD-ONLY
# must read 0). A copy of dqi_int_kcgpair.sh with the parent moved to today's line.
#   bash dqi_int_mhrpair.sh <tag> <cores> <workers> [salts]   lines -> ~/qcc/ws/dqi_int/pairmhr/<tag>_vs_abcA.tsv
set -u
N=$1; C=$2; W=${3:-2}; S=${4:-salt0,alt-salt-1..11}
D=~/qcc/ws/dqi_int; mkdir -p $D/pairmhr
until grep -a -q '^DONE' $D/f901_$N.txt 2>/dev/null; do sleep 30; done
cd ~/Quantum/tools && nice -n 10 taskset -c $C ~/qcc/.venv/bin/python -B dqi_int_rate.py $D/p901/$N/A/dqi.kmb \
  $D/p901/kcgABC/A/dqi.kmb $S $W $D/pairmhr/${N}_vs_abcA.tsv > $D/pairmhr/${N}_vs_abcA.log 2>&1
echo "DONE $(date -u +%H:%M:%SZ)" >> $D/pairmhr/${N}_vs_abcA.log
