#!/usr/bin/env bash
# dqi_int_kcgpair.sh -- B5 (20th holder, Nov 17; runs ON box 6, detached): waits for a 901 chain's f901 DONE, then the
# paired failure test of its file A against the published line's file A (t54chz, 880403d3) on twelve salts
# (tools/dqi_int_rate.py: the child's own official draws, every failing shot named; CHILD-ONLY must read 0).
#   bash dqi_int_kcgpair.sh <tag> <cores> <workers>      lines -> ~/qcc/ws/dqi_int/pairkcg/<tag>_vs_chzA.tsv
set -u
N=$1; C=$2; W=${3:-2}
D=~/qcc/ws/dqi_int; mkdir -p $D/pairkcg
until grep -a -q '^DONE' $D/f901_$N.txt 2>/dev/null; do sleep 30; done
cd ~/Quantum/tools && nice -n 10 taskset -c $C ~/qcc/.venv/bin/python -B dqi_int_rate.py $D/p901/$N/A/dqi.kmb \
  $D/p901/t54chz/A/dqi.kmb salt0,alt-salt-1..11 $W $D/pairkcg/${N}_vs_chzA.tsv > $D/pairkcg/${N}_vs_chzA.log 2>&1
echo "DONE $(date -u +%H:%M:%SZ)" >> $D/pairkcg/${N}_vs_chzA.log
