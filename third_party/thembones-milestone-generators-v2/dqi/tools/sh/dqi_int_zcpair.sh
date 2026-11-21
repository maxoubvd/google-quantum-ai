#!/usr/bin/env bash
# dqi_int_zcpair.sh -- B5 (22nd holder, Nov 18; runs ON box 6, detached): waits for a 901 chain's f901 DONE, then the
# paired failure test of its file A against the PUBLISHED line's file A (best_q30mhrE12 = ~/qcc/ws/dqi_int/p901/mhrE12/A,
# fbe6e5dd) on twelve salts (tools/dqi_int_rate.py; CHILD-ONLY must read 0). dqi_int_mhrpair.sh with the parent moved.
#   [B5_ZC_PARENT=<p901 tag> B5_ZC_PTAG=<label>] bash dqi_int_zcpair.sh <tag> <cores> <workers> [salts]   lines -> ~/qcc/ws/dqi_int/pairzc/<tag>_vs_e12A.tsv
set -u
N=$1; C=$2; W=${3:-2}; S=${4:-salt0,alt-salt-1..11}; PA=${B5_ZC_PARENT:-mhrE12}; PT=${B5_ZC_PTAG:-e12A}
D=~/qcc/ws/dqi_int; mkdir -p $D/pairzc
until grep -a -q '^DONE' $D/f901_$N.txt 2>/dev/null; do sleep 30; done
cd ~/Quantum/tools && nice -n 10 taskset -c $C ~/qcc/.venv/bin/python -B dqi_int_rate.py $D/p901/$N/A/dqi.kmb \
  $D/p901/$PA/A/dqi.kmb $S $W $D/pairzc/${N}_vs_$PT.tsv > $D/pairzc/${N}_vs_$PT.log 2>&1
echo "DONE $(date -u +%H:%M:%SZ)" >> $D/pairzc/${N}_vs_$PT.log
