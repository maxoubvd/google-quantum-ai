#!/usr/bin/env bash
# B5: after the pinned build ~/qcc/ws/dqi_int/pinvar/<name> has its raw quick check, run P1's dwx pass on its raw file
# into ~/qcc/ws/dqi_int/dwx/<name>/ and write one line to ~/qcc/ws/dqi_int/rebuild_<name>.txt.
#   bash dqi_int_dwxafter.sh <name> <cores>
set -u
N=$1; C=$2; D=~/qcc/ws/dqi_int; R=$D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u
until [ -s $R/quick.log ] || grep -a -q '^end_' $D/pinvar/$N/RESULT.txt 2>/dev/null; do sleep 30; done
cd ~/Quantum && nice -n 10 bash tools/sh/dqi_pp_dwx.sh $R/dqi.kmb $D/dwx/$N $C > $D/logs/dwx_$N.log 2>&1
echo "rebuild $N: quick $(head -1 $R/quick.log | cut -c1-80) raw $(sha256sum $R/dqi.kmb | cut -c1-16) final $(sha256sum $D/dwx/$N/final/dqi.kmb 2>/dev/null | cut -c1-16) $(date -u +%H:%MZ)" >> $D/rebuild_$N.txt
