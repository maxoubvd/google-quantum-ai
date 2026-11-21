#!/usr/bin/env bash
# dqi_int_c901.sh -- B5 (17th holder, Nov 16; runs ON box 6): the 902 census of one MARKER build for the 901 file:
# exact zero intervals (P1's dqi_pp_zero.py), real segments (r46_dqi_ext_segreal.py), idle holders of every run at
# live >= LEVEL with the greedy flag-park cover (RET2 re-pricer's r48_dqi_ret2_rp_idle.py), and the flag checks of the
# cover's parks (r48_dqi_ret2_rp_flagchk.py; census = dqi_int_idle901.py, a copy with B5_IDLE_EXCLUDE=<wires>: the tested word, writes to it, samples at open / close / room op).
#   bash dqi_int_c901.sh <marker build name | path to a marker dqi.kmb with segments.json beside it> <cores> [LEVEL=902]
# Output: ~/qcc/ws/dqi_int/p901/<name>/{zero.npz,segments_real.json,idle<LEVEL>.txt,flagchk<LEVEL>.txt}
set -u
N=${1:?name}; C=${2:?cores}; L=${3:-902}
D=~/qcc/ws/dqi_int; T=~/Quantum/tools; PX=~/qcc/.venv_work/bin/python
if [ -f "$N" ]; then K=$(readlink -f $N); NM=$(basename $(dirname $K)); else K=$D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u/dqi.kmb; NM=$N; fi
SJ=$(dirname $K)/segments.json
O=$D/p901/$NM; mkdir -p $O; cd $T
P="nice -n 10 taskset -c $C"
[ -s $O/zero.npz ] || $P $PX -B dqi_pp_zero.py $K $O/zero.npz > $O/zero.log 2>&1
[ -s $O/segments_real.json ] || $P $PX -B r46_dqi_ext_segreal.py $K $SJ $O/segments_real.json > $O/segreal.log 2>&1
$P $PX -B dqi_int_idle901.py $K $O/zero.npz extract,synthesis,convert,scan,unconvert,padded_unsynthesis,unextract $L $O/segments_real.json 80 > $O/idle$L.txt 2>&1
$P $PX -B r48_dqi_ret2_rp_flagchk.py $K $O/segments_real.json $O/idle$L.txt 16 7 $O/zero.npz $L > $O/flagchk$L.txt 2>&1
echo "c901 $NM level $L done $(date -u +%H:%M:%SZ)" >> $O/done.txt
