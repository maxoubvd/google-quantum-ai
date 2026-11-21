#!/usr/bin/env bash
# dqi_int_t54.sh -- B5 (18th holder, Nov 17; runs ON box 6): best_q27g100's whole 901 chain (dqi_int_r901repro.sh,
# fresh folders) with R54-DQI-TENSOR's even step (173 terms) as the r4jt tensor file (B5_R4JT, dqi_int_pinvar.sh),
# started detached. Output lines: ~/qcc/ws/dqi_int/repro901_<name>.txt; stdout ~/qcc/ws/dqi_int/bg/<name>.out.
#   bash dqi_int_t54.sh <fresh name> <c1,c2> <c1-c2> [MAXGAP=100] [tensor json] [extra VAR=V ...]
set -u
N=${1:?fresh name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; MG=${4:-100}; JT=${5:-$HOME/qcc/ws/r54_dqi_tensor/out/jt_tensors_r54.json}
shift; shift; shift; [ $# -gt 0 ] && shift; [ $# -gt 0 ] && shift
D=$HOME/qcc/ws/dqi_int; mkdir -p $D/bg
[ -s $JT ] || { echo "no tensor file $JT"; exit 1; }
echo "t54 $N tensor $JT $(sha256sum $JT | cut -c1-16) MAXGAP $MG extra $* start $(date -u +%H:%M:%SZ)" >> $D/bg/$N.out
B5_R4JT=$JT setsid nohup bash $HOME/Quantum/tools/sh/dqi_int_r901repro.sh $N $CC $PC l9r4oxg \
  B5_ZP_ORACLE=$D/zptrace/z9C_oracle.json B5_ZP_DEFER=1 B5_ZP_DOOM=1 B5_ZP_BELADY=1 B5_ZP_VICTIMS=1 B5_ZP_MAXGAP=$MG "$@" \
  >> $D/bg/$N.out 2>&1 < /dev/null &
echo "started $N on $CC (pid $!)"
