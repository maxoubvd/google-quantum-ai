#!/bin/bash
# dqi_kern_jtest.sh -- run a junction-harness script (workstream dqi_kern, builder B6) with the ip8 build environment.
#   bash dqi_kern_jtest.sh <script.py> [args...]
# Environment: the kit copy with ceiling C (default 903), the build line's patches as in ip8/ip9, my patches in
# DQI_JUNC_R8P (default ipjunc). DQI_JUNC_SHARE defaults to 1 (ip8's shared return).
T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
C=${C:-903}
export DQI_BUILD_SOURCE=~/qcc/ws/dqi_junc/r8/src/source_cap$C
export DQI_BUILD_PATCHES=${DQI_BUILD_PATCHES:-prefix5,t1,r4jt,r4k,r5w,r6n}
export DQI_BUILD_FIRST=3
export DQI_BUILD_R4JT=${DQI_BUILD_R4JT:-$HOME/qcc/ws/dqi_build/r8/jt_tensors_lens.json}
export DQI_JUNC_R8P=${DQI_JUNC_R8P-ipjunc}
export DQI_JUNC_SHARE=${DQI_JUNC_SHARE:-1}
export DQI_JUNC_R8_CEIL=${DQI_JUNC_R8_CEIL:-$C}
mkdir -p ~/qcc/ws/dqi_kern/logs
cd $T || exit 1
S=$1; shift
nice -n 10 taskset -c ${CORES:-8,9} $PY -B $T/$S "$@"
