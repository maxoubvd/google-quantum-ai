#!/bin/bash
# dqi_k1_jtest.sh -- K1 (round 3): run a junction-harness script with the ip8/ip9 build environment
# (a copy of B6's tools/sh/dqi_kern_jtest.sh with the tools path set to the workspace; nothing else changed).
#   CORES=6,7 bash dqi_k1_jtest.sh <script.py> [args...]
T=<HOME>/Quantum/tools; PY=~/qcc/.venv/bin/python
C=${C:-903}
export DQI_BUILD_SOURCE=~/qcc/ws/dqi_junc/r8/src/source_cap$C
export DQI_BUILD_PATCHES=${DQI_BUILD_PATCHES:-prefix5,t1,r4jt,r4k,r5w,r6n}
export DQI_BUILD_FIRST=3
export DQI_BUILD_R4JT=${DQI_BUILD_R4JT:-$HOME/qcc/ws/dqi_build/r8/jt_tensors_lens.json}
export DQI_JUNC_R8P=${DQI_JUNC_R8P-latem,div6}
export DQI_JUNC_SHARE=${DQI_JUNC_SHARE:-1}
export DQI_JUNC_R8_CEIL=${DQI_JUNC_R8_CEIL:-$C}
cd $T || exit 1
S=$1; shift
nice -n 10 taskset -c ${CORES:-6,7} $PY -B $T/$S "$@"
