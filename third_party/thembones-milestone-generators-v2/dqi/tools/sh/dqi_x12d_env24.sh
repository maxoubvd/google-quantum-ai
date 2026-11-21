#!/usr/bin/env bash
# Workstream dqi_x12d: run a command in the environment of one of builder J's 24-bit recipes (read-only use of J's
# kit copies and patches), niced, on one core.
#   bash tools/sh/dqi_x12d_env24.sh <ip4c907|ip8|ip9> <core> <command ...>
# My own patches are added with DQI_X12D_R8P (entry tools/dqi_x12d_r8_entry.py); set it before calling.
set -u
REC=${1:?recipe}; CORE=${2:?core}; shift 2
export DQI_BUILD_PATCHES=${DQI_BUILD_PATCHES:-prefix5,t1,r4jt,r4k,r5w,r6n,peak} DQI_BUILD_FIRST=3
export DQI_BUILD_R4JT=$HOME/qcc/ws/dqi_build/r8/jt_tensors_lens.json
export DQI_JUNC_FWD_LEVELS= DQI_JUNC_REV_KINDS=
case $REC in
  ip4c907) export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap907 DQI_JUNC_R8_CEIL=907
           export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem ;;
  ip8)     export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap905 DQI_JUNC_R8_CEIL=905 DQI_JUNC_SHARE=1
           export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b ;;
  ip9)     export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap903 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
           export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,div6,tailctl ;;
  *) echo "unknown recipe $REC"; exit 1 ;;
esac
[ -n "${DQI_X12D_SRC:-}" ] && export DQI_BUILD_SOURCE=$DQI_X12D_SRC     # a copy of mine instead (measurement builds)
cd <HOME>/Quantum/tools
exec nice -n 10 taskset -c "$CORE" "$@"
