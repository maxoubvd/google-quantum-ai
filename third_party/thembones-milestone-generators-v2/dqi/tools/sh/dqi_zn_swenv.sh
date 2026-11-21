#!/usr/bin/env bash
# Copy of S3's tools/sh/dqi_sw_env.sh for builder Z1 (read-only there): inputs and tools paths overridable (ZN_INPUTS, ZN_TOOLS).
# S3 sweep: run a command in the build environment of B5's recipe ip9nd (copied from tools/sh/dqi_int_env.sh,
# read-only there), with every input taken from MY folder ~/qcc/ws/dqi_sw/inputs (copies of B5's kit copy
# source_cap902, D1's frozen snapshots and J's fused-step tensors), so it runs on a cloud box as well.
#   bash tools/sh/dqi_sw_env.sh <cores> <command ...>      (niced; tools from this PC's repo or ~/Quantum/tools)
set -u
CORES=${1:?cores}; shift
I=${ZN_INPUTS:-$HOME/qcc/ws/dqi_sw/inputs}
export DQI_BUILD_PATCHES=${DQI_BUILD_PATCHES:-prefix5,t1,r4jt,r4k,r5w,r6n} DQI_BUILD_FIRST=3
export DQI_BUILD_R4JT=$I/jt_tensors_lens.json
export DQI_JUNC_FWD_LEVELS= DQI_JUNC_REV_KINDS=
export DQI_BUILD_SOURCE=${DQI_SW_SOURCE:-$I/source_cap902} DQI_JUNC_R8_CEIL=${DQI_JUNC_R8_CEIL:-903} DQI_JUNC_SHARE=${DQI_JUNC_SHARE:-1}
export DQI_JUNC_R8P=${DQI_JUNC_R8P:-ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl}
export DQI_X12D_R8P=${DQI_X12D_R8P-}
export DQI_INT_P=${DQI_INT_P-}
T=<HOME>/Quantum/tools
[ -d $T ] || T=$HOME/Quantum/tools; [ -n "${ZN_TOOLS:-}" ] && T=$ZN_TOOLS
cd $T
exec nice -n ${DQI_SW_NICE:-15} taskset -c "$CORES" "$@"
