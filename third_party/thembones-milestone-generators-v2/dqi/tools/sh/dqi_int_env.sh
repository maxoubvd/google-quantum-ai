#!/usr/bin/env bash
# Workstream dqi_int (builder B5): run a command in the environment of a 24-bit recipe (builder J's recipes, copied
# from tools/sh/dqi_x12d_env24.sh; J's kit copies and patches used read-only), niced, on the given cores.
#   bash tools/sh/dqi_int_env.sh <ip8|ip9|r1> <cores> <command ...>
# r1 = ip9 + x12d's scanword (DQI_X12D_R8P) + my patches (DQI_INT_P, comma list of tools/dqi_int_<name>.py), entry
# tools/dqi_int_entry.py. The `peak` diagnostic (no gate change) is in the build patches unless DQI_INT_NOPEAK=1.
set -u
REC=${1:?recipe}; CORES=${2:?cores}; shift 2
BASEP=prefix5,t1,r4jt,r4k,r5w,r6n
[ "${DQI_INT_NOPEAK:-0}" = 1 ] || BASEP=$BASEP,peak
export DQI_BUILD_PATCHES=${DQI_BUILD_PATCHES:-$BASEP} DQI_BUILD_FIRST=3
# on a cloud box the tensors are my copy (DQI_INT_R4JT=~/qcc/ws/dqi_int/src/jt_tensors_lens.json)
export DQI_BUILD_R4JT=${DQI_INT_R4JT:-$HOME/qcc/ws/dqi_build/r8/jt_tensors_lens.json}
export DQI_JUNC_FWD_LEVELS=${DQI_JUNC_FWD_LEVELS-} DQI_JUNC_REV_KINDS=${DQI_JUNC_REV_KINDS-}   # empty unless set (sweeps)
case $REC in
  ip8)  export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap905 DQI_JUNC_R8_CEIL=905 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b ;;
  ip9|r1)
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap903 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,div6,tailctl ;;
  ip9n) # ip9 without div6 (the un-synthesis division's narrow 6-bit inverse)
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap903 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl ;;
  ip9nc) # ip9n with the kit copy capped at 895 (the kit's ceiling-filled kernels pick narrower variants; J's
        # ceiling stays 903): a probe of which stages still need 903
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap895 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl ;;
  ip9nd) # ip9n with my kit copy capped at 902 (~/qcc/ws/dqi_int/src/source_cap902: J's cap903 copy with the
        # ceiling literal 903 -> 902 in its four files); J's ceiling stays 903
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_int/src/source_cap902 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl ;;
  ip9d) # ip9 (with div6) on my kit copy capped at 902: every stage but the converter at 902
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_int/src/source_cap902 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,div6,tailctl ;;
  ip9f) # rung 2 probe: ip9 (with div6) on my kit copy capped at 901 (source_cap901: the four ceiling literals 901)
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_int/src/source_cap901 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,div6,tailctl ;;
  ip9ne) # ip9nd with J's ceiling at 902 as well (J's kernels choose their variants against 902)
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_int/src/source_cap902 DQI_JUNC_R8_CEIL=902 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl ;;
  ip9nx) # ip9nd on a kit copy capped at DQI_INT_CAP (src/source_cap<N>, made by tools/sh/dqi_int_kitcap.sh): a
        # probe that measures a wider kernel's Toffolis before the qubits for it exist
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_int/src/source_cap${DQI_INT_CAP:?} DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,tailctl ;;
  ip9t) # ip9 without tailctl (the synthesis tail's factor-free correction)
        export DQI_BUILD_SOURCE=$HOME/qcc/ws/dqi_junc/r8/src/source_cap903 DQI_JUNC_R8_CEIL=903 DQI_JUNC_SHARE=1
        export DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b,div6 ;;
  *) echo "unknown recipe $REC"; exit 1 ;;
esac
[ -n "${DQI_INT_R8P:-}" ] && export DQI_JUNC_R8P=$DQI_INT_R8P   # sweeps: replace the recipe's J patch list
[ "${DQI_INT_NOPEAK:-0}" = 1 ] || export DQI_JUNC_R8P=$DQI_JUNC_R8P,truepeak     # J's no-gate diagnostic
if [ "$REC" = r1 ]; then export DQI_X12D_R8P=${DQI_X12D_R8P-scanword}; else export DQI_X12D_R8P=${DQI_X12D_R8P-}; fi
export DQI_INT_P=${DQI_INT_P-}
if [ -d <HOME>/Quantum/tools ]; then cd <HOME>/Quantum/tools
else cd $HOME/Quantum/tools; fi      # cloud box: tools pushed by box_sync.sh push
exec nice -n 10 taskset -c "$CORES" "$@"
