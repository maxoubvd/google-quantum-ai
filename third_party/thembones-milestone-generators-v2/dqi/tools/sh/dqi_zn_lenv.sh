#!/usr/bin/env bash
# dqi_zn_lenv.sh -- builder Z1: run a command in the LINE's build environment (r12hz's recipe: B5's
# ~/qcc/ws/dqi_int/r12/b_r12hz_prof/env.txt without the profile markers), through S3's tools/sh/dqi_sw_env.sh.
#   [ZN_SRC=<kit copy>] [ZN_TAIL=<extra DQI_INT_P tokens>] bash dqi_zn_lenv.sh <cores> <command ...>
set -u
CORES=${1:?cores}; shift
I=${ZN_INPUTS:-$HOME/qcc/ws/dqi_sw/inputs}
export DQI_INT_A1_FIX=1 DQI_INT_A1_MODE=kitco DQI_INT_A1_SNAP=$I/d1a1snap_9c05c258 DQI_INT_A3_SNAP=$I/d1snap_577a000b
export DQI_INT_FP=${ZN_FP:-a,b,c,h} DQI_INT_NOPEAK=1
P=retj,polepark2,flagpark,retl,div6f,a3,a1,qerase,var_sync3,var_keep,var_blockhold,var_gauge,pub_pfd,pub_zt,var_ztp
for t in $(echo "${ZN_DROP:-}" | tr ',' ' '); do P=$(echo ",$P," | sed "s/,$t,/,/" | sed 's/^,//; s/,$//'); done
export DQI_INT_P=$P${ZN_TAIL:+,$ZN_TAIL}
export DQI_PUB_A3_CORE=v5 DQI_PUB_A3_DIV=dr DQI_PUB_A3_SYNC=1 DQI_PUB_PFD_P0=1 DQI_PUB_PFD_RINF=1 DQI_PUB_ZJOIN=1
export DQI_SW_KEEP=1 DQI_VAR_BH_LEVELS= DQI_VAR_BH_PRE=13,14 DQI_VAR_GAUGE_MODE=sort DQI_VAR_KEEP_ADD=4,5,6,7
export DQI_BUILD_PATCHES=prefix5,t1,r4jt,r4k,r5w,r6n
L=$HOME/qcc/ws/dqi_zn/src/kv__noa16_cap900_noret
[ -d "$L" ] || L=$I/kv__noa16_cap900_noret
export DQI_SW_SOURCE=${ZN_SRC:-$L}
T=${ZN_TOOLS:-<HOME>/Quantum/tools}; [ -d $T ] || T=$HOME/Quantum/tools
export ZN_INPUTS=$I ZN_TOOLS=$T
tr -d '\r' < $T/sh/dqi_zn_swenv.sh > /tmp/dqi_zn_swenv_$$.sh
exec bash /tmp/dqi_zn_swenv_$$.sh "$CORES" "$@"
