#!/usr/bin/env bash
# dqi_int_zcgo.sh -- B5 (22nd holder, Nov 18; runs ON box 6): one zcache chain = best_q30mhrE12's recipe (dqi_int_mhr_go.sh
# with E at levels 11, 13) + tools/dqi_int_zcache.py installed last at ceiling CEIL, into FRESH folders, and its twelve-salt
# pair waiter against best_q30mhrE12's file A (dqi_int_zcpair.sh). Never start two of these in the same minute: the
# launcher rewrites the shared ~/qcc/ws/dqi_int/pinvar.sh.
#   [B5_ZC_FILES=a.py,b.py B5_ZC_EXT=a,b B5_ZC_VARS='X=1 Y=2'] bash dqi_int_zcgo.sh <fresh> <c1,c2> <c1-c2> <CEIL> <pair cores>
set -u
N=${1:?}; CC=${2:?}; PC=${3:?}; CE=${4:?}; PP=${5:?}
D=~/qcc/ws/dqi_int
B5_MHR_FILES=dqi_int_zcache.py${B5_ZC_FILES:+,$B5_ZC_FILES} B5_MHR_EXT=dqi_int_zcache${B5_ZC_EXT:+,$B5_ZC_EXT} bash ~/Quantum/tools/sh/dqi_int_mhr_go.sh $N $CC $PC \
  DQI_INT_MHR_E=1 DQI_INT_MHR_XOFF=1 DQI_INT_MHR_LEVELS=11,13 DQI_INT_ZC_CEIL=$CE DQI_INT_ZC_OUT=$D/mhr/${N}_zc.json DQI_INT_RC_OUT=$D/mhr/${N}_rc.json ${B5_ZC_VARS:-}
setsid nohup bash ~/Quantum/tools/sh/dqi_int_zcpair.sh $N $PP 2 > /dev/null 2>&1 < /dev/null &
echo "pair waiter for $N on $PP"
