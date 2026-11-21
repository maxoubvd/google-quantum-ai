#!/usr/bin/env bash
# dqi_int_zc5b.sh -- B5 (23rd holder, Nov 18; runs ON box 6). zc5b = best_q31zc4's recipe (zcache + rinvc + iicache on
# best_q30mhrE12's recipe) + DQ3's three rows installed LAST: npair (tools/dqi_dq3_npair.py), ff24 (tools/dqi_dq3_w901.py
# with DQI_DQ3_W901_NAMES=inverse24_forward_flagfree), hflag (tools/dqi_dq3_hflag.py). Fresh folders only; its twelve-salt
# pair waiter runs against zc4's file A (best_q31zc4). Never two chain launches in the same minute.
#   bash dqi_int_zc5b.sh <fresh> <c1,c2> <c1-c2> <pair cores> [extra files] [extra exts] [extra vars]
set -u
N=${1:?fresh}; CC=${2:-2,3}; PC=${3:-2-3}; PP=${4:-8,9}
XF=${5:-}; XE=${6:-}; XV=${7:-}
export B5_ZC_FILES=dqi_int_rinvc.py,dqi_int_iicache.py,dqi_dq3_npair.py,dqi_dq3_w901.py,dqi_dq3_hflag.py${XF:+,$XF}
export B5_ZC_EXT=dqi_int_rinvc,dqi_int_iicache,dqi_dq3_npair,dqi_dq3_w901,dqi_dq3_hflag${XE:+,$XE}
export B5_ZC_VARS="DQI_DQ3_W901_NAMES=inverse24_forward_flagfree${XV:+ $XV}"
# B5_INV17_MAP=dqi_int_inv17map3.py (env, from zc11 on): INV17 site map by frames; unset = map2 (the published files)
export B5_ZC_PARENT=${B5_ZC_PARENT:-zc4} B5_ZC_PTAG=${B5_ZC_PTAG:-zc4A}
echo "$(date -u +%H:%M:%SZ) launch $N files=$B5_ZC_FILES vars=$B5_ZC_VARS parent=$B5_ZC_PARENT" >> ~/qcc/ws/dqi_int/zc_launch.log
bash ~/Quantum/tools/sh/dqi_int_zcgo.sh $N $CC $PC 901 $PP
