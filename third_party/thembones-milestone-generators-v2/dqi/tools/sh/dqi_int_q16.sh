#!/usr/bin/env bash
# dqi_int_q16.sh -- workstream dqi_int (B5, Nov 14 15:50Z): q16 = best_q15b + R36-DQI-LOC's cross-word alias (Z1's
# tools/dqi_zn_alias.py: patch A at levels 6..16, ZN_ALIAS_ROOT=1 the root, ZN_ALIAS_B=3,4,5 patch B), installed after
# pdl2 / kpark and BEFORE Z2's pdr (Z1 frv11ac: after pdr it reaches no converter scope). Runs ON box 6, detached.
# The overlay snapshot SNAP (default snap_rv6a) = snap_rv6's modules + a LF copy of ~/Quantum/tools/dqi_zn_alias.py +
# my dqi_int_aliaslog.py (atexit print of the alias LOG), with SHA256SUMS.txt; made once, checked on every build.
#   bash dqi_int_q16.sh <fresh name> <c1,c2> <c1-c2> [VAR=value ...]     (extra settings are exported last)
# Everything else is REBUILD_q15b.sh's recipe (pin q8u, B7's route_3e862d11, rinv27, rv6 LATECONV, shiftscale <= 7).
# (16:05Z) SNAP=<new name> SNAP_FROM=<existing snapshot, default snap_rv6a> SNAP_ADD="<tools/*.py names>": a new
# snapshot = SNAP_FROM + LF copies of those files (e.g. SNAP=snap_rv6b SNAP_ADD=dqi_int_pdrfix.py for q17).
set -u
N=${1:?fresh name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; shift 3
D=~/qcc/ws/dqi_int; V=~/qcc/ws/dqi_var; F=$V/snap/route_3e862d11; S6=$D/zsnap/snap_rv6; S=$D/zsnap/${SNAP:-snap_rv6a}
if [ ! -d $S ] && [ -n "${SNAP_ADD:-}" ]; then
  SF=$D/zsnap/${SNAP_FROM:-snap_rv6a}
  ( cd $SF && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$SF differs"; exit 1; }
  rm -rf $S.tmp; mkdir -p $S.tmp; cp $SF/*.py $SF/USAGE.txt $S.tmp/
  for f in $SNAP_ADD; do tr -d '\r' < ~/Quantum/tools/$f > $S.tmp/$f; done
  echo "B5 $(date -u +%FT%H:%MZ): $(basename $SF) + $SNAP_ADD (LF copies of ~/Quantum/tools as of this time)" >> $S.tmp/USAGE.txt
  ( cd $S.tmp && sha256sum *.py > SHA256SUMS.txt ) && mv $S.tmp $S
fi
if [ ! -d $S ]; then
  rm -rf $S.tmp; mkdir -p $S.tmp
  cp $S6/*.py $S6/USAGE.txt $S.tmp/
  tr -d '\r' < ~/Quantum/tools/dqi_zn_alias.py > $S.tmp/dqi_zn_alias.py
  tr -d '\r' < ~/Quantum/tools/dqi_int_aliaslog.py > $S.tmp/dqi_int_aliaslog.py
  echo "B5 $(date -u +%FT%H:%MZ): snap_rv6 + Z1's dqi_zn_alias.py (as of this time, LF) + dqi_int_aliaslog.py" >> $S.tmp/USAGE.txt
  ( cd $S.tmp && sha256sum *.py > SHA256SUMS.txt ) && mv $S.tmp $S
fi
( cd $F && sha256sum -c --quiet MANIFEST.sha256 ) || { echo "route differs"; exit 1; }
( cd $S && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$S differs"; exit 1; }
( cd $D/pin/q8u && sha256sum -c --quiet MANIFEST.sha256 ) || { echo "pin q8u differs"; exit 1; }
sha256sum ~/Quantum/tools/dqi_var_rinv27.py | grep -q ^6954d6f0 || { echo "rinv27 differs"; exit 1; }
sha256sum $S/dqi_zn_mkunf5.py | grep -q ^38f96e37 || { echo "mkunf5 differs"; exit 1; }
tr -d '\r' < ~/Quantum/tools/sh/dqi_int_pinvar.sh > $D/pinvar.sh
SW="R12_FP=a,h DQI_INT_EXT=dqi_zn_pdl2,dqi_zn_kpark,dqi_z2_pdr ZN_PDL_BLOCKS=17,15,13,11,9 ZN_PDL2_FWDONLY=1 ZN_PDL2_SHIFT=16,14,12,10,8 ZN_PDL2_SHIFTPARK=1 ZN_PDL2_YSHARE=1 ZN_PDL2_HAND=7,8,9,10 Z2_PDR_SHARE=1 Z2_PDR_BLOCKS=15,13,11 Z2_PDR_RBLOCKS=15,13,11 Z2_PDR_REV=1 Z2_PDR_LOCAL=1 Z2_PDR_LIVEB=1 DQI_VAR_ROUTE_LOG=1"
EX="R12_FP=h ZN_PDL2_FWDONLY= ZN_PDL2_NPARK=1 DQI_INT_PP2_NOINPUT=1 ZN_PDL2_HAND=7,8,9,10,11,12,13 DQI_INT_SS_OFF=1 R12_SOURCE=$D/pinvar/$N/home/qcc/ws/dqi_zn/src/kv_unfused5"
Q15="ZN_PDL2_LATECONV=1 ZN_PDL2_SHIFTPARK= DQI_INT_SS_OFF= DQI_INT_SS_MAXM=7"
AL="DQI_INT_EXT=dqi_zn_pdl2,dqi_zn_kpark,dqi_zn_alias,dqi_int_aliaslog,dqi_z2_pdr ZN_ALIAS_ROOT=1 ZN_ALIAS_B=3,4,5"
mkdir -p $D/pinvar
cd $D && ROUTE_FROZEN=$F/dqi_var_route.py OVERLAY=$S KIT_MK=$S/dqi_zn_mkunf5.py KEEP_KMB=1 setsid nohup bash $D/pinvar.sh \
  q8u $N $CC $PC dqi_var_rinv27.py var_rinv27,var_route $SW $EX $Q15 $AL "$@" > $D/pinvar_$N.txt 2>&1 < /dev/null &
echo "q16 build $N started on $CC (P1 $PC), overlay $S"
