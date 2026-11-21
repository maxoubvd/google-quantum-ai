#!/usr/bin/env bash
# dqi_int_r48b.sh -- B5 (13th holder, Nov 16; runs ON box 6, detached): best_q24k3a's recipe (k3 of dqi_int_r46b.sh)
# + R48-DQI-RET2's late flag lf3 (r48_dqi_ret2_lateflag with R48_KPOINT=1 R48_N3FLAG=1, the copies RET2 built lf3 from)
# + my root block module dqi_int_pdr17 in Z2's pdr's place, into FRESH folders; then P1's dwx pass (dwxafter.sh).
#   bash dqi_int_r48b.sh <fresh name> <c1,c2> <c1-c2> <variant>[p] [VAR=V ...]
# variants: l3 = k3 + lf3 (pdr17 without block 17 = Z2's pdr: the control, raw must equal RET2's lf3 fa222822)
#           f  = l3 + block 17 forward (Z2_PDR_BLOCKS=17,15,13,11)
#           fr = l3 + block 17 both ways (Z2_PDR_BLOCKS=Z2_PDR_RBLOCKS=17,15,13,11)
#           suffix p = the same with r7prof markers (no P1, no dwx), for INV17's site map
# Overlay (made once per pdr17 version, sha-checked): zsnap/snap_b5r48_<pdr17 sha8> = snap_b5r46n + RET2's frozen
# lf3 copies of r48_dqi_ret2_lateflag.py / r48_dqi_ret2_forms.py + an LF copy of ~/Quantum/tools/dqi_int_pdr17.py.
# Result line: ~/qcc/ws/dqi_int/rebuild_<name>.txt; final file ~/qcc/ws/dqi_int/dwx/<name>/final/dqi.kmb.
set -u
N=${1:?fresh name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; VAR=${4:?variant}; shift 4
D=~/qcc/ws/dqi_int; V=~/qcc/ws/dqi_var; F=$V/snap/route_3e862d11; T=~/Quantum/tools
[ -e $D/pinvar/$N ] && { echo "$D/pinvar/$N exists: give a fresh name"; exit 1; }
[ -e $D/dwx/$N ] && { echo "$D/dwx/$N exists: give a fresh name"; exit 1; }
SN=$D/zsnap/snap_b5r46n; A=$D/zsnap/a3snap_k1j2; R2=~/qcc/ws/r48_dqi_ret2/pinvar/lf3/home/Quantum/tools
# 14th holder: variants l6 (RET2's lf6 = lf3 + P2': its frozen lf6 lateflag copy 2a96b6d3, R48_BETA0=1) and l6r4
# (l6 + my dqi_int_root4, CONV2 row 4); the overlay name then carries the lateflag and root4 shas.
case ${VAR%p} in l6*) R2=~/qcc/ws/r48_dqi_ret2/pinvar/lf6/home/Quantum/tools;; esac
R4=""; case ${VAR%p} in *r4) R4=$(tr -d '\r' < $T/dqi_int_root4.py | sha256sum | cut -c1-8);; esac
( cd $SN && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$SN differs"; exit 1; }
( cd $A && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$A differs"; exit 1; }
( cd $F && sha256sum -c --quiet MANIFEST.sha256 ) || { echo "route differs"; exit 1; }
( cd $D/pin/q8u && sha256sum -c --quiet MANIFEST.sha256 ) || { echo "pin q8u differs"; exit 1; }
sha256sum $T/dqi_var_rinv27.py | grep -q ^6954d6f0 || { echo "rinv27 differs"; exit 1; }
sha256sum $SN/dqi_zn_mkunf5.py | grep -q ^38f96e37 || { echo "mkunf5 differs"; exit 1; }
P17=$(tr -d '\r' < $T/dqi_int_pdr17.py | sha256sum | cut -c1-8); LF8=$(sha256sum $R2/r48_dqi_ret2_lateflag.py | cut -c1-8)
S=$D/zsnap/snap_b5r48_${P17}_${LF8}${R4:+_r4$R4}
if [ ! -d $S ]; then
  rm -rf $S.tmp; mkdir -p $S.tmp; cp $SN/*.py $SN/USAGE.txt $S.tmp/
  cp $R2/r48_dqi_ret2_lateflag.py $R2/r48_dqi_ret2_forms.py $S.tmp/
  tr -d '\r' < $T/dqi_int_pdr17.py > $S.tmp/dqi_int_pdr17.py
  [ -n "$R4" ] && tr -d '\r' < $T/dqi_int_root4.py > $S.tmp/dqi_int_root4.py
  echo "B5 $(date -u +%FT%H:%MZ): snap_b5r46n + RET2's copies from $R2 (lateflag $LF8, forms $(sha256sum $R2/r48_dqi_ret2_forms.py | cut -c1-8)) + dqi_int_pdr17.py $P17${R4:+ + dqi_int_root4.py $R4}" >> $S.tmp/USAGE.txt
  ( cd $S.tmp && sha256sum *.py > SHA256SUMS.txt ) && mv $S.tmp $S
fi
( cd $S && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$S differs"; exit 1; }
SW="R12_FP=a,h DQI_INT_EXT=dqi_zn_pdl2,dqi_zn_kpark,dqi_z2_pdr ZN_PDL_BLOCKS=17,15,13,11,9 ZN_PDL2_FWDONLY=1 ZN_PDL2_SHIFT=16,14,12,10,8 ZN_PDL2_SHIFTPARK=1 ZN_PDL2_YSHARE=1 ZN_PDL2_HAND=7,8,9,10 Z2_PDR_SHARE=1 Z2_PDR_BLOCKS=15,13,11 Z2_PDR_RBLOCKS=15,13,11 Z2_PDR_REV=1 Z2_PDR_LOCAL=1 Z2_PDR_LIVEB=1 DQI_VAR_ROUTE_LOG=1"
EX="R12_FP=h ZN_PDL2_FWDONLY= ZN_PDL2_NPARK=1 DQI_INT_PP2_NOINPUT=1 ZN_PDL2_HAND=7,8,9,10,11,12,13 DQI_INT_SS_OFF=1 R12_SOURCE=$D/pinvar/$N/home/qcc/ws/dqi_zn/src/kv_unfused5"
Q15="ZN_PDL2_LATECONV=1 ZN_PDL2_SHIFTPARK= DQI_INT_SS_OFF= DQI_INT_SS_MAXM=7"
AL="DQI_INT_EXT=dqi_zn_pdl2,dqi_zn_kpark,dqi_zn_alias,dqi_int_aliaslog,dqi_z2_pdr ZN_ALIAS_ROOT=1 ZN_ALIAS_B=3,4,5"
K3="ZN_PDL_BLOCKS=17,15,13,11 ZN_PDL2_SHIFT=16,14,12,10 ZN_ALIAS_C=16,14,12,10 ZN_ALIAS_CL=1 ZN_ALIAS_NPARK=1 DQI_INT_SS_LEVELS=6,7,8,9 DQI_INT_PFDP=1 DQI_INT_PFDT=1 DQI_INT_PFDC=1 DQI_INT_A3_SNAP_OVR=$A DQI_INT_FINH=1 DQI_INT_J6=1 DQI_INT_J6_ZSHARE=1 DQI_INT_J6_Z0SHARE=1"
E=dqi_zn_pdl2,dqi_zn_kpark,dqi_zn_alias,dqi_int_aliaslog,dqi_int_pdr17
X="DQI_INT_EXT=$E,r46_dqi_syn_tail,r46_dqi_ext_fgret,r48_dqi_ret2_lateflag,r46_dqi_conv_zpark R46_ZP_INCOND=1 DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b R46_TAIL_MODE=t4 R46_NOSU2=1 R48_KPOINT=1 R48_N3FLAG=1"
PROF=0; case $VAR in *p) PROF=1; VAR=${VAR%p};; esac
case $VAR in
  l6) B17="R48_BETA0=1";;
  l6r4) B17="R48_BETA0=1"; X="${X/r48_dqi_ret2_lateflag,r46_dqi_conv_zpark/r48_dqi_ret2_lateflag,dqi_int_root4,r46_dqi_conv_zpark}";;
  l3) B17="";;
  f) B17="Z2_PDR_BLOCKS=17,15,13,11";;
  fr) B17="Z2_PDR_BLOCKS=17,15,13,11 Z2_PDR_RBLOCKS=17,15,13,11";;
  *) echo "unknown variant $VAR"; exit 1;;
esac
[ $PROF = 1 ] && X="$X R12_PROF=1"
echo "$N variant $VAR prof $PROF overlay $(basename $S) extra: $B17 $*" > $D/r48b_$N.txt
mkdir -p $D/pinvar $D/logs
tr -d '\r' < $T/sh/dqi_int_pinvar.sh > $D/pinvar.sh
cd $D && ROUTE_FROZEN=$F/dqi_var_route.py OVERLAY=$S KIT_MK=$S/dqi_zn_mkunf5.py KEEP_KMB=1 setsid nohup bash $D/pinvar.sh \
  q8u $N $CC $PC dqi_var_rinv27.py var_rinv27,var_route $SW $EX $Q15 $AL $K3 $X $B17 "$@" > $D/pinvar_$N.txt 2>&1 < /dev/null &
if [ $PROF = 1 ]; then
  B=$D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u; mkdir -p $D/r46/$N
  setsid nohup bash -c "until [ -s $B/quick.log ] || grep -a -q '^end_' $D/pinvar/$N/RESULT.txt 2>/dev/null; do sleep 5; done; \
    ln -f $B/dqi.kmb $D/r46/$N/dqi.kmb; cp $B/segments.json $D/r46/$N/; echo \"marker \$(date -u +%H:%MZ)\" > $D/r46/$N/marker.txt" \
    > /dev/null 2>&1 < /dev/null &
  echo "r48b $N ($VAR, markers) started on $CC; marker file -> $D/r46/$N/dqi.kmb"; exit 0
fi
tr -d '\r' < $T/sh/dqi_int_dwxafter.sh > $D/dwxafter.sh
setsid nohup bash $D/dwxafter.sh $N $CC > /dev/null 2>&1 < /dev/null &
echo "r48b $N ($VAR) started (build on $CC, then dwx), overlay $(basename $S); result line in $D/rebuild_$N.txt"
