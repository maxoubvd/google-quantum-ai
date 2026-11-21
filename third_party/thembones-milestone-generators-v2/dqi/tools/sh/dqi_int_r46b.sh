#!/usr/bin/env bash
# dqi_int_r46b.sh -- B5 (Nov 15 18:30Z, twelfth holder): best_q23j5's recipe (= dqi_int_rebuild_j6j5.sh) + tonight's
# R46 pieces as EXT modules, from frozen copies, into FRESH folders on box 6; then P1's dwx pass (dwxafter.sh).
#   bash dqi_int_r46b.sh <fresh name> <c1,c2> <c1-c2> <variant>          (runs ON box 6, detached)
# variants: k1 = FG (r46_dqi_ext_fgret) + zpark v2 (r46_dqi_conv_zpark, R46_ZP_INCOND=1)
#           k2 = k1 + T4 (r46_dqi_syn_tail mode t4; tailctl taken out of DQI_JUNC_R8P)
#           k3 = k2 + sU2 dropped (R46_NOSU2=1 with R46-DQI-CONV's copies of pdl2 / pdr)   [gated piece]
#           t4 = T4 alone (control); fg = FG alone; ft = FG + T4; z = zpark v2 alone; s = sU2's drop alone
#           (single pieces for the lead's 18:22Z per-piece paired tests)
# Overlays (made once, sha-checked): zsnap/snap_b5r46a = snap_b5j6d + LF copies of the four R46 modules;
# zsnap/snap_b5r46n = snap_b5r46a with r46_dqi_conv_mknosu2.py's copies of dqi_zn_pdl2.py / dqi_z2_pdr.py.
# Published file: ~/qcc/ws/dqi_int/dwx/<name>/final/dqi.kmb; result line in ~/qcc/ws/dqi_int/rebuild_<name>.txt.
set -u
N=${1:?fresh name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; VAR=${4:?variant}; D=~/qcc/ws/dqi_int; T=~/Quantum/tools
[ -e $D/pinvar/$N ] && { echo "$D/pinvar/$N exists: give a fresh name"; exit 1; }
[ -e $D/dwx/$N ] && { echo "$D/dwx/$N exists: give a fresh name"; exit 1; }
S0=$D/zsnap/snap_b5j6d; A=$D/zsnap/a3snap_k1j2; SA=$D/zsnap/snap_b5r46a; SN=$D/zsnap/snap_b5r46n
chk() { sha256sum "$1" 2>/dev/null | grep -q "^$2" || { echo "differs or missing: $1 (want $2)"; exit 1; }; }
( cd $S0 && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$S0 sums differ"; exit 1; }
( cd $A && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$A sums differ"; exit 1; }
chk $D/pin/q8u/MANIFEST.sha256 a21a1ad8d1e2c936
chk ~/qcc/ws/dqi_var/snap/route_3e862d11/MANIFEST.sha256 8aac98c9412ef5a4
if [ ! -d $SA ]; then
  rm -rf $SA.tmp; mkdir -p $SA.tmp; cp $S0/*.py $S0/USAGE.txt $SA.tmp/
  for f in r46_dqi_ext_fgret.py r46_dqi_conv_zpark.py r46_dqi_syn_tail.py r46_dqi_syn_t4.py; do tr -d '\r' < $T/$f > $SA.tmp/$f; done
  chk $SA.tmp/r46_dqi_ext_fgret.py 96e3718c4d37e3ed; chk $SA.tmp/r46_dqi_conv_zpark.py 2350d3c6bbd21454
  chk $SA.tmp/r46_dqi_syn_tail.py a29f40481fbb34f3; chk $SA.tmp/r46_dqi_syn_t4.py f3dc04d467731c05
  echo "B5 $(date -u +%FT%H:%MZ): snap_b5j6d + R46 modules fgret 96e3718c, zpark 2350d3c6, syn_tail a29f4048, syn_t4 f3dc04d4" >> $SA.tmp/USAGE.txt
  ( cd $SA.tmp && sha256sum *.py > SHA256SUMS.txt ) && mv $SA.tmp $SA
fi
if [ ! -d $SN ]; then
  rm -rf $SN.tmp $SN.ovl; mkdir -p $SN.tmp; cp $SA/*.py $SA/USAGE.txt $SN.tmp/
  python3 $T/r46_dqi_conv_mknosu2.py $S0 $SN.ovl > /dev/null || { echo "mknosu2 failed"; exit 1; }
  cp $SN.ovl/dqi_zn_pdl2.py $SN.ovl/dqi_z2_pdr.py $SN.tmp/ && rm -rf $SN.ovl
  chk $SN.tmp/dqi_zn_pdl2.py def46af71a47db83; chk $SN.tmp/dqi_z2_pdr.py 690ce6dd687ca8f0
  echo "B5 $(date -u +%FT%H:%MZ): + r46_dqi_conv_mknosu2.py copies of dqi_zn_pdl2.py (def46af7) / dqi_z2_pdr.py (690ce6dd)" >> $SN.tmp/USAGE.txt
  ( cd $SN.tmp && sha256sum *.py > SHA256SUMS.txt ) && mv $SN.tmp $SN
fi
( cd $SA && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$SA sums differ"; exit 1; }
( cd $SN && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$SN sums differ"; exit 1; }
E=dqi_zn_pdl2,dqi_zn_kpark,dqi_zn_alias,dqi_int_aliaslog,dqi_z2_pdr
T4S="DQI_JUNC_R8P=ipjunc,nofullcache,cachelevels,latem,rinv24b R46_TAIL_MODE=t4"
PROF=0; case $VAR in *p) PROF=1; VAR=${VAR%p};; esac
case $VAR in
  k1) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_ext_fgret,r46_dqi_conv_zpark R46_ZP_INCOND=1";;
  k2) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_syn_tail,r46_dqi_ext_fgret,r46_dqi_conv_zpark R46_ZP_INCOND=1 $T4S";;
  k3) SN_=snap_b5r46n; X="DQI_INT_EXT=$E,r46_dqi_syn_tail,r46_dqi_ext_fgret,r46_dqi_conv_zpark R46_ZP_INCOND=1 $T4S R46_NOSU2=1";;
  t4) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_syn_tail $T4S";;
  fg) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_ext_fgret";;
  ft) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_syn_tail,r46_dqi_ext_fgret $T4S";;
  z) SN_=snap_b5r46a; X="DQI_INT_EXT=$E,r46_dqi_conv_zpark R46_ZP_INCOND=1";;
  s) SN_=snap_b5r46n; X="DQI_INT_EXT=$E R46_NOSU2=1";;
  *) echo "unknown variant $VAR"; exit 1;;
esac
# variant name + 'p' (k1p, k2p, k3p): the same recipe built with the r7prof markers (R12_PROF=1, no P1, no dwx); the
# marker file and its segments.json are hard-linked into ~/qcc/ws/dqi_int/r46/<name>/ as soon as the raw quick check
# is written (the driver deletes its own copy after the profile run). Used to map INV17's sites onto that recipe.
[ $PROF = 1 ] && X="$X R12_PROF=1"
echo "$N variant $VAR prof $PROF overlay $SN_ extra: $X" > $D/r46b_$N.txt
SNAP=$SN_ bash $D/dqi_int_q16.sh $N $CC $PC ZN_PDL_BLOCKS=17,15,13,11 ZN_PDL2_SHIFT=16,14,12,10 \
  ZN_ALIAS_C=16,14,12,10 ZN_ALIAS_CL=1 ZN_ALIAS_NPARK=1 DQI_INT_SS_LEVELS=6,7,8,9 DQI_INT_PFDP=1 DQI_INT_PFDT=1 DQI_INT_PFDC=1 \
  DQI_INT_A3_SNAP_OVR=$A DQI_INT_FINH=1 DQI_INT_J6=1 DQI_INT_J6_ZSHARE=1 DQI_INT_J6_Z0SHARE=1 $X
if [ $PROF = 1 ]; then
  B=$D/pinvar/$N/home/qcc/ws/dqi_int/r12/b_q8u; mkdir -p $D/r46/$N
  setsid nohup bash -c "until [ -s $B/quick.log ] || grep -a -q '^end_' $D/pinvar/$N/RESULT.txt 2>/dev/null; do sleep 5; done; \
    ln -f $B/dqi.kmb $D/r46/$N/dqi.kmb; cp $B/segments.json $D/r46/$N/; echo \"marker \$(date -u +%H:%MZ)\" > $D/r46/$N/marker.txt" \
    > /dev/null 2>&1 < /dev/null &
  echo "r46b $N ($VAR, markers) started on $CC; marker file -> $D/r46/$N/dqi.kmb"; exit 0
fi
tr -d '\r' < $T/sh/dqi_int_dwxafter.sh > $D/dwxafter.sh
setsid nohup bash $D/dwxafter.sh $N $CC > /dev/null 2>&1 < /dev/null &
echo "r46b $N ($VAR) started (build on $CC, then dwx); result line in $D/rebuild_$N.txt"
