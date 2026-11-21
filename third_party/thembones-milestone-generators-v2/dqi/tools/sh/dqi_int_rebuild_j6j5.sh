#!/usr/bin/env bash
# dqi_int_rebuild_j6j5.sh -- B5 (Nov 15): rebuild the J = 6 (B5-3) + J2 DQI file (pair-program delta hook skipped: no delta program) from frozen pieces, on box 6 or box 1, detached,
# into FRESH folders (refuses existing ones):
#   bash dqi_int_rebuild_j6j5.sh <fresh name> <c1,c2> <c1-c2>
# 1. pinned build ~/qcc/ws/dqi_int/pinvar/<name>/ (best_q22ch's recipe + DQI_INT_J6=1 DQI_INT_J6_ZSHARE=1 DQI_INT_J6_Z0SHARE=1 and K1's a3
#    snapshot with R40-DQI-RETURN's patch J2); 2. P1's post-pass tools/sh/dqi_pp_dwx.sh (dead-write -> links ->
#    dead-write) on its raw file into ~/qcc/ws/dqi_int/dwx/<name>/; the published file is dwx/<name>/final/dqi.kmb.
# Needs (each checked by hash): ~/qcc/ws/dqi_int/pin/q8u (MANIFEST a21a1ad8), B5's overlay zsnap/snap_b5j6d (= snap_b5ch + the R40-DQI-RETURN re-pricer's copy of
# dqi_pub_pair.py cb942adb (tools/r40_dqi_return_rp_mkpair.py on pin q8u's 7228f649) +
# tools/dqi_int_pfdp.py 341df422 as dqi_pub_pfd.py + dqi_int_j6blk.py 568067c8 + dqi_int_j6sched.py f48ec553 +
# r37_dqi_ext_chain.py b932ba46 + tools/dqi_int_zt6.py 1c434105 as dqi_pub_zt.py), K1's a3 snapshot zsnap/a3snap_k1j2
# (= tools/dqi_k1_a3snap_k1j2: finish 4954f174, a3 f08dda2d, kern_div d53547a8), B7's route_3e862d11 (MANIFEST 8aac98c9).
set -u
N=${1:?fresh name}; CC=${2:?c1,c2}; PC=${3:?c1-c2}; D=~/qcc/ws/dqi_int
[ -e $D/pinvar/$N ] && { echo "$D/pinvar/$N exists: give a fresh name"; exit 1; }
[ -e $D/dwx/$N ] && { echo "$D/dwx/$N exists: give a fresh name"; exit 1; }
S=$D/zsnap/snap_b5j6d; A=$D/zsnap/a3snap_k1j2
chk() { sha256sum "$1" 2>/dev/null | grep -q "^$2" || { echo "differs or missing: $1 (want $2)"; exit 1; }; }
chk $S/dqi_pub_pfd.py 341df422deac6e89
chk $S/dqi_int_j6blk.py 568067c8d3a0f7af
chk $S/dqi_int_j6sched.py f48ec55373bc4be0
chk $S/dqi_pub_pair.py cb942adbb493db59
chk $S/r37_dqi_ext_chain.py b932ba467fdf7cfd
chk $S/dqi_pub_zt.py 1c43410582822386
chk $S/dqi_int_shiftscale.py df1295159d120a43
chk $S/dqi_zn_alias.py 92cd29f7296db78a
chk $S/dqi_int_a3.py 318849fc1026de90
( cd $S && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$S sums differ"; exit 1; }
( cd $A && sha256sum -c --quiet SHA256SUMS.txt ) || { echo "$A sums differ"; exit 1; }
chk $A/dqi_pub_finish.py 4954f174
chk $A/dqi_pub_a3.py f08dda2d
chk $D/pin/q8u/MANIFEST.sha256 a21a1ad8d1e2c936
chk ~/qcc/ws/dqi_var/snap/route_3e862d11/MANIFEST.sha256 8aac98c9412ef5a4
SNAP=snap_b5j6d bash $D/dqi_int_q16.sh $N $CC $PC ZN_PDL_BLOCKS=17,15,13,11 ZN_PDL2_SHIFT=16,14,12,10 \
  ZN_ALIAS_C=16,14,12,10 ZN_ALIAS_CL=1 ZN_ALIAS_NPARK=1 DQI_INT_SS_LEVELS=6,7,8,9 DQI_INT_PFDP=1 DQI_INT_PFDT=1 DQI_INT_PFDC=1 \
  DQI_INT_A3_SNAP_OVR=$A DQI_INT_FINH=1 DQI_INT_J6=1 DQI_INT_J6_ZSHARE=1 DQI_INT_J6_Z0SHARE=1
tr -d '\r' < ~/Quantum/tools/sh/dqi_int_dwxafter.sh > $D/dwxafter.sh
setsid nohup bash $D/dwxafter.sh $N $CC > /dev/null 2>&1 < /dev/null &
echo "rebuild $N started (build on $CC, then dwx); result line in $D/rebuild_$N.txt"
