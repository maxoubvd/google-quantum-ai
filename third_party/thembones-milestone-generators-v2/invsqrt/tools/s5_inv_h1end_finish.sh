#!/usr/bin/env bash
# S5-INV-H1END: after tools/s5_inv_h1end_build.sh NAME: text of the final file, the EXACT all-input proof against the base
# fc19_b312 (tools/r63_inv_ghost_chk.py check BASE NEW 0: x restored, scratch zero, phase equal to the base's on every input,
# outcomes all 0 / all 1 / two random seeds), and the window where the two texts differ (bit names blanked).
#   bash tools/s5_inv_h1end_finish.sh NAME [CORES=11]
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
W=~/qcc/ws/s5/inv_h1end; NAME=$1; CORES=${2:-11}
N=$W/files/${NAME}_w.r91/invsqrt.kmb; B=$W/base/file.kmx; C=$W/check/$NAME
[ -f $N ] || { echo "no final file $N"; exit 1; }
mkdir -p $C
sha256sum $N | cut -c1-16 > $C/sha.txt; echo "final sha $(cat $C/sha.txt)"
(cd $Q/tools && $P -I -c "
import sys; sys.path.insert(0, '<HOME>/qcc/data/kickmix'); import kickmix as km
c = km.Circuit.from_file('$N', format='kmb'); c.to_file('$C/new.kmx', format='kmx'); print('text written, ops', len(c))")
echo "== window diff $(date -u +%H:%M:%SZ)"
(cd $Q/tools && $P -I s5_inv_h1end_ndiff.py $B $C/new.kmx) | tee $C/ndiff.log
echo "== exact proof against the base, ALL inputs $(date -u +%H:%M:%SZ)"
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $B $C/new.kmx 0 2>&1 | tail -8 | cut -c1-260) | tee $C/proof.log
echo "== done $(date -u +%H:%M:%SZ)"
