#!/bin/bash
# T10 seat 3 (ON a box, in ~/qcc/ws/chem_t10/pfin): THE TABLE STEP OF A RELATION TABLE UNDER R53's CARRIER RULE:
#   1 chem_t10_obcar.py (table7c with the 76 one-body rows rounded IN CARRIER FORM under the relations)  -> D/table.npz
#   2 r48_chem_q7c_obsigned.py (signed range: a gauge move)                                         -> D/table_s.npz
#   3 r48_chem_q7c_vtab.py (the hosts' evaluation)                                                  -> D/vtab_s.json
#   4 r48_chem_q7c_refit.py (all-gates refit of the weights, targets 19.9999)                       -> D/all.npz
#   5 C2's carrier gauge chem_ld2_e11car7.py (unchanged)                                            -> D/all_c.npz
#   6 chem_t10_pcheck.py: carriers on 75 columns, every relation on the STORED words (fitted targets) -> D/pcheck.json
#   SRC=<rounded.npz> D=<dir> REL=<layout for the one-body rows> RELG=<rel_used_g.txt> [OBT=24 OBW=12 OBF=3 T7ARGS=--quick] bash chem_t10_pfin.sh
cd ${W:-~/qcc/ws/chem_t10/pfin} || exit 1
: ${SRC:?} ${D:?} ${REL:?} ${RELG:?}
T=${T:-~/Quantum/tools}
if [ -z "$PY" ]; then
  for p in ~/qcc/.venv_work/bin/python ~/qcc/.venv/bin/python; do [ -x $p ] && { PY=$p; break; }; done
fi
export OMP_NUM_THREADS=${FTH:-1} OPENBLAS_NUM_THREADS=${FTH:-1} MKL_NUM_THREADS=${FTH:-1} CUDA_VISIBLE_DEVICES= JAX_PLATFORMS=cpu
mkdir -p $D
date -u +%H:%M:%SZ
T10_OB_RELAX=${OBRELAX:-0} $PY $T/chem_t10_obcar.py $SRC $D --rel @$REL --ob-tries ${OBT:-24} --ob-workers ${OBW:-12} --ob-fine ${OBF:-3} --ob-mode ${OBMODE:-ge} --brot 7 --rel-out $D/rel_ob.txt ${T7ARGS:-} > $D/t7.log 2>&1; echo "table7c (obcar) exit $?"
grep "\[table7" $D/t7.log | cut -c1-330 | tail -n 8
[ -s $D/table.npz ] || exit 2
if [ "${OBMODE:-ge}" = ge ]; then
  $PY $T/chem_t10_relmerge.py $RELG $D/rel_ob.txt $D/rel_final.txt && RELG=$D/rel_final.txt
fi
date -u +%H:%M:%SZ
$PY $T/r48_chem_q7c_obsigned.py $D/table.npz $D/table_s.npz --brot 7 | cut -c1-250 || exit 3
$PY $T/r48_chem_q7c_vtab.py $D/table_s.npz --brot 7 --loader ${LOADER:-3550} --json $D/vtab_s.json 2>&1 | tail -n 3 | cut -c1-300
date -u +%H:%M:%SZ
$PY $T/r48_chem_q7c_refit.py $D/table_s.npz $D/all.npz --brot 7 --start table --t2 ${RT2:-19.9999} --t1 ${RT1:-19.9999} --margin 1e-6 ${RFARGS:-} | grep -E 'BEST|NO FEAS|saved'
[ -s $D/all.npz ] || exit 4
$PY $T/chem_ld2_e11car7.py $D/all.npz $D/all_c.npz --brot 7 | tail -n 2
$PY $T/chem_t10_pcheck.py $D/all_c.npz --rel @$RELG --json $D/pcheck.json | tail -n 12
NB=$($PY -c "import json; d = json.load(open('$D/pcheck.json')); print(d['broken_one_body'] + d['broken_two_body'])")
if [ "${NB:-1}" != 0 ]; then     # relations a late pole broke: the cheapest exact repair on the stored words, then the refit
  date -u +%H:%M:%SZ; echo "repair: $NB broken (row, relation) pairs"
  $PY $T/chem_t10_prepair.py $D/all_c.npz $D/all_cr0.npz --rel @$RELG | tail -n 3
  $PY $T/r48_chem_q7c_refit.py $D/all_cr0.npz $D/all_r.npz --brot 7 --start table --t2 ${RT2:-19.9999} --t1 ${RT1:-19.9999} --margin 1e-6 | grep -E 'BEST|NO FEAS|saved'
  [ -s $D/all_r.npz ] && $PY $T/chem_ld2_e11car7.py $D/all_r.npz $D/all_c.npz --brot 7 | tail -n 2
  $PY $T/chem_t10_pcheck.py $D/all_c.npz --rel @$RELG --json $D/pcheck.json | tail -n 6
fi
$PY $T/r48_chem_q7c_vtab.py $D/all_c.npz --brot 7 --loader ${LOADER:-3550} --json $D/vall_c.json 2>&1 | tail -n 2 | cut -c1-300
cp $RELG $D/rel_used_g.txt 2>/dev/null
date -u +%H:%M:%SZ
