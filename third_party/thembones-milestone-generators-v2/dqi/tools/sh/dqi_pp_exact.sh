#!/usr/bin/env bash
# Workstream dqi_pp (P1): EXACT-BY-CONSTRUCTION retained-product post-pass on one RAW emitted DQI build, one command.
# Runs ON A BOX with ~/qcc/.venv_work (numba) and ~/qcc/.venv (verifier), tools under ~/Quantum/tools.
#   bash dqi_pp_exact.sh <raw build dqi.kmb> <fresh out dir> [cores=4-11] [paired batches=200] [extra zero.npz]
# 1 exact value numbering with CCX and CSWAP product forms (tools/dqi_pp_gvn.py)
# 2 exact zero intervals (tools/dqi_pp_zero.py); optionally merged with proved ones (tools/dqi_pp_zprove.py output)
# 3 links: consecutive free CCX / CSWAP whose products (a*b, or d = c*(x+y)) have the same EXACT form, chains over
#   3+ occurrences, cross-type allowed; roving scratch on exact intervals (tools/dqi_pp_ret.py, DQI_PP_EXACTALL=1)
# 4 sanity: every scratch checkpoint zero on 192 sample runs (tools/dqi_pp_zcheck.py); paired-seed test against the
#   input (tools/dqi_pp_paired.py, salt ""): CHILD-ONLY must read 0.
# Nothing sampled is relied on by the rewrite itself. Output: <out>/dqi.kmb, <out>/report.txt.
set -u
IN=$(readlink -f "${1:?raw build}"); OUT=${2:?out dir}; C=${3:-4-11}; NB=${4:-200}; XZ=${5:-}
T=~/Quantum/tools; PY=~/qcc/.venv_work/bin/python; PV=~/qcc/.venv/bin/python
[ -e "$OUT/dqi.kmb" ] && { echo "$OUT/dqi.kmb exists: give a fresh folder"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
P="nice -n 10 taskset -c $C"; export NUMBA_NUM_THREADS=8; t0=$(date +%s); R=report.txt; : > $R
log() { echo "$*" | tee -a $R; }
log "input $(sha256sum $IN | cut -c1-16) $(date -u +%H:%M)Z"
$P $PY $T/dqi_pp_gvn.py $IN gvn.npz 1 > /dev/null
$P $PY $T/dqi_pp_zero.py $IN zero.npz | grep intervals | tee -a $R
$P $PY $T/dqi_pp_pairs2.py $IN pairs.npz 2 5 > /dev/null          # interface only: exact mode ignores fingerprints
ZR=zero.npz
if [ -n "$XZ" ]; then
  $PY -c "
import numpy as np
a = np.load('zero.npz'); b = np.load('$XZ')
np.savez('zero_all.npz', wire=np.r_[a['wire'], b['wire']], start=np.r_[a['start'], b['start']],
         end=np.r_[a['end'], b['end']], structural=np.zeros(len(a['wire']) + len(b['wire']), bool))"
  ZR=zero_all.npz
fi
$PY -c "import numpy as np; np.savez('drop.npz', wire=np.zeros(0, np.int64), time=np.zeros(0, np.int64))"
DQI_PP_EXACTALL=1 DQI_PP_DROP=drop.npz $P $PY $T/dqi_pp_ret.py $IN pairs.npz gvn.npz $ZR dqi.kmb > ret.log 2>&1
grep -E "exact links|admitted|dropped by|hops|chains |static" ret.log | tee -a $R
$P $PY $T/dqi_pp_zcheck.py $IN dqi.kmb.assign.npz drop_it.npz 192 | tail -n 1 | tee -a $R
$P $PV $T/dqi_pp_paired.py dqi.kmb $IN salt0 8 $NB | head -n 3 | tee -a $R
log "output $(sha256sum dqi.kmb | cut -c1-16); runtime $(( $(date +%s) - t0 )) s"
# opt-in diagnostics (after the output is written; dqi.kmb is not touched): links kept with d fewer free wires per
# moment (the pass at a narrower layout) and the wire stretches worth most links (tools/dqi_pp_windows.py)
if [ "${DQI_PP_DIAG:-0}" = "1" ]; then
  $P $PY $T/dqi_pp_capshift.py $IN gvn.npz $ZR 0,2,4,8,12 2>&1 | sed 's/^/diag capacity: /' | tee -a $R
  DQI_PP_WIN_LMAX=1000000 $P $PY $T/dqi_pp_windows.py $IN gvn.npz $ZR 3 20 2>&1 | grep -E "^ ?[0-9]+\." | \
    sed 's/^/diag stretch <= 1M ops: /' | tee -a $R
fi
