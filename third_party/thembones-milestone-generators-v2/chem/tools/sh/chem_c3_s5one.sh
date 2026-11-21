#!/bin/bash
# Builder C3 (seventh holder, Nov 18 18:25Z): ON A BOX. ONE 30-output loader K rebuilt on a given table, the chain's way
# (chem_c3_chainbox.sh s30: r58_chem_kept_vpad.py -> chem_c3_bsrun.sh (S5, the fallback) -> R66's S5T through
# chem_c3_s5t.sh on the bsrun's stage-A plan). Circuits in OUT/b<K>/ (S5T: kmx_k<K>t_l*.kmx), log OUT/s5one.log, OUT/DONE.
#   S5T_SECS=2400 CORES="24 25 26 27" bash chem_c3_s5one.sh TABLE K OUT
set -u
TB=$1; K=$2; R=$3; CO=(${CORES:-24 25 26 27})
P=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
mkdir -p $R/vin; cd $R || exit 1; rm -f DONE
log() { echo "$(date -u +%H:%M:%SZ) k$K $*" >> $R/s5one.log; }
$P -I $T/r58_chem_kept_vpad.py $TB --ks $K --out vin/pad.npz > vin/vpad.log 2>&1
python3 -c "import json, sys; json.dump({'V': json.load(open(sys.argv[1]))['V']}, open(sys.argv[2], 'w'))" vin/pad.npz.json vin/proj.json \
  || { log "vpad failed"; echo fail > DONE; exit 1; }
COMBOS="0 1|1 1|2 1" CT=$T PY=$P nice -n 10 taskset -c ${CO[0]} bash $T/sh/chem_c3_bsrun.sh $R/vin $K $R/b$K 3 > $R/b$K.out 2>&1
TG=$(grep -h 'rc=0 .*_s5' $R/b$K/bsrun_k$K.txt 2>/dev/null | tail -n 1 | sed -n 's/.*rc=0 \(.*\)_s5[p0-9]*: EMITTED.*/\1/p')
log "S5: $(grep -h 'rc=0' $R/b$K/bsrun_k$K.txt 2>/dev/null | tail -n 1 | cut -c1-100)"
[ -n "$TG" ] && [ -s $R/b$K/planA_$TG.json ] || { log "no S5 plan"; echo "ok s5only" > DONE; exit 0; }
CT=$T PY=$P bash $T/sh/chem_c3_s5t.sh $R/vin/pad.npz $R/vin/proj.json $R/b$K/planA_$TG.json $K $R/b$K k${K}t "${CO[*]}" > /dev/null 2>&1
log "S5T: $(grep -h EMITTED $R/b$K/s5t_k$K.txt 2>/dev/null | head -n 1 | cut -c1-90)"
echo "ok" > DONE
