#!/bin/bash
# chem_ld2_e21keptjob.sh K CORE: builder C2 (twenty-first holder, Nov 17 13:10Z). Runs ON THE BOX in a kept run folder
# (in/ = pad.npz, proj.json, k<K>.json; shipped there by tools/sh/chem_ld2_e21keptbox2.sh as in/kjob.sh).
# One 30-output loader at 253 = a SEARCHED stage A at 149 (tools/chem_ld2_e21dual.py: four direct sums of pool-16
# blocks, D = 13,14,15 as block 4) + R62-CHEM-BILIN's point stage (r62_chem_bilin_mask.py: class by class, --cls 200,
# the relief --dc --eqw) + the re-pricer's emitter (r58_chem_kept_rp_emit68.py); all R58 / R62 tools read-only. Falls
# back to the in-folder bases (stage A 150: 254). Writes k<K>/angle_k_<K>.kmx and k<K>/RES; releases its core lock.
K=$1; C=$2
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1
PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
trap "rmdir ~/qcc/ws/chem_ld2/locks/c$C 2>/dev/null" EXIT
mkdir -p k$K
V=$(python3 -c "import json, sys; print(json.load(open(sys.argv[1])).get('V', {}).get(sys.argv[2], ''))" in/proj.json $K)
[ -n "$V" ] || { echo "k$K: no V" > k$K/RES; exit 0; }
NREAL=${NREAL:-1219}; PS=${PSEEDS:-12}; CLS=${CLS:-200}
best=""; bestf=""
try() {   # pfile rep D tag
  nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_c2a.py in/pad.npz --k $K --pfile in/$1 --rep $2 --D $3 \
    --out k$K/planA_$4.json > k$K/c2a_$4.log 2>&1 || return 1
  nice -n 10 taskset -c $C $PY -I $T/r62_chem_bilin_mask.py in/pad.npz --planA k$K/planA_$4.json --K 8 --onlyK --seed 0 \
    --pseeds $PS --dist 6,1,1 --V $V --dc --eqw --nreal $NREAL --cls $CLS --plan k$K/plan_$4.json > k$K/mask_$4.log 2>&1 || return 1
  nice -n 10 taskset -c $C $PY -I $T/r58_chem_kept_rp_emit68.py in/pad.npz k$K/plan_$4.json --V $V --nreal $NREAL \
    --out k$K/kmx_$4.kmx > k$K/emit_$4.log 2>&1 || return 1
  t=$(grep -h -o -E "[0-9]+ Toffolis" k$K/emit_$4.log | tail -n 1 | grep -o -E "^[0-9]+")
  [ -n "$t" ] || return 1
  if [ -z "$best" ] || [ "$t" -lt "$best" ]; then best=$t; bestf=$4; cp k$K/kmx_$4.kmx k$K/angle_k_$K.kmx; fi
  return 0
}
# 1. searched bases: stage A 149 on about three sets of four [M G8F k2]
nice -n 10 taskset -c $C $PY -I $T/chem_ld2_e21dual.py in/pad.npz --k $K --D 13,14,15 --n ${DUALN:-1000000} --seed 0 \
  --out in --pfx k${K}x_ > k$K/dual.log 2>&1
for pf in $(cd in && ls k${K}x_*.json 2>/dev/null); do
  try $pf 4 13,14,15 s${pf%.json}
  [ -n "$best" ] && [ "$best" -le 253 ] && break
done
# 1b. twenty-second holder (Nov 17 14:30Z; R62-CHEM-BILIN's 'For the builders': PATCH 2 into C2's kept step): S5 on
# the best stage-A plan (the point products read point 13's helpers, 7 helpers: 252 on a 149-product stage A [M R62,
# V C3]), by C3's tools/sh/chem_c3_pt5.sh (read-only; R62's mask2 / emit2), up to four partition seeds
if [ -n "$best" ] && [ "${S5:-1}" = 1 ]; then
  # 18:25Z: 12 partition seeds (a rewritten table can need seed 2 or 6 of 12 [M C3]; R65 saw 253 twice on k7 with 4)
  CT=$T PY=$PY nice -n 10 taskset -c $C bash $T/sh/chem_c3_pt5.sh in/pad.npz in/proj.json k$K/planA_$bestf.json $K k$K \
    ${bestf}_s5 0 ${S5SEEDS:-12} > k$K/s5.out 2>&1
  if [ -s k$K/${bestf}_s5.kmx ]; then
    t5=$(grep -c -E "^(CCX|CCZ|CSWAP) " k$K/${bestf}_s5.kmx)
    if [ "$t5" -gt 0 ] && [ "$t5" -lt "$best" ]; then
      best=$t5; bestf=${bestf}_s5; cp k$K/$bestf.kmx k$K/angle_k_$K.kmx
      echo "S5: $t5 Toffolis" > k$K/emit_$bestf.log
    fi
  fi
fi
# 2. the in-folder bases (150) when no searched basis emitted 253 or 254
if [ -z "$best" ] || [ "$best" -gt 254 ]; then
  for rd in 0:13,14,15 2:13,14,15 1:13,14,15 3:13,14,15 4:13,14,15 0:4,10,14 2:5,10,14; do
    R=${rd%%:*}; DD=${rd#*:}
    try k$K.json $R $DD r${R}_${DD//,/-}
    [ -n "$best" ] && [ "$best" -le 254 ] && break
  done
fi
if [ -n "$best" ]; then echo "k$K $bestf: $(tail -n 1 k$K/emit_$bestf.log | cut -c1-160)" > k$K/RES
else echo "k$K: none emitted" > k$K/RES; fi
