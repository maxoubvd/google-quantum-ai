#!/bin/bash
# Builder C3 (fourth holder, Nov 17 19:00Z): ON A BOX. Design seat R65-CHEM-WIRED221's owner route for ONE plain wired
# loader K of a table, automated (its 'For C2' steps 1-6; its tools read-only, run through tools/chem_c3_vrun.py so
# that WR7's 'S' wires emit):
#  1. prices (r65_chem_wired221_exact.py --pricesonly K) unless PRICES is given;
#  2. sweep: plan P1 = 6,6:1:3:0,6:1:3:1,7:2:0:2 (seeds 0-15 x tau 2e-3 / 1e-2) and P2 = 6,6:1:3:0,6:1:3:1,7:3:1:2 (the
#     three-kill last group, tau 5e-3): r65_chem_wired221_dc.py, each emitting its best seed (151 = 225 / 148 = 222);
#  3. P1's best seed dumped (r65_chem_wired221_build.py), group 0 annealed for 38 (builder C2's chem_ld2_e19grp.py),
#     plan shifts built (150 = 224), the best shift's last group annealed for 39 and built (149 = 223);
#  4. P2's best seed dumped, group 0 annealed for 38, plan shifts built (147 = 221);
#  5. pick: the emitted candidate with the least toffoli + price / 0.18 (R65's measured refit rule: a product pays
#     under about 0.18 mHa^2 of formula price) -> OUT/best/{angle_k_K.kmx, table_kK.npz}; OUT/cand.txt lists all.
#   CORES="10 11 12 13 14" bash chem_c3_r65own.sh TABLE K OUT [PRICES]      (writes OUT/DONE)
set -u
TB=$1; K=$2; OUT=$3; PRF=${4:-}
C=(${CORES:-10 11 12 13 14}); N=${#C[@]}
PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=${CHEM_HB:-77} OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $OUT; cd $OUT || exit 1; rm -f DONE
log() { echo "$(date -u +%H:%M:%SZ) k$K $*" >> $OUT/r65own.log; }
# Nov 18 12:05Z (seventh holder): VRUN=dead -> design seat R78-CHEM-DEAD's rule (d) launcher (its 'By-product for C3':
# owner units also at low points 0-2 whose block is outside both readers' block sets); default unchanged
VR="$PY -I $T/chem_c3_vrun.py"; [ "${VRUN:-}" = dead ] && VR="$PY -I $T/r78_chem_dead_lmin.py vrun"
P1="6,6:1:3:0,6:1:3:1,7:2:0:2"; O1="0,3,1,2"; P2="6,6:1:3:0,6:1:3:1,7:3:1:2"; O2="0,2,1,3"
GA=${GA_SECS:-90}; G3A=${G3_SECS:-120}; NS=${NSHIFT:-$N}
run1() {   # run1 CORE LOGFILE module args...  (from the tools' parent folder, nice 10, one core)
  local c=$1 lf=$2; shift 2
  (cd $T/.. && nice -n 10 taskset -c $c "$@" > $lf 2>&1)
}
log "start: table $TB, cores ${C[*]}"
# 1. prices
if [ -z "$PRF" ]; then
  run1 ${C[0]} $OUT/pr.log $PY -I $T/r65_chem_wired221_exact.py $TB --bits 1,2,3 --pricesonly $K --out $OUT/pr
  PRF=$OUT/pr/prices_k$K.npz
fi
[ -s "$PRF" ] || { log "prices FAILED: $(tail -n 2 $OUT/pr.log)"; echo fail > DONE; exit 1; }
log "prices $PRF"
# 19:25Z: a RELATION-BEARING owner (24 / 23 generic outputs) takes my owner step's plans and word orders with R65's
# live units through tools/chem_c3_vrun.py's relation mode (relation-keeping dead units + live units on relation-free
# wires); the sweep only (R65's group refinements are written for the 25-output plans)
NG=$(cd $T/.. && taskset -c ${C[0]} $PY -I $T/chem_c3_gencount.py $TB 2>/dev/null | python3 -c "import json,sys; print(json.load(sys.stdin)['gen'][sys.argv[1]])" $K)
log "generic outputs: ${NG:-?}"
if [ "${NG:-25}" != 25 ]; then
  # 20:12Z: R65's PLACED kill plans, shortened to 24 / 23 outputs (vrun appends the relation's dependent wire after the
  # placed ones): G52 k6 147 (39 | 36 | 36 | 36) against 149-150 with my tie plans, k14 (23) 146 [M local, exact]
  case "$NG" in
    24) RPL="6,6:1:3:0,6:1:3:1,6:1:3:1 6,6:1:3:0,6:1:3:1,6:1:2:1 6,6:1:3:0,6:1:2:1,6:1:3:1" ;;
    23) RPL="6,6:1:3:0,6:1:3:1,5:1:3:1 6,6:1:3:0,5:1:3:1,6:1:3:1 6,6:1:3:0,6:1:3:1,5:1:2:1" ;;
    *) log "no plan for $NG generic outputs"; echo none > DONE; exit 0 ;;
  esac
  i=0; pi=0
  for pl in $RPL; do
    for sd in ${SEEDSETS:-0-3 4-7 8-11 12-15}; do
      for oa in 0,2,1,3 0,3,1,2; do
        nm=c${pi}_${sd}_${oa//,/}
        while [ "$(jobs -r | wc -l)" -ge "$N" ]; do sleep 1; done
        run1 ${C[$((i % N))]} $OUT/$nm.log $VR r65_chem_wired221_dc $TB --k $K --plan "$pl" --oauto $oa --avoid 0 \
          --seeds $sd --tau 1e-2 --maxbit 3 --tries 100 --pad --prices $PRF --emit $OUT/$nm &
        i=$((i + 1))
      done
    done
    pi=$((pi + 1))
  done
  wait
fi
# 2. sweep (25 generic outputs)
JOBS=()
if [ "${NG:-25}" = 25 ]; then
  for s in ${SEEDSETS:-0-3 4-7 8-11 12-15}; do for t in 2e-3 1e-2; do JOBS+=("a_${s}_$t|$P1|$O1|$s|$t"); done; done
  for s in ${SEEDSETS:-0-3 4-7 8-11 12-15}; do JOBS+=("b_${s}_5e-3|$P2|$O2|$s|5e-3"); done
fi
i=0
for j in "${JOBS[@]}"; do
  IFS='|' read -r nm pl oa sd tu <<< "$j"
  while [ "$(jobs -r | wc -l)" -ge "$N" ]; do sleep 1; done
  run1 ${C[$((i % N))]} $OUT/$nm.log $VR r65_chem_wired221_dc $TB --k $K --plan "$pl" --oauto $oa --avoid 0 --seeds $sd \
    --tau $tu --maxbit 3 --tries 100 --pad --prices $PRF --emit $OUT/$nm &
  i=$((i + 1))
done
wait
# parse: one line per sweep run: name products seed price toffoli(0 = not emitted)
$PY - $OUT > $OUT/sweep.txt << 'PYEOF'
import glob, os, re, sys
d = sys.argv[1]
for f in sorted(glob.glob(os.path.join(d, "[abc]_*.log"))):
    t = open(f).read()
    m = re.search(r"BEST: (\d+) products, seed (\d+), formula price of the live flips ([-+0-9.e]+)", t)
    e = re.search(r"toffoli (\d+), qubits", t)
    if m:
        print(os.path.basename(f)[:-4], m.group(1), m.group(2), m.group(3), e.group(1) if e else 0)
PYEOF
log "sweep: $(awk '{print $1":"$2"/"$5"@"$4}' $OUT/sweep.txt | tr '\n' ' ')"
best() {   # best PREFIX -> "name seed price products" of the run with the fewest products, then the lowest price
  awk -v p="$1" 'index($1, p) == 1 && $5 > 0' $OUT/sweep.txt | sort -k2,2n -k4,4g | head -n 1 | awk '{print $1, $3, $4, $2}'
}
# g0 / g3 refinement of one dumped seed: refine NAME PLAN OAUTO SEED TAU PRICE WANT0 WANT1 WANT3
refine() {
  local nm=$1 pl=$2 oa=$3 sd=$4 tu=$5 pr=$6 w1=$7 w3=$8 D=$OUT/r_$1
  mkdir -p $D
  local mp=$(python3 -c "import sys; print(max(float(sys.argv[1]) + 0.03, 0.03))" $pr)
  local ARGS="--k $K --plan $pl --oauto $oa --avoid 0 --seeds $sd-$sd --tries 100 --pad --prices $PRF --tau $tu --maxbit 3"
  run1 ${C[0]} $D/dump.log $VR r65_chem_wired221_build $TB $ARGS --maxprice $mp --dump $D/d.json --emit $D/base
  [ -s $D/d.json ] || { log "$nm: dump FAILED: $(tail -n 1 $D/dump.log | cut -c1-150)"; return; }
  log "$nm: $(grep -h '^DUMP' $D/dump.log | cut -c1-160)"
  if [ "$w1" = auto ]; then       # 20:13Z: the dump's groups with the first one lower by one (the anneal's target)
    w1=$(sed -n 's/^DUMP.*groups \[\([0-9, ]*\)\].*/\1/p' $D/dump.log | head -n 1 | tr -d ' ' | awk -F, -v OFS=, '{$1 = $1 - 1; print}')
    log "$nm: want after the group-0 anneal: $w1"
  fi
  local r=1 c
  for c in "${C[@]}"; do
    (cd $T/.. && nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $D/d.json --g 0 --P 38 --secs $GA --multi \
       --rseed $((20 + r)) --save $D/g0sp_$r.json > $D/g0an_$r.log 2>&1) &
    r=$((r + 1))
  done
  wait
  $PY - $D << 'PYEOF'
import glob, json, sys
d = sys.argv[1]
out = []
for f in sorted(glob.glob(d + "/g0sp_*.json")):
    try:
        x = json.load(open(f))
    except Exception:
        continue
    out += x if isinstance(x, list) else [x]
n = len(out)
out = out[::max(1, n // 24)][:24]
json.dump(out, open(d + "/g0spaces.json", "w"))
print(n)
PYEOF
  local nsp=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" $D/g0spaces.json 2>/dev/null || echo 0)
  log "$nm: group 0 spaces at 38: $nsp kept"
  if [ "${nsp:-0}" -eq 0 ]; then
    # Nov 18 09:40Z (sixth holder): THE LAST-GROUP FIX AS DEFAULT. Group 0 has no space at 38 here, so the last-group anneal
    # below never ran; run it on this dump directly (tools/sh/chem_c3_lastg.sh; k8, k10, k12 of t10_g10k1: -1 each [M, m2])
    if [ "$w3" != "-" ] && [ "${LASTG:-1}" = 1 ]; then
      local gw=$(sed -n 's/^DUMP.*groups \[\([0-9, ]*\)\].*/\1/p' $D/dump.log | head -n 1 | tr -d ' ')
      local gl=${gw##*,}
      local wl=$(echo $gw | awk -F, -v OFS=, '{$NF = $NF - 1; print}')
      if [ -n "$gw" ] && [ "${gl:-0}" -ge 38 ]; then
        log "$nm: last-group fix on the dump (groups $gw -> $wl)"
        TB=$TB K=$K PRF=$PRF DUMP=$D/d.json PLAN=$pl OAUTO=$oa SEED=$sd TAU=$tu WANT=$wl P3=$((gl - 1)) PR3=$((gl - 3)) \
          MAXP=$mp CORES="${C[*]}" OUT=$D/lg bash $T/sh/chem_c3_lastg.sh
        log "$nm: last-group fix: $(cat $D/lg/RESULT 2>/dev/null)"
      fi
    fi
    # 09:58Z (sixth holder), C4OWN=1: builder C4's first-group screen on this dump's placement (another six wires at 38) +
    # the last group one lower (tools/sh/chem_c3_c4own.sh; C4 measured k10 of t10_g10k1 at 224 by the screen alone)
    if [ "$w3" != "-" ] && [ "${C4OWN:-0}" = 1 ]; then
      TB=$TB K=$K PRF=$PRF DD=$D PLAN=$pl SEED=$sd TAU=$tu MAXP=$mp CORES="${C[*]}" OUT=$D/c4o bash $T/sh/chem_c3_c4own.sh
      log "$nm: C4's first-group screen: $(cat $D/c4o/RESULT 2>/dev/null)"
    fi
    return
  fi
  run1 ${C[0]} $D/plan0.log $PY -I $T/chem_ld2_e19grp.py plan $D/d.json $D/g0spaces.json --g 0 --tries 1500 --out $D/plans_g0.json
  [ -s $D/plans_g0.json ] || { log "$nm: no group-0 plans"; return; }
  local i=0
  while [ $i -lt $NS ]; do
    mkdir -p $D/s2_$i
    run1 ${C[$((i % N))]} $D/s2_$i/build.log $VR r65_chem_wired221_build $TB $ARGS --want $w1 --maxprice $mp \
      --inj "0=$D/plans_g0.json" --pshift "0=$((i * 53))" --dump $D/s2_$i/d.json --emit $D/s2_$i/base &
    i=$((i + 1)); [ $((i % N)) = 0 ] && wait
  done
  wait
  local I=$(for j in $(seq 0 $((NS - 1))); do t=$(grep -h -o -E 'toffoli [0-9]+' $D/s2_$j/build.log 2>/dev/null | awk '{print $2}'); [ -n "$t" ] && echo "$t $j"; done | sort -n | head -n 1 | awk '{print $2}')
  log "$nm: shifts: $(for j in $(seq 0 $((NS - 1))); do echo -n "$j:$(grep -h -o -E 'toffoli [0-9]+' $D/s2_$j/build.log 2>/dev/null | awk '{print $2}') "; done)"
  [ -n "$I" ] || return
  [ "$w3" = "-" ] && return
  local S=$D/s2_$I; r=1
  for c in "${C[@]}"; do
    (cd $T/.. && nice -n 10 taskset -c $c $PY -I $T/chem_ld2_e19grp.py anneal $S/d.json --g 3 --P 39 --price 37 --secs $G3A --multi \
       --rseed $r --save $S/an_$r.json > $S/an_$r.log 2>&1) &
    r=$((r + 1))
  done
  wait
  $PY - $S << 'PYEOF'
import glob, json, sys
d = sys.argv[1]
out = []
for f in sorted(glob.glob(d + "/an_*.json")):
    try:
        x = json.load(open(f))
    except Exception:
        continue
    out += x if isinstance(x, list) else [x]
n = len(out)
out = out[::max(1, n // 30)][:30]
json.dump(out, open(d + "/spaces.json", "w"))
PYEOF
  local n3=$(python3 -c "import json,sys; print(len(json.load(open(sys.argv[1]))))" $S/spaces.json 2>/dev/null || echo 0)
  log "$nm: last-group spaces at 39 (shift $I): $n3"
  [ "${n3:-0}" -gt 0 ] || return
  run1 ${C[0]} $S/plan3.log $PY -I $T/chem_ld2_e19grp.py plan $S/d.json $S/spaces.json --g 3 --tries 3000 --out $S/plans_g3.json
  [ -s $S/plans_g3.json ] || return
  mkdir -p $D/s3_$I
  run1 ${C[0]} $D/s3_$I/build.log $VR r65_chem_wired221_build $TB $ARGS --want $w3 --maxprice $mp \
    --inj "0=$D/plans_g0.json;3=$S/plans_g3.json" --pshift "0=$((I * 53))" --emit $D/s3_$I
  log "$nm: last group: $(grep -h -o -E 'toffoli [0-9]+' $D/s3_$I/build.log | tail -n 1)"
}
# 20:45Z: the refinements gain at most about 3 on the sweep's best [M]: skipped when the sweep's best is CNT + 3 or more
# (CNT = the entry's count of this loader, given by chem_c3_r65box.sh)
SB=$(awk '$5 > 0 {print $5}' $OUT/sweep.txt | sort -n | head -n 1)
if [ -n "${CNT:-}" ] && [ -n "$SB" ] && [ "$SB" -ge $((CNT + 3)) ]; then
  log "sweep best $SB against the entry's $CNT: refinements skipped"
  REFINE=0
fi
[ "${REFINE:-1}" = 1 ] && { set -- $(best a_); [ $# -ge 4 ] && { A_NM=$1; A_SD=$2; A_PR=$3; A_TU=$(echo $1 | awk -F_ '{print $3}'); refine $A_NM "$P1" $O1 $A_SD $A_TU $A_PR 38,36,36,40 38,36,36,39; }; }
[ "${REFINE:-1}" = 1 ] && { set -- $(best b_); [ $# -ge 4 ] && { B_NM=$1; B_SD=$2; B_PR=$3; refine $B_NM "$P2" $O2 $B_SD 5e-3 $B_PR 38,36,36,37 -; }; }
# 20:13Z: the relation-bearing owner's best placed kill plan, group 0 annealed (no last-group step)
set -- ; [ "${REFINE:-1}" = 1 ] && set -- $(best c_)
if [ $# -ge 4 ]; then
  C_NM=$1; C_SD=$2; C_PR=$3
  CN=$(echo $C_NM | sed 's/^c\([0-9]*\)_.*/\1/'); PLC=$(echo $RPL | awk -v n=$((CN + 1)) '{print $n}')
  OAC=$(echo $C_NM | awk -F_ '{print $3}' | sed 's/./&,/g; s/,$//')
  refine $C_NM "$PLC" $OAC $C_SD 1e-2 $C_PR auto -
fi
# 5. pick
$PY - $OUT $K > $OUT/cand.txt << 'PYEOF'
import glob, os, re, sys
d, k = sys.argv[1], sys.argv[2]
for f in sorted(glob.glob(os.path.join(d, "**", f"angle_k_{k}.kmx"), recursive=True)):
    dd = os.path.dirname(f)
    if dd.endswith("/best") or not os.path.exists(os.path.join(dd, f"table_k{k}.npz")):
        continue
    logs = [dd + ".log", os.path.join(dd, "build.log"), os.path.join(os.path.dirname(dd), "build.log"),
            os.path.join(os.path.dirname(dd), "dump.log")]
    for lf in logs:
        if not os.path.exists(lf):
            continue
        t = open(lf).read()
        e = re.findall(r"emitted (\S+)/angle_k_%s\.kmx: toffoli (\d+)" % k, t)
        e = [x for x in e if os.path.abspath(x[0]) == os.path.abspath(dd)]
        p = re.findall(r"formula price (?:of the live flips )?([-+0-9.e]+)", t)
        if e and p:
            print(e[-1][1], p[-1], dd)
            break
PYEOF
SEL=$(awk '{printf "%.4f %s %s %s\n", $1 + $2 / 0.18, $1, $2, $3}' $OUT/cand.txt | sort -k1,1g | head -n 1)
log "candidates: $(awk '{print $1"@"$2}' $OUT/cand.txt | tr '\n' ' ')"
if [ -n "$SEL" ]; then
  set -- $SEL; mkdir -p $OUT/best; cp $4/angle_k_$K.kmx $4/table_k$K.npz $OUT/best/
  echo "$2 $3 $4" > $OUT/best/PICK; log "PICK toffoli $2 price $3 from $4"
else
  log "no candidate emitted"
fi
echo ok > DONE
