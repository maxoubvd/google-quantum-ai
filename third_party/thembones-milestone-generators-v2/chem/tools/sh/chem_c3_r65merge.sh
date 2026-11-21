#!/bin/bash
# Builder C3 (fifth holder, Nov 18 04:30Z): chem_c3_r65box.sh WITHOUT the R65 route step and without S5T+: the choice runs on
# k<K>/cand.txt files POOLED from earlier runs (candidate dirs copied into this run folder), then the transplant, the
# neighbours (g42, relrun, S5T on 30-output) and the refit chain exactly as there. Original header:
# Builder C3 (fourth holder, Nov 17 19:00Z): ON A BOX. The second stage on a verified entry's words (RUNDIR/w1.npz):
#  1. design seat R65-CHEM-WIRED221's owner route (tools/sh/chem_c3_r65own.sh) on each loader of KS, GS cores a loader;
#  2. the owners whose pick beats the entry's count (RUNDIR/counts.txt 'K toffoli') transplanted
#     (tools/chem_c3_transplant.py: batches K-1, K copied whole) -> final.npz;
#  3. the weights refit with every word fixed (tools/r48_chem_q7c_refit.py, R65's last step) -> refit.npz, one core,
#     while the neighbours build: wired (25 generic) by R59's g42 driver, 30-output by vpad + chem_c3_bsrun.sh (S5) +
#     R66-CHEM-PLAN244's S5T, relation-bearing (24 / 23) by tools/sh/chem_c3_relrun.sh.
#   KS="8 10 12 14" CORES="10 ... 19" GS=5 bash chem_c3_r65box.sh RUNDIR      (writes DONE)
set -u
cd "$1" || exit 1
R=$(pwd); CO=(${CORES:-10 11 12 13 14 15 16 17 18 19}); GS=${GS:-5}
PY=~/qcc/.venv_work/bin/python; T=~/Quantum/tools
export CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
log() { echo "$(date -u +%H:%M:%SZ) $*" >> run.log; }
rm -f DONE
# 19:58Z: R66-CHEM-PLAN244's S5T+ on k1 (tools/sh/chem_c3_s5tp.sh: batches 0, 1 are never an owner's) on the LAST three
# cores while R65's route runs on the others, when k1 is a kept loader and k2 a plain wired one; grafted before the refit
$PY -I $T/chem_c3_gencount.py w1.npz > gen1.json 2> /dev/null
g1() { python3 -c "import json,sys; print(json.load(open('gen1.json'))['gen'].get(sys.argv[1], 0))" "$1"; }
S5=0
case " ${OWNERS:-} ${KS:-} " in *" 2 "*) S5TP=0;; esac          # batch 1 is k2's when k2 owns
if [ "${S5TP:-0}" = 1 ] && [ "$(g1 1)" -ge 30 ] && [ "$(g1 2)" = 25 ] && [ ${#CO[@]} -ge 8 ]; then
  WL=$(python3 -c "import json; g=json.load(open('gen1.json'))['gen']; print(','.join(k for k in sorted(g, key=int) if 0 < int(g[k]) <= 25))")
  S5C="${CO[*]: -4}"; CO=(${CO[*]:0:$((${#CO[@]} - 4))}); S5=1          # 20:59Z: four processes
  rm -rf $R/s5tp
  ( CORES="$S5C" setsid nohup bash $T/sh/chem_c3_s5tp.sh $R $R/w1.npz 1 "$WL" > /dev/null 2>&1 < /dev/null & )
  log "S5T+ on k1 started, detached (cores $S5C, wired $WL)"
fi
[ "$S5" = 1 ] && GS=$(( ${#CO[@]} / 2 ))                      # two R65 groups on what S5T+ leaves
NG=$(( ${#CO[@]} / GS )); [ $NG -ge 1 ] || NG=1
log "start: loaders ${KS:-}, ${#CO[@]} cores in $NG groups of $GS"
# 19:20Z: the choice over all loaders together: per loader the candidate with the least toffoli + price / 0.18 (R65's
# rule), only where it beats the entry's count; then, while the summed formula price is over BUDGET mHa^2 (the refit
# absorbed 0.79 in the rehearsal [M]), one loader steps to its next cheaper-price candidate, the step that gives up the
# fewest Toffolis per mHa^2 freed first
$PY - "${KS:-}" "${BUDGET:-1.2}" > choice.txt << 'PYEOF'
import os, sys
ks, budget = [int(x) for x in sys.argv[1].split()], float(sys.argv[2])
cnt = {int(a): int(b) for a, b in (l.split()[:2] for l in open("counts.txt") if l.strip())}
opts = {}
for k in ks:
    f = f"k{k}/cand.txt"
    if not os.path.exists(f) or k not in cnt:
        continue
    c = []
    for l in open(f):
        t, p, d = l.split()[:3]
        # 21:24Z: also the same count at a NEGATIVE price (on a refit table many flips lower the error: the next refit
        # then lowers lambda_eff [M deep run on c3_g9s2rh65k1: prices -0.01 to -0.11 mHa^2])
        # 21:56Z: REVERTED to fewer Toffolis only: the same-count re-owns of k8 / k14 at -0.147 mHa^2 RAISED lambda_eff
        # after the refit (46.934728 -> 46.941249, c3_g9s2deep2 [V]): a formula price at fixed weights is not the refit's price
        if int(t) < cnt[k]:
            c.append((int(t), float(p), d))
    if c:
        # keep the Pareto front (fewer toffoli or lower price), ordered by toffoli + price / 0.18
        c.sort(key=lambda x: (x[0] + x[1] / 0.18, x[1]))
        opts[k] = c
pick = {k: 0 for k in opts}
def total():
    return sum(opts[k][pick[k]][1] for k in opts)
while opts and total() > budget:
    best = None
    for k in opts:
        cur = opts[k][pick[k]]
        for j in range(len(opts[k])):
            alt = opts[k][j]
            if alt[1] < cur[1]:
                r = (alt[0] - cur[0]) / (cur[1] - alt[1])
                if best is None or r < best[0]:
                    best = (r, k, j)
    if best is None:
        break
    pick[best[1]] = best[2]
for k in sorted(opts):
    t, p, d = opts[k][pick[k]]
    print(k, t, p, d)
print("TOTAL", sum(opts[k][pick[k]][1] for k in opts), file=sys.stderr)
PYEOF
TP=""; OWN=""
while read -r K t p d; do
  log "k$K: R65 choice $t (price $p) against the entry's $(awk -v k=$K '$1 == k {print $2}' counts.txt) ($d)"
  mkdir -p k$K/use; cp $d/angle_k_$K.kmx $d/table_k$K.npz k$K/use/
  TP="$TP $K=k$K/use/table_k$K.npz"; OWN="$OWN $K"
done < choice.txt
if [ -n "$TP" ]; then
  $PY -I $T/chem_c3_transplant.py --allowdead w1.npz final.npz $TP > tp.log 2>&1 || { log "transplant FAILED: $(tail -n 2 tp.log)"; echo fail > DONE; exit 1; }
  log "transplant: $(tail -n 2 tp.log | tr '\n' ' ' | cut -c1-300)"
  RB=$(sed -n 's/.*loaders whose delta changed: \[\(.*\)\]; written.*/\1/p' tp.log | tr -d ',')
else
  log "no R65 owner beats the entry"; cp w1.npz final.npz; RB=""
fi
NB=""; for K in $RB; do case " $OWN " in *" $K "*) ;; *) NB="$NB $K";; esac; done
# 21:02Z: the refit and the neighbours start AT ONCE on a copy (final_nb.npz); S5T+'s graft (dead words of batch 1 only:
# neither the refit's weights nor any neighbour but k2 depend on them) goes into final.npz and refit.npz at the end
[ -n "$TP" ] || [ "$S5" = 1 ] || { echo none > DONE; exit 0; }
cp final.npz final_nb.npz
# 00:12Z: the refit as a RESTART CHAIN (the refit is path-dependent: c3_g9s2rh65k1's weights 46.934728 -> 46.933103 ->
# 46.929584 -> 46.922016 from its own outputs [V c3_g9s2rh65k1n]): R65's targets, then 19.9995 / 19.998 twice from the last
# output, then 19.9995 / 19.9995; refit.npz = the feasible step with the lowest lambda_eff (REFITCHAIN=0: one step)
( c=${CO[0]}
  nice -n 10 taskset -c $c $PY -I $T/r48_chem_q7c_refit.py final_nb.npz refit_0.npz --brot 7 --t2 19.9989 --t1 19.99 --stages 1 > refit_0.log 2>&1
  src=refit_0.npz; i=1
  if [ "${REFITCHAIN:-1}" = 1 ]; then
    for tt in "19.9995 19.998" "19.9995 19.998" "19.9995 19.9995"; do
      set -- $tt
      [ -s $src ] || break
      nice -n 10 taskset -c $c $PY -I $T/r48_chem_q7c_refit.py $src refit_$i.npz --brot 7 --t2 $1 --t1 $2 --stages 1 > refit_$i.log 2>&1
      [ -s refit_$i.npz ] && src=refit_$i.npz; i=$((i + 1))
    done
  fi
  best=$(for f in refit_*.log; do l=$(grep -h BEST $f 2>/dev/null | tail -n 1 | sed -n 's/.*lambda_eff": \([0-9.]*\).*/\1/p'); \
         [ -n "$l" ] && [ -s ${f%.log}.npz ] && echo "$l ${f%.log}.npz"; done | sort -g | head -n 1)
  [ -n "$best" ] && cp ${best#* } refit.npz
  cat refit_*.log > refit.log 2>/dev/null
  log "refit chain: $(for f in refit_*.log; do echo -n "${f%.log} $(grep -h BEST $f 2>/dev/null | tail -n 1 | sed -n 's/.*lambda_eff": \([0-9.]*\).*/\1/p'); "; done) -> ${best:-none}" ) &
$PY -I $T/chem_c3_gencount.py final_nb.npz > genf.json 2> /dev/null
gen() { python3 -c "import json,sys; print(json.load(open('genf.json'))['gen'].get(sys.argv[1], 0))" "$1"; }
W=""; P=""; L=""
for K in $NB; do
  [ "$K" -ge 1 ] || continue
  n=$(gen $K)
  if [ "$n" -ge 30 ]; then P="$P $K"; elif [ "$n" -eq 25 ]; then W="$W $K"; else L="$L $K"; fi
done
log "owners$OWN; neighbours: 30-output$P | wired$W | relation-bearing$L"
i=1
for K in $W; do
  c1=${CO[$((i % ${#CO[@]}))]}; c2=${CO[$(((i + 1) % ${#CO[@]}))]}; i=$((i + 2))
  ( cd $T/.. && nice -n 10 $PY -I tools/r59_chem_groups_g42.py --table $R/final_nb.npz --k $K --seeds 0-63 --cores $c1,$c2 \
      --rounds 6 --secs 600 --out $R/n$K > $R/n$K.log 2>&1 ) &
done
for K in $L; do
  c=${CO[$((i % ${#CO[@]}))]}; i=$((i + 1))
  JOB=$(python3 -c "import json,sys; j=[x for x in json.load(open('genf.json'))['jobs'] if x['k']==int(sys.argv[1])]; j=j[0] if j else {}; print(j.get('plan','-'), j.get('modes','-'), j.get('stop',0), j.get('driver','tie'), j.get('wantg','-'))" $K)
  set -- $JOB
  ( CT=$T PY=$PY nice -n 10 bash $T/sh/chem_c3_relrun.sh $R/final_nb.npz $K "$1" "$2" "$3" $R/r$K "$4" "$5" $c > /dev/null 2>&1 ) &
done
if [ -n "$P" ]; then
  mkdir -p vin
  $PY -I $T/r58_chem_kept_vpad.py final_nb.npz --ks $(echo $P | tr ' ' ',') --out vin/pad.npz > vpad.log 2>&1
  python3 -c "import json, sys; json.dump({'V': json.load(open(sys.argv[1]))['V']}, open(sys.argv[2], 'w'))" vin/pad.npz.json vin/proj.json \
    || log "vpad FAILED: $(tail -n 2 vpad.log)"
  for K in $P; do
    # 19:17Z: S5T with FOUR processes sharing the lucky classes (two left k11 of the rehearsal at 8+ minutes)
    c=${CO[$((i % ${#CO[@]}))]}; S4="$c ${CO[$(((i + 1) % ${#CO[@]}))]} ${CO[$(((i + 2) % ${#CO[@]}))]} ${CO[$(((i + 3) % ${#CO[@]}))]}"
    i=$((i + 4))
    ( COMBOS="0 1|1 1|2 1" CT=$T PY=$PY nice -n 10 taskset -c $c bash $T/sh/chem_c3_bsrun.sh $R/vin $K $R/b$K 3 > b$K.out 2>&1
      TG=$(grep -h 'rc=0 .*_s5' b$K/bsrun_k$K.txt 2>/dev/null | tail -n 1 | sed -n 's/.*rc=0 \(.*\)_s5[p0-9]*: EMITTED.*/\1/p')
      PA=$R/b$K/planA_$TG.json
      [ -n "$TG" ] && [ -s "$PA" ] && CT=$T PY=$PY bash $T/sh/chem_c3_s5t.sh $R/vin/pad.npz $R/vin/proj.json "$PA" $K $R/b$K k${K}t \
          "$S4" > /dev/null 2>&1 ) &
  done
fi
# S5T+ (k1, and k2 rebuilt by it on the same words): wait up to S5WAIT s while the neighbours and the refit run
G=0
if [ "$S5" = 1 ]; then
  for i in $(seq 1 $(( ${S5WAIT:-900} / 15 ))); do [ -s s5tp/DONE ] && break; sleep 15; done
  if grep -q '^ok' s5tp/DONE 2>/dev/null; then
    $PY -I $T/chem_c3_transplant.py --nooutside final.npz final_g.npz 1=s5tp/words.npz > tp1.log 2>&1 && mv final_g.npz final.npz \
      && { mkdir -p k1/use && cp s5tp/angle_k_1.kmx k1/use/ && rm -rf n2 && cp -r s5tp/n2 n2 && OWN="$OWN 1"; G=1
           log "S5T+ k1 grafted: $(cat s5tp/DONE); $(tail -n 1 tp1.log | cut -c1-200)"; } \
      || log "S5T+ graft FAILED: $(tail -n 2 tp1.log | tr '\n' ' ' | cut -c1-200)"
  else
    log "S5T+ k1 not used: $(cat s5tp/DONE 2>/dev/null || echo 'not done in time')"; pkill -f "$R/s5tp" 2>/dev/null
  fi
fi
wait
if [ "$G" = 1 ] && [ -s refit.npz ]; then        # the same dead-word graft in the refit submission (weights unchanged)
  $PY -I $T/chem_c3_transplant.py --nooutside refit.npz refit_g.npz 1=s5tp/words.npz > tp1r.log 2>&1 && mv refit_g.npz refit.npz \
    || log "S5T+ graft into refit.npz FAILED: $(tail -n 1 tp1r.log | cut -c1-150)"
fi
[ -n "$TP" ] || [ "$G" = 1 ] || { log "nothing: no R65 owner and no S5T+"; echo none > DONE; exit 0; }
for K in $W; do log "wired k$K: $(tail -n 1 n$K/RESULT.txt 2>/dev/null | cut -c1-150)"; done
for K in $L; do log "relation-bearing k$K: $(tail -n 1 r$K/relrun_k$K.txt 2>/dev/null | cut -c1-150)"; done
for K in $P; do log "30-output k$K: $(grep -h 'rc=0' b$K/bsrun_k$K.txt 2>/dev/null | tail -n 1 | cut -c1-120); S5T: $(grep -h EMITTED b$K/s5t_k$K.txt 2>/dev/null | head -n 1 | cut -c1-80)"; done
echo "ok owners$OWN neighbours wired:$W relation:$L 30-output:$P" > DONE
log done
