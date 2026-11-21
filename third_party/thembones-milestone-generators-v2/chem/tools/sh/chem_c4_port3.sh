#!/bin/bash
# chem_c4_port3.sh: builder C4 (holder 5, Nov 19 00:03Z). ON BOX 6, DETACHED. My RIGHT CHAIN k5-k14 on a NEW table with
# today's levers, in build order (batch b is read by loaders k_b and k_b+1; each loader is built on the words the
# previous step leaves): k6 by design seat R86's group-4 form (batches 5 + 6) -> k7 S5T+ (builder C3's script, own
# basis; plain S5T if it fails) -> k8 dead23 -> k9 by R86's form on its OWN batch 9 (row24b if no layout closes) ->
# k10 row24b -> k12 by R86's form (batches 11 + 12, k11's V kept clean) -> k11 plain S5T || k13 / k14 dead23 -> k5
# plain S5T on [LEFT's batches 0-4 + my batches 5-14] once OUT/LEFT.npz exists (builder C3's entry words on the table).
# Writes OUT/port.log, OUT/parts.txt (5..14), OUT/T5.npz (= the entry's words), OUT/DONE.
#   OUT=... T0=<new table> WI=2,4,6,8,9,10,12,13,14 CO="c1 ... c13" setsid nohup bash chem_c4_port3.sh > /dev/null 2>&1 < /dev/null &
set -u
T=$HOME/Quantum/tools/sh; PY=$HOME/qcc/.venv_work/bin/python
mkdir -p $OUT; rm -f $OUT/DONE
CA=($CO); NC=${#CA[@]}
log() { echo "$(date -u +%H:%M:%SZ) $*" >> $OUT/port.log; }
vk() {  # vk K TABLE TAG -> prints K's V by its padding (R58's r58_chem_kept_vpad.py, read-only)
  mkdir -p $OUT/vk_$3
  $PY -I ~/Quantum/tools/r58_chem_kept_vpad.py $2 --ks $1 --out $OUT/vk_$3/pad.npz > $OUT/vk_$3/vpad.log 2>&1
  python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['V'][sys.argv[2]])" $OUT/vk_$3/pad.npz.json $1 2>/dev/null
}
r86() {  # r86 K TABLE own|both TAG: sets RK (circuit), RW (its table); fallback = the row24b build of the pairs run
  local k=$1 tb=$2 md=$3 w=$OUT/$4 i=0 j=0 p d n best
  RK=""; RW=""
  ( cd ~/Quantum/tools && export CHEM_HB=77 PY=$PY C4_DROP=1 MAXP=${R86_MAXP:-48} SAMPLES=150 CL="${CA[*]:0:8}" \
      C4_T1=148 C4_DEAD_ROW=3:2:72 C4_DEAD_ROW_TRIALS=100 C4_DEAD_OUT=EMIT && bash sh/chem_c4_loader2.sh $k 0 0 $tb $w/lp ) > /dev/null 2>&1
  log "k$k row24b: $(grep -h -E 'SUCCESS|end without' $w/lp/L$k/run.log 2>/dev/null | tail -n 1 | cut -c1-80); pairs $(ls -d $w/lp/L$k/pairs/p* 2>/dev/null | wc -l)"
  for p in $(ls -d $w/lp/L$k/pairs/p* 2>/dev/null); do
    ( R86_BATCHES=$md CENSUS_ONLY=1 bash $T/chem_c4_r86go.sh $k $tb $p $w/s/$(basename $p) ${CA[$((i % NC))]} > /dev/null 2>&1 ) &
    i=$((i + 1)); [ $((i % NC)) = 0 ] && wait
  done
  wait
  best=$(for d in $w/s/p*; do n=$(sed -n 's/.*census: \([0-9]*\) \([0-9]*\) .*/\2/p' $d/r86go.log 2>/dev/null | head -n 1); \
    [ -n "$n" ] && [ "$n" -gt 0 ] && echo "$n $(basename $d)"; done | sort -n -r | head -n $NC | awk '{print $2}')
  log "k$k deliverable pairs ($md): $(echo $best)"
  for p in $best; do
    ( R86_BATCHES=$md bash $T/chem_c4_r86go.sh $k $tb $w/lp/L$k/pairs/$p $w/b/$p ${CA[$((j % NC))]} > /dev/null 2>&1 ) &
    j=$((j + 1))
  done
  wait
  for p in $best; do
    [ -s $w/b/$p/angle_k_$k.kmx ] && grep -q "DONE emit_c4 toffoli 221" $w/b/$p/r86go.log && { RK=$w/b/$p/angle_k_$k.kmx; RW=$w/b/$p/table_dead.npz; break; }
  done
  if [ -z "$RK" ] && [ -s $w/lp/L$k/best/angle_k_$k.kmx ]; then RK=$w/lp/L$k/best/angle_k_$k.kmx; RW=$w/lp/L$k/best/table_dead.npz; fi
  log "k$k: ${RK:-none}"
}
chain1() {  # chain1 K MODE TABLE TAG: one loader of chem_c4_chain.sh; sets RK, RW
  CL="${CA[*]:0:8}" bash $T/chem_c4_chain.sh $3 $OUT/$4 "$1:$2" > $OUT/$4.log 2>&1
  RK=$OUT/$4/L$1/best/angle_k_$1.kmx; RW=$OUT/$4/table_final.npz
  log "k$1 $2: $(tail -n 2 $OUT/$4/chain.log 2>/dev/null | head -n 1 | cut -c1-120)"
  [ -s $RK ] || RK=""
}
vroute() {  # vroute K TABLE TAG: kept loader K by S5T+ on its own batch, builder C3's V route as it is (o29 probe, then its
  # V screen: variant B with the neighbour's relation, the first V that passes phase 1 gets the full plan); sets VK, VW
  local k=$1 tb=$2 w=$OUT/$3
  VK=""; VW=""; mkdir -p $w
  ( cd ~/Quantum/tools && CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 nice -n 10 $PY -I chem_c3_o29probe.py $tb --k $k --n 16 \
      --cores $(echo ${CA[*]:0:8} | tr ' ' ',') --out $w/probe > $w/probe.log 2>&1 )
  log "k$k probe: $(tail -n 1 $w/probe.log | cut -c1-120)"
  ( cd $T && K=$k TB=$tb PROBE=$w/probe CORES="${CA[*]:0:8}" RELNBR="0:1 1:1 2:1 3:1 4:1" WI=$WI ALLV=0 \
      OUT=$w/scr bash chem_c3_vscreen.sh > /dev/null 2>&1 )
  log "k$k V route: $(cat $w/scr/DONE 2>/dev/null) $(tail -n 1 $w/scr/vscreen.log 2>/dev/null | cut -c1-120)"
  [ "$(cat $w/scr/DONE 2>/dev/null)" = ok ] && { VK=$w/scr/full/s5tp/angle_k_$k.kmx; VW=$w/scr/B_table.npz; }
}
log "start: T0 $T0, cores $CO, K5FIRST ${K5FIRST:-0}"
: > $OUT/parts.txt
# holder 5 03:17Z: START=9 = T0 already carries the final batches 0-8 (e.g. a k8 built elsewhere): k9-k14 only, then the
# run ends with T.npz = W14 and COLS 45-74 (no LEFT / k5 step)
if [ "${START:-6}" = 9 ]; then W8=$T0; else
K6MODE=both
# holder 5 00:25Z: K5FIRST=1 (T0 must carry the final batches 0-4): k5 S5T+ on its own batch 5 FIRST, then k6 by the form
# on its OWN batch 6 (k5 250 + k6 221 in place of k5 plain 251 + k6 'both' 221); no k5 at the end
if [ "${K5FIRST:-0}" = 1 ]; then
  vroute 5 $T0 k5v
  [ -n "$VK" ] || { echo "fail k5 S5T+" > $OUT/DONE; exit 1; }
  echo "5=$VK" >> $OUT/parts.txt; T0=$VW; K6MODE=own
fi
r86 6 $T0 $K6MODE k6; [ -n "$RK" ] || { echo "fail k6" > $OUT/DONE; exit 1; }; echo "6=$RK" >> $OUT/parts.txt; W6=$RW
r=$OUT/k7/c_own; mkdir -p $r
DEADBITS=2,3,4,5 NONBR=1 S5TP_SECS=2400 CORES="${CA[*]:0:8}" bash $T/chem_c3_s5tp.sh $r $W6 7 $WI > /dev/null 2>&1
st=$(cat $r/s5tp/DONE 2>/dev/null); log "k7 S5T+ own: $st"
if [ "${st%% *}" = ok ]; then echo "7=$r/s5tp/angle_k_7.kmx" >> $OUT/parts.txt; W7=$r/s5tp/words.npz
else
  vroute 7 $W6 k7v
  if [ -n "$VK" ]; then echo "7=$VK" >> $OUT/parts.txt; W7=$VW
  else
    S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $W6 7 $OUT/k7/p > /dev/null 2>&1
    log "k7 plain: $(cat $OUT/k7/p/DONE 2>/dev/null)"; echo "7=dir:$OUT/k7/p/b7" >> $OUT/parts.txt; W7=$W6
  fi
fi
# holder 5 04:07Z (R89-CHEM-CASCADE's 03:35Z request): k8 by R89's form (full second column kill on group 4, 216 on pB1q)
# through its port (tools/sh/r89_chem_cascade_port.sh, read-only), my dead23 as its FALLBACK
FALLBACK=1 bash $T/r89_chem_cascade_port.sh 8 $W7 $OUT/R8/L8 "${CA[*]:0:8}" > $OUT/R8.log 2>&1
RK=$OUT/R8/L8/best/angle_k_8.kmx; RW=$OUT/R8/L8/best/table_dead.npz; [ -s $RW ] || RW=$W7; [ -s $RK ] || RK=""
log "k8 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R8/L8/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
[ -n "$RK" ] || { echo "fail k8" > $OUT/DONE; exit 1; }; echo "8=$RK" >> $OUT/parts.txt; W8=$RW
fi
r86 9 $W8 own k9; [ -n "$RK" ] || { echo "fail k9" > $OUT/DONE; exit 1; }; echo "9=$RK" >> $OUT/parts.txt; W9=$RW
VK11=$(vk 11 $W9 k11a)   # holder 5 00:05Z: k10 by the form on its own batch 10 too (pB1q: 221 on 1 of 6 pairs), k11's V kept clean
[ -n "$VK11" ] && export C4_DEAD_XV="10:$VK11"
r86 10 $W9 own k10; unset C4_DEAD_XV
[ -n "$RK" ] || { echo "fail k10" > $OUT/DONE; exit 1; }; echo "10=$RK" >> $OUT/parts.txt; W10=$RW
# holder 5 00:14Z: TAILV=1: k11 S5T+ by builder C3's V route (o29 probe + V screen, variant B with k12's relation) on W10
# (batch 11 still raw here), then k12 by the form on its OWN batch 12 (pB1q: 221 on 2 of 3 pairs); if k11 does not close:
# k12 'both' (k11's V kept clean) and k11 plain, as before.
K11=""
if [ "${TAILV:-1}" = 1 ]; then
  vroute 11 $W10 k11v
  [ -n "$VK" ] && { K11=$VK; W11=$VW; }
fi
if [ -n "$K11" ]; then
  echo "11=$K11" >> $OUT/parts.txt
  r86 12 $W11 own k12o
  # holder 5 03:43Z: the dry run's k12 found no own-batch layout on its 48 pairs: a second, wider pair set before the fallback
  if [ -n "$RK" ] && [ "${RK#*/k12o/b/}" = "$RK" ]; then RK0=$RK; RW0=$RW; R86_MAXP=96 r86 12 $W11 own k12o2
    [ -n "$RK" ] && [ "${RK#*/k12o2/b/}" != "$RK" ] || { RK=$RK0; RW=$RW0; }; fi
  [ -n "$RK" ] || { echo "fail k12" > $OUT/DONE; exit 1; }; echo "12=$RK" >> $OUT/parts.txt; W12=$RW
  CL="${CA[*]:0:8}" bash $T/chem_c4_chain.sh $W12 $OUT/R13 "13:dead23" > $OUT/R13.log 2>&1
  FALLBACK=1 bash $T/r89_chem_cascade_port.sh 14 $OUT/R13/table_final.npz $OUT/R13/L14 "${CA[*]:0:8}" > $OUT/R14.log 2>&1
  [ -s $OUT/R13/L14/best/table_dead.npz ] && cp $OUT/R13/L14/best/table_dead.npz $OUT/R13/table_final.npz
  log "k14 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R13/L14/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
  log "k13 / k14: $(grep -h -o 'k1[34] dead23: [0-9:]* SUCCESS.*toffoli [0-9]*' $OUT/R13/chain.log | sed 's/SUCCESS.*toffoli/=/' | tr '\n' ' ')"
else
  VK11=$(vk 11 $W10 k11)
  [ -n "$VK11" ] && export C4_DEAD_XV="11:$VK11"
  log "k11's V on W10: ${VK11:-none}"
  r86 12 $W10 both k12; unset C4_DEAD_XV
  [ -n "$RK" ] || { echo "fail k12" > $OUT/DONE; exit 1; }; echo "12=$RK" >> $OUT/parts.txt; W12=$RW
  S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $W12 11 $OUT/k11p > /dev/null 2>&1 &
  P11=$!
  CL="${CA[*]:4:8}" bash $T/chem_c4_chain.sh $W12 $OUT/R13 "13:dead23" > $OUT/R13.log 2>&1
  FALLBACK=1 bash $T/r89_chem_cascade_port.sh 14 $OUT/R13/table_final.npz $OUT/R13/L14 "${CA[*]:4:8}" > $OUT/R14.log 2>&1
  [ -s $OUT/R13/L14/best/table_dead.npz ] && cp $OUT/R13/L14/best/table_dead.npz $OUT/R13/table_final.npz
  log "k14 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R13/L14/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
  wait $P11
  log "k11 plain: $(cat $OUT/k11p/DONE 2>/dev/null) $(tail -n 1 $OUT/k11p/s5one.log 2>/dev/null | cut -c1-100)"
  log "k13 / k14: $(grep -h -o 'k1[34] dead23: [0-9:]* SUCCESS.*toffoli [0-9]*' $OUT/R13/chain.log | sed 's/SUCCESS.*toffoli/=/' | tr '\n' ' ')"
  grep -q ok $OUT/k11p/DONE 2>/dev/null || { echo "fail k11" > $OUT/DONE; exit 1; }
  echo "11=dir:$OUT/k11p/b11" >> $OUT/parts.txt
fi
for k in 13 14; do [ -s $OUT/R13/L$k/best/angle_k_$k.kmx ] && echo "$k=$OUT/R13/L$k/best/angle_k_$k.kmx" >> $OUT/parts.txt; done
W14=$OUT/R13/table_final.npz
[ -s $W14 ] || { echo "fail k13-14" > $OUT/DONE; exit 1; }
if [ "${K5FIRST:-0}" = 1 ]; then cp $W14 $OUT/T5.npz; echo ok > $OUT/DONE; log "done (K5FIRST): $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-400)"; exit 0; fi
if [ "${START:-6}" = 9 ]; then cp $W14 $OUT/T.npz; echo "45-74" > $OUT/COLS; echo ok > $OUT/DONE; log "done (START 9): $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-400)"; exit 0; fi
log "right chain done; waiting for OUT/LEFT.npz (C3's entry words) for k5"
while [ ! -s $OUT/LEFT.npz ]; do sleep 30; done
$PY - $OUT/LEFT.npz $W14 $OUT/T5.npz <<'EOP' >> $OUT/port.log 2>&1
import sys, numpy as np
l = dict(np.load(sys.argv[1])); r = np.load(sys.argv[2])["angle_words"]; wl = np.asarray(l["angle_words"])
w = wl.copy(); w[:, 25:75] = r[:, 25:75]
d = np.argwhere(w[:, :25] != r[:, :25]); print("LEFT vs my chain's raw batches 0-4: columns", sorted(set(d[:, 1].tolist())))
l["angle_words"] = w.astype(wl.dtype); np.savez(sys.argv[3], **l)
EOP
# holder 5 (Nov 19 02:44Z): CASCADE=1 (default): the LEFT CASCADE on LEFT's words (C3's k0 and weights kept): k2 by the
# form 'both' (batches 1 + 2, toggles off k1's and k3's V blocks at points 0-2) -> k1 plain S5T || k3 S5T+ by C3's V
# route on k2's batch 2 -> k4 row24b -> k5 plain on [k4's batch 4 + my batch 5]: -1 a pass against C3's k1 S5T+ / k2
# plain (on pB1q: c4_pB1q7). Any step missing -> the plain ending below (k5 plain on LEFT's batch 4), parts 5-14 only.
# holder 5 03:10Z: WIRED10=1 = the TEN-wired + P57x layout (wired 2,3,4,6,8,9,10,12,13,14: k3 wired): k3 is then built
# by the form on its OWN batch 3 (row24b if no layout closes) after k2, and only k1's V is kept clean (UNTESTED: no
# ten-wired P57x table existed at 03:10Z).
if [ "${CASCADE:-1}" = 1 ]; then
  V1=$(vk 1 $OUT/T5.npz c1)
  if [ "${WIRED10:-0}" = 1 ]; then V3=""; [ -n "$V1" ] && export C4_DEAD_XV="1:$V1"
  else V3=$(vk 3 $OUT/T5.npz c3); [ -n "$V1" ] && [ -n "$V3" ] && export C4_DEAD_XV="1:$V1,2:$V3"; fi
  log "cascade: k1's V ${V1:-none}, k3's V ${V3:-none} (WIRED10 ${WIRED10:-0})"
  r86 2 $OUT/T5.npz both k2c; unset C4_DEAD_XV
  if [ -n "$RK" ] && [ "${RK#*/k2c/b/}" != "$RK" ]; then
    K2=$RK; W2=$RW
    S5T_SECS=2400 CORES="${CA[*]:8:4}" bash $T/chem_c3_s5one.sh $W2 1 $OUT/k1c > /dev/null 2>&1 &
    P1=$!
    if [ "${WIRED10:-0}" = 1 ]; then r86 3 $W2 own k3c; K3=$RK; W3=$RW
    else vroute 3 $W2 k3c; K3=$VK; W3=$VW; fi
    wait $P1
    log "cascade: k2 $K2; k1 plain $(cat $OUT/k1c/DONE 2>/dev/null); k3 ${K3:-none}"
    if grep -q ok $OUT/k1c/DONE 2>/dev/null && [ -n "$K3" ]; then
      chain1 4 row24b $W3 R4c; K4=$RK; W4=$RW
      if [ -n "$K4" ]; then
        S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $W4 5 $OUT/k5c > /dev/null 2>&1
        log "cascade: k4 $K4; k5 plain $(cat $OUT/k5c/DONE 2>/dev/null)"
        if grep -q ok $OUT/k5c/DONE 2>/dev/null; then
          { echo "1=dir:$OUT/k1c/b1"; echo "2=$K2"; echo "3=$K3"; echo "4=$K4"; echo "5=dir:$OUT/k5c/b5"; } >> $OUT/parts.txt
          cp $W4 $OUT/T.npz; echo "5-74" > $OUT/COLS
          echo ok > $OUT/DONE; log "done (CASCADE): $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-500)"; exit 0
        fi
      fi
    fi
  fi
  log "cascade incomplete: the plain ending"
fi
S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $OUT/T5.npz 5 $OUT/k5p > /dev/null 2>&1
log "k5 plain: $(cat $OUT/k5p/DONE 2>/dev/null) $(tail -n 1 $OUT/k5p/s5one.log 2>/dev/null | cut -c1-100)"
grep -q ok $OUT/k5p/DONE 2>/dev/null || { echo "fail k5" > $OUT/DONE; exit 1; }
echo "5=dir:$OUT/k5p/b5" >> $OUT/parts.txt
cp $OUT/T5.npz $OUT/T.npz; echo "25-74" > $OUT/COLS
echo ok > $OUT/DONE; log "done: $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-400)"
