#!/bin/bash
# chem_c4_port5.sh: builder C4 (holder 7, Nov 19 08:00Z). Copy of chem_c4_port4.sh (holder 6) with ONE NEW MODE, FAST=1:
#   THE RIGHT HALF IN THREE SEGMENTS AT ONCE (right_split): loader k builds delta = batch k-1 XOR batch k, so a dead cell
#   toggled in batch k moves loader k+1 too and port4's k6..k14 chain is sequential (73.5-85 min). A loader that toggles
#   nothing breaks the chain: k7 and k11 (kept, 30 outputs) are built PLAIN (251 against S5T+ 250: +1 a pass each [M]),
#   so S1 = k6 'both' -> k7 plain, S2 = k8 -> k9 -> k10 (all on the raw words of batches 7-10) and S3 = k12 'both'
#   (k11's V kept clean) -> k13 -> k14 run at once, then k11 plain on [S2's batch 10 + S3's batch 11]. The LEFT half is
#   port4's left10 (k2 'both' || k4 'both', then k1 plain || k3 plain25: K6-1's measured path; +1 against left10f's k3
#   form) for LAYOUT 10, left9 for LAYOUT 9. Everything else (k5 early, merges, parts, DONE) is port4's.
#   Price by plan [P]: ten-wired about 3,470-3,471 a pass (port4's exact path 3,467); time [E] about 35-45 min on 24+ cores.
#   FAST=0 (default) = port4 exactly.
# --- port4's header follows ---
# chem_c4_port4.sh: builder C4 (holder 6, Nov 19 05:25Z). ON BOX 5 OR 6, DETACHED. THE UNATTENDED PORT OF A NEW TABLE for the
# NINE-wired (pB1q's P57x layout), TEN-wired (G10P57x: k3 wired with NO bit-1 relation = 25 generic outputs) and
# ELEVEN-wired (G11w11: wires only, every wired loader 25 outputs) layouts. Copy of chem_c4_port3.sh (holder 5) with:
#  - LEFT (k1-k4) and RIGHT (k6-k14) START TOGETHER on the raw table T0 (C3's k0 never changes batch 0 [M dry run:
#    LEFT vs raw = columns 5-9, 15-20 only]); k5 plain LAST on [LEFT's batches 0-4 + RIGHT's batches 5-14];
#  - LAYOUT=10: k2 by R86's form 'both' (k1's V clean) AND k4 by the form 'both' AT ONCE (disjoint batches 1-2 / 3-4:
#    k3 is wired and plain, it toggles nothing), then k1 plain S5T and k3 plain 25-output (loader2: 228) AT ONCE;
#  - LAYOUT=11: k1 / k5 / k7 S5T+ by C3's V route first (disjoint batches), then every wired loader by loader2 at once;
#  - R86 builds through R86GO (default design seat R89's strict r86go: a missed drop retries the group).
# Writes OUT/port.log, OUT/parts.txt (1..14), OUT/T.npz (the entry's words), OUT/COLS (5-74), OUT/DONE ('ok' or 'fail ...').
#   OUT=... T0=<table> REL=<its rel_used_g.txt> WI=<wired list> LAYOUT=9|10|11 RC="<right cores, 8+>" LC="<left cores, 8+>" \
#     setsid nohup bash <a copy of chem_c4_port4.sh> > /dev/null 2>&1 < /dev/null &
# REL gives each kept loader's V route its neighbour's real bit-1 relations (port3 passed a 24-output one to all).
set -u
T=$HOME/Quantum/tools/sh; PY=$HOME/qcc/.venv_work/bin/python
R86GO=${R86GO:-$T/r89_chem_cascade_r86go.sh}; [ -s "$R86GO" ] || R86GO=$T/chem_c4_r86go.sh
mkdir -p $OUT; rm -f $OUT/DONE $OUT/RIGHT.res $OUT/LEFT.res $OUT/k6.res $OUT/ng.txt $OUT/k4.res $OUT/k0w.part $OUT/k5c.part
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
  # holder 7: R86_CMIN=n (FAST mode default 3) ends the census after the first wave of NC pairs that leaves n deliverable
  # pairs (K7-1: R86's census is 1.4-1.8 min a pair on these words, so 48 pairs on 6 cores = 8 waves = 14.6 min; the
  # builds closed 221 on the first deliverable pair of k4, k6, k12). 0 = every pair (port4's rule).
  local nd
  for p in $(ls -d $w/lp/L$k/pairs/p* 2>/dev/null); do
    ( R86_BATCHES=$md CENSUS_ONLY=1 bash $R86GO $k $tb $p $w/s/$(basename $p) ${CA[$((i % NC))]} > /dev/null 2>&1 ) &
    i=$((i + 1))
    if [ $((i % NC)) = 0 ]; then
      wait
      if [ "${R86_CMIN:-0}" -gt 0 ]; then
        nd=$(for d in $w/s/p*; do sed -n 's/.*census: \([0-9]*\) \([0-9]*\) .*/\2/p' $d/r86go.log 2>/dev/null | head -n 1; done | awk '$1 > 0' | wc -l)
        [ "$nd" -ge "${R86_CMIN}" ] && { log "k$k census stopped after $i pairs: $nd deliverable (R86_CMIN ${R86_CMIN})"; break; }
      fi
    fi
  done
  wait
  best=$(for d in $w/s/p*; do n=$(sed -n 's/.*census: \([0-9]*\) \([0-9]*\) .*/\2/p' $d/r86go.log 2>/dev/null | head -n 1); \
    [ -n "$n" ] && [ "$n" -gt 0 ] && echo "$n $(basename $d)"; done | sort -n -r | head -n $NC | awk '{print $2}')
  log "k$k deliverable pairs ($md): $(echo $best)"
  for p in $best; do
    ( R86_BATCHES=$md bash $R86GO $k $tb $w/lp/L$k/pairs/$p $w/b/$p ${CA[$((j % NC))]} > /dev/null 2>&1 ) &
    j=$((j + 1))
  done
  wait
  for p in $best; do
    [ -s $w/b/$p/angle_k_$k.kmx ] && grep -q "DONE emit_c4 toffoli 221" $w/b/$p/r86go.log && { RK=$w/b/$p/angle_k_$k.kmx; RW=$w/b/$p/table_dead.npz; break; }
  done
  if [ -z "$RK" ] && [ -s $w/lp/L$k/best/angle_k_$k.kmx ]; then RK=$w/lp/L$k/best/angle_k_$k.kmx; RW=$w/lp/L$k/best/table_dead.npz; fi
  if [ -z "$RK" ]; then   # holder 7 (13:45Z; g2c's ten-wired tables: row24b closed only at 149 > its 148 target and the
    # form missed on every deliverable pair, so k2 had NO circuit and the left failed): the best CLOSED b1 of the pairs
    local bp; bp=$(for f in $w/lp/L$k/pairs/p*/b1.log; do n=$(grep -h -o -E 'CLOSED [0-9]+ products, emit_c4 toffoli [0-9]+ .checked.' $f 2>/dev/null | head -n 1 | awk '{print $2}'); \
      [ -n "$n" ] && [ -s $(dirname $f)/b1/angle_k_$k.kmx ] && echo "$n $(dirname $f)"; done | sort -n | head -n 1 | awk '{print $2}')
    [ -n "$bp" ] && { RK=$bp/b1/angle_k_$k.kmx; RW=$bp/b1/table_dead.npz; log "k$k: the best closed b1 ($bp)"; }
  fi
  [ -n "$RK" ] && [ ! -s "$RW" ] && RW=$tb   # a build without toggles leaves the words as they were
  log "k$k: ${RK:-none}"
}
plain25() {  # plain25 K TABLE TAG: a 25-output wired loader by loader2 (plan 6|6|7|6, group 3 annealed: 154 = 228); no toggles
  local k=$1 tb=$2 w=$OUT/$3
  RK=""
  ( cd ~/Quantum/tools && export CHEM_HB=77 PY=$PY MAXP=${P25_MAXP:-48} SAMPLES=150 CL="${CA[*]:0:8}" && \
      bash sh/chem_c4_loader2.sh $k 0 0 $tb $w ) > /dev/null 2>&1
  if [ -s $w/L$k/best/angle_k_$k.kmx ]; then RK=$w/L$k/best/angle_k_$k.kmx
  elif [ -s $w/L$k/best/angle_k_${k}_t1.kmx ]; then RK=$w/L$k/best/angle_k_${k}_t1.kmx; fi
  log "k$k plain25: $(grep -h -E 'SUCCESS|end without|kept|no plan' $w/L$k/run.log 2>/dev/null | tail -n 1 | cut -c1-110) -> ${RK:-none}"
}
chain1() {  # chain1 K MODE TABLE TAG: one loader of chem_c4_chain.sh; sets RK, RW
  CL="${CA[*]:0:8}" bash $T/chem_c4_chain.sh $3 $OUT/$4 "$1:$2" > $OUT/$4.log 2>&1
  RK=$OUT/$4/L$1/best/angle_k_$1.kmx; RW=$OUT/$4/table_final.npz
  log "k$1 $2: $(tail -n 2 $OUT/$4/chain.log 2>/dev/null | head -n 1 | cut -c1-120)"
  [ -s $RK ] || RK=""
}
relnbr() {  # relnbr K: the bit-1 relations of loader K+1 in the table's relation file REL ('j:b j:b ...;j:b ...'; empty if
  # none). holder 6: port3 passed '0:1 1:1 2:1 3:1 4:1' for every kept loader (right only for a 24-output neighbour)
  [ -s "${REL:-}" ] || { echo "0:1 1:1 2:1 3:1 4:1"; return; }
  awk -F'|' -v k=$(( $1 + 1 )) '($1 + 0) == k && $2 ~ /[0-4]:[1-9]/ {gsub(/^ +| +$/, "", $2); s = s (s == "" ? "" : ";") $2} END {print s}' $REL
}
vroute() {  # vroute K TABLE TAG: kept loader K by S5T+ on its own batch (C3's V route: o29 probe, then its V screen); sets VK, VW
  local k=$1 tb=$2 w=$OUT/$3 rn
  VK=""; VW=""; mkdir -p $w; rn=$(relnbr $k)
  ( cd ~/Quantum/tools && CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 nice -n 10 $PY -I chem_c3_o29probe.py $tb --k $k --n 16 \
      --cores $(echo ${CA[*]:0:8} | tr ' ' ',') --out $w/probe > $w/probe.log 2>&1 )
  log "k$k probe: $(tail -n 1 $w/probe.log | cut -c1-120); neighbour relations '${rn}'"
  ( cd $T && K=$k TB=$tb PROBE=$w/probe CORES="${CA[*]:0:8}" RELNBR="$rn" WI=$WI ALLV=0 \
      OUT=$w/scr bash chem_c3_vscreen.sh > /dev/null 2>&1 )
  log "k$k V route: $(cat $w/scr/DONE 2>/dev/null) $(tail -n 1 $w/scr/vscreen.log 2>/dev/null | cut -c1-120)"
  [ "$(cat $w/scr/DONE 2>/dev/null)" = ok ] && { VK=$w/scr/full/s5tp/angle_k_$k.kmx; VW=$w/scr/B_table.npz; }
}
cmerge() {  # cmerge OUTNPZ BASE "lo-hi:FILE ..." : BASE's words with the column ranges taken from the files; asserts that
  # each file differs from BASE only inside its own range (a lever that toggled a word outside it stops the port)
  $PY - "$@" <<'EOP' >> $OUT/port.log 2>&1
import sys, numpy as np
out, base = sys.argv[1], sys.argv[2]
z = dict(np.load(base)); w0 = np.asarray(z["angle_words"]); w = w0.copy()
for spec in sys.argv[3].split():
    rg, f = spec.split(":", 1)
    lo, hi = (int(x) for x in rg.split("-"))
    t = np.asarray(np.load(f)["angle_words"])
    d = sorted(set(np.argwhere(t != w0)[:, 1].tolist()))
    assert all(lo <= c <= hi for c in d), f"{f} changes columns {d} outside {rg}"
    w[:, lo:hi + 1] = t[:, lo:hi + 1]
    print(f"merge: {f} columns {d}")
z["angle_words"] = w.astype(w0.dtype); np.savez(out, **z)
EOP
}

right_chain() {  # k6-k14 as in port3 (holder 5's levers + R89's k8 / k14 port); writes OUT/parts_R.txt, OUT/RIGHT.res
  : > $OUT/parts_R.txt
  local W6 W7 W8 W9 W10 W11 W12 W14 r st K11
  r86 6 $T0 both k6; [ -n "$RK" ] || { echo "fail k6" > $OUT/RIGHT.res; return; }; echo "6=$RK" >> $OUT/parts_R.txt; W6=$RW
  echo "ok $W6" > $OUT/k6.res   # holder 6: k5 (batches 4 + 5) can be built as soon as LEFT and k6 are done
  if [ "$(ngof 7)" != 30 ]; then   # holder 6 (06:10Z): k7 WIRED (the twelve-wired rung and up): its form on its own batch
    wiredown 7 $W6 k7w; [ -n "$RK" ] || { plainw 7 $W6 k7p; RW=$W6; }
    [ -n "$RK" ] || { echo "fail k7" > $OUT/RIGHT.res; return; }; echo "7=$RK" >> $OUT/parts_R.txt; W7=$RW
  else
  r=$OUT/k7/c_own; mkdir -p $r
  DEADBITS=2,3,4,5 NONBR=1 S5TP_SECS=2400 CORES="${CA[*]:0:8}" bash $T/chem_c3_s5tp.sh $r $W6 7 $WI > /dev/null 2>&1
  st=$(cat $r/s5tp/DONE 2>/dev/null); log "k7 S5T+ own: $st"
  if [ "${st%% *}" = ok ]; then echo "7=$r/s5tp/angle_k_7.kmx" >> $OUT/parts_R.txt; W7=$r/s5tp/words.npz
  else
    vroute 7 $W6 k7v
    if [ -n "$VK" ]; then echo "7=$VK" >> $OUT/parts_R.txt; W7=$VW
    else
      S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $W6 7 $OUT/k7/p > /dev/null 2>&1
      log "k7 plain: $(cat $OUT/k7/p/DONE 2>/dev/null)"; echo "7=dir:$OUT/k7/p/b7" >> $OUT/parts_R.txt; W7=$W6
    fi
  fi
  fi
  FALLBACK=1 bash $T/r89_chem_cascade_port.sh 8 $W7 $OUT/R8/L8 "${CA[*]:0:8}" > $OUT/R8.log 2>&1
  RK=$OUT/R8/L8/best/angle_k_8.kmx; RW=$OUT/R8/L8/best/table_dead.npz; [ -s $RW ] || RW=$W7; [ -s $RK ] || RK=""
  log "k8 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R8/L8/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
  [ -n "$RK" ] || { echo "fail k8" > $OUT/RIGHT.res; return; }; echo "8=$RK" >> $OUT/parts_R.txt; W8=$RW
  r86 9 $W8 own k9; [ -n "$RK" ] || { echo "fail k9" > $OUT/RIGHT.res; return; }; echo "9=$RK" >> $OUT/parts_R.txt; W9=$RW
  local VK11 K10="" b; VK11=$(vk 11 $W9 k11a)
  [ -n "$VK11" ] && export C4_DEAD_XV="10:$VK11"
  # holder 6 (05:30Z): k10 FIRST by design seat R89's supply-aware strict pairs (its right chain's step, copied: pB1q 221
  # where my pairs ended row24b 222), kept only if C3's o29 probe still finds a stage-A-149 V for k11 on its words
  bash $T/r89_chem_cascade_r86s.sh 10 $W9 $OUT/k10s "${CA[*]:0:8}" 8 > /dev/null 2>&1
  for b in $(grep -l "DONE emit_c4 toffoli 221" $OUT/k10s/b/p*/r86go.log 2>/dev/null | xargs -r -n1 dirname); do
    ( cd ~/Quantum/tools && CHEM_BROT=7 CHEM_HB=77 OMP_NUM_THREADS=1 nice -n 10 $PY -I chem_c3_o29probe.py $b/table_dead.npz \
        --k 11 --n 16 --cores $(echo ${CA[*]:0:8} | tr ' ' ',') --out $b/probe11 > $b/probe11.log 2>&1 )
    log "k10 (R89 supply-aware) $b: k11 probe $(tail -n 1 $b/probe11.log | cut -c1-80)"
    grep -q "stage A min 149" $b/probe11.log && { K10=$b; break; }
  done
  if [ -n "$K10" ]; then RK=$K10/angle_k_10.kmx; RW=$K10/table_dead.npz; log "k10: $RK (R89 supply-aware)"
  else r86 10 $W9 own k10; fi
  unset C4_DEAD_XV
  [ -n "$RK" ] || { echo "fail k10" > $OUT/RIGHT.res; return; }; echo "10=$RK" >> $OUT/parts_R.txt; W10=$RW
  K11=""
  if [ "$(ngof 11)" != 30 ]; then   # holder 6 (06:10Z): k11 WIRED (the eleven-wired rung and up): its form, own batch
    wiredown 11 $W10 k11w
    if [ -n "$RK" ]; then K11=$RK; W11=$RW; else plainw 11 $W10 k11p; [ -n "$RK" ] && { K11=$RK; W11=$W10; }; fi
  else
    vroute 11 $W10 k11v
    [ -n "$VK" ] && { K11=$VK; W11=$VW; }
  fi
  if [ -n "$K11" ]; then
    echo "11=$K11" >> $OUT/parts_R.txt
    r86 12 $W11 own k12o
    if [ -n "$RK" ] && [ "${RK#*/k12o/b/}" = "$RK" ]; then local RK0=$RK RW0=$RW; R86_MAXP=96 r86 12 $W11 own k12o2
      [ -n "$RK" ] && [ "${RK#*/k12o2/b/}" != "$RK" ] || { RK=$RK0; RW=$RW0; }; fi
    [ -n "$RK" ] || { echo "fail k12" > $OUT/RIGHT.res; return; }; echo "12=$RK" >> $OUT/parts_R.txt; W12=$RW
    CL="${CA[*]:0:8}" bash $T/chem_c4_chain.sh $W12 $OUT/R13 "13:dead23" > $OUT/R13.log 2>&1
    FALLBACK=1 bash $T/r89_chem_cascade_port.sh 14 $OUT/R13/table_final.npz $OUT/R13/L14 "${CA[*]:0:8}" > $OUT/R14.log 2>&1
    [ -s $OUT/R13/L14/best/table_dead.npz ] && cp $OUT/R13/L14/best/table_dead.npz $OUT/R13/table_final.npz
  else
    VK11=$(vk 11 $W10 k11)
    [ -n "$VK11" ] && export C4_DEAD_XV="11:$VK11"
    r86 12 $W10 both k12; unset C4_DEAD_XV
    [ -n "$RK" ] || { echo "fail k12" > $OUT/RIGHT.res; return; }; echo "12=$RK" >> $OUT/parts_R.txt; W12=$RW
    S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $W12 11 $OUT/k11p > /dev/null 2>&1 &
    local P11=$!
    CL="${CA[*]:4:8}" bash $T/chem_c4_chain.sh $W12 $OUT/R13 "13:dead23" > $OUT/R13.log 2>&1
    FALLBACK=1 bash $T/r89_chem_cascade_port.sh 14 $OUT/R13/table_final.npz $OUT/R13/L14 "${CA[*]:4:8}" > $OUT/R14.log 2>&1
    [ -s $OUT/R13/L14/best/table_dead.npz ] && cp $OUT/R13/L14/best/table_dead.npz $OUT/R13/table_final.npz
    wait $P11
    log "k11 plain: $(cat $OUT/k11p/DONE 2>/dev/null)"
    grep -q ok $OUT/k11p/DONE 2>/dev/null || { echo "fail k11" > $OUT/RIGHT.res; return; }
    echo "11=dir:$OUT/k11p/b11" >> $OUT/parts_R.txt
  fi
  log "k14 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R13/L14/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
  log "k13 / k14: $(grep -h -o 'k1[34] dead23: [0-9:]* SUCCESS.*toffoli [0-9]*' $OUT/R13/chain.log | sed 's/SUCCESS.*toffoli/=/' | tr '\n' ' ')"
  for k in 13 14; do [ -s $OUT/R13/L$k/best/angle_k_$k.kmx ] && echo "$k=$OUT/R13/L$k/best/angle_k_$k.kmx" >> $OUT/parts_R.txt; done
  W14=$OUT/R13/table_final.npz
  [ -s $W14 ] || { echo "fail k13-14" > $OUT/RIGHT.res; return; }
  echo "ok $W14" > $OUT/RIGHT.res; log "RIGHT done: $(tr '\n' ' ' < $OUT/parts_R.txt | cut -c1-400)"
}

right_split() {  # holder 7 (FAST=1): k6-k14 in three segments at once on thirds of the right cores (see the header);
  # writes OUT/parts_R.txt, OUT/k6.res (as right_chain does: k5 starts early from it), OUT/RIGHT.res
  : > $OUT/parts_R.txt; rm -f $OUT/s1.res $OUT/s2.res $OUT/s3.res $OUT/s1.parts $OUT/s2.parts $OUT/s3.parts $OUT/s3k12.res
  local n=${#RA[@]} a b VK11 r1 r2 r3 W6 W10 W14
  a=$(( n / 3 )); b=$(( 2 * n / 3 ))
  local S1=(${RA[@]:0:$a}) S2=(${RA[@]:$a:$(( b - a ))}) S3=(${RA[@]:$b})
  VK11=$(vk 11 $T0 k11t0)
  log "right_split: S1 ${S1[*]} | S2 ${S2[*]} | S3 ${S3[*]}; k11's V on T0 ${VK11:-none} (kept clean by k10 and k12)"
  # S1: k6 by R86's form on batches 5-6, then k7 plain (no toggles: batch 7 stays raw for S2's k8)
  ( CA=(${S1[@]}); NC=${#CA[@]}
    r86 6 $T0 both k6; [ -n "$RK" ] || { echo "fail k6" > $OUT/s1.res; exit 0; }
    echo "6=$RK" > $OUT/s1.parts; echo "ok $RW" > $OUT/k6.res; local w6=$RW
    plaink 7 $w6 k7p; [ -s $OUT/k7p.part ] || { echo "fail k7 plain" > $OUT/s1.res; exit 0; }
    cat $OUT/k7p.part >> $OUT/s1.parts; echo "ok $w6" > $OUT/s1.res; log "S1 done: k6, k7 plain" ) &
  # S2: k8 (R89's port, own batch 8) -> k9 own -> k10 own (R89's supply-aware pairs first; k11 is plain, so no probe)
  ( CA=(${S2[@]}); NC=${#CA[@]}; local w8 w9 bb
    FALLBACK=1 bash $T/r89_chem_cascade_port.sh 8 $T0 $OUT/R8/L8 "${CA[*]:0:8}" > $OUT/R8.log 2>&1
    RK=$OUT/R8/L8/best/angle_k_8.kmx; RW=$OUT/R8/L8/best/table_dead.npz; [ -s $RW ] || RW=$T0; [ -s $RK ] || RK=""
    log "S2 k8 (R89 form, dead23 fallback): $(grep -h -E 'SUCCESS|FALLBACK' $OUT/R8/L8/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
    [ -n "$RK" ] || { echo "fail k8" > $OUT/s2.res; exit 0; }; echo "8=$RK" > $OUT/s2.parts; w8=$RW
    r86 9 $w8 own k9; [ -n "$RK" ] || { echo "fail k9" > $OUT/s2.res; exit 0; }; echo "9=$RK" >> $OUT/s2.parts; w9=$RW
    [ -n "$VK11" ] && export C4_DEAD_XV="10:$VK11"
    bash $T/r89_chem_cascade_r86s.sh 10 $w9 $OUT/k10s "${CA[*]:0:8}" 8 > /dev/null 2>&1
    bb=$(grep -o "SUCCESS [^ ]*" $OUT/k10s/r86s.log 2>/dev/null | head -n 1 | cut -d" " -f2)
    if [ -n "$bb" ] && [ -s $bb/angle_k_10.kmx ]; then RK=$bb/angle_k_10.kmx; RW=$bb/table_dead.npz; [ -s $RW ] || RW=$w9
      log "S2 k10: $RK (R89 supply-aware)"
    else r86 10 $w9 own k10; fi
    unset C4_DEAD_XV
    [ -n "$RK" ] || { echo "fail k10" > $OUT/s2.res; exit 0; }; echo "10=$RK" >> $OUT/s2.parts
    echo "ok $RW" > $OUT/s2.res; log "S2 done: k8, k9, k10" ) &
  local P2=$!
  # S3: k12 by R86's form on batches 11-12 (k11's V kept clean in batch 11) -> k13 dead23 -> k14 (R89's port)
  ( CA=(${S3[@]}); NC=${#CA[@]}; local w12 k
    [ -n "$VK11" ] && export C4_DEAD_XV="11:$VK11"
    r86 12 $T0 both k12; unset C4_DEAD_XV
    [ -n "$RK" ] || { echo "fail k12" > $OUT/s3k12.res; echo "fail k12" > $OUT/s3.res; exit 0; }; echo "12=$RK" > $OUT/s3.parts; w12=$RW
    echo "ok $w12" > $OUT/s3k12.res   # holder 7 (v3): batch 11 is final here (k13 / k14 toggle 13 / 14 only): k11 may start
    CL="${CA[*]:0:8}" bash $T/chem_c4_chain.sh $w12 $OUT/R13 "13:dead23" > $OUT/R13.log 2>&1
    # holder 7 (v4, K7-1: R89's k14 port ran 35+ min when its first pairs' forms missed, then its dead23 fallback):
    # R89's form (216) on 8 cores RACES dead23 (217, R89's own fallback recipe) on the others; R89's result is taken
    # if it lands before k11 plain and k5 plain are written (or K14CAP s), else its jobs are stopped and dead23 is used
    local T13=$OUT/R13/table_k13.npz P89 t89
    cp $OUT/R13/table_final.npz $T13
    ( FALLBACK=0 bash $T/r89_chem_cascade_port.sh 14 $T13 $OUT/R13/L14 "${CA[*]:0:8}" > $OUT/R14.log 2>&1 ) &
    P89=$!; t89=$(date +%s)
    local C14=(${CA[@]:8}); [ ${#C14[@]} -ge 2 ] || C14=(${CA[@]:0:4})
    CL="${C14[*]}" bash $T/chem_c4_chain.sh $T13 $OUT/R14d "14:dead23" > $OUT/R14d.log 2>&1   # with its own plain fallback
    log "S3 k14 dead23: $(tail -n 2 $OUT/R14d/chain.log 2>/dev/null | head -n 1 | cut -c1-110)"
    while kill -0 $P89 2>/dev/null; do
      [ -s $OUT/k11p.part ] && [ -s $OUT/k5c.part ] && break
      [ $(( $(date +%s) - t89 )) -ge ${K14CAP:-2700} ] && break
      sleep 20
    done
    if grep -q " SUCCESS " $OUT/R13/L14/port.log 2>/dev/null && [ -s $OUT/R13/L14/best/angle_k_14.kmx ]; then
      log "S3 k14: R89's form in time ($(( $(date +%s) - t89 )) s)"
    else
      pkill -f "$OUT/R13/L14" 2>/dev/null; sleep 2; pkill -9 -f "$OUT/R13/L14" 2>/dev/null
      rm -rf $OUT/R13/L14/best; mkdir -p $OUT/R13/L14/best
      if [ -s $OUT/R14d/L14/best/angle_k_14.kmx ]; then
        cp $OUT/R14d/L14/best/angle_k_14.kmx $OUT/R13/L14/best/
        [ -f $OUT/R14d/L14/best/table_dead.npz ] && cp $OUT/R14d/L14/best/table_dead.npz $OUT/R13/L14/best/
      fi
      log "S3 k14: dead23 taken (R89's form not in time after $(( $(date +%s) - t89 )) s)"
    fi
    [ -s $OUT/R13/L14/best/table_dead.npz ] && cp $OUT/R13/L14/best/table_dead.npz $OUT/R13/table_final.npz
    log "S3 k14 (R89 form | dead23 race): $(grep -h -E 'SUCCESS|end without' $OUT/R13/L14/port.log 2>/dev/null | tail -n 1 | cut -c1-120)"
    for k in 13 14; do [ -s $OUT/R13/L$k/best/angle_k_$k.kmx ] && echo "$k=$OUT/R13/L$k/best/angle_k_$k.kmx" >> $OUT/s3.parts; done
    [ -s $OUT/R13/table_final.npz ] && [ "$(grep -c '^1[34]=' $OUT/s3.parts)" = 2 ] || { echo "fail k13-14" > $OUT/s3.res; exit 0; }
    echo "ok $OUT/R13/table_final.npz" > $OUT/s3.res; log "S3 done: k12, k13, k14" ) &
  local P3=$! P11="" r12 W12
  # holder 7 (v3, K7-1: S2 ended 08:20Z while S3's k14 ran on): k11 reads batches 10 and 11 only, final once S2's k10
  # and S3's k12 are done, so k11 plain starts then on S2's cores while k13 / k14 go on (checked against the final words)
  while [ ! -s $OUT/s2.res ] && kill -0 $P2 2>/dev/null; do sleep 20; done
  r2=""; [ -s $OUT/s2.res ] && read -r r2 W10 < $OUT/s2.res
  if [ "${r2:-}" = ok ]; then
    while [ ! -s $OUT/s3k12.res ] && kill -0 $P3 2>/dev/null; do sleep 20; done
    r12=""; [ -s $OUT/s3k12.res ] && read -r r12 W12 < $OUT/s3k12.res
    if [ "${r12:-}" = ok ] && cmerge $OUT/K11IN.npz $T0 "40-54:$W10 55-64:$W12"; then
      ( CA=(${S2[@]}); NC=${#CA[@]}; plaink 11 $OUT/K11IN.npz k11p ) &
      P11=$!; log "k11 plain started early on [S2's batch 10 + S3's k12 batch 11]"
    fi
  fi
  wait
  read -r r1 W6 < $OUT/s1.res 2>/dev/null; read -r r2 W10 < $OUT/s2.res 2>/dev/null; read -r r3 W14 < $OUT/s3.res 2>/dev/null
  log "right_split segments: S1 $(cat $OUT/s1.res 2>/dev/null) | S2 $(cat $OUT/s2.res 2>/dev/null) | S3 $(cat $OUT/s3.res 2>/dev/null)"
  [ "${r1:-}" = ok ] && [ "${r2:-}" = ok ] && [ "${r3:-}" = ok ] || { echo "fail segments ${r1:-} ${r2:-} ${r3:-}" > $OUT/RIGHT.res; return; }
  # each segment's words differ from T0 only inside its own batches (cmerge asserts it): 5-6 | 8-10 | 11-14; batch 7 raw
  cmerge $OUT/WR.npz $T0 "25-34:$W6 40-54:$W10 55-74:$W14" || { echo "fail merge S1 / S2 / S3" > $OUT/RIGHT.res; return; }
  if [ -n "$P11" ] && [ -s $OUT/k11p.part ]; then   # the early k11 is valid only if batches 10-11 are the final ones
    $PY - $OUT/K11IN.npz $OUT/WR.npz <<'EOP' >> $OUT/port.log 2>&1 || { log "k11 early: batches 10-11 moved: rebuilt"; rm -rf $OUT/k11p $OUT/k11p.part; }
import sys, numpy as np
a = np.asarray(np.load(sys.argv[1])["angle_words"]); b = np.asarray(np.load(sys.argv[2])["angle_words"])
assert (a[:, 50:60] == b[:, 50:60]).all(), "batches 10-11 differ"
print("k11 early: batches 10-11 equal to the final words")
EOP
  fi
  if [ ! -s $OUT/k11p.part ]; then CA=(${RA[@]}); NC=${#CA[@]}; plaink 11 $OUT/WR.npz k11p; fi
  [ -s $OUT/k11p.part ] || { echo "fail k11 plain" > $OUT/RIGHT.res; return; }
  { cat $OUT/s1.parts; cat $OUT/s2.parts; cat $OUT/k11p.part; cat $OUT/s3.parts; } > $OUT/parts_R.txt
  echo "ok $OUT/WR.npz" > $OUT/RIGHT.res; log "RIGHT (split) done: $(tr '\n' ' ' < $OUT/parts_R.txt | cut -c1-400)"
}

left10() {  # the ten-wired G10P57x left: k2 'both' and k4 'both' at once, then k1 plain and k3 plain25 at once
  : > $OUT/parts_L.txt
  local V1 V5 K2 W2 K4 W4 n2 n4
  V1=$(vk 1 $T0 c1); V5=$(vk 5 $T0 c5); log "left10: k1's V ${V1:-none}, k5's V ${V5:-none}"
  n2=$(( ${#LA[@]} / 2 ))
  [ "${L10K2:-}" = "2/3" ] && n2=$(( 2 * ${#LA[@]} / 3 ))   # holder 7 (FAST): k2 is the slow one
  # k2 keeps k1's V blocks clean in batch 1; k4 keeps k5's in batch 4 (design seat R89's E1 rule); k3 is wired: no V
  ( CA=(${LA[@]:0:$n2}); NC=${#CA[@]}; [ -n "$V1" ] && export C4_DEAD_XV="1:$V1"; r86 2 $T0 both k2c; echo "$RK $RW" > $OUT/k2c.res ) &
  ( CA=(${LA[@]:$n2}); NC=${#CA[@]}; [ -n "$V5" ] && export C4_DEAD_XV="4:$V5"; r86 4 $T0 both k4c; echo "$RK $RW" > $OUT/k4c.res ) &
  wait
  read -r K2 W2 < $OUT/k2c.res; read -r K4 W4 < $OUT/k4c.res
  [ -n "${K2:-}" ] && [ -n "${K4:-}" ] || { echo "fail k2 / k4 ($K2 / $K4)" > $OUT/LEFT.res; return; }
  cmerge $OUT/W24.npz $T0 "5-14:$W2 15-24:$W4" || { echo "fail merge k2 / k4" > $OUT/LEFT.res; return; }
  echo "ok $OUT/W24.npz" > $OUT/k4.res   # holder 7 (v5): batch 4 is final here (k1 / k3 plain toggle nothing): k5 may start
  ( CA=(${LA[@]:0:4}); S5T_SECS=2400 CORES="${CA[*]}" bash $T/chem_c3_s5one.sh $OUT/W24.npz 1 $OUT/k1c > /dev/null 2>&1 ) &
  local P1=$!
  ( CA=(${LA[@]:4}); NC=${#CA[@]}; plain25 3 $OUT/W24.npz k3w; echo "$RK" > $OUT/k3w.res ) &
  wait
  local K3; K3=$(cat $OUT/k3w.res 2>/dev/null)
  log "left10: k2 $K2; k4 $K4; k1 plain $(cat $OUT/k1c/DONE 2>/dev/null); k3 ${K3:-none}"
  grep -q ok $OUT/k1c/DONE 2>/dev/null && [ -n "$K3" ] || { echo "fail k1 / k3" > $OUT/LEFT.res; return; }
  { echo "1=dir:$OUT/k1c/b1"; echo "2=$K2"; echo "3=$K3"; echo "4=$K4"; } >> $OUT/parts_L.txt
  echo "ok $OUT/W24.npz" > $OUT/LEFT.res
}

left9f() {  # holder 7 (v5, FAST nine-wired; y46a 08:28Z: left9's k2 -> k3 V route -> k4 chain took 32 min and k5 waited
  # for it, 49 min in all): k2 'both' (k1's and k3's V kept clean) || k4 'both' (k3's and k5's V kept clean) at once
  # (disjoint batches 1-2 / 3-4), then k1 plain || k3 plain (kept, no toggles: 251, +1 a pass against its S5T+);
  # batch 4 is final after the merge, so k5 starts then (OUT/k4.res)
  : > $OUT/parts_L.txt
  local V1 V3 V5 K2 W2 K4 W4 n2 x
  V1=$(vk 1 $T0 c1); V3=$(vk 3 $T0 c3); V5=$(vk 5 $T0 c5)
  log "left9f: k1's V ${V1:-none}, k3's V ${V3:-none}, k5's V ${V5:-none}"
  n2=$(( 2 * ${#LA[@]} / 3 ))
  ( CA=(${LA[@]:0:$n2}); NC=${#CA[@]}; x=""; [ -n "$V1" ] && x="1:$V1"; [ -n "$V3" ] && x="${x:+$x,}2:$V3"
    [ -n "$x" ] && export C4_DEAD_XV="$x"; r86 2 $T0 both k2c; echo "$RK $RW" > $OUT/k2c.res ) &
  ( CA=(${LA[@]:$n2}); NC=${#CA[@]}; x=""; [ -n "$V3" ] && x="3:$V3"; [ -n "$V5" ] && x="${x:+$x,}4:$V5"
    [ -n "$x" ] && export C4_DEAD_XV="$x"; r86 4 $T0 both k4c; echo "$RK $RW" > $OUT/k4c.res ) &
  wait
  read -r K2 W2 < $OUT/k2c.res; read -r K4 W4 < $OUT/k4c.res
  [ -n "${K2:-}" ] && [ -n "${K4:-}" ] || { echo "fail k2 / k4 ($K2 / $K4)" > $OUT/LEFT.res; return; }
  cmerge $OUT/W24.npz $T0 "5-14:$W2 15-24:$W4" || { echo "fail merge k2 / k4" > $OUT/LEFT.res; return; }
  echo "ok $OUT/W24.npz" > $OUT/k4.res
  ( CA=(${LA[@]:0:4}); NC=4; plaink 1 $OUT/W24.npz k1c ) &
  ( CA=(${LA[@]:4:4}); NC=4; plaink 3 $OUT/W24.npz k3c ) &
  wait
  log "left9f: k2 $K2; k4 $K4; k1 $(cat $OUT/k1c.part 2>/dev/null); k3 $(cat $OUT/k3c.part 2>/dev/null)"
  [ -s $OUT/k1c.part ] && [ -s $OUT/k3c.part ] || { echo "fail k1 / k3 plain" > $OUT/LEFT.res; return; }
  { cat $OUT/k1c.part; echo "2=$K2"; cat $OUT/k3c.part; echo "4=$K4"; } >> $OUT/parts_L.txt
  echo "ok $OUT/W24.npz" > $OUT/LEFT.res
}

ngof() {  # ngof K: generic outputs of loader K on T0 (30 kept, 25 wired, 24 / 23 wired with one / two bit-1 relations)
  [ -s $OUT/ng.txt ] || ( cd ~/Quantum/tools && $PY -I - $T0 "$HOME/Quantum/tools" > $OUT/ng.txt 2>> $OUT/port.log <<'EOP'
import sys
sys.path.insert(0, sys.argv[2])   # -I: the tools folder is not on the path by itself (loader2 does the same)
import numpy as np
import chem_loaders_core as core
import chem_ld2_e19grp as eg
core.J, core.BROT = 5, 7
w = np.asarray(np.load(sys.argv[1])["angle_words"], np.int64)
for k in range(0, 15):   # holder 7 (v5): k0 too (a k0-wired table: 25 generic outputs, built here by loader2)
    _, wp = eg.block_set(w, k)
    print(k, int(core.analyze(wp, k, 4).X0.shape[1]))
EOP
  )
  awk -v k=$1 '$1 == k {print $2}' $OUT/ng.txt
}
plainw() {  # plainw K TABLE TAG: any wired loader by loader2 with C4_DROP and NO dead words (no toggles): 25 outputs 228,
  # 24 about 223-224, 23 about 220; sets RK
  local k=$1 tb=$2 w=$OUT/$3
  RK=""
  ( cd ~/Quantum/tools && export CHEM_HB=77 PY=$PY C4_DROP=1 MAXP=${P25_MAXP:-48} SAMPLES=150 CL="${CA[*]:0:8}" && \
      bash sh/chem_c4_loader2.sh $k 0 0 $tb $w ) > /dev/null 2>&1
  if [ -s $w/L$k/best/angle_k_$k.kmx ]; then RK=$w/L$k/best/angle_k_$k.kmx
  elif [ -s $w/L$k/best/angle_k_${k}_t1.kmx ]; then RK=$w/L$k/best/angle_k_${k}_t1.kmx; fi
  log "k$k plain wired: $(grep -h -E 'SUCCESS|end without|kept|no plan' $w/L$k/run.log 2>/dev/null | tail -n 1 | cut -c1-110) -> ${RK:-none}"
}
plaink() {  # plaink K TABLE TAG: loader K with NO toggles, by its kind: kept (30 outputs) plain S5T, wired plainw.
  # Writes the part line to OUT/TAG.part (K=dir:... or K=<kmx>); empty on failure
  local k=$1 tb=$2 tag=$3 ng; ng=$(ngof $k); : > $OUT/$tag.part
  if [ "$ng" = 30 ]; then
    S5T_SECS=2400 CORES="${CA[*]:0:4}" bash $T/chem_c3_s5one.sh $tb $k $OUT/$tag > /dev/null 2>&1
    grep -q ok $OUT/$tag/DONE 2>/dev/null && echo "$k=dir:$OUT/$tag/b$k" > $OUT/$tag.part
  else
    plainw $k $tb $tag; [ -n "$RK" ] && echo "$k=$RK" > $OUT/$tag.part
  fi
  log "k$k plain ($ng outputs): $(cat $OUT/$tag.part)"
}
wiredown() {  # wiredown K TABLE TAG: a WIRED loader by its form on its OWN batch, by its kind: 25 outputs form25 (K6-4),
  # 24 R86's form own (row24b fallback), 23 R89's form (dead23 fallback); sets RK, RW
  local k=$1 tb=$2 tag=$3 ng; ng=$(ngof $k)
  case "$ng" in
    25) form25 $k $tb $tag;;
    24) r86 $k $tb own $tag;;
    23) FALLBACK=1 bash $T/r89_chem_cascade_port.sh $k $tb $OUT/$tag/L$k "${CA[*]:0:8}" > $OUT/$tag.log 2>&1
        RK=$OUT/$tag/L$k/best/angle_k_$k.kmx; RW=$OUT/$tag/L$k/best/table_dead.npz; [ -s $RW ] || RW=$tb; [ -s $RK ] || RK="";;
    *) RK=""; RW="";;
  esac
  log "k$k wired own ($ng outputs): ${RK:-none}"
}
anneal42() {  # anneal42 PAIRDIR: C2's 42-annealer on group 3 of that pair's b1 build (as loader2 runs it), 4 seeds x ASECS
  # on the first four cores; prints the first plan file written (empty if none)
  local p=$1 r
  [ -s $p/dumpb1.json ] || return 0
  ( cd ~/Quantum/tools && for r in 0 1 2 3; do nice -n 10 taskset -c ${CA[$((r % NC))]} $PY -I chem_ld2_e19grp.py anneal $p/dumpb1.json \
      --g 2 --P 42 --secs ${ASECS:-300} --every 20000 --price 40 --rseed $r --save $p/an_$r.json > $p/an_$r.log 2>&1 & done; wait
    for r in 0 1 2 3; do [ -s $p/an_$r.json ] || continue
      $PY -I chem_ld2_e19grp.py plan $p/dumpb1.json $p/an_$r.json --g 2 --tries 20000 --out $p/pl42_$r.json > $p/pl42_$r.log 2>&1
      done ) > /dev/null 2>&1
  ls $p/pl42_*.json 2>/dev/null | tr '\n' ' '
}
form25() {  # form25 K TABLE TAG (holder 6, v2 06:12Z): a 25-output wired loader by design seat R86's form on its last group,
  # OWN batch (K6-4: k3 227 against plain 228 [M]). v1 waited for loader2's annealing loop (5 min a pair, every pair until
  # one reached 154: K6-6's k3 stalled there) before the census. Now: (1) loader2's pairs and plain builds WITHOUT its
  # annealing loop (target 155); (2) R86's census on every pre-checked pair; (3) on up to three deliverable pairs:
  # C2's 42-annealer on group 3, then the form build (227); (4) else the plain loader annealed on one pair (154 = 228),
  # else the 155. Sets RK, RW (RW = its table_dead.npz: toggles in batch K only).
  local k=$1 tb=$2 w=$OUT/$3 p d n b pl i=0 nb=0 best RKP
  RK=""; RW=""
  ( cd ~/Quantum/tools && export CHEM_HB=77 PY=$PY C4_DROP=1 MAXP=${P25_MAXP:-48} SAMPLES=150 CL="${CA[*]:0:8}" C4_T2=155 && \
      bash sh/chem_c4_loader2.sh $k 0 0 $tb $w ) > /dev/null 2>&1
  RKP=$w/L$k/best/angle_k_$k.kmx; [ -s $RKP ] || RKP=""
  log "k$k form25: plain pairs $(grep -l WROTE $w/L$k/pairs/p*/pre.log 2>/dev/null | wc -l) pre-checked; plain $(grep -h -o 'CLOSED [0-9]* products' $w/L$k/run.log 2>/dev/null | tail -n 1)"
  for p in $(grep -l WROTE $w/L$k/pairs/p*/pre.log 2>/dev/null | xargs -r -n1 dirname); do
    ( R86_BATCHES=own CENSUS_ONLY=1 bash $T/chem_c4_r86go25.sh $k $tb $p $w/s/$(basename $p) ${CA[$((i % NC))]} > /dev/null 2>&1 ) &
    i=$((i + 1)); [ $((i % NC)) = 0 ] && wait
  done
  wait
  best=$(for d in $w/s/p*; do n=$(sed -n 's/.*census: \([0-9]*\) \([0-9]*\) .*/\2/p' $d/r86go.log 2>/dev/null | head -n 1); \
    [ -n "$n" ] && [ "$n" -gt 0 ] && echo "$n $(basename $d)"; done | sort -n -r | awk '{print $2}')
  log "k$k form25 deliverable pairs (own): $(echo $best | cut -c1-200)"
  # holder 6 (06:30Z, K6-6): an annealed 42-plan does not always close in the R86 build ('no closed solve for group ...'
  # on k2's words; and with the annealed 42 in place the census can read deficit 1 where the 43 read deliverable): every
  # plan of a pair is tried (R86's dump + census + loader with that plan; a failing census ends in about 30 s), up to
  # three pairs, so that the left half stays shorter than the right chain even in the worst case
  for b in $best; do
    p=$w/L$k/pairs/$b; pls=$(ls $p/pl42_*.json 2>/dev/null | tr '\n' ' ')
    [ -n "$pls" ] || pls=$(anneal42 $p)
    for pl in $pls; do
      d=$w/b/${b}_$(basename $pl .json)
      PL42=$pl R86_BATCHES=own bash $T/chem_c4_r86go25.sh $k $tb $p $d ${CA[0]} > /dev/null 2>&1
      log "k$k form25 $b $(basename $pl): $(tail -n 1 $d/r86go.log 2>/dev/null | cut -c1-100)"
      if grep -q "DONE emit_c4 toffoli 227" $d/r86go.log 2>/dev/null; then
        RK=$d/angle_k_$k.kmx; RW=$d/table_dead.npz; [ -s $RW ] || RW=$tb; break; fi
    done
    [ -n "$RK" ] && break
    nb=$((nb + 1)); [ $nb -ge 3 ] && break
  done
  if [ -z "$RK" ]; then   # holder 6 (06:35Z, K6-6): no form: the TESTED plain path, loader2 with its own annealing loop
    # over the pairs until one closes 154 (= 228: K6-1, K6-4); the 155 of step (1) only if that fails too
    plain25 $k $tb ${3}p; [ -n "$RK" ] || RK=$RKP; RW=$tb
  fi
  log "k$k form25: ${RK:-none}"
}

left10f() {  # the ten-wired G10P57x left with k3's form (holder 6, 05:55Z): k2 'both' -> k1 plain || k3 by form25 on its
  # own batch 3 -> k4 by R86's form on its own batch 4 through R89's supply-aware strict pairs (k5's V kept clean),
  # fallback my row24b (own batch). Batches: k2 1-2, k3 3, k4 4.
  : > $OUT/parts_L.txt
  local V1 V5 K2 W2 K3 W3 K4 W4 s
  V1=$(vk 1 $T0 c1); log "left10f: k1's V ${V1:-none}"
  if [ -s "${K2RES:-}" ]; then cp $K2RES $OUT/k2c.res; log "left10f: k2 taken from $K2RES (test)"   # test only
  else ( CA=(${LA[@]}); NC=${#CA[@]}; [ -n "$V1" ] && export C4_DEAD_XV="1:$V1"; r86 2 $T0 both k2c; echo "$RK $RW" > $OUT/k2c.res ); fi
  read -r K2 W2 < $OUT/k2c.res
  [ -n "${K2:-}" ] || { echo "fail k2" > $OUT/LEFT.res; return; }
  ( CA=(${LA[@]:0:4}); NC=4; plaink 1 $W2 k1c ) &
  ( CA=(${LA[@]:4}); NC=${#CA[@]}; wiredown 3 $W2 k3f; [ -n "$RK" ] || { plainw 3 $W2 k3p; RW=$W2; }; echo "$RK $RW" > $OUT/k3f.res ) &
  wait
  read -r K3 W3 < $OUT/k3f.res
  log "left10f: k2 $K2; k1 $(cat $OUT/k1c.part 2>/dev/null); k3 ${K3:-none}"
  [ -s $OUT/k1c.part ] && [ -n "${K3:-}" ] || { echo "fail k1 / k3" > $OUT/LEFT.res; return; }
  V5=$(vk 5 $W3 c5); [ -n "$V5" ] && export C4_DEAD_XV="4:$V5"
  bash $T/r89_chem_cascade_r86s.sh 4 $W3 $OUT/k4s "${LA[*]:0:8}" 8 > /dev/null 2>&1; unset C4_DEAD_XV
  s=$(grep -o "SUCCESS [^ ]*" $OUT/k4s/r86s.log 2>/dev/null | head -n 1 | cut -d" " -f2)
  if [ -n "$s" ]; then K4=$s/angle_k_4.kmx; W4=$s/table_dead.npz; log "left10f: k4 own (R89 supply-aware): $K4"
  else CA=(${LA[@]}); chain1 4 row24b $W3 R4c; K4=$RK; W4=$RW; fi
  [ -n "$K4" ] || { echo "fail k4" > $OUT/LEFT.res; return; }
  { cat $OUT/k1c.part; echo "2=$K2"; echo "3=$K3"; echo "4=$K4"; } >> $OUT/parts_L.txt
  echo "ok $W4" > $OUT/k4.res; echo "ok $W4" > $OUT/LEFT.res
}

left9() {  # the nine-wired left cascade of port3 (k2 'both' -> k1 plain || k3 S5T+ by the V route -> k4 row24b)
  : > $OUT/parts_L.txt
  local V1 V3 K2 W2 K3 W3 K4 W4
  V1=$(vk 1 $T0 c1); V3=$(vk 3 $T0 c3); [ -n "$V1" ] && [ -n "$V3" ] && export C4_DEAD_XV="1:$V1,2:$V3"
  log "left9: k1's V ${V1:-none}, k3's V ${V3:-none}"
  r86 2 $T0 both k2c; unset C4_DEAD_XV; K2=$RK; W2=$RW
  [ -n "$K2" ] || { echo "fail k2" > $OUT/LEFT.res; return; }
  ( CA=(${LA[@]:8:4}); [ ${#CA[@]} -ge 2 ] || CA=(${LA[@]:0:4}); NC=${#CA[@]}; plaink 1 $W2 k1c ) &
  local P1=$!
  vroute 3 $W2 k3c; K3=$VK; W3=$VW
  if [ -z "$K3" ]; then   # holder 7 (13:45Z, y46aR's exact port: no V passed phase 1, no fallback, the port failed): k3 plain
    plaink 3 $W2 k3p0; [ -s $OUT/k3p0.part ] && { K3=$(cut -d= -f2- $OUT/k3p0.part); W3=$W2; log "left9: k3 plain (its V route failed)"; }
  fi
  wait $P1
  log "left9: k2 $K2; k1 $(cat $OUT/k1c.part 2>/dev/null); k3 ${K3:-none}"
  [ -s $OUT/k1c.part ] && [ -n "$K3" ] || { echo "fail k1 / k3" > $OUT/LEFT.res; return; }
  # holder 6: k4 by R86's form on its OWN batch 4 through design seat R89's supply-aware strict pairs (its E1: 221 against
  # row24b's 222), k5's V kept clean; fallback my row24b
  local V5 s; V5=$(vk 5 $W3 c5); [ -n "$V5" ] && export C4_DEAD_XV="4:$V5"
  bash $T/r89_chem_cascade_r86s.sh 4 $W3 $OUT/k4s "${CA[*]:0:8}" 8 > /dev/null 2>&1; unset C4_DEAD_XV
  s=$(grep -o "SUCCESS [^ ]*" $OUT/k4s/r86s.log 2>/dev/null | head -n 1 | cut -d" " -f2)
  if [ -n "$s" ]; then K4=$s/angle_k_4.kmx; W4=$s/table_dead.npz; log "left9: k4 own (R89 supply-aware): $K4"
  else chain1 4 row24b $W3 R4c; K4=$RK; W4=$RW; fi
  [ -n "$K4" ] || { echo "fail k4" > $OUT/LEFT.res; return; }
  { cat $OUT/k1c.part; echo "2=$K2"; echo "3=$K3"; echo "4=$K4"; } >> $OUT/parts_L.txt
  echo "ok $W4" > $OUT/k4.res; echo "ok $W4" > $OUT/LEFT.res
}

eleven() {  # G11w11 (wires only): the kept k1 / k5 / k7 by S5T+ (C3's V route; batches 0-1 / 4-5 / 6-7 are disjoint, so
  # they run at once when there are 8 cores each, else one after another), then EVERY wired loader by loader2 (25 outputs,
  # no toggles: all on the same words, 4 cores each, as many at once as the cores allow); a kept loader whose V route
  # fails gets plain S5T last (no toggles). Every kept loader's neighbour carries no bit-1 relation here (relnbr: empty).
  : > $OUT/parts.txt
  local all=(${LA[@]} ${RA[@]}) k i=0 cur=$T0 m nv ns
  nv=$(( ${#all[@]} / 8 )); [ $nv -lt 1 ] && nv=1
  ns=$(( ${#all[@]} / 4 )); [ $ns -lt 1 ] && ns=1
  log "eleven: ${#all[@]} cores: $nv V routes at once, $ns wired loaders at once"
  for k in 1 5 7; do
    ( CA=(${all[@]:$((8 * (i % nv))):8}); NC=${#CA[@]}; vroute $k $T0 k${k}v; echo "$VK $VW" > $OUT/k${k}v.res ) &
    i=$((i + 1)); [ $((i % nv)) = 0 ] && wait
  done
  wait
  m=""
  for k in 1 5 7; do VK=""; VW=""; read -r VK VW < $OUT/k${k}v.res
    if [ -n "${VK:-}" ]; then echo "$k=$VK" >> $OUT/parts.txt; m="$m $((5 * k))-$((5 * k + 4)):$VW"
    else echo "$k=PLAIN" >> $OUT/parts.txt; fi; done
  if [ -n "$m" ]; then cmerge $OUT/W11.npz $T0 "${m# }" || { echo "fail merge 1 / 5 / 7" > $OUT/DONE; exit 1; }; cur=$OUT/W11.npz; fi
  log "eleven: kept loaders $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-300)"
  i=0
  for k in $(echo $WI | tr ',' ' '); do
    ( CA=(${all[@]:$((4 * (i % ns))):4}); NC=${#CA[@]}; plain25 $k $cur w$k; echo "$k=$RK" > $OUT/w$k.res ) &
    i=$((i + 1)); [ $((i % ns)) = 0 ] && wait
  done
  wait
  for k in $(echo $WI | tr ',' ' '); do grep -q '=/' $OUT/w$k.res 2>/dev/null && cat $OUT/w$k.res >> $OUT/parts.txt; done
  for k in 1 5 7; do grep -q "^$k=PLAIN" $OUT/parts.txt || continue
    sed -i "/^$k=PLAIN/d" $OUT/parts.txt
    ( CA=(${all[@]:0:4}); S5T_SECS=2400 CORES="${CA[*]}" bash $T/chem_c3_s5one.sh $cur $k $OUT/k${k}p > /dev/null 2>&1 )
    log "eleven: k$k plain S5T $(cat $OUT/k${k}p/DONE 2>/dev/null)"
    grep -q ok $OUT/k${k}p/DONE 2>/dev/null && echo "$k=dir:$OUT/k${k}p/b$k" >> $OUT/parts.txt
  done
  cp $cur $OUT/T.npz; echo "5-74" > $OUT/COLS
  n=$(grep -c '=' $OUT/parts.txt)
  if [ $n -ge 14 ]; then echo ok > $OUT/DONE; else echo "fail: $n of 14 parts" > $OUT/DONE; fi
  log "done (ELEVEN): $(tr '\n' ' ' < $OUT/parts.txt | cut -c1-600)"
}

RA=($RC); LA=($LC)
log "start: T0 $T0, LAYOUT ${LAYOUT:-9}, WI $WI, right cores $RC, left cores $LC, R86GO $(basename $R86GO)"
ngof 1 > /dev/null; log "generic outputs k1..k14: $(awk '{print $2}' $OUT/ng.txt | tr '
' ' ')"
if [ "${LAYOUT:-9}" = 11 ]; then eleven; exit 0; fi
if [ "${LEFTONLY:-0}" = 1 ]; then   # test only: the left half alone
  CA=(${LA[@]}); NC=${#CA[@]}; if [ "$(ngof 3)" != 30 ] && [ "${K3FORM:-1}" = 1 ]; then left10f; elif [ "$(ngof 3)" != 30 ]; then left10; else left9; fi
  log "LEFT (LEFTONLY): $(cat $OUT/LEFT.res 2>/dev/null)"; echo "left $(cat $OUT/LEFT.res 2>/dev/null)" > $OUT/DONE; exit 0
fi
P5=""; P0=""
if [ "${FAST:-0}" = 1 ]; then   # holder 7: the split right half and the parallel lefts (left10 ten-wired, left9f nine-wired)
  # K7-1 (08:15Z): R86's census over all 48 pairs took 14.6 min (k9, 6 cores) and 15 min (k2, 4 cores) = the long pole
  # of every R86 step; FAST stops it after the first wave with 3 deliverable pairs (R86_CMIN=0 restores port4's rule),
  # and k2 (the slow one of left10's pair: 28 min against k4's 7.3 on 4 cores each) gets two thirds of the left cores
  export R86_CMIN=${R86_CMIN:-3} L10K2=${L10K2:-2/3}
  LF=left9f; [ "${LEFT9F:-1}" = 1 ] || LF=left9; [ "$(ngof 3)" != 30 ] && LF=left10
  # holder 7 (15:20Z, the lead's 15:16Z: T11's r10ag2ca = ten-wired + ONE bit-1 relation on k3, so k3 has 24 generic
  # outputs and is worth R86's form (221) where left10's plain build gives about 223): left10f, whose wiredown builds
  # a 24-output k3 by R86's form on its own batch, then k4 own (R89 supply-aware) - sequential, inside the right's time
  [ "$(ngof 3)" = 24 ] && LF=left10f
  log "FAST: right_split on ${RA[*]}; $LF on ${LA[*]}; R86_CMIN $R86_CMIN"
  # holder 7 (v5): a k0-WIRED table (k0 with 25 generic outputs: T11's layout) gets k0 by loader2 HERE, at once, on the
  # cores K0C (batch 0 is read by k0 and k1 only and toggled by nobody: k1 is plain), so its entry needs no k0 pipeline
  if [ "$(ngof 0)" != 30 ]; then
    K0A=(${K0C:-${LA[@]: -8}})
    ( CA=(${K0A[@]}); NC=${#CA[@]}; plain25 0 $T0 k0w; [ -n "$RK" ] && echo "0=$RK" > $OUT/k0w.part ) &
    P0=$!; log "k0 WIRED ($(ngof 0) generic outputs): loader2 on ${K0A[*]}"
  fi
  ( CA=(${RA[@]}); NC=${#CA[@]}; right_split ) &
  PR=$!
  ( CA=(${LA[@]}); NC=${#CA[@]}; $LF ) &
  PL=$!
  # holder 7 (v5, y46a: k5 waited for the whole left, 17 min after it): k5 the minute batches 4 and 5 are final (k4.res
  # from the left, k6.res from S1), on the last 4 left cores (left10 / left9f leave them idle after their k2 / k4 phase)
  while { [ ! -s $OUT/k4.res ] || [ ! -s $OUT/k6.res ]; } && kill -0 $PL 2>/dev/null; do sleep 20; done
  if [ -s $OUT/k4.res ] && [ -s $OUT/k6.res ]; then
    read -r s4 W4F < $OUT/k4.res; read -r s6 W6F < $OUT/k6.res
    if [ "$s4" = ok ] && [ "$s6" = ok ] && cmerge $OUT/K5IN.npz $T0 "5-24:$W4F 25-34:$W6F"; then
      ( CA=(${LA[@]: -4}); NC=4; plaink 5 $OUT/K5IN.npz k5c ) &
      P5=$!; log "k5 plain started early on [k4's batch 4 + k6's batch 5]"
    fi
  fi
else
( CA=(${RA[@]}); NC=${#CA[@]}; right_chain ) &
PR=$!
# holder 6 (06:20Z): if the form left (left10f) fails, the ten-wired left of test K6-1 (left10, measured on the projection) runs
( CA=(${LA[@]}); NC=${#CA[@]}
  if [ "$(ngof 3)" != 30 ] && [ "${K3FORM:-1}" = 1 ]; then left10f
    read -r st10 rest10 < $OUT/LEFT.res; [ "$st10" = ok ] || { log "left10f failed ($(cat $OUT/LEFT.res)): the plain ten-wired left"; left10; }
  elif [ "$(ngof 3)" != 30 ]; then left10; else left9; fi ) &
PL=$!
fi
wait $PL; log "LEFT: $(cat $OUT/LEFT.res 2>/dev/null)"
# holder 6 (06:00Z): k5 plain EARLY, on the left cores, the minute LEFT and k6 are done (k5 reads batches 4 and 5 only;
# batch 5 is final after k6 'both': no later loader of the right chain touches it), while the right chain goes on
read -r sl WL < $OUT/LEFT.res
if [ "$sl" = ok ] && [ -z "$P5" ]; then
  while [ ! -s $OUT/k6.res ] && [ ! -s $OUT/RIGHT.res ]; do sleep 30; done
  if [ -s $OUT/k6.res ]; then
    read -r s6 W6F < $OUT/k6.res
    if cmerge $OUT/K5IN.npz $T0 "5-24:$WL 25-34:$W6F"; then
      ( CA=(${LA[@]}); NC=${#CA[@]}; plaink 5 $OUT/K5IN.npz k5c ) &
      P5=$!; log "k5 plain started early on [LEFT's batch 4 + k6's batch 5]"
    fi
  fi
fi
wait $PR; log "RIGHT: $(cat $OUT/RIGHT.res 2>/dev/null)"
read -r sr WR < $OUT/RIGHT.res
[ "$sl" = ok ] && [ "$sr" = ok ] || { [ -n "$P5" ] && kill $P5 2>/dev/null; echo "fail: LEFT $(cat $OUT/LEFT.res) RIGHT $(cat $OUT/RIGHT.res)" > $OUT/DONE; exit 1; }
cmerge $OUT/T5.npz $T0 "5-24:$WL 25-74:$WR" || { echo "fail merge LEFT / RIGHT" > $OUT/DONE; exit 1; }
if [ -n "$P5" ]; then
  wait $P5
  [ -s $OUT/k5c.part ] || { log "k5 early failed: rebuilt on the final words"; rm -rf $OUT/k5c $OUT/k5c.part; P5=""; }
fi
if [ -n "$P5" ]; then
  # the early k5 is valid only if the final batches 4-5 equal its input's (they must: nothing after k6 touches them)
  $PY - $OUT/K5IN.npz $OUT/T5.npz <<'EOP' >> $OUT/port.log 2>&1 || { log "k5 early: batches 4-5 moved: rebuilt"; rm -rf $OUT/k5c $OUT/k5c.part; P5=""; }
import sys, numpy as np
a = np.asarray(np.load(sys.argv[1])["angle_words"]); b = np.asarray(np.load(sys.argv[2])["angle_words"])
assert (a[:, 20:30] == b[:, 20:30]).all(), "batches 4-5 differ"
print("k5 early: batches 4-5 equal to the final words")
EOP
fi
if [ -z "$P5" ]; then CA=(${RA[@]}); NC=${#CA[@]}; plaink 5 $OUT/T5.npz k5c; fi
log "k5: $(cat $OUT/k5c.part 2>/dev/null)"
[ -s $OUT/k5c.part ] || { echo "fail k5" > $OUT/DONE; exit 1; }
if [ -n "$P0" ]; then   # holder 7 (v5): the k0-wired table's own k0 (batch 0 untouched by every lever: checked below)
  wait $P0; log "k0 wired: $(cat $OUT/k0w.part 2>/dev/null)"
  [ -s $OUT/k0w.part ] || { echo "fail k0 wired" > $OUT/DONE; exit 1; }
  $PY - $T0 $OUT/T5.npz <<'EOP' >> $OUT/port.log 2>&1 || { echo "fail k0 wired: batch 0 moved" > $OUT/DONE; exit 1; }
import sys, numpy as np
a = np.asarray(np.load(sys.argv[1])["angle_words"]); b = np.asarray(np.load(sys.argv[2])["angle_words"])
assert (a[:, 0:5] == b[:, 0:5]).all(), "batch 0 differs"
print("k0 wired: batch 0 equal to the final words")
EOP
fi
{ [ -s $OUT/k0w.part ] && cat $OUT/k0w.part; cat $OUT/parts_L.txt; cat $OUT/k5c.part; cat $OUT/parts_R.txt; } > $OUT/parts.txt
cp $OUT/T5.npz $OUT/T.npz; echo "5-74" > $OUT/COLS
echo ok > $OUT/DONE; log "done (LAYOUT ${LAYOUT:-9}): $(tr '
' ' ' < $OUT/parts.txt | cut -c1-600)"
