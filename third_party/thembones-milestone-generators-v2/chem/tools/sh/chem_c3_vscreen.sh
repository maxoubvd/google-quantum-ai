#!/bin/bash
# Builder C3 (EIGHTH holder, Nov 18 21:00Z). ON BOX 5, DETACHED. S5T+ VARIANT B ON ANOTHER V: waits for a V probe
# (tools/chem_c3_o29probe.py: PROBE/summary.json), screens its 8 largest-padding V's at staircase 149 for phase 1 (one core
# each, P1ONLY, 400 draws x seeds 0-1), builds the first that passes on all CORES (tools/sh/chem_c3_s5tpb.sh), then
# optionally builder C4's chain tool on its words (CHAIN, e.g. "4:row24", read-only) and a plain S5T of loader PLAIN.
# Writes OUT/vscreen.log, OUT/full/s5tp/{angle_k_K.kmx, words.npz}, OUT/parts.txt, OUT/B_table.npz, OUT/DONE.
#   K=3 TB=... PROBE=... CORES="c ..." RELNBR="0:1 ..." WI=2,4,... CHAIN="4:row24" PLAIN=5 OUT=... setsid nohup bash chem_c3_vscreen.sh &
set -u
T=$HOME/Quantum/tools/sh; CO=($CORES)
mkdir -p $OUT; rm -f $OUT/DONE
log() { echo "$(date -u +%H:%M:%SZ) $*" >> $OUT/vscreen.log; }
i=0; while [ ! -s $PROBE/summary.json ] && [ $i -lt 60 ]; do sleep 20; i=$((i + 1)); done
# 21:05Z: ALLV=1 takes the V's whatever the probe's own-basis staircase read (with PFILE = a chain basis that reads 149,
# exported to the s5tpb calls below: pB1pp's k3 reads 150 on its own basis for every V, 149 on the chain's bs2)
VS=$(python3 -c "import json,sys; s=[r for r in json.load(open(sys.argv[1])) if r.get('stageA') == 149 or sys.argv[2] == '1']; s.sort(key=lambda r: -r['size']); print(' '.join(r['V'] for r in s[:8]))" $PROBE/summary.json "${ALLV:-0}")
log "start: k$K, V's $VS"
n=0
for V in $VS; do
  R=$OUT/scr_$V; mkdir -p $R; echo "{\"V\": {\"$K\": \"$V\"}}" > $R/vproj.json
  ( VPROJ=$R/vproj.json P1ONLY=1 DEADBITS=1,2,3,4,5 RELNBR="$RELNBR" P1TRIES=400 P1SEEDS="0 1" P1SECS=900 NONBR=1 \
      CORES="${CO[$((n % ${#CO[@]}))]}" bash $T/chem_c3_s5tpb.sh $R $TB $K $WI > /dev/null 2>&1 ) &
  n=$((n + 1)); [ $((n % ${#CO[@]})) = 0 ] && wait
done
wait
PV=""; PS=""
for V in $VS; do
  st=$(cat $OUT/scr_$V/s5tp/DONE 2>/dev/null); log "V $V: $st"
  [ -z "$PV" ] && [ "${st%% *}" = p1ok ] && { PV=$V; PS=$(echo $st | awk '{print $3}'); }
done
[ -n "$PV" ] || { log "no V passed phase 1"; echo "fail screen" > $OUT/DONE; exit 0; }
log "full plan on V $PV (seed $PS), cores $CORES"
VPROJ=$OUT/scr_$PV/vproj.json DEADBITS=1,2,3,4,5 RELNBR="$RELNBR" P1TRIES=400 P1SEEDS="$PS" P1SECS=1500 NONBR=1 S5TP_SECS=2400 \
  CORES="$CORES" bash $T/chem_c3_s5tpb.sh $OUT/full $TB $K $WI > /dev/null 2>&1
st=$(cat $OUT/full/s5tp/DONE 2>/dev/null); log "k$K full: $st"
[ "${st%% *}" = ok ] || { echo "fail full $st" > $OUT/DONE; exit 0; }
echo "$K=$OUT/full/s5tp/angle_k_$K.kmx" > $OUT/parts.txt; W=$OUT/full/s5tp/words.npz
if [ -n "${CHAIN:-}" ]; then
  CL="$CORES" bash $T/chem_c4_chain.sh $W $OUT/L1 "$CHAIN" > $OUT/L1.log 2>&1
  [ -s $OUT/L1/table_final.npz ] || { log "chain $CHAIN gave no table"; echo "fail chain" > $OUT/DONE; exit 0; }
  log "chain: $(grep -h -o 'k[0-9]* [a-z0-9]*: [0-9:]* SUCCESS.*toffoli [0-9]*' $OUT/L1/chain.log | sed 's/SUCCESS.*toffoli/=/' | tr '\n' ' ')"
  for km in $CHAIN; do k=${km%%:*}; echo "$k=$OUT/L1/L$k/best/angle_k_$k.kmx" >> $OUT/parts.txt; done
  W=$OUT/L1/table_final.npz
fi
if [ -n "${PLAIN:-}" ]; then
  S5T_SECS=2400 CORES="$CORES" bash $T/chem_c3_s5one.sh $W $PLAIN $OUT/k${PLAIN}p > /dev/null 2>&1
  log "k$PLAIN plain: $(cat $OUT/k${PLAIN}p/DONE 2>/dev/null)"; echo "$PLAIN=dir:$OUT/k${PLAIN}p/b$PLAIN" >> $OUT/parts.txt
fi
cp $W $OUT/B_table.npz; echo ok > $OUT/DONE; log "done"
