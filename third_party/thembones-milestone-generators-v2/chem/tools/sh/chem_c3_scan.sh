#!/usr/bin/env bash
# Builder C3: scan (basis file x replaced block x point set D) for ONE loader until the 68-AND kept-product loader
# emits at or under TARGET Toffolis (default 256); every try is one line of OUT/scan_k<K>.txt.
#   bash chem_c3_scan.sh IN K OUT [TARGET] [MAXOK]
# Run it pinned: taskset -c <core> nice -n 10 bash chem_c3_scan.sh ...
set -u
IN=$1; K=$2; OUT=$3; TARGET=${4:-256}; MAXOK=${5:-1}
S=${CT:-<HOME>/Quantum/tools}/sh
mkdir -p "$OUT"
SC=$OUT/scan_k$K.txt
DS="13,14,15 4,10,14 5,10,14 12,13,14 3,9,14 6,11,15 7,11,13 2,9,15"
PFS=$(cd "$IN" && ls k${K}.json k${K}b*.json 2>/dev/null)
nok=0
for PF in $PFS; do
  for D in $DS; do
    for REP in 0 1 2 3 4; do
      L=$(bash $S/chem_c3_one.sh "$IN" "$K" "$PF" "$REP" "$D" "$OUT" 2)
      rc=$?
      echo "$(date -u +%H:%M:%SZ) rc=$rc $L" >> "$SC"
      if [ $rc -eq 0 ]; then
        n=$(echo "$L" | sed -n 's/.*: \([0-9]*\) Toffolis.*/\1/p'); q=$(echo "$L" | sed -n 's/.*qubits \([0-9]*\).*/\1/p')
        if [ -n "$n" ] && [ "$n" -le "$TARGET" ] && [ "$q" -le 128 ]; then
          nok=$((nok + 1))
          [ $nok -ge $MAXOK ] && { echo "$(date -u +%H:%M:%SZ) DONE k$K: $n" >> "$SC"; exit 0; }
        fi
      fi
    done
  done
done
echo "$(date -u +%H:%M:%SZ) END k$K: $nok at or under $TARGET" >> "$SC"
