#!/bin/bash
# P1 (sixth holder): post-pass ORDER variant on one RAW DQI build: dead-write (c3) -> P1's exact links pass
# (tools/sh/dqi_pp_exact.sh, unchanged) -> dead-write again.  Dead writes removed first free holder capacity for links.
# Runs ON A BOX.  Usage: bash dqi_pp_dwx.sh <raw dqi.kmb> <fresh out dir> [cores=4-11]
# Output: <out>/final/dqi.kmb (+ paired.txt at each dead-write step; the links step's own report in <out>/links)
set -u
IN=$(readlink -f "$1"); OUT=$2; C=${3:-4-11}; T=~/Quantum/tools
[ -e "$OUT" ] && { echo "$OUT exists"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
bash $T/sh/dqi_pp_deadw.sh $IN dw1 $C 200 || exit 1
tr -d '\r' < $T/sh/dqi_pp_exact.sh > exact_frozen.sh
bash exact_frozen.sh $(readlink -f dw1/dqi.kmb) links $C | tail -n 5
bash $T/sh/dqi_pp_deadw.sh $(readlink -f links/dqi.kmb) final $C 200
echo "dwx done $(date -u +%H:%M)Z"
