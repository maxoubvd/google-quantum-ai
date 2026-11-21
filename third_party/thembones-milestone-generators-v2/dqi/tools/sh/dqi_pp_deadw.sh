#!/bin/bash
# P1 (sixth holder): the DEAD-WRITE pass (public rewrite c3) on one DQI file, with its paired-seed check.
# Runs ON A BOX (venvs as in dqi_pp_exact.sh).  Usage: bash dqi_pp_deadw.sh <in dqi.kmb> <fresh out dir> [cores] [batches]
# Output: <out>/dqi.kmb (child), dqi.kmb.sites.txt, paired.txt (CHILD-ONLY must read 0), dqi.kmb.xz
set -u
IN=$(readlink -f "$1"); OUT=$2; C=${3:-4-11}; NB=${4:-200}
T=~/Quantum/tools; PY=~/qcc/.venv_work/bin/python; PV=~/qcc/.venv/bin/python; P="nice -n 10 taskset -c $C"
[ -e "$OUT/dqi.kmb" ] && { echo "$OUT/dqi.kmb exists"; exit 1; }
mkdir -p "$OUT" && cd "$OUT" || exit 1
echo "input $(sha256sum $IN | cut -c1-16) $(date -u +%H:%M)Z"
$P $PY $T/dqi_pp_deadw.py $IN dqi.kmb || exit 1
echo "child $(sha256sum dqi.kmb | cut -c1-16) $(date -u +%H:%M)Z"
$P $PV $T/dqi_pp_paired.py dqi.kmb $IN salt0 8 $NB > paired.txt 2>&1; tail -6 paired.txt
$P xz -T4 -6 -k dqi.kmb && ls -la dqi.kmb.xz
echo "done $(date -u +%H:%M)Z"
