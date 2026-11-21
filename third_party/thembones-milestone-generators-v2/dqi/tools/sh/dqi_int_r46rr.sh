#!/usr/bin/env bash
# dqi_int_r46rr.sh -- B5 (Nov 15 18:55Z; runs ON box 6, detached): re-rolls (tools/dqi_int_pptag.py: one INFO_TAG, no
# gate change) of ~/qcc/ws/dqi_int/r46/<name>/dqi.kmb until <want> of them read 0 / 0 on both official salts
# (tools/qcc_verify.py, 102,400 shots each); an unclean re-roll's kmb is deleted at once. Lines in rr_<name>.txt.
#   bash dqi_int_r46rr.sh <name> <cores> <workers> <want> <tag> [<tag> ...]
set -u
N=$1; C=$2; W=$3; WANT=$4; shift 4
D=~/qcc/ws/dqi_int; O=$D/r46/$N; T=~/Quantum/tools; PW=~/qcc/.venv/bin/python; G=$D/rr_$N.txt
clean=0
for tg in "$@"; do
  [ $clean -ge $WANT ] && break
  R=$O/rr_$tg; mkdir -p $R
  (cd $T && nice -n 10 taskset -c $C $PW dqi_int_pptag.py $O/dqi.kmb $R/dqi.kmb r46$N$tg > $R/tag.log 2>&1)
  line="$tg $(sha256sum $R/dqi.kmb | cut -c1-16)"
  ok=1
  for salt in "" alt-salt-1; do
    tag=${salt:-salt0}
    r=$(cd $T && nice -n 10 taskset -c $C $PW qcc_verify.py $R --tasks dqi --salt "$salt" --workers $W --json $R/verify_$tag.json 2>&1 | grep -a "cost" | tail -1)
    line="$line | $tag $r"
    echo "$r" | grep -q "errors 0/0" || { ok=0; break; }
  done
  echo "$line" >> $G
  if [ $ok = 1 ]; then clean=$((clean+1)); else rm -f $R/dqi.kmb; fi
done
echo "clean $clean $(date -u +%H:%MZ)" >> $G
