#!/usr/bin/env bash
# f9_i1_h6vloc.sh [DIR] (seat I1, sixth holder, LOCAL WSL): THE LOCAL FALLBACK FOR THE VERIFICATION LANE (box 4 does not
# answer since about 10:55Z). Runs the team's wrapper of the official verifier (tools/qcc_verify.py: the kit in
# ~/qcc/data/verifiers, the same script box_verify.sh runs on box 4) on THIS PC, both salts one after the other, 2 workers on
# the last two cores, nice 10. Default DIR = a CONTROL: a copy of b196lend's invsqrt.kmb (sha c639a045...) in
# ~/qcc/ws/f9_i1/vtest/b196copy (b196lend's own folder is only read). Expected on the control: 63,055.04 / 63,055.06.
# JSON files land in DIR as vloc_salt0.json / vloc_alt-salt-1.json. Log ~/qcc/ws/f9_i1/logs/vloc_<dir>.log.
Q=<HOME>/Quantum; W=~/qcc/ws/f9_i1; P=~/qcc/.venv/bin/python
cd $Q || exit 1
D=${1:-$W/vtest/b196copy}
if [ -z "${1:-}" ]; then
  mkdir -p $D; cp ~/qcc/ws/s5/inv_h1end/best/b196lend/invsqrt.kmb $D/invsqrt.kmb
fi
D=$(readlink -f $D); L=$W/logs/vloc_$(basename $D).log
N=$(nproc); C="$((N - 2)),$((N - 1))"
{
echo "== vloc $(date -u +%H:%M:%SZ): $D, invsqrt.kmb sha $(sha256sum < $D/invsqrt.kmb | cut -c1-16); nproc $N, cores $C; load $(cut -d' ' -f1-3 /proc/loadavg)"
for salt in "" alt-salt-1; do
  tag=${salt:-salt0}; s=$(date +%s)
  nice -n 10 taskset -c $C $P tools/qcc_verify.py $D --tasks invsqrt --salt "$salt" --workers 2 --json $D/vloc_$tag.json 2>&1 | grep -v '^\s*$' | tail -4
  echo "-- $tag: $(( $(date +%s) - s )) s, ended $(date -u +%H:%M:%SZ)"
done
echo "== vloc end $(date -u +%H:%M:%SZ)"
} > $L 2>&1
cat $L
