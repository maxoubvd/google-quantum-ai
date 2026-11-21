#!/bin/bash
# chem_own: pick the cheapest stored plan of every circuit among several search runs, then run the ONE-COMMAND
# build (tools/sh/chem_own_build.sh) on it: circuits -> every address checked -> official verifier on salt ""
# and alt-salt-1 -> published in ~/qcc/ws/chem_own/best_<LABEL>/ with MANIFEST.json and plan_spec.json.
#   LABEL    entry label (required)
#   DIRS     plan directories, e.g. "plans/r1 box_out/s1 box_out/s2" (relative to ~/qcc/ws/chem_own)
#   TABLE    table npz (default: ref_s007c/submission.npz = public words + our three-stage refit, 46.911221)
#   WEIGHTS  keep (default: the table's weights are ours already) | refit | scratch
#   WORKERS  processes for the build (default 3)       LOG  run in the background with this log name
# The build command printed in MANIFEST.json, run on a clean work directory, recomputes the same plans
# (search is deterministic in (words, seed, effort, rounds); each circuit's seed and rounds are in plan_spec.json).
set -e
cd ~/qcc/ws/chem_own
T=<HOME>/Quantum/tools
[ -n "$LABEL" ] || { echo "LABEL required"; exit 2; }
TABLE=${TABLE:-$HOME/qcc/ws/chem_own/ref_s007c/submission.npz}
mkdir -p build/$LABEL/plans
venv/bin/python $T/chem_own_pick.py build/$LABEL/plans ${DIRS:-plans/r1}
export ARGS="$TABLE $LABEL --weights ${WEIGHTS:-keep} --spec $HOME/qcc/ws/chem_own/build/$LABEL/plans/plan_spec.json --workers ${WORKERS:-3} $EXTRA"
tr -d '\r' < $T/sh/chem_own_build.sh | bash
