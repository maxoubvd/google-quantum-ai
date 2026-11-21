#!/usr/bin/env bash
# New file added by the team (Maxime Bendavid), 2026-11-21, plan step S2.26.
# Rebuilds ThemBones' invsqrt (b196lend, sha256 c639a045...) with their one-command pipeline (tools/f9_i1_stack.sh),
# in the sandbox made by setup_sandbox.sh. Run inside WSL:  bash generators/thembones/build_invsqrt.sh [CORE]
#
# Two inputs the pipeline needs are not in their dataset, so this script supplies them as follows:
#  - fc18_M's emit (~/qcc/ws/invsqrt_fc/s18/files/M/invsqrt.kmx), used only by the pipeline's first all-input check:
#    made here by the same emit command on FC's srcM with cfg_M and R59's angles, without the S5 generator patch.
#  - the FC18-M9 device file s18/cub/m9dev.json (config key "_fc18m9") is missing, as is the srcR tree that cfg_M names:
#    the key is removed (the generator then emits its built-in m9 branch, +0.5 Toffoli a shot) and srcR is replaced by
#    srcM, as the pipeline's own q18 recipe does. The output is therefore NOT byte-identical to c639a045; the pipeline's
#    final step proves it exactly against c639a045 on all 3,133,452 inputs (same phase, x restored, clean scratch).
# Time: about 5 min on one core, 1.8 GB peak (measured 2026-11-21: 301 s, output e2787831...).
set -u
TBH=${TBH:-$HOME/tb}; CORE=${1:-2}
export HOME=$TBH            # their scripts use ~ for ~/qcc and ~/Quantum
Q=$TBH/Quantum; P=$TBH/qcc/.venv/bin/python; S18=$TBH/qcc/ws/invsqrt_fc/s18; R59=$TBH/qcc/ws/invsqrt_fc/best/r59b0
LOG=$TBH/inv_build.log; : > $LOG; T0=$(date +%s)
$P - "$S18/cfg/cfg_M.json" "$S18/srcM" <<'PY'
import json, sys
p, src = sys.argv[1:3]
c = json.load(open(p)); c.pop("_fc18m9", None)
c["fc2_src"] = src; c["h1c"]["cubprog"] = src + "/cub_rot2e6e_acc.json"
json.dump(c, open(p, "w"), indent=1)
print("cfg_M: key _fc18m9 removed (m9dev.json not published); fc2_src ->", src)
PY
rm -rf $TBH/qcc/ws/f9_i1 $S18/files/M; mkdir -p $S18/files/M
O=$S18/files/M; cp $S18/cfg/cfg_M.json $O/cfg.json; cp $R59/angles.json $O/angles.json
echo "== M emit $(date -u +%H:%M:%SZ)" | tee -a $LOG
( export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=7,3,0;9,3,1" R41_XC=DIV1 FXACCCAP=30 \
    NUMBA_CACHE_DIR=/tmp/tbnb_M R55_NET=$R59/router.json I39_CAPS= PYTHONDONTWRITEBYTECODE=1 S5LEND=0
  cd /tmp && $P -I $Q/tools/r64_inv_stack_xcrun.py emit $S18/srcM $O $O/angles.json $O/cfg.json >> $LOG 2>&1 )
[ -f $O/invsqrt.kmx ] || { echo "M emit failed, see $LOG"; exit 1; }
cd $TBH/ds/invsqrt
echo "== f9_i1_stack $(date -u +%H:%M:%SZ)" | tee -a $LOG
bash $Q/tools/f9_i1_stack.sh q18 dq18c prog/dq18c.staged.json $CORE 2>&1 | tee -a $LOG
B=$TBH/qcc/ws/f9_i1/best/dq18c/invsqrt.kmb
echo "== wall $(( $(date +%s) - T0 )) s; output $(sha256sum $B 2>/dev/null | cut -c1-64)" | tee -a $LOG
