#!/usr/bin/env bash
# FC eighteenth holder, LOCAL WSL: a found q18 END-exit program (GHOST's staged format, inputs q17, w20..w25, at most 3
# temporaries) -> a verified file ON R25c, in one command (I39's tools/sh/inv_i39_q18find.sh, re-based on my srcR).
#   bash tools/sh/invsqrt_fc_s18q18.sh NAME PROG.json
# Steps: copy of srcR with fcsup/sc8_q18.prog.json replaced -> emit -> ALL-input check against R25c's emit (four outcome
# modes) -> the hold pipeline (relabel to 39, P1, post-chain) -> pack fc18_NAME with SOURCE.txt -> verifier, both salts.
set -u
Q=<HOME>/Quantum; S18=~/qcc/ws/invsqrt_fc/s18; P=~/qcc/.venv/bin/python
NAME=$1; PROG=$(readlink -f $2); SRC=$S18/src_$NAME; CFG=$S18/cfg/cfg_$NAME.json
rm -rf $SRC && cp -r $S18/srcR $SRC && rm -rf $SRC/__pycache__
python3 - "$PROG" "$SRC/fcsup/sc8_q18.prog.json" "$S18/cfg/cfg_R.json" "$CFG" "$SRC" <<'PY'
import json, sys
prog, dst, cbase, cout, src = sys.argv[1:6]
P = json.load(open(prog))
assert "stages" in P and P["inputs"] == ["q17", "w20", "w21", "w22", "w23", "w24", "w25"], P.get("inputs")
P.setdefault("internal", 0); P["target_name"] = "q18"
P["note"] = "FC18: q18 at %d ANDs (exhaustive gate-0 search), replaces the 6-AND exit" % sum(len(s["gates"]) for s in P["stages"])
json.dump(P, open(dst, "w"), indent=1)
c = json.load(open(cbase)); c["fc2_src"] = src; c["h1c"]["cubprog"] = src + "/cub_rot2e6e_acc.json"
json.dump(c, open(cout, "w"), indent=1)
print("q18 program:", sum(len(s["gates"]) for s in P["stages"]), "ANDs; cfg", cout)
PY
cd $Q
bash tools/sh/invsqrt_fc_s18emit.sh $NAME $SRC $CFG 10 | grep -E "emit /|EMIT|Error" | cut -c1-160
[ -f $S18/files/$NAME/invsqrt.kmx ] || { echo "EMIT FAILED"; exit 1; }
(cd tools && nice -n 10 taskset -c 10 $P r63_inv_ghost_chk.py check $S18/files/R/invsqrt.kmx $S18/files/$NAME/invsqrt.kmx 0 2>&1 \
  | tail -6 | cut -c1-200) | tee $S18/logs/$NAME.chk.log
grep -q "CHECK PASS" $S18/logs/$NAME.chk.log || { echo "ALL-INPUT CHECK FAILED: stop"; exit 1; }
bash tools/sh/invsqrt_fc_s18hold.sh $NAME $SRC $CFG 10,11 relabel > $S18/logs/hold_$NAME.log 2>&1
tail -n 3 $S18/logs/hold_$NAME.log | cut -c1-200
bash tools/sh/invsqrt_fc_s18pack.sh $NAME ~/qcc/ws/invsqrt_fc/best/fc18_R25c \
  "fc18_R25c with q18's END exit re-made by $PROG (FC18's exhaustive gate-0 search on I39's staged SAT); exact on ALL inputs against R25c's emit" \
  "SRC = srcR (I39's srcN + tools/invsqrt_fc18_patch_r25c.py) with fcsup/sc8_q18.prog.json = the program; CFG = $CFG; bash tools/sh/invsqrt_fc_s18hold.sh $NAME SRC CFG 10,11 relabel"
