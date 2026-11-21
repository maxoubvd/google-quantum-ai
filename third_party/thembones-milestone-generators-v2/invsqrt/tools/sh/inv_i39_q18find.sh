#!/usr/bin/env bash
# I39 (third holder), LOCAL WSL: a found q18 END-exit program (GHOST's staged format, inputs q17, w20..w25, at most 3
# temporaries) -> a verified file in one command, the way Egq6's 6-AND q18 went in.
#   bash tools/sh/inv_i39_q18find.sh NAME PROG.json
# Steps: copy of Egq6n25's source (hold/srcN) with fcsup/sc8_q18.prog.json replaced -> emit + ALL-input check against
# Egq6n25's emit (four outcome modes) -> the hold pipeline (relabel to 39, P1, post-chain) -> SOURCE.txt -> verifier.
# Safe to run: writes only my folders.
set -u
Q=<HOME>/Quantum; H=~/qcc/ws/inv_i39/hold; F=$H/files17
NAME=$1; PROG=$2; S=$H/src_$NAME; CFG=$H/cfg_$NAME.json
rm -rf $S && cp -r $H/srcN $S && rm -rf $S/__pycache__
cp $S/fcsup/sc8_q18.prog.json $S/fcsup/sc8_q18.prog.json.pre_$NAME
python3 - "$PROG" "$S/fcsup/sc8_q18.prog.json" "$H/cfg_N.json" "$CFG" "$S" <<'PY'
import json, sys
prog, dst, cbase, cout, src = sys.argv[1:6]
P = json.load(open(prog))
assert "stages" in P and P["inputs"] == ["q17", "w20", "w21", "w22", "w23", "w24", "w25"], P.get("inputs")
P.setdefault("internal", 0); P["target_name"] = "q18"
P["note"] = "I39 third holder: q18 at %d ANDs (gate-0 portfolio SAT), replaces the 6-AND exit" % sum(len(s["gates"]) for s in P["stages"])
json.dump(P, open(dst, "w"), indent=1)
c = json.load(open(cbase)); c["fc2_src"] = src; c["h1c"]["cubprog"] = src + "/cub_rot2e6e_acc.json"
json.dump(c, open(cout, "w"), indent=1)
print("q18 program:", sum(len(s["gates"]) for s in P["stages"]), "ANDs; cfg", cout)
PY
bash $Q/tools/sh/inv_i39_emitchk.sh $NAME $S $CFG 4,5 $H/files17/N25/invsqrt.kmx 2>&1 | tee $H/files18/$NAME.chk.log | tail -n 8
if ! grep -q "CHECK PASS" $H/files18/$NAME.chk.log || grep -q -E "CHECK FAIL|EMIT FAILED|differs from the base on [1-9]|restored False|dirty wires \[[0-9]|Traceback" $H/files18/$NAME.chk.log; then
  echo "ALL-INPUT CHECK FAILED: stop"; exit 1
fi
sed -e "s/^invsqrt: 39-QUBIT FILE Egq6n25 = /invsqrt: 39-QUBIT FILE $NAME = Egq6n25 with q18's END exit re-made at fewer ANDs (I39 third holder, gate-0 portfolio SAT); Egq6n25 = /" \
    ~/qcc/ws/inv_i39/best/hold39egq6n25/SOURCE.txt > $H/SOURCE_$NAME.txt
echo "- $NAME: fcsup/sc8_q18.prog.json replaced by $PROG (tools/inv_i39_ssat.py with g0in; checked on ALL inputs against Egq6n25's emit)." >> $H/SOURCE_$NAME.txt
bash $Q/tools/sh/inv_i39_s17hold.sh $NAME $S $CFG 4,5 relabel > $F/$NAME.log 2>&1
tail -n 4 $F/$NAME.log
bash $Q/tools/sh/inv_i39_autoverify.sh $NAME hold39$(echo $NAME | tr 'A-Z' 'a-z') $H/SOURCE_$NAME.txt $S/fcsup/sc8_q18.prog.json 2>&1 | tail -n 8
