#!/usr/bin/env bash
# f9_i1_stack.sh (seat I1, lead quantum-f9, LOCAL WSL): ONE COMMAND from a found program to a file STACKED ON b196lend.
#   bash tools/f9_i1_stack.sh KIND NAME PROG.json [CORES=10]
#     KIND d20  : digit 20 lent in MAIN's undo  (cfg "_w21ret" = PROG; source = copy of FC's s18/srcP13M)   = s18digitm's recipe
#     KIND d21  : digit 21 re-made by PROG      (cfg "_w22ret" = PROG; source = copy of FC's s18/srcP13M)   = s18digitm's recipe
#     KIND q18  : q18's END exit re-made by PROG (source = copy of FC's s18/srcM with fcsup/sc8_q18.prog.json) = s18q18m's recipe
#     PROG: for d20 / d21 the decode of I40's tool (F1 flat) or a staged file, model convention; for q18 GHOST's staged file.
#   env: S5LEND=1 (default; 0 = the same pipeline WITHOUT b196lend's change), VERIFY=1 (official verifier, both salts, at the
#        end: tools/sh/box_verify.sh, box 4 lane; since 11:45Z Nov 20 box 7, then this PC, when box 4 does not answer),
#        BASE=<text the final file is proved against> (default b196lend's text).
# Pipeline = S5's tools/s5_inv_h1end_build.sh (itself a copy of FC's s19hold) with SRC / CFG / output made parameters:
#   copy of the source -> S5's generator patch (tools/s5_inv_h1end_genpatch.py: H1's word un-add lends q11) -> cfg (FC's cfg
#   + the stub key + cfg_b312g's "r91_splitdebt") -> emit -> ALL-input check of the emit against fc18_M's emit -> width step
#   (relabel to 39) -> all-input check -> P1 pass -> post-chain -> R91 split-debt pass (the b312 carry) -> all-input check
#   -> EXACT PROOF of the final file against BASE on ALL 3,133,452 inputs (tools/r63_inv_ghost_chk.py check BASE NEW 0:
#   x restored, every other wire zero, phase equal to the base's, outcomes 0 / 1 / two random seeds).
# Writes ONLY under ~/qcc/ws/f9_i1/ (src/NAME, cfg/, prog/, files/NAME*, check/NAME, best/NAME) and /tmp/f9i1nb_NAME.
set -u
Q=<HOME>/Quantum; P=~/qcc/.venv/bin/python
S18=~/qcc/ws/invsqrt_fc/s18; R59=~/qcc/ws/invsqrt_fc/best/r59b0
W=~/qcc/ws/f9_i1; F=$W/files
KIND=$1; NAME=$2; PROG=$(readlink -f "$3"); CORES=${4:-10}
BASE=$(readlink -f ${BASE:-~/qcc/ws/s5/inv_h1end/best/b196lend/invsqrt.kmx})
export S5LEND=${S5LEND:-1}
case "$KIND" in
  d20|d21) PARENT=$S18/srcP13M; CB=$S18/cfg/cfg_P13M.json;;
  q18) PARENT=$S18/srcM; CB=$S18/cfg/cfg_M.json;;
  *) echo "KIND must be d20, d21 or q18"; exit 2;;
esac
[ -f "$PROG" ] || { echo "no program file $3"; exit 2; }
mkdir -p $W/src $W/cfg $W/prog $W/check/$NAME $F
T0=$(date -u +%s); el() { echo "$(( $(date -u +%s) - T0 )) s"; }
SRC=$W/src/$NAME; CFG=$W/cfg/cfg_$NAME.json; PG=$W/prog/$NAME.prog.json; C=$W/check/$NAME
echo "== f9_i1_stack $KIND $NAME S5LEND=$S5LEND cores $CORES start $(date -u +%H:%M:%SZ); program $PROG; base $BASE"
cp "$PROG" $PG
rm -rf $SRC; cp -r $PARENT $SRC; rm -rf $SRC/__pycache__
$P -I $Q/tools/s5_inv_h1end_genpatch.py $SRC || { echo "GENERATOR PATCH FAILED"; exit 1; }
mkcfg() {   # $1 = flip (0 / 1): digit 20 only, the output constant flipped (FC's s19digfin retry, CH1 08:36Z)
python3 - "$KIND" "$PG" "$CB" "$CFG" "$SRC" "$S18/cfg/cfg_b312g.json" "$1" <<'PY'
import json, sys
kind, pg, cb, cout, src, cb312, flip = sys.argv[1:8]
P = json.load(open(pg)); c = json.load(open(cb))
if kind in ("d20", "d21"):
    if "stages" not in P and "gates" in P and "outputs" in P:      # s18digitm's conversion: F1's flat output -> one stage
        if flip == "1":
            out = P["outputs"][0]
            P["outputs"][0] = [j for j in out if j != -1] if -1 in out else [-1] + out
        n, k = len(P["inputs"]), len(P["gates"])
        out = list(P["outputs"][0])
        if out.count(n + k - 1) == 1:
            st = {"gates": P["gates"], "internal": k - 1, "out": [o for o in out if o != n + k - 1]}
        else:
            st = {"gates": P["gates"], "internal": k, "out": out}
        P = {"inputs": P["inputs"], "internal": 0, "stages": [st], "outputs": [[]], "target_name": "d%s" % kind[1:]}
        json.dump(P, open(pg + ".staged.json", "w"), indent=1)
        pg = pg + ".staged.json"
        print("converted F1's flat program to one stage: %d gates, %d temporaries, flip %s" % (k, st["internal"], flip))
    c["_w21ret" if kind == "d20" else "_w22ret"] = pg
    print("cfg key", "_w21ret" if kind == "d20" else "_w22ret", "=", pg)
else:
    assert "stages" in P, "q18 needs GHOST's staged format (tools/invsqrt_fc_s19stage.py)"
    if P["inputs"] != ["q17", "w20", "w21", "w22", "w23", "w24", "w25"]:
        print("NOTE: q18 inputs in another order than the search's:", P["inputs"])
    P.setdefault("internal", 0); P["target_name"] = "q18"
    json.dump(P, open(src + "/fcsup/sc8_q18.prog.json", "w"), indent=1)
    c["fc2_src"] = src; c["h1c"]["cubprog"] = src + "/cub_rot2e6e_acc.json"
    print("q18 program:", sum(len(s["gates"]) for s in P["stages"]), "ANDs")
c["r91_splitdebt"] = json.load(open(cb312))["r91_splitdebt"]
json.dump(c, open(cout, "w"), indent=1)
PY
}
O=$F/$NAME; WD=${O}_w
export R41_XT=MAIN C1ADD=1 C1PRED=ctop C1UNDO=1 C1MIX=1 C1CCZ=1 C1MIXBEST=1 "FCGEN=${FCGEN:-7,3,0;9,3,1}" R41_XC=DIV1 FXACCCAP=30
unset R56_LEV
export NUMBA_CACHE_DIR=/tmp/f9i1nb_$NAME R55_NET=$R59/router.json I39_CAPS=${I39_CAPS:-} PYTHONDONTWRITEBYTECODE=1
emit() {
  mkdir -p $O $WD; rm -f $O/invsqrt.kmb $O/invsqrt.kmx $WD/invsqrt.kmx
  cp $CFG $O/cfg.json; cp $R59/angles.json $O/angles.json
  echo "== emit $(date -u +%H:%M:%SZ) ($(el)) S5LEND=$S5LEND"
  (cd /tmp && nice -n 10 taskset -c $CORES $P -I $Q/tools/r64_inv_stack_xcrun.py emit $SRC $O $O/angles.json $O/cfg.json > $C/emit.log 2>&1)
  grep -E "S5LEND|R55 net|emit |per-tag|r63:|Traceback|Error|assert" $C/emit.log | cut -c1-260
  [ -f $O/invsqrt.kmx ] || { echo "EMIT FAILED (no .kmx); log $C/emit.log"; return 1; }
  echo "== all-input check of the emit against fc18_M's emit $(date -u +%H:%M:%SZ) ($(el))"
  (cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $S18/files/M/invsqrt.kmx $O/invsqrt.kmx 0 2>&1 | tail -6 | cut -c1-220) | tee $C/emit_chk.log
  grep -q "CHECK PASS" $C/emit_chk.log
}
mkcfg 0 || exit 1
if ! emit; then
  [ "$KIND" = d20 ] && [ ! -f $PG.staged.json ] && { echo "EMIT / CHECK FAILED: stop"; exit 1; }
  [ "$KIND" = d20 ] || { echo "EMIT / ALL-INPUT CHECK FAILED: stop"; exit 1; }
  echo "digit 20: retrying with the output constant flipped (complement convention)"
  mkcfg 1 || exit 1
  emit || { echo "EMIT / ALL-INPUT CHECK FAILED with the flipped constant too: stop"; exit 1; }
fi
if [ "$S5LEND" = 1 ]; then grep -q "S5LEND: lent" $C/emit.log || { echo "S5LEND DID NOT ACT in the emit: stop"; exit 1; }; fi
echo "== width step (relabel) $(date -u +%H:%M:%SZ) ($(el))"
mkdir -p ${O}_e ${O}_x
(cd $Q/tools && $P r63_inv_ghost_early.py $O/invsqrt.kmx ${O}_e/invsqrt.kmx q18,q19,q20,q21 2>&1 | tail -1)
(cd $Q/tools && $P inv_i39_xpair.py ${O}_e/invsqrt.kmx ${O}_x/invsqrt.kmx 2>&1 | tail -1 | cut -c1-80)
(cd $Q/tools && nice -n 10 taskset -c $CORES $P inv_i39_relabel.py ${O}_x/invsqrt.kmx $WD/invsqrt.kmx 2>&1 | tail -1)
[ -f $WD/invsqrt.kmx ] || { echo "WIDTH STEP FAILED"; exit 1; }
echo "== all-input check (emit vs relabelled) $(date -u +%H:%M:%SZ) ($(el))"
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $O/invsqrt.kmx $WD/invsqrt.kmx 0 2>&1 | tail -2 | cut -c1-200) | tee $C/width_chk.log
grep -q "CHECK PASS" $C/width_chk.log || { echo "WIDTH ALL-INPUT CHECK FAILED: stop"; exit 1; }
(cd $Q/tools && $P -c "import invsqrt_fw_rw as R; R.to_kmb('$WD/invsqrt.kmx', '$WD/invsqrt.kmb')" 2>&1 | tail -2)
cp $O/cfg.json $O/angles.json $WD/
echo "== P1 pass + post-chain $(date -u +%H:%M:%SZ) ($(el))"
cd $Q && bash tools/sh/invsqrt_pp_affpass.sh $WD/invsqrt.kmb ${WD}.aff $CORES 2>&1 | grep -E "expected Toffoli-class|sites [0-9]+, rewritten|CHECK" | tr '\n' ' '; echo
A=${WD}.aff/invsqrt.kmb; [ -f $A ] || A=$WD/invsqrt.kmb
bash tools/sh/invsqrt_fc_s7pubchain.sh ${WD}.pc $CORES $A 2>&1 | tail -2
FIN=${WD}.pc/$(basename $(dirname $A))/final.kmb
[ -f $FIN ] || { echo "POST-CHAIN FAILED (no final.kmb)"; exit 1; }
echo "post-chain final: $FIN  sha $(sha256sum $FIN | cut -c1-16)  $(date -u +%H:%M:%SZ) ($(el))"
echo "== R91 split-debt pass (the b312 carry) $(date -u +%H:%M:%SZ) ($(el))"
rm -rf ${WD}.r91
(cd $Q/tools && nice -n 10 taskset -c $CORES $P invsqrt_fc_s19r91pass.py $FIN $O/cfg.json ${WD}.r91 2>&1 | tail -6 | cut -c1-240) | tee $C/r91.log
CARRY=1
if [ -f ${WD}.r91/invsqrt.kmb ]; then
  (cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check ${WD}.r91/x.kmx ${WD}.r91/invsqrt.kmx 0 2>&1 | tail -2 | cut -c1-200) | tee ${WD}.r91/chk.log
  grep -q "CHECK PASS" ${WD}.r91/chk.log || CARRY=0
else CARRY=0; fi
if [ "$CARRY" = 1 ]; then LAST=${WD}.r91/invsqrt.kmb; else LAST=$FIN; echo "!! b312 CARRY NOT INSTALLED (block not found or check failed): the final file is the post-chain file (+0.5 a shot); re-derive with tools/r91_inv_step_build_exact.py"; fi
B=$W/best/$NAME; mkdir -p $B; cp $LAST $B/invsqrt.kmb; cp $O/cfg.json $O/angles.json $B/; cp $PG $B/program.json
(cd $Q/tools && $P -I -c "
import sys; sys.path.insert(0, '<HOME>/qcc/data/kickmix'); import kickmix as km
c = km.Circuit.from_file('$B/invsqrt.kmb', format='kmb'); c.to_file('$B/invsqrt.kmx', format='kmx'); print('text written, ops', len(c), 'qubits', c.num_qubits)")
sha256sum $B/invsqrt.kmb | cut -c1-64 > $B/sha256.txt
echo "final: $B/invsqrt.kmb  sha $(cut -c1-16 $B/sha256.txt)  carry $CARRY  $(date -u +%H:%M:%SZ) ($(el))"
NQ=$(grep -oE '\bq[0-9]+\b' $B/invsqrt.kmx | sort -u | wc -l)   # 05:19Z Nov 20: the 39-wire gate
[ "$NQ" = 39 ] && echo "wires: 39" || echo "!! THE FILE USES $NQ WIRES, NOT 39: its cost is $NQ x (Toffolis a shot + 237.42); above 39 wires it is no candidate (b196lend is 39 x 1,616.80)"
echo "== EXACT PROOF against the base, ALL inputs $(date -u +%H:%M:%SZ) ($(el))"
(cd $Q/tools && nice -n 10 taskset -c $CORES $P r63_inv_ghost_chk.py check $BASE $B/invsqrt.kmx 0 2>&1 | tail -8 | cut -c1-260) | tee $C/proof.log
{ echo "seat I1 (session quantum-f9), $(date -u +%Y-%m-%dT%H:%MZ): bash tools/f9_i1_stack.sh $KIND $NAME $PROG  (S5LEND=$S5LEND, b312 carry $CARRY)"
  echo "  = FC's hold pipeline on a copy of $PARENT + S5's generator patch (q11 lent in H1's word un-add) + the program $PG"
  echo "  exact proof against $BASE: $(grep -c 'CHECK PASS' $C/proof.log) CHECK PASS line(s), log $C/proof.log"
  echo "  sha256 $(cat $B/sha256.txt)"; } > $B/SOURCE.txt
grep -q "CHECK PASS" $C/proof.log || { echo "EXACT PROOF FAILED: NOT a file. $(date -u +%H:%M:%SZ) ($(el))"; exit 1; }
echo "== PROVED $(date -u +%H:%M:%SZ): wall-clock $(el) from the program to the proved file $B/invsqrt.kmb"
if [ "${VERIFY:-0}" = 1 ]; then
  # 11:45Z Nov 20 (sixth holder): AWS box 4 (the lane) stopped answering at about 10:55Z. Probe box 4, then box 7 (our own
  # AWS box, verifier kit present, checked 11:39Z), with 20 s ssh connects; if neither answers, verify ON THIS PC with
  # tools/f9_i1_h6vloc.sh (the same tools/qcc_verify.py, 2 workers, both salts). VBOX=n in the env forces a box.
  VB=${VBOX:-}
  if [ -z "$VB" ]; then
    for b in 4 7; do BOX=$b timeout 40 bash $Q/tools/sh/box_sync.sh ssh 'true' < /dev/null > /dev/null 2>&1 && { VB=$b; break; }; done
  fi
  if [ -n "$VB" ]; then
    echo "== official verifier, both salts (box $VB) $(date -u +%H:%M:%SZ)"
    (cd $Q && VBOX=$VB bash tools/sh/box_verify.sh $B invsqrt f9_i1_$NAME 2>&1 | cut -c1-200)
  else
    echo "== official verifier, both salts, ON THIS PC (boxes 4 and 7 do not answer) $(date -u +%H:%M:%SZ)"
    (cd $Q && bash tools/f9_i1_h6vloc.sh $B 2>&1 | cut -c1-200)
  fi
  echo "== verified $(date -u +%H:%M:%SZ): wall-clock $(el)"
fi
