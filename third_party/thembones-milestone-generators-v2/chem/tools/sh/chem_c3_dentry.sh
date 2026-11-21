#!/usr/bin/env bash
# Builder C3: a candidate chemistry entry = a verified entry with named loaders replaced by kept-product circuits
# (design seat R58-CHEM-KEPT's entry tool, read-only), SOURCE.txt with the chemistry credits block of the lead's
# s061 word for word, the local verifier (salt ""), then the official verifier on both salts on the box lane.
#   bash chem_c3_dentry.sh BASE TAG WORDS.npz "K=file.kmx ..." [NOTE]  (LOCAL WSL; third holder Nov 17 15:00Z: the
# base entry with its angle words replaced by WORDS (a table rewritten at dead words only: R63-CHEM-WIRED entry tool,
# --equal so that rebuilt neighbours of the same count go in)
# Writes ~/qcc/ws/chem_c3/best/TAG/{dfthcsos,SOURCE.txt,verify_salt0.json,box_verify_*.json}.
set -u
BASE=$1; TAG=$2; WORDS=$3; PAIRS=$4; NOTE=${5:-}
Q=<HOME>/Quantum
OUT=$HOME/qcc/ws/chem_c3/best/$TAG
mkdir -p "$OUT"
cd "$Q" || exit 2
~/qcc/.venv/bin/python -I ${ETOOL:-tools/r63_chem_wired_entry.py} "$BASE" "$OUT" --equal "WORDS=$WORDS" $PAIRS || exit 3
# third holder (18:05Z): SUBMISSION=<npz> = a complete submission whose weights were refit with the words FIXED
# (tools/r48_chem_q7c_refit.py): it replaces the entry's submission.npz when its angle words and tree are WORDS's
if [ -n "${SUBMISSION:-}" ]; then
  ~/qcc/.venv/bin/python - "$SUBMISSION" "$WORDS" <<'EOF' || exit 4
import sys
import numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
assert (np.asarray(a["angle_words"])[:1219] == np.asarray(b["angle_words"])[:1219]).all(), "refit changed the words"
assert (np.asarray(a["tree_edges"]) == np.asarray(b["tree_edges"])).all(), "refit changed the tree"
print("refit submission: same words and tree")
EOF
  cp "$SUBMISSION" "$OUT/dfthcsos/submission.npz"
  echo "submission.npz = $SUBMISSION (weights refit, words fixed)" >> "$OUT/SOURCE.txt"
fi
python3 - "$OUT/SOURCE.txt" "$BASE" "$NOTE" <<'EOF'
import sys
src, base, note = sys.argv[1], sys.argv[2], sys.argv[3]
# since Nov 17 12:05Z (second holder): the chemistry credits block of the lead's s065, word for word, from the
# chemistry section of its SOURCE.txt (the line block that starts FIRST after the 'dfthcsos/' line)
S65 = "<HOME>/Quantum/submissions/s065_famC_r59b_q28t54chz_c3g8f6rp/SOURCE.txt"
FIRST = "  CREDIT (added 2026-11-16 12:40Z, research/midprize/CREDITS_and_lineage.md section 3)"
LAST = "loader compiler are our own."
lines = open(S65, encoding="utf-8").read().splitlines()
c0 = next(n for n, s in enumerate(lines) if s.startswith("dfthcsos/"))
i = next(n for n in range(c0, len(lines)) if lines[n].startswith(FIRST))
j = next(n for n in range(i, len(lines)) if lines[n].rstrip().endswith(LAST))
with open(src, "a", encoding="utf-8") as fh:
    fh.write("\nBuilder C3 (research/ws_chem_c3.md): the replaced 30-output loaders are THE KEPT-PRODUCT LOADERS of design\n"
             "seat R58-CHEM-KEPT (research/panel/R58-CHEM-KEPT.md: the point stage with class-mask helpers whose wires\n"
             "are read by later products), on builder C2's 68-AND high cache (8 helper wires), with its re-pricer's\n"
             "relief (research/panel/R58-CHEM-KEPT-reprice.md: the last address block's never-run rows left free at the\n"
             "helper points; helper classes of equal size); stage A by builder C2's engine on a searched basis; tools:\n"
             "chem_c3_bsearch.py, r58_chem_kept_c2a.py, r58_chem_kept_rp_mask.py (--dc --eqw), r58_chem_kept_rp_emit68.py,\n"
             "tools/sh/chem_c3_bsrun.sh. Base entry: " + base + "\n")
    if note:
        fh.write(note + "\n")
    fh.write("\n".join(lines[i:j + 1]) + "\n")
print("SOURCE.txt: credits block appended (%d lines)" % (j - i + 1))
EOF
~/qcc/.venv/bin/python tools/qcc_verify.py "$OUT" --tasks chem --workers 2 --json "$OUT/verify_salt0.json" 2>&1 \
  | grep -E "cost|FAIL|fail|rror" | head -8
[ "${NOBOX:-0}" = "1" ] && exit 0
bash tools/sh/box_verify.sh "$OUT" chem "chem_c3_$TAG" 2>&1 | tail -6
python3 - "$OUT" <<'EOF'
import json, sys, os
d = sys.argv[1]
for s in ("salt0", "alt-salt-1"):
    f = os.path.join(d, f"box_verify_{s}.json")
    if os.path.isfile(f):
        j = json.load(open(f))
        print(s, json.dumps(j)[:400])
EOF
