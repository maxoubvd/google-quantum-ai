"""invsqrt_gen: after a publication, bring the short facts of the state block of research/ws_invsqrt_gen.md in
line with the file at the top of best/ (own tool; the long text of the block is rewritten by hand when the
circuit changes in kind).

usage: invsqrt_gen_newbest.py "<what the build adds, for the table of published builds>" [key=old=>new ...]
Reads best/STATE.txt, best/ENTRY and the two verifier files and
  - rewrites the paragraph "**Best verified build: ..." up to its sha256 and place (tag, scores, Q, Toffolis,
    rotations, worst input, sha256, folder, time) and the entry line that follows "Entry (one line";
  - appends one row to the table of published builds (if the tag has none yet);
  - sets the "last edit" time of the header.
Optional arguments of the form  old text=>new text  are plain replacements inside the state block (for the
stage table's numbers). Everything else is left as it is.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time

LOG = "<HOME>/Quantum/research/ws_invsqrt_gen.md"
W = os.path.expanduser("~/qcc/ws/invsqrt_gen")


def main():
    what = sys.argv[1]
    ent = open(os.path.join(W, "best", "ENTRY")).read().split()
    tag = ent[0]
    v = [json.load(open(os.path.join(W, "best", n)))["invsqrt"] for n in ("verify_salt_empty.json", "verify_salt_alt1.json")]
    assert all(x["status"] == "PASS" for x in v) and v[0]["circuit_sha256"] == v[1]["circuit_sha256"]
    sha = v[0]["circuit_sha256"]
    assert open(os.path.join(W, "best", "invsqrt.kmb.sha256")).read().split()[0] == sha
    st = open(os.path.join(W, "best", "STATE.txt")).read()
    when = re.search(r"published (\S+) (\d\d:\d\d)", st)
    s = open(LOG, encoding="utf-8").read()
    i0 = s.index("**Best verified build:")
    i1 = s.index("`~/qcc/ws/invsqrt_gen/best/STATE.txt`", i0)
    pub = 93863.0
    para = (f"**Best verified build: {tag}, {v[0]['score']:,.2f} / {v[1]['score']:,.2f} (salt \"\" / alt-salt-1, both PASS at "
            f"{v[0]['verified_shots']:,} shots on\nthe verification box; Q {v[0]['num_qubits']}, "
            f"{v[0]['toffoli_count']:,.2f} Toffolis, {v[0]['rotation_count']:.0f} rotations + {v[0]['t_count']:.0f} T, worst "
            f"{v[0]['worst_error_over_tolerance']:.5f} at x = {v[0]['worst_x']:,}), sha256\n{sha}, in "
            f"`~/qcc/ws/invsqrt_gen/best/{tag}/` (circuit,\nsha256, both verifier files, build.json, ENTRY, the table "
            f"file, and in `lut/` the lookup programs and specs\nit used) and at the top of `best/` since {when.group(2)}Z "
            f"({when.group(1)}). The public file is at 93,863 with 62 qubits:\nour own circuit is "
            f"{pub - max(x['score'] for x in v):,.0f} points ({100 * (pub - max(x['score'] for x in v)) / pub:.1f}%) below it, on "
            f"{62 - v[0]['num_qubits']} wires less [M].\n")
    s = s[:i0] + para + s[i1:]
    j0 = s.index("`invsqrt_gen_build.sh ", s.index("Entry (one line"))
    j1 = s.index("`", j0 + 1)
    words, lines, cur = ["invsqrt_gen_build.sh"] + ent, [], ""
    for w in words:
        if len(cur) + len(w) + 1 > 108:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    lines.append(cur)
    s = s[:j0] + "`" + "\n".join(lines) + s[j1:]
    row = (f"| {tag} | {what} | {v[0]['num_qubits']} | {v[0]['toffoli_count']:,.2f} | {v[0]['rotation_count']:.0f} | "
           f"{v[0]['score']:,.2f} / {v[1]['score']:,.2f} |")
    if f"\n| {tag} |" not in s:
        k0 = s.index("| tag | what it adds |")
        k1 = s.index("\n\n", k0)
        s = s[:k1] + "\n" + row + s[k1:]
    for a in sys.argv[2:]:
        old, new = a.split("=>", 1)
        assert old in s, f"not found: {old}"
        s = s.replace(old, new, 1)
    s = re.sub(r"^## Goals and state \(kept current; last edit [0-9-]+ [0-9:]+Z\)",
               "## Goals and state (kept current; last edit " + time.strftime("%Y-%m-%d %H:%MZ", time.gmtime()) + ")",
               s, count=1, flags=re.M)
    open(LOG, "w", encoding="utf-8", newline="\n").write(s)
    print(para.split(" (salt")[0] + " ... sha " + sha[:16] + "; row: " + row[:90])


if __name__ == "__main__":
    main()
