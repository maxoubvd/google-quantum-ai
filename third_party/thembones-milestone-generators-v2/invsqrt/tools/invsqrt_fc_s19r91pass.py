#!/usr/bin/env python3
"""invsqrt_fc_s19r91pass.py (FC nineteenth holder): the generator's R91 split-debt pass (cfg key "r91_splitdebt").

After FC's post-chain, each entry of cfg["r91_splitdebt"] pays one coin's phase debt by Cliffords in two pieces (R91 seat
R91-INV-C2-SPLITDEBT's method and patch tool) instead of the coin's own repair block. An entry names the debt by TEXT, coin id
free, so it survives re-numbering of coins and ops by other changes:
  {"name": "b312", "hmr": "HMR q38 b{K}", "block": ["CCX q12 q13 q29 if b{K}", "CZ q12 q14 if b{K}", "PUSH_CONDITION if b{K}"],
   "t_gate": "CCX q11 q29 q33", "g_terms": [[10], [10, 13], ...], "h_terms": [[11], [31]]}
The pass finds the HMR whose next 200 ops hold every block line (exactly one hit), the first t_gate after the block, writes
C2's spec (b, h_op, t_op, g_terms, h_terms) and runs C2's tools/r91_inv_step_build_patch.py (read-only use): it deletes the
coin's conditioned ops and PUSH/POP span and inserts the Cliffords under the coin. No key -> the input is copied unchanged.
  python invsqrt_fc_s19r91pass.py FINAL.kmb CFG.json OUTDIR      (writes OUTDIR/x.kmx, spec.json, invsqrt.kmx, invsqrt.kmb)
The caller then runs the R63 all-input check of OUTDIR/invsqrt.kmx against OUTDIR/x.kmx (phase equal on every input)."""
import json
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.expanduser("~/qcc/data/kickmix"))
import kickmix as km  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    fin, cfgp, out = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(out, exist_ok=True)
    entries = json.load(open(cfgp)).get("r91_splitdebt", [])
    c = km.Circuit.from_file(fin, format="kmb")
    text = str(c).strip()
    open(os.path.join(out, "x.kmx"), "w").write(text + "\n")
    if not entries:
        shutil.copy(fin, os.path.join(out, "invsqrt.kmb"))
        shutil.copy(os.path.join(out, "x.kmx"), os.path.join(out, "invsqrt.kmx"))
        print("r91pass: no r91_splitdebt key: copied unchanged")
        return
    lines = text.split("\n")
    nhead = 0
    while nhead < len(lines) and (lines[nhead].startswith("REGISTER") or lines[nhead].startswith("APPEND_TO_REGISTER")):
        nhead += 1
    ops = [s.strip() for s in lines[nhead:]]
    spec = []
    for e in entries:
        pat = re.compile("^" + re.escape(e["hmr"]).replace(re.escape("{K}"), r"(\d+)") + "$")
        hits = []
        for i, s in enumerate(ops):
            m = pat.match(s)
            if not m:
                continue
            k = m.group(1)
            need = [b.replace("{K}", k) for b in e["block"]]
            win = ops[i + 1:i + 200]
            if all(n in win for n in need):
                last = i + 1 + max(win.index(n) for n in need)
                t = next((j for j in range(last, min(len(ops), last + 400)) if ops[j] == e["t_gate"]), None)
                hits.append((int(k), i, t))
        if len(hits) != 1 or hits[0][2] is None:
            raise SystemExit("r91pass: %s: block found %d times (%s): re-derive with tools/r91_inv_step_build_exact.py" % (
                e["name"], len(hits), hits))
        k, h_op, t_op = hits[0]
        spec.append({"b": k, "h_op": h_op, "t_op": t_op, "g_terms": e["g_terms"], "h_terms": e["h_terms"]})
        print("r91pass: %s -> coin b%d, HMR op %d, t op %d" % (e["name"], k, h_op, t_op))
    json.dump(spec, open(os.path.join(out, "spec.json"), "w"), indent=1)
    r = subprocess.run([sys.executable, os.path.join(HERE, "r91_inv_step_build_patch.py"), os.path.join(out, "x.kmx"),
                        os.path.join(out, "spec.json"), out], capture_output=True, text=True)
    print(r.stdout.strip()[-600:])
    if r.returncode != 0 or not os.path.exists(os.path.join(out, "invsqrt.kmb")):
        raise SystemExit("r91pass: C2's patch tool failed: %s" % r.stderr[-400:])


if __name__ == "__main__":
    main()
