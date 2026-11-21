"""R58-CHEM-KEPT: a padding of the table with a block set V for MORE loaders (builder C2's vset and pad_table,
imported read-only), so that the kept-product loader can use the 68-AND cache on them.

  CHEM_BROT=7 CHEM_HB=77 python -I r58_chem_kept_vpad.py TABLE.npz --ks 1,2,3,5,7,11 --old in/pad.npz --out pad2.npz
Prints the V of every loader asked for, checks the padded table (check_padded) and that the real rows are the
table's own; --old: the builders' padding, to report whether its V's are kept.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_c4 as c4  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--ks", default="1,2,3,5,7,11")
    ap.add_argument("--old", default=None)
    ap.add_argument("--proj", default=None, help="keep the V's of this proj.json where it has them")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J = 5
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    A = words.shape[0]
    words = words[:min(A, 1219)] if A > 1219 else words
    print(f"table rows {A} -> real rows {words.shape[0]}")
    keep = json.load(open(a.proj)).get("V", {}) if a.proj else {}
    Vs = {}
    for k in [int(x) for x in a.ks.split(",")]:
        V = int(keep[str(k)], 16) if str(k) in keep else c4.vset(words, k)
        Vs[k] = V
        print(f"loader {k}: V {hex(V) if V is not None else None}" + (f" (size {bin(V).count('1')})" if V else ""),
              flush=True)
    wp = c4.pad_table(words, Vs)
    bad = c4.check_padded(wp, Vs)
    print("check_padded (bits whose V-sum is not constant, per loader):", bad)
    assert np.array_equal(wp[:words.shape[0]], words)
    if a.old:
        old = np.asarray(np.load(a.old)["angle_words"], np.int64)
        print("padding rows equal to the builders' padding:", bool(np.array_equal(old, wp)))
    np.savez(a.out, angle_words=wp, tree_edges=z["tree_edges"] if "tree_edges" in z.files else np.zeros((0, 2)))
    json.dump({"V": {str(k): hex(v) for k, v in Vs.items() if v is not None}}, open(a.out + ".json", "w"))
    print("written", a.out)


if __name__ == "__main__":
    main()
