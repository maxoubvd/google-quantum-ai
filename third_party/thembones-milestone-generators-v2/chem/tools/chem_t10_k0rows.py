"""chem_t10_k0rows (T10 seat 5, Nov 19; kill test K0 (b)): WHERE the relations of a finished relation table break.

Reads a run's END (rounded.npz: words of the placed rows, rows in address order 76 + i; rel_used_g.txt: the relations
with their fitted block constants) and prints, per row with a broken relation on C2's STORED words (R53's carrier form,
tools/chem_t10_relge.stored_words), the broken relations, their loader k and whether they hold a bit-0 output (a wire,
which eps cannot change) or only bit-1 outputs (which a late eps change can flip).

    python chem_t10_k0rows.py RUN_DIR [--json OUT.json]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_t10_relge import stored_words  # noqa: E402  (read-only imports)
from r65_chem_all15_relg import parse_g, ttable_g  # noqa: E402


def per_row(words, addrs, rels):
    """bad[row, relation] = 1 where the relation fails on the given words (rows in address order)."""
    w = np.asarray(words, np.int64)
    bad = np.zeros((len(w), len(rels)), np.int64)
    for i, (k, o, t) in enumerate(rels):
        acc = np.zeros(len(w), np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        bad[:, i] = (acc != ttable_g(t, addrs)).astype(np.int64)
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    z = np.load(os.path.join(a.run, "rounded.npz"))
    words = np.asarray(z["words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    bits = int(z["bits"])
    rels = parse_g("@" + os.path.join(a.run, "rel_used_g.txt"))
    n = words.shape[0]
    addrs = 76 + np.arange(n)
    st = stored_words(np.mod(words, 1 << bits), edges, addrs, bits)
    bad_st = per_row(st, addrs, rels)
    bad_z = per_row(np.mod(words, 1 << bits), addrs, rels)
    wire = np.array([any(b == 0 for (_j, b) in o) for (_k, o, _t) in rels])
    rows = np.flatnonzero(bad_st.sum(1))
    print(f"{a.run}: rows {n}, relations a row {len(rels)}; broken on the stored words {int(bad_st.sum())} in {len(rows)} rows; "
          f"wires among the broken {int(bad_st[:, wire].sum())}; broken on the signed-range words {int(bad_z.sum())}")
    out = []
    for r in rows:
        idx = np.flatnonzero(bad_st[r])
        desc = [f"k{rels[i][0]}{'w' if wire[i] else 'b1'}" for i in idx]
        print(f"  row {r} (address {addrs[r]}, group {r // 127}): {len(idx)} broken: {' '.join(desc)}")
        out.append({"row": int(r), "addr": int(addrs[r]), "broken": [int(i) for i in idx]})
    if a.json:
        json.dump({"rows": out, "total": int(bad_st.sum())}, open(a.json, "w"))


if __name__ == "__main__":
    main()
