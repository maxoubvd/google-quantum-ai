#!/usr/bin/env python3
"""r89_chem_cascade_cells.py (seat R89-CHEM-CASCADE): usable dead cells of a loader's batch(es) for EVERY generic output
(C4's rules, tools/chem_c4_dead.usable_cells, read-only), with the output -> wire map (wire basis).
  python -I r89_chem_cascade_cells.py <table.npz> --k 14 [--batches own|prev|both]
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_c4_wb  # noqa: E402
import chem_c4_dead as cd  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, default=14)
    ap.add_argument("--batches", default="own")
    ap.add_argument("--json", default="", help="write {output: usable cells} here")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    import chem_ld2_e19grp as eg
    z = np.load(a.table)
    w0 = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    V, words = eg.block_set(w0, a.k)
    an = core.analyze(words, a.k, 4)
    wires = chem_c4_wb.wire_basis(an)
    gen = [o for o in range(an.X0.shape[1]) if an.out_level[o] > 2]
    db = {"own": None, "prev": [a.k - 1], "both": [a.k, a.k - 1]}[a.batches]
    cells = cd.usable_cells(w0, edges, a.k, an, wires, gen, V, db)
    cnt = {o: 0 for o in gen}
    for c in cells:
        cnt[c[0]] += 1
    print(f"k{a.k} ({a.batches}): {len(cells)} usable cells over {len(gen)} generic outputs")
    for o in sorted(gen, key=lambda o: -cnt[o]):
        print(f"  output {o:2d} wire {an.out_wires[wires[o]]}: {cnt[o]}")
    if a.json:
        import json
        json.dump({str(o): cnt[o] for o in gen}, open(a.json, "w"))


if __name__ == "__main__":
    main()
