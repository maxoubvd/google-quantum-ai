"""chem_t10_pcheck (T10 seat 3): THE ACCEPTANCE of a relation table after C2's carrier gauge (all_c.npz, rows in address
order 0..1218): carriers exact on all 75 columns (C2's slots 3, 5, 6, 9, 10), bit-0 wires and bit-1 relations held on the
STORED words with the run's fitted targets (rel_used_g.txt; G/ targets), per relation and split one-body / two-body.
    python chem_t10_pcheck.py all_c.npz --rel @rel_used_g.txt [--json OUT.json]
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r65_chem_all15_relg import parse_g, ttable_g  # noqa: E402
from r53_chem_gauge_car import MONO  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--rel", required=True)
    ap.add_argument("--json", default="")
    a = ap.parse_args()
    w = np.mod(np.asarray(np.load(a.table)["angle_words"], np.int64), 128)
    n, T = w.shape
    addr = np.arange(n)
    low = addr % 16
    want = np.stack([((low & MONO[t % 5]) == MONO[t % 5]).astype(np.int64) for t in range(T)], 1)
    car = ((w >> 6) ^ (w >> 5)) & 1
    cols_ok = int((car == want).all(axis=0).sum())
    rels = parse_g(a.rel)
    out = dict(carriers_exact_columns=cols_ok, wrong_carrier_entries=int((car != want).sum()), rels=[])
    tot1 = tot2 = 0
    for (k, o, t) in rels:
        acc = np.zeros(n, np.int64)
        for (j, b) in o:
            acc ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        bad = acc != ttable_g(t, addr)
        b1, b2 = int(bad[:76].sum()), int(bad[76:].sum())
        tot1 += b1
        tot2 += b2
        out["rels"].append(dict(k=k, outs=[list(x) for x in o], one_body=b1, two_body=b2))
        if b1 or b2 or any(b == 1 for (_j, b) in o):
            print(f"  k{k} {o}: broken one-body {b1}, two-body {b2}")
    out.update(broken_one_body=tot1, broken_two_body=tot2)
    print(f"{a.table}: carriers exact on {cols_ok}/75 columns ({out['wrong_carrier_entries']} wrong entries); "
          f"{len(rels)} relations: BROKEN one-body {tot1}, two-body {tot2}")
    if a.json:
        json.dump(out, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
