"""chem_t10_prepair (T10 seat 3, Nov 18): THE END REPAIR OF A RELATION TABLE IN CARRIER FORM (all_c.npz after C2's gauge).

A relation broken on the stored words (a late pole changed an eps after its words froze) is repaired on the STORED words:
per row, the cheapest set of 1 to 3 +-2 moves on bit-1 relation columns (bit 0 kept, carriers bit6^bit5 kept) after which
EVERY relation of the row holds (all loaders: a move of a batch-k word toggles loaders k and k+1), priced EXACTLY at fixed
weights (two-body rows: Delta error_2^2 by chem_t10_k1post.Pricer; one-body rows: Delta error_1^2). The weights are then
re-solved by the caller's refit. Rows in address order 0..1218; targets from the relation file (R65's G/ targets).
    python chem_t10_prepair.py all_c.npz OUT.npz --rel @rel_final.txt [--brot 7]
"""
import argparse
import itertools
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_t10_k1post import Pricer  # noqa: E402
from r65_chem_all15_relg import parse_g, ttable_g  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("out")
    ap.add_argument("--rel", required=True)
    ap.add_argument("--brot", type=int, default=7)
    a = ap.parse_args()
    z = np.load(a.table)
    d = {k: np.array(z[k]) for k in z.files}
    w0 = np.mod(np.asarray(d["angle_words"], np.int64), 1 << a.brot)
    edges = np.asarray(d["tree_edges"], np.int64)
    wt = np.asarray(d["W"], np.float64)
    shape = wt.shape
    Woff = np.asarray(d["H2_I"], np.float64) - wt.sum(axis=1)
    n = w0.shape[0]
    N1 = edges.shape[0] + 1
    addr = np.arange(n)
    rels = parse_g(a.rel)
    T = np.stack([ttable_g(t, addr) for (_k, _o, t) in rels], axis=1)          # (rows, relations)

    def parities(w):
        out = np.zeros((len(w), len(rels)), np.int64)
        for i, (k, o, _t) in enumerate(rels):
            for (j, b) in o:
                out[:, i] ^= ((w[:, 5 * k + j] ^ w[:, 5 * k - 5 + j]) >> b) & 1
        return out

    bad = parities(w0) != T
    rows = np.flatnonzero(bad.any(axis=1))
    print(f"{a.table}: {int(bad.sum())} (row, relation) pairs broken on {len(rows)} rows "
          f"(one-body {int((rows < N1).sum())}, two-body {int((rows >= N1).sum())})", flush=True)
    if not len(rows):
        np.savez(a.out, **d)
        print("nothing to repair")
        return
    b1cols = sorted({5 * k + j for (k, o, _t) in rels for (j, b) in o if b == 1} |
                    {5 * k - 5 + j for (k, o, _t) in rels for (j, b) in o if b == 1})
    car0 = ((w0 >> 6) ^ (w0 >> 5)) & 1
    pr = Pricer(w0, edges, shape, a.brot, wt, Woff)
    w1 = w0.copy()
    tot2 = tot1 = 0.0
    left = 0
    for r in rows:
        moves = []                                   # (col, new word)
        for c in b1cols:
            for dlt in (2, -2):
                nw = (w1[r, c] + dlt) % (1 << a.brot)
                if (((nw >> 6) ^ (nw >> 5)) & 1) == car0[r, c]:
                    moves.append((c, nw))
        # price every single move exactly
        if r >= N1:
            rows2 = np.repeat(w1[r][None, :], len(moves), axis=0)
            for i, (c, nw) in enumerate(moves):
                rows2[i, c] = nw
            price = pr.dE2(np.full(len(moves), r - N1, np.int64), pr.vec(rows2)) * 1e6
        else:
            price = np.array([pr.dE1(r, pr.vec(np.where(np.arange(w1.shape[1]) == c, nw, w1[r])[None, :])[0]) * 1e6
                              for (c, nw) in moves])
        order = np.argsort(price)
        best = None
        for size in (1, 2, 3):
            pool = order[: (len(order) if size < 3 else min(len(order), 24))]
            for comb in itertools.combinations(pool, size):
                if len({moves[i][0] for i in comb}) < size:
                    continue
                v = w1[r].copy()
                for i in comb:
                    v[moves[i][0]] = moves[i][1]
                if np.array_equal(parities(v[None, :])[0], T[r]):
                    cst = float(sum(price[i] for i in comb))
                    if best is None or cst < best[0]:
                        best = (cst, comb)
            if best is not None:
                break
        if best is None:
            left += 1
            continue
        for i in best[1]:
            w1[r, moves[i][0]] = moves[i][1]
        if r >= N1:
            tot2 += best[0]
        else:
            tot1 += best[0]
    assert np.array_equal(car0, ((w1 >> 6) ^ (w1 >> 5)) & 1) and np.array_equal(w0 & 1, w1 & 1)
    nb = int((parities(w1) != T).sum())
    print(f"REPAIRED: rows left {left}; relations broken now {nb}; first-order price two-body {tot2:.3f} mHa^2 of error_2^2, "
          f"one-body {tot1:.3f} mHa^2 of error_1^2", flush=True)
    d["angle_words"] = w1
    np.savez(a.out, **d)
    print("saved", a.out)


if __name__ == "__main__":
    main()
