#!/usr/bin/env python
"""chem_own_pairs: twins for the strongest vectors of a continuous table (panel design CHEM-1, step 1; CPU).

    chem_own_pairs.py IN.npz OUT.npz [--k 6] [--p 0.3333]

IN.npz: U (R, B', 76), W (R, B', C): a continuous point whose weakest k rows per group have already been removed
and the rest refitted (tools/chem_own_fit.py, phases "dropb:k;dlbfgs:...;adam:..."). OUT.npz: U (R, B' + k, 76),
W (R, B' + k, C) and `pairs` (flat row indices master, twin): in every group the k strongest vectors
(largest |w_rb|^2) get a twin with the same vector; the master keeps the share p of the weights, the twin 1 - p.
The table is the same function as IN (the two rows add up to the old one), so its error is IN's error.
tools/chem_own_round7.py keeps the angles of a pair tied and freezes them jointly ("antithetic split copies":
with p = 1/3 the weighted mean of two neighbouring words moves in thirds of a step). The mechanism and the
numbers (4% of the vectors, p near 1/3) are the panel's (research/panel/CHEM-1_best_chemistry.md); the code is ours.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_own_round as RD  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inp")
    ap.add_argument("out")
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--p", type=float, default=1.0 / 3.0)
    a = ap.parse_args()
    z = np.load(a.inp)
    U, w = np.array(z["U"], np.float64), np.array(z["W"], np.float64)
    U = U / np.linalg.norm(U, axis=2, keepdims=True)
    R, B0, C = w.shape
    B = B0 + a.k
    U2, w2 = np.zeros((R, B, RD.N)), np.zeros((R, B, C))
    U2[:, :B0], w2[:, :B0] = U, w
    pairs = []
    s = np.sum(w ** 2, axis=2)
    for r in range(R):
        top = np.argsort(-s[r], kind="stable")[:a.k]
        for j, b in enumerate(top):
            U2[r, B0 + j] = U[r, b]
            w2[r, B0 + j] = (1.0 - a.p) * w[r, b]
            w2[r, b] = a.p * w[r, b]
            pairs.append((r * B + int(b), r * B + B0 + j))
    tgt = RD.Target()
    e0, e1 = tgt.error2_np(U, w), tgt.error2_np(U2, w2)
    med = np.median(s)
    print(f"{a.inp}: {(R, B0, C)} -> {(R, B, C)}, {len(pairs)} pairs (k {a.k} per group, p {a.p:.4f}); error_2 "
          f"{1e3 * e0:.4f} -> {1e3 * e1:.4f} mHa (the same table); strength of the paired vectors {s[np.arange(R)[:, None], np.argsort(-s, axis=1)[:, :a.k]].min() / med:.1f}"
          f" to {s.max() / med:.1f} x the median, their share of sum |w|^2: {np.sort(s, axis=1)[:, ::-1][:, :a.k].sum() / s.sum():.3f}")
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out, U=U2, W=w2, pairs=np.array(pairs, np.int64), pair_p=a.p)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
