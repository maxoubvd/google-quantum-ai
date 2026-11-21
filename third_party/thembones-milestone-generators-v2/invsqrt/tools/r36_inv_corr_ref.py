#!/usr/bin/env python3
"""r36_inv_corr_ref.py -- R36 designer seat INV-CORR (own file): integer reference of the higher-order stage
(research/panel/R36-INV-CORR.md, "For the builders"), for a builder's stage check.

Stage input: y (19-bit two's complement word left by plain DIV1 + DIV2).  T = y >> 9, L = y & 511,
s = [y < 0], Tf = T XOR (-s) restricted to 9 bits (the sign-folded top), in [0, 511].
  A   = 3 Tf - sum_{i<j<=8} w_ij (Tf_i XOR Tf_j)       (half units)
        w_ij = 3 * 2^(i+j-9) for i + j >= 9, except w_45 = 2;  w_ij = 1 for i + j = 8;  0 below
  h   = (A >> 1) + kappa                 for y >= 0        kappa  = [Tf >= 320]
      = (A >> 1) - kappa - kappa2        for y <  0        kappa2 = [Tf >= 448]   (variant v2; v1: kappa2 = 0)
  U   = y + h                    (on the same 19 wires)
  C   = (L + h) >> 9 in {0, 1, 2};  c0 = C & 1 (the one held wire);  c1 = C >> 1 (erased:
        c1 = NOT c0 AND Tf'8 AND Tf'7 AND NOT L'8 with Tf', L' taken from U)
  guard wire (half unit) ends as Tf_5 and is cleared by one CNOT.

  python -I tools/r36_inv_corr_ref.py [gtab5.json]      self-test (and comparison with the simulated tables)
"""
import json
import sys

import numpy as np


def h_table(kappa_thr=320, neg2=448):
    """tab[sign][Tf], sign 0 for y >= 0.  neg2 = None gives variant v1."""
    Tf = np.arange(512, dtype=np.int64)
    A = 3 * Tf
    for j in range(9):
        for i in range(j):
            if (i, j) == (4, 5):
                w = 2
            elif i + j >= 9:
                w = 3 << (i + j - 9)
            elif i + j == 8:
                w = 1
            else:
                continue
            A = A - w * (((Tf >> i) ^ (Tf >> j)) & 1)
    h = A >> 1
    k = (Tf >= kappa_thr).astype(np.int64)
    k2 = (Tf >= neg2).astype(np.int64) if neg2 else 0
    return np.array([h + k, h - k - k2])


def stage_ref(y, kappa_thr=320, neg2=448):
    """y: int64 array.  Returns U, c0, c1 (before its erase), C."""
    y = np.asarray(y, dtype=np.int64)
    s = (y < 0).astype(np.int64)
    T = y >> 9
    Tf = np.where(s == 1, ~T, T)
    h = h_table(kappa_thr, neg2)[s, Tf]
    C = ((y & 511) + h) >> 9
    return y + h, C & 1, C >> 1, C


def main():
    y = np.arange(-261824, 261379, dtype=np.int64)
    out = {}
    for nm, neg2 in (("v2", 448), ("v1", None)):
        U, c0, c1, C = stage_ref(y, 320, neg2)
        T2 = U >> 9
        s2 = (U < 0).astype(np.int64)
        Tf2 = np.where(s2 == 1, ~T2, T2)
        L8 = (U >> 8) & 1
        f = (1 - c0) & ((Tf2 >> 8) & 1) & ((Tf2 >> 7) & 1) & (1 - L8)
        out[nm] = dict(inputs=len(y), h_rng=[int((U - y).min()), int((U - y).max())],
                       U_rng=[int(U.min()), int(U.max())], C_vals=sorted(int(t) for t in np.unique(C)),
                       c1_formula_ok=bool((f == c1).all()),
                       injective_with_c0=bool(len(np.unique(U * 2 + c0)) == len(y)),
                       injective_without=bool(len(np.unique(U)) == len(y)))
    if len(sys.argv) > 1:
        g = json.load(open(sys.argv[1]))
        for nm, key, neg2 in (("v2", "v2_jtop7", 448), ("v1", "v1_jtop7", None)):
            if key in g:
                t = np.array(g[key]["tab"])
                seen = t > -10 ** 5
                out[nm]["equals_simulated_table"] = bool((t[seen] == h_table(320, neg2)[seen]).all())
                out[nm]["table_entries"] = int(seen.sum())
    print(json.dumps(out))


if __name__ == "__main__":
    main()
