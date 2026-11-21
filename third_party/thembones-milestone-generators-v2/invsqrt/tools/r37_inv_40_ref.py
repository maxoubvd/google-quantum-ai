#!/usr/bin/env python3
"""r37_inv_40_ref.py -- R37 designer INV-40 (own code): integer reference of the higher-order stage H1' for a root
one digit shorter (m = 7; nb = 10 folded bits of T) -- the R36 / FC2 stage (tools/r36_inv_corr_ref.py, read as a
model) generalised to nb bits, with a cubic term made of pieces that are free in the circuit.

Input y: (27 - m)-bit two's complement word left by plain DIV1 + DIV2 (20 wires at m = 7).
  T = y >> 9 (nb + 1 bits signed), L = y & 511, s = [y < 0], Tf = T xor (-s) on nb bits (sign-folded top).
  A   = lin Tf - sum_{i<j<nb} w_ij (Tf_i xor Tf_j)     half units;  lin = 3 * 2^(nb-9)   (3 at nb 9, 6 at nb 10)
        w_ij = 3 * 2^(i+j-9) for i + j >= 9 (w_45 = 2);  w_ij = 1 for i + j = 8 (unless in `drop8`);  0 below
  cub = sum_k g_k P_k(Tf)          units; P_k an AND of folded bits (spec "cub": list of [weight, [bits]]);
  h   = (A + K2) >> 1 + cub   for y >= 0;   (A + K2) >> 1 - cub   for y < 0      (K2: constant, half units)
  U   = y + h;   C = (L + h) >> 9;   c0 = C & 1 is the one held wire, the upper bits of C are rebuilt from U.
nb = 9, cub = [[1, [8]]], K2 = 0 is R36's v00 (h_table(320, None) of the R36 file with kappa at 256).

  r37_inv_40_ref.py selftest                      R36 v00 reproduced; m 7 table ranges
  r37_inv_40_ref.py carry LO HI 'SPEC'            carry analysis on y in [LO, HI]: values of C, injectivity with
                                                   c0, and the smallest supports of C >> 1 in the output bits
"""
import itertools
import json
import sys

import numpy as np

DEF7 = dict(nb=10, cub=[[4, [9]], [8, [9, 8]], [2, [9, 7]], [1, [9, 6]]], K2=0, drop8=[])


def h_table(spec):
    """tab[sign][Tf] (units), sign 0 for y >= 0."""
    nb = int(spec["nb"])
    Tf = np.arange(1 << nb, dtype=np.int64)
    A = (3 << (nb - 9)) * Tf
    drop8 = {tuple(p) for p in spec.get("drop8", [])}
    for j in range(nb):
        for i in range(j):
            if (i, j) == (4, 5):
                w = 2
            elif i + j >= 9:
                w = 3 << (i + j - 9)
            elif i + j == 8 and (i, j) not in drop8:
                w = 1
            elif (i, j) in {tuple(q) for q in spec.get("cin_pairs", [])}:
                w = 1                              # a dropped pair put back as a half-unit borrow-in
            else:
                continue
            A = A - w * (((Tf >> i) ^ (Tf >> j)) & 1)
    for b in spec.get("cin_bits", []):               # borrow-ins that are folded bits themselves (half a unit each)
        A = A - ((Tf >> int(b)) & 1)
    cub = np.zeros_like(Tf)
    for g, bits in spec.get("cub", []):
        p = np.ones_like(Tf)
        for b in bits:
            p = p & ((Tf >> b) & 1)
        cub = cub + int(g) * p
    if spec.get("cubtab"):                           # cubic as a table of the top k folded bits: [k, [2^k values]]
        k, vals = spec["cubtab"]
        cub = cub + np.array(vals, dtype=np.int64)[Tf >> (nb - int(k))]
    h = (A + int(spec.get("K2", 0))) >> 1
    return np.array([h + cub, h - cub])


def stage_ref(y, spec):
    """y: int64 array. Returns U, C."""
    y = np.asarray(y, dtype=np.int64)
    s = (y < 0).astype(np.int64)
    T = y >> 9
    Tf = np.where(s == 1, ~T, T)
    h = h_table(spec)[s, Tf]
    C = ((y & 511) + h) >> 9
    return y + h, C


def fold_out(U, nb):
    T2 = U >> 9
    s2 = (U < 0).astype(np.int64)
    return np.where(s2 == 1, ~T2, T2) & ((1 << nb) - 1), U & 511, s2


def carry_report(lo, hi, spec):
    y = np.arange(lo, hi + 1, dtype=np.int64)
    nb = int(spec["nb"])
    U, C = stage_ref(y, spec)
    c0 = C & 1
    out = dict(inputs=len(y), h_rng=[int((U - y).min()), int((U - y).max())], U_rng=[int(U.min()), int(U.max())],
               C_vals=sorted(int(t) for t in np.unique(C)),
               injective_with_c0=bool(len(np.unique(U * 2 + c0)) == len(y)),
               injective_without=bool(len(np.unique(U)) == len(y)))
    Tf2, L2, s2 = fold_out(U, nb)
    hi_c = C >> 1
    sup = []
    for kt in range(nb, 2, -1):                # keep the top kt bits of Tf'
        for kl in range(0, 5):                 # and the top kl bits of L'
            key = ((Tf2 >> (nb - kt)) << 8) | ((L2 >> (9 - kl)) << 2 if kl else 0)
            key = (key << 2) | (c0 << 1) | s2
            o = np.argsort(key, kind="stable")
            ks, vs = key[o], hi_c[o]
            brk = np.flatnonzero(np.diff(ks)) + 1
            mn = np.minimum.reduceat(vs, np.r_[0, brk])
            mx = np.maximum.reduceat(vs, np.r_[0, brk])
            if (mn == mx).all():
                sup.append(dict(Tf_top_bits=kt, L_top_bits=kl, with_sign=True, classes=int(len(mn))))
    out["supports_of_C_hi(Tf' top bits, L' top bits, c0, sign)"] = sup[-6:]
    return out


def main():
    if sys.argv[1] == "selftest":
        import importlib.util
        import os
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "r36_inv_corr_ref.py")
        sp = importlib.util.spec_from_file_location("r36ref", p)
        R36 = importlib.util.module_from_spec(sp)
        sp.loader.exec_module(R36)
        mine = h_table(dict(nb=9, cub=[[1, [8]]], K2=0))
        print(json.dumps(dict(r36_v00_equal=bool((mine == R36.h_table(256, None)).all()))))
        t7 = h_table(DEF7)
        Tf = np.arange(1024)
        ideal = 3 * (Tf + 0.5) ** 2 / 1024.0
        cube = 5 * ((Tf + 0.5) * 512.0) ** 3 / 2.0 ** 55
        print(json.dumps(dict(m7_h_pos=[int(t7[0].min()), int(t7[0].max())], m7_h_neg=[int(t7[1].min()), int(t7[1].max())],
                              err_pos=[float((t7[0] - ideal - cube).min()), float((t7[0] - ideal - cube).max())],
                              err_neg=[float((t7[1] - ideal + cube).min()), float((t7[1] - ideal + cube).max())])))
    elif sys.argv[1] == "carry":
        spec = dict(DEF7)
        if len(sys.argv) > 4:
            spec.update(json.loads(sys.argv[4]))
        print(json.dumps(carry_report(int(sys.argv[2]), int(sys.argv[3]), spec)))


if __name__ == "__main__":
    main()
