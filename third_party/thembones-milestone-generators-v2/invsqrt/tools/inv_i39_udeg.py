#!/usr/bin/env python3
"""inv_i39_udeg.py -- I39: the LOWEST ANF degree of the slot's parity bit u(T) = hn(T) mod 2 (or b = u xor t) over the
T that occur (the unseen T are don't-cares): GF(2) solve with all monomials of degree <= d in the 11 bits of T.
(07:28Z's 'degree 10' filled the 7 unseen T with neighbours: an artefact.)   inv_i39_udeg.py HTAB.npz [TABLE.json]
Safe: prints only."""
import itertools
import json
import sys

import numpy as np

hz = np.load(sys.argv[1])
h0 = hz["h"].astype(np.int64)
seen = hz["seen"].astype(bool)
Tt = np.arange(-1024, 1024)
s = (Tt < 0).astype(np.int64)
Tf = np.where(s == 1, ~Tt, Tt) & 1023
hn = h0.copy()
if len(sys.argv) > 2:
    D = json.load(open(sys.argv[2]))["S1"]["tables"]
    hn = h0 - np.array(D[0])[(Tf >> 1) & 127] + np.where(s == 1, -1, 1) * np.array(D[1])[Tf >> 4]
bits = [((Tt >> i) & 1) for i in range(11)]           # T0..T10 (T10 = sign)
idx = np.flatnonzero(seen)
for name, f in (("u", hn & 1), ("b (par 0)", (hn & 1) ^ (Tt & 1))):
    for d in range(1, 7):
        mons = [()] + [c for k in range(1, d + 1) for c in itertools.combinations(range(11), k)]
        A = np.zeros((len(idx), len(mons) + 1), dtype=np.uint8)
        for j, mset in enumerate(mons):
            col = np.ones(len(Tt), dtype=np.int64)
            for q in mset:
                col = col & bits[q]
            A[:, j] = col[idx]
        A[:, -1] = f[idx]
        M = A.copy()
        piv = 0
        for col in range(M.shape[1] - 1):
            sel = np.flatnonzero(M[piv:, col]) + piv
            if len(sel) == 0:
                continue
            M[[piv, sel[0]]] = M[[sel[0], piv]]
            oth = np.flatnonzero(M[:, col])
            oth = oth[oth != piv]
            M[oth] ^= M[piv]
            piv += 1
            if piv == M.shape[0]:
                break
        bad = int(((M[:, :-1].sum(axis=1) == 0) & (M[:, -1] == 1)).sum())
        print("%s: degree <= %d over %d occurring T: %s" % (name, d, len(idx), "FITS" if bad == 0 else
                                                             "no (%d inconsistent)" % bad), flush=True)
        if bad == 0:                                   # one solution (free monomials 0): its monomials by degree,
            sol = np.zeros(len(mons), dtype=np.uint8)  # in the folded variables Tf0..9 and s as well
            for r_ in range(piv):
                nz = np.flatnonzero(M[r_, :-1])
                if len(nz):
                    sol[nz[0]] = M[r_, -1]
            cnt = {}
            for j, mset in enumerate(mons):
                if sol[j]:
                    cnt[len(mset)] = cnt.get(len(mset), 0) + 1
            print("   one solution in T0..T10: monomials by degree %s" % dict(sorted(cnt.items())), flush=True)
            break
