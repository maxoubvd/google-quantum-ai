#!/usr/bin/env python3
"""invsqrt_pp_exlinks.py -- workstream invsqrt_pp (P1): EXACT check of candidate links (from invsqrt_pp_links.py)
on all 3,133,452 inputs, plus the weight of every product involved.

For each candidate (i, j): reduce d_j + d_i modulo the exact affine span of the 46 wire vectors just before j
(Gaussian elimination on the full vectors, padded lanes masked). Writes the links file again with a 4th field
exact (1 / 0) and prints how many are exact and how many of those keep a zero wire throughout.
Usage: python invsqrt_pp_exlinks.py <circuit.kmb> <links.json> <out.json>
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_links import asap_order, permute  # noqa: E402
from invsqrt_pp_hold import prod_of  # noqa: E402

LASTMASK = np.uint64((1 << (sim.NX - 64 * (sim.W - 1))) - 1)


def masked(v):
    v = v.copy()
    v[-1] &= LASTMASK
    return v


def exact_reduce(S, T):
    """Reduce rows of T (k, W) modulo span{1, rows of S}; returns reduced copy (padded lanes masked)."""
    M = np.concatenate([np.full((1, sim.W), sim.ONES, np.uint64), S], axis=0).copy()
    M[:, -1] &= LASTMASK
    T = T.copy()
    T[:, -1] &= LASTMASK
    nr = M.shape[0]
    for r in range(nr):
        v = M[r]
        nz = np.flatnonzero(v)
        if len(nz) == 0:
            continue
        w = int(nz[0])
        bb = np.uint64(int(v[w]) & -int(v[w]))
        for s in range(r + 1, nr):
            if M[s, w] & bb:
                M[s] ^= v
        hit = (T[:, w] & bb) != 0
        if hit.any():
            T[hit] ^= v
    return T


def popcount(v):
    return int(np.unpackbits(masked(v).view(np.uint8)).sum())


def main():
    t0 = time.time()
    path, lpath, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    L = json.load(open(lpath))
    C0 = sim.load(path)
    C = permute(C0, asap_order(C0)[0]) if L['asap'] else C0
    links = L['links']
    need_i = set(a for a, b, f in links)
    by_j = {}
    for k, (a, b, f) in enumerate(links):
        by_j.setdefault(b, []).append(k)
    D = {}
    wt = {}
    exact = [0] * len(links)

    def hook(i, S, P, bits, act):
        if i in need_i or i in by_j:
            d = prod_of(C, i, S)
            wt[i] = popcount(d)
            if i in need_i:
                D[i] = d
            if i in by_j:
                ks = by_j[i]
                T = np.stack([d ^ D[links[k][0]] for k in ks], axis=0)
                R = exact_reduce(S, T)
                for t, k in enumerate(ks):
                    exact[k] = int(not R[t].any())
    sim.run(C, hook=hook)
    out = []
    for k, (a, b, f) in enumerate(links):
        out.append((a, b, f, exact[k], wt.get(a, -1), wt.get(b, -1)))
    L2 = dict(L)
    L2['links'] = out
    json.dump(L2, open(outp, 'w'))
    ex = [x for x in out if x[3]]
    print('candidates', len(out), 'exact', len(ex), 'exact and fitting', sum(1 for x in ex if x[2] >= 1),
          'distinct j exact', len(set(x[1] for x in ex)), round(time.time() - t0, 1), 's')
    sparse = [x for x in out if not x[3]]
    if sparse:
        w = np.array([min(x[4], x[5]) for x in sparse])
        print('inexact candidates: product weight min/median/max', w.min(), int(np.median(w)), w.max())
    if ex:
        w = np.array([min(x[4], x[5]) for x in ex])
        print('exact links: product weight min/median/max', w.min(), int(np.median(w)), w.max())
        ln = np.array([x[1] - x[0] for x in ex])
        print('exact links: interval length pct 10/50/90', np.percentile(ln, [10, 50, 90]))


if __name__ == '__main__':
    main()
