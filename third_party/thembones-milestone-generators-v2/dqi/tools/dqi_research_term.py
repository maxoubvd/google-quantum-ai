#!/usr/bin/env python
"""Early termination in the homogeneous extraction (workstream dqi_research).

Runs the zero-safe affine extraction (kit convention: inverse of 0 is 0, jet (0, v) -> (0, q v)) and the
homogeneous chain with the kit's pole swap on random instances of every weight, and prints per weight
which groups see D_f = 0 (pole), whether L_f is 0 there, and whether the two extractions give the same heads.
"""
import random
import sys
from dqi_research_model import *


def einv0(u):
    return einv(u) if u else 0


def affine_zero_safe(J):
    heads = []
    for g in range(G):
        u, v = J[g]
        for i in range(2 * g):
            a, b = heads[i]
            q = KAP[i // 2] ^ KAP[g]
            u ^= head_value(a, b, CS[g])
            v ^= b
            z = int(u == 0)
            ui = einv0(u)
            v = emul(v, esq(ui ^ z))
            u = ui
            v = u ^ escale(v, q)
            u = escale(u, q)
        a, b = coords(u, CS[g])
        heads.append((a, b))
        v ^= b
        heads.append(coords(einv0(v), CS[g]))
    return heads


def homogeneous(J):
    """Kit convention: chain from (N0, 1, L0); pole = [D_f == 0] swaps N and D; normalise by D."""
    heads = []
    info = []
    for g in range(G):
        N, D, L = chain_forward(heads, g, J[g][0], 1, J[g][1])
        pole = int(D == 0)
        if pole:
            N, D = D, N
        di = einv(D)
        u = emul(N, di)
        v = emul(L, esq(di))
        info.append((pole, D, L))
        a, b = coords(u, CS[g])
        heads.append((a, b))
        v ^= b
        heads.append(coords(einv0(v), CS[g]))
    return heads, info


if __name__ == '__main__':
    rng = random.Random(7)
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    for weight in range(0, 36):
        for rep in range(reps):
            roots, weights = instance(rng, weight)
            J = jets(roots, weights)
            ha = affine_zero_safe(J)
            hh, info = homogeneous(J)
            poles = ''.join(str(p) for p, _, _ in info)
            lz = ''.join('z' if L == 0 else '.' for _, _, L in info)
            nzh = sum(1 for h in ha if h != (0, 0))
            first_diff = next((i for i in range(34) if ha[i] != hh[i]), None)
            if rep == 0:
                print(f'w={weight:2d} nonzero heads {nzh:2d} poles {poles} Lzero {lz} same_heads {ha == hh} first_diff {first_diff}')
