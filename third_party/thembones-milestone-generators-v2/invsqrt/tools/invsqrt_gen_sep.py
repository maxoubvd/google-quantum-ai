"""invsqrt_gen: choose a few GF(2)-linear forms of the lo table bits that, together with the c2 table, tell
all 256 addresses apart (own tool). Needed by variant g4: before A^-1 all lo wires are erased except these
forms, and the erased wires must be functions of what stays.

separating_forms(C2, C0, nlo, m) -> list of m bit masks over lo wire indices 0..nlo-1 (lo wire i = bit i+1 of C0),
or None. Deterministic: forms 1..m-1 are chosen greedily among single wires and seeded random masks (most
colliding pairs split), the last form solves the linear system "f . delta = 1" for every pair still colliding.
"""
from __future__ import annotations

import random


def _pairs(keys):
    groups = {}
    for a, k in enumerate(keys):
        groups.setdefault(k, []).append(a)
    out = []
    for g in groups.values():
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                out.append((g[i], g[j]))
    return out


def _solve_all_ones(deltas, nbits):
    """f with parity(f & d) == 1 for every d in deltas, or None."""
    rows = [(d, 1) for d in deltas]
    piv = {}
    for d, r in rows:
        while d:
            h = d.bit_length() - 1
            if h in piv:
                d ^= piv[h][0]
                r ^= piv[h][1]
            else:
                piv[h] = (d, r)
                break
        else:
            if r:
                return None
    f = 0
    for h in sorted(piv):                      # back-substitution, free variables = 0
        d, r = piv[h]
        low = d & ~(1 << h)
        if bin(low & f).count("1") % 2 != r:
            f |= 1 << h
    for d in deltas:
        assert bin(d & f).count("1") % 2 == 1
    return f


def separating_forms(C2, C0, nlo, m, seed=1, tries=4000):
    vec = [(C0[a] >> 1) & ((1 << nlo) - 1) for a in range(256)]
    rng = random.Random(f"{seed}:sep")
    cands = [1 << i for i in range(nlo)] + [rng.getrandbits(nlo) or 1 for _ in range(tries)]
    forms = []
    keys = [(C2[a],) for a in range(256)]
    for step in range(m):
        pairs = _pairs(keys)
        if not pairs:
            break
        if step == m - 1:
            f = _solve_all_ones([vec[a] ^ vec[b] for a, b in pairs], nlo)
            if f is None:
                return None
        else:
            f = max(cands, key=lambda g: sum(bin((vec[a] ^ vec[b]) & g).count("1") & 1 for a, b in pairs))
        forms.append(f)
        keys = [k + (bin(vec[a] & f).count("1") & 1,) for a, k in enumerate(keys)]
    if _pairs(keys):
        return None
    return forms


if __name__ == "__main__":
    import sys
    sys.path.insert(0, "<HOME>/Quantum/tools")
    import invsqrt_gen_model as M
    P = M.Params(glo=2, gu=1, fc2=1, virt=1, tie_c1=1)
    tabs = M.fit_tables(P)
    print("colliding pairs under c2 alone:", len(_pairs([(v,) for v in tabs["C2"]])),
          "distinct c2 values:", len(set(tabs["C2"])))
    for m in (1, 2, 3, 4):
        f = separating_forms(tabs["C2"], tabs["C0"], P.wlo - 1, m)
        print(m, "forms:", None if f is None else [bin(x).count("1") for x in f])
