"""chem_ld2_rel: the relation space R = F-perp of a dissected loader plan, and whether it is block-diagonal.

Usage: python -I chem_ld2_rel.py <dissect.json>
F = span of the features y (x) s (15 low monomials x 40 word bits). Prints dim F, the outputs used, dim R inside
the used outputs, and for every 3-subset S of low monomials dim(R restricted to rows S); then the best disjoint
cover by five 3-subsets.
"""
import itertools
import json
import sys

NOUT = 40


def echelon(vecs):
    piv = {}
    for v in vecs:
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                v ^= piv[hb]
            else:
                piv[hb] = v
                break
    return piv


def reduce(piv, v):
    while v:
        hb = v.bit_length() - 1
        if hb in piv:
            v ^= piv[hb]
        else:
            return v
    return 0


def nullspace(vecs, n):
    """Basis of {x in GF(2)^n : <x, v> = 0 for all v in vecs} (Gaussian elimination on bit columns)."""
    rows = [v for v in vecs]
    piv_cols = []
    r = 0
    rows = list(rows)
    m = len(rows)
    where = {}
    for col in range(n):
        sel = None
        for i in range(r, m):
            if (rows[i] >> col) & 1:
                sel = i
                break
        if sel is None:
            continue
        rows[r], rows[sel] = rows[sel], rows[r]
        for i in range(m):
            if i != r and (rows[i] >> col) & 1:
                rows[i] ^= rows[r]
        where[col] = r
        piv_cols.append(col)
        r += 1
    free = [c for c in range(n) if c not in where]
    basis = []
    for f in free:
        x = 1 << f
        for col, i in where.items():
            if (rows[i] >> f) & 1:
                x |= 1 << col
        basis.append(x)
    return basis


def main():
    d = json.load(open(sys.argv[1]))
    prods, fan = d["prods"], d["fan"]
    used = 0
    for f in fan:
        used |= f
    outs = [t for t in range(NOUT) if (used >> t) & 1]
    no = len(outs)
    oi = {t: i for i, t in enumerate(outs)}
    print("outputs used:", no)

    def feat(y, s):
        v = 0
        sc = 0
        for t in outs:
            if (s >> t) & 1:
                sc |= 1 << oi[t]
        for m in range(1, 16):
            if (y >> m) & 1:
                v |= sc << (no * (m - 1))
        return v
    F = []
    for p, s in zip(prods, fan):
        if p["l2"] and p["h1"]:
            F.append(feat(p["l2"], s))
        if p["l1"] and p["h2"]:
            F.append(feat(p["l1"], s))
    n = 15 * no
    piv = echelon(F)
    print("features", len(F), "dim F", len(piv), "cells", n)
    R = nullspace(list(piv.values()), n)
    print("dim R", len(R))
    # dim of R restricted to rows in S: relations supported on rows S = R intersect span(e_(m,t), m in S)
    # = vectors of R with zero on rows outside S: nullspace computation on coordinates outside S.
    rowmask = [((1 << no) - 1) << (no * (m - 1)) for m in range(1, 16)]

    def dim_supported(S):
        outside = 0
        for m in range(1, 16):
            if m not in S:
                outside |= rowmask[m - 1]
        # R vectors as combinations; find combos with zero outside: rank(R) - rank(R projected outside)
        proj = [r & outside for r in R]
        return len(R) - len(echelon(proj))
    best = {}
    for S in itertools.combinations(range(1, 16), 3):
        best[S] = dim_supported(set(S))
    top = sorted(best.items(), key=lambda kv: -kv[1])[:20]
    print("top 3-subsets by supported relation dim:", [([format(m, '04b') for m in S], v) for S, v in top])
    # best disjoint cover (greedy over sorted + exhaustive over partitions would be large; do a DP over masks)
    from functools import lru_cache
    subs = [(sum(1 << (m - 1) for m in S), v) for S, v in best.items()]

    @lru_cache(maxsize=None)
    def dp(mask):
        if mask == (1 << 15) - 1:
            return 0, ()
        first = next(i for i in range(15) if not (mask >> i) & 1)
        res = (-1, ())
        for sm, v in subs:
            if (sm >> first) & 1 and not (sm & mask):
                sub = dp(mask | sm)
                if sub[0] + v > res[0]:
                    res = (sub[0] + v, (sm,) + sub[1])
        return res
    tot, parts = dp(0)
    print("best partition into 3-subsets: sum of supported dims", tot,
          [[format(m + 1, '04b') for m in range(15) if (sm >> m) & 1] for sm in parts])


if __name__ == "__main__":
    main()
