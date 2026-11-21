#!/usr/bin/env python3
"""chem_c4_flips.py: builder C4. Census of a full loader decomposition (tools/chem_c4_dump.py output):
  * exactness: every output's mixed table == XOR of the R_p of the products that fan into it;
  * the flips available among the products of chosen groups (Kauers-Moosbauer moves in this model):
      F1  rank(R_p + R_q) <= 2           -> R_p += R_q, s_q += s_p   (keeps the sum; a product-level flip)
      F2  shared low vector (planes meet) -> half-product flip (a (x) y (x) s + a' (x) y (x) s')
      F3  shared high vector (A's meet)   -> half-product flip
      F4  equal fan-outs                  -> the two products' R may be redistributed (merge if rank <= 2)
  * --grp g: writes group g in the e19grp dump format (wb, a, b, prods) for r59_chem_groups_lab.py.

  python -I chem_c4_flips.py d_k2.json --groups 0,1 [--grp 0 --grpout g0.json]
"""
import argparse
import json
import sys

import numpy as np


def bitsv(x, n):
    return np.array([(x >> i) & 1 for i in range(n)], np.uint8)


def gf2_rank(M):
    M = (M.copy() & 1).astype(np.uint8)
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if M[i, c]:
                piv = i
                break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def load(path):
    d = json.load(open(path))
    nr, nl = d["nr"], d["nl"]
    R = []
    for (y1, A1, y2, A2, s, ph) in d["prods"]:
        m = np.outer(bitsv(A1, nr), bitsv(y1, nl))
        if y2:
            m ^= np.outer(bitsv(A2, nr), bitsv(y2, nl))
        R.append(m & 1)
    X = np.array([[bitsv(v, nl) for v in row] for row in d["X0"]], np.uint8)   # nr x nb x nl
    return d, R, X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dump")
    ap.add_argument("--groups", default="0,1", help="indices among the generic groups")
    ap.add_argument("--grp", type=int, default=-1)
    ap.add_argument("--grpout", default="")
    a = ap.parse_args()
    d, R, X = load(a.dump)
    nr, nl, nb = d["nr"], d["nl"], d["nb"]
    # exactness of the whole decomposition
    acc = np.zeros((nr, nb, nl), np.uint8)
    for (y1, A1, y2, A2, s, ph), m in zip(d["prods"], R):
        for o in range(nb):
            if (s >> o) & 1:
                acc[:, o, :] ^= m
    bad = int((acc != X).sum())
    print(f"k {d['k']}: {len(R)} products, nr {nr}, nl {nl}, outputs {nb}; exact: {'YES' if bad == 0 else 'NO (' + str(bad) + ' cells)'}")
    gens = [b for b in d["bounds"] if b["tag"] == "generic"]
    for i, b in enumerate(gens):
        print(f"  generic group {i}: outs {b['outs']} products {b['start']}..{b['end'] - 1} ({b['made']}) {b['log']}")
    sel = [int(x) for x in a.groups.split(",")]
    idx = [p for g in sel for p in range(gens[g]["start"], gens[g]["end"])]
    outs = [o for g in sel for o in gens[g]["outs"]]
    omask = sum(1 << o for o in outs)
    P = d["prods"]
    # fan-out bits outside the selected outputs (re-fan-outs into later groups: the pool use)
    ext = sum(1 for p in idx if P[p][4] & ~omask)
    print(f"selected groups {sel}: {len(idx)} products on {len(outs)} outputs; {ext} of them also fan out into later "
          f"groups (pool use)")
    # flips
    f1 = f2 = f3 = f4 = 0
    low = {}
    for p in idx:
        y1, A1, y2, A2, s, ph = P[p]
        ys = {y1, y2, y1 ^ y2} if y2 else {y1}
        for y in ys:
            low.setdefault(y, []).append(p)
    f2 = sum(len(v) * (len(v) - 1) // 2 for v in low.values())
    high = {}
    for p in idx:
        y1, A1, y2, A2, s, ph = P[p]
        As = {A1, A2, A1 ^ A2} if y2 else {A1}
        for A in As:
            high.setdefault(A, []).append(p)
    f3 = sum(len(v) * (len(v) - 1) // 2 for v in high.values())
    fan = {}
    for p in idx:
        fan.setdefault(P[p][4] & omask, []).append(p)
    f4 = sum(len(v) * (len(v) - 1) // 2 for v in fan.values())
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            if gf2_rank(R[idx[i]] ^ R[idx[j]]) <= 2:
                f1 += 1
    rank1 = sum(1 for p in idx if not P[p][2])
    print(f"FLIPS among {len(idx)} products ({len(idx) * (len(idx) - 1) // 2} pairs): F1 rank(Rp+Rq)<=2: {f1}; "
          f"F2 shared low vector: {f2}; F3 shared high vector: {f3}; F4 equal fan-out (on the selected outputs): {f4}; "
          f"rank-1 products: {rank1}")
    if a.grp >= 0 and a.grpout:
        g = gens[a.grp]
        if g["start"] == g["end"]:
            sys.exit("empty group")
        lc = list(range(nl))
        b = len(g["outs"])
        rows = []
        for M in range(nr):
            v = 0
            for gi, o in enumerate(g["outs"]):
                for ci in range(nl):
                    if X[M, o, ci]:
                        v |= 1 << (gi * nl + ci)
            rows.append(v)
        piv, wb = {}, []
        for v in rows:
            w = v
            while w:
                h = w.bit_length() - 1
                if h in piv:
                    w ^= piv[h]
                else:
                    piv[h] = w
                    wb.append(v)
                    break
        oi = {o: i for i, o in enumerate(g["outs"])}
        prods = [[P[p][0], P[p][2], sum(1 << oi[o] for o in oi if (P[p][4] >> o) & 1)] for p in range(g["start"], g["end"])]
        json.dump(dict(groups=[dict(outs=g["outs"], coords=lc, a=nl, b=b, dim_w=len(wb), wb=[hex(v) for v in wb],
                                    prods=prods, made=g["made"])]), open(a.grpout, "w"))
        print(f"group {a.grp} written to {a.grpout}: b {b} a {nl} dim_w {len(wb)} (residual = X0: valid for a group "
              f"without kills only)")


if __name__ == "__main__":
    main()
