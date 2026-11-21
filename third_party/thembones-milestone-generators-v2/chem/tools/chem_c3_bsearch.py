#!/usr/bin/env python3
"""chem_c3_bsearch.py: builder C3. A BASIS SEARCH FOR STAGE A of the kept-product loader (design R58-CHEM-KEPT,
section 6: 'stage A at the generic law 150; C2's basis search on the 12-coordinate table would reopen it').

Stage A plans the table on the 13 low points outside D: its coordinates are 12 functionals c (15-bit masks over the
loader's low basis coordinates) that vanish on the three point indicators of D. A block is a 3-dim space U of such
functionals; its cells are (u, o) for u in a basis of U and the 30 outputs, its rows are the 76 block monomials:
X_(c,o) = the 76-bit column 'bit c of X0[M, o]' over M. The engine keeps at most 90 - rank relations in a block
(codim); the table's C-4 condition makes the generic rank 75 (codim 15, stage A (360 - 60) / 2 = 150). A block of
rank 74 allows codim 16: two such blocks give 149, four give 148 [P on the engine's count].
This tool samples random U, keeps those of rank <= 74, picks four blocks forming a direct sum (as many low-rank ones
as it finds), and writes a pfile for tools/r58_chem_kept_c2a.py (P = the new basis in old coordinates, the point
block = block 4, run with --rep 4 --D <the same D>).
  CHEM_BROT=7 python -I chem_c3_bsearch.py PAD.npz --k K --D 13,14,15 --tries 200000 --seed 0 --out pf.json
"""
import argparse
import json
import os
import random
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401

bits = core.bits


@njit(cache=True)
def rank_of(vlo, vhi):
    """GF(2) rank of n vectors of up to 128 bits (lo, hi uint64 halves)."""
    n = vlo.shape[0]
    plo = np.zeros(128, np.uint64)
    phi = np.zeros(128, np.uint64)
    used = np.zeros(128, np.bool_)
    r = 0
    for i in range(n):
        a, b = vlo[i], vhi[i]
        while a != 0 or b != 0:
            if b != 0:
                hb = 64
                x = b
                while x > 1:
                    x >>= np.uint64(1)
                    hb += 1
            else:
                hb = 0
                x = a
                while x > 1:
                    x >>= np.uint64(1)
                    hb += 1
            if used[hb]:
                a ^= plo[hb]
                b ^= phi[hb]
            else:
                used[hb] = True
                plo[hb] = a
                phi[hb] = b
                r += 1
                break
    return r


@njit(cache=True)
def batch_ranks(Ylo, Yhi, C, nout):
    """C: (T, 3) functional masks (15 bits); Y*: (15, nout) column halves. Returns the rank of every block."""
    T = C.shape[0]
    out = np.zeros(T, np.int64)
    vlo = np.zeros(3 * nout, np.uint64)
    vhi = np.zeros(3 * nout, np.uint64)
    for t in range(T):
        for j in range(3):
            c = C[t, j]
            for o in range(nout):
                a = np.uint64(0)
                b = np.uint64(0)
                for i in range(15):
                    if (c >> i) & 1:
                        a ^= Ylo[i, o]
                        b ^= Yhi[i, o]
                vlo[j * nout + o] = a
                vhi[j * nout + o] = b
        out[t] = rank_of(vlo, vhi)
    return out


def gf2_rank(vecs):
    piv = {}
    for v in vecs:
        while v:
            hb = v.bit_length() - 1
            if hb in piv:
                v ^= piv[hb]
            else:
                piv[hb] = v
                break
    return len(piv)


def solve_right_inverse(Cm, n=15):
    """Cm: 12 row masks (functionals). Return 12 column masks p_j with parity(c_i & p_j) = delta_ij."""
    m = len(Cm)
    rows = [(Cm[i], 1 << i) for i in range(m)]           # augmented: functional | identity tag
    piv = []
    for col in range(n):
        sel = next((r for r in range(len(piv), m) if rows[r][0] >> col & 1), None)
        if sel is None:
            continue
        k = len(piv)
        rows[k], rows[sel] = rows[sel], rows[k]
        for r in range(m):
            if r != k and rows[r][0] >> col & 1:
                rows[r] = (rows[r][0] ^ rows[k][0], rows[r][1] ^ rows[k][1])
        piv.append(col)
    assert len(piv) == m, "functionals dependent"
    # rows[k] = (e_{piv[k]} + free cols, combination tag): p_j = sum over k with tag bit j of e_{piv[k]}
    P = [0] * m
    for k, col in enumerate(piv):
        for j in range(m):
            if rows[k][1] >> j & 1:
                P[j] |= 1 << col
    for i in range(m):
        for j in range(m):
            assert bin(Cm[i] & P[j]).count("1") & 1 == (i == j)
    return P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--D", default="13,14,15")
    ap.add_argument("--tries", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--maxhits", type=int, default=64)
    ap.add_argument("--climb", action="store_true", help="also test every single-functional neighbour of each hit")
    a = ap.parse_args()
    t0 = time.time()
    core.BROT, core.J = int(os.environ.get("CHEM_BROT", "7")), 5
    w = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    an = core.analyze(w, a.k, 4)
    X0 = np.asarray(an.X0, np.int64)
    nr, nout = X0.shape
    nl = an.nl
    assert nl == 15 and nr <= 128
    # columns over the block monomials: Y[i][o] = bit i of X0[M, o] over M
    Y = [[0] * nout for _ in range(nl)]
    for M in range(nr):
        for o in range(nout):
            v = int(X0[M, o])
            for i in bits(v):
                Y[i][o] |= 1 << M
    Ylo = np.array([[Y[i][o] & ((1 << 64) - 1) for o in range(nout)] for i in range(nl)], np.uint64)
    Yhi = np.array([[Y[i][o] >> 64 for o in range(nout)] for i in range(nl)], np.uint64)
    # point indicators in basis coordinates; their annihilator (12-dim)
    LB = []
    for bs in an.basis:
        tt = 0
        for l in range(16):
            v = 0
            for m in bits(bs["anf"]):
                if (l & m) == m:
                    v ^= 1
            tt |= v << l
        LB.append(tt)
    sp = core.Span()
    for t in range(nl):
        sp.add(LB[t], 1 << t)
    D = [int(x) for x in a.D.split(",")]
    mD = []
    for p in D:
        r, c = sp.reduce(1 << p)
        assert r == 0
        mD.append(c)
    ann = [c for c in range(1, 1 << nl) if all(bin(c & d).count("1") % 2 == 0 for d in mD)]
    abasis = []
    for c in ann:
        if gf2_rank(abasis + [c]) > len(abasis):
            abasis.append(c)
        if len(abasis) == 12:
            break
    assert len(abasis) == 12
    rng = np.random.default_rng(a.seed)

    def rand_func():
        m = int(rng.integers(1, 1 << 12))
        c = 0
        for i in range(12):
            if m >> i & 1:
                c ^= abasis[i]
        return c
    hits = {}
    batch = 4096
    done = 0
    hist = {}
    while done < a.tries and len(hits) < a.maxhits:
        C = np.zeros((batch, 3), np.int64)
        for t in range(batch):
            while True:
                tri = [rand_func(), rand_func(), rand_func()]
                if gf2_rank(tri) == 3:
                    break
            C[t] = tri
        R = batch_ranks(Ylo, Yhi, C, nout)
        for t in range(batch):
            hist[int(R[t])] = hist.get(int(R[t]), 0) + 1
            if R[t] <= 74:
                # canonical form of U: reduced echelon of the three masks
                rows = sorted(int(x) for x in C[t])
                ech = []
                for v in rows:
                    for e in ech:
                        if v ^ e < v:
                            v ^= e
                    if v:
                        ech.append(v)
                key = tuple(sorted(ech))
                hits[key] = int(R[t])
        done += batch
    print(f"k {a.k} D {D}: {done} blocks sampled, rank histogram {dict(sorted(hist.items()))}; "
          f"{len(hits)} distinct blocks of rank <= 74; {time.time() - t0:.0f} s", flush=True)
    if a.climb:
        # every single-functional replacement of every hit (3 x 4,095 neighbours): blocks of rank <= 73
        allf = []
        for m in range(1, 1 << 12):
            c = 0
            for i in range(12):
                if m >> i & 1:
                    c ^= abasis[i]
            allf.append(c)
        low = {}
        hh = {}
        for key in list(hits):
            tri = list(key)
            cand = []
            for j in range(3):
                keep = [tri[x] for x in range(3) if x != j]
                for c in allf:
                    if gf2_rank(keep + [c]) == 3:
                        cand.append(keep + [c])
            C = np.array(cand, np.int64)
            R = batch_ranks(Ylo, Yhi, C, nout)
            for t in range(len(cand)):
                hh[int(R[t])] = hh.get(int(R[t]), 0) + 1
                if R[t] <= 73:
                    low[tuple(sorted(int(x) for x in C[t]))] = int(R[t])
        print(f"CLIMB k {a.k}: {sum(hh.values())} neighbours of {len(hits)} hits, rank histogram "
              f"{dict(sorted(hh.items()))}; {len(low)} blocks of rank <= 73; {time.time() - t0:.0f} s", flush=True)
        for key, r in low.items():
            hits[key] = r
    # choose blocks: greedy over hits (lowest rank first), keep a direct sum, then fill with random blocks
    hl = sorted(hits.items(), key=lambda kv: kv[1])
    best = None
    pyr = random.Random(a.seed)
    for attempt in range(400):
        order = list(hl)
        pyr.shuffle(order)
        order.sort(key=lambda kv: kv[1])
        chosen, span = [], []
        for key, r in order:
            if len(chosen) == 4:
                break
            if gf2_rank(span + list(key)) == len(span) + 3:
                chosen.append((list(key), r))
                span += list(key)
        good = len(chosen)
        while len(chosen) < 4:
            tri = [rand_func(), rand_func(), rand_func()]
            if gf2_rank(span + tri) == len(span) + 3:
                chosen.append((tri, None))
                span += tri
        if best is None or good > best[0]:
            best = (good, chosen)
        if good == 4:
            break
    good, chosen = best
    Cm = [c for blk, _ in chosen for c in blk]
    # exact ranks of the chosen blocks
    rk = []
    for blk, _ in chosen:
        vecs = []
        for c in blk:
            for o in range(nout):
                v = 0
                for i in bits(c):
                    v ^= Y[i][o]
                vecs.append(v)
        rk.append(gf2_rank(vecs))
    P12 = solve_right_inverse(Cm)
    P = P12 + mD                                  # columns: 12 stage-A basis vectors, then the point block
    json.dump(dict(k=a.k, D=D, P=[hex(x) for x in P], blocks=[[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11],
                                                               [12, 13, 14]],
                   ranks=rk, codim_bound=[90 - r for r in rk], hits=len(hits), sampled=done,
                   mode="C3 bsearch (stage-A blocks of low rank; point block = block 4)"), open(a.out, "w"))
    print(f"pfile {a.out}: block ranks {rk} (codim bounds {[90 - r for r in rk]}), {good} low-rank blocks in a direct "
          f"sum; run r58_chem_kept_c2a.py --pfile {a.out} --rep 4 --D {a.D}", flush=True)


if __name__ == "__main__":
    main()
