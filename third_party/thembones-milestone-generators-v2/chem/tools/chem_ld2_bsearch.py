"""chem_ld2_bsearch: choose the low basis of one ANGLE_k circuit so that the table becomes DEGENERATE for the block
engine (ours; criterion read from the public loaders' bases, see research/ws_chem_ld2.md).

Facts used [P, ws_chem_ld2]: in the new basis, block b's relation pool is the set of table relations whose low
support lies in the dual block Z_b (3-dim subspace of functionals on the old low coordinates). A pool relation
that vanishes on a point y has its low support in the 2-dim subspace Z_y = Z_b meet ann(y), i.e. it is a
rank-<=2 relation of the whole table r = z1 (x) t1 + z2 (x) t2. Its multiplicity K_Z = 2 nb - rank of the table
projected on Z (x) outputs depends only on the 2-dim subspace Z. A block is good when many of its seven 2-dim
subspaces have K_Z >= 1; the public k1 plan has sum K = 35 over its five blocks, random bases about 1.
Here: hill climbing over five dual blocks forming a direct sum, objective sum_b sum_{Z in Z_b} K_Z.
Only for circuits with one low-cache phase and outputs that see every low coordinate (k = 1..3); coded and
two-phase circuits need the filtration / phase constraints (not done here).

Usage (numba venv): chem_ld2_bsearch.py <table.npz> <k> [--iters N] [--restarts R] [--seed S] [--out basis.json]
The output holds P (rows = new low functions over the old coordinates) for chem_ld2_try --pfile.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_own_compile as cc  # noqa: E402
import chem_ld2_basis as lb2  # noqa: E402


def projections(an, free_top=True):
    """proj[z, M] = the nb-bit vector (parity(X0[M, o] & z))_o for every functional z in [0, 2^nl)."""
    X = an.X0.astype(np.int64)
    nr, nb = X.shape
    rows = nr - 1 if free_top else nr
    n = an.nl
    Z = np.arange(1 << n, dtype=np.int64)
    proj = np.zeros((1 << n, rows), np.int64)
    for o in range(nb):
        for M in range(rows):
            x = int(X[M, o])
            v = Z & x
            # parity of v
            v ^= v >> 8
            v ^= v >> 4
            v ^= v >> 2
            v ^= v >> 1
            proj[:, M] |= (v & 1) << o
    return proj


@njit(cache=True)
def kz(proj, z1, z2, nb):
    """K_Z for Z = span(z1, z2): 2 nb - rank of the rows (proj[z1, M] | proj[z2, M] << nb)."""
    rows = proj.shape[1]
    lo = np.empty(rows, np.int64)
    hi = np.empty(rows, np.int64)
    for M in range(rows):
        lo[M] = proj[z1, M]
        hi[M] = proj[z2, M]
    rank = 0
    # pivot over the 2 nb columns: hi word first then lo word
    for half in range(2):
        for bit in range(nb - 1, -1, -1):
            sel = -1
            for M in range(rank, rows):
                w = hi[M] if half == 0 else lo[M]
                if half == 1 and hi[M] != 0:
                    continue
                if (w >> bit) & 1:
                    sel = M
                    break
            if sel < 0:
                continue
            lo[rank], lo[sel] = lo[sel], lo[rank]
            hi[rank], hi[sel] = hi[sel], hi[rank]
            for M in range(rows):
                if M != rank:
                    w = hi[M] if half == 0 else lo[M]
                    if (w >> bit) & 1:
                        lo[M] ^= lo[rank]
                        hi[M] ^= hi[rank]
            rank += 1
    return 2 * nb - rank


@njit(cache=True)
def block_score(proj, a, b, c, nb):
    s = 0
    s += kz(proj, a, b, nb)
    s += kz(proj, a, c, nb)
    s += kz(proj, b, c, nb)
    s += kz(proj, a, b ^ c, nb)
    s += kz(proj, b, a ^ c, nb)
    s += kz(proj, c, a ^ b, nb)
    s += kz(proj, a ^ b, a ^ c, nb)
    return s


def rank_of(vecs):
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv:
                v ^= piv[h]
            else:
                piv[h] = v
                break
    return len(piv)


def climb(proj, n, nb, rng, iters, start=None):
    """Hill climbing on 5 x 3 dual vectors (direct sum). Returns (score, blocks)."""
    if start is None:
        while True:
            vecs = [rng.randrange(1, 1 << n) for _ in range(15)]
            if rank_of(vecs) == 15:
                break
        blocks = [vecs[3 * i:3 * i + 3] for i in range(5)]
    else:
        blocks = [list(b) for b in start]
    scores = [block_score(proj, *b, nb) for b in blocks]
    tot = sum(scores)
    for it in range(iters):
        i = rng.randrange(5)
        j = rng.randrange(3)
        old = blocks[i][j]
        r = rng.random()
        if r < 0.5:
            new = rng.randrange(1, 1 << n)
        else:                                   # add a vector of another block or of the same block
            i2 = rng.randrange(5)
            j2 = rng.randrange(3)
            if i2 == i and j2 == j:
                continue
            new = old ^ blocks[i2][j2]
        blocks[i][j] = new
        if rank_of([v for b in blocks for v in b]) < 15:
            blocks[i][j] = old
            continue
        s = block_score(proj, *blocks[i], nb)
        if s >= scores[i]:
            tot += s - scores[i]
            scores[i] = s
        else:
            blocks[i][j] = old
    return tot, blocks, scores


@njit(cache=True)
def through(proj, a, nb, n):
    """K of span(a, b) for every b (0 where b in span(a))."""
    out = np.zeros(1 << n, np.int64)
    for b in range(1, 1 << n):
        if b == a:
            continue
        out[b] = kz(proj, a, b, nb)
    return out


@njit(cache=True)
def best_blocks_through(proj, a, cand, nb, keep):
    """Score every block span(a, b, c) for b, c in cand (degenerate partners of a); return the best `keep`
    as rows (score, b, c)."""
    m = len(cand)
    best = np.zeros((keep, 3), np.int64)
    best[:, 0] = -1
    for i in range(m):
        b = cand[i]
        for j in range(i + 1, m):
            c = cand[j]
            if c == (a ^ b):
                continue
            s = block_score(proj, a, b, c, nb)
            # insert if better than the worst kept
            w = 0
            for q in range(1, keep):
                if best[q, 0] < best[w, 0]:
                    w = q
            if s > best[w, 0]:
                best[w, 0] = s
                best[w, 1] = b
                best[w, 2] = c
    return best


def collect_blocks(proj, n, nb, rng, n_a, keep=30):
    """Candidate blocks (score, a, b, c) from n_a random dual vectors a."""
    out = []
    for _ in range(n_a):
        a = rng.randrange(1, 1 << n)
        K = through(proj, a, nb, n)
        cand = np.array(sorted({min(b, b ^ a) for b in np.flatnonzero(K > 0).tolist()}), np.int64)
        if len(cand) < 2:
            continue
        best = best_blocks_through(proj, a, cand, nb, keep)
        for s, b, c in best.tolist():
            if s > 0:
                out.append((s, a, b, c))
    return out


def select_five(blocks, rng, rounds=2000):
    """Five blocks forming a direct sum, maximal total score: greedy, then random single-block exchanges."""
    blocks = sorted(blocks, reverse=True)
    chosen = []
    for blk in blocks:
        vecs = [v for c in chosen for v in c[1:]] + list(blk[1:])
        if rank_of(vecs) == len(vecs):
            chosen.append(blk)
            if len(chosen) == 5:
                break
    if len(chosen) < 5:
        return None
    tot = sum(c[0] for c in chosen)
    top = blocks[:2000]
    for _ in range(rounds):
        i = rng.randrange(5)
        cand = top[rng.randrange(len(top))]
        trial = chosen[:i] + [cand] + chosen[i + 1:]
        vecs = [v for c in trial for v in c[1:]]
        if rank_of(vecs) == 15:
            t = sum(c[0] for c in trial)
            if t > tot:
                chosen, tot = trial, t
    return tot, chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("k", type=int)
    ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--restarts", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--through", type=int, default=0, help="candidate blocks from this many random dual vectors")
    args = ap.parse_args()
    words, _ = cc.load_words(args.table)
    an, lb = cc.analyze_first(words, args.k, cc.feasible_lbs(words), False)
    assert len(an.ymasks) == 1, "one low-cache phase only"
    nb = an.X0.shape[1]
    t0 = time.time()
    proj = projections(an, True)
    n = an.nl
    ident = [[1 << (3 * i), 1 << (3 * i + 1), 1 << (3 * i + 2)] for i in range(5)]
    print(f"k={args.k} nb={nb} proj {proj.shape} ({time.time() - t0:.1f} s); identity-like score "
          f"{sum(block_score(proj, *b, nb) for b in ident)}", flush=True)
    rng = random.Random(args.seed)
    best = None
    if args.through:
        cands = collect_blocks(proj, n, nb, rng, args.through)
        print(f"{len(cands)} candidate blocks from {args.through} dual vectors; best scores "
              f"{sorted((c[0] for c in cands), reverse=True)[:12]} ({time.time() - t0:.1f} s)", flush=True)
        sel = select_five(cands, rng)
        if sel is not None:
            tot, chosen = sel
            blocks = [list(c[1:]) for c in chosen]
            scores = [c[0] for c in chosen]
            print(f"selected: {tot} {scores} ({time.time() - t0:.1f} s)", flush=True)
            if args.iters:
                tot, blocks, scores = climb(proj, n, nb, rng, args.iters, start=blocks)
                print(f"after climb: {tot} {scores} ({time.time() - t0:.1f} s)", flush=True)
            best = (tot, blocks, scores)
    for rs in range(args.restarts):
        tot, blocks, scores = climb(proj, n, nb, rng, args.iters)
        print(f"restart {rs}: score {tot} per block {scores} ({time.time() - t0:.1f} s)", flush=True)
        if best is None or tot > best[0]:
            best = (tot, blocks, scores)
    tot, blocks, scores = best
    # dual blocks -> primal basis: Q has columns = dual vectors (u = 3 i + j), P = Q^-1
    cols = [v for b in blocks for v in b]
    # Q as rows over old t: Q[t] = mask over u with bit u = bit t of cols[u]
    Qrows = [sum(((cols[u] >> t) & 1) << u for u in range(15)) for t in range(15)]
    P = lb2.mat_inv(Qrows, 15)
    out = dict(k=args.k, score=int(tot), scores=[int(s) for s in scores], dual=[[hex(v) for v in b] for b in blocks],
               P=[hex(v) for v in P], seed=args.seed, iters=args.iters, restarts=args.restarts,
               words_sha256=cc.words_sha(words))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh)
    print("best score", tot, scores, flush=True)


if __name__ == "__main__":
    main()
