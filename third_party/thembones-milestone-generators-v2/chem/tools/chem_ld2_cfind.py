r"""chem_ld2_cfind: dual blocks of a CODED circuit (k >= 4) with large relation pools (ours).

Measured (ws_chem_ld2.md 15:15Z): the public coded blocks hold no planar relations; their pools are larger than
ours ([17,17,19,18,14] against [17,17,17,17,11] on k10, same 454 cells) because the table rows restricted to a
block have rank below 75: every unit of rank deficiency is one more relation to choose from.

Model.  Dual vectors d_u (functionals on the old low coordinates) give the new coordinates; output o (old support
S_o) has a cell at u iff d_u & S_o != 0.  The supports are nested, 0xf (linear) < 0x3ff (degree <= 2) < all, so the
cell count is minimal (as in the identity basis) iff the duals fall into three classes: class 1 touches 0xf
(4 of them, independent on 0xf), class 2 touches 0x3ff but not 0xf (6, independent on 0x3f0), class 3 lives on
0x7c00 (5, a basis there).  Then Q is block lower triangular and P = Q^-1 is a valid filtered basis.
A block is a class pattern such as (1,2,3) or (2,2,2); its pool = cells - rank of the 75 (free top row) rows.
Search: coordinate scans (all admissible values of one dual vector, the block's pool computed for each), from
random starts, for a list of block patterns.

Usage (numba venv): chem_ld2_cfind.py <table.npz> <k> [--patterns "123,123,123,123,223;133,133,122,123,222"]
                    [--rounds 40] [--restarts 3] [--seed S] [--out basis.json]
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
import chem_loaders_core as core  # noqa: E402
from chem_ld2_bsearch import rank_of  # noqa: E402
# fifteenth holder (17:30Z): the carrier-free build (env C2_NOCAR=1) needs my width module's hook here too; every
# other run is unchanged (this module never imported it)
if os.environ.get("C2_NOCAR") == "1":
    import chem_ld2_wwidth  # noqa: E402,F401

M1, M2, M3 = 0x000F, 0x03F0, 0x7C00


def colt_table(an, free_top=True):
    """colt[t, o] = the row vector (bit M = bit t of X0[M, o]) over the rows used, as 2 uint64 words."""
    X = np.asarray(an.X0, np.int64)
    nr, nb = X.shape
    rows = nr - 1 if free_top else nr
    out = np.zeros((an.nl, nb, 2), np.uint64)
    for t in range(an.nl):
        for o in range(nb):
            lo = 0
            hi = 0
            for M in range(rows):
                if (int(X[M, o]) >> t) & 1:
                    if M < 64:
                        lo |= 1 << M
                    else:
                        hi |= 1 << (M - 64)
            out[t, o, 0] = lo
            out[t, o, 1] = hi
    return out, rows


@njit(cache=True)
def block_pool(colt, supp, d0, d1, d2, rows):
    """(cells, rank) of the block with duals d0, d1, d2."""
    nl, nb, _ = colt.shape
    lo = np.zeros(3 * nb, np.uint64)
    hi = np.zeros(3 * nb, np.uint64)
    n = 0
    for d in (d0, d1, d2):
        for o in range(nb):
            if d & supp[o]:
                a = np.uint64(0)
                b = np.uint64(0)
                for t in range(nl):
                    if (d >> t) & 1:
                        a ^= colt[t, o, 0]
                        b ^= colt[t, o, 1]
                lo[n] = a
                hi[n] = b
                n += 1
    # rank of n vectors of `rows` bits
    r = 0
    for col in range(rows):
        if col < 64:
            m = np.uint64(1) << np.uint64(col)
        else:
            m = np.uint64(1) << np.uint64(col - 64)
        sel = -1
        for i in range(r, n):
            w = lo[i] if col < 64 else hi[i]
            if w & m:
                sel = i
                break
        if sel < 0:
            continue
        a, b = lo[sel], hi[sel]
        lo[sel], hi[sel] = lo[r], hi[r]
        lo[r], hi[r] = a, b
        for i in range(r + 1, n):
            w = lo[i] if col < 64 else hi[i]
            if w & m:
                lo[i] ^= a
                hi[i] ^= b
        r += 1
        if r == n:
            break
    return n, r


@njit(cache=True)
def scan_slot(colt, supp, cand, others0, others1, rows):
    """Pool of the block (c, others0, others1) for every candidate c."""
    out = np.zeros(len(cand), np.int64)
    for i in range(len(cand)):
        n, r = block_pool(colt, supp, cand[i], others0, others1, rows)
        out[i] = n - r
    return out


def cls_of(d):
    if d & M1:
        return 1
    if d & M2:
        return 2
    return 3


def valid(duals):
    """Classes have the right counts and the diagonal parts are independent."""
    c1 = [d for d in duals if cls_of(d) == 1]
    c2 = [d for d in duals if cls_of(d) == 2]
    c3 = [d for d in duals if cls_of(d) == 3]
    return (len(c1) == 4 and len(c2) == 6 and len(c3) == 5 and rank_of([d & M1 for d in c1]) == 4
            and rank_of([d & M2 for d in c2]) == 6 and rank_of([d & M3 for d in c3]) == 5)


ALL = np.arange(1, 1 << 15, dtype=np.int64)
CLASS_CANDS = {1: ALL[(ALL & M1) != 0], 2: ALL[((ALL & M1) == 0) & ((ALL & M2) != 0)],
               3: ALL[((ALL & (M1 | M2)) == 0) & ((ALL & M3) != 0)]}


def random_start(pattern, rng):
    while True:
        duals = []
        for blk in pattern:
            for c in blk:
                pool = CLASS_CANDS[c]
                duals.append(int(pool[rng.randrange(len(pool))]))
        if valid(duals):
            return duals


def total(colt, supp, duals, rows):
    pools = []
    for b in range(5):
        n, r = block_pool(colt, supp, duals[3 * b], duals[3 * b + 1], duals[3 * b + 2], rows)
        pools.append(n - r)
    return pools


def exhaustive_133(colt, supp, rows):
    """Every block (d1 of class 1, plane of two class-3 duals) whose pool exceeds the generic one: all 155 planes
    of the class-3 space times all class-1 vectors. Returns [(pool, d1, d2, d3)], best first."""
    c3 = [int(x) for x in CLASS_CANDS[3]]
    planes = {}
    for i in c3:
        for j in c3:
            if i != j:
                key = tuple(sorted({i, j, i ^ j}))
                planes.setdefault(key, (i, j))
    c1 = CLASS_CANDS[1]
    out = []
    for d2, d3 in planes.values():
        sc = scan_slot(colt, supp, c1, d2, d3, rows)
        base = int(np.median(sc))
        for h in np.flatnonzero(sc > base):
            out.append((int(sc[h]), int(c1[h]), d2, d3))
    out.sort(reverse=True)
    return out


def start_with(pattern, fixed, rng, tries=200000):
    """A valid start whose blocks listed in `fixed` ({block: (d0, d1, d2)}) are given."""
    for _ in range(tries):
        duals = []
        for b, blk in enumerate(pattern):
            if b in fixed:
                duals += list(fixed[b])
            else:
                for c in blk:
                    pool = CLASS_CANDS[c]
                    duals.append(int(pool[rng.randrange(len(pool))]))
        if valid(duals):
            return duals
    return None


def climb(colt, supp, pattern, rows, rng, rounds, start=None, sub=40000, frozen=()):
    duals = list(start) if start else random_start(pattern, rng)
    pools = total(colt, supp, duals, rows)
    movable = [b for b in range(5) if b not in frozen]
    for it in range(rounds):
        b = movable[rng.randrange(len(movable))]
        j = rng.randrange(3)
        u = 3 * b + j
        c = pattern[b][j]
        cand = CLASS_CANDS[c]
        if len(cand) > sub:
            cand = cand[np.array(rng.sample(range(len(cand)), sub))]
        o = [duals[3 * b + q] for q in range(3) if q != j]
        sc = scan_slot(colt, supp, cand, o[0], o[1], rows)
        order = np.argsort(-sc, kind="stable")
        for i in order[:64]:
            if sc[i] < pools[b]:
                break
            trial = list(duals)
            trial[u] = int(cand[i])
            if valid(trial):
                if sc[i] > pools[b] or rng.random() < 0.5:
                    duals = trial
                    pools[b] = int(sc[i])
                break
    return duals, pools


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("k", type=int)
    ap.add_argument("--patterns", default="123,123,123,123,223;133,133,122,123,222;133,122,122,123,233")
    ap.add_argument("--rounds", type=int, default=60)
    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--sub", type=int, default=40000, help="candidates scanned per slot (all if larger)")
    ap.add_argument("--cfind2", action="store_true", help="candidate lists per pattern + compatible selection")
    ap.add_argument("--seed133", action="store_true", help="fix two deficient (1,3,3) blocks from an exhaustive scan")
    ap.add_argument("--fix222", default="", help="chem_ld2_c222 output: fix its best (2,2,2) blocks (+ seed133 pairs)")
    ap.add_argument("--n222", type=int, default=6, help="--fix222: how many of the best (2,2,2) blocks to try")
    args = ap.parse_args()
    words, _ = cc.load_words(args.table)
    an, lb = cc.analyze_first(words, args.k, cc.feasible_lbs(words), False)
    bd = core.BlockDecomp(an, 0, 1)
    bd.prepare()
    supp = np.array([int(s) for s in bd.supp], np.int64)
    assert set(int(s) for s in supp) <= {M1, M1 | M2, 0x7FFF}, "unexpected supports"
    t0 = time.time()
    colt, rows = colt_table(an, True)
    ident = [1 << t for t in (0, 6, 11, 1, 7, 12, 2, 8, 13, 3, 9, 14, 4, 5, 10)]
    print(f"k={args.k}: identity blocks pools {total(colt, supp, ident, rows)} ({time.time() - t0:.1f} s)", flush=True)
    rng = random.Random(args.seed)
    best = None
    if args.cfind2:
        extra = None
        if args.fix222:                    # with --cfind2: exhaustive (2,2,2) blocks join the candidate list
            extra = [tuple(x) for x in json.load(open(args.fix222))["blocks"]]
        res = cfind2(colt, supp, rows, rng, log=lambda m: print(m, flush=True), extra222=extra)
        if res is not None:
            duals, pools = res
            pat = "133,133,122,123,222"
            pattern = [tuple(int(ch) for ch in blk) for blk in pat.split(",")]
            if args.rounds:
                duals, pools = climb(colt, supp, pattern, rows, rng, args.rounds, start=duals, sub=args.sub)
                print(f"cfind2 + climb: pools {pools} total {sum(pools)} ({time.time() - t0:.1f} s)", flush=True)
            best = (duals, pools, pat)
    if args.fix222:
        # builder C2 (third holder): best exhaustive (2,2,2) blocks fixed as block 4, two compatible deficient (1,3,3)
        # blocks as blocks 0-1, blocks 2-3 ((1,2,2), (1,2,3)) climbed
        b222 = [tuple(x) for x in json.load(open(args.fix222))["blocks"]]
        ex = exhaustive_133(colt, supp, rows)
        print(f"exhaustive (1,3,3): {len(ex)} blocks above the median pool; best {[e[0] for e in ex[:10]]} "
              f"({time.time() - t0:.1f} s)", flush=True)
        pat = "133,133,122,123,222"
        pattern = [tuple(int(ch) for ch in blk) for blk in pat.split(",")]
        used = set()
        n_try = 0
        for blk in b222:
            key = tuple(sorted(blk[1:]))
            if key in used:
                continue
            used.add(key)
            n_try += 1
            if n_try > args.n222:
                break
            m2 = [d & M2 for d in blk[1:]]
            pairs = []
            for i in range(min(len(ex), 300)):
                for j in range(i + 1, min(len(ex), 300)):
                    a, b_ = ex[i], ex[j]
                    if rank_of([a[2], a[3], b_[2], b_[3]]) == 4 and rank_of([a[1] & M1, b_[1] & M1]) == 2:
                        pairs.append((a[0] + b_[0], a, b_))
            pairs.sort(reverse=True)
            for tot2, a, b_ in pairs[:max(1, args.restarts)]:
                st = start_with(pattern, {0: a[1:], 1: b_[1:], 4: tuple(blk[1:])}, rng)
                if st is None:
                    continue
                duals, pools = climb(colt, supp, pattern, rows, rng, args.rounds, start=st, sub=args.sub,
                                     frozen=(0, 1, 4))
                print(f"fix222 {blk[0]} {[hex(d) for d in blk[1:]]} + pair {a[0]}+{b_[0]}: pools {pools} total "
                      f"{sum(pools)} ({time.time() - t0:.1f} s)", flush=True)
                if best is None or sum(pools) > sum(best[1]):
                    best = (duals, pools, pat)
            _ = m2
    if args.seed133 and best is None:
        # two deficient (1,3,3) blocks from the exhaustive scan, then the other three blocks climbed
        ex = exhaustive_133(colt, supp, rows)
        print(f"exhaustive (1,3,3): {len(ex)} blocks above the median pool; best {[e[0] for e in ex[:10]]} "
              f"({time.time() - t0:.1f} s)", flush=True)
        pat = "133,133,122,123,222"
        pattern = [tuple(int(ch) for ch in blk) for blk in pat.split(",")]
        pairs = []
        for i in range(min(len(ex), 300)):
            for j in range(i + 1, min(len(ex), 300)):
                a, b_ = ex[i], ex[j]
                if rank_of([a[2], a[3], b_[2], b_[3]]) == 4 and rank_of([a[1] & M1, b_[1] & M1]) == 2:
                    pairs.append((a[0] + b_[0], a, b_))
        pairs.sort(reverse=True)
        for tot2, a, b_ in pairs[:args.restarts]:
            st = start_with(pattern, {0: a[1:], 1: b_[1:]}, rng)
            if st is None:
                continue
            duals, pools = climb(colt, supp, pattern, rows, rng, args.rounds, start=st, sub=args.sub, frozen=(0, 1))
            print(f"seed133 pair {a[0]}+{b_[0]}: pools {pools} total {sum(pools)} ({time.time() - t0:.1f} s)",
                  flush=True)
            if best is None or sum(pools) > sum(best[1]):
                best = (duals, pools, pat)
    for pat in ([] if (args.seed133 or args.cfind2 or args.fix222) else args.patterns.split(";")):
        pattern = [tuple(int(ch) for ch in blk) for blk in pat.split(",")]
        for rs in range(args.restarts):
            duals, pools = climb(colt, supp, pattern, rows, rng, args.rounds, sub=args.sub)
            print(f"pattern {pat} restart {rs}: pools {pools} total {sum(pools)} ({time.time() - t0:.1f} s)",
                  flush=True)
            if best is None or sum(pools) > sum(best[1]):
                best = (duals, pools, pat)
    duals, pools, pat = best
    assert valid(duals)
    Qrows = [sum(((duals[u] >> t) & 1) << u for u in range(15)) for t in range(15)]
    P = lb2.mat_inv(Qrows, 15)
    print("best", pat, pools, sum(pools), flush=True)
    out = dict(k=args.k, pattern=pat, pools=[int(p) for p in pools], total=int(sum(pools)),
               dual=[hex(d) for d in duals], P=[hex(v) for v in P], seed=args.seed, rounds=args.rounds,
               words_sha256=cc.words_sha(words))
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out, fh)



# ---- candidate lists per block pattern, then a compatible selection (cfind2) --------------------------------------
def sample_deficient(colt, supp, rows, pattern, n_pairs, rng):
    """Deficient blocks of one pattern (c1, c2, c3) (slot 0 scanned over all its class candidates, slots 1-2 drawn
    at random n_pairs times). Returns [(pool, d0, d1, d2)] with pool > the median of each scan."""
    c0, c1, c2 = pattern
    cand = CLASS_CANDS[c0]
    out = []
    for _ in range(n_pairs):
        a = int(CLASS_CANDS[c1][rng.randrange(len(CLASS_CANDS[c1]))])
        b = int(CLASS_CANDS[c2][rng.randrange(len(CLASS_CANDS[c2]))])
        if a == b or rank_of([a, b]) < 2:
            continue
        sc = scan_slot(colt, supp, cand, a, b, rows)
        base = int(np.median(sc))
        for h in np.flatnonzero(sc > base):
            d0 = int(cand[h])
            if rank_of([d0, a, b]) == 3:
                out.append((int(sc[h]), d0, a, b))
    out.sort(reverse=True)
    return out


def select_set(lists, pattern_order, tries=200000, rng=None):
    """One block from each list (lists[i] for pattern_order[i]) with valid() duals and the largest total pool.
    Greedy over the best few of each list with random restarts."""
    best = None
    tops = [lst[:60] for lst in lists]
    if any(not t for t in tops):
        return None
    for _ in range(tries):
        pick = [t[min(len(t) - 1, int(rng.expovariate(0.15)))] for t in tops]
        duals = [v for blk in pick for v in blk[1:]]
        if len(set(duals)) < 15 or not valid(duals):
            continue
        tot = sum(blk[0] for blk in pick)
        if best is None or tot > best[0]:
            best = (tot, pick)
    return best


def generic_222(colt, supp, rows, rng, n=200):
    """(2,2,2) blocks at whatever pool they have (used when no deficient one is known): random class-2 triples
    with independent class-2 parts. Returns [(pool, d0, d1, d2)], best first. (builder C2, third holder)"""
    c2 = CLASS_CANDS[2]
    out = []
    while len(out) < n:
        ds = [int(c2[rng.randrange(len(c2))]) for _ in range(3)]
        if rank_of([d & M2 for d in ds]) < 3:
            continue
        cells, r = block_pool(colt, supp, ds[0], ds[1], ds[2], rows)
        out.append((int(cells - r), ds[0], ds[1], ds[2]))
    out.sort(reverse=True)
    return out


def cfind2(colt, supp, rows, rng, budget_pairs=(120, 200, 60), log=print, extra222=None):
    """Pattern 133,133,122,123,222: exhaustive (1,3,3); sampled (1,2,2), (1,2,3), (2,2,2); compatible selection.
    extra222: (2,2,2) blocks from chem_ld2_c222 (exhaustive), merged into the sampled list; if the list stays
    empty, generic (2,2,2) blocks are used so that the deficient blocks of the other patterns still combine."""
    import time as _t
    t0 = _t.time()
    l133 = exhaustive_133(colt, supp, rows)
    l123 = sample_deficient(colt, supp, rows, (1, 2, 3), budget_pairs[0], rng)
    l122 = sample_deficient(colt, supp, rows, (1, 2, 2), budget_pairs[1], rng)
    l222 = sample_deficient(colt, supp, rows, (2, 2, 2), budget_pairs[2], rng)
    if extra222:
        l222 = sorted(set(l222) | set(tuple(x) for x in extra222), reverse=True)
    if not l222:
        l222 = generic_222(colt, supp, rows, rng)
        log(f"cfind2: no deficient (2,2,2) block; {len(l222)} generic ones (pool {l222[0][0]})")
    log(f"cfind2 lists: 133 {len(l133)} (best {[x[0] for x in l133[:4]]}), 123 {len(l123)} "
        f"({[x[0] for x in l123[:4]]}), 122 {len(l122)} ({[x[0] for x in l122[:4]]}), 222 {len(l222)} "
        f"({[x[0] for x in l222[:4]]}) ({_t.time() - t0:.0f} s)")
    sel = select_set([l133, l133, l122, l123, l222], "133,133,122,123,222", rng=rng)
    if sel is None:
        return None
    tot, pick = sel
    duals = [v for blk in pick for v in blk[1:]]
    pools = total(colt, supp, duals, rows)
    log(f"cfind2 selected: pools {pools} total {sum(pools)} ({_t.time() - t0:.0f} s)")
    return duals, pools


if __name__ == "__main__":
    main()
