"""chem_ld2_e21dual: BASIS SEARCH FOR THE FOUR-BLOCK STAGE A of the kept-product loader (twenty-first holder C2).

R58-CHEM-KEPT's stage A (tools/r58_chem_kept_c2a.py) plans the table on the 13 low points outside D = {a, b, c}
in four blocks of three new low coordinates; block `rep` of the basis is replaced by the indicators of D.
A block's relation pool (its codim) = cells - rank of the rows-M x cells-(u, o) matrix, and it depends only on
the block's DUAL subspace Phi_b (span of the dual functionals of its three coordinates = the columns of
Q = P^-1). Every stage-A dual must annihilate the three point indicators of D, so the four Phi_b are a direct sum
of a 12-dimensional space of duals. Generic pool 15 (rank 75: the V-padded table's columns lie in a hyperplane);
a pool of 16 is one free cell less (two such blocks = one product off stage A if the pairing stays perfect).

This tool (1) checks its pool against the engine on a given pfile (--check), (2) samples random 3-dim dual
subspaces of the 12-dim space (numba), keeps those with pool >= --minpool, (3) picks four that form a direct sum
with as many deficient blocks as possible, (4) writes pfiles {P, blocks} with block 4 = the D block (use rep 4).

  CHEM_BROT=7 CHEM_HB=77 python -I chem_ld2_e21dual.py PAD.npz --k 2 --D 13,14,15 --n 2000000 --seed 0 --out DIR
        [--check PFILE --rep R] [--npf 4]
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
import chem_ld2_basis as lb2  # noqa: E402
from chem_ld2_cfind import colt_table, block_pool  # noqa: E402

bits = core.bits


def rank_of(vs):
    sp = core.Span()
    return sum(1 for v in vs if v and sp.add(v))


def point_coords(an, D):
    """mD[i] = the indicator of low point D[i] in the analysis' low basis coordinates (as c2a computes it)."""
    LB = []
    for bs in an.basis:
        tt = 0
        for lpt in range(16):
            v = 0
            for m in bits(bs["anf"]):
                if (lpt & m) == m:
                    v ^= 1
            tt |= v << lpt
        LB.append(tt)
    sp = core.Span()
    for t in range(an.nl):
        sp.add(LB[t], 1 << t)
    out = []
    for p in D:
        r, c = sp.reduce(1 << p)
        assert r == 0
        out.append(c)
    return out


def dual_space(mD, n=15):
    """Basis of the duals d (masks over the n old coordinates) with parity(d & m) = 0 for every m in mD."""
    rows = list(mD)
    # nullspace of the matrix with rows mD (as functionals d -> parity(d & m))
    return core.nullspace(rows, n)


def pools_of_pfile(colt, supp, rows, P2, blocks):
    """Pools of every block of basis P2 (rows = new vectors in old coordinates)."""
    n = len(P2)
    Q = lb2.mat_inv(P2, n)
    duals = [sum(((Q[t] >> u) & 1) << t for t in range(n)) for u in range(n)]
    out = []
    for b in blocks:
        cn, r = block_pool(colt, supp, np.int64(duals[b[0]]), np.int64(duals[b[1]]), np.int64(duals[b[2]]), rows)
        out.append(int(cn - r))
    return out, duals


@njit(cache=True)
def _sample(colt, supp, dsub, rows, n, seed, minpool, out):
    np.random.seed(seed)
    cnt = 0
    nd = dsub.shape[0]
    full = (1 << nd) - 1
    for it in range(n):
        c0 = np.random.randint(1, full + 1)
        c1 = np.random.randint(1, full + 1)
        c2 = np.random.randint(1, full + 1)
        if c1 == c0 or c2 == c0 or c2 == c1 or c2 == (c0 ^ c1):
            continue
        d0 = 0
        d1 = 0
        d2 = 0
        for i in range(nd):
            if (c0 >> i) & 1:
                d0 ^= dsub[i]
            if (c1 >> i) & 1:
                d1 ^= dsub[i]
            if (c2 >> i) & 1:
                d2 ^= dsub[i]
        cn, r = block_pool(colt, supp, d0, d1, d2, rows)
        if cn - r >= minpool and cnt < out.shape[0]:
            out[cnt, 0] = cn - r
            out[cnt, 1] = d0
            out[cnt, 2] = d1
            out[cnt, 3] = d2
            cnt += 1
    return cnt


def canon(ds):
    """Canonical key of the span of the given masks (reduced row echelon form, distinct leading bits)."""
    basis = []
    for v in ds:
        for b in sorted(basis, reverse=True):
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    basis.sort(reverse=True)
    for j in range(len(basis)):
        hb = 1 << (basis[j].bit_length() - 1)
        for i in range(len(basis)):
            if i != j and basis[i] & hb:
                basis[i] ^= basis[j]
    return tuple(sorted(basis, reverse=True))


def pick_sets(hits, need, rng, tries=20000, base_blocks=None, nset=4):
    """Up to `need` distinct sets of nset blocks in direct sum maximising the deficient count.
    hits: [(pool, (d0, d1, d2))]; generic blocks to complete a set are taken from base_blocks (list of triples)."""
    sets = {}
    hits = sorted(hits, key=lambda h: -h[0])
    for _ in range(tries):
        order = list(range(len(hits)))
        rng.shuffle(order)
        chosen = []
        vecs = []
        for i in order:
            ds = hits[i][1]
            if rank_of(vecs + list(ds)) == len(vecs) + 3:
                chosen.append(i)
                vecs += list(ds)
                if len(chosen) == nset:
                    break
        if len(chosen) < nset and base_blocks:
            for blk in base_blocks:
                if len(chosen) >= nset:
                    break
                if rank_of(vecs + list(blk)) == len(vecs) + 3:
                    chosen.append(("g", blk))
                    vecs += list(blk)
        if len(vecs) < 3 * nset:
            continue
        key = tuple(sorted(str(c) for c in chosen))
        score = sum(hits[c][0] for c in chosen if not isinstance(c, tuple))
        sets[key] = (score, chosen)
        if len(sets) >= 4 * need:
            break
    best = sorted(sets.values(), key=lambda s: -s[0])[:need]
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pad")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--D", default="13,14,15")
    ap.add_argument("--n", type=int, default=2000000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minpool", type=int, default=16)
    ap.add_argument("--check", default="")
    ap.add_argument("--rep", type=int, default=-1)
    ap.add_argument("--npf", type=int, default=4)
    ap.add_argument("--out", default="")
    ap.add_argument("--pfx", default="", help="pfile name prefix (default k<K>d)")
    ap.add_argument("--nset", type=int, default=4, help="deficient blocks a set (4: kept stage A with D; 3: the "
                    "staircase's stage A, --D none: blocks 0-2 = stage A, blocks 3-4 = the complement)")
    a = ap.parse_args()
    t0 = time.time()
    core.J = 5
    words = np.asarray(np.load(a.pad)["angle_words"], np.int64)
    an = core.analyze(words, a.k, 4)
    nl = an.nl
    X = np.asarray(an.X0, np.int64)
    nr, nb = X.shape
    supp = np.array([int(np.bitwise_or.reduce(X[:, o])) for o in range(nb)], np.int64)
    colt, rows = colt_table(an, False)
    D = [int(x) for x in a.D.split(",")] if a.D != "none" else []
    mD = point_coords(an, D)
    dsub = dual_space(mD, nl) if D else [1 << t for t in range(nl)]
    print(f"k {a.k}: nr {nr} nb {nb} rows {rows} ymasks {[hex(y) for y in an.ymasks]} supp {sorted(set(hex(s) for s in supp))}"
          f" D {D} dual space dim {len(dsub)}", flush=True)
    if a.check:
        p0 = json.load(open(a.check))
        P = [int(x, 16) for x in p0["P"]]
        five = [list(b) for b in p0["blocks"]]
        pools5, _ = pools_of_pfile(colt, supp, rows, P, five)
        print(f"check {os.path.basename(a.check)}: five-block pools (rows {rows}) {pools5}", flush=True)
        if a.rep >= 0:
            P2 = list(P)
            for i, t in enumerate(five[a.rep]):
                P2[t] = mD[i]
            try:
                pools4, duals = pools_of_pfile(colt, supp, rows, P2, five)
                ok = all(((d & m).bit_count() % 2 == 0) for b in range(5) if b != a.rep for t in five[b]
                         for d in [duals[t]] for m in mD)
                print(f"check rep {a.rep}: pools {pools4} (stage-A duals annihilate D: {ok})", flush=True)
            except ValueError:
                print("check: singular", flush=True)
    if not a.out:
        return
    os.makedirs(a.out, exist_ok=True)
    dsa = np.array(dsub, np.int64)
    out = np.zeros((200000, 4), np.int64)
    t1 = time.time()
    cnt = _sample(colt, supp, dsa, rows, a.n, a.seed, a.minpool, out)
    dt = time.time() - t1
    hits = {}
    for i in range(cnt):
        ds = (int(out[i, 1]), int(out[i, 2]), int(out[i, 3]))
        key = canon(ds)
        if key not in hits or hits[key][0] < int(out[i, 0]):
            hits[key] = (int(out[i, 0]), ds)
    pc = {}
    for p, _ in hits.values():
        pc[p] = pc.get(p, 0) + 1
    print(f"sampled {a.n} subspaces in {dt:.0f} s: {cnt} hits, {len(hits)} distinct, pools {sorted(pc.items())}",
          flush=True)
    hl = list(hits.values())
    json.dump(dict(k=a.k, D=D, n=a.n, seed=a.seed, hits=[[p, [hex(d) for d in ds]] for p, ds in hl]),
              open(os.path.join(a.out, f"hits_k{a.k}_D{'-'.join(map(str, D))}_s{a.seed}.json"), "w"))
    # generic completion blocks: random triples of the dual space
    rng = random.Random(a.seed)
    gen = []
    for _ in range(400):
        cs = [0, 0, 0]
        while rank_of(cs) < 3:
            cs = [rng.randrange(1, 1 << len(dsub)) for _ in range(3)]
        gen.append(tuple(lb2.vec_mul(c, dsub) for c in cs))
    sets = pick_sets(hl, a.npf, rng, base_blocks=gen, nset=a.nset)
    for si, (score, chosen) in enumerate(sets):
        blocks4 = [hl[c][1] if not isinstance(c, tuple) else c[1] for c in chosen]
        duals12 = [d for blk in blocks4 for d in blk]
        # complete to a basis of all 15 (with D: any completion, the D block's primal span is fixed; without D:
        # random completions, the one whose two complement blocks have the largest pools)
        comp = []
        if D:
            for t in range(nl):
                if rank_of(duals12 + comp + [1 << t]) == len(duals12) + len(comp) + 1:
                    comp.append(1 << t)
                if len(comp) == nl - len(duals12):
                    break
        else:
            bestc = None
            for _ in range(200):
                cc = []
                while len(cc) < nl - len(duals12):
                    v = rng.randrange(1, 1 << nl)
                    if rank_of(duals12 + cc + [v]) == len(duals12) + len(cc) + 1:
                        cc.append(v)
                pc = sum(int(block_pool(colt, supp, np.int64(cc[3 * i]), np.int64(cc[3 * i + 1]),
                                        np.int64(cc[3 * i + 2]), rows)[0] - block_pool(
                    colt, supp, np.int64(cc[3 * i]), np.int64(cc[3 * i + 1]), np.int64(cc[3 * i + 2]), rows)[1])
                         for i in range(len(cc) // 3))
                if bestc is None or pc > bestc[0]:
                    bestc = (pc, cc)
            comp = bestc[1]
        duals = duals12 + comp
        Qrows = [sum(((duals[u] >> t) & 1) << u for u in range(nl)) for t in range(nl)]
        P = lb2.mat_inv(Qrows, nl)
        five = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 14]]
        P2 = list(P)
        if D:
            for i, t in enumerate(five[4]):
                P2[t] = mD[i]
        pools4, _ = pools_of_pfile(colt, supp, rows, P2, five)
        fn = os.path.join(a.out, f"{a.pfx or f'k{a.k}d'}{si}.json")
        json.dump(dict(k=a.k, D=D, rep=4, pools=pools4, P=[hex(v) for v in P], blocks=five,
                       dual=[hex(d) for d in duals], seed=a.seed, src="chem_ld2_e21dual"), open(fn, "w"))
        print(f"set {si}: pools (D block replaced) {pools4} -> {fn}", flush=True)
    print(f"done {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
