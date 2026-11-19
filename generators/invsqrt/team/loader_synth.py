#!/usr/bin/env python3
# New file added by the team (Maxime Bendavid), 2026-11-18, for idea INV-A02 (plan step S2.8).
"""XOR-AND synthesis of a table loader: few AND gates for many outputs of a 6-bit address.

Problem: given T output bits (truth tables over the 2^n addresses, as integers), find a program that computes
them using only XOR/NOT (free) and as few two-input ANDs of affine functions as possible.

Method (bilinear split): split the n address bits into a high group H and a low group L (|H| + |L| = n). All
monomials of each group (2^|H| and 2^|L| of them, the constant and the variables included) cost
2^|H| - |H| - 1 + 2^|L| - |L| - 1 ANDs (each monomial of degree >= 2 is one AND from smaller ones). Every output
is then a bilinear form f_t(h, l) = m_H(h)^T M_t m_L(l) over GF(2). Terms with the constant monomial on either
side are free (XOR of already computed monomials). The rest needs products r_k = (u_k . m_H)(w_k . m_L) (one AND
each, the two factors being XORs of monomials): find few rank-1 matrices u_k w_k^T whose span, together with the
free terms, contains every M_t. Lower bound: the GF(2) rank of the M_t modulo the free terms. Upper bound: a
randomized greedy that always keeps the remaining target rank as low as possible, with restarts.

    python loader_synth.py MODEL.json [--restarts 200]

Output: per split, monomial cost + products (= total ANDs), the lower bound, and the best program as JSON.
Complexity: per greedy step, all (2^a - 1)(2^b - 1) rank-1 candidates are scored with a GF(2) rank update
(bit-packed integers), i.e. O(candidates x T) word operations per step.
"""
from __future__ import annotations

import argparse
import itertools
import json
import random
import sys
import time


def gf2_reduce(basis: dict, v: int) -> int:
    """Reduce v by a GF(2) basis stored as {pivot_bit: vector}."""
    while v:
        p = v.bit_length() - 1
        b = basis.get(p)
        if b is None:
            return v
        v ^= b
    return 0


def gf2_insert(basis: dict, v: int) -> bool:
    v = gf2_reduce(basis, v)
    if v:
        basis[v.bit_length() - 1] = v
        return True
    return False


def monomials(bits: tuple) -> list:
    """All subsets of `bits` (as tuples), constant first, by degree."""
    out = []
    for d in range(len(bits) + 1):
        out += list(itertools.combinations(bits, d))
    return out


def matrices(words: list, out_bits: list, n: int, H: tuple, L: tuple):
    """Coefficient matrices M_t (as flat bitmasks over (i, j) with i in H-monomials, j in L-monomials)."""
    mh, ml = monomials(H), monomials(L)
    # ANF of each output via the Moebius transform over the n address bits.
    mats = []
    for (w, bit) in out_bits:
        f = [(words[a] >> bit) & 1 for a in range(1 << n)]
        anf = f[:]
        for i in range(n):
            for a in range(1 << n):
                if a >> i & 1:
                    anf[a] ^= anf[a ^ (1 << i)]
        m = 0
        for a in range(1 << n):
            if anf[a]:
                vars_ = tuple(i for i in range(n) if a >> i & 1)
                hi = tuple(i for i in vars_ if i in H)
                lo = tuple(i for i in vars_ if i in L)
                m |= 1 << (mh.index(hi) * len(ml) + ml.index(lo))
        mats.append(m)
    return mats, mh, ml


def free_space(nh: int, nl: int) -> dict:
    """Matrices supported on the constant row or the constant column (free XORs of monomials)."""
    basis = {}
    for j in range(nl):
        gf2_insert(basis, 1 << (0 * nl + j))
    for i in range(nh):
        gf2_insert(basis, 1 << (i * nl + 0))
    return basis


def rank1(u: int, w: int, nh: int, nl: int) -> int:
    m = 0
    for i in range(nh):
        if u >> i & 1:
            m |= w << (i * nl)
    return m


def residual_rank(basis: dict, targets: list) -> int:
    b = dict(basis)
    r = 0
    for t in targets:
        if gf2_insert(b, t):
            r += 1
    return r


_CANDS = {}


def candidates(nh, nl):
    """All rank-1 products whose two factors both involve a non-constant monomial, as (u, w, matrix)."""
    key = (nh, nl)
    if key not in _CANDS:
        _CANDS[key] = [(u, w, rank1(u, w, nh, nl)) for u in range(1, 1 << nh) for w in range(1, 1 << nl)
                       if (u & ~1) and (w & ~1)]
    return _CANDS[key]


def greedy(mats, nh, nl, rng, max_products=64):
    """Randomized greedy. Adding a product m lowers the remaining target rank by one exactly when m, reduced by
    the current basis, lies in the span of the reduced targets; that test is one reduction per candidate.
    Ties are broken by how many targets become directly covered, then at random."""
    basis = free_space(nh, nl)
    cands = candidates(nh, nl)
    chosen = []
    while True:
        tbasis = dict(basis)
        res = [t for t in (gf2_reduce(basis, t) for t in mats) if t]
        if not res:
            return chosen
        if len(chosen) >= max_products:
            return None
        for t in res:
            gf2_insert(tbasis, t)
        useful = [(u, w, m) for (u, w, m) in cands if gf2_reduce(basis, m) and not gf2_reduce(tbasis, m)]
        if not useful:
            return None
        # Tie-break: prefer products that make some target land in the span immediately.
        scored = []
        for (u, w, m) in useful:
            b = dict(basis)
            gf2_insert(b, m)
            covered = sum(1 for t in res if not gf2_reduce(b, t))
            scored.append((covered, u, w, m))
        top = max(s[0] for s in scored)
        _, u, w, m = rng.choice([s for s in scored if s[0] == top])
        gf2_insert(basis, m)
        chosen.append((u, w))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--a-bits", type=int, default=13)
    ap.add_argument("--s-bits", type=int, default=7)
    ap.add_argument("--json")
    a = ap.parse_args()
    model = json.load(open(a.model))
    A, S = model["A"], model["S"]
    n = 6
    assert len(A) == len(S) == 1 << n
    words = [A[i] | (S[i] << a.a_bits) for i in range(1 << n)]
    out_bits = [(0, b) for b in range(a.a_bits + a.s_bits)]
    results = []
    t0 = time.perf_counter()
    for hsize in (3, 2):
        for H in itertools.combinations(range(n), hsize):
            L = tuple(i for i in range(n) if i not in H)
            mats, mh, ml = matrices(words, out_bits, n, H, L)
            nh, nl = len(mh), len(ml)
            mono_cost = (nh - len(H) - 1) + (nl - len(L) - 1)
            lb = residual_rank(free_space(nh, nl), mats)
            rng = random.Random(1)
            best = None
            for _ in range(a.restarts):
                c = greedy(mats, nh, nl, rng)
                if c is not None and (best is None or len(c) < len(best)):
                    best = c
            results.append(dict(H=H, L=L, monomial_ANDs=mono_cost, lower_bound_products=lb,
                                products=len(best) if best else None,
                                total_ANDs=mono_cost + len(best) if best else None,
                                program=[[u, w] for u, w in best] if best else None))
            print(f"H={H} L={L}: monomials {mono_cost} + products {results[-1]['products']} (lower bound {lb}) "
                  f"= {results[-1]['total_ANDs']}  [{time.perf_counter() - t0:.0f} s]", flush=True)
    best = min((r for r in results if r["total_ANDs"]), key=lambda r: r["total_ANDs"])
    print("best:", {k: v for k, v in best.items() if k != "program"})
    if a.json:
        json.dump(dict(model=a.model, results=results, best=best), open(a.json, "w"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
