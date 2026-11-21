#!/usr/bin/env python3
"""invsqrt_pp_smallfn.py -- workstream invsqrt_pp (P1): rebuild a dropped value that is an exact function of a few
live wires with as few ANDs as possible (patterns that never occur are don't-cares).

support_tt(S, w, cands): greedy smallest support (<= 6 wires) on sample lanes, then the exact truth table over ALL
  inputs with its care mask (patterns that occur).
synth(tt, care, n): the cheapest program found, as ('aff', L) | ('one', L, A, B) | ('two', L, A, B, Cf, D) |
  ('nest', L, L1, A, B, Cf): forms are (mask over the n support variables, constant):
    aff   f = L
    one   f = L + A B
    two   f = L + A B + C D
    nest  f = L + (L1 + A B) C            (needs one scratch wire at the rebuild point)
emit(prog, sup, w, scratch): the gate list that writes f into w (w zero before), every support wire restored.
All tests on truth tables of the n <= 6 support variables (64 patterns); exactness on all inputs follows from the
exact truth table (a value that is a function of the support on every input).
"""
import itertools

import numpy as np

import invsqrt_pp_sim as sim

X_, CX, CCX, HMR, CZ = 0x11, 0x21, 0x31, 0x50, 0x20
LANES = np.sort(np.random.default_rng(11).choice(sim.NX, size=4096, replace=False)).astype(np.int64)
LW, LB = LANES // 64, (LANES % 64).astype(np.uint64)


def _bits(v):
    return np.unpackbits(v.view(np.uint8), bitorder='little')[:sim.NX]


def support_tt(S, w, cands, kmax=6):
    vs = ((S[w][LW] >> LB) & np.uint64(1)).astype(np.int64)
    X = ((S[cands][:, LW] >> LB[None, :]) & np.uint64(1)).astype(np.int64)
    sup = []
    bc = None
    for _ in range(kmax):
        best, bc = None, None
        for t in range(len(cands)):
            if t in sup:
                continue
            key = np.zeros(len(vs), np.int64)
            for u in sup + [t]:
                key = key * 2 + X[u]
            c = len(set(key[vs == 0].tolist()) & set(key[vs == 1].tolist()))
            if bc is None or c < bc:
                best, bc = t, c
        sup.append(best)
        if bc == 0:
            break
    if bc != 0:
        return None
    U = [cands[t] for t in sup]
    key = np.zeros(sim.NX, np.int64)
    for j, u in enumerate(U):
        key |= _bits(S[u]).astype(np.int64) << j
    v = _bits(S[w]).astype(np.int64)
    n = len(U)
    ones = np.bincount(key, weights=v, minlength=1 << n)
    tot = np.bincount(key, minlength=1 << n)
    if np.any((ones > 0) & (ones < tot)):
        return None                                  # not a function of the support on all inputs
    tt = (ones > 0).astype(np.int64)
    care = (tot > 0).astype(np.int64)
    return U, tt, care


def _forms(n):
    """truth tables (as int arrays of 0/1 over 2^n patterns) of all affine forms, keyed by (mask, const)."""
    pats = np.arange(1 << n)
    out = {}
    for mask in range(1 << n):
        lin = np.zeros(1 << n, np.int64)
        for j in range(n):
            if mask >> j & 1:
                lin ^= (pats >> j) & 1
        out[(mask, 0)] = lin
        out[(mask, 1)] = lin ^ 1
    return out


def _pack(arr):
    return int(sum(int(b) << i for i, b in enumerate(arr.tolist())))


def synth(tt, care, n, allow_nest=True):
    F = _forms(n)
    keys = list(F.keys())
    cm = _pack(care)
    f = _pack(tt) & cm
    packed = {k: _pack(F[k]) & cm for k in keys}
    by_val = {}
    for k in keys:
        by_val.setdefault(packed[k], k)
    if f in by_val:
        return ('aff', by_val[f])
    nonconst = [k for k in keys if k[0] != 0]
    prods = {}
    for a, b in itertools.combinations(nonconst, 2):
        if a[0] == b[0]:
            continue
        p = packed[a] & packed[b]
        prods.setdefault(p, (a, b))
    for p, (a, b) in prods.items():
        r = f ^ p
        if r in by_val:
            return ('one', by_val[r], a, b)
    # two independent products: f + L = P1 + P2  ->  for each P1, is f + P1 + P2 affine for some P2?
    aff_vals = set(by_val)
    plist = list(prods.items())
    for p1, (a, b) in plist:
        r = f ^ p1
        for av in aff_vals:
            p2 = r ^ av
            if p2 in prods:
                c, d = prods[p2]
                return ('two', by_val[av], a, b, c, d)
    if allow_nest and n <= 4:
        for p1, (a, b) in plist:
            for l1 in keys:
                t1 = packed[l1] ^ p1
                for cf in nonconst:
                    r = f ^ (t1 & packed[cf])
                    if r in aff_vals:
                        return ('nest', by_val[r], l1, a, b, cf)
    return None


def _expose(seq, sup, cur, form, avoid=None):
    """put affine form (mask, const) over the ORIGINAL support values on one support wire; cur[j] = mask of original
    variables the wire j now holds (plus const in cur_c). Returns the pivot j."""
    mask, const = form
    # express form in the current wire basis by elimination: find a combination of current wires equal to mask
    n = len(sup)
    rows = [(cur[j][0], 1 << j) for j in range(n)]
    piv = []
    for v, tg in rows:
        for (pb, pv, pt) in piv:
            if v >> pb & 1:
                v ^= pv
                tg ^= pt
        if v:
            pb = v.bit_length() - 1
            piv.append((pb, v, tg))
    t, comb = mask, 0
    for (pb, pv, pt) in piv:
        if t >> pb & 1:
            t ^= pv
            comb ^= pt
    assert t == 0
    js = [j for j in range(n) if comb >> j & 1]
    js.sort(key=lambda j: j == avoid)
    p = js[0]
    assert p != avoid
    for j in js[1:]:
        seq.append((CX, (sup[j], sup[p]), 0, 0))
        cur[p] = (cur[p][0] ^ cur[j][0], cur[p][1] ^ cur[j][1])
    # now wire p holds mask with constant cur[p][1]; flip to the wanted constant
    if cur[p][1] != const:
        seq.append((X_, (sup[p],), 0, 0))
        cur[p] = (cur[p][0], cur[p][1] ^ 1)
    return p


def emit(prog, sup, w, scratch=None):
    """gates writing f into zero wire w; support wires restored. Returns (gates, number of CCX)."""
    n = len(sup)
    seq = []
    kind = prog[0]

    def add_affine(form, target):
        mask, const = form
        for j in range(n):
            if mask >> j & 1:
                seq.append((CX, (sup[j], target), 0, 0))
        if const:
            seq.append((X_, (target,), 0, 0))

    def product_into(a, b, target):
        cur = [(1 << j, 0) for j in range(n)]
        mark = len(seq)
        pa = _expose(seq, sup, cur, a)
        save_a = list(seq[mark:])
        mark2 = len(seq)
        pb = _expose(seq, sup, cur, b, avoid=pa)
        assert pb != pa
        save_b = list(seq[mark2:])
        seq.append((CCX, (sup[pa], sup[pb], target), 0, 0))
        for g in reversed(save_b):
            seq.append(g)
        for g in reversed(save_a):
            seq.append(g)
    if kind == 'aff':
        add_affine(prog[1], w)
        return seq, 0
    if kind == 'one':
        add_affine(prog[1], w)
        product_into(prog[2], prog[3], w)
        return seq, 1
    if kind == 'two':
        add_affine(prog[1], w)
        product_into(prog[2], prog[3], w)
        product_into(prog[4], prog[5], w)
        return seq, 2
    if kind == 'nest':
        assert scratch is not None
        L, L1, a, b, cf = prog[1:]
        # s = L1 + a b ; w += s * cf ; erase s by HMR + phase repair (CZ on the exposed a, b; Z on L1's wires)
        add_affine(L1, scratch)
        product_into(a, b, scratch)
        cur = [(1 << j, 0) for j in range(n)]
        mark = len(seq)
        pc = _expose(seq, sup, cur, cf)
        save_c = list(seq[mark:])
        seq.append((CCX, (scratch, sup[pc], w), 0, 0))
        for g in reversed(save_c):
            seq.append(g)
        # erase the scratch: recompute it away (one more CCX) -- kept simple and exact
        product_into(a, b, scratch)
        add_affine(L1, scratch)
        add_affine(L, w)
        return seq, 3
    raise ValueError(kind)
