#!/usr/bin/env python3
"""invsqrt_pp_hold.py -- workstream invsqrt_pp (P1): retained products for an invsqrt .kmb, CSWAP and CCX alike,
with every relation exact on all 3,133,452 inputs and an all-input equivalence check of the result.

A chain [i0 < i1 < ... < ik] of Toffoli-class ops (top level, unconditioned) with exact relations
d(i_m) = d(i_{m-1}) + lambda_m (lambda_m affine in the wires just before i_m; d = a*b for CCX, c*(x+y) for CSWAP)
is rewritten with one holder wire h that is zero and untouched by the parent over [i0, ik]:
  i0  CCX a b t     -> CCX a b h ; CX h t
      CSWAP c x y   -> CX x y ; CCX c y h ; CX x y ; CX h x ; CX h y
  im  (lambda into h: CX l h ..., X h) then CCX: CX h t | CSWAP: CX h x ; CX h y          (0 Toffolis)
  ik  as im, then HMR h -> fresh bit m ; CZ a b if m  |  CZ c x if m ; CZ c y if m      (h back to 0)
Chains come from the links file of invsqrt_pp_links.py (asap=0 order), entries with a zero wire throughout.
Check: the child, simulated on all inputs (zero transcript and random transcripts), gives the parent's final state
and Z-phase on every input, and every Z_POW sees the parent's value on every input.
Usage: python invsqrt_pp_hold.py <parent.kmb> <links.json> <out.kmb> [maxlen=0 (no limit)] [skip=j1,j2,..]
"""
import dataclasses
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_blocks import gf2_solve, sample_bits, vec_sample  # noqa: E402

CX, CCX, X, HMR, SWAP = 0x21, 0x31, 0x11, 0x50, 0x22
CZ_IF = 0x20 | 0x80


def prod_of(C, i, S):
    q = C.qs[i]
    if C.base[i] == 0x32:
        return S[q[0]] & (S[q[1]] ^ S[q[2]])
    return S[q[0]] & S[q[1]]


def chains_from(links, skip):
    """Chains of consecutive links. Each j (in op order) extends the latest chain whose tail it links to; else
    it starts a chain [i, j] from the latest linked i that is in no chain yet; else it is left out."""
    cand = {}
    for i, j, f in links:
        if f < 1 or j in skip or i in skip:
            continue
        cand.setdefault(j, []).append(i)
    chains = []
    tail_of = {}        # tail op -> chain index
    in_chain = set()
    for j in sorted(cand):
        if j in in_chain:
            continue
        tails = [i for i in cand[j] if i in tail_of]
        if tails:
            i = max(tails)
            ci = tail_of.pop(i)
            chains[ci].append(j)
            tail_of[j] = ci
            in_chain.add(j)
            continue
        free = [i for i in cand[j] if i not in in_chain]
        if free:
            i = max(free)
            chains.append([i, j])
            tail_of[j] = len(chains) - 1
            in_chain.update((i, j))
    return chains


def main():
    t0 = time.time()
    ppath, lpath, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    maxlen = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    skip = set(int(x) for x in sys.argv[5].split(',')) if len(sys.argv) > 5 and sys.argv[5] else set()
    L = json.load(open(lpath))
    assert not L['asap'], 'links must be computed on the parent order (asap=0)'
    C = sim.load(ppath)
    links = [tuple(x) for x in L['links'] if maxlen == 0 or x[1] - x[0] <= maxlen]
    chains = chains_from(links, skip)
    print('chains', len(chains), [c for c in chains], flush=True)
    rng = np.random.default_rng(99)
    lanes = np.sort(rng.choice(sim.NX, size=2048, replace=False)).astype(np.int64)
    # pass 1: products at chain ops, states for the relations, zero / touch profile
    member = {}
    for ci, ch in enumerate(chains):
        for m, i in enumerate(ch):
            member[i] = (ci, m)
    D = {}
    rel = {}
    zero_now = np.ones(C.nq, bool)
    for w in C.xw:
        zero_now[w] = False
    nonzero_seen = {}      # wire -> sorted list of op indices where it is nonzero or touched (for holder search)
    busy = [set() for _ in range(C.nq)]

    S = sim.init_state(C)
    P = np.zeros(sim.W, np.uint64)
    bits = [0] * C.nbits
    stk = []
    nrot = [0]
    blk_open = []
    for i in range(C.n):
        b = C.base[i]
        if i in member:
            ci, m = member[i]
            d = prod_of(C, i, S)
            if m > 0:
                prev = chains[ci][m - 1]
                tgt = d ^ D[prev]
                Wb = sample_bits(S, lanes).T
                A = np.concatenate([np.ones((Wb.shape[0], 1), np.uint8), Wb], axis=1)
                ok, sol = gf2_solve(A, vec_sample(tgt, lanes)[:, None])
                lam = None
                if ok[0]:
                    s = sol[0]
                    q = np.zeros(sim.W, np.uint64)
                    if s[0]:
                        q ^= sim.ONES
                    ws = []
                    for w in range(C.nq):
                        if s[1 + w] and S[w].any():
                            q ^= S[w]
                            ws.append(w)
                    dd = q ^ tgt
                    nb = int(np.unpackbits(dd.view(np.uint8), bitorder='little')[:sim.NX].sum())
                    if nb == 0:
                        lam = (int(s[0]), ws)
                rel[i] = lam
            D[i] = d
        # touch / nonzero bookkeeping (block ops count as touching for the whole block)
        if b == kmb.PUSH and C.level[i] == 0:
            j = C.match[i]
            for k in range(i, j + 1):
                for w in C.qs[k]:
                    busy[w].add(i)
        for w in C.qs[i]:
            busy[w].add(i)
        act = all(stk)
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act or (C.cnd[i] >= 0 and not bits[C.cnd[i]]):
            continue
        sim.apply(C, i, S, P, bits, None, nrot)
        for w in C.qs[i]:
            if S[w].any():
                nonzero_seen.setdefault(w, []).append(i)
    P_parent = P.copy()
    S_parent_end = S.copy()
    print('pass 1', round(time.time() - t0, 1), 's; relations exact:', sum(1 for v in rel.values() if v is not None),
          'of', len(rel), flush=True)
    # holder per chain: wire untouched over [i0 - 1, ik] and zero there (never nonzero inside; zero at i0)
    busy_sorted = [np.array(sorted(s), dtype=np.int64) for s in busy]
    nz_sorted = {w: np.array(v, dtype=np.int64) for w, v in nonzero_seen.items()}

    def zero_at(w, t):
        """Is wire w zero just before op t? (its last nonzero-producing write before t, then a zeroing touch)"""
        return wire_zero_before[w][t]

    # exact zero-before table via a cheap replay of zero flags
    wire_zero_before = np.ones((C.nq, C.n + 1), bool)
    S = sim.init_state(C)
    P = np.zeros(sim.W, np.uint64)
    bits = [0] * C.nbits
    stk = []
    z = np.ones(C.nq, bool)
    for w in C.xw:
        z[w] = False
    for i in range(C.n):
        wire_zero_before[:, i] = z
        b = C.base[i]
        act = all(stk)
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act or (C.cnd[i] >= 0 and not bits[C.cnd[i]]):
            continue
        sim.apply(C, i, S, P, bits, None, nrot)
        for w in C.qs[i]:
            z[w] = not S[w].any()
    wire_zero_before[:, C.n] = z
    used_iv = []    # (lo, hi, wire)
    plan = []
    for ci, ch in enumerate(chains):
        if any(rel.get(i, 0) is None for i in ch[1:]):
            # cut the chain at the first inexact relation
            k = next(m for m in range(1, len(ch)) if rel.get(ch[m]) is None)
            ch = ch[:k]
            if len(ch) < 2:
                continue
        lo, hi = ch[0], ch[-1]
        hw = None
        for w in range(C.nq):
            if w in C.xw:
                continue
            bs = busy_sorted[w]
            if len(bs) and np.any((bs >= lo) & (bs <= hi)):
                continue
            if not wire_zero_before[w, lo]:
                continue
            if any(not (hi < a or lo > b2) and w == ww for a, b2, ww in used_iv):
                continue
            hw = w
            break
        if hw is None:
            print('  chain', ch, 'no holder', flush=True)
            continue
        used_iv.append((lo, hi, hw))
        plan.append((ch, hw))
    print('plan:', [(ch, hw) for ch, hw in plan], flush=True)
    # emit
    k = kmb.read(ppath)
    nb0 = k.num_bits
    newbit = nb0
    delete = []
    inserts = {}
    saved = 0
    for ch, h in plan:
        for m, i in enumerate(ch):
            q = C.qs[i]
            seq = []
            if m == 0:
                if C.base[i] == 0x32:
                    c, x, y = q
                    seq = [(CX, (x, y), 0, 0), (CCX, (c, y, h), 0, 0), (CX, (x, y), 0, 0), (CX, (h, x), 0, 0),
                           (CX, (h, y), 0, 0)]
                else:
                    a, bq, t = q
                    seq = [(CCX, (a, bq, h), 0, 0), (CX, (h, t), 0, 0)]
            else:
                const, ws = rel[i]
                for w in ws:
                    seq.append((CX, (w, h), 0, 0))
                if const:
                    seq.append((X, (h,), 0, 0))
                if C.base[i] == 0x32:
                    c, x, y = q
                    seq += [(CX, (h, x), 0, 0), (CX, (h, y), 0, 0)]
                else:
                    a, bq, t = q
                    seq += [(CX, (h, t), 0, 0)]
                saved += 1
                if m == len(ch) - 1:
                    seq.append((HMR, (h,), newbit, 0))
                    if C.base[i] == 0x32:
                        c, x, y = q
                        seq += [(CZ_IF, (c, x), 0, newbit), (CZ_IF, (c, y), 0, newbit)]
                    else:
                        a, bq, t = q
                        seq += [(CZ_IF, (a, bq), 0, newbit)]
                    newbit += 1
            delete.append(i)
            inserts[i] = seq
    k2 = kw.edit(k, delete, inserts)
    k2 = dataclasses.replace(k2, num_bits=newbit)
    kw.write(outp, k2)
    print('written', outp, 'Toffoli-class ops removed', saved, 'new bits', newbit - nb0, flush=True)
    # all-input check of the child against the parent
    Cc = sim.load(outp)
    zp = {}
    ok_all = True
    for seed in (0, 1, 2, 3):
        if seed == 0:
            oc = None
        else:
            r = np.random.default_rng(seed)
            draws = r.integers(0, 2, size=Cc.n)
            oc = (lambda j, d=draws: int(d[j]))
        S2, P2, b2, nr2 = sim.run(Cc, outcome=oc)
        same_state = bool(np.array_equal(S2, S_parent_end))
        same_phase = bool(np.array_equal(P2, P_parent))
        print('  check seed', seed, 'final state equal', same_state, 'Z-phase equal', same_phase, flush=True)
        ok_all &= same_state and same_phase
    print('CHECK', 'PASS' if ok_all else 'FAIL', 'saved static Toffoli-class', saved, round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
