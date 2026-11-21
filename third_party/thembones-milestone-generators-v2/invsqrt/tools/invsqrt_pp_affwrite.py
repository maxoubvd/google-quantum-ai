#!/usr/bin/env python3
"""invsqrt_pp_affwrite.py -- workstream invsqrt_pp (P1): rewrite the affine sites of an invsqrt .kmb (0 Toffolis).

For every top-level, unconditioned Toffoli-class op named in the audit file (invsqrt_pp_affine.py) the relation
  product = c0 + sum_{w in R} w        (exact on all 3,133,452 inputs, wire values just before the op)
is solved again with the wires that make it (exact Gaussian elimination with tracked combinations), and the op is
replaced by Clifford gates:
  CCX a b t, t not in R   -> X t (if c0); CX w t (w in R)                                   t ^= product
  CCX a b t, t in R       -> HMR t -> new bit m; if m: CZ a b, Z w (w in R \\ t), NEG (if c0);
                             then CX w t (w in R \\ t), X t (if c0)                          t <- t + product
  CSWAP c x y, x, y not in R -> (X x; X y if c0); CX w x; CX w y (w in R)                    x, y ^= d
Sites whose relation needs anything else are skipped and listed. The child is checked on all inputs: final state
and Z-phase equal to the parent's for the zero transcript and three random transcripts, every Z_POW sees the
parent's value.
Usage: python invsqrt_pp_affwrite.py <parent.kmb> <affine.json> <out.kmb>
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
from invsqrt_pp_asap import zpow_values  # noqa: E402

X_, CX, CCX, CSWAP, CZ, Z, HMR, NEG, PUSH, POP = 0x11, 0x21, 0x31, 0x32, 0x20, 0x10, 0x50, 0x00, 0x81, 0x01
LASTMASK = np.uint64((1 << (sim.NX - 64 * (sim.W - 1))) - 1)


def solve_affine(S, p):
    """p = c0 + sum_{w in R} S[w] on all real lanes? Returns (c0, R) or None. Exact elimination with tracking."""
    nq = S.shape[0]
    rows = [np.full(sim.W, sim.ONES, np.uint64)] + [S[w].copy() for w in range(nq)]
    for r in rows:
        r[-1] &= LASTMASK
    tags = [1 << 0] + [1 << (w + 1) for w in range(nq)]          # combination of (1, w0, w1, ...) per row
    piv = []                                                    # (word, bit, row vector, tag)
    for v, tg in zip(rows, tags):
        v = v.copy()
        for (w, bb, pv, pt) in piv:
            if v[w] & bb:
                v ^= pv
                tg ^= pt
        nz = np.flatnonzero(v)
        if len(nz) == 0:
            continue
        w = int(nz[0])
        bb = np.uint64(int(v[w]) & -int(v[w]))
        piv.append((w, bb, v, tg))
    t = p.copy()
    t[-1] &= LASTMASK
    comb = 0
    for (w, bb, pv, pt) in piv:
        if t[w] & bb:
            t ^= pv
            comb ^= pt
    if t.any():
        return None
    c0 = comb & 1
    R = [w for w in range(nq) if (comb >> (w + 1)) & 1]
    return c0, R


def main():
    t0 = time.time()
    src, apath, dst = sys.argv[1], sys.argv[2], sys.argv[3]
    C = sim.load(src)
    k0 = kmb.read(src)
    sites = [s for s in json.load(open(apath)) if s.get('ctx') == 'top' and s['kind'] in ('CCX', 'CSWAP')]
    want = set(s['op'] for s in sites)
    rel = {}

    def hook(i, S, P, bits, act):
        if i in want and act:
            q = C.qs[i]
            if C.base[i] == CCX:
                p = S[q[0]] & S[q[1]]
            else:
                p = S[q[0]] & (S[q[1]] ^ S[q[2]])
            rel[i] = solve_affine(S, p)
    sim.run(C, hook=hook)
    nb = k0.num_bits
    delete, inserts, skipped, done = [], {}, [], []
    for i in sorted(want):
        r = rel.get(i)
        q = C.qs[i]
        if r is None:
            skipped.append((i, 'no relation'))
            continue
        c0, R = r
        if C.base[i] == CCX:
            a, b, t = q
            if t not in R:
                seq = ([(X_, (t,), 0, 0)] if c0 else []) + [(CX, (w, t), 0, 0) for w in R]
            else:
                m = nb
                nb += 1
                Rt = [w for w in R if w != t]
                seq = [(HMR, (t,), m, 0), (PUSH, (), 0, m), (CZ, (a, b), 0, 0)] + \
                      [(Z, (w,), 0, 0) for w in Rt] + ([(NEG, (), 0, 0)] if c0 else []) + [(POP, (), 0, 0)] + \
                      [(CX, (w, t), 0, 0) for w in Rt] + ([(X_, (t,), 0, 0)] if c0 else [])
        else:
            c, x, y = q
            if x in R or y in R:
                skipped.append((i, 'CSWAP target in relation'))
                continue
            seq = []
            for tt in (x, y):
                seq += ([(X_, (tt,), 0, 0)] if c0 else []) + [(CX, (w, tt), 0, 0) for w in R]
        delete.append(i)
        inserts[i] = seq
        done.append((i, C.base[i], c0, R))
    print('sites %d, rewritten %d, skipped %s' % (len(want), len(done), skipped), flush=True)
    for d in done[:40]:
        print('  op %d %s c0 %d R %s' % (d[0], 'CCX' if d[1] == CCX else 'CSWAP', d[2], d[3]))
    kB = dataclasses.replace(kw.edit(k0, delete, inserts), num_bits=nb)
    kw.write(dst, kB)
    D = sim.load(dst)
    texp = lambda Cc: sum(0.5 ** (Cc.level[j] + (1 if Cc.cnd[j] >= 0 else 0)) for j in range(Cc.n)  # noqa: E731
                          if Cc.base[j] in (0x30, 0x31, 0x32))
    print('expected Toffoli-class: %.3f -> %.3f' % (texp(C), texp(D)), flush=True)
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=D.n)
            oc = (lambda j, d=dr: int(d[j]))
        S1_, P1_, _, _ = sim.run(D, outcome=oc)
        e = bool(np.array_equal(S1_, S0)) and bool(np.array_equal(P1_, P0))
        print('  seed', seed, 'final state and Z-phase equal to the parent:', e, flush=True)
        ok &= e
    zp = zpow_values(C) == zpow_values(D)
    print('  every Z_POW sees the parent value:', zp)
    print('CHECK', 'PASS' if ok and zp else 'FAIL', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
