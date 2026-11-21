#!/usr/bin/env python3
"""invsqrt_pp_pubc3q.py -- workstream invsqrt_pp (P1, sixth holder): dead writes whose operands changed (c3, extended).

Site: a top-level unconditioned CCX a b t whose target is only XOR-written until an unconditioned top-level HMR t
(the c3 structure of tools/invsqrt_pp_pubrw.py), where a*b at the HMR is NOT a*b of the CCX's wires any more.  The
kick m * p (p = a*b of the CCX time) is paid at the HMR by a QUADRATIC FORM of the wires there (t excluded):
  p = c0 + sum Z-wires w + sum pairs x*y           (solved on 4,096 sample lanes, then checked on all 3,133,452)
Rewrite: delete the CCX; after the HMR: PUSH m; CZ x y ...; Z w ...; NEG (c0); POP.  All Clifford: -1 per site.
Sites in one round do not read a wire another site of the round changes at that moment.  The child is checked on
all inputs (tools/invsqrt_pp_pubpass.check).
Usage: python invsqrt_pp_pubc3q.py <parent.kmb> <out.kmb>
"""
import dataclasses
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
import invsqrt_pp_pubrw as pr  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
from invsqrt_pp_blocks import gf2_solve, quad_basis, sample_bits, vec_sample  # noqa: E402
from invsqrt_pp_pubpass import check  # noqa: E402

CZ, Z, NEG, PUSH, POP = 0x20, 0x10, 0x00, 0x81, 0x01


def find(src):
    R = pr.main(src, quiet=True)
    C = sim.load(src)
    C.path = src
    ctx = pr.contexts(C)
    cand = [(i, h) for (i, h, ok, w, same) in R['c3_list'] if not ok and same and len(ctx[i]) == 0 and len(ctx[h]) == 0]
    at_i = {i: h for (i, h) in cand}
    at_h = {}
    for (i, h) in cand:
        at_h.setdefault(h, []).append(i)
    rng = np.random.default_rng(7)
    lanes = np.sort(rng.choice(sim.NX, size=4096, replace=False)).astype(np.int64)
    S, kind, regs = pr.lanes(C)
    bits = [0] * C.nbits
    stk = []
    P = {}
    found = []
    for j in range(C.n):
        b = C.base[j]
        act = all(stk)
        if j in at_h:
            for i in at_h[j]:
                t = C.qs[i][2]
                wires = [w for w in range(C.nq) if w != t]
                Wb = sample_bits(S[wires], lanes).T
                A, cols = quad_basis(Wb)
                g = vec_sample(P[i], lanes).reshape(-1, 1)
                okk, sol = gf2_solve(A, g)
                if not okk[0]:
                    continue
                terms = [cols[c] for c in np.flatnonzero(sol[0])]
                # exact check on every lane
                v = np.zeros_like(P[i])
                c0 = 0
                zs, cz = [], []
                for tm in terms:
                    if tm[0] == '1':
                        c0 ^= 1
                    elif tm[0] == 'w':
                        zs.append(wires[tm[1]])
                        v ^= S[wires[tm[1]]]
                    else:
                        x, y = wires[tm[1]], wires[tm[2]]
                        cz.append((x, y))
                        v ^= S[x] & S[y]
                if c0:
                    v ^= sim.ONES
                if np.array_equal(v, P[i]):
                    found.append((i, j, c0, zs, cz))
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[j]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act or (C.cnd[j] >= 0 and not bits[C.cnd[j]]):
            continue
        if j in at_i:
            P[j] = S[C.qs[j][0]] & S[C.qs[j][1]]
        if b == 0x50:
            S[C.qs[j][0]][:] = 0
            bits[C.bit[j]] = 1
        else:
            sim.apply(C, j, S, np.zeros(S.shape[1], np.uint64), bits, None, [0])
    return C, cand, found


def main():
    t0 = time.time()
    src, dst = sys.argv[1], sys.argv[2]
    C, cand, found = find(src)
    print('candidates %d, quadratic-form sites %d (%.0fs)' % (len(cand), len(found), time.time() - t0), flush=True)
    changed = [(C.qs[i][2], i, h) for (i, h, c0, zs, cz) in found]
    k0 = kmb.read(src)
    delete, inserts = [], {}
    n = 0
    for (i, h, c0, zs, cz) in found:
        rd = set(zs) | {w for p in cz for w in p}
        if any(w in rd and lo < h < hi for (w, lo, hi) in changed if not (lo == i and hi == h)):
            continue
        m = C.bit[h]
        seq = [(PUSH, (), -1, m)] + [(CZ, p, -1, -1) for p in cz] + [(Z, (w,), -1, -1) for w in zs] + \
              ([(NEG, (), -1, -1)] if c0 else []) + [(POP, (), -1, -1)]
        delete.append(i)
        inserts.setdefault(h, []).extend(seq)
        n += 1
    if not n:
        print('no change')
        return
    kw.write(dst, kw.edit(k0, sorted(delete), inserts))
    ok, a, b = check(src, dst)
    print('rewritten %d; expected Toffoli-class %.3f -> %.3f; CHECK %s (%.0fs)' % (n, a, b, 'PASS' if ok else 'FAIL',
                                                                                   time.time() - t0))
    if not ok:
        os.remove(dst)


if __name__ == '__main__':
    main()
