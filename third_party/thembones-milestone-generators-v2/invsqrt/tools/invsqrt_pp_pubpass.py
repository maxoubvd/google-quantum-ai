#!/usr/bin/env python3
"""invsqrt_pp_pubpass.py -- workstream invsqrt_pp (P1, sixth holder): the PUBLIC exact rewrites on an invsqrt .kmb.

One round on the census of tools/invsqrt_pp_pubrw.py (exact on all 3,133,452 inputs):
  P2  CCX a b t ... CCX a' b' t, equal products, same block and condition, t untouched between: both deleted (-2).
  c3  CCX a b t, t only XOR-written until an unconditioned HMR t in the same block, a*b equal at the HMR on every
      lane (wires a, b of the CCX): CCX deleted, CZ a b conditioned on the HMR bit after the HMR (-1).
  c2  CCX a b t ... CCX a' b' t (equal products, same block and condition, t touched between) and a wire z zero
      and untouched over [i, k] (each z holds one pair at a time): CCX a b z; CX z t at i; CX z t; HMR z -> new bit m;
      CZ a' b' if m at k (-1).
Each op takes part in one rewrite per round; rounds repeat until nothing is found (or MAXR rounds).
Check of every round's child on all inputs: final state and Z-phase equal to the parent's for the zero transcript
and three random transcripts, every Z_POW sees the parent's value (tools/invsqrt_pp_affwrite.zpow_values).
Usage: python invsqrt_pp_pubpass.py <parent.kmb> <out.kmb> [maxrounds=4]
"""
import dataclasses
import os
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
import invsqrt_pp_pubrw as pr  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
from invsqrt_pp_affwrite import zpow_values  # noqa: E402

CCX, CX, CZ, HMR = 0x31, 0x21, 0x20, 0x50


def texp(Cc):
    return sum(0.5 ** (Cc.level[j] + (1 if Cc.cnd[j] >= 0 else 0)) for j in range(Cc.n) if Cc.base[j] in (0x30, 0x31, 0x32))


def one_round(src, dst):
    R = pr.main(src, quiet=True)
    C = sim.load(src)
    k0 = kmb.read(src)
    nb = k0.num_bits
    used = set()
    delete, inserts = [], {}

    def code(op, cnd):
        return op | (0x80 if cnd >= 0 else 0)
    # P2
    n2 = 0
    changed = []                                     # (wire, lo, hi): value differs from the parent inside (lo, hi)
    for (i, k, w) in R['P2_list']:
        if i in used or k in used:
            continue
        used.update((i, k))
        delete += [i, k]
        changed.append((C.qs[i][2], i, k))
        n2 += 1
    # c3 (same block, unconditioned CCX): candidates, then drop those whose new CZ reads a changed wire
    cand = []
    for (i, h, ok, w, same) in R['c3_list']:
        if ok and same and i not in used:
            cand.append((i, h))
            used.add(i)
    ch3 = changed + [(C.qs[i][2], i, h) for (i, h) in cand]
    n3 = 0
    for (i, h) in cand:
        a, b, t = C.qs[i]
        if any(w in (a, b) and lo < h < hi for (w, lo, hi) in ch3 if not (lo == i and hi == h)):
            used.discard(i)
            continue
        delete.append(i)
        inserts.setdefault(h, []).append((code(CZ, 0), (a, b), 0, C.bit[h]))
        n3 += 1
    # c2 (free wires assigned without overlap; ops not used above)
    rows = [(r[0], r[1], r[2], r[3], r[4]) for r in R['c2_list']]
    fit = pr.assign_c2(rows, taken_ops=used)
    # guard: a c3 CZ inserted inside [i, k] must not read the pair's z
    cz_at = {h: [q for (_, q, _, _) in v] for h, v in inserts.items()}
    n1 = 0
    for (i, k, z, w) in fit:
        if any(i < h < k and any(z in q for q in qs) for h, qs in cz_at.items()):
            continue
        a, b, t = C.qs[i]
        a2, b2, t2 = C.qs[k]
        cb = C.cnd[i]
        m = nb
        nb += 1
        delete += [i, k]
        inserts.setdefault(i, []).extend([(code(CCX, cb), (a, b, z), 0, cb), (code(CX, cb), (z, t), 0, cb)])
        inserts.setdefault(k, []).extend([(code(CX, cb), (z, t), 0, cb), (code(HMR, cb), (z,), m, cb),
                                          (code(CZ, 0) | 0x80, (a2, b2), 0, m)])
        used.update((i, k))
        n1 += 1
    if not delete:
        return 0, (0, 0, 0)
    # kw.edit takes (opcode, qubits, bit_t, cond); drop the unused fields' placeholders
    ins = {p: [(op, q, bt if (op & 0x40) else -1, cd if (op & 0x80) else -1) for (op, q, bt, cd) in v]
           for p, v in inserts.items()}
    kB = dataclasses.replace(kw.edit(k0, sorted(set(delete)), ins), num_bits=nb)
    kw.write(dst, kB)
    return len(delete), (n2, n3, n1)


def check(parent, child):
    C, D = sim.load(parent), sim.load(child)
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=D.n)
            oc = (lambda j, d=dr: int(d[j]))
        S1, P1, _, _ = sim.run(D, outcome=oc)
        e = bool(np.array_equal(S1, S0)) and bool(np.array_equal(P1, P0))
        ok &= e
    zp = zpow_values(C) == zpow_values(D)
    return ok and zp, texp(C), texp(D)


def main():
    t0 = time.time()
    src, dst = sys.argv[1], sys.argv[2]
    maxr = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    cur = src
    for r in range(maxr):
        nxt = dst + '.r%d' % r
        nd, (n2, n3, n1) = one_round(cur, nxt)
        if nd == 0:
            print('round %d: no sites' % r, flush=True)
            break
        ok, a, b = check(cur, nxt)
        print('round %d: P2 %d, c3 %d, c2 %d; expected Toffoli-class %.3f -> %.3f; CHECK %s (%.0fs)' % (
            r, n2, n3, n1, a, b, 'PASS' if ok else 'FAIL', time.time() - t0), flush=True)
        if not ok:
            print('round %d FAILED its check: stopping at the previous file' % r)
            break
        cur = nxt
    if cur == src:
        print('no change')
        return
    shutil.copy(cur, dst)
    ok, a, b = check(src, dst)
    print('FINAL %s: expected Toffoli-class %.3f -> %.3f; CHECK against the input %s (%.0fs)' % (
        dst, a, b, 'PASS' if ok else 'FAIL', time.time() - t0))


if __name__ == '__main__':
    main()
