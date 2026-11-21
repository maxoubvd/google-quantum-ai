#!/usr/bin/env python3
"""dqi_int_ppfp.py -- workstream dqi_int (B5), gate-level post-pass census by Boolean fingerprints (no rewrite).

Simulates the VALUES of an emitted .kmb bit-parallel on 64 random inputs per run (syndrome qubits and mask bits random,
HMR outcomes random per sample; phases ignored), two independent runs. For every unconditioned CCX(a, b, t) outside
every PUSH_CONDITION block it records the fingerprints of a and b, and looks for a wire w != t that holds exactly
the product a*b at that moment in BOTH runs (then CCX -> CX(w, t) is exact with high confidence: -1 Toffoli, no
scratch). It also groups the CCX by their unordered pair of control functions (repeated products, for a retained
copy). Needs numba (box venv ~/qcc/.venv_work).
Usage: python dqi_int_ppfp.py <circuit.kmb> <out.npz> [seed=1]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def _rng(s):
    s ^= (s << np.uint64(13)) & ALL
    s ^= s >> np.uint64(7)
    s ^= (s << np.uint64(17)) & ALL
    return s


@njit(cache=True)
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, nq, nbits, seed, nfree):
    v = np.zeros(nq, np.uint64)
    bits = np.zeros(nbits, np.uint64)
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    for i in range(840):                          # syndrome qubits q0..q839
        s = _rng(s)
        v[i] = s
    for i in range(min(49140, nbits)):            # mask bits b0..b49139
        s = _rng(s)
        bits[i] = s
    fa = np.zeros(nfree, np.uint64)
    fb = np.zeros(nfree, np.uint64)
    hit = np.full(nfree, -1, np.int64)            # a wire holding a*b at that moment, -1 if none
    stack = np.zeros(256, np.uint64)
    sp = 0
    active = ALL
    i3 = 0
    i2 = 0
    i1 = 0
    ib = 0
    ic = 0
    k = 0
    for oi in range(ops.shape[0]):
        op = ops[oi]
        if op == 0x81:                            # PUSH_CONDITION
            c = cd[ic]
            ic += 1
            stack[sp] = active
            sp += 1
            active = active & bits[c]
            continue
        base = op & 0x7F
        if base == 0x01:                          # POP_CONDITION
            sp -= 1
            active = stack[sp]
            continue
        nqk = (op & 0x30) >> 4
        a = 0
        b = 0
        t = 0
        if nqk == 3:
            a = q3a[i3]
            b = q3b[i3]
            t = q3c[i3]
            i3 += 1
        elif nqk == 2:
            a = q2a[i2]
            b = q2b[i2]
            i2 += 1
        elif nqk == 1:
            a = q1[i1]
            i1 += 1
        bit = 0
        if (op & 0x40) != 0 and base != 0x06 and base != 0x07:
            bit = bt[ib]
            ib += 1
        m = active
        conditioned = (op & 0x80) != 0
        if conditioned:
            m = m & bits[cd[ic]]
            ic += 1
        if base == 0x11:                          # X
            v[a] ^= m
        elif base == 0x12:                        # R
            v[a] &= ~m
        elif base == 0x21:                        # CX a -> b
            v[b] ^= v[a] & m
        elif base == 0x22:                        # SWAP
            d = (v[a] ^ v[b]) & m
            v[a] ^= d
            v[b] ^= d
        elif base == 0x31:                        # CCX a b -> t
            if sp == 0 and not conditioned:
                p = v[a] & v[b]
                fa[k] = v[a]
                fb[k] = v[b]
                for w in range(nq):
                    if w != t and v[w] == p:
                        hit[k] = w
                        break
                k += 1
            v[t] ^= v[a] & v[b] & m
        elif base == 0x32:                        # CSWAP c x y
            d = (v[b] ^ v[t]) & v[a] & m
            v[b] ^= d
            v[t] ^= d
        elif base == 0x50:                        # HMR q -> bit
            s = _rng(s)
            bits[bit] = (bits[bit] & ~m) | (s & m)
            v[a] &= ~m
        elif base == 0x40:                        # BIT_STORE0
            bits[bit] &= ~m
        elif base == 0x41:                        # BIT_STORE1
            bits[bit] |= m
        elif base == 0x42:                        # BIT_INVERT
            bits[bit] ^= m
    return fa, fb, hit, k


def popc(x):
    x = x.astype(np.uint64)
    c = np.zeros(len(x), np.int64)
    for i in range(64):
        c += ((x >> np.uint64(i)) & np.uint64(1)).astype(np.int64)
    return c


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    ops = np.ascontiguousarray(np.asarray(c.ops))
    base = ops & 0x7F
    push = ops == 0x81
    pop = ops == 0x01
    depth = np.cumsum(push.astype(np.int64) - pop.astype(np.int64))
    free_ccx = np.flatnonzero((base == 0x31) & ((ops & 0x80) == 0) & (depth == 0))
    nfree = len(free_ccx)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    runs = []
    for r in range(2):
        fa, fb, hit, k = simulate(ops, *arr, c.num_qubits, c.num_bits, seed * 1000 + r, nfree)
        assert k == nfree, (k, nfree)
        runs.append((fa, fb, hit))
        print('run', r, 'free CCX', nfree, 'hits', int((hit >= 0).sum()), '%.1fs' % (time.time() - t0), flush=True)
    (fa0, fb0, h0), (fa1, fb1, h1) = runs
    both = (h0 >= 0) & (h0 == h1)
    dense = (popc(fa0 & fb0) >= 6) & (popc(fa1 & fb1) >= 6)
    print('CCX -> CX(w, t) sites, same wire in both runs: %d (product dense: %d)' % (int(both.sum()),
          int((both & dense).sum())), flush=True)
    # repeated products by control-function pair (both runs)
    lo0, hi0 = np.minimum(fa0, fb0), np.maximum(fa0, fb0)
    lo1, hi1 = np.minimum(fa1, fb1), np.maximum(fa1, fb1)
    keys = np.stack([lo0, hi0, lo1, hi1], axis=1)
    order = np.lexsort(keys.T[::-1])
    ks = keys[order]
    same = np.all(ks[1:] == ks[:-1], axis=1)
    nontriv = (popc(lo0 & hi0)[order[1:]] >= 6)
    print('repeated (dense) products: %d consecutive pairs' % int((same & nontriv).sum()), flush=True)
    np.savez(out, free_ccx=free_ccx, fa0=fa0, fb0=fb0, fa1=fa1, fb1=fb1, h0=h0, h1=h1, both=both, dense=dense,
             pair_first=free_ccx[order[:-1][same & nontriv]], pair_second=free_ccx[order[1:][same & nontriv]])
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
