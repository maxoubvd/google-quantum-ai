#!/usr/bin/env python3
"""dqi_int_pprep.py -- workstream dqi_int (B5): exactness filter for repeated-product pairs (no rewrite).

For candidate pairs (first, second) of unconditioned free CCX ops, checks over R runs of 64 samples whether the two
control functions at `first` equal the two at `second` (as an unordered pair) in EVERY sample of every run. Runs cycle
through three input kinds: valid inputs by the DQI weight law, valid inputs with the error weight uniform on 18..35
(policy W's supported range), and uniformly random inputs. Pairs are read from an npz (pair_first, pair_second) and
optionally restricted to a mask (`--covered` uses tools/dqi_int_ppgap.py's rule inline).
Usage: python dqi_int_pprep.py <circuit.kmb> <fp.npz> <out.npz> [runs=192] [seed=5]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_int_ppfp2 import _rng, pack, v, km  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def controls(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, nfree):
    """Fingerprints (a, b) of the controls of every free CCX, in order."""
    v = vinit.copy()
    bits = binit.copy()
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    fa = np.zeros(nfree, np.uint64)
    fb = np.zeros(nfree, np.uint64)
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
        if op == 0x81:
            c = cd[ic]
            ic += 1
            stack[sp] = active
            sp += 1
            active = active & bits[c]
            continue
        base = op & 0x7F
        if base == 0x01:
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
        if (op & 0x40) != 0:
            bit = bt[ib]
            ib += 1
        m = active
        conditioned = (op & 0x80) != 0
        if conditioned:
            m = m & bits[cd[ic]]
            ic += 1
        if base == 0x11:
            v[a] ^= m
        elif base == 0x12:
            v[a] &= ~m
        elif base == 0x21:
            v[b] ^= v[a] & m
        elif base == 0x22:
            d = (v[a] ^ v[b]) & m
            v[a] ^= d
            v[b] ^= d
        elif base == 0x31:
            if sp == 0 and not conditioned:
                fa[k] = v[a]
                fb[k] = v[b]
                k += 1
            v[t] ^= v[a] & v[b] & m
        elif base == 0x32:
            d = (v[b] ^ v[t]) & v[a] & m
            v[b] ^= d
            v[t] ^= d
        elif base == 0x50:
            s = _rng(s)
            bits[bit] = (bits[bit] & ~m) | (s & m)
            v[a] &= ~m
        elif base == 0x40:
            bits[bit] &= ~m
        elif base == 0x41:
            bits[bit] |= m
        elif base == 0x42:
            bits[bit] ^= m
    return fa, fb, k


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    fp = np.load(sys.argv[2])
    out = sys.argv[3]
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 192
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    ops = np.ascontiguousarray(np.asarray(c.ops))
    base = ops & 0x7F
    depth = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
    free_ccx = np.flatnonzero((base == 0x31) & ((ops & 0x80) == 0) & (depth == 0))
    pos = {int(o): i for i, o in enumerate(free_ccx)}
    pf = np.array([pos[int(x)] for x in fp['pair_first']], dtype=np.int64)
    ps = np.array([pos[int(x)] for x in fp['pair_second']], dtype=np.int64)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wlaw = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    wuni = np.zeros_like(wlaw)
    wuni[18:36] = 1.0
    wuni /= wuni.sum()
    rng = np.random.default_rng(seed)
    alive = np.ones(len(pf), dtype=bool)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        kind = r % 3
        if kind < 2:
            syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wlaw if kind == 0 else wuni)
            vinit[:840] = pack(v.to_bits(syn, spec.degree))
            mb = pack(v.to_bits(mbuc, spec.degree))
            binit[:len(mb)] = mb
        else:
            vinit[:840] = rng.integers(0, 2 ** 62, size=840).astype(np.uint64) * np.uint64(4) + \
                rng.integers(0, 4, size=840).astype(np.uint64)
            binit[:49140] = rng.integers(0, 2 ** 62, size=49140).astype(np.uint64) * np.uint64(4) + \
                rng.integers(0, 4, size=49140).astype(np.uint64)
        fa, fb, k = controls(ops, *arr, vinit, binit, seed * 1000 + r, len(free_ccx))
        lo, hi = np.minimum(fa, fb), np.maximum(fa, fb)
        alive &= (lo[pf] == lo[ps]) & (hi[pf] == hi[ps])
        if r % 16 == 15:
            print('run %d: pairs still exact %d of %d %.1fs' % (r, int(alive.sum()), len(pf), time.time() - t0),
                  flush=True)
    np.savez(out, pair_first=fp['pair_first'][alive], pair_second=fp['pair_second'][alive])
    print('exact pairs kept: %d' % int(alive.sum()))


if __name__ == '__main__':
    main()
