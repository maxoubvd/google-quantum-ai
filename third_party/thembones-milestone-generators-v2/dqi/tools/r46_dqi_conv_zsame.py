#!/usr/bin/env python3
"""r46_dqi_conv_zsame.py -- R46-DQI-CONV: do all zero tests of one word slot give the SAME answer on every sample?

Input: a kmb, the list made by `r46_dqi_conv_ztests.py ... <stages> 0 <tests.npy>` (static index of the last AND of
each 12-bit zero test = the predicate [word == 0], slot id). Bit-parallel value model on valid verifier samples
(P1's, as tools/r46_dqi_conv_census.py). For each slot: samples on which the predicate was 1 at its first test, and
samples on which any later test of the slot disagreed with the first.
Usage: python r46_dqi_conv_zsame.py <dqi.kmb> <tests.npy> [runs=256] [seed=7]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km, _rng, popc64  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, grp, ns, first, differ):
    vv = vinit.copy()
    bits = binit.copy()
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    stack = np.zeros(256, np.uint64)
    ref = np.zeros(ns, np.uint64)
    seen = np.zeros(ns, np.uint8)
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
            if sp > 0:
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
        if (op & 0x80) != 0:
            m = m & bits[cd[ic]]
            ic += 1
        if base == 0x11:
            vv[a] ^= m
        elif base == 0x12:
            vv[a] &= ~m
        elif base == 0x21:
            vv[b] ^= vv[a] & m
        elif base == 0x22:
            d = (vv[a] ^ vv[b]) & m
            vv[a] ^= d
            vv[b] ^= d
        elif base == 0x31:
            d = vv[a] & vv[b] & m
            g = grp[k]
            if g >= 0:
                if seen[g] == 0:
                    seen[g] = 1
                    ref[g] = d
                    first[g] += popc64(d)
                else:
                    differ[g] |= d ^ ref[g]
            vv[t] ^= d
            k += 1
        elif base == 0x32:
            d = (vv[b] ^ vv[t]) & vv[a] & m
            vv[b] ^= d
            vv[t] ^= d
            k += 1
        elif base == 0x30:
            k += 1
        elif base == 0x50:
            s = _rng(s)
            bits[bit] = (bits[bit] & ~m) | (s & m)
            vv[a] &= ~m
        elif base == 0x40:
            bits[bit] &= ~m
        elif base == 0x41:
            bits[bit] |= m
        elif base == 0x42:
            bits[bit] ^= m
    return k


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    tests = np.load(sys.argv[2])
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    ops = np.ascontiguousarray(np.asarray(c.ops))
    nt = int((((ops & 0x7C) == 0x30) & (ops != 0x81)).sum())
    grp = -np.ones(nt, np.int64)
    grp[tests[:, 0]] = tests[:, 1]
    ns = int(tests[:, 1].max()) + 1
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    first = np.zeros(ns, np.int64)
    bad = np.zeros(ns, np.int64)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        differ = np.zeros(ns, np.uint64)
        k = simulate(ops, *arr, vinit, binit, seed * 1000 + r, grp, ns, first, differ)
        assert k == nt
        for g in range(ns):
            bad[g] += bin(int(differ[g])).count('1')
        if (r + 1) % 128 == 0 or r == R - 1:
            print('run %d (%d samples, %.0fs): predicate true at the first test per slot %s | samples with a later '
                  'test disagreeing %s' % (r + 1, 64 * (r + 1), time.time() - t0, first.tolist(), bad.tolist()),
                  flush=True)
    wires = {int(s): int(w) for _, s, w in tests}
    print('slots (first wire): ' + ', '.join('%d:q%d' % (s, wires[s]) for s in sorted(wires)))


if __name__ == '__main__':
    main()
