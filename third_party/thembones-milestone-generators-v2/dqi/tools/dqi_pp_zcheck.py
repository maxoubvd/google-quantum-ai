#!/usr/bin/env python3
"""dqi_pp_zcheck.py -- workstream dqi_pp (P1): sample check that every scratch wire of a retained pass is zero.

Reads the checkpoints (after-op time t, wire w) written by tools/dqi_pp_ret.py (<out>.assign.npz) and simulates the
INPUT circuit on R runs of 64 samples (even runs valid verifier samples, odd runs uniform inputs, HMR outcomes
random): a checkpoint fails when w is nonzero right after op t on some sample. Writes the failing (wire, time) for
DQI_PP_DROP. Necessary condition only (the official verifier and the paired-seed test follow).
Usage: python dqi_pp_zcheck.py <input circuit.kmb> <assign.npz> <out drop.npz> [R=192]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_fp as F  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def sim_check(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, ct, cw, bad):
    v = vinit.copy()
    bits = binit.copy()
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    stack = np.zeros(256, np.uint64)
    sp = 0
    active = ALL
    i3 = 0
    i2 = 0
    i1 = 0
    ib = 0
    ic = 0
    q = 0
    nc = ct.shape[0]
    while q < nc and ct[q] < 0:
        if v[cw[q]] != 0:
            bad[q] = True
        q += 1
    for oi in range(ops.shape[0]):
        op = ops[oi]
        base = op & 0x7F
        if op == 0x81:
            c = cd[ic]
            ic += 1
            stack[sp] = active
            sp += 1
            active = active & bits[c]
        elif base == 0x01:
            if sp > 0:
                sp -= 1
                active = stack[sp]
        else:
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
                v[t] ^= v[a] & v[b] & m
            elif base == 0x32:
                d = (v[b] ^ v[t]) & v[a] & m
                v[b] ^= d
                v[t] ^= d
            elif base == 0x50:
                s = F._rng(s)
                bits[bit] = (bits[bit] & ~m) | (s & m)
                v[a] &= ~m
            elif base == 0x40:
                bits[bit] &= ~m
            elif base == 0x41:
                bits[bit] |= m
            elif base == 0x42:
                bits[bit] ^= m
        while q < nc and ct[q] == oi:
            if v[cw[q]] != 0:
                bad[q] = True
            q += 1
    return q


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    A = np.load(sys.argv[2])
    out = sys.argv[3]
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 192
    chk = A['chk']
    o = np.argsort(chk[:, 0], kind='stable')
    ct, cw = chk[o, 0].astype(np.int64), chk[o, 1].astype(np.int64)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = F.v.ProblemSpecDQI.benchmark_instance(num_shots=F.v.DEFAULT_NUM_SHOTS)
    tables = F.v.FieldTables(F.km.GF2Field(spec.degree))
    wprob = F.v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(29)
    bad = np.zeros(len(ct), np.bool_)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        if r % 2 == 0:
            syn, mbuc, _ = F.v.sample_shots(spec, tables, rng, 64, wprob)
            vinit[:840] = F.pack(F.v.to_bits(syn, spec.degree))
            mb = F.pack(F.v.to_bits(mbuc, spec.degree))
            binit[:len(mb)] = mb
        else:
            vinit[:840] = rng.integers(0, 2 ** 63, size=840, dtype=np.int64).astype(np.uint64) * np.uint64(2) + \
                rng.integers(0, 2, size=840).astype(np.uint64)
            binit[:49140] = rng.integers(0, 2 ** 63, size=49140, dtype=np.int64).astype(np.uint64) * np.uint64(2) \
                + rng.integers(0, 2, size=49140).astype(np.uint64)
        sim_check(ops, *arr, vinit, binit, 900 + r, ct, cw, bad)
    print('checkpoints %d, nonzero on some sample %d (%.1fs)' % (len(ct), int(bad.sum()), time.time() - t0))
    np.savez(out, wire=cw[bad], time=ct[bad])


if __name__ == '__main__':
    main()
