#!/usr/bin/env python3
"""dqi_pp_act.py -- workstream dqi_pp (P1): EXECUTED activity of every Toffoli-class op and every PUSH block.

Bit-parallel value simulation (phases ignored) of an emitted .kmb on R runs of 64 VALID verifier samples (syndromes
and masks from verify_submission_dqi.sample_shots with the verifier's weight law, HMR outcomes random), the same model
as tools/dqi_pp_fp.py without the fingerprint tests. For every Toffoli-class op (CCX, CCZ, CSWAP) it counts the
samples on which the op executes (its block conditions and its own condition all 1); for every PUSH_CONDITION it
counts the samples on which the block body runs. Mean executed Toffolis per shot = act.sum() / (64 R), which must
match the verifier's toffoli_count (calibration line printed at the end).
Usage: python dqi_pp_act.py <circuit.kmb> <out.npz> [runs=32] [seed=5]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

DATA = os.path.expanduser('~/qcc/data')
sys.path[:0] = [os.path.join(DATA, 'kickmix'), os.path.join(DATA, 'verifiers')]
import kickmix as km  # noqa: E402
from dqi import verify_submission_dqi as v  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def _rng(s):
    s ^= (s << np.uint64(13)) & ALL
    s ^= s >> np.uint64(7)
    s ^= (s << np.uint64(17)) & ALL
    return s


@njit(cache=True)
def popc64(x):
    x = x - ((x >> np.uint64(1)) & np.uint64(0x5555555555555555))
    x = (x & np.uint64(0x3333333333333333)) + ((x >> np.uint64(2)) & np.uint64(0x3333333333333333))
    x = (x + (x >> np.uint64(4))) & np.uint64(0x0F0F0F0F0F0F0F0F)
    return (x * np.uint64(0x0101010101010101)) >> np.uint64(56)


@njit(cache=True)
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, act, pact):
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
    k = 0
    kp = 0
    for oi in range(ops.shape[0]):
        op = ops[oi]
        if op == 0x81:
            c = cd[ic]
            ic += 1
            stack[sp] = active
            sp += 1
            active = active & bits[c]
            pact[kp] += popc64(active)
            kp += 1
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
            act[k] += popc64(m)
            v[t] ^= v[a] & v[b] & m
            k += 1
        elif base == 0x32:
            act[k] += popc64(m)
            d = (v[b] ^ v[t]) & v[a] & m
            v[b] ^= d
            v[t] ^= d
            k += 1
        elif base == 0x30:
            act[k] += popc64(m)
            k += 1
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
    return k, kp


def pack(boolmat):
    w = (np.uint64(1) << np.arange(64, dtype=np.uint64))
    return (boolmat.astype(np.uint64) * w[None, :]).sum(axis=1).astype(np.uint64)


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 32
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    ops = np.ascontiguousarray(np.asarray(c.ops))
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    t_op = np.flatnonzero(toff)
    p_op = np.flatnonzero(ops == 0x81)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    act = np.zeros(len(t_op), np.int64)
    pact = np.zeros(len(p_op), np.int64)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        k, kp = simulate(ops, *arr, vinit, binit, seed * 1000 + r, act, pact)
        assert k == len(t_op) and kp == len(p_op)
        if r % 8 == 7 or r == R - 1:
            print('run %d: mean executed Toffoli-class %.1f (%.1fs)' % (r, act.sum() / (64.0 * (r + 1)),
                                                                          time.time() - t0), flush=True)
    np.savez(out + '.tmp.npz', t_op=t_op, act=act, p_op=p_op, pact=pact, runs=R)
    os.replace(out + '.tmp.npz', out)
    print('calibration: mean executed %.2f per shot over %d shots; static %d' % (act.sum() / (64.0 * R), 64 * R,
                                                                                  len(t_op)))
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
