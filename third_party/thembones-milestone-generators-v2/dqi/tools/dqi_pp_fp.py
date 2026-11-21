#!/usr/bin/env python3
"""dqi_pp_fp.py -- workstream dqi_pp (P1): fingerprint CANDIDATES for the exact prover (tools/dqi_pp_prove.py).

Bit-parallel value simulation (phases ignored) of an emitted .kmb on R runs of 64 VALID verifier samples (syndromes
and masks from verify_submission_dqi.sample_shots, HMR outcomes random). For every Toffoli-class op (CCX, CCZ,
CSWAP; free or inside PUSH_CONDITION blocks), over the samples where the op is active:
  Z    the product is never 1 (CCX a*b, CCZ a*b*c, CSWAP c*(x+y))      -> delete, if proved
  EA   CCX: a*b == a on every sample (a implies b); EB: a*b == b       -> CX(a, t) / CX(b, t), if proved
  U    CCX: t == a*b before the op                                     -> HMR(t) + CZ(a, b) if m, if proved
  H    CCX: some wire w != t equals a*b (lowest such w, same in every run) -> CX(w, t), if proved
  ones number of samples where the product was 1; act = samples where the op was active
These are candidates only: a near-identity passes them. Nothing is rewritten without the prover.
Usage: python dqi_pp_fp.py <circuit.kmb> <out.npz> [runs=64] [seed=11]
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
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, first,
             fz, fea, feb, fu, fh, ones, act):
    v = vinit.copy()
    bits = binit.copy()
    nq = v.shape[0]
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
            p = v[a] & v[b]
            act[k] += popc64(m)
            ones[k] += popc64(p & m)
            if (p & m) != 0:
                fz[k] = False
            if ((p ^ v[a]) & m) != 0:
                fea[k] = False
            if ((p ^ v[b]) & m) != 0:
                feb[k] = False
            if ((p ^ v[t]) & m) != 0:
                fu[k] = False
            if first:
                fh[k] = -1
                for w in range(nq):
                    if w != t and ((v[w] ^ p) & m) == 0:
                        fh[k] = w
                        break
            elif fh[k] >= 0:
                w = fh[k]
                if ((v[w] ^ p) & m) != 0:
                    fh[k] = -1
            v[t] ^= p & m
            k += 1
        elif base == 0x32:
            d = (v[b] ^ v[t]) & v[a]
            act[k] += popc64(m)
            ones[k] += popc64(d & m)
            if (d & m) != 0:
                fz[k] = False
            d &= m
            v[b] ^= d
            v[t] ^= d
            k += 1
        elif base == 0x30:
            p = v[a] & v[b] & v[t]
            act[k] += popc64(m)
            ones[k] += popc64(p & m)
            if (p & m) != 0:
                fz[k] = False
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
    return k


def pack(boolmat):
    w = (np.uint64(1) << np.arange(64, dtype=np.uint64))
    return (boolmat.astype(np.uint64) * w[None, :]).sum(axis=1).astype(np.uint64)


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 11
    ops = np.ascontiguousarray(np.asarray(c.ops))
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    t_op = np.flatnonzero(toff)
    nt = len(t_op)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    fz = np.ones(nt, np.bool_)
    fea = np.ones(nt, np.bool_)
    feb = np.ones(nt, np.bool_)
    fu = np.ones(nt, np.bool_)
    fh = np.full(nt, -1, np.int64)
    ones = np.zeros(nt, np.int64)
    act = np.zeros(nt, np.int64)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        k = simulate(ops, *arr, vinit, binit, seed * 1000 + r, r == 0, fz, fea, feb, fu, fh, ones, act)
        assert k == nt
        if r % 8 == 7 or r == R - 1:
            print('run %d: Z %d  EA %d  EB %d  U %d  H %d  (%.1fs)' % (
                r, int(fz.sum()), int(fea.sum()), int(feb.sum()), int(fu.sum()), int((fh >= 0).sum()),
                time.time() - t0), flush=True)
    np.savez(out, t_op=t_op, fz=fz, fea=fea, feb=feb, fu=fu, fh=fh, ones=ones, act=act, runs=R)
    base = ops[t_op] & 0x7F
    S = R * 64.0
    for code, nm in ((0x31, 'CCX'), (0x30, 'CCZ'), (0x32, 'CSWAP')):
        m = base == code
        for fl, fn in ((fz, 'Z'), (fea & ~fz, 'EA'), (feb & ~fz & ~fea, 'EB'), (fu & ~fz, 'U'),
                       ((fh >= 0) & ~fz & ~fu, 'H')):
            mm = m & fl
            if mm.any():
                print('%-5s %-2s sites %7d  expected (measured activity) %9.1f  of which free %d' % (
                    nm, fn, int(mm.sum()), float(act[mm].sum()) / S, int((mm & (act == R * 64)).sum())))
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
