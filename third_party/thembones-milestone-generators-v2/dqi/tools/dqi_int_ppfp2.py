#!/usr/bin/env python3
"""dqi_int_ppfp2.py -- workstream dqi_int (B5): strict filter for the CCX -> CX(w, t) sites (no rewrite).

As tools/dqi_int_ppfp.py, but over R runs of 64 samples each: even runs take VALID inputs from the official
verifier's own sampler (syndromes of errors drawn by the DQI weight law, masks, as the verifier writes them), odd runs
take uniformly random inputs. A site is kept when the same wire w != t holds the product a*b at that moment in EVERY
run, and the product is dense (popcount >= 6 of 64) in every random run.
Usage: python dqi_int_ppfp2.py <circuit.kmb> <out.npz> [runs=32] [seed=7]
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
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, nfree, sets, first, unc):
    v = vinit.copy()
    bits = binit.copy()
    nq = v.shape[0]
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    pv = np.zeros(nfree, np.uint64)
    hit = np.zeros(nfree, np.int64)
    nw = sets.shape[1]
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
                p = v[a] & v[b]
                pv[k] = p
                if v[t] == p:
                    unc[k] = True
                for j in range(nw):
                    word = np.uint64(0)
                    for u in range(64):
                        w = j * 64 + u
                        if w < nq and w != t and v[w] == p:
                            word |= np.uint64(1) << np.uint64(u)
                    if first:
                        sets[k, j] = word
                    else:
                        sets[k, j] &= word
                    if sets[k, j] != 0:
                        hit[k] = 1
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
    return pv, hit, k


def pack(boolmat):
    """(rows, 64) bool -> (rows,) uint64 with sample s in bit s."""
    w = (np.uint64(1) << np.arange(64, dtype=np.uint64))
    return (boolmat.astype(np.uint64) * w[None, :]).sum(axis=1).astype(np.uint64)


def popc(x):
    c = np.zeros(len(x), np.int64)
    for i in range(64):
        c += ((x >> np.uint64(i)) & np.uint64(1)).astype(np.int64)
    return c


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 32
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 7
    ops = np.ascontiguousarray(np.asarray(c.ops))
    base = ops & 0x7F
    depth = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
    free_ccx = np.flatnonzero((base == 0x31) & ((ops & 0x80) == 0) & (depth == 0))
    nfree = len(free_ccx)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    nw = (c.num_qubits + 63) // 64
    sets = np.zeros((nfree, nw), np.uint64)
    dense = np.ones(nfree, dtype=bool)
    uncall = np.ones(nfree, dtype=bool)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        if r % 2 == 0:
            syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
            vinit[:840] = pack(v.to_bits(syn, spec.degree))
            mb = pack(v.to_bits(mbuc, spec.degree))
            binit[:len(mb)] = mb
        else:
            vinit[:840] = rng.integers(0, 1 << 63, size=840, dtype=np.int64).astype(np.uint64) * np.uint64(2) + \
                rng.integers(0, 2, size=840).astype(np.uint64)
            binit[:49140] = rng.integers(0, 1 << 63, size=49140, dtype=np.int64).astype(np.uint64) * np.uint64(2) \
                + rng.integers(0, 2, size=49140).astype(np.uint64)
        unc_r = np.zeros(nfree, dtype=np.bool_)
        pv, hit, k = simulate(ops, *arr, vinit, binit, seed * 1000 + r, nfree, sets, r == 0, unc_r)
        uncall &= unc_r
        assert k == nfree
        if r % 2 == 1:
            dense &= popc(pv) >= 6
        alive = hit > 0
        print('run %d (%s): sites with a common matching wire %d (dense %d); uncompute sites %d (dense %d) %.1fs' % (
            r, 'valid' if r % 2 == 0 else 'random', int(alive.sum()), int((alive & dense).sum()),
            int(uncall.sum()), int((uncall & dense).sum()), time.time() - t0), flush=True)
    alive = sets.any(axis=1)
    # lowest common wire per site
    w_all = np.full(nfree, -1, np.int64)
    for j in range(nw - 1, -1, -1):
        col = sets[:, j]
        nz = col != 0
        low = np.zeros(nfree, np.int64)
        x = col.copy()
        for u in range(63, -1, -1):
            bit = ((x >> np.uint64(u)) & np.uint64(1)) == 1
            low = np.where(bit, u, low)
        w_all = np.where(nz, j * 64 + low, w_all)
    keep = alive & dense
    np.savez(out, free_ccx=free_ccx, w=w_all, keep=keep, sparse_consistent=alive & ~dense, unc=uncall,
             unc_dense=uncall & dense)
    print('kept %d sites (a common matching wire in every run, dense); common but sparse %d' % (int(keep.sum()),
          int((alive & ~dense).sum())))


if __name__ == '__main__':
    main()
