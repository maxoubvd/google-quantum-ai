#!/usr/bin/env python3
"""dqi_pp_pairs2.py -- workstream dqi_pp (P1): retained SUMS on fingerprint pairs (no rewrite).

Product key of every free CCX = rolling hash of its product value over R runs of 64 samples (even runs VALID
verifier samples, odd runs uniform inputs, HMR outcomes random); popcount of the product over the uniform runs kept
as `dens`. Pairs: consecutive occurrences of one key. XOR segment: a maximal run of events on a wire in which every
event only XORs INTO the wire (target of X, CX, CCX, conditioned or not): such writes commute and the wire is not read.
Group = pairs with first CCX in XOR segment F and second in XOR segment M (F != M); one zero-idle scratch wire from the
group's first forward CCX to the start of M saves n Toffolis (see tools/dqi_pp_pairs.py). F == M: cancels in place.
Usage: python dqi_pp_pairs2.py <circuit.kmb> <out.npz> [runs=64] [seed=5]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_prove as P  # noqa: E402
import dqi_pp_fp as F  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def smix(x):
    x = x + np.uint64(0x9E3779B97F4A7C15)
    x = (x ^ (x >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
    return x ^ (x >> np.uint64(31))


@njit(cache=True)
def sim(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, isfree, key, dens, count_dens):
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
            if isfree[oi]:
                key[k] = smix(key[k] ^ p)
                if count_dens:
                    dens[k] += F.popc64(p)
                k += 1
            v[t] ^= p & m
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
    return k


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 64
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    ops, W0, W1, W2, BT, CB, DEP = P.decode(c)
    n = len(ops)
    base = ops & 0x7F
    isfree = (base == 0x31) & (CB < 0) & (DEP == 0)
    if os.environ.get('DQI_PP_ALLCCX'):                 # measurement: in-block CCX too
        isfree = (base == 0x31)
    fo = np.flatnonzero(isfree)
    nf = len(fo)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = F.v.ProblemSpecDQI.benchmark_instance(num_shots=F.v.DEFAULT_NUM_SHOTS)
    tables = F.v.FieldTables(F.km.GF2Field(spec.degree))
    wprob = F.v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    key = np.zeros(nf, np.uint64)
    dens = np.zeros(nf, np.int64)
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
        k = sim(ops, *arr, vinit, binit, seed * 1000 + r, isfree, key, dens, r % 2 == 1)
        assert k == nf
    print('fingerprints: %d free CCX, %d runs (%.1fs)' % (nf, R, time.time() - t0), flush=True)
    # events and XOR segments
    Wm = np.stack([W0, W1, W2], axis=1)
    ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
    ev_w = Wm.reshape(-1).astype(np.int64)
    ev_slot = np.tile(np.arange(3), n)
    kk = ev_w >= 0
    ev_op, ev_w, ev_slot = ev_op[kk], ev_w[kk], ev_slot[kk]
    b_ = base[ev_op]
    xorw = ((b_ == 0x11) & (ev_slot == 0)) | ((b_ == 0x21) & (ev_slot == 1)) | ((b_ == 0x31) & (ev_slot == 2))
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, xorw, ev_slot = ev_op[o], ev_w[o], xorw[o], ev_slot[o]
    ne = len(ev_op)
    newseg = np.ones(ne, bool)
    newseg[1:] = (ev_w[1:] != ev_w[:-1]) | ~xorw[1:] | ~xorw[:-1]
    segid = np.cumsum(newseg) - 1
    seg_start = ev_op[newseg]
    tsel = np.flatnonzero(isfree[ev_op] & (ev_slot == 2))
    ccx_seg = np.full(n, -1, np.int64)
    ccx_seg[ev_op[tsel]] = segid[tsel]
    sg = ccx_seg[fo]
    print('events %d, XOR segments %d (%.1fs)' % (ne, int(newseg.sum()), time.time() - t0), flush=True)
    # pairs by key
    order = np.lexsort((fo, key))
    same = key[order][1:] == key[order][:-1]
    pi, pj = order[:-1][same], order[1:][same]
    half = (R // 2) * 64
    dense = (dens[pi] >= 0.25 * half) & (dens[pi] <= 0.75 * half)
    fi, fj, si, sj = fo[pi], fo[pj], sg[pi], sg[pj]
    cancel = si == sj
    print('fingerprint pairs %d (dense %d); cancelling in one XOR segment %d' % (len(pi), int(dense.sum()),
          int(cancel.sum())), flush=True)
    keyg = si * np.int64(1 << 32) + sj
    keyu, inv, cnt = np.unique(keyg[~cancel], return_inverse=True, return_counts=True)
    hist = np.bincount(np.minimum(cnt, 40))
    print('groups %d; size histogram 1..40+: %s' % (len(keyu), hist[1:].tolist()))
    print('pairs in groups of size >= 2: %d, >= 4: %d, >= 8: %d' % (int(cnt[cnt >= 2].sum()), int(cnt[cnt >= 4].sum()),
          int(cnt[cnt >= 8].sum())))
    gfirst = np.full(len(keyu), np.iinfo(np.int64).max)
    np.minimum.at(gfirst, inv, fi[~cancel])
    gms = seg_start[(keyu & ((1 << 32) - 1)).astype(np.int64)]
    np.savez(out, fo=fo, key=key, dens=dens, pi=fi, pj=fj, si=si, sj=sj, cancel=cancel, dense=dense, keyu=keyu,
             cnt=cnt, gfirst=gfirst, gms=gms, seg_start=seg_start, runs=R)
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
