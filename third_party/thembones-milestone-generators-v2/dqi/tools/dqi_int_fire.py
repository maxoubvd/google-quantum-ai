#!/usr/bin/env python3
"""dqi_int_fire.py -- B5 (16th holder, Nov 16): which Toffoli-class ops of a decoder file EXECUTE (active lane) but
never CHANGE a value on N samples? Read only. (SCAN's value simulator, 64 lanes per run, the verifier's law and uniform
weights alternately; HMR outcomes random; phases not tracked.) A CCX "fires" when a AND b AND active is nonzero in a
lane; a CSWAP when a AND (x XOR y) AND active is. CCZ ops are phase ops: listed apart (executed count only).
Prints: executed Toffoli-class ops per shot (mean), and the ops that executed but never fired, with their mean
executed count (= what deleting each would save if it is a true no-op), grouped by stage.
Usage: python dqi_int_fire.py <dqi.kmb> [R=64] [seed=31] [out.npz]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km  # noqa: E402
from r46_dqi_scan_rank import _rng  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


@njit(cache=True)
def popc(x):
    c = 0
    while x:
        x &= x - np.uint64(1)
        c += 1
    return c


@njit(cache=True)
def simulate_fire(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, ex, fi):
    vv = vinit.copy()
    bits = binit.copy()
    s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
    stack = np.zeros(256, np.uint64)
    sp = 0
    active = np.uint64(0xFFFFFFFFFFFFFFFF)
    i3 = 0
    i2 = 0
    i1 = 0
    ib = 0
    ic = 0
    kt = 0
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
        if base == 0x30 or base == 0x31 or base == 0x32:
            ex[kt] += popc(m)
            if base == 0x31:
                f = vv[a] & vv[b] & m
                if f != 0:
                    fi[kt] = 1
            elif base == 0x32:
                f = (vv[b] ^ vv[t]) & vv[a] & m
                if f != 0:
                    fi[kt] = 1
            else:
                fi[kt] = 2
            kt += 1
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
            vv[t] ^= vv[a] & vv[b] & m
        elif base == 0x32:
            d = (vv[b] ^ vv[t]) & vv[a] & m
            vv[b] ^= d
            vv[t] ^= d
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
    return kt


def main():
    t0 = time.time()
    path = sys.argv[1]
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 31
    out = sys.argv[4] if len(sys.argv) > 4 else None
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    nq = c.num_qubits
    tofmask = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    tpos = np.flatnonzero(tofmask)
    nt = len(tpos)
    ex = np.zeros(nt, np.int64)
    fi = np.zeros(nt, np.int8)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    law = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    uni = np.full(spec.max_errors + 1, 1.0 / (spec.max_errors + 1))
    rng = np.random.default_rng(seed)
    for r in range(R):
        vinit = np.zeros(nq, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, law if r % 2 == 0 else uni)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        kt = simulate_fire(ops, *arr, vinit, binit, seed * 1000 + r, ex, fi)
        assert kt == nt
    N = 64 * R
    mean_ex = ex / N
    never = (fi == 0) & (ex > 0)
    print('%s: Toffoli-class ops %d, executed per shot (mean over %d samples, law/uniform mix) %.1f (%.0fs)' % (
        path, nt, N, mean_ex.sum(), time.time() - t0))
    print('executed but NEVER fired: %d ops, %.1f executed per shot; CCZ %d ops (%.1f per shot)' % (
        int(never.sum()), mean_ex[never].sum(), int((fi == 2).sum()), mean_ex[fi == 2].sum()))
    st = stages(c)
    for nm, a, b in st:
        sel = never & (tpos >= a) & (tpos <= b)
        if sel.any():
            print('  %-20s %6d ops, %8.1f per shot' % (nm, int(sel.sum()), mean_ex[sel].sum()))
    if out:
        np.savez(out, pos=tpos, ex=ex, fi=fi, N=N)
    segp = os.environ.get('B5_FIRE_SEGS')               # marker builds: never-fired ops by kernel key
    if segp:
        import collections
        import json
        S = json.load(open(segp))
        stacks, segs = S['stacks'], S['segments']
        base = ops & 0x7F
        mark = np.flatnonzero(base == 0x07)
        js = np.searchsorted(mark, tpos[never], 'right') - 1
        cnt = collections.Counter()
        for j, e in zip(js.tolist(), mean_ex[never].tolist()):
            if 0 <= j < len(segs):
                names = [x.split('@')[0] for x in stacks[segs[j][0]]]
                dd = [x for i, x in enumerate(names) if i == 0 or names[i - 1] != x]
                cnt[' > '.join(dd[-3:])] += e
        for key, e in cnt.most_common(30):
            print('  %8.1f  %s' % (e, key))


if __name__ == '__main__':
    main()
