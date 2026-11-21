#!/usr/bin/env python3
"""dqi_int_ppzero.py -- workstream dqi_int (B5): idle-zero intervals of every wire of an emitted .kmb (no rewrite).

A wire is idle between two consecutive touches; it is ZERO there when its value right after the earlier touch is 0.
Structural zeros (after HMR / R, before the first touch) are taken as proved; the others are accepted when the value
is 0 in EVERY sample of R runs (valid inputs by the law, valid with weight uniform on 18..35, random inputs).
Output npz: wire, start (earlier touch op, -1 at the start), end (next touch op, n at the end), structural flag.
Usage: python dqi_int_ppzero.py <circuit.kmb> <out.npz> [runs=192] [seed=31] [safe]
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
SAFE = len(sys.argv) > 5 and sys.argv[5] == 'safe'   # argv[5] = 'safe': P1's two audit fixes (default off: the
                                                      # published files were made without them)


@njit(cache=True)
def zeros_after(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, nz):
    """nz[3 * op + slot] &= (wire of that operand is 0 in all 64 samples right after the op)."""
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
            s = _rng(s)
            bits[bit] = (bits[bit] & ~m) | (s & m)
            v[a] &= ~m
        elif base == 0x40:
            bits[bit] &= ~m
        elif base == 0x41:
            bits[bit] |= m
        elif base == 0x42:
            bits[bit] ^= m
        if nqk >= 1:
            nz[3 * oi] &= v[a] == 0
        if nqk >= 2:
            nz[3 * oi + 1] &= v[b] == 0
        if nqk == 3:
            nz[3 * oi + 2] &= v[t] == 0


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 192
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 31
    ops = np.ascontiguousarray(np.asarray(c.ops))
    n = len(ops)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wlaw = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    wuni = np.zeros_like(wlaw)
    wuni[18:36] = 1.0
    wuni /= wuni.sum()
    rng = np.random.default_rng(seed)
    nz = np.ones(3 * n, dtype=np.bool_)
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
        zeros_after(ops, *arr, vinit, binit, seed * 1000 + r, nz)
        if r % 32 == 31:
            print('run %d %.1fs' % (r, time.time() - t0), flush=True)
    base = ops & 0x7F
    nq = (ops & 0x30) >> 4
    W = np.full((n, 3), -1, dtype=np.int64)
    for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
        sel = np.flatnonzero(nq == k)
        for j in range(k):
            W[sel, j] = np.asarray(arrs[j], dtype=np.int64)
    ev = np.flatnonzero(W.reshape(-1) >= 0)
    ev_op, ev_w = ev // 3, W.reshape(-1)[ev]
    reset = np.isin(base, (0x50, 0x12))
    if SAFE:
        # P1's audit (Nov 13): a reset is a structural zero only when it is unconditioned and outside every
        # PUSH_CONDITION block; the others must be zero on the samples like any other interval
        depth = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
        depth_before = np.concatenate([[0], depth[:-1]])
        reset = reset & ((ops & 0x80) == 0) & (depth_before == 0)
    zafter = nz[ev] | reset[ev_op]
    struct_after = reset[ev_op]
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, zafter, struct_after = ev_op[o], ev_w[o], zafter[o], struct_after[o]
    same = np.zeros(len(ev_op), dtype=bool)
    same[:-1] = ev_w[:-1] == ev_w[1:]
    nxt = np.where(same, np.concatenate([ev_op[1:], [n]]), n)
    selz = zafter
    wires = np.concatenate([np.arange(c.num_qubits), ev_w[selz]])
    firsts = np.searchsorted(ev_w, np.arange(c.num_qubits), side='left')
    first_op = np.where(firsts < len(ev_w), ev_op[np.minimum(firsts, len(ev_w) - 1)], n)
    first_op = np.where((firsts < len(ev_w)) & (ev_w[np.minimum(firsts, len(ev_w) - 1)] == np.arange(c.num_qubits)),
                        first_op, n)
    starts = np.concatenate([np.full(c.num_qubits, -1), ev_op[selz]])
    ends = np.concatenate([first_op, nxt[selz]])
    structural = np.concatenate([np.ones(c.num_qubits, dtype=bool), struct_after[selz]])
    if SAFE:   # the 840 syndrome wires hold the input before their first touch: no interval there
        keep = ~((starts == -1) & (wires < 840))
        wires, starts, ends, structural = wires[keep], starts[keep], ends[keep], structural[keep]
    np.savez(out, wire=wires, start=starts, end=ends, structural=structural)
    print('idle-zero intervals %d (structural %d, by samples %d) %.1fs' % (
        len(wires), int(structural.sum()), int((~structural).sum()), time.time() - t0))


if __name__ == '__main__':
    main()
