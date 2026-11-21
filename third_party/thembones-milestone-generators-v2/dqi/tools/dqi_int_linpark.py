#!/usr/bin/env python3
"""dqi_int_linpark.py -- B5 (16th holder, Nov 16): FREE parks by linear dependency? Read only.

For each named run [ra, rb) of a decoder file: the idle wires I (no touch inside the run, live through it, as
tools/dqi_int_runwho.py). If one idle wire w is an affine function of other idle wires (w = c + sum of u_j), then
CNOTs (and an X) zero w before the run and the same CNOTs restore it after: a park that costs 0 Toffolis.
Test on samples: the values of every wire of I just before ra on N shots (SCAN's simulate_snapop, the verifier's law
and uniform weights alternately), Gaussian elimination over GF(2) with a constant row; every dependency found is a
candidate (true ones always show; a false one survives N shots with probability about 2^-(N - |I|)).
Prints per run: |I|, rank, the number of dependencies and the smallest few (support size, wires).
Usage: python dqi_int_linpark.py <dqi.kmb> <zero.npz> <ra:rb,...> [R=48 runs of 64 shots] [seed=5]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km  # noqa: E402
from r46_dqi_scan_rank import simulate_snapop  # noqa: E402

SH = 27


def idle_wires(c, ops, zw, zs, ze, tk, ra, rb):
    n = len(ops)
    nq = c.num_qubits
    W = np.arange(nq, dtype=np.int64)
    p = np.searchsorted(tk, (W << SH) | ra, 'left')
    pc = np.minimum(p, len(tk) - 1)
    same = (p < len(tk)) & ((tk[pc] >> SH) == W)
    nxt = np.where(same, tk[pc] & ((1 << SH) - 1), n)
    idle = nxt >= rb
    zke = (zw << SH) | ze
    q = np.searchsorted(zke, (W << SH) | ra, 'right')
    qc = np.minimum(q, len(zke) - 1)
    hit = (q < len(zke)) & (zw[qc] == W) & (zs[qc] < rb)
    idle &= ~hit
    return np.flatnonzero(idle)


def deps(rows):
    """rows: list of python ints (bit vectors). Returns the list of dependencies as sets of row indices
    (each XORs to 0), from a reduced echelon pass."""
    piv = {}            # pivot bit -> (vector, combination as python int over row indices)
    out = []
    for i, r in enumerate(rows):
        comb = 1 << i
        x = r
        while x:
            hb = x.bit_length() - 1
            if hb in piv:
                pv, pc = piv[hb]
                x ^= pv
                comb ^= pc
            else:
                piv[hb] = (x, comb)
                break
        if x == 0:
            out.append(comb)
    return out


def main():
    t0 = time.time()
    path, zp, runs = sys.argv[1:4]
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 48
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 5
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    n = len(ops)
    nq = c.num_qubits
    M = kmb.ID_MASK
    base = ops & 0x7F
    k = ((ops & 0x30) >> 4).astype(np.int8)
    Z = np.load(zp)
    zw = Z['wire'].astype(np.int64)
    zs = np.clip(Z['start'].astype(np.int64), 0, n)
    ze = np.clip(Z['end'].astype(np.int64), 0, n)
    ok = ze > zs
    zw, zs, ze = zw[ok], zs[ok], ze[ok]
    o = np.lexsort((ze, zw))
    zw, zs, ze = zw[o], zs[o], ze[o]
    real = (ops != 0x81) & (base != 0x17) & (base != 0x07)
    i1 = np.cumsum(k == 1, dtype=np.int64) - 1
    i2 = np.cumsum(k == 2, dtype=np.int64) - 1
    i3 = np.cumsum(k == 3, dtype=np.int64) - 1
    parts = []
    for cls, arrs, ci in ((1, (c.q1,), i1), (2, c.q2, i2), (3, c.q3, i3)):
        idx = np.flatnonzero((k == cls) & real)
        for a in arrs:
            w = np.asarray(a).astype(np.int64)[ci[idx]] & M
            parts.append((w << SH) | idx)
    tk = np.sort(np.concatenate(parts))
    del parts
    RUNS = [tuple(int(x) for x in rr.split(':')) for rr in runs.split(',')]
    I = [idle_wires(c, ops, zw, zs, ze, tk, ra, rb) for ra, rb in RUNS]
    print('%s: %d runs; idle sizes %s (%.0fs)' % (path, len(RUNS), [len(x) for x in I], time.time() - t0), flush=True)
    sites = np.unique(np.array([ra for ra, rb in RUNS], np.int64))
    snap = np.zeros((len(sites), nq), np.uint64)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    law = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    uni = np.full(spec.max_errors + 1, 1.0 / (spec.max_errors + 1))
    rng = np.random.default_rng(seed)
    vals = [[0] * len(x) for x in I]     # per run, per idle wire: python int of N bits
    for r in range(R):
        vinit = np.zeros(nq, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, law if r % 2 == 0 else uni)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        simulate_snapop(ops, *arr, vinit, binit, seed * 1000 + r, sites, snap)
        for j, (ra, rb) in enumerate(RUNS):
            row = snap[int(np.searchsorted(sites, ra))]
            for i, w in enumerate(I[j].tolist()):
                vals[j][i] |= int(row[w]) << (64 * r)
    N = 64 * R
    for j, (ra, rb) in enumerate(RUNS):
        rows = vals[j] + [(1 << N) - 1]          # the constant row last
        D = deps(rows)
        sup = sorted((bin(d).count('1'), d) for d in D)
        print('\nRUN %d..%d: idle %d, samples %d, dependencies %d' % (ra, rb, len(I[j]), N, len(D)), flush=True)
        for s, d in sup[:8]:
            mem = [i for i in range(len(rows)) if d >> i & 1]
            names = [('1' if i == len(rows) - 1 else str(int(I[j][i]))) for i in mem]
            print('  support %d: %s' % (s, ' + '.join(names[:16]) + (' ...' if len(names) > 16 else '')))
    print('runtime %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
