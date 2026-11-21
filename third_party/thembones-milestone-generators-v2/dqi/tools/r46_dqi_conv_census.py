#!/usr/bin/env python3
"""r46_dqi_conv_census.py -- R46-DQI-CONV: per static Toffoli-class op of a DQI decoder, on valid verifier samples
(P1's bit-parallel value model, tools/dqi_pp_act.py; a copy of R44-XFER's tools/r44_xfer_dqiones.py with the two
controls counted separately and the wires stored):
  act = samples on which the op executes; ones = executes AND its product is 1;
  ca  = executes AND control a is 1; cb = executes AND the second factor is 1 (CCX: wire b; CSWAP: b xor t;
  CCZ: b and c). Wires a, b, t of every op are stored, so a rare control wire (a flag) can be read off per op.
Optional watch list (npy of static indices k): `fired` = samples on which any listed op has product 1.
Usage: python r46_dqi_conv_census.py <dqi.kmb> <out.npz> [runs=256] [seed=5] [watch.npy]
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
def simulate(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, act, ones, ca, cb, cand):
    vv = vinit.copy()
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
    anyf = np.uint64(0)
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
            act[k] += popc64(m)
            ones[k] += popc64(d)
            ca[k] += popc64(vv[a] & m)
            cb[k] += popc64(vv[b] & m)
            if cand[k]:
                anyf |= d
            vv[t] ^= d
            k += 1
        elif base == 0x32:
            x = (vv[b] ^ vv[t]) & m
            d = x & vv[a]
            act[k] += popc64(m)
            ones[k] += popc64(d)
            ca[k] += popc64(vv[a] & m)
            cb[k] += popc64(x)
            if cand[k]:
                anyf |= d
            vv[b] ^= d
            vv[t] ^= d
            k += 1
        elif base == 0x30:
            x = vv[b] & vv[t] & m
            d = vv[a] & x
            act[k] += popc64(m)
            ones[k] += popc64(d)
            ca[k] += popc64(vv[a] & m)
            cb[k] += popc64(x)
            if cand[k]:
                anyf |= d
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
    return k, anyf


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    R = int(sys.argv[3]) if len(sys.argv) > 3 else 256
    seed = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    ops = np.ascontiguousarray(np.asarray(c.ops))
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    t_op = np.flatnonzero(toff)
    nt = len(t_op)
    cand = np.zeros(nt, np.uint8)
    if len(sys.argv) > 5:
        cand[np.load(sys.argv[5])] = 1
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    # wires of the Toffoli-class ops: the j-th 3-qubit op among all 3-qubit ops
    is3 = ((ops & 0x30) == 0x30) & (ops != 0x81) & ((ops & 0x7F) != 0x01)
    idx3 = np.cumsum(is3) - 1
    j = idx3[t_op]
    wa, wb, wt = arr[0][j].astype(np.int32), arr[1][j].astype(np.int32), arr[2][j].astype(np.int32)
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    act = np.zeros(nt, np.int64)
    ones = np.zeros(nt, np.int64)
    ca = np.zeros(nt, np.int64)
    cb = np.zeros(nt, np.int64)
    fired = 0
    print("read %.0fs: %d ops, %d Toffoli-class static, %d watched" % (time.time() - t0, len(ops), nt, cand.sum()),
          flush=True)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        k, anyf = simulate(ops, *arr, vinit, binit, seed * 1000 + r, act, ones, ca, cb, cand)
        assert k == nt
        fired += bin(int(anyf)).count("1")
        if (r + 1) % 64 == 0 or r == R - 1:
            n = 64 * (r + 1)
            print("run %d: %d samples, mean executed %.1f, watched fired on %d samples (%.0fs)"
                  % (r + 1, n, act.sum() / n, fired, time.time() - t0), flush=True)
            np.savez(out + ".tmp.npz", t_op=t_op, base=(ops[t_op] & 0x7F), act=act, ones=ones, ca=ca, cb=cb,
                     wa=wa, wb=wb, wt=wt, samples=n, fired=fired)
            os.replace(out + ".tmp.npz", out)
    print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
