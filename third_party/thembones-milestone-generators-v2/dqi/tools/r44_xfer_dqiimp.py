#!/usr/bin/env python3
"""r44_xfer_dqiimp.py -- R44-XFER: the two-wire affine relations of every Toffoli-class op of a DQI decoder on valid
inputs (the invsqrt affine pass restricted to the forms the stride-200 census found: all but one of its relations
have support 1 or 3).

Per static op, over R runs of 64 valid verifier samples (P1's model), on the samples where the op executes:
  col 0 act, col 1 ones (product 1),
  CCX a b t / CSWAP a (x = b + t):  col 2 = a and not x (violations of "a implies x": product == a),
                                    col 3 = x and not a (product == x), col 4 = neither (a or x always: product ==
                                    1 + a + x)
  CCZ a b t: cols 2..7 = a&~b, b&~a, a&~t, t&~a, b&~t, t&~b (an implication turns the CCZ into a CZ)
With a candidate file (npz: k, mode) it also counts the samples on which ANY candidate relation is violated:
mode 1 product never 1, 2 product == a, 3 product == x, 4 product == 1 + a + x, 12..17 the CCZ implications.
Usage: python r44_xfer_dqiimp.py <dqi.kmb> <out.npz> [runs=256] [seed=5] [candidates.npz]
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
def simulate_imp(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, cnt, mode):
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
    w = np.zeros(8, np.uint64)
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
        elif base == 0x31 or base == 0x32:
            va = vv[a]
            x = vv[b]
            if base == 0x32:
                x = vv[b] ^ vv[t]
            d = va & x & m
            w[1] = d
            w[2] = va & ~x & m
            w[3] = x & ~va & m
            w[4] = ~va & ~x & m
            cnt[k, 0] += popc64(m)
            for j in range(1, 5):
                cnt[k, j] += popc64(w[j])
            if mode[k] > 0:
                anyf |= w[mode[k]]
            if base == 0x31:
                vv[t] ^= d
            else:
                vv[b] ^= d
                vv[t] ^= d
            k += 1
        elif base == 0x30:
            va = vv[a]
            vb = vv[b]
            vt = vv[t]
            w[1] = va & vb & vt & m
            w[2] = va & ~vb & m
            w[3] = vb & ~va & m
            w[4] = va & ~vt & m
            w[5] = vt & ~va & m
            w[6] = vb & ~vt & m
            w[7] = vt & ~vb & m
            cnt[k, 0] += popc64(m)
            for j in range(1, 8):
                cnt[k, j] += popc64(w[j])
            if mode[k] == 1:
                anyf |= w[1]
            elif mode[k] >= 12:
                anyf |= w[mode[k] - 10]
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
    mode = np.zeros(nt, np.int64)
    if len(sys.argv) > 5:
        d = np.load(sys.argv[5])
        mode[d["k"]] = d["mode"]
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    rng = np.random.default_rng(seed)
    cnt = np.zeros((nt, 8), np.int64)
    fired = 0
    print("read %.0fs: %d ops, %d Toffoli-class static, %d candidates" % (time.time() - t0, len(ops), nt,
                                                                          int((mode > 0).sum())), flush=True)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        k, anyf = simulate_imp(ops, *arr, vinit, binit, seed * 1000 + r, cnt, mode)
        assert k == nt
        fired += bin(int(anyf)).count("1")
        if (r + 1) % 64 == 0 or r == R - 1:
            n = 64 * (r + 1)
            print("run %d: %d samples, mean executed %.1f, samples with a candidate relation violated %d (%.0fs)"
                  % (r + 1, n, cnt[:, 0].sum() / n, fired, time.time() - t0), flush=True)
            np.savez(out + ".tmp.npz", t_op=t_op, base=(ops[t_op] & 0x7F), raw=ops[t_op], cnt=cnt, samples=n,
                     fired=fired)
            os.replace(out + ".tmp.npz", out)
    print("done %.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
