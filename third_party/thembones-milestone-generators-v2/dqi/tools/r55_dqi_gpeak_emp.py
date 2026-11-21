#!/usr/bin/env python3
"""r55_dqi_gpeak_emp.py -- R55-DQI-GPEAK (designer): EMPIRICAL zero intervals of every wire, to set against P1's
EXACT ones (tools/dqi_pp_zero.py). Read only.

The value simulator is B5's tools/dqi_int_fire.py loop (SCAN's simulator: 64 lanes per run, the verifier's law and
uniform weights alternately, HMR outcomes random, phases not tracked; credit B5 / R46-DQI-SCAN; originals unchanged).
Added: after every op, for each operand wire, whether it is nonzero in some lane (OR over all runs).
A wire is EMPIRICALLY zero between two consecutive events on it when it was zero in every lane of every run after the
first event. P1's exact intervals are a subset (value numbering is sound and incomplete); the difference is
  (a) wires written and erased inside PUSH_CONDITION blocks, where value numbering loses the identity
      K a + K b = K (a + b): they hold nothing, and stay "live" in the exact file until a later unconditional HMR;
  (b) values that are rarely nonzero (gap flags at 2^-12): empirically zero here by sampling only.
Class (a) is separated structurally: every write to the wire since its last exact-zero moment is a conditional op
(inside a PUSH block or with the _IF bit).
Output npz: wire, start, end (the same fields as P1's zero.npz: free at after-op times start .. end-1), cls (1 = also
exact, 2 = conditional-only writes, 3 = other).
Usage: python r55_dqi_gpeak_emp.py <dqi.kmb> <exact zero.npz> <out.npz> [R=32] [seed=31]
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

SH = 27


@njit(cache=True)
def simulate_nz(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, nz, cnd):
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
        if sp > 0 or (op & 0x80) != 0:
            cnd[oi] = 1
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
        if nqk >= 1:
            z = nz[oi]
            if vv[a] != 0:
                z |= 1
            if nqk >= 2 and vv[b] != 0:
                z |= 2
            if nqk == 3 and vv[t] != 0:
                z |= 4
            nz[oi] = z
    return 0


@njit(cache=True)
def cond_only(ev_w, exact, ev_wr, condw):
    """Per event (sorted by wire, op): True while every write to the wire since its last exact-zero event (or since
    the start, for a non-input wire) was a conditional op."""
    ne = ev_w.shape[0]
    allc = np.zeros(ne, np.bool_)
    cur_w = -1
    flag = True
    for i in range(ne):
        if ev_w[i] != cur_w:
            cur_w = ev_w[i]
            flag = cur_w >= 840
        if exact[i]:
            flag = True
        elif ev_wr[i] and not condw[i]:
            flag = False
        allc[i] = flag
    return allc


def main():
    t0 = time.time()
    path, zp, out = sys.argv[1:4]
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 32
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 31
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    n = len(ops)
    nq = c.num_qubits
    nz = np.zeros(n, np.uint8)
    cnd = np.zeros(n, np.uint8)
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
        simulate_nz(ops, *arr, vinit, binit, seed * 1000 + r, nz, cnd)
        if r % 8 == 7 or r == R - 1:
            print('run %d of %d (%.0fs)' % (r + 1, R, time.time() - t0), flush=True)
    # events per wire
    M = kmb.ID_MASK
    base = ops & 0x7F
    k = ((ops & 0x30) >> 4).astype(np.int8)
    real = (ops != 0x81)
    i1 = np.cumsum(k == 1, dtype=np.int64) - 1
    i2 = np.cumsum(k == 2, dtype=np.int64) - 1
    i3 = np.cumsum(k == 3, dtype=np.int64) - 1
    q1 = np.asarray(c.q1).astype(np.int64) & M
    q2 = [np.asarray(x).astype(np.int64) & M for x in c.q2]
    q3 = [np.asarray(x).astype(np.int64) & M for x in c.q3]
    ev_w, ev_op, ev_z, ev_wr = [], [], [], []
    for cls, arrs, ii in ((1, (q1,), i1), (2, q2, i2), (3, q3, i3)):
        idx = np.flatnonzero((k == cls) & real)
        b_ = base[idx]
        for slot, a in enumerate(arrs):
            ev_w.append(a[ii[idx]])
            ev_op.append(idx)
            ev_z.append(((nz[idx] >> slot) & 1) == 0)
            if cls == 3:
                wr = ((b_ == 0x31) & (slot == 2)) | ((b_ == 0x32) & (slot >= 1))
            elif cls == 2:
                wr = ((b_ == 0x21) & (slot == 1)) | (b_ == 0x22)
            else:
                wr = (b_ == 0x11) | (b_ == 0x12) | (b_ == 0x50)
            ev_wr.append(wr)
    ev_w = np.concatenate(ev_w)
    ev_op = np.concatenate(ev_op)
    ev_z = np.concatenate(ev_z)
    ev_wr = np.concatenate(ev_wr)
    o = np.lexsort((ev_op, ev_w))
    ev_w, ev_op, ev_z, ev_wr = ev_w[o], ev_op[o], ev_z[o], ev_wr[o]
    ne = len(ev_op)
    same = ev_w[1:] == ev_w[:-1]
    nxt = np.full(ne, n, np.int64)
    nxt[:-1][same] = ev_op[1:][same]
    # exact zero: is the wire exactly zero after event e? (inside an exact interval [s, e) with s <= op < e)
    Z = np.load(zp)
    zw = Z['wire'].astype(np.int64)
    zs = np.clip(Z['start'].astype(np.int64), 0, n)
    ze = np.clip(Z['end'].astype(np.int64), 0, n)
    okz = ze > zs
    zw, zs, ze = zw[okz], zs[okz], ze[okz]
    oz = np.lexsort((zs, zw))
    zw, zs, ze = zw[oz], zs[oz], ze[oz]
    zks = (zw << SH) | zs
    p = np.searchsorted(zks, (ev_w << SH) | ev_op, 'right') - 1
    pc = np.maximum(p, 0)
    exact = (p >= 0) & (zw[pc] == ev_w) & (ze[pc] > ev_op)
    # conditional-only writes since the last exact-zero event: running flag per wire
    condw = cnd[ev_op].astype(bool)
    allc = cond_only(ev_w, exact, ev_wr, condw)
    sel = np.flatnonzero(ev_z)
    cls_ = np.where(exact[sel], 1, np.where(allc[sel], 2, 3)).astype(np.int8)
    wire, start, end = ev_w[sel], ev_op[sel], nxt[sel]
    # untouched non-input wires and the span before a non-input wire's first event
    first = np.ones(ne, bool)
    first[1:] = ~same
    fw, fo = ev_w[first], ev_op[first]
    mfw = fw >= 840
    touched = np.zeros(nq, bool)
    touched[fw] = True
    never = np.flatnonzero(~touched & (np.arange(nq) >= 840))
    wire = np.concatenate([wire, fw[mfw], never])
    start = np.concatenate([start, np.zeros(int(mfw.sum()), np.int64), np.zeros(len(never), np.int64)])
    end = np.concatenate([end, fo[mfw], np.full(len(never), n, np.int64)])
    cls_ = np.concatenate([cls_, np.ones(int(mfw.sum()), np.int8), np.ones(len(never), np.int8)])
    np.savez(out, wire=wire, start=start, end=end, cls=cls_, samples=64 * R)
    # profile: tight live under exact, exact + class 2, and all empirical intervals
    tof = ((ops & 0x7C) == 0x30) & (ops != 0x81)

    def livecount(mask):
        d = np.zeros(n + 2, np.int32)
        s_ = np.clip(start[mask], 0, n)
        e_ = np.clip(end[mask], 0, n)
        okk = e_ > s_
        np.add.at(d, s_[okk], 1)
        np.add.at(d, e_[okk], -1)
        return (nq - np.cumsum(d[:n])).astype(np.int32)
    Lx = livecount(cls_ == 1)
    Lc = livecount(cls_ <= 2)
    Le = livecount(cls_ <= 3)
    print('%s: Q %d, ops %d, samples %d; intervals: exact %d, conditional-only %d, other empirical %d (%.0fs)' % (
        path, nq, n, 64 * R, int((cls_ == 1).sum()), int((cls_ == 2).sum()), int((cls_ == 3).sum()), time.time() - t0))
    levels = (864, 877, 889, 896, 899, 900)
    print('tight peak of the file: exact %d | exact + conditional-only %d | all empirical %d' % (
        int(Lx.max()), int(Lc.max()), int(Le.max())))
    print('%-20s %-22s %6s  %s' % ('stage', 'zero intervals', 'peak', 'Toffoli-class ops at live > ' + ' / '.join(map(str, levels))))
    for nm, a, b in stages(c):
        T = tof[a:b]
        for lab_, L in (('exact (P1)', Lx), ('exact + conditional', Lc), ('all empirical', Le)):
            LL = L[a:b]
            print('%-20s %-22s %6d  %s' % (nm, lab_, int(LL.max()), ' / '.join('%d' % int(T[LL > lv].sum()) for lv in levels)))
    if os.environ.get('R55_EMP_LIVE') == '1':      # 1 GB per decoder file: only on request (box 6's disk is tight)
        np.savez(out.replace('.npz', '_live.npz'), Lx=Lx, Lc=Lc, Le=Le)
    print('runtime %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
