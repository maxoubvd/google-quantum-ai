#!/usr/bin/env python3
"""dqi_pp_deadw.py -- workstream dqi_pp (P1, sixth holder): the DEAD-WRITE pass (public rewrite c3), exact by structure.

Site: an unconditioned CCX a b t at block depth 0 whose target t is next touched only by XOR-writes (X t, CX . t,
CCX . . t: t is never read) and then by an unconditioned HMR t at depth 0, with a and b not modified between the CCX
and the HMR (input wires included: the public scans skipped them).  The HMR's phase kick is m * (value of t), and the
CCX adds a*b to that value only.  Rewrite: delete the CCX and put CZ a b conditioned on the HMR's bit right after
the HMR: the kick m*a*b is paid by the CZ, every wire value, outcome distribution and later correction is unchanged.
One Toffoli less per site on every shot.  (Several sites may share one HMR; each gets its own CZ.)
Usage: python dqi_pp_deadw.py <in.kmb> <out.kmb>
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402


class P:
    """Per-op decode (copy of tools/dqi_pp_prove.decode, numpy only: this pass needs no numba)."""
    @staticmethod
    def decode(c):
        ops = np.asarray(c.ops)
        n = len(ops)
        nq = (ops & 0x30) >> 4
        W = np.full((3, n), -1, dtype=np.int32)
        for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
            sel = np.flatnonzero(nq == k)
            for j in range(k):
                W[j, sel] = np.asarray(arrs[j], dtype=np.int32)
        base = ops & 0x7F
        hb = ((ops & 0x40) != 0) & (base != 0x06) & (base != 0x07)
        BT = np.full(n, -1, dtype=np.int32)
        BT[hb] = np.asarray(c.bit_t, dtype=np.int32)
        hc = (ops & 0x80) != 0
        CB = np.full(n, -1, dtype=np.int32)
        CB[hc] = np.asarray(c.cond, dtype=np.int32)
        push = ops == 0x81
        pop = ops == 0x01
        DEP = (np.cumsum(push.astype(np.int64) - pop.astype(np.int64)) - push).astype(np.int8)
        return np.ascontiguousarray(ops), W[0].copy(), W[1].copy(), W[2].copy(), BT, CB, DEP


def sites(c):
    """List of (ccx op, hmr op, a, b, hmr bit) for every exact dead-write site."""
    ops, W0, W1, W2, BT, CB, DEP = P.decode(c)
    n = len(ops)
    base = ops & 0x7F
    Wm = np.stack([W0, W1, W2], axis=1)
    ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
    ev_w = Wm.reshape(-1).astype(np.int64)
    ev_slot = np.tile(np.arange(3, dtype=np.int8), n)
    kk = ev_w >= 0
    ev_op, ev_w, ev_slot = ev_op[kk], ev_w[kk], ev_slot[kk]
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, ev_slot = ev_op[o], ev_w[o], ev_slot[o]
    b_ = base[ev_op]
    xorw = ((b_ == 0x11) & (ev_slot == 0)) | ((b_ == 0x21) & (ev_slot == 1)) | ((b_ == 0x31) & (ev_slot == 2))
    ne = len(ev_op)
    key = ev_w * np.int64(1 << 40) + ev_op
    cc = np.flatnonzero((base == 0x31) & (CB < 0) & (DEP == 0))
    et = np.searchsorted(key, W2[cc].astype(np.int64) * np.int64(1 << 40) + cc)
    notx = np.flatnonzero(~xorw)
    q = np.searchsorted(notx, et + 1)
    nx = np.where(q < len(notx), notx[np.minimum(q, len(notx) - 1)], ne)
    nxc = np.minimum(nx, ne - 1)
    rop = ev_op[nxc]
    ok = (nx < ne) & (ev_w[nxc] == W2[cc]) & (base[rop] == 0x50) & (CB[rop] < 0) & (DEP[rop] == 0)
    mod = xorw | np.isin(b_, (0x12, 0x50, 0x22)) | ((b_ == 0x32) & (ev_slot >= 1))
    mkey = key[mod]

    def modified(wire, lo, hi):
        i1 = np.searchsorted(mkey, wire * np.int64(1 << 40) + lo, side='right')
        i2 = np.searchsorted(mkey, wire * np.int64(1 << 40) + hi, side='left')
        return i2 > i1
    s = np.flatnonzero(ok)
    a, b = W0[cc[s]].astype(np.int64), W1[cc[s]].astype(np.int64)
    clean = ~modified(a, cc[s], rop[s]) & ~modified(b, cc[s], rop[s])
    s = s[clean]
    out = []
    for j in s:
        i, r = int(cc[j]), int(rop[j])
        # independent re-check by a direct scan of the ops between (slices; no event bookkeeping)
        t, aa, bb = int(W2[i]), int(W0[i]), int(W1[i])
        sl = slice(i + 1, r)
        bs, w0, w1, w2 = base[sl], W0[sl], W1[sl], W2[sl]
        reads_t = ((w0 == t) & ~(bs == 0x11)) | ((w1 == t) & ~(bs == 0x21)) | ((w2 == t) & ~(bs == 0x31))
        mod_ab = np.zeros(len(bs), bool)
        for w in (aa, bb):
            mod_ab |= ((w0 == w) & np.isin(bs, (0x11, 0x12, 0x50, 0x22))) | ((w1 == w) & np.isin(bs, (0x21, 0x22, 0x32))) \
                | ((w2 == w) & np.isin(bs, (0x31, 0x32)))
        assert not reads_t.any() and not mod_ab.any(), ('site fails the direct check', i, r)
        assert W0[r] == t and BT[r] >= 0
        out.append((i, r, aa, bb, int(BT[r])))
    return out


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    S = sites(c)
    print('dead-write sites: %d (%.0fs)' % (len(S), time.time() - t0), flush=True)
    delete = [s[0] for s in S]
    inserts = {}
    for (i, r, a, b, bit) in S:
        inserts.setdefault(r, []).append((0xA0, (a, b), -1, bit))      # CZ a b, conditioned on the HMR's bit
    k1 = kw.edit(c, sorted(delete), inserts)
    kw.write(sys.argv[2], k1)
    print('written %s: -%d CCX, +%d conditioned CZ (%.0fs)' % (sys.argv[2], len(delete), len(delete), time.time() - t0))
    with open(sys.argv[2] + '.sites.txt', 'w') as f:
        f.write('ccx_op hmr_op a b bit\n')
        for s in S:
            f.write('%d %d %d %d %d\n' % s)


if __name__ == '__main__':
    main()
