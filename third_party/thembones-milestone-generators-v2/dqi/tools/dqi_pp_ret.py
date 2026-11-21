#!/usr/bin/env python3
"""dqi_pp_ret.py -- workstream dqi_pp (P1): retained products with CHAINS and roving scratch (rewrite).

A product computed m times (free CCX o_1 < ... < o_m, same product function) is held on one zero-idle wire from o_s
to o_e for any run of consecutive occurrences, saving e - s Toffolis:
  o_s: CCX(a, b, t) -> CCX(a, b, w); CX(w, t)     middle o_k: CCX(..., t_k) -> CX(w, t_k)
  o_e: CCX(a', b', t') -> CX(w, t'); HMR(w) -> m; CZ(a', b') if m
w hops between wires that are zero and untouched in the original (SWAP before the original touches it).
Links (o_k, o_k+1) are admitted greedily, shortest first, against the free-wire profile (segment tree); adjacent
admitted links of one product merge into one held copy. Sound products: product function equal on every sample of
192 runs and dense (>= 1/16 of uniform samples), or links whose two occurrences have the same exact form
(tools/dqi_pp_gvn.py). Zero intervals: exact ones (tools/dqi_pp_zero.py / dqi_pp_zprove.py).
Usage: python dqi_pp_ret.py <in.kmb> <pairs2.npz (192 runs)> <gvn.npz> <zero.npz> <out.kmb> [dens_div=16]
"""
import dataclasses
import heapq
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402


@njit(cache=False)   # cached recursive code crashed on reload (Nov 13)
def st_build(vals, size):
    mn = np.full(2 * size, 1 << 40, np.int64)
    lz = np.zeros(2 * size, np.int64)
    mn[size:size + vals.shape[0]] = vals
    for i in range(size - 1, 0, -1):
        mn[i] = min(mn[2 * i], mn[2 * i + 1])
    return mn, lz


@njit(cache=False)   # cached recursive code crashed on reload (Nov 13)
def st_query(mn, lz, size, l, r):
    # iterative range min with lazy adds on the path (top-down recursion emulated by a stack)
    res = 1 << 40
    stack = np.zeros((64, 4), np.int64)        # node, nl, nr, acc
    sp = 0
    stack[0, 0] = 1
    stack[0, 1] = 0
    stack[0, 2] = size - 1
    stack[0, 3] = 0
    sp = 1
    while sp > 0:
        sp -= 1
        node, nl, nr, acc = stack[sp, 0], stack[sp, 1], stack[sp, 2], stack[sp, 3]
        if r < nl or nr < l:
            continue
        if l <= nl and nr <= r:
            v = mn[node] + acc
            if v < res:
                res = v
            continue
        mid = (nl + nr) // 2
        a2 = acc + lz[node]
        stack[sp, 0] = 2 * node
        stack[sp, 1] = nl
        stack[sp, 2] = mid
        stack[sp, 3] = a2
        sp += 1
        stack[sp, 0] = 2 * node + 1
        stack[sp, 1] = mid + 1
        stack[sp, 2] = nr
        stack[sp, 3] = a2
        sp += 1
    return res


@njit(cache=False)   # cached recursive code crashed on reload (Nov 13)
def st_add(mn, lz, node, nl, nr, l, r, v):
    if r < nl or nr < l:
        return
    if l <= nl and nr <= r:
        mn[node] += v
        lz[node] += v
        return
    mid = (nl + nr) // 2
    st_add(mn, lz, 2 * node, nl, mid, l, r, v)
    st_add(mn, lz, 2 * node + 1, mid + 1, nr, l, r, v)
    mn[node] = lz[node] + min(mn[2 * node], mn[2 * node + 1])


@njit(cache=False)   # cached recursive code crashed on reload (Nov 13)
def admit(lk_lo, lk_hi, lk_prod, lk_k, order, prod_off, nocc, pts, mn, lz, size):
    """Greedy admission of links in `order`. A link (product p, index k) spans occurrences o_k .. o_k+1 given as
    after-op ranges: own range [lk_lo, lk_hi] when standalone. Returns admitted flags."""
    nl = lk_lo.shape[0]
    adm = np.zeros(nl, np.bool_)
    # link index of (product, k): prod_off[p] + k
    for q in range(order.shape[0]):
        L = order[q]
        p = lk_prod[L]
        k = lk_k[L]
        base_ = prod_off[p]
        left_adm = k > 0 and adm[base_ + k - 1]
        right_adm = k + 1 < nocc[p] - 1 and adm[base_ + k + 1]
        lo = lk_lo[L]
        hi = lk_hi[L]
        if left_adm:
            lo = lk_lo[L] + 2                  # previous link already holds through o_k
        if right_adm:
            hi = lk_hi[L] - 2                  # next link already holds from o_k+1 - 1
        if hi < lo:
            adm[L] = True
            continue
        a = np.searchsorted(pts, lo, side='right') - 1
        b = np.searchsorted(pts, hi, side='right') - 1
        if st_query(mn, lz, size, a, b) >= 1:
            st_add(mn, lz, 1, 0, size - 1, a, b, -1)
            adm[L] = True
    return adm


def main():
    t0 = time.time()
    src, pairs_npz, gvn_npz, zero_npz, dst = sys.argv[1:6]
    ddiv = int(sys.argv[6]) if len(sys.argv) > 6 else 16
    c = kmb.read(src)
    ops = np.asarray(c.ops)
    n = len(ops)
    nqo = (ops & 0x30) >> 4
    W = np.full((n, 3), -1, dtype=np.int64)
    for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
        sel = np.flatnonzero(nqo == k)
        for j in range(k):
            W[sel, j] = np.asarray(arrs[j], dtype=np.int64)
    Pp = np.load(pairs_npz)
    G = np.load(gvn_npz)
    fo, key, dens = Pp['fo'], Pp['key'], Pp['dens']
    half = (int(Pp['runs']) // 2) * 64
    # exact form per free CCX
    gi = np.searchsorted(G['t_op'], fo)
    assert (G['t_op'][gi] == fo).all()
    e0, e1, ec = G['t_p0'][gi], G['t_p1'][gi], G['t_pc'][gi]
    LINKS = None
    if os.environ.get('DQI_PP_LINKS'):
        Lk = np.load(os.environ['DQI_PP_LINKS'])
        LINKS = set(zip(Lk['i'].tolist(), Lk['j'].tolist()))
        print('proved links offered: %d' % len(LINKS), flush=True)
    # products: group by key, occurrences in time order
    order = np.lexsort((fo, key))
    ks = key[order]
    brk = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1], True])
    lk_lo, lk_hi, lk_prod, lk_k, lk_i, lk_j = [], [], [], [], [], []
    prod_occ = []
    for g_ in range(len(brk) - 1):
        occ = order[brk[g_]:brk[g_ + 1]]
        if len(occ) < 2:
            continue
        dense = dens[occ[0]] >= half // ddiv
        ok_links = []
        for k in range(len(occ) - 1):
            x, y = occ[k], occ[k + 1]
            exact = e0[x] == e0[y] and e1[x] == e1[y] and ec[x] == ec[y]
            if LINKS is not None:                        # only links proved exact (tools/dqi_pp_pprove.py)
                ok_links.append(bool(exact or (int(fo[x]), int(fo[y])) in LINKS))
            else:
                ok_links.append(bool(dense or exact))
        # split into runs of sound links: each run becomes its own "product" chain domain
        k = 0
        while k < len(ok_links):
            if not ok_links[k]:
                k += 1
                continue
            s = k
            while k < len(ok_links) and ok_links[k]:
                k += 1
            seg = occ[s:k + 1]
            p = len(prod_occ)
            prod_occ.append(fo[seg])
            for kk in range(len(seg) - 1):
                i, j = int(fo[seg[kk]]), int(fo[seg[kk + 1]])
                lk_lo.append(i - 1)
                lk_hi.append(j)
                lk_prod.append(p)
                lk_k.append(kk)
                lk_i.append(i)
                lk_j.append(j)
    nccx_prod = len(prod_occ)
    exact_all = os.environ.get('DQI_PP_EXACTALL') == '1'
    if exact_all:                                       # links only from exact forms, CCX and CSWAP together
        prod_occ, lk_lo, lk_hi, lk_prod, lk_k, lk_i, lk_j = [], [], [], [], [], [], []
    if (os.environ.get('DQI_PP_CSWAP', '1') == '1' or exact_all) and 't_p0' in G.files:
        # links: consecutive free CSWAP (and, in exact_all, free CCX) with the same EXACT product form
        # (CCX a*b, CSWAP d = c * (x + y); value numbering): exact by construction
        bt_ = ops[G['t_op']] & 0x7F
        msw = ((bt_ == 0x32) | ((bt_ == 0x31) & exact_all)) & (G['t_lvl'] == 0)
        so = G['t_op'][msw]
        k0, k1, kc = G['t_p0'][msw], G['t_p1'][msw], G['t_pc'][msw]
        srt = np.lexsort((so, kc, k1, k0))
        ks = np.stack([k0[srt], k1[srt], kc[srt].astype(np.uint64)], 1)
        brk2 = np.flatnonzero(np.r_[True, np.any(ks[1:] != ks[:-1], axis=1), True])
        nsw = 0
        for g_ in range(len(brk2) - 1):
            occ = so[srt[brk2[g_]:brk2[g_ + 1]]]
            if len(occ) < 2:
                continue
            p = len(prod_occ)
            prod_occ.append(occ.astype(np.int64))
            for kk in range(len(occ) - 1):
                i, j = int(occ[kk]), int(occ[kk + 1])
                lk_lo.append(i - 1)
                lk_hi.append(j)
                lk_prod.append(p)
                lk_k.append(kk)
                lk_i.append(i)
                lk_j.append(j)
                nsw += 1
        print('CSWAP exact links %d' % nsw, flush=True)
    lk_lo, lk_hi = np.array(lk_lo, np.int64), np.array(lk_hi, np.int64)
    lk_prod, lk_k = np.array(lk_prod, np.int64), np.array(lk_k, np.int64)
    nocc = np.array([len(x) for x in prod_occ], np.int64)
    prod_off = np.r_[0, np.cumsum(nocc - 1)[:-1]].astype(np.int64)
    print('sound products %d, links %d (%.1fs)' % (len(prod_occ), len(lk_lo), time.time() - t0), flush=True)
    # free profile
    Z = np.load(zero_npz)
    zw, zs, ze = Z['wire'].astype(np.int64), Z['start'].astype(np.int64), Z['end'].astype(np.int64)
    kx = int(os.environ.get('DQI_PP_EXTRA_WIRES', '0'))          # measurement only: k always-free extra wires
    if kx:
        zw = np.r_[zw, np.arange(c.num_qubits, c.num_qubits + kx)]
        zs = np.r_[zs, np.full(kx, -1)]
        ze = np.r_[ze, np.full(kx, n)]
    lo_t, hi_t = zs + 1, ze - 2      # conservative on both sides: a hop needs a wire free one op before and after
    # every inserted op must run OUTSIDE PUSH_CONDITION blocks: a hop happens right after op `until`, so the usable
    # end of an interval is the last position at or before it where the block depth after the op is 0
    dep_after = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
    last0 = np.maximum.accumulate(np.where(dep_after == 0, np.arange(n), -1))
    hi_t = np.where(hi_t >= 0, last0[np.clip(hi_t, 0, n - 1)], hi_t)
    okz = hi_t >= lo_t
    zw, lo_t, hi_t = zw[okz], lo_t[okz], hi_t[okz]
    # merge overlapping intervals of one wire (sets from several sources may repeat or nest)
    o_ = np.lexsort((lo_t, zw))
    zw, lo_t, hi_t = zw[o_], lo_t[o_], hi_t[o_]
    mw, ml, mh = [], [], []
    for w_, l_, h_ in zip(zw.tolist(), lo_t.tolist(), hi_t.tolist()):
        if mw and mw[-1] == w_ and l_ <= mh[-1]:      # share an after-op time; NOT merely adjacent (op between)
            mh[-1] = max(mh[-1], h_)
        else:
            mw.append(w_)
            ml.append(l_)
            mh.append(h_)
    zw, lo_t, hi_t = np.array(mw, np.int64), np.array(ml, np.int64), np.array(mh, np.int64)
    drop = os.environ.get('DQI_PP_DROP')                # npz (wire, time): drop the interval holding that moment
    if drop:
        D = np.load(drop)
        bad = np.zeros(len(zw), bool)
        for w_, t_ in zip(D['wire'].tolist(), D['time'].tolist()):
            bad |= (zw == w_) & (lo_t <= t_) & (hi_t >= t_)
        zw, lo_t, hi_t = zw[~bad], lo_t[~bad], hi_t[~bad]
        print('dropped %d intervals refuted by samples' % int(bad.sum()), flush=True)
    print('zero intervals after merging per wire: %d' % len(zw), flush=True)
    pts = np.unique(np.concatenate([lo_t, hi_t + 1, lk_lo, lk_hi + 1, lk_lo + 2, lk_hi - 1, [-1, n]]))
    seg_of = lambda t: np.searchsorted(pts, t, side='right') - 1  # noqa: E731
    diff = np.zeros(len(pts) + 1, np.int64)
    np.add.at(diff, seg_of(lo_t), 1)
    np.add.at(diff, seg_of(hi_t + 1), -1)
    free = np.cumsum(diff)[:len(pts)]
    size = 1
    while size < len(free):
        size *= 2
    mn, lz = st_build(free.astype(np.int64), size)
    order_mode = os.environ.get('DQI_PP_ORDER', 'gap')      # gap: shortest first; end: earliest end first
    if order_mode == 'end':
        lord = np.lexsort((lk_hi - lk_lo, lk_hi))
    else:
        lord = np.argsort(lk_hi - lk_lo, kind='stable')
    adm = admit(lk_lo, lk_hi, lk_prod, lk_k, lord, prod_off, nocc, pts, mn, lz, size)
    print('links admitted %d of %d (%.1fs)' % (int(adm.sum()), len(adm), time.time() - t0), flush=True)
    # chains: maximal runs of admitted links per product
    chains = []
    for p in range(len(prod_occ)):
        b0 = prod_off[p]
        k = 0
        m_ = nocc[p] - 1
        while k < m_:
            if not adm[b0 + k]:
                k += 1
                continue
            s = k
            while k < m_ and adm[b0 + k]:
                k += 1
            chains.append(prod_occ[p][s:k + 1].astype(np.int64))
    keep_rng = os.environ.get('DQI_PP_CHAINS')            # 'lo:hi' bisection subset of the admitted chains
    if keep_rng:
        lo_, hi_ = (int(x) for x in keep_rng.split(':'))
        chains = chains[lo_:hi_]
    excl = os.environ.get('DQI_PP_EXCLUDE')               # comma list of first-occurrence ops to leave out
    if excl:
        ex = set(int(x) for x in excl.split(','))
        chains = [ch for ch in chains if int(ch[0]) not in ex]
    saved = sum(len(ch) - 1 for ch in chains)
    print('chains %d, Toffolis saved %d (chain lengths: %s)' % (len(chains), saved,
          np.bincount([len(ch) for ch in chains]).tolist()), flush=True)
    # roving wire assignment; a chain that cannot get a wire at some moment is dropped and the sweep redone
    starts = {}
    for w, l, h in zip(zw.tolist(), lo_t.tolist(), hi_t.tolist()):
        starts.setdefault(l, []).append((w, h))

    def sweep(chains):
        cstart, cend, cmid = {}, {}, {}
        for idx, ch in enumerate(chains):
            cstart.setdefault(int(ch[0]) - 1, []).append(idx)
            cend.setdefault(int(ch[-1]), []).append(idx)
            for o in ch[1:-1].tolist():
                cmid.setdefault(int(o), []).append(idx)
        free_heap = []
        holding = {}
        pwire, start_wire, end_wire, mid_wire, hops = {}, {}, {}, {}, {}

        def take_free(t):
            while free_heap:
                nu, w = heapq.heappop(free_heap)
                if -nu >= t and w not in holding:
                    return w, -nu
            return None, None

        queue = sorted(set(starts) | set(cstart) | set(cend) | set(cmid))
        heapq.heapify(queue)
        seen = set()
        while queue:
            t = heapq.heappop(queue)
            if t in seen:
                continue
            seen.add(t)
            for w, h in starts.get(t, []):
                heapq.heappush(free_heap, (-h, w))
            for idx in cmid.get(t, []):
                mid_wire[(idx, t)] = pwire[idx]
            for idx in cend.get(t, []):
                w = pwire.pop(idx)
                end_wire[idx] = w
                until = holding.pop(w)[1]
                if until > t:
                    heapq.heappush(free_heap, (-until, w))
            for idx in cstart.get(t, []):
                w, until = take_free(t + 1)
                if w is None:
                    return idx, None
                holding[w] = (idx, until)
                pwire[idx] = w
                start_wire[idx] = w
                heapq.heappush(queue, until)
            for w in [w for w, (idx, until) in holding.items() if until == t]:
                idx, _ = holding.pop(w)
                w2, until2 = take_free(t + 1)
                if w2 is None:
                    return idx, None
                holding[w2] = (idx, until2)
                pwire[idx] = w2
                hops.setdefault(t, []).append((w, w2))
                heapq.heappush(queue, until2)
        return None, (start_wire, end_wire, mid_wire, hops)

    ndrop = 0
    while True:
        bad, res = sweep(chains)
        if bad is None:
            break
        chains = chains[:bad] + chains[bad + 1:]
        ndrop += 1
    start_wire, end_wire, mid_wire, hops = res
    if ndrop:
        print('chains dropped by the wire assignment: %d; kept %d saving %d' % (
            ndrop, len(chains), sum(len(ch) - 1 for ch in chains)), flush=True)
    print('hops %d (%.1fs)' % (sum(len(v) for v in hops.values()), time.time() - t0), flush=True)
    nb = c.num_bits
    inserts = {}
    delete = set()
    for t, lst in hops.items():
        inserts.setdefault(t, []).extend((0x22, (a, b), 0, 0) for a, b in lst)
    def first_ops(o, w):
        a_, b_, t_ = (int(x) for x in W[o])
        if (ops[o] & 0x7F) == 0x32:   # CSWAP(c=a_, x=b_, y=t_): w = c (x + y), then x ^= w, y ^= w
            return [(0x21, (b_, t_), 0, 0), (0x31, (a_, t_, w), 0, 0), (0x21, (b_, t_), 0, 0),
                    (0x21, (w, b_), 0, 0), (0x21, (w, t_), 0, 0)]
        return [(0x31, (a_, b_, w), 0, 0), (0x21, (w, t_), 0, 0)]

    def use_ops(o, w):
        a_, b_, t_ = (int(x) for x in W[o])
        if (ops[o] & 0x7F) == 0x32:
            return [(0x21, (w, b_), 0, 0), (0x21, (w, t_), 0, 0)]
        return [(0x21, (w, t_), 0, 0)]

    def erase_ops(o, w, m_):
        a_, b_, t_ = (int(x) for x in W[o])
        if (ops[o] & 0x7F) == 0x32:   # phase m c (x + y): two conditioned CZ (x + y is unchanged by the swap)
            return [(0x50, (w,), m_, 0), (0xA0, (a_, b_), 0, m_), (0xA0, (a_, t_), 0, m_)]
        return [(0x50, (w,), m_, 0), (0xA0, (a_, b_), 0, m_)]

    for idx, ch in enumerate(chains):
        o_s, o_e = int(ch[0]), int(ch[-1])
        a, b, t = (int(x) for x in W[o_s])
        w0 = start_wire[idx]
        m = nb + idx
        a2, b2, t2 = (int(x) for x in W[o_e])
        w1 = end_wire[idx]
        if exact_all:
            inserts[o_s] = first_ops(o_s, w0) + inserts.get(o_s, [])
            delete.add(o_s)
            for o in ch[1:-1].tolist():
                inserts[o] = use_ops(o, mid_wire[(idx, o)]) + inserts.get(o, [])
                delete.add(o)
            inserts[o_e] = use_ops(o_e, w1) + erase_ops(o_e, w1, m) + inserts.get(o_e, [])
            delete.add(o_e)
            continue
        if (ops[o_s] & 0x7F) == 0x32:
            # CSWAP(c=a, x=b, y=t): w = c * (x + y), then x ^= w, y ^= w (one Toffoli, same as the CSWAP)
            # CX(x, y): y = x + y; CCX(c, y, w): w = d; CX(x, y): y restored; then the swap as two CX from w
            inserts[o_s] = [(0x21, (b, t), 0, 0), (0x31, (a, t, w0), 0, 0), (0x21, (b, t), 0, 0),
                            (0x21, (w0, b), 0, 0), (0x21, (w0, t), 0, 0)] + inserts.get(o_s, [])
            delete.add(o_s)
            for o in ch[1:-1].tolist():
                wk = mid_wire[(idx, o)]
                inserts[o] = [(0x21, (wk, int(W[o, 1])), 0, 0), (0x21, (wk, int(W[o, 2])), 0, 0)] + inserts.get(o, [])
                delete.add(o)
            inserts[o_e] = [(0x21, (w1, b2), 0, 0), (0x21, (w1, t2), 0, 0), (0x50, (w1,), m, 0),
                            (0xA0, (a2, b2), 0, m), (0xA0, (a2, t2), 0, m)] + inserts.get(o_e, [])
            delete.add(o_e)
            continue
        inserts[o_s] = [(0x31, (a, b, w0), 0, 0), (0x21, (w0, t), 0, 0)] + inserts.get(o_s, [])
        delete.add(o_s)
        for o in ch[1:-1].tolist():
            tk = int(W[o, 2])
            inserts[o] = [(0x21, (mid_wire[(idx, o)], tk), 0, 0)] + inserts.get(o, [])
            delete.add(o)
        inserts[o_e] = [(0x21, (w1, t2), 0, 0), (0x50, (w1,), m, 0), (0xA0, (a2, b2), 0, m)] + inserts.get(o_e, [])
        delete.add(o_e)
    out = kw.edit(c, delete, inserts)
    out = dataclasses.replace(out, num_bits=nb + len(chains))
    nqub = kw.write(dst, out)
    print('qubits %d, bits %d -> %d, static Toffoli-like %d -> %d (%.1fs)' % (
        nqub, nb, out.num_bits, int(((ops & 0x7C) == 0x30).sum()), int(((out.ops & 0x7C) == 0x30).sum()),
        time.time() - t0))
    # scratch checkpoints (after-op time, wire) for an independent zero check: each chain start and each hop target
    chk = [(int(ch[0]) - 1, start_wire[idx]) for idx, ch in enumerate(chains)]
    chk += [(t, b) for t, lst in hops.items() for (a, b) in lst]
    np.savez(dst + '.assign.npz', chk=np.array(chk, dtype=np.int64).reshape(-1, 2),
             chain_first=np.array([int(ch[0]) for ch in chains], np.int64),
             chain_last=np.array([int(ch[-1]) for ch in chains], np.int64),
             chain_len=np.array([len(ch) for ch in chains], np.int64),
             start_wire=np.array([start_wire[i] for i in range(len(chains))], np.int64))


if __name__ == '__main__':
    main()
