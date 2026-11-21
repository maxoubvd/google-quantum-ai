#!/usr/bin/env python3
"""dqi_pp_windows.py -- workstream dqi_pp (P1): which wire, between which gates, would have to be free for the
exact retained-product pass to admit more links (the "wire-time" list for B5 and D1).

Rebuilds exactly what tools/dqi_pp_ret.py does in DQI_PP_EXACTALL=1 mode (links of CCX / CSWAP with equal exact
forms, exact zero intervals, conservative hop margins, shortest-first admission) and then:
 1. residual capacity after the baseline admission; windows = maximal runs of segments with residual 0;
 2. for each rejected link, the windows inside its span; windows ranked by links blocked by that window alone;
 3. for the top windows, the TRUE gain: admission rerun with one more free wire over the window;
 4. per top window, the wires untouched over it that hold a value (not zero there), with the gates that bracket
    their idle stretch: freeing such a wire there (parking or recomputing its value) gives the gain.
Usage: python dqi_pp_windows.py <raw.kmb> <gvn.npz> <zero.npz> [top=10] [cand=40]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_ret import admit, st_build, st_query  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402
from numba import njit  # noqa: E402


def build(c, G, Z, ops, n):
    bt_ = ops[G['t_op']] & 0x7F
    msw = ((bt_ == 0x32) | (bt_ == 0x31)) & (G['t_lvl'] == 0)
    so = G['t_op'][msw]
    k0, k1, kc = G['t_p0'][msw], G['t_p1'][msw], G['t_pc'][msw]
    srt = np.lexsort((so, kc, k1, k0))
    ks = np.stack([k0[srt], k1[srt], kc[srt].astype(np.uint64)], 1)
    brk = np.flatnonzero(np.r_[True, np.any(ks[1:] != ks[:-1], axis=1), True])
    prod_occ, lk_lo, lk_hi, lk_prod, lk_k = [], [], [], [], []
    for g_ in range(len(brk) - 1):
        occ = so[srt[brk[g_]:brk[g_ + 1]]]
        if len(occ) < 2:
            continue
        p = len(prod_occ)
        prod_occ.append(occ.astype(np.int64))
        for kk in range(len(occ) - 1):
            lk_lo.append(int(occ[kk]) - 1)
            lk_hi.append(int(occ[kk + 1]))
            lk_prod.append(p)
            lk_k.append(kk)
    lk_lo, lk_hi = np.array(lk_lo, np.int64), np.array(lk_hi, np.int64)
    lk_prod, lk_k = np.array(lk_prod, np.int64), np.array(lk_k, np.int64)
    nocc = np.array([len(x) for x in prod_occ], np.int64)
    prod_off = np.r_[0, np.cumsum(nocc - 1)[:-1]].astype(np.int64)
    zw, zs, ze = Z['wire'].astype(np.int64), Z['start'].astype(np.int64), Z['end'].astype(np.int64)
    mlo, mhi = (int(x) for x in os.environ.get('DQI_PP_MARGINS', '1,2').split(','))   # measurement knob
    lo_t, hi_t = zs + mlo, ze - mhi
    dep_after = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
    last0 = np.maximum.accumulate(np.where(dep_after == 0, np.arange(n), -1))
    hi_t = np.where(hi_t >= 0, last0[np.clip(hi_t, 0, n - 1)], hi_t)
    okz = hi_t >= lo_t
    zw, lo_t, hi_t = zw[okz], lo_t[okz], hi_t[okz]
    o_ = np.lexsort((lo_t, zw))
    zw, lo_t, hi_t = zw[o_], lo_t[o_], hi_t[o_]
    mw, ml, mh = [], [], []
    for w_, l_, h_ in zip(zw.tolist(), lo_t.tolist(), hi_t.tolist()):
        if mw and mw[-1] == w_ and l_ <= mh[-1]:
            mh[-1] = max(mh[-1], h_)
        else:
            mw.append(w_)
            ml.append(l_)
            mh.append(h_)
    zw, lo_t, hi_t = np.array(mw, np.int64), np.array(ml, np.int64), np.array(mh, np.int64)
    pts = np.unique(np.concatenate([lo_t, hi_t + 1, lk_lo, lk_hi + 1, lk_lo + 2, lk_hi - 1, [-1, n]]))
    seg_of = lambda t: np.searchsorted(pts, t, side='right') - 1  # noqa: E731
    diff = np.zeros(len(pts) + 1, np.int64)
    np.add.at(diff, seg_of(lo_t), 1)
    np.add.at(diff, seg_of(hi_t + 1), -1)
    free = np.cumsum(diff)[:len(pts)]
    lord = np.argsort(lk_hi - lk_lo, kind='stable')
    return dict(lk_lo=lk_lo, lk_hi=lk_hi, lk_prod=lk_prod, lk_k=lk_k, nocc=nocc, prod_off=prod_off, pts=pts,
                free=free, lord=lord, zw=zw, zlo=lo_t, zhi=hi_t, seg_of=seg_of)


def touch_scan(Wq, at, nq, back, chunk=1 << 20):
    """Per wire: last op < at touching it (back) or first op >= at touching it (forward); -1 if none."""
    n = Wq.shape[0]
    res = np.full(nq, -1, np.int64)
    found = np.zeros(nq, bool)
    pos = at
    while (pos > 0 if back else pos < n) and not found.all():
        lo, hi = (max(0, pos - chunk), pos) if back else (pos, min(n, pos + chunk))
        blk = Wq[lo:hi]
        idx = np.broadcast_to(np.arange(lo, hi)[:, None], blk.shape)
        v = blk >= 0
        wq, ix = blk[v], idx[v]
        tmp = np.full(nq, -1 if back else n, np.int64)
        if back:
            np.maximum.at(tmp, wq, ix)
            new = (~found) & (tmp >= 0)
        else:
            np.minimum.at(tmp, wq, ix)
            new = (~found) & (tmp < n)
        res[new] = tmp[new]
        found |= new
        pos = lo if back else hi
    return res


@njit(cache=False)
def proxy_nb(sa, sb, bmin, bmax):
    """Per stretch: rejected links (sorted by bmax) whose zero-residual span lies inside, packed disjointly."""
    out = np.zeros(sa.shape[0], np.int64)
    for k in range(sa.shape[0]):
        last = -1
        cnt = 0
        for t in range(bmax.shape[0]):
            if bmax[t] > sb[k]:
                break
            if bmin[t] >= sa[k] and bmin[t] > last:
                cnt += 1
                last = bmax[t]
        out[k] = cnt
    return out


def wire_mode(c, ops, n, B, adm, resid, zero_seg, top, ncand, t0, lmin=500):
    """Rank (wire, idle stretch) pairs: the wire is untouched between two consecutive touches lo < hi and holds a
    value there (no exact zero interval covers the stretch). Proxy: rejected links whose zero-residual segments all
    lie inside the stretch, packed with disjoint zero-residual spans; true gain: admission rerun with +1 over it."""
    from dqi_pp_gvn import stages
    seg_of = B['seg_of']
    nseg = len(B['free'])
    nqo = (ops & 0x30) >> 4
    Wq = np.full((n, 3), -1, dtype=np.int64)
    for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
        sel = np.flatnonzero(nqo == k)
        for j in range(k):
            Wq[sel, j] = np.asarray(arrs[j], dtype=np.int64)
    flat = Wq.ravel()
    v = flat >= 0
    wires = flat[v]
    opi = np.repeat(np.arange(n, dtype=np.int64), 3)[v]
    del flat, v
    o_ = np.lexsort((opi, wires))
    wires, opi = wires[o_], opi[o_]
    del o_
    same = wires[1:] == wires[:-1]
    g_w, g_lo, g_hi = wires[:-1][same], opi[:-1][same], opi[1:][same]
    del wires, opi, same
    lmax = int(os.environ.get('DQI_PP_WIN_LMAX', '0')) or (1 << 62)
    keep = ((g_hi - g_lo) >= lmin) & ((g_hi - g_lo) <= lmax)
    g_w, g_lo, g_hi = g_w[keep], g_lo[keep], g_hi[keep]
    # keep only stretches that hold at least one zero-residual segment (else +1 there admits nothing)
    cz0 = np.r_[0, np.cumsum(zero_seg.astype(np.int64))]
    s0a, s0b = seg_of(g_lo + 1), seg_of(g_hi - 1)
    h0 = (cz0[s0b + 1] - cz0[s0a]) > 0
    g_w, g_lo, g_hi = g_w[h0], g_lo[h0], g_hi[h0]
    print('stretches >= %d ops over a zero-residual segment: %d (%.1fs)' % (lmin, len(g_w), time.time() - t0),
          flush=True)
    o2 = np.argsort(g_w, kind='stable')
    g_w, g_lo, g_hi = g_w[o2], g_lo[o2], g_hi[o2]
    # drop stretches where the wire is zero (exact zero interval of that wire covering the stretch)
    zw, zlo, zhi = B['zw'], B['zlo'], B['zhi']
    isz = np.zeros(len(g_w), bool)
    bnd = np.flatnonzero(np.r_[True, g_w[1:] != g_w[:-1], True])
    for q_ in range(len(bnd) - 1):
        m = np.arange(bnd[q_], bnd[q_ + 1])
        w = int(g_w[m[0]])
        zm = zw == w
        if not zm.any():
            continue
        zl, zh = zlo[zm], zhi[zm]
        # the wire is untouched over the stretch, so its value is constant there: zero iff some exact zero
        # interval of the wire meets the open stretch (lo, hi)
        k = np.searchsorted(zl, g_hi[m] - 1, side='right') - 1
        okk = k >= 0
        cov = np.zeros(len(m), bool)
        cov[okk] = zh[k[okk]] >= g_lo[m][okk] + 1
        isz[m] = cov
    g_w, g_lo, g_hi = g_w[~isz], g_lo[~isz], g_hi[~isz]
    print('idle stretches with a value, >= %d ops: %d (%.1fs)' % (lmin, len(g_w), time.time() - t0), flush=True)
    # rejected links: zero-residual span inside each
    rej = np.flatnonzero(~adm)
    bmin = np.full(len(rej), -1, np.int64)
    bmax = np.full(len(rej), -1, np.int64)
    zidx = np.flatnonzero(zero_seg)
    for t, L in enumerate(rej.tolist()):
        a, b = seg_of(B['lk_lo'][L]), seg_of(B['lk_hi'][L])
        lo_i, hi_i = np.searchsorted(zidx, a), np.searchsorted(zidx, b, side='right')
        if hi_i > lo_i:
            bmin[t], bmax[t] = zidx[lo_i], zidx[hi_i - 1]
    have = bmin >= 0
    bmin, bmax = bmin[have], bmax[have]
    srt = np.argsort(bmax, kind='stable')
    bmin, bmax = bmin[srt], bmax[srt]
    sa = seg_of(g_lo + 1)
    sb = seg_of(g_hi - 1)
    cz = np.r_[0, np.cumsum(zero_seg.astype(np.int64))]
    hasz = (cz[sb + 1] - cz[sa]) > 0
    g_w, g_lo, g_hi, sa, sb = g_w[hasz], g_lo[hasz], g_hi[hasz], sa[hasz], sb[hasz]
    print('stretches holding a zero-residual segment: %d' % len(g_w), flush=True)
    proxy = proxy_nb(sa.astype(np.int64), sb.astype(np.int64), bmin, bmax)
    order = np.argsort(-proxy, kind='stable')[:ncand]
    print('proxy computed (%.1fs); best proxies %s' % (time.time() - t0, proxy[order[:10]].tolist()), flush=True)
    base = int(adm.sum())
    res = []
    for k in order.tolist():
        if proxy[k] == 0:
            break
        f2 = B['free'].copy()
        f2[sa[k]:sb[k] + 1] += 1
        a2, _, _, _ = run_admit(B, f2)
        res.append((int(a2.sum()) - base, int(proxy[k]), k))
    res.sort(reverse=True)
    st = stages(c)

    def stage_at(op):
        for name, a, b in st:
            if a <= op <= b:
                return name
        return '-'
    regname = {}
    for name, ids, isq in c.registers:
        for pos, (i_, q_) in enumerate(zip(np.asarray(ids).tolist(), np.asarray(isq).tolist())):
            if q_:
                regname[i_] = '%s[%d]' % (name, pos)
    print('\nTOP (wire, idle stretch) by TRUE gain (wire emptied between its two touches; admission rerun):')
    for rank, (g, p, k) in enumerate(res[:top]):
        w, lo, hi = int(g_w[k]), int(g_lo[k]), int(g_hi[k])
        print('%2d. +%d links (proxy %d): wire %d %s, idle between op %d (%s) and op %d (%s), %d ops, stage %s -> %s' %
              (rank + 1, g, p, w, regname.get(w, ''), lo, kmb.NAMES.get(int(ops[lo]) & 0x7F, '?'), hi,
               kmb.NAMES.get(int(ops[hi]) & 0x7F, '?'), hi - lo, stage_at(lo), stage_at(hi)))
    print('done %.1fs' % (time.time() - t0))


def run_admit(B, free):
    size = 1
    while size < len(free):
        size *= 2
    mn, lz = st_build(free.astype(np.int64), size)
    adm = admit(B['lk_lo'], B['lk_hi'], B['lk_prod'], B['lk_k'], B['lord'], B['prod_off'], B['nocc'], B['pts'],
                mn, lz, size)
    return adm, mn, lz, size


def main():
    t0 = time.time()
    raw, gvn, zero = sys.argv[1:4]
    top = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    ncand = int(sys.argv[5]) if len(sys.argv) > 5 else 40
    c = kmb.read(raw)
    ops = np.asarray(c.ops)
    n = len(ops)
    G = np.load(gvn)
    Z = np.load(zero)
    B = build(c, G, Z, ops, n)
    adm, mn, lz, size = run_admit(B, B['free'])
    base = int(adm.sum())
    print('links %d, admitted %d (baseline, %.1fs)' % (len(adm), base, time.time() - t0), flush=True)
    a1, _, _, _ = run_admit(B, B['free'] + 1)
    print('reference: one always-free extra wire: +%d' % (int(a1.sum()) - base), flush=True)
    nseg = len(B['free'])
    resid = np.array([st_query(mn, lz, size, s, s) for s in range(nseg)], np.int64)
    zero_seg = resid <= 0
    # windows: maximal runs of residual-0 segments
    idx = np.flatnonzero(zero_seg)
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1) if len(idx) else []
    win_of_seg = np.full(nseg, -1, np.int64)
    for w, r in enumerate(runs):
        win_of_seg[r] = w
    seg_of = B['seg_of']
    single = np.zeros(len(runs), np.int64)
    anyc = np.zeros(len(runs), np.int64)
    for L in np.flatnonzero(~adm).tolist():
        a = seg_of(B['lk_lo'][L])
        b = seg_of(B['lk_hi'][L])
        ws = np.unique(win_of_seg[a:b + 1])
        ws = ws[ws >= 0]
        anyc[ws] += 1
        if len(ws) == 1:
            single[ws[0]] += 1
    order = np.argsort(-single, kind='stable')[:ncand]
    print('windows %d; rejected links %d; ranked by links blocked by that window alone (%.1fs)' %
          (len(runs), int((~adm).sum()), time.time() - t0), flush=True)
    occ = os.environ.get('DQI_PP_WIN_OCC')       # cost mode: "lo:hi,lo:hi" or a file of "lo hi" lines (op indices)
    if occ:
        if os.path.exists(occ):
            iv = [tuple(int(x) for x in ln.split()[:2]) for ln in open(occ) if ln.strip()]
        else:
            iv = [tuple(int(x) for x in t.split(':')) for t in occ.split(',')]
        tot = B['free'].copy()
        print('\nCOST of held intervals (one wire occupied over each; links lost from the baseline %d):' % base)
        for lo, hi in iv:
            f2 = B['free'].copy()
            a_, b_ = seg_of(lo), seg_of(hi)
            f2[a_:b_ + 1] -= 1
            tot[a_:b_ + 1] -= 1
            a2, _, _, _ = run_admit(B, np.maximum(f2, 0))
            print('  ops %d-%d: -%d links' % (lo, hi, base - int(a2.sum())))
        a2, _, _, _ = run_admit(B, np.maximum(tot, 0))
        print('  all %d together: -%d links' % (len(iv), base - int(a2.sum())))
        return
    if os.environ.get('DQI_PP_WIN_MODE', 'wires') == 'wires':
        wire_mode(c, ops, n, B, adm, resid, zero_seg, top, ncand, t0)
        return
    st = stages(c)

    def stage_at(op):
        for name, a, b in st:
            if a <= op <= b:
                return name
        return '-'
    gains = []
    for w in order.tolist():
        if single[w] == 0:
            break
        f2 = B['free'].copy()
        f2[runs[w]] += 1
        a2, _, _, _ = run_admit(B, f2)
        gains.append((int(a2.sum()) - base, w))
    gains.sort(reverse=True)
    # operands per op for idle-wire search
    nqo = (ops & 0x30) >> 4
    Wq = np.full((n, 3), -1, dtype=np.int64)
    for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
        sel = np.flatnonzero(nqo == k)
        for j in range(k):
            Wq[sel, j] = np.asarray(arrs[j], dtype=np.int64)
    regname = {}
    for name, ids, isq in c.registers:
        for pos, (i_, q_) in enumerate(zip(np.asarray(ids).tolist(), np.asarray(isq).tolist())):
            if q_:
                regname[i_] = '%s[%d]' % (name, pos)
    pts = B['pts']
    print('\nTOP windows by TRUE gain (one more free wire over the window; admission rerun):')
    for rank, (g, w) in enumerate(gains[:top]):
        r = runs[w]
        lo_op = int(pts[r[0]]) + 1
        hi_op = int(pts[r[-1] + 1]) if r[-1] + 1 < len(pts) else n - 1
        seg = Wq[lo_op:hi_op + 1]
        touched = np.unique(seg[seg >= 0])
        busy_zero = set()
        zsel = (B['zlo'] <= lo_op) & (B['zhi'] >= hi_op)
        busy_zero.update(B['zw'][zsel].tolist())
        idle = [q for q in range(c.num_qubits) if q not in set(touched.tolist()) and q not in busy_zero]
        lt_all = touch_scan(Wq, lo_op, c.num_qubits, back=True)
        nt_all = touch_scan(Wq, hi_op + 1, c.num_qubits, back=False)
        cands = []
        for q in idle:
            lt, nt = int(lt_all[q]), int(nt_all[q])
            cands.append((nt - lt if lt >= 0 and nt >= 0 else 10 ** 9, q, lt, nt))
        cands.sort()
        print('%2d. +%d links (blocked by it alone %d, crossing it %d): ops %d-%d (%d ops), stage %s; idle wires with a '
              'value there %d' % (rank + 1, g, single[w], anyc[w], lo_op, hi_op, hi_op - lo_op + 1, stage_at(lo_op),
                                  len(idle)))
        for span, q, lt, nt in cands[:3]:
            gl = kmb.NAMES.get(int(ops[lt]) & 0x7F, '?') if lt >= 0 else '-'
            gn = kmb.NAMES.get(int(ops[nt]) & 0x7F, '?') if nt >= 0 else '-'
            print('      wire %d %s: idle from op %d (%s, %s) to op %d (%s, %s), %d ops' %
                  (q, regname.get(q, ''), lt, gl, stage_at(lt), nt, gn, stage_at(nt), span))
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
