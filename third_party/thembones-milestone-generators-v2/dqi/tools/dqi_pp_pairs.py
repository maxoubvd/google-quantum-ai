#!/usr/bin/env python3
"""dqi_pp_pairs.py -- workstream dqi_pp (P1): retained SUMS, measured exactly (no rewrite).

Pairs: consecutive occurrences of the same EXACT product form (tools/dqi_pp_gvn.py hashes) among free CCX.
Segment: a maximal run of events on one wire that are all free-CCX-target events (the wire is neither read nor
otherwise written inside it). Group: the pairs whose first CCX lies in forward segment F and second in segment M.
A group of n pairs needs ONE zero-idle scratch wire s from its first forward CCX to the start of M:
  forward members CCX(a, b, t) -> CCX(a, b, s); CX(s, t) after the last member (t is not read inside F)
  before M: CX(s, t'); HMR(s) -> m;  mirror members CCX(a', b', t') -> CZ(a', b') if m
saving n Toffolis (n = 1 is B5's retained product). F == M: the product cancels in place, both CCX deleted.
Zero-idle wires: structural (after an unconditioned R / HMR outside blocks until the next event; untouched
non-input wires from the start). Scheduling: greedy, largest n first, one group at a time per wire interval.
Usage: python dqi_pp_pairs.py <circuit.kmb> <gvn.npz> <out.npz>
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_prove as P  # noqa: E402


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    G = np.load(sys.argv[2])
    out = sys.argv[3]
    ops, W0, W1, W2, BT, CB, DEP = P.decode(c)
    n = len(ops)
    base = ops & 0x7F
    t_op = G['t_op']
    free = (G['t_lvl'] == 0) & (base[t_op] == 0x31)
    fo = t_op[free]
    p0, p1, pc = G['t_p0'][free], G['t_p1'][free], G['t_pc'][free]
    print('free CCX %d (%.1fs)' % (len(fo), time.time() - t0), flush=True)
    # events (op, wire) for every operand; role 1 = target of a free CCX
    Wm = np.stack([W0, W1, W2], axis=1)
    ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
    ev_w = Wm.reshape(-1).astype(np.int64)
    ev_slot = np.tile(np.arange(3), n)
    k = ev_w >= 0
    ev_op, ev_w, ev_slot = ev_op[k], ev_w[k], ev_slot[k]
    isfree = np.zeros(n, bool)
    isfree[fo] = True
    role = isfree[ev_op] & (ev_slot == 2)
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, role = ev_op[o], ev_w[o], role[o]
    ne = len(ev_op)
    print('events %d (%.1fs)' % (ne, time.time() - t0), flush=True)
    # segments: new segment where wire changes, or role changes, or previous event not a free-CCX target
    newseg = np.ones(ne, bool)
    newseg[1:] = (ev_w[1:] != ev_w[:-1]) | ~role[1:] | ~role[:-1]
    segid = np.cumsum(newseg) - 1
    seg_start = ev_op[newseg]                       # first op of each segment
    # previous event op on the same wire before each segment start (for placement), and next after its end
    # segment of each free CCX (its target event)
    tgt_ev = np.flatnonzero(role)
    ccx_seg = np.full(n, -1, np.int64)
    ccx_seg[ev_op[tgt_ev]] = segid[tgt_ev]
    sg = ccx_seg[fo]
    # exact pairs: sort by product form then op
    order = np.lexsort((fo, pc, p1, p0))
    same = (p0[order][1:] == p0[order][:-1]) & (p1[order][1:] == p1[order][:-1]) & (pc[order][1:] == pc[order][:-1])
    pi = order[:-1][same]
    pj = order[1:][same]
    print('exact repeated pairs (consecutive) %d (%.1fs)' % (len(pi), time.time() - t0), flush=True)
    fi, fj = fo[pi], fo[pj]
    si, sj = sg[pi], sg[pj]
    cancel = si == sj
    print('pairs cancelling inside one segment: %d' % int(cancel.sum()))
    # groups by (forward segment, mirror segment)
    key = si * np.int64(1 << 32) + sj
    keyu, inv, cnt = np.unique(key[~cancel], return_inverse=True, return_counts=True)
    gfirst = np.full(len(keyu), np.iinfo(np.int64).max)
    np.minimum.at(gfirst, inv, fi[~cancel])
    gms = seg_start[(keyu & ((1 << 32) - 1)).astype(np.int64)]
    hist = np.bincount(np.minimum(cnt, 20))
    print('groups %d, pairs in groups %d; size histogram (1..20+): %s' % (len(keyu), int(cnt.sum()),
          hist[1:].tolist()), flush=True)
    gap = gms - gfirst
    print('group span ops: median %d, p10 %d, p90 %d' % (int(np.median(gap)), int(np.percentile(gap, 10)),
          int(np.percentile(gap, 90))))
    # zero-idle intervals per wire (structural)
    reset = np.isin(base[ev_op], (0x50, 0x12)) & (CB[ev_op] < 0) & (DEP[ev_op] == 0)
    zs, ze, zw = [], [], []
    wstart = np.searchsorted(ev_w, np.arange(c.num_qubits), side='left')
    wend = np.searchsorted(ev_w, np.arange(c.num_qubits), side='right')
    for w in range(c.num_qubits):
        e = ev_op[wstart[w]:wend[w]]
        r = reset[wstart[w]:wend[w]]
        if w >= 840:
            zs.append(-1)
            ze.append(e[0] if len(e) else n)
            zw.append(w)
        idx = np.flatnonzero(r)
        nxt = np.append(e[1:], n)
        zs.extend(e[idx].tolist())
        ze.extend(nxt[idx].tolist())
        zw.extend([w] * len(idx))
    zs, ze, zw = np.array(zs), np.array(ze), np.array(zw)
    print('zero-idle intervals %d (%.1fs)' % (len(zs), time.time() - t0), flush=True)
    # greedy scheduling: largest groups first; a group needs an interval with start < first and end > gms
    zo = np.argsort(zs)
    zs, ze, zw = zs[zo], ze[zo], zw[zo]
    used = {}
    chosen = []
    saved = 0
    for g in np.lexsort((gap, -cnt)):
        a_, b_ = int(gfirst[g]), int(gms[g])
        hi = np.searchsorted(zs, a_, side='left')     # intervals starting before the first member
        cand = np.flatnonzero(ze[:hi] > b_)
        for q in cand[::-1]:
            L = used.setdefault(int(q), [])
            if any(not (b_ < x or a_ > y) for x, y in L):
                continue
            L.append((a_, b_))
            chosen.append((g, int(zw[q])))
            saved += int(cnt[g])
            break
    print('scheduled groups %d, Toffolis saved %d (+ %d cancelling pairs x 2)' % (len(chosen), saved,
          int(cancel.sum())))
    np.savez(out, pi=fi, pj=fj, si=si, sj=sj, cancel=cancel, keyu=keyu, cnt=cnt, gfirst=gfirst, gms=gms,
             chosen=np.array(chosen, dtype=np.int64).reshape(-1, 2))
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
