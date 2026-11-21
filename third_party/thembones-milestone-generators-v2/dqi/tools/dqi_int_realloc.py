#!/usr/bin/env python3
"""dqi_int_realloc.py -- B5 (16th holder, Nov 16): WIRE RE-ASSIGNMENT of a DQI .kmb into Q wires (the 901 file).

Input: a marker-free .kmb and its exact zero intervals (tools/dqi_pp_zero.py on THAT file: a wire is |0> at after-op
times start .. end-1, end = its next touch). Every wire's timeline is cut into BUSY segments [a, b] (closed op
ranges between zero intervals: a = the touch that ends a zero interval, b = the touch after which the next one
starts). Segments are coloured with Q colours by one sweep in order of a (greedy interval colouring, optimal for
intervals): the 840 syndrome wires (registers S0..S69) keep their own colour in every segment (input and output
stay in place); any other segment takes a colour that is free over [a, b], preferring its own wire id, else the free
colour whose next syndrome segment starts last-but-after b (best fit). Ops are rewritten with the colour of the
segment holding each operand; an operand touch that lies in no segment (ops that never execute) gets colour 0.
Exact by construction: two segments share a colour only if they share no op, and a colour is |0> between them.
Usage: python dqi_int_realloc.py <in.kmb> <zero.npz> <out.kmb> <Q>   (B5_REALLOC_PIN=ends: pin only the first and last
segment of each syndrome wire)
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

SH = 27
NIN = 840


def main():
    t0 = time.time()
    src, zp, dst, Q = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    k = kmb.read(src)
    ops = np.asarray(k.ops)
    n = len(ops)
    nq = k.num_qubits
    M = kmb.ID_MASK
    base = ops & 0x7F
    assert not np.isin(base[ops != 0x81], (0x07, 0x17, 0x47)).any(), 'strip the markers first'
    kk = ((ops & 0x30) >> 4).astype(np.int8)
    i1 = np.cumsum(kk == 1, dtype=np.int64) - 1
    i2 = np.cumsum(kk == 2, dtype=np.int64) - 1
    i3 = np.cumsum(kk == 3, dtype=np.int64) - 1
    # touches: (wire, op) for every qubit operand
    tw, to = [], []
    for cls, arrs, ix in ((1, (k.q1,), i1), (2, k.q2, i2), (3, k.q3, i3)):
        idx = np.flatnonzero(kk == cls)
        for a in arrs:
            tw.append(np.asarray(a).astype(np.int64)[ix[idx]] & M)
            to.append(idx)
    tw = np.concatenate(tw)
    to = np.concatenate(to)
    # zero intervals per wire, merged
    Z = np.load(zp)
    zw, zs, ze = Z['wire'].astype(np.int64), Z['start'].astype(np.int64), Z['end'].astype(np.int64)
    o = np.lexsort((zs, zw))
    zw, zs, ze = zw[o], zs[o], ze[o]
    mw, ms, me = [], [], []
    for w_, s_, e_ in zip(zw.tolist(), zs.tolist(), ze.tolist()):
        if mw and mw[-1] == w_ and s_ < me[-1]:       # strict overlap only: s == previous end is a touch
            me[-1] = max(me[-1], e_)
        else:
            mw.append(w_)
            ms.append(s_)
            me.append(e_)
    zw, zs, ze = np.array(mw, np.int64), np.array(ms, np.int64), np.array(me, np.int64)
    print('%s: Q %d -> %d, ops %d, touches %d, zero intervals %d (%.0fs)' % (src, nq, Q, n, len(tw), len(zw),
                                                                            time.time() - t0), flush=True)
    # busy segments per wire: complement of the zero intervals on [-1, n]
    sw, sa, sb = [], [], []
    for w in range(nq):
        lo, hi = np.searchsorted(zw, w, 'left'), np.searchsorted(zw, w, 'right')
        s_, e_ = zs[lo:hi], ze[lo:hi]
        cur = -1                                 # the wire is busy from cur (inclusive) until the next zero start
        if len(s_) and s_[0] <= -1:
            cur = int(e_[0])
            s_, e_ = s_[1:], e_[1:]
        elif w >= NIN:
            # never zero from the start (should not happen for ancillas): busy from -1
            pass
        for a_, b_ in zip(s_.tolist(), e_.tolist()):
            if a_ >= cur:
                sw.append(w)
                sa.append(cur)
                sb.append(a_)                    # busy through the touch a_ (after which it is zero)
            cur = b_
        if cur < n:
            sw.append(w)
            sa.append(cur)
            sb.append(n)                         # busy to the end
    sw, sa, sb = np.array(sw, np.int64), np.array(sa, np.int64), np.array(sb, np.int64)
    # drop empty segments of wires never touched after a zero start at n
    keep = sb >= sa
    sw, sa, sb = sw[keep], sa[keep], sb[keep]
    # B5_REALLOC_PIN=ends (17th holder): only a syndrome wire's FIRST segment (the input, from -1) and LAST segment (the
    # output, to the end) keep the wire's own colour; its middle segments are coloured like any other (default 'all':
    # every syndrome segment pinned, the 16th holder's rule, kept for reproducibility)
    if os.environ.get('B5_REALLOC_PIN', 'all') == 'ends':
        pinned = (sw < NIN) & ((sa == -1) | (sb == n))
    else:
        pinned = sw < NIN
    print('busy segments %d (syndrome wires %d) (%.0fs)' % (len(sw), int(pinned.sum()), time.time() - t0), flush=True)
    # sweep
    order = np.lexsort((~pinned, sa))            # by start; pinned first at equal start
    busy_until = np.full(Q, -2, np.int64)
    pin_starts = [[] for _ in range(Q)]
    for j in np.flatnonzero(pinned).tolist():
        pin_starts[int(sw[j])].append((int(sa[j]), int(sb[j])))
    for c in range(Q):
        pin_starts[c].sort()
    pin_ptr = np.zeros(Q, np.int64)
    INF = np.int64(1) << 62
    nps = np.array([pin_starts[c][0][0] if pin_starts[c] else INF for c in range(Q)], np.int64)
    color = np.full(len(sw), -1, np.int64)
    for j in order.tolist():
        a, b, w = int(sa[j]), int(sb[j]), int(sw[j])
        if pinned[j]:
            c = w
            assert busy_until[c] < a, ('pinned colour busy', w, a, int(busy_until[c]))
            busy_until[c] = b
            pin_ptr[c] += 1
            nps[c] = pin_starts[c][pin_ptr[c]][0] if pin_ptr[c] < len(pin_starts[c]) else INF
            color[j] = c
            continue
        ok = (busy_until < a) & (nps > b)
        if w < Q and ok[w]:
            c = w
        else:
            cand = np.flatnonzero(ok)
            if len(cand) == 0:
                live = int((busy_until >= a).sum())
                print('FAILED at op %d (segment of wire %d, %d..%d): %d colours busy, %d blocked by syndrome segments'
                      % (a, w, a, b, live, int(((busy_until < a) & (nps <= b)).sum())), flush=True)
                sys.exit(2)
            c = int(cand[np.argmin(nps[cand])])
        busy_until[c] = b
        color[j] = c
    print('coloured (%.0fs)' % (time.time() - t0), flush=True)
    # remap every touch: the segment of its wire with a <= op <= b
    key = (sw << SH) | (sa + 1)                  # +1: segments that start at -1 (inputs)
    o2 = np.argsort(key)
    key_s, b_s, c_s, w_s = key[o2], sb[o2], color[o2], sw[o2]
    pos = np.searchsorted(key_s, (tw << SH) | (to + 1), 'right') - 1
    okp = (pos >= 0)
    posc = np.maximum(pos, 0)
    inseg = okp & (w_s[posc] == tw) & (b_s[posc] >= to)
    newc = np.where(inseg, c_s[posc], 0)
    print('touches outside any segment (never-executed ops): %d' % int((~inseg).sum()), flush=True)
    # write back per operand array, in the same order the touches were built
    pos_ = 0
    out_q1 = np.asarray(k.q1).copy()
    out_q2 = [np.asarray(x).copy() for x in k.q2]
    out_q3 = [np.asarray(x).copy() for x in k.q3]
    for cls, arrs, ix in ((1, [out_q1], i1), (2, out_q2, i2), (3, out_q3, i3)):
        idx = np.flatnonzero(kk == cls)
        for a in arrs:
            m = len(idx)
            seg = newc[pos_:pos_ + m]
            hi_bits = a[ix[idx]] & ~np.uint32(M)
            a[ix[idx]] = (seg.astype(np.uint32) | hi_bits).astype(a.dtype)
            pos_ += m
    out = kmb.Kmb(Q, k.num_bits, k.registers, k.ops, tuple(out_q3), tuple(out_q2), out_q1, k.bit_t, k.cond,
                  k.angles, k.tags)
    nqw = kw.write(dst, out)
    print('written %s: Q %d (%.0fs)' % (dst, nqw, time.time() - t0))


if __name__ == '__main__':
    main()
