#!/usr/bin/env python3
"""dqi_int_zdiff.py -- B5 (24th holder, Nov 18): WHERE two builds of the same recipe differ in P1's exact free
wire-time (the zero file P1's links pass reads: zero_all.npz = P1's zeros + R55's). The two circuits are aligned by
their INFO_TAG marks (the same tag sequence is required; checked); the free wire-ops of every inter-tag segment are
summed in both, and the segments where the second file has the most LESS free time are listed with their tags.
Why: the EARLYF hold (zd5 / zd6 / zd7) costs P1 exactly 37 links at either level and on two bases (slope 20:34),
with R55's free wire-ops 1.95M lower; the cause was never found (Nov 17).
Usage: python dqi_int_zdiff.py <kmb A> <zero A> <kmb B> <zero B> [top = 25]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402


def profile(path, zpath):
    c = kmb.read(path)
    n = len(c.ops)
    Z = np.load(zpath)
    w, s, e = Z['wire'].astype(np.int64), Z['start'].astype(np.int64), Z['end'].astype(np.int64)
    ok = e >= s
    w, s, e = w[ok], s[ok], e[ok]
    # merge overlapping intervals per wire (P1 merges them too)
    o = np.lexsort((s, w))
    w, s, e = w[o], s[o], e[o]
    mw, ms, me = [], [], []
    for a, b_, c_ in zip(w.tolist(), s.tolist(), e.tolist()):
        if mw and mw[-1] == a and b_ <= me[-1] + 1:
            me[-1] = max(me[-1], c_)
        else:
            mw.append(a)
            ms.append(b_)
            me.append(c_)
    s, e = np.array(ms, np.int64), np.array(me, np.int64)
    d = np.zeros(n + 2, np.int64)
    np.add.at(d, np.clip(s, 0, n), 1)
    np.add.at(d, np.clip(e + 1, 0, n + 1), -1)
    free = np.cumsum(d)[:n]
    cum = np.concatenate([[0], np.cumsum(free)])          # cum[t] = free wire-ops before op t
    pos = np.asarray(c.tag_positions()).astype(np.int64)
    tags = list(c.tags)
    return n, cum, pos, tags, int(free.sum())


def main():
    ka, za, kb, zb = sys.argv[1:5]
    top = int(sys.argv[5]) if len(sys.argv) > 5 else 25
    na, ca, pa, ta, fa = profile(ka, za)
    nb, cb, pb, tb, fb = profile(kb, zb)
    print('A: ops %d, tags %d, free wire-ops %d; B: ops %d, tags %d, free wire-ops %d; B - A = %d' % (
        na, len(ta), fa, nb, len(tb), fb, fb - fa), flush=True)
    assert ta == tb, 'tag sequences differ (%d / %d)' % (len(ta), len(tb))
    edges_a = np.concatenate([[0], pa, [na]])
    edges_b = np.concatenate([[0], pb, [nb]])
    fa_seg = ca[edges_a[1:]] - ca[edges_a[:-1]]
    fb_seg = cb[edges_b[1:]] - cb[edges_b[:-1]]
    diff = fb_seg - fa_seg
    print('segments %d; sum of (B - A) %d; segments with B less: %d, more: %d' % (
        len(diff), int(diff.sum()), int((diff < 0).sum()), int((diff > 0).sum())))
    names = ['<start>'] + ta
    order = np.argsort(diff)[:top]
    for k in order.tolist():
        print('  seg %5d after %-42s A ops %10d-%10d  B - A free wire-ops %10d  (len A %d, B %d)' % (
            k, names[k][:42], edges_a[k], edges_a[k + 1], int(diff[k]), int(edges_a[k + 1] - edges_a[k]),
            int(edges_b[k + 1] - edges_b[k])))
    # running sum: where the deficit accumulates
    run = np.cumsum(diff)
    for q in (0.1, 0.25, 0.5, 0.75, 0.9, 1.0):
        k = int(np.argmax(run <= run[-1] * q)) if run[-1] < 0 else 0
        print('  %3d%% of the deficit reached at segment %d (%s, A op %d)' % (q * 100, k, names[min(k, len(names) - 1)][:40],
                                                                     edges_a[min(k, len(edges_a) - 1)]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
