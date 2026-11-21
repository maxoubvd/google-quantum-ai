#!/usr/bin/env python3
"""dqi_int_inv17map.py -- workstream dqi_int (B5, Nov 15 twelfth holder): R46-DQI-EXT's INV17 map
(tools/r46_dqi_ext_inv17.py map, credit R46-DQI-EXT) made to work on files whose CONVERTER differs from the marker
file's (zpark keeps flags allocated, so every later allocation of convert / unconvert gets other wire numbers and the
one-way greedy alignment stops there: 'no resync' at the start of convert).

The same selection (ops of inverse6_joint18_12extra in extract / unextract with a or b = 1 on every sample of the
stat files; class junction = under external_node_chain), aligned to the target file TWICE: forward from the first
Toffoli (reaches extract) and backward from the last (reaches unextract). An op mapped both ways must map to the
same target op (else it is dropped and reported). Writes <prefix>_A.npz (junction) and <prefix>_B.npz (all) in the
format of tools/r44_xfer_dqidrop.py build2 (k, mode = 4).
Usage: python dqi_int_inv17map.py <marker.kmb> <segments_real.json> <target.kmb> <out prefix> <stat.npz> [...]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r46_dqi_ext_inv17 import toffolis  # noqa: E402


def align(ms, ls):
    """R46-DQI-EXT's greedy alignment with resync; returns mp (marker index -> target index or -1) and the stop."""
    mp = np.full(len(ms), -1, np.int64)
    i = j = 0
    W = 6
    res = 0
    stop = None
    while i < len(ms) and j < len(ls):
        if ms[i] == ls[j]:
            m = min(len(ms) - i, len(ls) - j)
            neq = np.flatnonzero(ms[i:i + m] != ls[j:j + m])
            r = int(neq[0]) if len(neq) else m
            mp[i:i + r] = np.arange(j, j + r)
            i += r
            j += r
            continue
        found = None
        for tot in range(1, 6000):
            cands = [(tot, 0), (0, tot)]
            if tot <= 60:
                cands += [(di, tot - di) for di in range(1, tot)]
            for di, dj in cands:
                a, b = i + di, j + dj
                w = min(W, len(ms) - a, len(ls) - b)
                if w > 0 and (ms[a:a + w] == ls[b:b + w]).all():
                    found = (di, dj)
                    break
            if found:
                break
        if not found:
            stop = (i, j)
            break
        res += 1
        i += found[0]
        j += found[1]
    return mp, res, stop


def main():
    mk, segf, line, pre = sys.argv[1:5]
    S = json.load(open(segf))
    marks, st, segs = S['marks'], S['stacks'], S['segments']
    act = None
    n = 0
    for f in sys.argv[5:]:
        d = np.load(f)
        if act is None:
            act, ones, pa, pb = (d[k].astype(np.int64).copy() for k in ('act', 'ones', 'pa', 'pb'))
            seg = d['seg']
        else:
            act += d['act']
            ones += d['ones']
            pa += d['pa']
            pb += d['pb']
        n += int(d['samples'])
    rel = np.flatnonzero((act == n) & (pa + pb - ones == act) & (ones > 0))
    pick = []
    for k in rel.tolist():
        sid, mkid, _n = segs[int(seg[k])]
        names = [x.split('@')[0] for x in st[sid]]
        if not marks[mkid] or marks[mkid][0] not in ('extract', 'unextract'):
            continue
        if 'inverse6_joint18_12extra' not in names:
            continue
        pick.append((k, marks[mkid][0], 'junction' if 'external_node_chain' in names else 'other'))
    print('marker file: %d samples; picked %d (junction %d, other %d)' % (
        n, len(pick), sum(1 for p in pick if p[2] == 'junction'), sum(1 for p in pick if p[2] == 'other')))
    mo, mt, ms = toffolis(mk)
    lo, lt, ls = toffolis(line)
    mpf, rf, sf = align(ms, ls)
    mpr, rr, sr = align(ms[::-1].copy(), ls[::-1].copy())
    mpb = np.full(len(ms), -1, np.int64)
    ok = mpr >= 0
    mpb[len(ms) - 1 - np.flatnonzero(ok)] = len(ls) - 1 - mpr[ok]
    print('forward: %d aligned, %d resyncs, stop %s; backward: %d aligned, %d resyncs, stop %s (marker %d, target %d)'
          % (int((mpf >= 0).sum()), rf, sf, int((mpb >= 0).sum()), rr,
             None if sr is None else (len(ms) - 1 - sr[0], len(ls) - 1 - sr[1]), len(ms), len(ls)))
    out = {'A': [], 'B': []}
    lost = conflict = 0
    for k, stg, cls in pick:
        a, b = int(mpf[k]), int(mpb[k])
        if a >= 0 and b >= 0 and a != b:
            conflict += 1
            continue
        t = a if a >= 0 else b
        if t < 0:
            lost += 1
            continue
        out['B'].append(t)
        if cls == 'junction':
            out['A'].append(t)
    print('picked ops mapped: %d of %d (junction %d); lost %d, conflicts %d' % (
        len(out['B']), len(pick), len(out['A']), lost, conflict))
    for tag in ('A', 'B'):
        ks = np.array(sorted(set(out[tag])), np.int64)
        assert (lo[lt[ks]] & 0x7F == 0x31).all(), 'not all CCX'
        np.savez(pre + '_' + tag + '.npz', k=ks, mode=np.full(len(ks), 4, np.int64))
        print('written %s_%s.npz: %d CCX' % (pre, tag, len(ks)))


if __name__ == '__main__':
    main()
