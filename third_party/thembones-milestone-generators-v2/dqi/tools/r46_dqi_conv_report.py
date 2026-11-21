#!/usr/bin/env python3
"""r46_dqi_conv_report.py -- R46-DQI-CONV: the profile table of convert / unconvert and the flag census.

Inputs: an r7prof profile folder (segments_real.json + executed.json: executed Toffolis per call-stack segment) and,
optionally, census files of tools/r46_dqi_conv_census.py made on the SAME static Toffoli sequence (the raw build).
  prof  <profdir> [stage,...]            kernel table: executed per decode by innermost named frames, share of block
  flags <profdir> <npz,...> [stage,...]  per kernel: executed per decode split by the rate of the op's RARER control
                                         (never, <1e-6, <1e-5, <1e-4, <1e-3, <1e-2, rest)
  wires <profdir> <npz,...> <stage> [maxrate=1e-3] [top=60]  rare control wires: ops, executed, rate, stacks
  stack <profdir> <stage> <substring> [depth=6]   call stacks containing a frame name, executed per decode
"""
import json
import sys
from collections import defaultdict

import numpy as np

REC = ('look_solve', 'handoff_direct_solve', 'handoff_generic_solve', 'look_root_cached', 'look_pd')
BINS = [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
BN = ['never', '<1e-6', '<1e-5', '<1e-4', '<1e-3', '<1e-2', 'rest']


def load(pdir):
    d = json.load(open(pdir + '/segments_real.json'))
    e = json.load(open(pdir + '/executed.json'))
    return d['marks'], d['stacks'], d['segments'], np.array(e['executed'], float) / e['shots'], e['shots']


def names(st):
    return [x.split('@')[0] for x in st]


def key(st, n=2):
    nm = names(st)
    depth = sum(1 for x in nm if x in REC)
    return '%s' % '>'.join(nm[-n:]), depth


def census(files):
    act = ones = ca = cb = None
    n = 0
    for f in files:
        z = np.load(f)
        if act is None:
            act, ones, ca, cb = (z[k].astype(np.int64) for k in ('act', 'ones', 'ca', 'cb'))
            wa, wb, wt, base = z['wa'], z['wb'], z['wt'], z['base']
        else:
            act += z['act']; ones += z['ones']; ca += z['ca']; cb += z['cb']
        n += int(z['samples'])
    return act, ones, ca, cb, wa, wb, wt, base, n


def seg_of_op(seg):
    c = np.array([s[2] for s in seg], np.int64)
    return np.repeat(np.arange(len(seg)), c), int(c.sum())


def main():
    mode, pdir = sys.argv[1], sys.argv[2]
    marks, stacks, seg, ex, shots = load(pdir)
    if mode == 'prof':
        stages = sys.argv[3].split(',') if len(sys.argv) > 3 else ['convert', 'unconvert']
        n = int(sys.argv[4]) if len(sys.argv) > 4 else 2
        thr = float(sys.argv[5]) if len(sys.argv) > 5 else 20
        tot = defaultdict(float)
        for i, (s, mk, c) in enumerate(seg):
            tot[marks[mk][0] if marks[mk] else '-'] += ex[i]
        print('stages: ' + ', '.join('%s %.1f' % (k, v) for k, v in tot.items()) + ' | all %.1f (%d shots)'
              % (sum(tot.values()), shots))
        block = sum(tot[s] for s in stages)
        for stg in stages:
            by = defaultdict(float)
            calls = defaultdict(int)
            for i, (s, mk, c) in enumerate(seg):
                if marks[mk] and marks[mk][0] == stg:
                    k, _ = key(stacks[s], n)
                    tag = '/'.join(marks[mk][1:])
                    by[k] += ex[i]
                    calls[k] += 1
            print('== %s %.1f' % (stg, tot[stg]))
            for k in sorted(by, key=lambda k: -by[k]):
                if by[k] >= thr:
                    print('  %-88s %9.1f %5.1f%% segs %d' % (k[-88:], by[k], 100 * by[k] / block, calls[k]))
        return
    if mode == 'stack':
        stg, sub = sys.argv[3], sys.argv[4]
        dep = int(sys.argv[5]) if len(sys.argv) > 5 else 6
        by = defaultdict(float)
        for i, (s, mk, c) in enumerate(seg):
            if marks[mk] and marks[mk][0] == stg and any(sub in x for x in stacks[s]):
                by['>'.join(names(stacks[s])[-dep:]) + ' [' + '/'.join(marks[mk][1:]) + ']'] += ex[i]
        for k in sorted(by, key=lambda k: -by[k])[:60]:
            print('%9.1f  %s' % (by[k], k))
        print('total %.1f' % sum(by.values()))
        return
    files = sys.argv[3].split(',')
    act, ones, ca, cb, wa, wb, wt, base, n = census(files)
    so, tot = seg_of_op(seg)
    assert tot == len(act), (tot, len(act))
    rate = np.minimum(ca, cb) / np.maximum(act, 1)
    binid = np.searchsorted(np.array(BINS), rate, side='left')      # 0: never ... 6: rest
    binid[(rate > 0) & (binid == 0)] = 1
    exo = act / n
    if mode == 'flags':
        stages = sys.argv[4].split(',') if len(sys.argv) > 4 else ['convert', 'unconvert']
        dep = int(sys.argv[5]) if len(sys.argv) > 5 else 2
        print('%d samples; executed per decode by the rate of the rarer control' % n)
        for stg in stages:
            tab = defaultdict(lambda: np.zeros(7))
            for i, (s, mk, c) in enumerate(seg):
                if c and marks[mk] and marks[mk][0] == stg:
                    sel = np.flatnonzero(so == i) if False else None
            # vectorised: per segment sums by bin
            segstage = np.array([1 if (marks[mk] and marks[mk][0] == stg) else 0 for (s, mk, c) in seg])
            m = segstage[so] == 1
            idx = np.flatnonzero(m)
            for b in range(7):
                sel = idx[binid[idx] == b]
                if len(sel) == 0:
                    continue
                sums = np.bincount(so[sel], weights=exo[sel], minlength=len(seg))
                for i in np.flatnonzero(sums):
                    k, _ = key(stacks[seg[i][0]], dep)
                    tab[k][b] += sums[i]
            allb = sum(tab.values())
            print('== %s: ' % stg + ' '.join('%s %.1f' % (BN[b], allb[b]) for b in range(7)) + ' | %.1f' % allb.sum())
            for k in sorted(tab, key=lambda k: -tab[k][:5].sum()):
                if tab[k][:5].sum() >= 5:
                    print('  %-50s ' % k[-50:] + ' '.join('%8.1f' % x for x in tab[k]))
        return
    if mode == 'wires':
        stg = sys.argv[4]
        mx = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-3
        top = int(sys.argv[6]) if len(sys.argv) > 6 else 60
        dep = int(sys.argv[7]) if len(sys.argv) > 7 else 3
        segstage = np.array([1 if (marks[mk] and marks[mk][0] == stg) else 0 for (s, mk, c) in seg])
        idx = np.flatnonzero((segstage[so] == 1) & (rate <= mx))
        rare = np.where(ca[idx] <= cb[idx], wa[idx], wb[idx])
        grp = defaultdict(lambda: [0, 0.0, 0, 0, None, None])
        for j, k in enumerate(idx):
            i = so[k]
            kk, depth = key(stacks[seg[i][0]], dep)
            g = grp[(kk, depth, int(rare[j]) if base[k] != 0x30 else -1)]
            g[0] += 1
            g[1] += exo[k]
            g[2] += int(min(ca[k], cb[k]))
            g[3] += int(ones[k])
            g[4] = k if g[4] is None else g[4]
            g[5] = k
        print('%d samples; stage %s; ops whose rarer control has rate <= %g: %d ops, %.1f executed per decode'
              % (n, stg, mx, len(idx), exo[idx].sum()))
        print('  ops   exec  flag hits/op  product hits/op  wire depth  first..last k   stack')
        for (kk, depth, w), g in sorted(grp.items(), key=lambda t: -t[1][1])[:top]:
            print('%5d %7.1f %9.1f %9.1f  q%-4d d%-2d %7d..%-7d %s'
                  % (g[0], g[1], g[2] / g[0], g[3] / g[0], w, depth, g[4], g[5], kk[-70:]))
        return


if __name__ == '__main__':
    main()
