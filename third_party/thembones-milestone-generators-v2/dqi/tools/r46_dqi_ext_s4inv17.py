#!/usr/bin/env python3
"""r46_dqi_ext_s4inv17.py -- R46-DQI-EXT session 4: INV17 (the redundant AND of the GF(64) inverse at call sites
whose input is nonzero) over the WHOLE decoder: r46_dqi_ext_inv17.py's map without the stage filter.

  map <marker.kmb> <segments_real.json> <file.kmb> <out prefix> <stat1.npz> [...]
      the ops of the marker file inside inverse6_joint18_12extra that execute on every sample and have a or b = 1
      on every sample of ALL the stat files given (pass the verifier-law files and the uniform-law files together);
      printed by stage and innermost frames, with the calls of the same kernel that do NOT hold the relation;
      aligned to <file.kmb> as in r46_dqi_ext_inv17.py; writes <prefix>_B.npz (extract + unextract) and
      <prefix>_C.npz (all stages) in the format of tools/r44_xfer_dqidrop.py build2 (k, mode = 4).
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r46_dqi_ext_inv17 import toffolis  # noqa: E402


def main():
    assert sys.argv[1] == 'map'
    mk, segf, line, pre = sys.argv[2:6]
    S = json.load(open(segf))
    marks, st, segs = S['marks'], S['stacks'], S['segments']
    act = None
    n = 0
    for f in sys.argv[6:]:
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
    first = {}
    pick = []
    groups = {}
    for k in range(len(seg)):
        s = int(seg[k])
        if s not in first:
            first[s] = k
        sid, mkid, cnt = segs[s]
        names = [x.split('@')[0] for x in st[sid]]
        if cnt != 9 or len(names) < 2 or names[-1] != 'rankone_xor' or names[-2] != 'inverse6_joint18_12extra':
            continue
        if k - first[s] != 6:
            continue
        stage = marks[mkid][0] if marks[mkid] else '-'
        hold = bool(act[k] == n and pa[k] + pb[k] - ones[k] == act[k] and ones[k] > 0)
        p0 = 1 - (pa[k] + pb[k] - ones[k]) / max(1, act[k])
        g = groups.setdefault((stage, hold, ' > '.join(names[-6:-2])), [0, 1.0, 0.0, 1.0])
        g[0] += 1
        g[1] = min(g[1], p0)
        g[2] = max(g[2], p0)
        g[3] = min(g[3], act[k] / n)
        if hold:
            pick.append((k, stage))
    print('marker file: %d samples; calls of inverse6_joint18_12extra %d; position-6 AND with a or b = 1 on every '
          'sample: %d (extract + unextract %d)' % (n, sum(g[0] for g in groups.values()), len(pick),
                                                   sum(1 for p in pick if p[1] in ('extract', 'unextract'))))
    print('stage / relation / calls / P(a = b = 0) min..max / least executed share / frames above the kernel')
    for (stage, hold, name), g in sorted(groups.items(), key=lambda x: (x[0][0], not x[0][1], -x[1][0])):
        print('  %-18s %-5s %3d  %.4f..%.4f  %.3f  %s' % (stage, 'holds' if hold else 'FAILS', g[0], g[1], g[2], g[3], name))
    mo, mt, ms = toffolis(mk)
    lo, lt, ls = toffolis(line)
    print('Toffoli-class static: marker file %d, target file %d' % (len(ms), len(ls)))
    mp = np.full(len(ms), -1, np.int64)
    i = j = 0
    W = 6
    res = 0
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
            print('no resync at marker Toffoli %d / target %d' % (i, j))
            break
        res += 1
        i += found[0]
        j += found[1]
    print('aligned %d of %d marker-file Toffolis to the target file (%d resyncs)' % (int((mp >= 0).sum()), len(ms), res))
    out = {'B': [], 'C': []}
    lost = {}
    for k, stage in pick:
        if mp[k] < 0:
            lost[stage] = lost.get(stage, 0) + 1
            continue
        out['C'].append(int(mp[k]))
        if stage in ('extract', 'unextract'):
            out['B'].append(int(mp[k]))
    print('picked ops present in the target file: %d of %d (block %d); lost by stage %s'
          % (len(out['C']), len(pick), len(out['B']), lost))
    for tag in ('B', 'C'):
        ks = np.array(sorted(out[tag]), np.int64)
        assert len(set(ks.tolist())) == len(ks)
        assert (lo[lt[ks]] & 0x7F == 0x31).all(), 'not all CCX'
        np.savez(pre + '_' + tag + '.npz', k=ks, mode=np.full(len(ks), 4, np.int64))
        print('written %s_%s.npz: %d CCX' % (pre, tag, len(ks)))


if __name__ == '__main__':
    main()
