#!/usr/bin/env python3
"""r46_dqi_ext_inv17.py -- R46-DQI-EXT: the redundant AND of the GF(64) inverse at nonzero-input call sites of
extract / unextract, found on the marker file and carried to the line's published file.

In every call of the kit's `inverse6_joint18_12extra` (18 ANDs; hosts `inverse6_fixed18`, `inverse6_init18`) one CCX
has controls that are never both 0 when the inverse's input is nonzero: its product is a + b + 1 (CX a t, CX b t,
X t: zero Toffolis). The 12-bit inverse's input at the junctions is Norm(n) with n the chain's final denominator
(forward) or the carried scale (return): nonzero on every shot the line decodes (see research/panel/R46-DQI-EXT.md).

  map <marker.kmb> <segments_real.json> <line.kmb> <out prefix> <stat1.npz> [stat2.npz ...]
      the ops of the marker file with "a or b always 1" inside those kernels (stat files of r46_dqi_ext_scan.py),
      by class (junction = under external_node_chain; other); their ordinals among the Toffoli-class ops of the
      line's file, by aligning the two Toffoli sequences (base, a, b, t); writes <prefix>_A.npz (junction) and
      <prefix>_B.npz (all) in the format of tools/r44_xfer_dqidrop.py build2 (k, mode = 4).
  check <line.kmb> <cand.npz> <stat of the line file.npz>
      on the line file's own statistics: the listed ops execute on every sample and have a or b = 1 on all.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402


def toffolis(path):
    c = kmb.read(path)
    ops = np.asarray(c.ops)
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    t_op = np.flatnonzero(toff)
    base = (ops[t_op] & 0x7F).astype(np.int64)
    q = [np.asarray(x).astype(np.int64) for x in c.q3]
    n3 = np.cumsum(((ops & 0x30) >> 4) == 3) - 1
    i3 = n3[t_op]
    sig = (base << 48) | (q[0][i3] << 32) | (q[1][i3] << 16) | q[2][i3]
    return ops, t_op, sig


def main():
    mode = sys.argv[1]
    if mode == 'map':
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
        rel = np.flatnonzero((act == n) & (pa + pb - ones == act) & (ones > 0))
        pick = []
        for k in rel.tolist():
            sid, mkid, _n = segs[int(seg[k])]
            names = [x.split('@')[0] for x in st[sid]]
            if not marks[mkid] or marks[mkid][0] not in ('extract', 'unextract'):
                continue
            if 'inverse6_joint18_12extra' not in names:
                continue
            pick.append((k, marks[mkid][0], 'junction' if 'external_node_chain' in names else 'other',
                         ' > '.join(names[-4:])))
        print('marker file: %d samples; ops with a or b always 1 in inverse6_joint18_12extra of the block: %d '
              '(junction %d, other %d)' % (n, len(pick), sum(1 for p in pick if p[2] == 'junction'),
                                           sum(1 for p in pick if p[2] == 'other')))
        mo, mt, ms = toffolis(mk)
        lo, lt, ls = toffolis(line)
        print('Toffoli-class static: marker file %d, line file %d' % (len(ms), len(ls)))
        # greedy alignment with resync
        mp = np.full(len(ms), -1, np.int64)
        i = j = 0
        W = 6
        res = 0
        while i < len(ms) and j < len(ls):
            if ms[i] == ls[j]:
                # run of equal entries
                m = min(len(ms) - i, len(ls) - j)
                neq = np.flatnonzero(ms[i:i + m] != ls[j:j + m])
                r = int(neq[0]) if len(neq) else m
                mp[i:i + r] = np.arange(j, j + r)
                i += r
                j += r
                continue
            found = None
            for tot in range(1, 6000):
                cands = [(tot, 0), (0, tot)]              # a run deleted from / added to the line file
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
                print('no resync at marker Toffoli %d / line %d' % (i, j))
                break
            res += 1
            i += found[0]
            j += found[1]
        print('aligned %d of %d marker-file Toffolis to the line file (%d resyncs); line Toffolis matched %d of %d'
              % (int((mp >= 0).sum()), len(ms), res, int((mp >= 0).sum()), len(ls)))
        out = {'A': [], 'B': []}
        lost = 0
        for k, stg, cls, name in pick:
            if mp[k] < 0:
                lost += 1
                continue
            out['B'].append(int(mp[k]))
            if cls == 'junction':
                out['A'].append(int(mp[k]))
        print('picked ops present in the line file: %d of %d (junction %d)' % (len(out['B']), len(pick), len(out['A'])))
        for tag in ('A', 'B'):
            ks = np.array(sorted(out[tag]), np.int64)
            assert (lo[lt[ks]] & 0x7F == 0x31).all(), 'not all CCX'
            np.savez(pre + '_' + tag + '.npz', k=ks, mode=np.full(len(ks), 4, np.int64))
            print('written %s_%s.npz: %d CCX' % (pre, tag, len(ks)))
    elif mode == 'check':
        line, cand, statf = sys.argv[2:5]
        d = np.load(statf)
        ks = np.load(cand)['k']
        n = int(d['samples'])
        act, ones, pa, pb = (d[k][ks] for k in ('act', 'ones', 'pa', 'pb'))
        ok = (act == n) & (pa + pb - ones == act)
        print('%s on %d samples of the line file: %d ops; execute on every sample %d; a or b = 1 on every sample %d; '
              'P(product) %.3f..%.3f' % (os.path.basename(cand), n, len(ks), int((act == n).sum()), int(ok.sum()),
                                         (ones / np.maximum(act, 1)).min(), (ones / np.maximum(act, 1)).max()))


if __name__ == '__main__':
    main()
