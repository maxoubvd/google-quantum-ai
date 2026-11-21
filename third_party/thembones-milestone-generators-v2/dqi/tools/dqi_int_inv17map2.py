#!/usr/bin/env python3
"""dqi_int_inv17map2.py -- workstream dqi_int (B5, Nov 15 twelfth holder): INV17's sites (R46-DQI-EXT, credit) carried
to a file of ANOTHER recipe through that recipe's own marker file.

R46-DQI-EXT found the sites on its marker file of the line (prof1): position 6 of the 9-AND `rankone_xor` segment of
`inverse6_joint18_12extra` holds a or b = 1 on every sample of the stat files (tools/r46_dqi_ext_s4inv17.py's rule).
On a recipe whose converter allocates differently (zpark, sU2's drop) the line's marker file does not align to the
new file past convert. Here the calls are matched by ORDER: the k-th call of the kernel in stage s of prof1 is the
k-th call in stage s of the new recipe's marker file (the counts per stage and the frames above the kernel must agree,
else the stage is refused); the new marker file is then aligned to the new target file (same allocation: the
alignment runs through) as in R46-DQI-EXT's map.
Usage: python dqi_int_inv17map2.py <prof1 marker.kmb> <prof1 segments_real.json> <new marker.kmb>
           <new segments_real.json> <target.kmb> <out prefix> <stat.npz of prof1> [...]
Writes <prefix>_A.npz (junction sites of extract / unextract: the 28 proved), <prefix>_B.npz (extract / unextract),
<prefix>_C.npz (all stages), in tools/r44_xfer_dqidrop.py build2's format (k, mode = 4).
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r46_dqi_ext_inv17 import toffolis  # noqa: E402
from dqi_int_inv17map import align  # noqa: E402


def calls(path, segf):
    """Per stage, in order: (static Toffoli index of position 6, frames above the kernel, junction?) of every call."""
    S = json.load(open(segf))
    marks, st, segs = S['marks'], S['stacks'], S['segments']
    ops, t_op, sig = toffolis(path)
    mark = np.flatnonzero(ops == 0x07)
    seg = np.searchsorted(mark, t_op, side='right') - 1
    first = {}
    out = {}
    for k in range(len(seg)):
        s = int(seg[k])
        if s not in first:
            first[s] = k
        if s < 0 or k - first[s] != 6:
            continue
        sid, mkid, cnt = segs[s]
        names = [x.split('@')[0] for x in st[sid]]
        if cnt != 9 or len(names) < 2 or names[-1] != 'rankone_xor' or names[-2] != 'inverse6_joint18_12extra':
            continue
        stage = marks[mkid][0] if marks[mkid] else '-'
        out.setdefault(stage, []).append((k, ' > '.join(names[-6:-2]), 'external_node_chain' in names))
    return out, sig


def main():
    pm, pseg, nm, nseg, tgt, pre = sys.argv[1:7]
    act = None
    n = 0
    for f in sys.argv[7:]:
        d = np.load(f)
        if act is None:
            act, ones, pa, pb = (d[k].astype(np.int64).copy() for k in ('act', 'ones', 'pa', 'pb'))
        else:
            act += d['act']
            ones += d['ones']
            pa += d['pa']
            pb += d['pb']
        n += int(d['samples'])
    hold = (act == n) & (pa + pb - ones == act) & (ones > 0)
    pc, _ = calls(pm, pseg)
    nc, nsig = calls(nm, nseg)
    print('stat samples %d; calls per stage prof1 / new: %s' % (
        n, ', '.join('%s %d/%d' % (s, len(pc[s]), len(nc.get(s, []))) for s in pc)))
    pick = []
    for s in pc:
        if len(pc[s]) != len(nc.get(s, [])):
            print('stage %s refused: call counts differ' % s)
            continue
        bad = sum(1 for a, b in zip(pc[s], nc[s]) if a[1] != b[1] or a[2] != b[2])
        if bad:
            print('stage %s refused: %d calls with other frames' % (s, bad))
            continue
        for a, b in zip(pc[s], nc[s]):
            if hold[a[0]]:
                pick.append((b[0], s, a[2]))
    blk = [p for p in pick if p[1] in ('extract', 'unextract')]
    print('picked on the new marker file: %d (extract + unextract %d, junction %d)' % (
        len(pick), len(blk), sum(1 for p in blk if p[2])))
    to, tt, ts = toffolis(tgt)
    mp, res, stop = align(nsig, ts)
    print('new marker file -> target: %d of %d aligned (%d resyncs, stop %s); target Toffolis %d' % (
        int((mp >= 0).sum()), len(nsig), res, stop, len(ts)))
    out = {'A': [], 'B': [], 'C': []}
    lost = {}
    for k, s, junc in pick:
        t = int(mp[k])
        if t < 0:
            lost[s] = lost.get(s, 0) + 1
            continue
        out['C'].append(t)
        if s in ('extract', 'unextract'):
            out['B'].append(t)
            if junc:
                out['A'].append(t)
    print('mapped: C %d, B %d, A %d; lost by stage %s' % (len(out['C']), len(out['B']), len(out['A']), lost))
    for tag in ('A', 'B', 'C'):
        ks = np.array(sorted(set(out[tag])), np.int64)
        assert (to[tt[ks]] & 0x7F == 0x31).all(), 'not all CCX'
        np.savez(pre + '_' + tag + '.npz', k=ks, mode=np.full(len(ks), 4, np.int64))
        print('written %s_%s.npz: %d CCX' % (pre, tag, len(ks)))


if __name__ == '__main__':
    main()
