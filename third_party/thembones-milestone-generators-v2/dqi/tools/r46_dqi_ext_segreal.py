#!/usr/bin/env python3
"""r46_dqi_ext_segreal.py -- R46-DQI-EXT: the real decode's segments of an r7prof segments.json, found from the
marker file itself (tools/dqi_sw_segslice.py needs the tail's static sum to equal the built count; on best_q23j5's
recipe several scratch decodes run in the same process, and inside the real decode some recorded segments have no
marker in the file, so no tail has that sum).

The built file holds M DEBUG_PRINT markers; per marker the Toffoli-class ops up to the next marker are counted in
the file. The json's segments are aligned to the markers from the head (the offset with the longest exact prefix)
and resynchronised at every mismatch by the smallest skip (json segments without a marker, or markers without a
json segment) after which the next 8 counts agree again. Writes segments_real.json: one entry per marker of the
file [stack id, marks id, Toffoli-class ops in the file]; markers that could not be paired get the stack and marks
of the nearest paired marker before them and are listed in 'r46_unpaired'.
Usage: python r46_dqi_ext_segreal.py <dqi.kmb> <segments.json> <segments_real.json>
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

c = kmb.read(sys.argv[1])
ops = np.asarray(c.ops)
mark = np.flatnonzero((ops & 0x7F) == 0x07)
toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
ct = np.concatenate(([0], np.cumsum(toff)))
M = len(mark)
edges = np.concatenate((mark, [len(ops)]))
per = (ct[edges[1:]] - ct[edges[:-1]]).astype(np.int64)
print('file: %d ops, %d Toffoli-class, %d markers, %d Toffoli-class before the first marker'
      % (len(ops), int(toff.sum()), M, int(ct[mark[0]]) if M else 0))
S = json.load(open(sys.argv[2]))
seg = S['segments']
st = np.array([x[2] for x in seg], np.int64)
print('segments.json: %d segments, static total %d' % (len(seg), int(st.sum())))
# head offset: longest exact prefix
best = (0, -1)
for off in range(max(0, len(seg) - M - 400), len(seg) - 50):
    n = min(M, len(seg) - off)
    neq = np.flatnonzero(st[off:off + n] != per[:n])
    pre = int(neq[0]) if len(neq) else n
    if pre > best[1]:
        best = (off, pre)
off = best[0]
print('head offset %d: exact prefix of %d markers' % best)
W = 8
i, j = off, 0
pair = [-1] * M
skips = []
while j < M and i < len(seg):
    if st[i] == per[j]:
        pair[j] = i
        i += 1
        j += 1
        continue
    found = None
    for tot in range(1, 400):
        for di in range(tot + 1):
            dj = tot - di
            a, b = i + di, j + dj
            n = min(W, len(seg) - a, M - b)
            if n > 0 and (st[a:a + n] == per[b:b + n]).all():
                found = (di, dj)
                break
        if found:
            break
    if not found:
        print('no resync at marker %d (json %d)' % (j, i))
        break
    skips.append((j, i, found[0], found[1]))
    i += found[0]
    j += found[1]
unp = [j for j in range(M) if pair[j] < 0]
print('paired %d of %d markers; resyncs %d; json segments skipped %d; unpaired markers %d holding %d Toffoli-class ops'
      % (M - len(unp), M, len(skips), sum(x[2] for x in skips), len(unp), int(per[unp].sum()) if unp else 0))
for j, i, di, dj in skips[:12]:
    s, mk, n = seg[i]
    print('  resync at marker %d: skipped %d json segments, %d markers; marks %s; leaf %s'
          % (j, di, dj, '/'.join(S['marks'][mk]), ' > '.join(x.split('@')[0] for x in S['stacks'][s][-3:])))
print('json index of the last marker %d (json has %d)' % (pair[-1], len(seg)))
out = []
last = off
for j in range(M):
    if pair[j] >= 0:
        last = pair[j]
    s, mk, _n = seg[last]
    out.append([s, mk, int(per[j])])
S['segments'] = out
S['r46_unpaired'] = unp
S['r46_note'] = 'one entry per marker of the file; head offset %d; static counts from the file' % off
json.dump(S, open(sys.argv[3], 'w'))
by = {}
for s, mk, n in out:
    k = S['marks'][mk][0] if S['marks'][mk] else '-'
    by[k] = by.get(k, 0) + n
print('static by first mark: ' + ', '.join('%s %d' % kv for kv in by.items()))
print('written %s' % sys.argv[3])
