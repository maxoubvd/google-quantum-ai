#!/usr/bin/env python3
"""dqi_sw_segslice.py -- S3 sweep: keep the real decode's segments of an r7prof segments.json.

With D1's a3 and DQI_PUB_A3_SYNC=1 a scratch decode of the unpatched build runs first in the same process, so
the r7prof patch records its segments too (its markers go into a discarded builder). The real decode's segments
are the tail whose static Toffolis sum to the built circuit's static count.
Usage: python dqi_sw_segslice.py <segments.json> <out.json> <static Toffolis of the built circuit>
"""
import json
import sys

src, out, static = sys.argv[1], sys.argv[2], int(sys.argv[3])
S = json.load(open(src))
seg = S['segments']
total = sum(x[2] for x in seg)
cum, cut = 0, None
for i, x in enumerate(seg):
    if total - cum == static:
        cut = i
        break
    cum += x[2]
assert cut is not None, 'no cut with %d static Toffolis after it (total %d)' % (static, total)
# move the cut back over zero-static segments that belong to the real decode's first marks (same stage tag)
S['segments'] = seg[cut:]
json.dump(S, open(out, 'w'))
print('kept segments %d of %d from %d; static %d of %d' % (len(seg) - cut, len(seg), cut, static, total))
