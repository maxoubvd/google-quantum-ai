#!/usr/bin/env python3
"""dqi_int_pptag.py -- workstream dqi_int (B5): a re-roll of an emitted .kmb: the same gates plus one INFO_TAG op
at the end (no gate change; the file's bytes and so the verifier's draw change).
Usage: python dqi_int_pptag.py <in.kmb> <out.kmb> <label>
"""
import dataclasses
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

src, dst, label = sys.argv[1], sys.argv[2], sys.argv[3]
k = kmb.read(src)
out = dataclasses.replace(k, ops=np.append(np.asarray(k.ops, dtype=np.uint8), np.uint8(0x06)),
                          tags=list(k.tags) + ['reroll:' + label])
kw.write(dst, out)
print('tag appended: reroll:%s' % label)
