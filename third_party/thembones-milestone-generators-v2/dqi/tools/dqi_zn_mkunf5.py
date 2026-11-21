#!/usr/bin/env python3
"""dqi_zn_mkunf5.py -- builder Z1 (Nov 14, fourth holder): make kv_unfused5 from a copy of Z2's kv_unfused4.
kv_unfused5 = kv_unfused4 without a14's reverse fusion (research/agents/a14/reverse_fusion.py); a12's
reverse_root_lifetime (which a14's install also installs) is kept. One edit in research/production/selected.py.
Usage: python3 dqi_zn_mkunf5.py <copy of kv_unfused4 to patch in place>
Measured on the line's converter (tools/dqi_zn_lct.py): +599.8 per pair against kv_unfused4 alone; with flagpark part
h and tools/dqi_zn_kpark.py it is 878 / 878, exact; with tools/dqi_zn_pdl2.py's reverse pairing see research/ws_dqi_zn.md.
"""
import os
import sys

root = sys.argv[1]
p = os.path.join(root, 'research/production/selected.py')
s = open(p).read()
old = "reverse_fusion.install(cv)\n"
new = ("import reverse_root_lifetime  # Z1 kv_unfused5: a14 reverse fusion NOT installed, a12 kept\n"
       "reverse_root_lifetime.install(cv)\n")
assert s.count(old) == 1, 'a14 install line not found once'
s = s.replace(old, new)
open(p, 'w').write(s)
print('kv_unfused5 patched in', root)
