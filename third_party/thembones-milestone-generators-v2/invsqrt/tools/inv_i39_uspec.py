#!/usr/bin/env python3
"""inv_i39_uspec.py -- I39: the slot's parity bit u(T) = hn(T) mod 2 as a family-SAT spec (FC's famg format) over the
stage's own wires at the rounding point: inputs x0..x9 = Tf0..9 (T folded by its sign), x10 = s (the sign); one output;
the T that never occur (2,048 - 2,041) are don't-cares. For the 'exact clear' rounding: slot ^= u(T) (computed), then
slot ^= t: the register becomes 2 y19 + hn - u + t (no carry, no held bit).   inv_i39_uspec.py HTAB.npz OUT.json [TABLE]
Safe: writes OUT only."""
import json
import sys

import numpy as np

hz = np.load(sys.argv[1])
h0 = hz["h"].astype(np.int64)
seen = hz["seen"].astype(bool)
Tt = np.arange(-1024, 1024)
s = (Tt < 0).astype(np.int64)
Tf = np.where(s == 1, ~Tt, Tt) & 1023
hn = h0.copy()
if len(sys.argv) > 3:
    D = json.load(open(sys.argv[3]))["S1"]["tables"]
    hn = h0 - np.array(D[0])[(Tf >> 1) & 127] + np.where(s == 1, -1, 1) * np.array(D[1])[Tf >> 4]
vals = [0] * 2048
alt = {}
for i in range(2048):
    r = int(Tf[i]) | (int(s[i]) << 10)
    vals[r] = int(hn[i] & 1)
    if not seen[i]:
        alt[str(r)] = [1 - vals[r]]
json.dump({"vals": vals, "alt": alt, "nb": 1, "note": "u(T) over Tf0..9, s; unseen T free"}, open(sys.argv[2], "w"))
print("%s: 2048 rows, %d free (unseen T), ones %d" % (sys.argv[2], len(alt), sum(vals)))
