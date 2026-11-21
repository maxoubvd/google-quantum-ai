#!/usr/bin/env python3
"""invsqrt_fc2_chi.py -- FC2 (own code): the upper carry bits Chi = C >> 1 of H1' (m 7) as functions of a few bits,
before (pre: folded T, s) and after (post: folded T', s') the increment T += C; prints the smallest supports and,
for the chosen one, the truth table classes.  invsqrt_fc2_chi.py 'SPEC' [Y_LO Y_HI]"""
import json
import sys
import os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r37_inv_40_ref as R37

spec = dict(R37.DEF7)
if len(sys.argv) > 1 and sys.argv[1].startswith("{"):
    spec.update(json.loads(sys.argv[1]))
lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (-523202, 521230)
nb = int(spec["nb"])
y = np.arange(lo, hi + 1, dtype=np.int64)
U, C = R37.stage_ref(y, spec)
s = (y < 0).astype(np.int64)
T = y >> 9
Tf = np.where(s == 1, ~T, T) & ((1 << nb) - 1)
L2 = U & 511
c0 = C & 1
chi = C >> 1
T2 = U >> 9
s2 = (U < 0).astype(np.int64)
Tf2 = np.where(s2 == 1, ~T2, T2) & ((1 << nb) - 1)
print(json.dumps(dict(C_vals=sorted(set(C.tolist())), chi_vals=sorted(set(chi.tolist())))))


def supports(tf, sg, name):
    out = []
    for kt in range(2, nb + 1):
        for kl in range(0, 6):
            for use_s in (0, 1):
                key = (tf >> (nb - kt)) * 1024 + ((L2 >> (9 - kl)) if kl else 0) * 4 + c0 * 2 + (sg if use_s else 0)
                o = np.argsort(key, kind="stable")
                ks, vs = key[o], chi[o]
                brk = np.flatnonzero(np.diff(ks)) + 1
                mn = np.minimum.reduceat(vs, np.r_[0, brk])
                mx = np.maximum.reduceat(vs, np.r_[0, brk])
                if (mn == mx).all():
                    out.append((kt + kl + 1 + use_s, kt, kl, use_s))
    out.sort()
    print(name, "smallest supports (total bits, Tf top, L'' top, sign used):", out[:6])
    return out


sp = supports(Tf, s, "PRE ")
so = supports(Tf2, s2, "POST")
for nm, tf, sg, best in (("PRE", Tf, s, sp[0]), ("POST", Tf2, s2, so[0])):
    _, kt, kl, us = best
    z = tf >> (nb - kt)
    l = (L2 >> (9 - kl)) if kl else np.zeros_like(L2)
    tab = {}
    for zz, ll, cc, ss, v in zip(z.tolist(), l.tolist(), c0.tolist(), (sg if us else 0 * sg).tolist(), chi.tolist()):
        tab[(cc, ss, ll, zz)] = v
    print(nm, "kt", kt, "kl", kl, "sign", us, "entries", len(tab))
    for cc in (0, 1):
        for ss in ((0, 1) if us else (0,)):
            for ll in range(1 << kl):
                row = "".join(str(tab.get((cc, ss, ll, zz), ".")) for zz in range(1 << kt))
                print("  c0 %d s %d l %d : %s" % (cc, ss, ll, row))


def thresholds(c0, l1, l0):
    """FC2's reading of the tables (19:30Z): Chi = [z >= t0] + [z >= t1] + [z >= t2], z = top 4 folded bits."""
    p = l1 & l0
    t0 = 6 + 3 * c0 + p
    t1 = 11 + 2 * c0 + np.where(c0 == 1, p, l1)
    t2 = np.where(c0 == 1, 16, 14 + l1)
    return t0, t1, t2


for nm, tf in (("PRE", Tf), ("POST", Tf2)):
    z = tf >> (nb - 4)
    l1, l0 = (L2 >> 8) & 1, (L2 >> 7) & 1
    t0, t1, t2 = thresholds(c0, l1, l0)
    got = (z >= t0).astype(np.int64) + (z >= t1) + (z >= t2)
    print(nm, "threshold rule mismatches on all y:", int((got != chi).sum()), "of", len(y))
