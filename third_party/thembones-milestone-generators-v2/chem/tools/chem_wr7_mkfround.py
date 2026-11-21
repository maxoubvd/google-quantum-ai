"""chem_wr7_mkfround (WR7): writes tools/chem_wr7_froundl.py = R50-CHEM-LBL's third-hold run file
tools/r50_chem_lbl_fround.py (FreezePlacer + the bit-0 block code, --fp --fp-code 13; it carries WR7's guard) + R8's
in-stage LM lever: the same four exact-string insertions that chem_wr7_mkround.py puts into chem_wr7_roundl.py.
Asserts (1) WR7's guard is in fround.py byte for byte, (2) every lever anchor is found exactly once.
    python chem_wr7_mkfround.py
"""
import os
import re

T = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(T, "r50_chem_lbl_fround.py")).read()

# (1) the guard of chem_wr7_roundg.py, copied from chem_wr7_mkround.py's GUARD (new text)
mk = open(os.path.join(T, "chem_wr7_mkround.py")).read()
g0 = mk.index('GUARD = [')
g1 = mk.index(']\nfor old, new in GUARD:')
ns = {}
exec(mk[g0:g1 + 1], ns)
guard_new = ns["GUARD"][0][1]
assert src.count(guard_new) == 1, "WR7's guard is NOT in r50_chem_lbl_fround.py byte for byte"
print("guard: found once in r50_chem_lbl_fround.py")

# (2) the lever's four insertions, taken from chem_wr7_mkround.py's LEVER list (same text)
l0 = mk.index('LEVER = [')
l1 = mk.index(']\nsrc_l = src_g')
ns2 = {}
exec(mk[l0:l1 + 1], ns2)
for old, new in ns2["LEVER"]:
    n = src.count(old)
    assert n == 1, ("lever anchor found %d times" % n, old[:90])
    src = src.replace(old, new)
head = '"""r50_chem_lbl_fround = '
assert src.count(head) == 1, "docstring head not found once"
src = src.replace(head, '"""chem_wr7_froundl (WR7) = r50_chem_lbl_fround.py + R8\'s in-stage LM lever (R8 holder 6; '
                  'chem_wr7_mkfround.py).\nr50_chem_lbl_fround = ', 1)
open(os.path.join(T, "chem_wr7_froundl.py"), "w").write(src)
print("written tools/chem_wr7_froundl.py (%d lever insertions)" % len(ns2["LEVER"]))
