#!/usr/bin/env python3
"""r46_dqi_conv_mknosu2.py -- R46-DQI-CONV: copies of Z1's dqi_zn_pdl2.py and Z2's dqi_z2_pdr.py (the FROZEN ones of
B5's overlay) in which the paired blocks' carry routes leave out the sU2 term when R46_NOSU2=1:
    sU2 = (not B) and C and C2 and (not F)   with B, C, F, C2 = [q1 == 0] at levels U-1, U-2, U-3, U-4
i.e. the term that acts only when TWO quotient heads vanish in one block window (zero pivots at levels U-1 and U-3).
Removed with it: the flag's own AND chain and its erasure. Nothing else changes; with R46_NOSU2 unset the copies
emit the frozen circuit.
  python r46_dqi_conv_mknosu2.py <folder with the frozen dqi_zn_pdl2.py, dqi_z2_pdr.py> <out folder>
"""
import os
import sys

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)


def sub(s, old, new, n=1):
    assert s.count(old) == n, (old, s.count(old))
    return s.replace(old, new)


s = open(os.path.join(src, 'dqi_zn_pdl2.py')).read()
s = sub(s, "        sU2, lU2 = PC.flag_and(b, [(h1, False), (C2q, False), (Fq, True)])\n",
        "        NOSU2 = os.environ.get('R46_NOSU2') == '1'\n"
        "        if not NOSU2:\n"
        "            sU2, lU2 = PC.flag_and(b, [(h1, False), (C2q, False), (Fq, True)])\n")
s = sub(s, "            for j in range(12):\n                b.ccx(sU2, X[k][j], t2[k][j])\n",
        "            for j in range(0 if NOSU2 else 12):\n                b.ccx(sU2, X[k][j], t2[k][j])\n")
s = sub(s, "        PC.erase_flag_and(b, sU2, lU2)\n",
        "        if not NOSU2:\n            PC.erase_flag_and(b, sU2, lU2)\n")
open(os.path.join(out, 'dqi_zn_pdl2.py'), 'w').write(s)

s = open(os.path.join(src, 'dqi_z2_pdr.py')).read()
s = sub(s, "    sL, sU2 = flag_and(b, lL), flag_and(b, lU2)\n",
        "    NOSU2 = os.environ.get('R46_NOSU2') == '1'\n"
        "    sL = flag_and(b, lL)\n"
        "    sU2 = None if NOSU2 else flag_and(b, lU2)\n")
s = sub(s, "        if 'su2' not in RDROP:\n", "        if 'su2' not in RDROP and not NOSU2:\n")
s = sub(s, "    erase_flag_and(b, sU2, lU2)\n    erase_flag_and(b, sL, lL)\n",
        "    if not NOSU2:\n        erase_flag_and(b, sU2, lU2)\n    erase_flag_and(b, sL, lL)\n")
open(os.path.join(out, 'dqi_z2_pdr.py'), 'w').write(s)
print('written', out)
