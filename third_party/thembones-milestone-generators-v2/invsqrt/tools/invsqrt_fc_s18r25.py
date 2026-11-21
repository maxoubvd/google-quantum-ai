#!/usr/bin/env python3
"""invsqrt_fc_s18r25.py -- FC eighteenth holder: THE TOP OF R25 IN ROW 24's UNDO (exact degrees on ALL inputs).
Row 24's undo (skip E) ends by re-making R25's bits 7 and 8 (w7', w8': 1 AND each) so that the chain block's mirror
(ROW9's chain, block 37) can erase bits 9..4 by their one-AND forms.  If a top bit is NOT re-made, the chain mirror
must erase the bits above it by X-measurement and a phase over the wires live then.  This tool prints the exact GF(2)
degree (1, 2 = free; 3 = >= 0.5 a shot; >= 4 = >= 1.0) of R25's bit p (p = 9 = digit 24, 8, 7) as a function of
   A: v0..v6 and R25's bits 0..(p-1)          (every lower bit present)
   B: v0..v6 and R25's bits 0..6               (bits 7, 8 never re-made)
   C: v0..v6 and U's three top bits u1..u3     (after the chain mirror turned bits 0..3 back into U's bits)
(R54's integer model through the designer's tools/r65_inv_chainall_carry.chain, read-only import.)
  invsqrt_fc_s18r25.py SRC HTAB.npz
Safe to run: reads only, prints only."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import invsqrt_fc_s18debt as DB  # noqa: E402


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    hpath = os.path.expanduser(sys.argv[2])
    import r65_inv_chainall_carry as CT
    d, U = CT.chain(src, hpath)
    v = (d - 257) >> 1
    vb = [(v >> j) & 1 for j in range(7)]
    N = (np.int64(1) << 35) + (U << 8)
    R = N.copy()
    for b in (26, 25):
        neg = (R < 0)
        R = R - np.where(neg, -1, 1) * (d << b)
    W = (R >> 25) & 1023                         # R25's window, positions 25..34 (true value)
    wb = [(W >> p) & 1 for p in range(10)]
    ub = [((N >> (24 + j)) & 1) for j in range(4)]
    sg = (R < 0).astype(np.int64)
    print("R25: %d inputs; bit 9 == sign of R25: %s" % (len(d), bool(np.array_equal(wb[9], sg))))
    for p in (9, 8, 7):
        out = ["bit %d" % p]
        for nm, cols in (("A", vb + wb[:p]), ("B", vb + wb[:7] if p > 7 else vb + wb[:p]), ("C", vb + ub[1:])):
            rows = DB.fn_rows(cols, wb[p])
            if rows is None:
                out.append("%s: not a function" % nm)
            else:
                dg = DB.degree(rows, len(cols))
                out.append("%s: deg %s (%d patterns)" % (nm, ">=4" if dg > 3 else dg, len(rows)))
        print("  " + "  |  ".join(out), flush=True)


if __name__ == "__main__":
    main()
