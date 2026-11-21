#!/usr/bin/env python3
"""inv_i39_slotmodel.py -- I39 (builder): BIT-EXACT CLASSICAL MODEL of the slot stage's arithmetic (step 2 of the
recipe), checked on ALL inputs: the register U = [slot, y19 bits] = 2 y19; the passes add hn(T) (T = y19 >> 8, a
table of T: today's built h0 plus the table's moves); then the ROUNDING: u = slot (= hn mod 2), cy = u AND NOT t
(t = T_0 xor par), the slot is set to t (u leaves by measurement: its phase is repaid in the mirror), the part above
the slot gets + cy; then the CARRY RULE (emit_chi's thresholds) gives C >> 1, c0 = C & 1, T += C; then c0 is CLEARED
by c0 xor U9 xor U0 xor par. MIRROR: T -= C (c0 re-made from U), the passes SUBTRACTED with the slot = t: the register
must then hold 2 y19 + b (b = u xor t on the slot, y19 above, no borrow).  Checks, each on all 3,133,452 inputs:
  (1) U_forward == 2 y19 + hp(T) (hp = hn + ((hn xor t) & 1));  (2) the carry rule == (low + hp) >> 9 >> 1;
  (3) c0 == U9 xor U0 xor par after T += C;  (4) the mirror's register == 2 y19 + b and its slot == b;
  (5) U injective (no two inputs share U).
  inv_i39_slotmodel.py SRC WD HTAB 'SCHED1' 'SCHED2' PAR TABLE.json|zero [CLO]
Safe: prints only."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r37_inv_40_model as RM          # noqa: E402  (reads sys.argv[1] = SRC)
from inv_i39_chicheck import chi      # noqa: E402


def main():
    hpath = os.path.expanduser(sys.argv[3])
    s1, s2, par = json.loads(sys.argv[4]), json.loads(sys.argv[5]), int(sys.argv[6])
    tab = sys.argv[7]
    clo = int(sys.argv[8]) if len(sys.argv) > 8 else 4
    hz = np.load(hpath)
    h0 = hz["h"].astype(np.int64)
    seen = hz["seen"].astype(bool)
    idx = np.arange(2048)
    good = idx[seen]
    near = good[np.argmin(np.abs(idx[:, None] - good[None, :]), axis=1)]
    h0 = np.where(seen, h0, h0[near])
    Tt = np.arange(-1024, 1024)
    s = (Tt < 0).astype(np.int64)
    Tf = np.where(s == 1, ~Tt, Tt) & 1023
    if tab == "zero":
        hn = h0
    else:
        D = json.load(open(tab))["S1"]["tables"]
        dw, dc = np.array(D[0]), np.array(D[1])
        hn = h0 - dw[(Tf >> 1) & 127] + np.where(s == 1, -1, 1) * dc[Tf >> clo]
    t_of = (Tt & 1) ^ par
    hp_tab = hn + ((hn ^ t_of) & 1)
    _, _, y, _ = RM.chain(7, 4, None, s1, s2)
    T = y >> 8
    ti = T + 1024
    hnv, t, hpv = hn[ti], t_of[ti], hp_tab[ti]
    # forward: register R = 2 y + hn (the passes), as a 20-bit two's-complement word
    R = 2 * y + hnv
    u = R & 1
    cy = u & (1 - t)
    R = (R - u + t) + 2 * cy                          # slot := t, the part above the slot + cy
    ok1 = np.array_equal(R, 2 * y + hpv)
    # carry rule on the low 9 bits (before T += C): z from T (folded), l1 l0 from the low part after the table
    low = (2 * y) & 511
    Ls = low + hpv
    C = Ls >> 9
    Lpp = Ls & 511
    Tf_in = np.where(T < 0, ~T, T) & 1023
    z = (Tf_in >> 6) & 15
    ch = chi(z, (Lpp >> 7) & 1, (Lpp >> 8) & 1, C & 1)
    ok2 = bool(np.array_equal(ch, C >> 1))
    # U after T += C equals R (R already holds 2 y + hp: its bits 9.. are T + C by construction)
    U = R
    U9 = (U >> 9) & 1
    U0 = U & 1
    ok3 = bool(np.array_equal(C & 1, U9 ^ U0 ^ par))
    # mirror: subtract hn (the passes) from U with the slot = t
    M = U - hnv
    b = (hnv ^ t) & 1
    ok4 = bool(np.array_equal(M, 2 * y + b)) and bool(np.array_equal(M & 1, b))
    ok5 = len(np.unique(U)) == len(np.unique(y))        # the stage maps y -> U: injective on the distinct y
    print("table %s par %d: (1) forward == 2y + hp: %s; (2) carry rule exact: %s; (3) c0 = U9^U0^par: %s; "
          "(4) mirror leaves 2y + b with b on the slot: %s; (5) U injective: %s; rounding up on %d inputs, slot "
          "carried (cy = 1) on %d" % (os.path.basename(tab), par, ok1, ok2, ok3, ok4, ok5, int((hpv != hnv).sum()),
                                      int(cy.sum())), flush=True)


if __name__ == "__main__":
    main()
