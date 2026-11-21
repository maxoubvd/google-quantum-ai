#!/usr/bin/env python3
"""inv_i39_chicheck.py -- I39: does the stage's CARRY RULE (emit_chi: C >> 1 = [z >= t0] + [z >= t1] + [z >= t2],
t0 = 6 + 3 c0 + l1 l0, t1 = 11 + 2 c0 + (c0 ? l1 l0 : l1), t2 = c0 ? 16 : 14 + l1; z = T's top 4 folded bits, l1 l0 =
bits 8, 7 of the low part after the table, c0 = C & 1) still hold for a SLOT table? Model level, all inputs:
  control: today's 20-digit chain (F1 5, today's schedules), U = y + h0(T), T = y >> 9 (must be exact);
  slot:    the 19-digit chain (F1 4, SCHED1 / SCHED2), U = 2 y19 + hp(T), T = y19 >> 8, hp = the table + round-up.
  inv_i39_chicheck.py SRC WD HTAB 'SCHED1' 'SCHED2' PAR TABLE.json|zero [CLO]
Safe: prints only."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r37_inv_40_model as RM          # noqa: E402  (reads sys.argv[1] = SRC)


def chi(z, l0, l1, c0):
    t0 = 6 + 3 * c0 + (l1 & l0)
    t1 = 11 + 2 * c0 + np.where(c0 == 1, l1 & l0, l1)
    t2 = np.where(c0 == 1, 16, 14 + l1)
    return (z >= t0).astype(np.int64) + (z >= t1) + (z >= t2)


def check(y, hp_of_T, shift, tag):
    T = y >> shift
    s = (T < 0).astype(np.int64)
    Tf = np.where(s == 1, ~T, T) & 1023
    z = (Tf >> 6) & 15
    h = hp_of_T[T + 1024]
    low = (y if shift == 9 else 2 * y) & 511
    Ls = low + h
    C = Ls >> 9
    Lpp = Ls & 511
    l0, l1, c0 = (Lpp >> 7) & 1, (Lpp >> 8) & 1, C & 1
    ch = chi(z, l0, l1, c0)
    bad = int((ch != (C >> 1)).sum())
    print("%s: inputs %d, C range %d..%d, carry rule wrong on %d inputs" % (tag, len(y), C.min(), C.max(), bad),
          flush=True)
    return bad


def main():
    wd, hpath = os.path.expanduser(sys.argv[2]), os.path.expanduser(sys.argv[3])
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
    # control: today's chain (F1 5, k9f's schedules as in R49's calibration: sched1 [0,0,0,2,3,4,5,6], sched2 [0,0,1,2,3,4,6,6])
    _, _, y20, _ = RM.chain(7, 5, None, [0, 0, 0, 2, 3, 4, 5, 6], [0, 0, 1, 2, 3, 4, 6, 6])
    check(y20, h0, 9, "control (today's 20-digit chain, today's table)")
    _, _, y19, _ = RM.chain(7, 4, None, s1, s2)
    T = np.arange(-1024, 1024)
    s = (T < 0).astype(np.int64)
    Tf = np.where(s == 1, ~T, T) & 1023
    if tab == "zero":
        hn = h0
    else:
        D = json.load(open(tab))["S1"]["tables"]
        dw, dc = np.array(D[0]), np.array(D[1])
        hn = h0 - dw[(Tf >> 1) & 127] + np.where(s == 1, -1, 1) * dc[Tf >> clo]
    t = (T & 1) ^ par
    hp = hn + ((hn ^ t) & 1)
    check(y19, hp, 8, "slot (19-digit chain, %s, par %d)" % (os.path.basename(tab), par))


if __name__ == "__main__":
    main()
