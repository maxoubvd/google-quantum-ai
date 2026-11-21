#!/usr/bin/env python3
"""r85_chem_budget_t1_gate: MAXEND (pre-hold END, mHa^2) for a ticket's count a pass, from MEASURED fixed-word prices in
pB1's END frame (END 399.34, count 3,511). Plain python3, no numpy: drop-in for chem_t10_s4fin.sh L13:
  MAXEND=${MAXEND:-$(python3 <HOME>/Quantum/tools/r85_chem_budget_t1_gate.py ${CNT#E})}
Slopes, loader Toffolis a pass a mHa^2 [M box 5 refits, R85-CHEM-BUDGET + -T1]: room side D30 3.16 on [398.255, 399.34],
then the mean of D30 and E30 (2.43) to 397.15 (g10k1's place in this frame), then E30 1.70; over side B30 6.29 on
[399.34, 400.294], J30 8.63 on [400.294, 401.369], 8.63 again to 401.87 (unmeasured: the lenient end), H30 27.8 beyond."""
import sys
E0, BASE = 399.34, 3511
OVER = [(0.954, 6.29), (1.075, 8.63), (0.501, 8.63), (1e9, 27.8)]       # (width, slope) going up from E0
ROOM = [(1.085, 3.16), (1.105, 2.43), (1e9, 1.70)]                      # going down from E0
def line(count):
    s = BASE - count; x = E0
    for w, k in (OVER if s > 0 else ROOM):
        need = abs(s)
        if need <= w * k:
            return x + (need / k if s > 0 else -need / k)
        s = (abs(s) - w * k) * (1 if s > 0 else -1); x += w if s > 0 else -w
    return x
if __name__ == '__main__':
    print('%.3f' % line(int(sys.argv[1])))
