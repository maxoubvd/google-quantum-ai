#!/usr/bin/env python3
"""invsqrt_fw_m7chi.py -- F3 (twelfth holder, 2026-11-15): the m 7 rule's measured clear in the OTHER order (cfg h1c
"chiclr": 1): c2 is X-measured first and its phase repaid by C2G1 (c2 given c1, which is still live: 2 gates, 1 final
= 0.5 expected), then c1 by C1 (c1 alone). Programs: F1 format over [z0..z3 = Tf6..9 folded, l0 = L'7, l1 = L'8, c0]
(+ c1 as input 7 for C2G1), exact on the rule's 104 forward care patterns [M tools/invsqrt_fw_minand_fin.py, box 4
~/qcc/ws/invsqrt_fw/m7box/rule/]. Today's order (CHI6_C1G 1.0 + CHI6_C2 1.5) = 2.5."""

# c2 given c1: k 2, last gate final (k 2 r 2 UNSAT): 0.5
C2G1 = {"gates": [{"left": [1, 3, 7], "right": [2, 7]}, {"left": [1, 8], "right": [3]}], "outputs": [[9]]}
# c1 alone: k 5 with 2 final gates (k 4 UNSAT; k 5 r 4 UNSAT): 1.5
C1 = {"gates": [{"left": [0, 1, 2, 5], "right": [3]}, {"left": [1, 2, 3, 4, 6], "right": [7]},
                {"left": [2, 5, 7, 8], "right": [3, 5, 6, 8]}, {"left": [1, 5, 6, 7, 9], "right": [2, 3, 5, 7, 9]},
                {"left": [2, 5, 7, 9], "right": [3, 5, 7, 9]}], "outputs": [[3, 10, 11]]}


def check(care_c1, care_c2):
    """exactness of both programs on the care rows: care_* = (rows, inputs tt list, target tt) as in the SAT files."""
    import json
    out = {}
    for nm, prog, pf in (("C2G1", C2G1, care_c2), ("C1", C1, care_c1)):
        d = json.load(open(pf))
        rows = int(d["rows"])
        ins = [int(x["tt"], 16) for x in d["inputs"]]
        tg = int(d["target"], 16)
        n = len(ins)
        bad = 0
        for r in range(rows):
            v = [(t >> r) & 1 for t in ins]
            for g in prog["gates"]:
                a = 0
                for x in g["left"]:
                    a ^= v[x]
                b = 0
                for x in g["right"]:
                    b ^= v[x]
                v.append(a & b)
            o = 0
            for x in prog["outputs"][0]:
                o ^= 1 if x == -1 else v[x]
            bad += o != ((tg >> r) & 1)
        out[nm] = (rows, bad, n)
    return out


if __name__ == "__main__":
    import sys
    print(check(sys.argv[1], sys.argv[2]))
