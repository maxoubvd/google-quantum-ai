#!/usr/bin/env python3
"""r89_chem_cascade_pair.py (seat R89-CHEM-CASCADE): SUPPLY-AWARE pairs for a 23-output loader.

Copy of builder C4's tools/chem_c4_pair.py (read-only) with one change: the pairs of disjoint pool-free 6-groups are
ranked by the DEAD-CELL SUPPLY that their rest leaves to group 4 (the sum of the usable dead cells of the four
dead-richest outputs of the rest, from r89_chem_cascade_cells.py --json), and the output order puts those four LAST
(group 4), the other seven as group 3. The kill order is 't,X' with X = group 4's killed column (scanned later).
  python -I r89_chem_cascade_pair.py ns1.json ... --cells cells.json --nb 23 --max 8 --outdir pairs/ [--min 300]
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--cells", required=True, help="json {output: usable dead cells}")
    ap.add_argument("--nb", type=int, default=23)
    ap.add_argument("--g4", type=int, default=4)
    ap.add_argument("--g3first", type=int, default=0, help="richest outputs of the rest kept for group 3 first (a payer)")
    ap.add_argument("--max", type=int, default=8)
    ap.add_argument("--min", type=int, default=0, help="least group-4 supply written")
    ap.add_argument("--x", type=int, default=13, help="group 4's killed column written into args (scanned later)")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    cnt = {int(k): v for k, v in json.load(open(a.cells)).items()}
    finds = []
    for f in a.files:
        finds += json.load(open(f))
    byt = {}
    for x in finds:
        byt.setdefault(x["t"], []).append(x)
    os.makedirs(a.outdir, exist_ok=True)
    pairs = []
    for t, lst in byt.items():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                A, B = lst[i], lst[j]
                if set(A["outs"]) & set(B["outs"]):
                    continue
                rest = [o for o in range(a.nb) if o not in A["outs"] and o not in B["outs"]]
                srt = sorted(rest, key=lambda o: (-cnt.get(o, 0), o))
                g3f = srt[:a.g3first]                     # e.g. group 3's row-drop payer
                g4 = sorted(srt[a.g3first:a.g3first + a.g4])
                sup = sum(cnt.get(o, 0) for o in g4) + sum(cnt.get(o, 0) for o in g3f)
                pairs.append((sup, len(A["plans"]) * len(B["plans"]), t, A, B, rest, g4, g3f))
    pairs.sort(key=lambda p: (-p[0], -p[1]))
    seen, n = set(), 0
    for sup, _, t, A, B, rest, g4, g3f in pairs:
        if sup < a.min:
            break
        key = (t, tuple(A["outs"]), tuple(B["outs"]))
        if key in seen:
            continue
        seen.add(key)
        g3 = g3f + [o for o in rest if o not in g4 and o not in g3f]
        order = A["outs"] + B["outs"] + g3 + g4
        x = a.x if a.x != t else (a.x - 1)
        d = os.path.join(a.outdir, f"p{n}")
        os.makedirs(d, exist_ok=True)
        json.dump(A["plans"], open(os.path.join(d, "g0.json"), "w"))
        json.dump(B["plans"], open(os.path.join(d, "g1.json"), "w"))
        open(os.path.join(d, "args.txt"), "w").write(f"--order {','.join(map(str, order))} --kills {t},{x}\n")
        json.dump(dict(t=t, g4=g4, supply=sup), open(os.path.join(d, "r89.json"), "w"))
        print(f"p{n}: t {t} supply {sup} group 4 {g4} S1 {A['outs']} ({len(A['plans'])}) S2 {B['outs']} ({len(B['plans'])})")
        n += 1
        if n >= a.max:
            break
    print(f"{len(pairs)} disjoint pairs; {n} written (supply >= {a.min})")


if __name__ == "__main__":
    main()
