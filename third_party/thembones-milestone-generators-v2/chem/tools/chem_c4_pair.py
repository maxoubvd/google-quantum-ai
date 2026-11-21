#!/usr/bin/env python3
"""chem_c4_pair.py: builder C4. From chem_c4_ns.py finds (pool-free 6-groups at 38, every product non-silent on a
coordinate t), pick pairs of DISJOINT output subsets on the same t and write build inputs: plans for group 1 and
group 2, the output order (S1, S2, then the other outputs: 7 for group 3, 6 for group 4), the kill order (t first:
group 3's column kill; then the engine's next coordinate for group 4).

  python -I chem_c4_pair.py ns1.json ns2.json ... --nb 25 --t2 14 --max 8 --outdir pairs/
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--nb", type=int, default=25)
    ap.add_argument("--t2", default="14,13", help="group 4's kill coordinate: the first of these that is not t")
    ap.add_argument("--max", type=int, default=8)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    finds = []
    for f in a.files:
        finds += json.load(open(f))
    byt = {}
    for x in finds:
        byt.setdefault(x["t"], []).append(x)
    os.makedirs(a.outdir, exist_ok=True)
    t2s = [int(v) for v in a.t2.split(",")]
    pairs = []
    for t, lst in sorted(byt.items(), key=lambda kv: -len(kv[1])):
        lst = sorted(lst, key=lambda x: -len(x["plans"]))
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                A, B = lst[i], lst[j]
                if set(A["outs"]) & set(B["outs"]):
                    continue
                pairs.append((len(A["plans"]) * len(B["plans"]), t, A, B))
    pairs.sort(key=lambda p: -p[0])
    seen, n = set(), 0
    for _, t, A, B in pairs:
        key = (t, tuple(A["outs"]), tuple(B["outs"]))
        if key in seen:
            continue
        seen.add(key)
        rest = [o for o in range(a.nb) if o not in A["outs"] and o not in B["outs"]]
        order = A["outs"] + B["outs"] + rest
        t2 = [v for v in t2s if v != t][0]
        d = os.path.join(a.outdir, f"p{n}")
        os.makedirs(d, exist_ok=True)
        json.dump(A["plans"], open(os.path.join(d, "g0.json"), "w"))
        json.dump(B["plans"], open(os.path.join(d, "g1.json"), "w"))
        open(os.path.join(d, "args.txt"), "w").write(f"--order {','.join(map(str, order))} --kills {t},{t2}\n")
        print(f"p{n}: t {t} t2 {t2} S1 {A['outs']} ({len(A['plans'])} plans) S2 {B['outs']} ({len(B['plans'])} plans)")
        n += 1
        if n >= a.max:
            break
    print(f"{len(pairs)} disjoint pairs; {n} written; finds by t: { {t: len(v) for t, v in sorted(byt.items())} }")


if __name__ == "__main__":
    main()
