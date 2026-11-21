#!/usr/bin/env python3
"""invsqrt_fc_s19stage.py (FC nineteenth holder, Nov 19): an F1 program (I40's decode / I40's mink output: {"gates": [{"left",
"right"}], "outputs": [[...]], "inputs": [...]}, indices = inputs 0..n-1 then gates n..n+K-1) written in GHOST's STAGED shape
(one stage: its first I gates on temporaries, the rest CCX straight into the target), the form that
tools/sh/invsqrt_fc_s18q18m.sh integrates on fc18_M. Temporaries = the gates another gate reads; terminal gates = the rest
(each must be in the output; a terminal gate outside it is dead and dropped). Checked by replay on EVERY row of the problem.
  invsqrt_fc_s19stage.py F1PROG.json PROBLEM.json OUT.json [maxtmp=3] [name=q18]
Safe to run: reads two files, writes OUT only."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def rows_of(prob):
    """the problem's rows as (input bits, target) from its truth tables (I40 / F1 problem format)."""
    n = len(prob["inputs"])
    R = prob["rows"]
    tts = [int(x["tt"], 16) for x in prob["inputs"]]
    tg = int(prob["target"], 16)
    return [([(tts[j] >> r) & 1 for j in range(n)], (tg >> r) & 1) for r in range(R)]


def main():
    F = json.load(open(sys.argv[1]))
    prob = json.load(open(sys.argv[2]))
    outp = sys.argv[3]
    opt = dict(a.split("=", 1) for a in sys.argv[4:] if "=" in a)
    maxtmp = int(opt.get("maxtmp", 3))
    n = len(F["inputs"])
    assert F["inputs"] == [x["name"] for x in prob["inputs"]], (F["inputs"], [x["name"] for x in prob["inputs"]])
    G = F["gates"]
    K = len(G)
    outs = set(F["outputs"][0])
    read = set()
    for g in G:
        read.update(i for i in g["left"] + g["right"] if i >= n)
    internal = [j for j in range(K) if n + j in read]
    terminal = [j for j in range(K) if n + j not in read and n + j in outs]
    dead = [j for j in range(K) if n + j not in read and n + j not in outs]
    print("F1 program: %d gates; temporaries %s, terminal %s, dead %s" % (K, internal, terminal, dead))
    if len(internal) > maxtmp:
        print("NOT STAGEABLE in %d temporaries: %d gates are read by others" % (maxtmp, len(internal)))
        sys.exit(2)
    order = internal + terminal
    new = {n + j: n + p for p, j in enumerate(order)}
    rm = lambda lst: sorted(new.get(i, i) for i in lst)  # noqa: E731
    gates = [{"left": rm(G[j]["left"]), "right": rm(G[j]["right"])} for j in order]
    out = sorted(i for i in outs if i < n) + sorted(new[i] for i in outs if i >= n and (i - n) in internal)
    P = {"inputs": F["inputs"], "internal": 0, "stages": [{"gates": gates, "internal": len(internal), "out": out}],
         "outputs": [[]], "target_name": opt.get("name", "q18"),
         "note": "FC19: %d-AND program staged from an F1 program (%s)" % (len(order), os.path.basename(sys.argv[1]))}
    bad = 0
    for bits, y in rows_of(prob):                      # replay: gates in order; target = out terms ^ terminal gates
        v = list(bits)
        for g in gates:
            a = 0
            for i in g["left"]:
                a ^= v[i]
            b = 0
            for i in g["right"]:
                b ^= v[i]
            v.append(a & b)
        t = 0
        for i in out:
            t ^= v[i]
        for p in range(len(internal), len(order)):
            t ^= v[n + p]
        bad += (t != y)
    print("staged replay: wrong on %d of %d rows" % (bad, prob["rows"]))
    if bad:
        sys.exit(1)
    json.dump(P, open(outp, "w"), indent=1)
    print("wrote %s: %d ANDs, %d temporaries" % (outp, len(order), len(internal)))


if __name__ == "__main__":
    main()
