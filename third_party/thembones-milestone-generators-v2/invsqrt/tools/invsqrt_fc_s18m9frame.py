#!/usr/bin/env python3
"""invsqrt_fc_s18m9frame.py -- FC eighteenth holder: digit 23's X-measurement (bit m9 of row 24's undo) read off the
EMITTED circuit: the kit's exact simulator (invsqrt_gen_circ.Sim, read-only) runs the .kmx on ALL inputs; at every HMR
whose wire holds digit 23 (R24's w9, or its complement) on every input, the other live non-x wires are named against
the arithmetic bits (v0..v6, sg, w0..w9, framed by 1 and by sg) and the spec of the measured wire as a function of
them is written in the CIRCUIT'S OWN FRAME (input of tools/invsqrt_fc_s18cubrank.py).
  invsqrt_fc_s18m9frame.py KMX SRC HTAB OUTPREFIX         (writes OUTPREFIX_<op>.json)
Safe to run: reads only; writes the spec files."""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    kmx, src, hpath, outp = [os.path.expanduser(a) for a in sys.argv[1:5]]
    import r65_inv_chainall_carry as CT
    d, U = CT.chain(os.path.abspath(src), hpath)
    import invsqrt_gen_circ as C
    import invsqrt_lib as L
    v = (d - 257) >> 1
    R = (np.int64(1) << 35) + (U << 8)
    for b in (26, 25, 24):
        neg = R < 0
        Rn = R - np.where(neg, -1, 1) * (d << b)
        if b == 24:
            sg = neg.astype(np.uint8)
            Wa = (Rn >> 24) & 1023
        R = Rn
    ar = {"v%d" % j: ((v >> j) & 1).astype(np.uint8) for j in range(7)}
    ar["sg"] = sg
    for p in range(10):
        ar["w%d" % p] = ((Wa >> p) & 1).astype(np.uint8)
    n = len(d)
    xs = np.arange(1, n + 1, dtype=np.uint64)
    pk = {k: C.pack_bits(a.astype(np.uint64), 0) for k, a in ar.items()}
    header, ops = L.parse_kmx(kmx)
    ops = [op for op in ops if op[0] != "INFO_TAG"]
    nq = 1 + max(q for op in ops for q in op[1])
    sim = C.Sim(xs, nq, hmr="zero")
    full = sim.full
    w9 = pk["w9"]
    nhit = 0
    at = set(int(a) for a in os.environ.get("FC18_AT", "").split(",") if a)    # explicit op indices (HMRs)
    for i, op in enumerate(ops):
        if op[0] == "HMR" and (i in at or (not at and sim.q[op[1][0]] in (w9, w9 ^ full))):
            t = op[1][0]
            nhit += 1
            print("op %d: HMR on q%d holds %s" % (i, t, "digit 23" if sim.q[t] == w9 else "digit 23 ^ 1" if
                                                    sim.q[t] == w9 ^ full else "(not digit 23)"), flush=True)
            live = [w for w in range(nq) if w != t and sim.q[w]]
            names, keep = [], []
            for w in live:
                lab = []
                for k, p in pk.items():
                    for fr, s in (("", 0), ("^1", full), ("^sg", pk["sg"]), ("^sg^1", pk["sg"] ^ full)):
                        if sim.q[w] == p ^ s:
                            lab.append(k + fr)
                if lab:
                    names.append("q%d" % w)
                    keep.append(w)
                    print("   q%d = %s" % (w, " = ".join(lab)))
            print("   %d live wires, %d of them arithmetic bits (named above)" % (len(live), len(keep)))
            allc = {w: C.unpack(sim.q[w], n).astype(np.int64) for w in live}
            tg = C.unpack(sim.q[t], n).astype(np.int64)

            def isfun(ws):
                if len(ws) > 62:
                    return True
                key = np.zeros(n, dtype=np.int64)
                for j, w in enumerate(ws):
                    key |= allc[w] << j
                o = np.argsort(key, kind="stable")
                k2, t2 = key[o], tg[o]
                s2 = np.flatnonzero(np.r_[True, k2[1:] != k2[:-1]])
                return bool(np.array_equal(np.maximum.reduceat(t2, s2), np.minimum.reduceat(t2, s2)))
            sup = list(live)                       # greedy: drop wires (x register first) while it stays a function
            for w in sorted(live, key=lambda w: (w >= 22, w)):
                if w in keep:
                    continue
                trial = [u for u in sup if u != w]
                if isfun(trial):
                    sup = trial
            print("   a minimal support (greedy, x register dropped first): %s" % ["q%d" % w for w in sup], flush=True)
            if os.environ.get("FC18_WIDE"):        # the support plus every live non-x wire (AND results can lower the degree)
                sup = sorted(set(sup) | set(w for w in live if w >= 22))
                print("   WIDE support: %s" % ["q%d" % w for w in sup], flush=True)
            if os.environ.get("FC18_HELP"):        # the support plus every live wire that is itself a function of it
                helpers = []
                for w in live:
                    if w in sup:
                        continue
                    key = np.zeros(n, dtype=np.int64)
                    for j, u in enumerate(sup):
                        key |= allc[u] << j
                    o = np.argsort(key, kind="stable")
                    k2, t2 = key[o], allc[w][o]
                    s2 = np.flatnonzero(np.r_[True, k2[1:] != k2[:-1]])
                    if np.array_equal(np.maximum.reduceat(t2, s2), np.minimum.reduceat(t2, s2)):
                        helpers.append(w)
                print("   HELPERS (live wires that are functions of the support): %s" % ["q%d" % w for w in helpers],
                      flush=True)
                sup = sorted(set(sup) | set(helpers))
            names = ["q%d" % w for w in sup]
            keep = sup
            live = keep
            cols = [C.unpack(sim.q[w], n).astype(np.int64) for w in live]
            tgt = C.unpack(sim.q[t], n).astype(np.int64)
            key = np.zeros(n, dtype=np.int64)
            for j, c in enumerate(cols):
                key |= c << j
            order = np.argsort(key, kind="stable")
            ks, ts = key[order], tgt[order]
            st = np.flatnonzero(np.r_[True, ks[1:] != ks[:-1]])
            mx, mn = np.maximum.reduceat(ts, st), np.minimum.reduceat(ts, st)
            func = bool(np.array_equal(mx, mn))
            print("   %d named wires, %d patterns; the measured wire IS %sa function of them" % (
                len(live), len(st), "" if func else "NOT "), flush=True)
            if func:
                out = "%s_%d.json" % (outp, i)
                json.dump(dict(inputs=names, rows=list(zip(ks[st].tolist(), mn.tolist())), target="q%d" % t, op=i),
                          open(out, "w"))
                print("   spec", out, flush=True)
        sim.run(ops, i, i + 1)
    print("%d HMRs examined" % nhit)


if __name__ == "__main__":
    main()
