"""invsqrt_gen: take a table set made by another workstream of our team (e.g. a tie set of the table-freedom
workstream, research/ws_invsqrt_ties.md) into OUR folder and check it (own tool).

usage: invsqrt_gen_tabs.py import <source.json> <name> <Params key=value ...>
         copies the lists H, C0, C1, C2 of the source file to ~/qcc/ws/invsqrt_gen/tabs/<name>.json (with the path
         and sha256 of the source), then checks them as `check` does.
       invsqrt_gen_tabs.py check <name> <Params key=value ...>
         the ties the builder relies on (C0 even, C1 bit 0 = C2 bit 7, window safety, parity of C0 & tie_mask =
         bit 0 of T, range of the hi word) and the nominal pipeline on ALL inputs against the verifier's
         reference (worst error / tolerance, margin, T per rotation), and how the tables differ from our fit.
The source is only read. A build takes the tables with the option tabs=<name> (invsqrt_gen_build.py).
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_model as M


def check(name, P):
    tabs = M.load_tables(name)
    C0, C1, C2, H = (np.array(tabs[k], dtype=np.int64) for k in ("C0", "C1", "C2", "H"))
    T = np.array([P.top_table(int(v))[1] for v in C1])
    par = np.zeros(256, dtype=np.int64)
    for bit in range(1, P.wlo):
        if (P.tie_mask >> bit) & 1:
            par ^= (C0 >> bit) & 1
    qmax = sum(M.rr(C2, P.a_shift(i)) for i in P.a_rows)
    ties = dict(c0_even=int((C0 & 1).sum()) == 0,
                tie_c1=int(((C1 ^ (C2 >> P.a_shift(P.a_lo))) & 1).sum()) == 0,
                window=bool((qmax <= (1 << (P.win - 1))).all() and (C1 >= (1 << P.win)).all()),
                tie_t0=int((par ^ (T & 1)).sum()) == 0,
                hi=(int(H.min()), int(H.max())), c2=(int(C2.min()), int(C2.max())),
                c1=(int(C1.min()), int(C1.max())))
    if P.zone:
        # zone format: an address that does not occur must hold a copy of the even address below it (the
        # builder keeps every function of the address in that form), and u must stay in range on every input
        src = [(a & ~1) if M.zone1(a >> 7, a & 127) else a for a in range(256)]
        ties["zone_copies"] = all(tabs[k][a] == tabs[k][src[a]] for k in ("H", "C0", "C1", "C2") for a in range(256))
        ties["used"] = len(set(src))
    n_rot = 40 if P.neg else 41
    out, ev, _ = M.check(P, tabs, n_rot=n_rot)
    out.pop("by_level")
    if P.zone:
        _, _, _, _, a8p, tp = M.decompose(ev["x"])
        a8z, tz = M.zone_apply(a8p, tp)
        uu, _ = M.rows(P, np.asarray(tabs["C1"], dtype=np.int64)[a8z], np.asarray(tabs["C2"], dtype=np.int64)[a8z], tz)
        ties["u_in_range"] = bool((uu >= 0).all() and (uu < (1 << P.wu)).all())
        ties["zone1_inputs"] = int((a8z != a8p).sum() + ((a8p & 1) == 0)[M.zone1(a8p >> 7, a8p & 127) == 1].sum())
    print(f"tables {name} sha {tabs['sha']}  ties {ties}")
    print(f"  all inputs: {'PASS' if out['ok'] else 'FAIL'}  failing {out['n_fail']}  worst {out['worst']:.5f} at x = "
          f"{out['worst_x']}  margin {out['min_margin']} ({100 * out['margin_frac']:.2f}% of the smallest tolerance)"
          f"  T per rotation at {n_rot}: {out.get('t_per_rot', float('nan')):.3f}")
    own = os.path.join(M.WS, "scratch", f"tables_{P.key()}.json")
    if os.path.exists(own):
        o = json.load(open(own))
        print("  segments that differ from our fit: " + "  ".join(
            f"{k}: {sum(1 for a, b in zip(tabs[k], o[k]) if a != b)}" for k in ("H", "C0", "C1", "C2")))


def main():
    cmd = sys.argv[1]
    if cmd == "import":
        src, name, rest = sys.argv[2], sys.argv[3], sys.argv[4:]
        raw = open(os.path.expanduser(src), "rb").read()
        d = json.loads(raw)
        out = {k: [int(v) for v in d[k]] for k in ("H", "C0", "C1", "C2")}
        out["source"] = src
        out["source_sha256"] = hashlib.sha256(raw).hexdigest()
        os.makedirs(os.path.join(M.WS, "tabs"), exist_ok=True)
        with open(os.path.join(M.WS, "tabs", f"{name}.json"), "w") as fh:
            json.dump(out, fh)
        print(f"imported {src} (sha256 {out['source_sha256'][:16]}...) as tabs/{name}.json; other fields of the "
              f"source: {sorted(k for k in d if k not in out)}")
    else:
        name, rest = sys.argv[2], sys.argv[3:]
    P = M.Params(**{k: int(v) for k, v in (a.split("=") for a in rest)})
    check(name, P)


if __name__ == "__main__":
    main()
