#!/usr/bin/env python3
"""r63_chem_wired_rp_dump.py: re-pricer of R63-CHEM-WIRED. Runs the DESIGNER'S engine (tools/r63_chem_wired_dc.py:
dead words as free kill unknowns, over R59's TieKill; imported read-only) and writes a GROUP DUMP in the format of
builder C2's group tool (tools/chem_ld2_e19grp.py cmd_dump: outs, coords, a, b, dim_w, wb, prods, made), so that
C2's annealer / stats / lk can be run on the designer's groups, in particular the LAST 7-output group with two
column kills (b 7, a 13, 91 cells: the engine takes the bonus a - 2 = 11 relations = 40 products; the designer's
counted next step asks 13 relations = 39 products, plan 152 = 226).

    CHEM_BROT=7 CHEM_HB=77 python -I tools/r63_chem_wired_rp_dump.py <table.npz> --k 14 --own cur,prev \
        --plan "6,6,6:1:3:1,7:2:0:2" --modes reg,tie3,-,- --avoid 0,1 --seeds 0-0 --pad [--worder auto] \
        [--nofree] --out dump.json

Then (C2's tool, unchanged):  python -I tools/chem_ld2_e19grp.py anneal dump.json --g 3 --P 39 --secs 600 ...
Nothing of the designer's or of the builders' is edited; the probe is a subclass.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_e17kill as ek  # noqa: E402
import r63_chem_wired_dc as dc  # noqa: E402
import r59_chem_groups_tie as tie  # noqa: E402

bits = core.bits


class ProbeDcTie(dc.DcMixin, tie.TieKill):
    """the designer's DcTie with C2's dump probe (records every generic group as solved)."""

    def _probe_init(self):
        self.groups_dump = []
        self._live, self._wb = None, None

    def apply_kills(self, outs, coords, dl, d, shift=0):
        live = super().apply_kills(outs, coords, dl, d, shift)
        if live is not None:
            self._live = live
        return live

    def _select(self, keys, wb, dim_w):
        self._wb = list(wb)
        return super()._select(keys, wb, dim_w)

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        self._live, self._wb = coords, None
        n0 = len(self.prods)
        made = super().solve_group(outs, coords, dl, d, tag)
        if tag == "generic" and self._wb is not None:
            lc = bits(self._live)
            loc = {t: i for i, t in enumerate(lc)}
            oi = {o: i for i, o in enumerate(outs)}
            prods = []
            for p in self.prods[n0:]:
                if any(t not in loc for t in bits(p[0])) or (p[2] and any(t not in loc for t in bits(p[2]))):
                    continue
                y1 = sum(1 << loc[t] for t in bits(p[0]))
                y2 = sum(1 << loc[t] for t in bits(p[2])) if p[2] else 0
                s = sum(1 << oi[o] for o in bits(p[4]) if o in oi)
                prods.append([y1, y2, s])
            self.groups_dump.append(dict(outs=list(outs), coords=lc, a=len(lc), b=len(outs), dim_w=len(self._wb),
                                         wb=[hex(v) for v in self._wb], prods=prods, made=made,
                                         log=self.log[-1] if self.log else ""))
        return made


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--modes", default="reg,tie3,-,-")
    ap.add_argument("--avoid", default="")
    ap.add_argument("--seeds", default="0-0")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--own", default="cur,prev")
    ap.add_argument("--nofree", action="store_true")
    ap.add_argument("--lmin", type=int, default=3)
    ap.add_argument("--worder", default="")
    ap.add_argument("--want", default="", help="accept the first seed whose group counts read e.g. 39,38,36,40")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    V, kap, wp = None, 0, words
    if a.pad:
        import chem_ld2_c4 as c4
        V = c4.vset(words, a.k)
        wp = c4.pad_table(words, {a.k: V})
    an = core.analyze(wp, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    if a.worder:
        dpw = dc.dead_per_word(words, edges, a.k, a.own.split(","), a.lmin)
        wo = sorted(range(core.J), key=lambda j: (dpw[j], j)) if a.worder == "auto" else [int(x) for x in a.worder.split(",")]
        order = sorted(range(len(an.out_wires)), key=lambda q: (wo.index(an.out_wires[q][0]), an.out_wires[q][1]))
        dc.rebase(an, order)
        print(f"k{a.k}: dead addresses a word {dpw}; word order {wo}", flush=True)
    free, nfree, ndead = ({}, 0, 0)
    if not a.nofree:
        free, nfree, ndead = dc.free_bits(words, edges, a.k, an, V, a.own.split(","), a.lmin, None)
    print(f"k{a.k}: generic outputs {an.X0.shape[1]}, infra {an.infra}; dead words owned {ndead}; free bits {nfree}",
          flush=True)
    plan = ek.parse_plan(a.plan)
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    modes = a.modes.split(",")
    for s in ek.rng_list(a.seeds):
        t0 = time.time()
        kd = ProbeDcTie(an, s, max(2, a.effort), plan, avoid, modes, V, kap, a.tries if a.pad else 1, 40)
        kd._dc_init(free)
        kd._probe_init()
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {exc}", flush=True)
            continue
        made = [g["made"] for g in kd.groups_dump]
        line = " | ".join(f"b{g['b']} a{g['a']} dimW{g['dim_w']} -> {g['made']}" for g in kd.groups_dump)
        print(f"seed {s}: {len(prods)} products; {line} ({time.time() - t0:.1f} s)", flush=True)
        if a.want and ",".join(str(m) for m in made) != a.want:
            continue
        out = dict(table=a.table, k=a.k, plan=a.plan, seed=s, pad=a.pad, products=len(prods), infra=an.infra,
                   nb=int(an.X0.shape[1]), nr=int(an.X0.shape[0]), free=not a.nofree, groups=kd.groups_dump)
        with open(a.out, "w") as fh:
            json.dump(out, fh, indent=0)
        print(f"DUMP {a.out}: seed {s}, {len(prods)} products, groups {made}", flush=True)
        return
    print("no seed matched; no dump written", flush=True)


if __name__ == "__main__":
    main()
