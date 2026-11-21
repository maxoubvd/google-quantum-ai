#!/usr/bin/env python3
"""chem_c4_dump.py: builder C4. Run the production wired-loader engine (designer R59's TieKill over builder C2's
placed-kill engine; both imported read-only) on one loader of a table and save the FULL decomposition:
every product [y1, A1, y2, A2, s, phase] (low masks over the 15 low coordinates, high vectors over the nr rows as
ints, fan-out over the output basis), the group boundaries, and the mixed table X0 (nr x outputs, uint32 low masks).

  CHEM_BROT=7 CHEM_HB=77 python -I chem_c4_dump.py <submission.npz> --k 2 --seeds 0-15 --out d.json
The plan and modes are the production ones for a 25-output wired loader (R59 'For C2' item 4).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_e17kill as kl  # noqa: E402
import r59_chem_groups_tie as rt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", default="6,6,7:1:0:0,6:1:3:1")
    ap.add_argument("--modes", default="reg,tie3,-,-")
    ap.add_argument("--seeds", default="0-15")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--stop", type=int, default=156)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    import chem_ld2_c4 as c4
    import chem_ld2_e19grp as eg
    words0 = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = eg.block_set(words0, a.k)
    an = core.analyze(words, a.k, 4)
    kap = c4.kappa_wires(an, V)
    plan = kl.parse_plan(a.plan)
    modes = a.modes.split(",")
    best = None

    class Probe(rt.TieKill):
        def solve_group(self, outs, coords, dl=0, d=0, tag=""):
            n0 = len(self.prods)
            made = super().solve_group(outs, coords, dl, d, tag)
            self.bounds.append(dict(tag=tag, outs=list(outs), start=n0, end=len(self.prods), made=made,
                                    log=self.log[-1] if self.log else ""))
            return made

    for s in kl.rng_list(a.seeds):
        t0 = time.time()
        kd = Probe(an, s, 3, plan, [0, 1], modes, V, kap, a.tries, 40)
        kd.bounds = []
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {type(exc).__name__} ({time.time() - t0:.1f} s)", flush=True)
            continue
        n = len(prods)
        print(f"seed {s}: {n} products = {an.infra + n} ({time.time() - t0:.1f} s); "
              + " | ".join(f"{b['tag']} {b['made']}" for b in kd.bounds if b['tag'] == 'generic'), flush=True)
        if best is None or n < best[0]:
            best = (n, s, [list(map(int, p)) for p in prods], kd.bounds)
        if a.stop and n <= a.stop:
            break
    if best is None:
        sys.exit("no seed closed")
    n, s, prods, bounds = best
    out = dict(table=a.table, k=a.k, seed=s, products=n, infra=an.infra, nr=int(an.nr), nl=int(an.nl),
               nb=int(an.X0.shape[1]), level=[int(x) for x in an.level], out_level=[int(x) for x in an.out_level],
               X0=[[int(v) for v in row] for row in an.X0], prods=prods, bounds=bounds,
               V=int(V) if V is not None else None)
    json.dump(out, open(a.out, "w"))
    print(f"BEST {n} products (seed {s}) -> {a.out}", flush=True)


if __name__ == "__main__":
    main()
