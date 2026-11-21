#!/usr/bin/env python3
"""chem_ld2_e17kill.py: builder C2 (seventeenth holder, Nov 17). PLACED KILLS IN THE PLAIN GROUP ENGINE (toward 232).

The group engine (chem_loaders_core.Decomp) reads 160 products on a wired loader: groups 6/6/6/7 at 39/39/39/43, one
column kill on the last group [M]. A kill of low coordinate t for a group needs, per output, its 76 high coefficients
on t written by the earlier products; a product whose low plane lies in t-perp (both low vectors 0 on t: one in four)
gives nothing, so the planner's rule allows a kill only from about 104 products on (3/4 n >= 76 + margin), and the
first three groups never get one. Here the earlier groups AVOID silent planes: a fan-out key whose plane lies in
t-perp for a planned kill coordinate t is not used (the greedy covers the group from the other keys). Then a kill is
possible from about 78 products on, and the leftover pool rank buys row kills (14 equations each, -1/2 product each
on a 6-group with a column kill). Paper plan [P]: 6 / 6 / 7 + kill t0 / 6 + kill t1 + 2 row kills = 39 + 39 + 43 +
37 = 158 products (233; 232 with the C-4 closure of tools/chem_ld2_e17gc4.py).

    python chem_ld2_e17kill.py <table.npz> --k 14 --plan "6,6,7:1:0:0,6:1:2:1" [--avoid 0,1] [--seeds 0-3] [--pad]

--plan: groups in order, each b[:dl:d:shift] (b outputs, dl column kills, d row kills, kill-order shift);
--avoid i,j: the kill-order positions whose coordinates the groups WITHOUT kills avoid; --pad: on the C-4 padded table
(V = vset), and the closure (tau == kappa per group, as chem_ld2_e17gc4) is required too.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_e17gc4 as gc  # noqa: E402


def rng_list(s):
    out = []
    for part in s.split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out


class KillDecomp(gc.TauDecomp):
    """Group engine with a fixed group plan, placed kills and silent-plane avoidance (closure optional)."""

    def __init__(self, an, seed, effort, plan, avoid_pos, V=None, kap=0, tries=1):
        super().__init__(an, seed, effort, V, kap, tries)
        self.plan_list = plan
        self.order = core.kill_order((1 << an.nl) - 1, an)
        self.avoid_t = [self.order[i] for i in avoid_pos]
        self.cur_avoid = []
        self.shift0 = 0
        self.closure = V is not None

    def _tau_bad(self, wmask):
        if not self.closure:
            return 0
        return super()._tau_bad(wmask)

    def apply_kills(self, outs, coords, dl, d, shift=0):
        return super().apply_kills(outs, coords, dl, d, (shift + self.shift0) % 6)

    def _select(self, keys, wb, dim_w):
        if self.cur_avoid:
            keep = {}
            for key, lst in keys.items():
                if all(any((y >> t) & 1 for y, _ in lst) for t in self.cur_avoid):
                    keep[key] = lst
            res = super()._select(keep, wb, dim_w)
            if res is not None:
                return res
        return super()._select(keys, wb, dim_w)

    def run(self):
        an = self.an
        nb = self.X.shape[1]
        allc = (1 << an.nl) - 1
        classes = {}
        for o in range(nb):
            lv = an.out_level[o]
            classes.setdefault(lv if lv <= 2 else 3, []).append(o)
        for lv in (1, 2):
            outs = classes.get(lv, [])
            coords = sum(1 << t for t in range(an.nl) if an.level[t] <= lv)
            size = 10 if lv == 1 else 2
            for i in range(0, len(outs), size):
                self.solve_group(outs[i:i + size], coords, tag=f"level{lv}")
        self.n_low = len(self.prods)             # products of the level-1 / level-2 groups (coded bits)
        gen = list(classes.get(3, []))
        for (b, dl, d, sh) in self.plan_list:
            if not gen:
                break
            self.cur_avoid = self.avoid_t if (dl == 0 and d == 0) else []
            self.shift0 = sh
            self.solve_group(gen[:b], allc, dl, d, tag="generic")
            gen = gen[b:]
        self.cur_avoid = []
        if gen:                                  # leftovers: the engine's own planner
            plan = core.make_planner(an, allc)
            while gen:
                b, dl, d = plan(len(gen), len(self.prods))
                self.solve_group(gen[:b], allc, dl, d, tag="generic")
                gen = gen[b:]
        assert not self.X.any()
        return self.prods


# 02:12Z: the model-optimal placed-kill plan per number of generic outputs (exhaustive search under the measured kill
# rule, avoid 0,1 in the groups without kills) and its model floor in products; 26 outputs: 164 measured on k4 of
# rk7_3 (the 25-output plan + a leftover: 165) [M]
PLANS_BY_N = {
    18: ("6,6,6:1:0:0", 116),
    19: ("6,6,7:1:0:0", 121),
    20: ("6,6,1:1:0:0,7:1:0:1", 128),
    21: ("6,6,1:1:0:0,1:1:0:1,7:1:0:0", 135),
    22: ("6,6,1:1:0:0,1:1:0:1,1:1:0:0,7:1:0:1", 142),
    23: ("6,6,1:1:0:0,1:1:0:1,1:1:0:0,1:1:0:1,7:1:0:0", 149),
    24: ("6,6,6:1:0:0,6:1:2:1", 153),
    25: ("6,6,7:1:0:0,6:1:3:1", 157),   # 07:05Z Nov 17 (nineteenth holder): 3 row kills: 157 = 231 closed on k4 of
                                        # wr7_gclfq, seeds 5, 7, 9 of 0-9 [M]; a failed third kill falls back to 2
    26: ("6,6,7:1:0:0,7:1:0:1", 164),
    27: ("6,6,7:1:0:0,8:1:0:1", 171),
    28: ("6,6,7:1:0:0,7:1:0:1,2:2:0:2", 177),
}


def parse_plan(s):
    out = []
    for g in s.split(","):
        f = [int(x) for x in g.split(":")] + [0, 0, 0]
        out.append((f[0], f[1], f[2], f[3]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--avoid", default="")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--emit", default=None)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, kap = None, 0
    if a.pad:
        import chem_ld2_c4 as c4
        V = c4.vset(words, a.k)
        words = c4.pad_table(words, {a.k: V})
    an = core.analyze(words, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    plan = parse_plan(a.plan)
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    best = None
    for s in rng_list(a.seeds):
        t0 = time.time()
        kd = KillDecomp(an, s, max(2, a.effort), plan, avoid, V, kap, a.tries if a.pad else 1)
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {exc}", flush=True)
            continue
        n = len(prods)
        line = " | ".join(x.split("generic ")[-1] for x in kd.log if "generic" in x)
        print(f"seed {s}: {n} products = {an.infra + n} ({an.infra} infra){' -1 with C-4' if a.pad else ''}; "
              f"{line} ({time.time() - t0:.1f} s)", flush=True)
        if best is None or n < best[0]:
            best = (n, s, prods)
    if best is None:
        return
    print(f"BEST: {best[0]} products, seed {best[1]}", flush=True)
    if a.emit and a.pad:
        import chem_ld2_c4 as c4
        text, st = c4.emit_c4(an, best[2], V)
        os.makedirs(a.emit, exist_ok=True)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx: toffoli {st['toffoli']} (68-AND cache, checked)", flush=True)


if __name__ == "__main__":
    main()
