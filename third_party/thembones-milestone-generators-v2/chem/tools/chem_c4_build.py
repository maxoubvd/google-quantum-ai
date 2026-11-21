#!/usr/bin/env python3
"""chem_c4_build.py: builder C4. Close a wired loader with chosen generic groups planned from plan files and the rest
by the production engine (designer R59's TieKill over builder C2's placed-kill engine with the C-4 closure; all imported
read-only), then emit with C2's emit_c4 (68-AND cache, every wire checked on every address).

New against C2's `chem_ld2_e19grp.py build`:
  * the injected group's plan advances on EVERY attempt of the closure's retry loop (C2's advances only on kills);
  * --order: the generic outputs in any order (a group = any subset of outputs);
  * --kills: the kill coordinates first in the kill order (group 3's column kill on t1, group 4's on t2);
  * --dumpg: every generic group's W basis (after its kills) and plan, in C2's e19grp dump format (for its annealer).

  CHEM_BROT=7 CHEM_HB=77 python -I chem_c4_build.py <table.npz> --k 2 --order 2,4,5,... --kills 4,13 \
      --inj 0=plansA.json --inj 1=plansB.json --modes=-,-,-,- [--seed 0] [--emit DIR] [--dumpg g.json]
plans files: a list of {planes: [[y1, y2, s], ...]} or of plans [[y1, y2, s], ...] (s over the group's outputs in
the given order, y over the live low coordinates).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_e17kill as kl  # noqa: E402
import r59_chem_groups_tie as rt  # noqa: E402

bits = core.bits


def load_plans(path):
    raw = json.load(open(path))
    out = []
    for pl in raw:
        planes = pl["planes"] if isinstance(pl, dict) else pl
        out.append([tuple(x) for x in planes])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", default="6,6,7:1:0:0,6:1:3:1")
    ap.add_argument("--modes", default="reg,tie3,-,-")
    ap.add_argument("--inj", action="append", default=[], help="g=plans.json (generic group index, 0-based)")
    ap.add_argument("--order", default="", help="generic outputs, comma list (default: index order)")
    ap.add_argument("--kills", default="", help="coordinates put first in the kill order, comma list")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--emit", default=None)
    ap.add_argument("--dumpg", default="")
    ap.add_argument("--dumpall", default="", help="JSON with X0 and every product (y1, A1, y2, A2, s)")
    ap.add_argument("--wirebasis", action="store_true", help="the non-free wires as the output basis (chem_c4_wb)")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    if a.kills:
        first = [int(x) for x in a.kills.split(",")]
        orig = core.kill_order

        def kill_order(coords, an):           # the given coordinates first, then the engine's own order
            base = orig(coords, an)
            return [t for t in first if (coords >> t) & 1] + [t for t in base if t not in first]
        core.kill_order = kill_order           # the engine reads it through the module (e17kill, Decomp)
    import chem_ld2_c4 as c4
    import chem_ld2_e19grp as eg
    if os.environ.get("C4_DROP") == "1":          # holder 3: the kill system's free bits buy one more row dependency
        import chem_c4_drop
        chem_c4_drop.install()
    cubic = None
    if os.environ.get("C4_CUBIC", "") != "":       # holder 4: a rank-additive cubic block in generic group C4_CUBIC
        import chem_c4_cubic as cubic
        cubic.install(int(os.environ["C4_CUBIC"]))
    if os.environ.get("C4_DEBUG"):                 # log why attempts of the engine's group solve fail
        inner = core.Decomp.solve_group
        cnt = {}

        def traced(self, outs, coords, dl=0, d=0, tag=""):
            try:
                r = inner(self, outs, coords, dl, d, tag)
                key = "ok " + (self.log[-1].split("generic ")[-1] if self.log else "")
            except (RuntimeError, AssertionError) as exc:
                key = f"{type(exc).__name__}: {str(exc)[:90]}"
                cnt[key] = cnt.get(key, 0) + 1
                if cnt[key] <= 2:
                    print(f"  attempt on {outs}: {key}", flush=True)
                raise
            cnt[key] = cnt.get(key, 0) + 1
            if cnt[key] <= 2:
                print(f"  attempt on {outs}: {key}", flush=True)
            return r
        core.Decomp.solve_group = traced
    words0 = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, words = eg.block_set(words0, a.k)
    an = core.analyze(words, a.k, 4)
    wires = None
    if a.wirebasis:
        import chem_c4_wb
        wires = chem_c4_wb.wire_basis(an)
    # holder 4, 19:13Z: C4_DEAD_BATCHES = own (default) | prev | both: which batches' dead words pay (prev = batch k - 1,
    # whose other reader k - 1 must then be rebuilt on the toggled words)
    _db = os.environ.get("C4_DEAD_BATCHES", "own")
    dbatches = {"own": None, "prev": [a.k - 1], "both": [a.k, a.k - 1]}[_db]
    deadm = None
    if os.environ.get("C4_DEAD", "") != "":       # holder 4: batch k's dead words as free unknowns (chem_c4_dead)
        import chem_c4_dead as deadm
        edges = np.asarray(np.load(a.table)["tree_edges"], np.int64)
        gen_outs = [o for o in range(an.X0.shape[1]) if an.out_level[o] > 2]
        dcells = deadm.usable_cells(words0, edges, a.k, an, wires, gen_outs, V, dbatches)
        conj = {c["wire"][0] for c in an.carriers if c["conj"]}
        deadm.install(dcells, int(os.environ["C4_DEAD"]), conj)
        print(f"dead cells usable: {len(dcells)}; relations asked {os.environ['C4_DEAD']}", flush=True)
    if os.environ.get("C4_DEAD_ROW", "") != "":   # holder 4: dead-paid left-kernel vectors in a row-bound group (m:gi)
        import chem_c4_dead as deadm
        specs = {}
        for sp in os.environ["C4_DEAD_ROW"].split(","):   # holder 4, 18:42Z: several groups "m:gi:rank,m:gi:rank"
            fr = [int(x) for x in sp.split(":")]
            specs[fr[1]] = (fr[0], fr[2] if len(fr) > 2 else None)
        fr = [int(x) for x in os.environ["C4_DEAD_ROW"].split(",")[0].split(":")]
        mrow, grow = fr[0], fr[1]
        rkrow = fr[2] if len(fr) > 2 else None    # holder 4, 18:33Z: m:gi:rank (the rank the group must reach)
        edges = np.asarray(np.load(a.table)["tree_edges"], np.int64)
        gen_outs = [o for o in range(an.X0.shape[1]) if an.out_level[o] > 2]
        dcells_r = deadm.usable_cells(words0, edges, a.k, an, wires, gen_outs, V, dbatches)
        conj = {c["wire"][0] for c in an.carriers if c["conj"]}
        deadm.install_row(dcells_r, mrow, grow, conj, rkrow, specs if len(specs) > 1 else None)
        print(f"dead cells usable (row hook): {len(dcells_r)}; left-kernel vectors {mrow} in group {grow}", flush=True)
    kap = c4.kappa_wires(an, V)
    injs = {}
    for s in a.inj:
        g, f = s.split("=")
        injs[int(g)] = load_plans(f)
    modes = a.modes.split(",")
    for g in injs:
        if g < len(modes):
            modes[g] = "-"
    stats = {g: dict(att=0, refused=0) for g in injs}
    order = [int(x) for x in a.order.split(",")] if a.order else None
    dumped = []

    class PlanTie(rt.TieKill):
        def run(self):
            if order is None:
                return super().run()
            an_ = self.an
            gen = [o for o in range(self.X.shape[1]) if an_.out_level[o] > 2]
            assert sorted(order) == sorted(gen), "--order must list every generic output once"
            # level-1 / level-2 outputs as the engine does them (wired loaders have none)
            assert all(an_.out_level[o] > 2 for o in range(self.X.shape[1])), "coded outputs: not supported"
            self.n_low = len(self.prods)
            allc = (1 << an_.nl) - 1
            gen = list(order)
            for gi, (b, dl, d, sh) in enumerate(self.plan_list):
                if not gen:
                    break
                self.cur_avoid = self.avoid_t if (dl == 0 and d == 0) else []
                self.shift0 = sh
                self._mode = self.modes[gi] if gi < len(self.modes) and self.modes[gi] != "-" else ""
                self.solve_group(gen[:b], allc, dl, d, tag="generic")
                gen = gen[b:]
            self._mode = ""
            self.cur_avoid = []
            assert not gen, "plan does not cover every output"
            assert not self.X.any()
            return self.prods

        def _tau_bad(self, wmask):
            bad = super()._tau_bad(wmask)
            if os.environ.get("C4_DEBUG"):
                self._taulog = getattr(self, "_taulog", {})
                k = (getattr(self, "_gi", -1), wmask, bad)
                self._taulog[k] = self._taulog.get(k, 0) + 1
                if self._taulog[k] == 1:
                    print(f"  tau group {k[0]}: wmask {wmask:#x} mismatch {bad:#x} ({self.log[-1].split('generic ')[-1] if self.log else ''})",
                          flush=True)
            return bad

        def apply_kills(self, outs, coords, dl, d, shift=0):
            live = super().apply_kills(outs, coords, dl, d, shift)
            if live is not None:
                self._live = live
            return live

        def solve_group(self, outs, coords, dl=0, d=0, tag=""):
            self._cur, self._live, self._wb = None, coords, None
            if tag == "generic":
                self._gi = getattr(self, "_gi", -1) + 1
                if self._gi in injs:
                    self._cur = self._gi
            n0 = len(self.prods)
            made = super().solve_group(outs, coords, dl, d, tag)
            if tag == "generic":
                self.final_log = getattr(self, "final_log", []) + [self.log[-1].split("generic ")[-1]]
                if self._wb is not None:
                    lc = bits(self._live)
                    loc = {t: i for i, t in enumerate(lc)}
                    oi = {o: i for i, o in enumerate(outs)}
                    pr = []
                    for p in self.prods[n0:]:
                        y1 = sum(1 << loc[t] for t in bits(p[0]))
                        y2 = sum(1 << loc[t] for t in bits(p[2])) if p[2] else 0
                        pr.append([y1, y2, sum(1 << oi[o] for o in bits(p[4]) if o in oi)])
                    dumped.append(dict(outs=list(outs), coords=lc, a=len(lc), b=len(outs), dim_w=len(self._wb),
                                       wb=[hex(v) for v in self._wb], prods=pr, made=made, log=self.log[-1]))
            self._cur = None
            return made

        def _select(self, keys, wb, dim_w):
            self._wb = list(wb)
            if self._cur is None:
                return super()._select(keys, wb, dim_w)
            st = stats[self._cur]
            lst = injs[self._cur]
            plan = lst[st["att"] % len(lst)]
            st["att"] += 1
            lc = bits(self._live)
            A = len(lc)
            ylocal = []
            for ym in self.an.ymasks:
                m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
                if m and m not in ylocal:
                    ylocal.append(m)
            chosen, F = [], eg.Span()
            for (y1, y2, s) in plan:
                ps = [p for p, ym in enumerate(ylocal) if not ((y1 | y2) & ~ym)]
                if not ps:
                    st["refused"] += 1
                    return None
                key = (s, ps[0])
                for y in (y1, y2):
                    flat = sum(y << (gi * A) for gi in bits(s))
                    if not F.add(flat):
                        st["refused"] += 1
                        return None
                    chosen.append((key, y, flat))
            if any(F.reduce(w) for w in wb):
                st["refused"] += 1
                return None
            st["plan"] = (st["att"] - 1) % len(lst)
            return len(plan), chosen

    kd = PlanTie(an, a.seed, 3, kl.parse_plan(a.plan), [0, 1], modes, V, kap, a.tries, 40)
    try:
        prods = kd.run()
    except (RuntimeError, AssertionError) as exc:
        print(f"FAILED {type(exc).__name__}: {str(exc)[:150]}; injected {stats}; groups so far "
              f"{getattr(kd, 'final_log', [])}", flush=True)
        if os.environ.get("C4_DROP") == "1":
            print(f"drop stats {chem_c4_drop.STATS}", flush=True)
        if a.dumpg:
            json.dump(dict(table=a.table, k=a.k, groups=dumped), open(a.dumpg, "w"))
        sys.exit(1)
    if deadm is not None:
        toggles = getattr(kd, "dead_toggles", [])
        print(f"dead toggles chosen: {len(toggles)}; hook stats {deadm.HOOK['stats']} row {deadm.HOOK_ROW['stats']}", flush=True)
        if toggles:
            w_new = deadm.write_toggles(words0, a.k, an, wires, toggles, deadm.HOOK["conj"])
            words2 = c4.pad_table(w_new, {a.k: V})   # the OLD block set V stays valid (toggles avoid its rows l <= 2)
            an2 = core.analyze(words2, a.k, 4)
            import chem_c4_wb
            wires2 = chem_c4_wb.wire_basis(an2)
            assert wires2 == wires and an2.G == an.G, "the wire map moved"
            an = an2
            dout = os.environ.get("C4_DEAD_OUT", "")
            if dout == "EMIT" and a.emit:                 # one table per build, next to its circuit
                os.makedirs(a.emit, exist_ok=True)
                dout = os.path.join(a.emit, "table_dead.npz")
            if dout and dout != "EMIT":
                z = dict(np.load(a.table))
                z["angle_words"] = w_new
                np.savez(dout, **z)
                print(f"table with the dead toggles -> {dout}", flush=True)
    blocks = getattr(kd, "cubic_blocks", []) if cubic is not None else []
    if blocks:
        text, st = cubic.emit_c4x(an, prods, V, blocks)
        print(f"cubic blocks {len(blocks)}: " + "; ".join(f"rank {b['rank0']} -> {b['rank1']} (attempt {b['attempt']})"
                                                         for b in blocks) + f"; stats {cubic.STATS}", flush=True)
    else:
        text, st = c4.emit_c4(an, prods, V)
        if cubic is not None:
            print(f"cubic: no block in the closed build; stats {cubic.STATS}", flush=True)
    line = " | ".join(getattr(kd, "final_log", []))
    print(f"CLOSED {len(prods) + 2 * len(blocks)} products, emit_c4 toffoli {st['toffoli']} (checked); injected {stats}; {line}",
          flush=True)
    if os.environ.get("C4_DROP") == "1":
        print(f"drop stats {chem_c4_drop.STATS}", flush=True)
    if a.dumpg:
        json.dump(dict(table=a.table, k=a.k, groups=dumped), open(a.dumpg, "w"))
    if a.dumpall:                                  # the whole decomposition: X0 = sum of products (holder 2)
        json.dump(dict(table=a.table, k=a.k, nl=int(an.nl), nr=int(an.nr), order=order, plan=a.plan,
                       kills=a.kills, X0=[[int(v) for v in row] for row in np.asarray(an.X0)],
                       ymasks=[int(m) for m in an.ymasks], groups=dumped,
                       prods=[[int(p[0]), hex(int(p[1])), int(p[2]), hex(int(p[3])), int(p[4])] for p in prods]),
                  open(a.dumpall, "w"))
        print(f"dumped all {len(prods)} products -> {a.dumpall}", flush=True)
    if a.emit:
        os.makedirs(a.emit, exist_ok=True)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx", flush=True)


if __name__ == "__main__":
    main()
