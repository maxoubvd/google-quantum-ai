#!/usr/bin/env python3
"""r65_chem_wired221_build.py: designer R65-CHEM-WIRED221. The owner engine with LIVE UNITS AT EXACT PRICE
(tools/r65_chem_wired221_dc.py) + the annealed last group (builder C2's group tool chem_ld2_e19grp.py; the probe and
the injection are copies of R63's re-pricer's tools r63_chem_wired_rp_dump.py / _rp_build.py over my free bits).

    1. dump   python -I tools/r65_chem_wired221_build.py T --k 14 --plan P --worder W --avoid 0 --seeds 0-7 --pad \
                  --prices prices_k14.npz --tau 1e-3 --maxbit 2 [--want 39,36,36,40] --dump d.json
    2. anneal python -I tools/chem_ld2_e19grp.py anneal d.json --g 3 --P 39 --price 37 --secs 300 --multi --rseed R --save an_R.json
    3. plans  python -I tools/chem_ld2_e19grp.py plan d.json spaces.json --g 3 --tries 4000 --out plans.json
    4. build  python -I tools/r65_chem_wired221_build.py T --k 14 --plan P --worder W --avoid 0 --seeds S --pad \
                  --prices prices_k14.npz --tau 1e-3 --maxbit 2 --plans plans.json --g 3 --a 13 --b 7 --emit DIR
       (same table, seed and options as the dump). Writes DIR/angle_k_K.kmx and DIR/table_kK.npz (the new words).
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("CHEM_BROT", "7")
os.environ.setdefault("CHEM_HB", "77")
import chem_loaders_core as core  # noqa: E402
import chem_ld2_wwidth  # noqa: E402,F401
import chem_ld2_e17kill as ek  # noqa: E402
import r63_chem_wired_dc as r63  # noqa: E402
import r59_chem_groups_tie as tie  # noqa: E402
import chem_ld2_e19grp as eg  # noqa: E402
import r65_chem_wired221_dc as my  # noqa: E402
from r62_chem_cache_dead import dead_mask  # noqa: E402

bits = core.bits


PSHIFT = {}       # SECOND HOLDER: generic group -> index of the first plan offered (a scan over plans = runs with shifts)


def make_cls(injs, pa, pb):
    class Kd(r63.DcMixin, tie.TieKill):
        """probe (dump of every generic group) + injection of annealed planes into chosen groups."""

        def _x_init(self):
            self.groups_dump = []
            self._live, self._wb, self._cur, self._inject, self._pidx = None, None, None, None, {}
            self._gi = -1

        def apply_kills(self, outs, coords, dl, d, shift=0):
            live = super().apply_kills(outs, coords, dl, d, shift)
            if live is not None:
                self._live = live
            return live

        def solve_group(self, outs, coords, dl=0, d=0, tag=""):
            self._live, self._wb, self._cur, self._inject = coords, None, None, None
            n0 = len(self.prods)
            if tag == "generic":
                self._gi += 1
                if self._gi in injs:
                    self._cur = self._gi
                    self._inject = (outs, coords)
                    self._planes = injs[self._gi][0]
            made = super().solve_group(outs, coords, dl, d, tag)
            self._cur = None
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

        def _select(self, keys, wb, dim_w):
            self._wb = list(wb)
            if self._cur is None:
                return super()._select(keys, wb, dim_w)
            outs, coords = self._inject
            lc = bits(self._live)
            A = len(lc)
            lst = injs[self._cur]
            pi = self._pidx.get(self._cur, PSHIFT.get(self._cur, 0) - 1) + 1     # SECOND HOLDER: next plan a call
            self._pidx[self._cur] = pi
            self._planes = lst[pi % len(lst)]
            ylocal = []
            for ym in self.an.ymasks:
                m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
                if m and m not in ylocal:
                    ylocal.append(m)
            chosen, F, why = [], eg.Span(), None
            for (y1, y2, s) in self._planes:
                ps = [p for p, ym in enumerate(ylocal) if not ((y1 | y2) & ~ym)]
                if not ps:
                    why = "no phase"
                    break
                key = (s, ps[0])
                for y in (y1, y2):
                    flat = sum(y << (gi * A) for gi in bits(s))
                    if not F.add(flat):
                        why = "dependent features"
                    chosen.append((key, y, flat))
            if why is None and any(F.reduce(w) for w in wb):
                why = "W not inside the plan's span"
            if why is not None:
                self._refused = getattr(self, "_refused", {})
                self._refused[why] = self._refused.get(why, 0) + 1
                return None
            self._accepted = getattr(self, "_accepted", 0) + 1
            return len(self._planes), chosen

    return Kd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--modes", default="-,-,-,-")
    ap.add_argument("--avoid", default="")
    ap.add_argument("--seeds", default="0-0")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--own", default="cur,prev")
    ap.add_argument("--lmin", type=int, default=3)
    ap.add_argument("--worder", default="")
    ap.add_argument("--wpat", default="", help="word order as a pattern over the words ranked by dead addresses (0 = poorest), e.g. 0,3,1,2,4")
    ap.add_argument("--prices", default=None)
    ap.add_argument("--oauto", default="")
    ap.add_argument("--oorder", default="")
    ap.add_argument("--tau", type=float, default=0.0)
    ap.add_argument("--maxbit", type=int, default=1)
    ap.add_argument("--want", default="")
    ap.add_argument("--maxprice", type=float, default=1e9, help="dump / build only a seed whose formula price is at most this")
    ap.add_argument("--dump", default=None)
    ap.add_argument("--plans", default=None)
    ap.add_argument("--g", type=int, default=3)
    ap.add_argument("--a", type=int, default=13)
    ap.add_argument("--b", type=int, default=7)
    ap.add_argument("--maxplans", type=int, default=5000)
    ap.add_argument("--inj", default="", help="SECOND HOLDER: several groups: G=plans.json;G=plans.json (replaces --plans / --g)")
    ap.add_argument("--pshift", default="", help="SECOND HOLDER: G=N;G=N: the first plan offered to group G")
    ap.add_argument("--emit", default=None)
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    z = np.load(a.table)
    words = np.asarray(z["angle_words"], np.int64)
    edges = np.asarray(z["tree_edges"], np.int64)
    import chem_ld2_c4 as c4
    V, kap, wp = None, 0, words
    if a.pad:
        V = c4.vset(words, a.k)
        wp = c4.pad_table(words, {a.k: V})
    an = core.analyze(wp, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    own = a.own.split(",")
    if a.worder or a.wpat:
        dpw = r63.dead_per_word(words, edges, a.k, own, a.lmin)
        if a.wpat:
            rk = sorted(range(core.J), key=lambda j: (dpw[j], j))
            wo = [rk[int(x)] for x in a.wpat.split(",")]
        else:
            wo = sorted(range(core.J), key=lambda j: (dpw[j], j)) if a.worder == "auto" else [int(x) for x in a.worder.split(",")]
        print(f"k{a.k}: dead addresses a word {dpw}; word order {wo}", flush=True)
        order = sorted(range(len(an.out_wires)), key=lambda q: (wo.index(an.out_wires[q][0]), an.out_wires[q][1]))
        r63.rebase(an, order)
    live = None
    if a.tau > 0:
        live = my.load_prices(a.prices, a.tau, a.maxbit)
    if a.oauto or a.oorder:
        sup = my.supply(an, words, edges, a.k, own, a.lmin, live)
        sizes = [g[0] for g in ek.parse_plan(a.plan)]
        order, seq = my.placed_order(an, sup, sizes, a.oauto, a.oorder)
        r63.rebase(an, order)
        print(f"k{a.k}: outputs placed: " + ",".join(f"{w[0]}:{w[1]}" for w in seq), flush=True)
    free, nfree, nlive, ndead = my.free_bits(words, edges, a.k, an, V, own, a.lmin, live)
    injs = {}
    if a.plans:
        plans = json.load(open(a.plans))
        injs = {a.g: [[tuple(x) for x in pl["planes"]] for pl in plans[:a.maxplans]]}
    for item in [x for x in a.inj.split(";") if x.strip()]:
        gi, f = item.split("=", 1)
        injs[int(gi)] = [[tuple(x) for x in pl["planes"]] for pl in json.load(open(f))[:a.maxplans]]
    for item in [x for x in a.pshift.split(";") if x.strip()]:
        gi, n = item.split("=", 1)
        PSHIFT[int(gi)] = int(n)
    print(f"k{a.k}: generic outputs {an.X0.shape[1]}, infra {an.infra}; dead bits {nfree}; live units {nlive}; "
          f"plans offered { {g: len(v) for g, v in injs.items()} }", flush=True)
    Cls = make_cls(injs, a.a, a.b)
    plan = ek.parse_plan(a.plan)
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    ntry = max([a.tries] + [len(v) for v in injs.values()])
    done = None
    for s in ek.rng_list(a.seeds):
        t0 = time.time()
        kd = Cls(an, s, max(2, a.effort), plan, avoid, a.modes.split(","), V, kap, ntry if a.pad else 1, 40)
        kd._dc_init(free)
        kd._x_init()
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError, ValueError) as exc:
            print(f"seed {s}: FAILED {type(exc).__name__}: {str(exc)[:160]}; accepted {getattr(kd, '_accepted', 0)}, "
                  f"refused {getattr(kd, '_refused', {})} ({time.time() - t0:.1f} s)", flush=True)
            continue
        n = len(prods)
        used = {o: list(v) for o, v in kd.used.items()}
        lcost = sum(f["cost"] for lst in used.values() for f in lst)
        nl = sum(1 for lst in used.values() for f in lst if f["cost"] != 0.0)
        made = [g["made"] for g in kd.groups_dump]
        line = " | ".join(f"{b}->{m}" + (f" [kills {st[0]},{st[1]} pool {st[2]} free {st[3]}]" if st else "")
                          for b, m, st in kd.gcount)
        print(f"seed {s}: {n} products = {an.infra + n} (-1 with C-4); {line}; live {nl}, formula price {lcost:+.5f} "
              f"mHa^2; accepted {getattr(kd, '_accepted', 0)}, refused {getattr(kd, '_refused', {})} "
              f"({time.time() - t0:.1f} s)", flush=True)
        if a.want and ",".join(str(m) for m in made) != a.want:
            continue
        if lcost > a.maxprice:
            continue
        done = (s, kd, prods, used, lcost, made)
        break
    if done is None:
        print("no seed matched", flush=True)
        return
    s, kd, prods, used, lcost, made = done
    if a.dump:
        out = dict(table=a.table, k=a.k, plan=a.plan, seed=s, pad=a.pad, products=len(prods), infra=an.infra,
                   nb=int(an.X0.shape[1]), nr=int(an.X0.shape[0]), free=True, groups=kd.groups_dump)
        with open(a.dump, "w") as fh:
            json.dump(out, fh, indent=0)
        print(f"DUMP {a.dump}: seed {s}, {len(prods)} products, groups {made}, formula price {lcost:+.5f}", flush=True)
    X0 = an.X0.copy()
    for o, lst in used.items():
        for f in lst:
            for o2 in bits(f["g"]):
                X0[:, o2] ^= f["vec"]
    w2 = r63.modified_words(words, used)
    dead = dead_mask(words, edges, core.BROT)
    nch = int((w2 != words).sum())
    live_changed = int(((w2 != words) & ~dead).sum())
    w2p = w2
    if a.pad:
        if c4.vset(w2, a.k) != V:
            print("the block set V moved: no emission", flush=True)
            return
        w2p = c4.pad_table(w2, {a.k: V})
    an2 = core.analyze(w2p, a.k, 4)
    t_new = r63.wire_tensor(an2, an2.X0)
    an2.G, an2.X0, an2.out_level = an.G, X0, an.out_level
    t_mod = r63.wire_tensor(an2, X0)
    bad = sum(1 for x, y in zip(t_new, t_mod) if not np.array_equal(x, y))
    print(f"words changed {nch} (live {live_changed}); wires whose tensor differs from the model: {bad} of {len(t_new)}",
          flush=True)
    if bad or (a.pad and c4.kappa_wires(an2, V) != kap):
        print("model / kappa mismatch: no emission", flush=True)
        return
    if a.emit:
        os.makedirs(a.emit, exist_ok=True)
        text, st = c4.emit_c4(an2, prods, V) if a.pad else core.emit(an2, prods)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        np.savez(os.path.join(a.emit, f"table_k{a.k}.npz"), angle_words=w2, tree_edges=edges)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx: toffoli {st['toffoli']}, qubits {st.get('qubits')} (checked on "
              f"every address by the emitter); products {len(prods)}; formula price {lcost:+.5f} mHa^2; new words in "
              f"{a.emit}/table_k{a.k}.npz", flush=True)


if __name__ == "__main__":
    main()
