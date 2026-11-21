# C2 (22nd holder): PINNED COPY of tools/r59_chem_groups_tie.py made by tools/chem_ld2_e22mktie.py (C2_TIE_EXTRA)
#!/usr/bin/env python3
"""r59_chem_groups_tie.py: designer R59-CHEM-GROUPS (Nov 17). COLUMN TIES: FORCED RELATIONS OF RANK 2 FROM A SMALL POOL.

The placed-kill group engine (builder C2: tools/chem_ld2_e17kill.py on tools/chem_loaders_core.py; imported read-only,
nothing of theirs is edited) forces relations on a group only as COLUMN KILLS: the rank-one relations t (x) e_o,
76 equations in n pool unknowns per output, so a group needs a pool n >= 76 and the SECOND group (pool 39) gets
nothing: its 6 x 39 = 234 re-fan-out bits stay unused.

A TIE is the rank-two relation  column ta of output oa  ==  column tb of output ob  (as vectors over the 76 high rows).
It is forced by re-fanning pool products into oa and ob: 76 equations (75 independent: the padded table has one row
dependency) in 2 n unknowns, solvable from n >= 38 IF the pool planes are REGULAR for (ta, tb): the 2 x 2 matrix
[y_i[t]] (i = the product's two low vectors, t = ta, tb) invertible, so that a product gives both its high vectors.
Three ties on the six outputs of group 2 are three relations whose restriction to any fan-out has rank <= 2 (two low
forms), so 11 generic relations of W-perp can be added and EVERY fan-out keeps a plane: c = 14, products (90-14)/2 = 38
where the engine's law reads 39.

  group modes (--modes, one per plan group):  reg = planes regular for (ta, tb), 12 relations (3-dim kernels);
                                               tieN = N ties from the pool, then N + 11 relations;  - = the engine.
  ta, tb = the LAST two coordinates of the engine's kill order (never a kill coordinate: a tie makes the group's own
  planes silent on ta or tb for the fan-outs that hold one side of a pair); the planes of both groups are also kept
  off the silent side of the next column kill's coordinate, so the later kills keep their pool.

    python -I r59_chem_groups_tie.py <table.npz> --k 4 --plan "6,6,7:1:0:0,6:1:3:1" --modes reg,tie3,-,- \
        [--seeds 0-7] [--pad] [--tries 400] [--emit DIR] [--dump g.json]
Run in WSL with ~/qcc/ws/chem_own/venv/bin/python or on a box with ~/qcc/.venv_work/bin/python.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_ld2_e17kill as kl  # noqa: E402  (read-only: KillDecomp, parse_plan)

bits = core.bits


def span_points(basis):
    out = []
    for m in range(1, 1 << len(basis)):
        v = 0
        for q in bits(m):
            v ^= basis[q]
        out.append(v)
    return out


def pick_plane(ys, reg, crit, pref, rng):
    """two vectors (u, v) of span(ys): REGULAR for the local coordinates reg = (c1, c2) when reg is given (the 2 x 2
    matrix [[u_c1, u_c2], [v_c1, v_c2]] invertible), not silent on every coordinate of crit (u or v has the bit), and
    not silent on as many coordinates of pref as possible. None if no plane fits."""
    pts = span_points(ys)
    rng.shuffle(pts)
    best, bscore = None, -1
    for i, u in enumerate(pts):
        if reg:
            su = ((u >> reg[0]) & 1, (u >> reg[1]) & 1)
            if su == (0, 0):
                continue
        for v in pts[i + 1:]:
            if reg:
                sv = ((v >> reg[0]) & 1, (v >> reg[1]) & 1)
                if sv == (0, 0) or sv == su:
                    continue
            w = u | v
            if any(not (w >> t) & 1 for t in crit):
                continue
            score = sum((w >> t) & 1 for t in pref)
            if score > bscore:
                best, bscore = [u, v], score
                if score == len(pref):
                    return best
    return best


class TieCore(core.Decomp):
    """core.Decomp.solve_group with two more group modes (reg, tie). Sits under TauDecomp in the MRO of TieKill."""

    _mode = ""

    def _sel(self, keys, wb, dim_w):
        """the engine's (virtual) _select without the silent-plane filter: the keys given here are filtered already;
        going through self._select keeps the builder's dump / build hooks (Probe, PlanKill) alive."""
        sav = getattr(self, "cur_avoid", [])
        self.cur_avoid = []
        try:
            return self._select(keys, wb, dim_w)
        finally:
            self.cur_avoid = sav

    def apply_ties(self, outs, nt, ta, tb):
        n = len(self.prods)
        if n == 0:
            return None
        P = self.PM[:n]
        ca = core.rows_to_ints((P >> np.uint32(ta)) & 1)
        cb = core.rows_to_ints((P >> np.uint32(tb)) & 1)
        perm = list(outs)
        self.rng.shuffle(perm)
        forced = []
        for i in range(nt):
            oa, ob = perm[2 * i], perm[2 * i + 1]
            tgt = core.to_int((self.X[:, oa] >> np.uint32(ta)) & 1) ^ core.to_int((self.X[:, ob] >> np.uint32(tb)) & 1)
            sp = core.Span()
            idxs = list(range(2 * n))
            self.rng.shuffle(idxs)
            for j in idxs:
                sp.add(ca[j] if j < n else cb[j - n], 1 << j)
            r, c = sp.reduce(tgt)
            if r:
                self.tie_fail = getattr(self, "tie_fail", 0) + 1
                self.tie_rank = len(sp)
                return None
            for j in bits(c):
                p, o = (j, oa) if j < n else (j - n, ob)
                self.X[:, o] ^= P[p]
                self.prods[p][4] ^= 1 << o
            assert np.array_equal((self.X[:, oa] >> np.uint32(ta)) & 1, (self.X[:, ob] >> np.uint32(tb)) & 1)
            forced.append((oa, ta, ob, tb))
            self.tie_rank = len(sp)
        return forced

    def solve_group(self, outs, coords, dl=0, d=0, tag=""):
        mode = self._mode if tag == "generic" else ""
        if not mode:
            return core.Decomp.solve_group(self, outs, coords, dl, d, tag)
        an, nr, rng = self.an, self.nr, self.rng
        b = len(outs)
        ta, tb = self.tie_t
        live = None
        while dl or d:
            for shift in range(6):
                live = self.apply_kills(outs, coords, dl, d, shift)
                if live is not None:
                    break
            if live is not None:
                break
            if d:
                d -= 1
            else:
                dl -= 1
        if live is None:
            live = coords
        forced = []
        if mode.startswith("tie"):
            forced = self.apply_ties(outs, int(mode[3:]), ta, tb)
            if forced is None:
                raise RuntimeError("ties not solvable from the pool")
        lc = bits(live)
        a = len(lc)
        loc = {t: i for i, t in enumerate(lc)}
        gidx = {o: i for i, o in enumerate(outs)}
        Xg = self.X[:, outs]
        assert not (Xg & ~np.uint32(live)).any(), "group has content outside its live coordinates"
        rows = []
        for M in range(nr):
            v = 0
            for gi in range(b):
                x = int(Xg[M, gi])
                if x:
                    for ci, t in enumerate(lc):
                        if x >> t & 1:
                            v |= 1 << (gi * a + ci)
            rows.append(v)
        W = core.Span()
        wb = [v for v in rows if W.add(v)]
        dim_w, amb = len(wb), a * b
        rel = core.nullspace(wb, amb)
        fmask = []
        for (oa, t1, ob, t2) in forced:
            f = (1 << (gidx[oa] * a + loc[t1])) | (1 << (gidx[ob] * a + loc[t2]))
            assert all(bin(f & w).count("1") % 2 == 0 for w in wb), "a tie is not a relation of the residual"
            fmask.append(f)
        ylocal, yphase = [], []
        for q, ym in enumerate(an.ymasks):
            m = sum(1 << ci for ci, t in enumerate(lc) if ym >> t & 1)
            if m and m not in ylocal:
                ylocal.append(m)
                yphase.append(q)
        full_a = (1 << a) - 1
        bon = core._bonus(live, an.ymasks)
        if mode == "reg":
            c_tot, target = min(len(rel), bon - 1), (amb - min(len(rel), bon - 1) + 1) // 2
        else:
            c_tot = min(len(rel), len(fmask) + bon - 2 + int(os.environ.get("C2_TIE_EXTRA", "0")))   # C2 copy: + EXTRA
            target = max((amb - c_tot + 1) // 2, (dim_w + 1) // 2)
        ca, cb = loc[ta], loc[tb]
        crit = [loc[t] for t in self.kill_t[:1] if t in loc]       # the next column kill's coordinate
        pref = [loc[t] for t in self.kill_t[1:] if t in loc]
        best = None
        for trial in range(self.gtrials):
            ind = core.Span()
            chosen_rel = []
            for f in fmask:
                assert ind.add(f)
                chosen_rel.append(f)
            while len(chosen_rel) < c_tot:
                v = 0
                for r in rel:
                    if rng.random() < 0.5:
                        v ^= r
                if ind.add(v):
                    chosen_rel.append(v)
            keys, keys_all = {}, {}
            for s in range(1, 1 << b):
                gis = bits(s)
                cons = []
                for r in chosen_rel:
                    v = 0
                    for gi in gis:
                        v ^= (r >> (gi * a)) & full_a
                    if v:
                        cons.append(v)
                for p, ym in enumerate(ylocal):
                    if ym == full_a:
                        ys = core.nullspace(cons, a)
                    else:
                        ys = core.nullspace(cons + [1 << ci for ci in range(a) if not ym >> ci & 1], a)
                    if len(ys) < 2:
                        continue
                    if mode == "reg":                      # one regular plane of the kernel, nothing else
                        pl = pick_plane(ys, (ca, cb), crit, pref, rng)
                        if pl:
                            keys[(s, p)] = [(y, sum(y << (gi * a) for gi in gis)) for y in pl]
                        continue
                    # tie mode: the whole kernel, a plane that is not silent on ta first (the later column kill)
                    pl = pick_plane(ys, None, crit, pref, rng)
                    full = [(y, sum(y << (gi * a) for gi in gis)) for y in ys]
                    if pl:
                        sp2 = core.Span()
                        lst = []
                        for y in pl + ys:
                            if sp2.add(y):
                                lst.append((y, sum(y << (gi * a) for gi in gis)))
                        keys[(s, p)] = lst
                        keys_all[(s, p)] = lst
                    else:
                        keys_all[(s, p)] = full
            res = self._sel(keys, wb, dim_w)
            if (res is None or res[0] > target) and keys_all:
                res2 = self._sel(keys_all, wb, dim_w)
                if res2 is not None and (res is None or res2[0] < res[0]):
                    res = res2
            if res is not None and (best is None or res[0] < best[0]):
                best = res
            if best is not None and best[0] <= target:
                break
        if best is None:
            raise RuntimeError("no feature set spans the group")
        _, chosen = best
        fs = core.Span()
        for idx, (_, _, flat) in enumerate(chosen):
            ok = fs.add(flat, 1 << idx)
            assert ok
        Acoef = [0] * len(chosen)
        for M, v in enumerate(rows):
            r, c = fs.reduce(v)
            assert r == 0
            for idx in bits(c):
                Acoef[idx] |= 1 << M
        bykey = {}
        for idx, (key, y, _) in enumerate(chosen):
            if Acoef[idx]:
                bykey.setdefault(key, []).append(idx)
        n0 = len(self.prods)
        for (s, p), idxs in bykey.items():
            phase = yphase[p]
            sfull = sum(1 << outs[gi] for gi in bits(s))
            for i in range(0, len(idxs), 2):
                i1 = idxs[i]
                y1 = sum(1 << lc[ci] for ci in bits(chosen[i1][1]))
                if i + 1 < len(idxs):
                    i2 = idxs[i + 1]
                    y2 = sum(1 << lc[ci] for ci in bits(chosen[i2][1]))
                    self.add_product(y1, Acoef[i1], y2, Acoef[i2], sfull, phase)
                else:
                    self.add_product(y1, Acoef[i1], 0, 0, sfull, phase)
        assert not self.X[:, outs].any(), "group not cleared"
        made = len(self.prods) - n0
        self.log.append(f"{tag} b={b} a={a} kills=({dl},{d}) mode={mode} dimW={dim_w} amb={amb} c={c_tot} "
                        f"ties={len(fmask)} -> {made}")
        return made


class TieKill(kl.KillDecomp, TieCore):
    """KillDecomp (plan, placed kills, closure) over TieCore: MRO TieKill, KillDecomp, TauDecomp, TieCore, Decomp."""

    def __init__(self, an, seed, effort, plan, avoid_pos, modes, V=None, kap=0, tries=1, gtrials=40):
        super().__init__(an, seed, effort, plan, avoid_pos, V, kap, tries)
        self.modes = modes
        self.kill_t = [self.order[0], self.order[1]]     # the later groups' column kills (shift 0, shift 1)
        self.tie_t = (self.order[-1], self.order[-2])    # the ties' two columns: never a kill coordinate
        self.gtrials = gtrials

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
        self.n_low = len(self.prods)
        gen = list(classes.get(3, []))
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
        if gen:
            plan = core.make_planner(an, allc)
            while gen:
                b, dl, d = plan(len(gen), len(self.prods))
                self.solve_group(gen[:b], allc, dl, d, tag="generic")
                gen = gen[b:]
        assert not self.X.any()
        return self.prods


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--modes", required=True)
    ap.add_argument("--avoid", default="0,1")
    ap.add_argument("--seeds", default="0-3")
    ap.add_argument("--effort", type=int, default=3)
    ap.add_argument("--pad", action="store_true")
    ap.add_argument("--tries", type=int, default=400)
    ap.add_argument("--gtrials", type=int, default=40)
    ap.add_argument("--emit", default=None)
    ap.add_argument("--stop", type=int, default=0, help="stop at the first seed with at most this many products")
    ap.add_argument("--noclose", action="store_true", help="with --pad: the padded table without the closure (counts)")
    a = ap.parse_args()
    core.J, core.BROT = 5, 7
    words = np.asarray(np.load(a.table)["angle_words"], np.int64)
    V, kap = None, 0
    if a.pad:
        import chem_ld2_c4 as c4
        import chem_ld2_e19grp as eg
        V, words = eg.block_set(words, a.k)
    an = core.analyze(words, a.k, 4)
    if a.pad:
        kap = c4.kappa_wires(an, V)
    plan = kl.parse_plan(a.plan)
    modes = a.modes.split(",")
    avoid = [int(x) for x in a.avoid.split(",") if x.strip()]
    print(f"k {a.k}: {an.X0.shape[1]} outputs ({sum(1 for x in an.out_level if x > 2)} generic), nr {an.nr}, "
          f"phases {len(an.ymasks)}, infra {an.infra}, kill order {core.kill_order((1 << an.nl) - 1, an)[:3]}",
          flush=True)
    best = None
    for s in kl.rng_list(a.seeds):
        t0 = time.time()
        if a.noclose:
            kd = TieKill(an, s, max(2, a.effort), plan, avoid, modes, None, 0, 1, a.gtrials)
        else:
            kd = TieKill(an, s, max(2, a.effort), plan, avoid, modes, V, kap, a.tries if a.pad else 1, a.gtrials)
        try:
            prods = kd.run()
        except (RuntimeError, AssertionError) as exc:
            print(f"seed {s}: FAILED {type(exc).__name__} {exc} (tie fails {getattr(kd, 'tie_fail', 0)}, tie rank "
                  f"{getattr(kd, 'tie_rank', None)}); log {[x.split('generic ')[-1] for x in kd.log]} "
                  f"({time.time() - t0:.1f} s)", flush=True)
            continue
        n = len(prods)
        line = " | ".join(x.split("generic ")[-1] for x in kd.log if "generic" in x)
        t_a = kd.kill_t[0]
        sil = sum(1 for p in prods[:77] if not ((p[0] | p[2]) >> t_a) & 1)
        print(f"seed {s}: {n} products = {an.infra + n} ({an.infra} infra){' -1 with C-4' if a.pad else ''}; {line} "
              f"(tie rank {getattr(kd, 'tie_rank', None)}; silent on the first kill coordinate among the first 77: "
              f"{sil}; {time.time() - t0:.1f} s)", flush=True)
        if best is None or n < best[0]:
            best = (n, s, prods)
        if a.stop and n <= a.stop:
            break
    if best is None:
        return
    print(f"BEST: {best[0]} products, seed {best[1]}", flush=True)
    if a.emit and a.pad and not a.noclose:
        import chem_ld2_c4 as c4
        text, st = c4.emit_c4(an, best[2], V)
        os.makedirs(a.emit, exist_ok=True)
        with open(os.path.join(a.emit, f"angle_k_{a.k}.kmx"), "w") as fh:
            fh.write(text)
        print(f"emitted {a.emit}/angle_k_{a.k}.kmx: toffoli {st['toffoli']} (68-AND cache, checked by emit_c4)",
              flush=True)


if __name__ == "__main__":
    main()
