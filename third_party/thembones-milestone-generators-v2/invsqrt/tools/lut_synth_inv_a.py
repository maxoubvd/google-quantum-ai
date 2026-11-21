"""lut_synth design A ("consumed algebra") for the incumbent's 35-wire invsqrt table.

  lut_synth_inv_a.py [--seed S] [--split 12] [--cubic 0|1] [--late-limit 24] [--b3-outer 64] [--extra-b 1]
                     [--smooth 1] [--out file.json]

Phase 0: helper ANDs build W = all functions of address bits (7,6,5,4) [11 monomial wires] + bit3 * indicator(a) for
         `split` of the 16 classes a [one wire each], plus the local monomials x2x1, x2x0, x1x0 (+ x2x1x0 if --cubic).
Phase B: targets go into empty / fill wires with two-sided products (u0 + u')(v0 + v'), u', v' unknown in W.
Phase A: remaining targets consume the algebra wires with one-sided products (fixed u, free v), Gray search on the
         last multiplier. The last four wires hold x7x6, x7x5, x6x5, x7x6x5 (free erasure).
Every program is replayed from scratch and compared with the 35 truth tables (all 256 addresses)."""
import argparse
import json
import random
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
from lut_synth_core import (BASE, FULL, NB, Mach, Proj, Span, _search2, bits, gray_search, replay, solve_cols,
                            t_rem)
from lut_synth_fast import fsearch, fsearch2, split2

import os

TAB = os.path.expanduser("~/qcc/ws/invsqrt_new/scratch/inc_tab35.json")
SL = lambda bit: 1 + bit  # slot of an address bit


def all_subspaces(n, k):
    """Every k-dimensional subspace of GF(2)^n once, as a reduced echelon basis (list of k ints)."""
    import itertools
    for piv in itertools.combinations(range(n), k):
        free = [[q for q in range(p + 1, n) if q not in piv] for p in piv]
        for vals in itertools.product(*[range(1 << len(f)) for f in free]):
            basis = []
            for p, f, v in zip(piv, free, vals):
                row = 1 << p
                for j, q in enumerate(f):
                    if (v >> j) & 1:
                        row |= 1 << q
                basis.append(row)
            yield basis


def load_targets(path=TAB):
    D = json.load(open(path))
    return [int(t, 16) for t in D["tt"]], D["wires"]


_JOB = {}


def _b_worker(seed):
    S, ctx, r, n_outer, n_att = _JOB["job"]
    return S._b_attempts(ctx, r, n_outer, n_att, seed)


def _a_worker(rng_):
    return _JOB["scan"](rng_[0], rng_[1])


class SynthA:
    def __init__(self, T, seed=0, split=12, cubic=0, nw=39, log=print, merged=None):
        self.T = T
        self.rng = random.Random(seed)
        self.M = Mach(nw)
        self.log = log
        self.role = {}
        self.stats = {"helpers": 0, "B": [], "A": [], "tries": 0}
        M = self.M
        nxt = [NB]

        def new_wire(role):
            s = nxt[0]
            nxt[0] += 1
            self.role[s] = role
            return s

        self.new_wire = new_wire
        # W_16: monomials of bits 7,6,5,4
        hb = [7, 6, 5, 4]
        self.mono = {}
        for i in range(4):
            for j in range(i + 1, 4):
                w = new_wire("alg")
                M.AND(1 << SL(hb[i]), 1 << SL(hb[j]), w, "mono")
                self.mono[frozenset((hb[i], hb[j]))] = w
        for tri in ((7, 6, 5), (7, 6, 4), (7, 5, 4), (6, 5, 4)):
            w = new_wire("alg")
            M.AND(1 << self.mono[frozenset(tri[:2])], 1 << SL(tri[2]), w, "mono")
            self.mono[frozenset(tri)] = w
        w = new_wire("alg")
        M.AND(1 << self.mono[frozenset((7, 6, 5))], 1 << SL(4), w, "mono")
        self.mono[frozenset(hb)] = w
        self.w16_slots = [0] + [SL(b) for b in hb] + list(self.mono.values())
        self.d16 = {}
        for a in range(16):
            tt = sum(1 << x for x in range(256) if (x >> 4) == a)
            m = M.express(tt, self.w16_slots)
            assert m is not None
            self.d16[a] = (tt, m)
        order = list(range(16))
        self.rng.shuffle(order)
        if merged is not None:
            order = [a for a in order if a not in merged] + [a for a in order if a in merged]
        self.S = sorted(order[:split])
        self.x3wire = {}
        for a in self.S:
            w = new_wire("x3")
            M.AND(self.d16[a][1], 1 << SL(3), w, "x3d")
            self.x3wire[a] = w
        # local monomials (consumable fills)
        self.loc = []
        for (p, q) in ((2, 1), (2, 0), (1, 0)):
            w = new_wire("fill")
            M.AND(1 << SL(p), 1 << SL(q), w, "loc")
            self.loc.append(w)
        if cubic:
            w = new_wire("fill")
            M.AND(1 << self.loc[0], 1 << SL(0), w, "loc")
            self.loc.append(w)
        self.stats["helpers"] = M.ands
        self.first_free = nxt[0]
        self.smooth_slots = set()
        self.reserve = 0
        self.bias = 0
        self.slack = 4
        self.workers = 1
        self.smooth_fills = True

    # ---- bookkeeping
    def empties(self):
        return [s for s in range(self.first_free, self.M.ns) if s not in self.role]

    def fills(self):
        return [s for s, r in self.role.items() if r == "fill"]

    def classes(self):
        out = []
        x3 = BASE[SL(3)]
        for a in range(16):
            tt, m = self.d16[a]
            if a in self.x3wire:
                w = self.x3wire[a]
                out.append((tt & x3, 1 << w))
                out.append((tt & ~x3 & FULL, m ^ (1 << w)))
            else:
                out.append((tt, m))
        return out

    def rem(self):
        return t_rem(self.M.val, self.T)

    def settle(self, tw, tau, Trem):
        """After the products of a step: XOR the linear correction into the target wire so that it holds the
        target combination itself. The correction uses pool wires that may be consumed later (fills, local
        monomials, algebra wires), so it must be applied now."""
        M = self.M
        t = 0
        for j, v in enumerate(Trem):
            if (tau >> j) & 1:
                t ^= v
        d = M.val[tw] ^ t
        if d:
            mask = M.express(d, [s for s in range(M.ns) if s != tw and M.val[s] != 0])
            assert mask is not None, "linear correction is not in the pool span"
            M.CX(mask, tw)
        assert M.val[tw] == t

    # ---- smooth targets: fixed parts are local functions (closed local algebra, no junk outside the pool)
    def step_smooth(self, r, samples):
        M, rng = self.M, self.rng
        emp = self.empties()
        if not emp or any(self.role.get(w) != "fill" for w in self.loc):
            return False
        Trem = self.rem()
        R = len(Trem)
        tw = emp[0]
        pool_slots = [s for s in range(M.ns) if s != tw and M.val[s] != 0]
        pool = [(M.val[s], 1 << s) for s in pool_slots]
        locs = [SL(0), SL(1), SL(2)] + list(self.loc)
        cls = self.classes()
        NS = M.ns
        t0 = time.time()
        seen = set()
        subs = list(all_subspaces(len(locs), 2 * r))
        rng.shuffle(subs)
        for it, basis in enumerate(subs[:samples]):
            masks = [sum(1 << locs[b] for b in range(len(locs)) if (v >> b) & 1) for v in basis]
            seen.add(it)
            ftt = [M.form(m) for m in masks]
            cols = []
            cm = {}
            cid = NS
            J = 0
            for i in range(r):
                for ci, (dtt, dm) in enumerate(cls):
                    cols.append((dtt & ftt[2 * i], 1 << cid))
                    cm[(i, "B", ci)] = cid
                    cid += 1
                    cols.append((dtt & ftt[2 * i + 1], 1 << cid))
                    cm[(i, "A", ci)] = cid
                    cid += 1
                J ^= ftt[2 * i] & ftt[2 * i + 1]
            sol = solve_cols(pool, cols, J, Trem)
            if sol is None:
                continue
            tau, comb = sol
            n0 = M.ands
            for i in range(r):
                U, V = masks[2 * i], masks[2 * i + 1]
                used = False
                for ci, (dtt, dm) in enumerate(cls):
                    if (comb >> cm[(i, "A", ci)]) & 1:
                        U ^= dm
                        used = True
                    if (comb >> cm[(i, "B", ci)]) & 1:
                        V ^= dm
                        used = True
                if not used:
                    continue    # the junk of an unused product is a local monomial, i.e. a pool element
                M.AND(U, V, tw, "S")
                if i < r - 1 and self.smooth_fills:
                    e = [s for s in self.empties() if s != tw]
                    if e:
                        M.CX(1 << tw, e[0])
                        self.role[e[0]] = "fill"
                        self.smooth_slots.add(e[0])
            self.settle(tw, tau, Trem)
            self.role[tw] = "tgt"
            self.smooth_slots.add(tw)
            R2 = len(self.rem())
            assert R2 == R - 1, f"smooth step did not finish a target ({R} -> {R2})"
            self.stats["B"].append((r, len(cls), R, 7, 0, it + 1, round(time.time() - t0, 1)))
            self.log(f"  S r={r} R={R}->{R2} sample={it + 1} ands={M.ands} t={time.time() - t0:.1f}s")
            return True
        self.log(f"  S r={r} R={R}: none in {samples} samples ({len(seen)} subspaces) t={time.time() - t0:.1f}s")
        return False

    # ---- two-sided step, both fixed parts of the last product searched (numba two-level Gray walk)
    def usable_empties(self):
        e = self.empties()
        return e[:max(0, len(e) - self.reserve)]

    def _b_attempts(self, ctx, r, n_outer, attempts, seed):
        """Search part of a two-sided step; does not touch the machine (safe in forked workers).
        Returns a dict with the fixed parts and the solution, or None."""
        M = self.M
        rng = random.Random(seed)
        Trem, pool, nonalg, cls, J0 = ctx["Trem"], ctx["pool"], ctx["nonalg"], ctx["cls"], ctx["J0"]
        H, K, NS = len(cls), len(nonalg), M.ns
        tries = 0
        for att in range(attempts):
            fixed = []
            for i in range(r - 1):
                um = vm = 0
                while um == 0 or vm == 0 or um == vm:
                    um = sum(1 << s for s in nonalg if rng.random() < 0.5)
                    vm = sum(1 << s for s in nonalg if rng.random() < 0.5)
                fixed.append((um, vm))
            ftt = [(M.form(um), M.form(vm)) for um, vm in fixed]
            cols = []
            cid = NS
            colmap = {}
            Jfix = J0
            for i in range(r - 1):
                for ci, (dtt, dm) in enumerate(cls):
                    cols.append((dtt & ftt[i][0], 1 << cid))
                    colmap[(i, "B", ci)] = cid
                    cid += 1
                    cols.append((dtt & ftt[i][1], 1 << cid))
                    colmap[(i, "A", ci)] = cid
                    cid += 1
                Jfix ^= ftt[i][0] & ftt[i][1]
            P = Proj(Trem + [v for v, _ in pool] + [v for v, _ in cols])
            tabD = [dtt & M.val[s] for s in nonalg for (dtt, dm) in cls]
            tabQ = [M.val[p] & M.val[q] for p in nonalg for q in nonalg]
            D = P.pack(tabD).reshape(K, H, P.wd)
            Q = P.pack(tabQ).reshape(K, K, P.wd)
            bb = P.pack([Jfix])[0]
            fast = P.wd <= 2
            if fast:
                Dl, Dh = split2(D)
                Ql, Qh = split2(Q)
                bl, bh = split2(bb[None, :])
                P0 = Proj([v for v, _ in pool] + [v for v, _ in cols])      # quotient without the targets
                use0 = P0.wd <= 2
                if use0:
                    D0l, D0h = split2(P0.pack(tabD).reshape(K, H, P0.wd))
                    Q0l, Q0h = split2(P0.pack(tabQ).reshape(K, K, P0.wd))
                    b0l, b0h = split2(P0.pack([Jfix]))
                else:
                    D0l, D0h, Q0l, Q0h, b0l, b0h = Dl, Dh, Ql, Qh, bl, bh
            y, x = 1, 1
            while y <= n_outer:
                if fast:
                    yi, xi = fsearch2(Dl, Dh, Ql, Qh, bl[0], bh[0], P.c, use0, D0l, D0h, Q0l, Q0h, b0l[0], b0h[0],
                                      y, n_outer - y + 1, x, 1 << K)
                else:
                    yi, xi = _search2(D, Q, bb, y, n_outer - y + 1, x, 1 << K)
                yi, xi = int(yi), int(xi)
                if yi < 0:
                    tries += (n_outer - y + 1) << K
                    break
                tries += ((yi - y) << K) + xi
                gy, gx = yi ^ (yi >> 1), xi ^ (xi >> 1)
                vm = sum(1 << nonalg[k] for k in range(K) if (gy >> k) & 1)
                um = sum(1 << nonalg[k] for k in range(K) if (gx >> k) & 1)
                if um and vm and um != vm:
                    utt, vtt = M.form(um), M.form(vm)
                    cols2 = list(cols)
                    cm2 = dict(colmap)
                    c2 = cid
                    for ci, (dtt, dm) in enumerate(cls):
                        cols2.append((dtt & utt, 1 << c2))
                        cm2[(r - 1, "B", ci)] = c2
                        c2 += 1
                        cols2.append((dtt & vtt, 1 << c2))
                        cm2[(r - 1, "A", ci)] = c2
                        c2 += 1
                    sol = solve_cols(pool, cols2, Jfix ^ (utt & vtt), Trem)
                    if sol is not None:
                        return {"fixed": fixed + [(um, vm)], "tau": sol[0], "comb": sol[1], "cm": cm2,
                                "att": att + 1, "y": yi, "codim": P.c, "tries": tries}
                y, x = yi, xi + 1
                if x >= (1 << K):
                    y, x = yi + 1, 1
        return None

    def step_B(self, r, log2_tries, attempts=1, consume=False, allow_loc=True):
        M, rng = self.M, self.rng
        Trem = self.rem()
        R = len(Trem)
        consumed = None
        emp = self.usable_empties()
        if consume:
            a = rng.choice(sorted(self.x3wire))
            tw = self.x3wire.pop(a)
            consumed = a
        elif emp:
            tw = emp[0]
        else:
            fl = [s for s in self.fills() if s not in self.loc and s not in self.smooth_slots]
            if not fl and allow_loc:
                fl = self.fills()
            if not fl:
                return False
            tw = rng.choice(fl)
        pool_slots = [s for s in range(M.ns) if s != tw and M.val[s] != 0]
        pool = [(M.val[s], 1 << s) for s in pool_slots]
        allna = [s for s in pool_slots if s in (1, 2, 3, 4) or self.role.get(s) in ("fill", "tgt")]
        # enumeration basis of the fixed parts: generic (locally rich) wires first, then the address bits and the
        # local monomials; the smooth-phase wires are left out (they add no local richness, only duplicates)
        basic = set((1, 2, 3, 4)) | set(self.loc)
        genl = [s for s in allna if s not in basic and s not in self.smooth_slots]
        nonalg = genl + [s for s in allna if s in basic]
        gen = len(genl)
        cls = self.classes()
        H = len(cls)
        K = len(nonalg)
        ctx = {"Trem": Trem, "pool": pool, "nonalg": nonalg, "cls": cls, "J0": M.val[tw]}
        t0 = time.time()
        # measured log2 of the tries this step needs: codimension left by random complete products
        needs = []
        for _ in range(5):
            cols = []
            for i in range(r):
                um = sum(1 << s for s in nonalg if rng.random() < 0.5) or (1 << nonalg[0])
                vm = sum(1 << s for s in nonalg if rng.random() < 0.5) or (1 << nonalg[-1])
                ut, vt = M.form(um), M.form(vm)
                cols += [dtt & ut for (dtt, dm) in cls] + [dtt & vt for (dtt, dm) in cls]
            sp = Span()
            for v in Trem + [v for v, _ in pool] + cols:
                sp.add(v)
            needs.append(256 - sp.rank())
        need = max(0, sorted(needs)[2] - self.bias)
        if r == 3 and ((gen == 0 and need > 6) or need > log2_tries - 1):
            self.log(f"  B r={r} H={H} R={R} K={K} gen={gen}: skipped (needs 2^{need}, cap 2^{log2_tries})")
            if consumed is not None:
                self.x3wire[consumed] = tw
            return False
        budget = 1 << min(log2_tries, max(need + self.slack, 18))
        per_att = min(2 * K - 2, max(K, need - 2))      # log2 tries per attempt (fresh products each time)
        n_outer = max(1, min((1 << K) - 1, 1 << max(0, per_att - K)))
        n_att = max(1, min(attempts, -(-budget // (n_outer << K))))
        W = self.workers if n_att >= 2 * self.workers and n_att * (n_outer << K) >= (1 << 24) else 1
        self.log(f"  B r={r} H={H} R={R} K={K} gen={gen}: need 2^{need}, {n_att} attempts of "
                 f"2^{K + int(np.log2(n_outer))}, workers {W}")
        if W > 1:
            import multiprocessing as mp
            _JOB["job"] = (self, ctx, r, n_outer, -(-n_att // W))
            with mp.get_context("fork").Pool(W) as pool_:
                sol = None
                for res in pool_.imap_unordered(_b_worker, [rng.getrandbits(48) for _ in range(W)]):
                    if res is not None:
                        sol = res
                        pool_.terminate()
                        break
        else:
            sol = self._b_attempts(ctx, r, n_outer, n_att, rng.getrandbits(48))
        if sol is None:
            if consumed is not None:
                self.x3wire[consumed] = tw
            self.log(f"  B r={r} H={H} R={R} K={K} failed ({n_att} attempts) t={time.time() - t0:.1f}s")
            return False
        fixed, tau, comb, cm = sol["fixed"], sol["tau"], sol["comb"], sol["cm"]
        self.stats["tries"] += sol["tries"]
        for i in range(r):
            U, V = fixed[i]
            for ci, (dtt, dm) in enumerate(cls):
                if (comb >> cm[(i, "A", ci)]) & 1:
                    U ^= dm
                if (comb >> cm[(i, "B", ci)]) & 1:
                    V ^= dm
            M.AND(U, V, tw, "B")
            if i < r - 1:
                e = [s for s in self.usable_empties() if s != tw]
                if e:
                    M.CX(1 << tw, e[0])
                    self.role[e[0]] = "fill"
        self.settle(tw, tau, Trem)
        self.role[tw] = "tgt"
        R2 = len(self.rem())
        assert R2 == R - 1, f"step B did not finish a target ({R} -> {R2})"
        self.stats["B"].append((r, H, R, K, sol["codim"], sol["att"], round(time.time() - t0, 1)))
        self.log(f"  B r={r} H={H} R={R}->{R2} K={K} codim={sol['codim']} att={sol['att']} y={sol['y']} "
                 f"ands={M.ands} t={time.time() - t0:.1f}s")
        return True

    def step_B_old(self, r, outer_max, consume=False):
        M, rng = self.M, self.rng
        Trem = self.rem()
        R = len(Trem)
        consumed = None
        emp = self.empties()
        if consume:
            a = rng.choice(sorted(self.x3wire))
            tw = self.x3wire.pop(a)
            consumed = a
        elif emp:
            tw = emp[0]
        else:
            fl = self.fills()
            if not fl:
                return False
            tw = rng.choice(fl)
        pool_slots = [s for s in range(M.ns) if s != tw and M.val[s] != 0]
        pool = [(M.val[s], 1 << s) for s in pool_slots]
        nonalg = [s for s in pool_slots if s in (1, 2, 3, 4) or self.role.get(s) in ("fill", "tgt")]
        cls = self.classes()
        H = len(cls)
        K = len(nonalg)
        J0 = M.val[tw]
        NS = M.ns
        t0 = time.time()
        for outer in range(outer_max):
            fixed = []
            for i in range(r):
                um = vm = 0
                while um == 0 or vm == 0 or um == vm:
                    um = sum(1 << s for s in nonalg if rng.random() < 0.5)
                    vm = sum(1 << s for s in nonalg if rng.random() < 0.5)
                fixed.append((um, vm))
            ftt = [(M.form(um), M.form(vm)) for um, vm in fixed]
            cols = []
            cid = NS
            colmap = {}
            for i in range(r):
                for ci, (dtt, dm) in enumerate(cls):
                    if i < r - 1:
                        cols.append((dtt & ftt[i][0], 1 << cid))
                        colmap[(i, "B", ci)] = cid
                        cid += 1
                    cols.append((dtt & ftt[i][1], 1 << cid))
                    colmap[(i, "A", ci)] = cid
                    cid += 1
            Jfix = J0
            for i in range(r - 1):
                Jfix ^= ftt[i][0] & ftt[i][1]
            v0r = ftt[r - 1][1]
            P = Proj(Trem + [v for v, _ in pool] + [v for v, _ in cols])
            G = P.pack([dtt & M.val[s] for s in nonalg for (dtt, dm) in cls]).reshape(K, H, P.wd)
            Bm = P.pack([M.val[s] & v0r for s in nonalg])
            g0 = np.zeros((H, P.wd), dtype=np.uint64)
            b0 = P.pack([Jfix])[0]

            def verify(x):
                um = 0
                for k in range(K):
                    if (x >> k) & 1:
                        um |= 1 << nonalg[k]
                if um == 0 or um == fixed[r - 1][1]:
                    return None
                utt = M.form(um)
                cols2 = list(cols)
                cm2 = dict(colmap)
                c2 = cid
                for ci, (dtt, dm) in enumerate(cls):
                    cols2.append((dtt & utt, 1 << c2))
                    cm2[(r - 1, "B", ci)] = c2
                    c2 += 1
                sol = solve_cols(pool, cols2, Jfix ^ (utt & v0r), Trem)
                if sol is None:
                    return None
                return um, sol, cm2

            st = {}
            res = gray_search(G, Bm, g0, b0, K, verify, stats=st)
            self.stats["tries"] += st.get("tries", 0)
            if res is None:
                continue
            um_last, (tau, comb), cm = res
            fixed[r - 1] = (um_last, fixed[r - 1][1])
            for i in range(r):
                U, V = fixed[i]
                for ci, (dtt, dm) in enumerate(cls):
                    if (comb >> cm[(i, "A", ci)]) & 1:
                        U ^= dm
                    if (comb >> cm[(i, "B", ci)]) & 1:
                        V ^= dm
                M.AND(U, V, tw, "B")
                if i < r - 1:
                    e = [s for s in self.empties() if s != tw]
                    if e:
                        M.CX(1 << tw, e[0])
                        self.role[e[0]] = "fill"
            self.role[tw] = "tgt"
            R2 = len(self.rem())
            assert R2 == R - 1, f"step B did not finish a target ({R} -> {R2})"
            self.stats["B"].append((r, H, R, K, P.c, outer + 1, round(time.time() - t0, 1)))
            self.log(f"  B r={r} H={H} R={R}->{R2} K={K} codim={P.c} outer={outer + 1} "
                     f"ands={M.ands} t={time.time() - t0:.1f}s")
            return True
        if consumed is not None:
            self.x3wire[consumed] = tw
        return False

    # ---- one-sided step (consumes an algebra wire)
    def consumable(self):
        if self.x3wire:
            a = sorted(self.x3wire)[0]
            return self.x3wire[a], ("x3", a)
        for key, w in self.mono.items():
            if 4 in key and self.role.get(w) == "alg":
                return w, ("mono", key)
        # blocks without scratch: the last algebra wires (monomials of bits 7, 6, 5) are consumed too, highest degree first
        for key, w in sorted(self.mono.items(), key=lambda kv: -len(kv[0])):
            if self.role.get(w) == "alg":
                return w, ("mono", key)
        fl = self.fills()
        if fl:
            return fl[0], ("fill", None)
        return None, None

    def step_A(self, r, limit_log2):
        M, rng = self.M, self.rng
        Trem = self.rem()
        R = len(Trem)
        tw, what = self.consumable()
        if tw is None:
            return False
        pool_slots = [s for s in range(M.ns) if s != tw and M.val[s] != 0]
        pool = [(M.val[s], 1 << s) for s in pool_slots]
        pn = [s for s in pool_slots if s != 0]
        K = len(pn)
        NS = M.ns
        J = M.val[tw]
        t0 = time.time()
        us = []
        for i in range(r - 1):
            um = 0
            while um == 0 or um in us:
                um = sum(1 << s for s in pn if rng.random() < 0.5)
            us.append(um)
        cols = []
        colmap = {}
        cid = NS
        for i, um in enumerate(us):
            utt = M.form(um)
            for s in pn:
                cols.append((utt & M.val[s], 1 << cid))
                colmap[(i, s)] = cid
                cid += 1
        P = Proj(Trem + [v for v, _ in pool] + [v for v, _ in cols])
        # measured need: codimension left by a random last multiplier (median of 3); no search if it is hopeless
        needs = []
        for _ in range(3):
            um = sum(1 << s for s in pn if rng.random() < 0.5) or (1 << pn[0])
            utt = M.form(um)
            sq = Span()
            for v in P.pack([utt & M.val[s] for s in pn]):
                sq.add(sum(int(v[w]) << (64 * w) for w in range(P.wd)))
            needs.append(P.c - sq.rank())
        need = sorted(needs)[1]
        if need - 3 > limit_log2:
            self.log(f"  A r={r} R={R} codim={P.c}: skipped (needs about 2^{need}, cap 2^{limit_log2})")
            return False
        G = P.pack([M.val[a] & M.val[b] for a in pn for b in pn]).reshape(K, K, P.wd)
        Bm = np.zeros((K, P.wd), dtype=np.uint64)
        g0 = np.zeros((K, P.wd), dtype=np.uint64)
        b0 = P.pack([J])[0]

        def verify(x):
            um = 0
            for k in range(K):
                if (x >> k) & 1:
                    um |= 1 << pn[k]
            if um == 0 or um in us:
                return None
            utt = M.form(um)
            cols2 = list(cols)
            cm2 = dict(colmap)
            c2 = cid
            for s in pn:
                cols2.append((utt & M.val[s], 1 << c2))
                cm2[(r - 1, s)] = c2
                c2 += 1
            sol = solve_cols(pool, cols2, J, Trem)
            if sol is None:
                return None
            return um, sol, cm2

        st = {}
        if P.wd <= 2:
            Gl, Gh = split2(G)
            bl, bh = split2(b0[None, :])
            P0 = Proj([v for v, _ in pool] + [v for v, _ in cols])      # quotient without the targets
            if P0.wd <= 2:
                G0l, G0h = split2(P0.pack([M.val[a] & M.val[b] for a in pn for b in pn]).reshape(K, K, P0.wd))
                c0l, c0h = split2(P0.pack([J]))
            else:
                G0l = G0h = None
                c0l = c0h = [0]
            total = 1 << min(K, limit_log2)

            def scan(start, stop):
                tries = 0
                while start < stop:
                    cnt = min(1 << 24, stop - start)
                    idx = int(fsearch(Gl, Gh, bl[0], bh[0], start, cnt, G0l, G0h, c0l[0], c0h[0]))
                    if idx < 0:
                        tries += cnt
                        start += cnt
                        continue
                    tries += idx - start + 1
                    r_ = verify(idx ^ (idx >> 1))
                    if r_ is not None:
                        return r_, tries
                    start = idx + 1
                return None, tries

            expect = P.c - (K - 1)          # rough log2 of the tries needed
            res = None
            if self.workers > 1 and expect >= 24:
                import multiprocessing as mp
                chunk = 1 << 25
                _JOB["scan"] = scan
                with mp.get_context("fork").Pool(self.workers) as pool_:
                    for r_, tries in pool_.imap_unordered(_a_worker, ((a0, min(total, a0 + chunk))
                                                                       for a0 in range(1, total, chunk))):
                        st["tries"] = st.get("tries", 0) + tries
                        if r_ is not None:
                            res = r_
                            pool_.terminate()
                            break
            else:
                res, tries = scan(1, total)
                st["tries"] = tries
        else:
            res = gray_search(G, Bm, g0, b0, limit_log2, verify, stats=st)
        self.stats["tries"] += st.get("tries", 0)
        if res is None:
            self.log(f"  A r={r} R={R} codim={P.c} FAILED after 2^{limit_log2} ({time.time() - t0:.0f}s)")
            return False
        um_last, (tau, comb), cm = res
        us.append(um_last)
        for i, um in enumerate(us):
            V = 0
            for s in pn:
                if (comb >> cm[(i, s)]) & 1:
                    V |= 1 << s
            if V == 0:
                continue
            M.AND(um, V, tw, "A")
        if what[0] == "x3":
            del self.x3wire[what[1]]
        self.settle(tw, tau, Trem)
        self.role[tw] = "tgt"
        R2 = len(self.rem())
        assert R2 == R - 1, f"step A did not finish a target ({R} -> {R2})"
        self.stats["A"].append((r, R, P.c, st.get("tries", 0), round(time.time() - t0, 1)))
        self.log(f"  A r={r} R={R}->{R2} codim={P.c} tries=2^{np.log2(max(1, st.get('tries', 1))):.1f} "
                 f"ands={M.ands} t={time.time() - t0:.1f}s")
        return True


def finalize(M, role, T, log=print):
    """CX network that puts T[o] exactly on the o-th target wire; returns (target wires, leftover wires, erasable)."""
    tw = sorted(s for s, r in role.items() if r == "tgt")
    assert len(tw) == len(T)
    live = [s for s in range(M.ns) if M.val[s] != 0]
    masks = [M.express(t, live) for t in T]
    assert all(m is not None for m in masks), "a target is not in the span of the final wires"
    idx = {s: j for j, s in enumerate(tw)}
    rows = []
    rest = []
    for m in masks:
        a = 0
        for s in bits(m):
            if s in idx:
                a |= 1 << idx[s]
        rows.append(a)
        rest.append(m & ~sum(1 << s for s in tw))
    n = len(tw)
    ops = []
    A = list(rows)
    for c in range(n):
        if not (A[c] >> c) & 1:
            p = next(i for i in range(c + 1, n) if (A[i] >> c) & 1)
            A[c] ^= A[p]
            ops.append((c, p))
        for i in range(n):
            if i != c and (A[i] >> c) & 1:
                A[i] ^= A[c]
                ops.append((i, c))
    for (i, j) in reversed(ops):
        M.CX(1 << tw[j], tw[i])
    for o in range(n):
        if rest[o]:
            M.CX(rest[o], tw[o])
    assert all(M.val[tw[o]] == T[o] for o in range(n)), "final CX network does not reproduce the table"
    left = [s for s in range(NB, M.ns) if s not in tw and M.val[s] != 0]
    # free erasure: content must be a quadratic form over the wires still alive
    alive = [s for s in range(M.ns) if M.val[s] != 0]
    order = []
    todo = list(left)
    ok = True
    while todo:
        done = False
        for s in todo:
            others = [x for x in alive if x != s]
            sp = Span()
            for i, p in enumerate(others):
                for q in others[i:]:
                    sp.add(M.val[p] & M.val[q])
            if sp.contains(M.val[s]):
                order.append(s)
                alive.remove(s)
                todo.remove(s)
                done = True
                break
        if not done:
            ok = False
            break
    return tw, left, ok, order


def choose_merged(sm, split, nsmooth=6):
    """16-address blocks that stay merged: avoid the transition blocks of the smoothest target elements."""
    import itertools
    sm = sm[:nsmooth]
    best = None
    for comb in itertools.combinations(range(16), 16 - split):
        cs = set(comb)
        cost = sum(1.0 / (1 + j) for j, o in enumerate(sm) if cs & set(o["blocks"]))
        if best is None or cost < best[0]:
            best = (cost, comb)
    return list(best[1])


def drive(S, a, t0=None):
    """Phase sequence of design A on a prepared SynthA (helpers built). Returns True when every target is done."""
    t0 = t0 or time.time()
    nT = len(S.rem())
    import lut_synth_fast
    print(f"helpers {S.M.ands}, algebra H = {len(S.classes())}, empty wires {len(S.empties())}, "
          f"fills {len(S.fills())}, kernel {lut_synth_fast.KERNEL}", flush=True)
    if a.order == "smooth-first":
        if a.smooth and a.cubic:
            while S.step_smooth(1, a.s1):
                sys.stdout.flush()
            while S.step_smooth(2, a.s2):
                sys.stdout.flush()
    else:
        # generic targets first (they need the large target freedom), smooth targets afterwards in reserved wires
        S.reserve = a.reserve if (a.smooth and a.cubic) else 0
        for r in range(3, 7):
            if S.step_B(r, a.b_log, attempts=a.b_att, allow_loc=False):
                break
        sys.stdout.flush()
        while S.rem() and S.step_B(3, a.b_log, attempts=a.b_att, allow_loc=False):
            sys.stdout.flush()
        S.reserve = 0
        if a.smooth and a.cubic:
            while S.rem() and S.step_smooth(1, a.s1):
                sys.stdout.flush()
            while S.rem() and S.step_smooth(2, a.s2):
                sys.stdout.flush()
    rmin = 3
    while S.rem() and (S.empties() or S.fills()):
        done = False
        for r in range(rmin, 7):
            ok = S.step_B(r, a.b_log, attempts=a.b_att)
            if ok:
                done = True
                break
        if not done:
            print("phase B stuck", flush=True)
            break
        sys.stdout.flush()
    for _ in range(a.extra_b):
        if not S.rem() or not S.x3wire or not S.step_B(3, a.extra_log, attempts=1, consume=True):
            break
        sys.stdout.flush()
    print(f"after phase B: {S.M.ands} ANDs, targets done {nT - len(S.rem())}, t = {time.time() - t0:.0f}s", flush=True)
    if getattr(a, "save_b", ""):
        import pickle
        lg = S.log
        S.log = None
        S.new_wire = None
        pickle.dump(S, open(a.save_b, "wb"))
        S.log = lg
        print("saved phase-B state to", a.save_b, flush=True)
    return late(S, a)


def late(S, a):
    """Remaining targets consume the algebra wires. Per target, the cheapest step that is feasible within the
    search caps: two-sided with 4 products while the algebra still has late_b_min classes or more (little search
    at 26, growing below), one-sided with 4 products (search up to 2^late_limit; given up for the rest of the run
    after its first failure, because the need only grows as targets run out), then one-sided with 5 or 6."""
    a4_dead = False
    b4_dead = False
    bmin = getattr(a, "late_b_min", 24)
    while S.rem():
        ok = False
        if not b4_dead and S.x3wire and len(S.classes()) >= bmin:
            ok = S.step_B(4, min(a.b_log, 30), attempts=a.b_att, consume=True)
            if not ok:
                b4_dead = True
        if not ok and not a4_dead:
            ok = S.step_A(4, a.late_limit)
            if not ok:
                a4_dead = True
        if not ok:
            ok = S.step_A(5, a.late_limit) or S.step_A(6, a.late_limit)
        if not ok:
            print("phase A stuck", flush=True)
            break
        sys.stdout.flush()
    return not S.rem()


def add_design_args(ap):
    ap.add_argument("--split", type=int, default=12)
    ap.add_argument("--cubic", type=int, default=1)
    ap.add_argument("--late-limit", type=int, default=24)
    ap.add_argument("--extra-b", type=int, default=0, help="two-sided targets that consume an algebra wire")
    ap.add_argument("--smooth", type=int, default=1)
    ap.add_argument("--s1", type=int, default=3000)
    ap.add_argument("--s2", type=int, default=12000)
    ap.add_argument("--b-log", type=int, default=31)
    ap.add_argument("--b-att", type=int, default=4096)
    ap.add_argument("--nsmooth", type=int, default=6)
    ap.add_argument("--order", default="generic-first")
    ap.add_argument("--reserve", type=int, default=5)
    ap.add_argument("--need-bias", type=int, default=6)
    ap.add_argument("--b-slack", type=int, default=8, help="search budget = 2^(need + slack)")
    ap.add_argument("--extra-log", type=int, default=28)
    ap.add_argument("--workers", type=int, default=1, help="processes for the heavy two-sided searches")
    ap.add_argument("--late-b-min", type=int, default=24,
                    help="late phase: two-sided 4-product steps while the algebra has at least this many classes")
    ap.add_argument("--smooth-fills", type=int, default=0,
                    help="0: smooth steps leave their partial sums out of the pool (keeps empty wires for generic fills)")
    ap.add_argument("--save-b", default="", help="pickle the state after phase B")
    ap.add_argument("--resume", default="", help="continue a pickled phase-B state with the late phase")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    add_design_args(ap)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    T, wires = load_targets()
    t0 = time.time()
    merged = None
    sm_path = os.path.expanduser("~/qcc/ws/lut_synth/smooth16.json")
    if a.smooth and os.path.exists(sm_path) and a.split < 16:
        sm = json.load(open(sm_path))
        merged = choose_merged(sm, a.split, a.nsmooth)
        print("merged 16-address blocks:", merged, "smooth elements hit:",
              [j for j, o in enumerate(sm[:a.nsmooth]) if set(merged) & set(o["blocks"])], flush=True)
    if a.resume:
        import pickle
        S = pickle.load(open(a.resume, "rb"))
        S.log = print
        S.rng = random.Random(a.seed + 1000)
        S.workers = a.workers
        print(f"resumed {a.resume}: {S.M.ands} ANDs, remaining targets {len(S.rem())}", flush=True)
        late(S, a)
    else:
        S = SynthA(T, seed=a.seed, split=a.split, cubic=a.cubic, merged=merged)
        S.bias = a.need_bias
        S.slack = a.b_slack
        S.workers = a.workers
        S.smooth_fills = bool(a.smooth_fills)
        drive(S, a, t0)
    M = S.M
    ands = M.ands
    if not S.rem():
        tw, left, ok, order = finalize(M, S.role, T)
        M2 = replay(M.prog, M.nw)
        good = all(M2.val[tw[o]] == T[o] for o in range(len(T))) and M2.ands == ands
        print(f"RESULT seed {a.seed} split {a.split} cubic {a.cubic}: ANDs {ands} (helpers {S.stats['helpers']}, "
              f"B {sum(x[0] for x in S.stats['B'])} for {len(S.stats['B'])} targets, "
              f"A {sum(x[0] for x in S.stats['A'])} for {len(S.stats['A'])} targets); wires {M.nw}; "
              f"replay truth-table check on 256 addresses x 35 wires: {good}; leftover wires {len(left)} "
              f"free erasure {ok}; search tries 2^{np.log2(max(1, S.stats['tries'])):.1f}; "
              f"time {time.time() - t0:.0f}s", flush=True)
        if a.out:
            json.dump({"prog": [[s[0]] + [int(x) if not isinstance(x, str) else x for x in s[1:]] for s in M.prog],
                       "nw": M.nw, "ands": ands, "target_wires": tw, "leftover": left, "erase_order": order,
                       "inc_wires": wires, "stats": S.stats, "args": vars(a)}, open(a.out, "w"))
    else:
        print(f"INCOMPLETE: ANDs so far {ands}, remaining {len(S.rem())}", flush=True)


if __name__ == "__main__":
    main()
