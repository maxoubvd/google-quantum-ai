"""chem_ld2_engine: chem_own_engine.Engine with the plan extraction extended to a free top row.

With chem_ld2_basis.transform(..., free_top=True) the engine searches relation spans that need not annihilate the
top row of the mixed tensor (the row of the top high monomial, whose coefficients on the low monomials >= nv
are don't-cares because the table ends inside the last address block). At extraction the don't-cares are
chosen so that the true top row plus a don't-care vector lies in the span of the extracted features (one
linear solve); if no such choice exists the extraction fails and says so. Everything else (search, pricing,
lines) is the parent class unchanged.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
if os.environ.get("C2_ENGX") == "1":    # copy with CMAX = 6 (merged blocks), chem_ld2_engx
    import chem_ld2_engx as eng  # noqa: E402
else:
    import chem_own_engine as eng  # noqa: E402
from chem_loaders_core import Span, bits  # noqa: E402

_Base = eng.Engine
# sixteenth holder (Nov 16): env C2_SCANCAP = at most this many candidates per exact scan (zero-loss cap and ntop);
# on a WIRED loader (25 generic outputs) the uncapped exact scan stalls for 10+ minutes after level 48 [M]; unset =
# the parent's caps (256 / 96), unchanged
SCANCAP = int(os.environ.get("C2_SCANCAP", "0"))


SCANLINES = int(os.environ.get("C2_SCANLINES", "0"))


QGROW = int(os.environ.get("C2_QGROW", "0"))


class Engine2(_Base):
    def level_up(self, cur, quick=60, trials=150, pick=0, lateral=0.15, walk=60, verbose=False):
        # sixteenth holder: env C2_QGROW = n: a level step WITHOUT the exact scan, for loaders whose line table is too
        # large (wired loaders: 4-coordinate blocks): n random grows per block, the best one, then the parent's climb
        # and its swap2 / trial repairs until the pairing is perfect again (or `trials` steps without a gain)
        if QGROW <= 0:
            return super().level_up(cur, quick=quick, trials=trials, pick=pick, lateral=lateral, walk=walk,
                                    verbose=verbose)
        rng = self.rng
        base = self.snapshot()
        self._sync()
        cands = []
        for b in range(self.nblk):
            for _ in range(QGROW):
                new = self.grow(b)
                if new is None:
                    continue
                old = self.R[b]
                self.set_block(b, new)
                p, _ = self.price()
                cands.append((p, rng.random(), b, new))
                self.set_block(b, old)
        if not cands:
            return cur, False
        cands.sort()
        p, _, b, new = cands[min(pick, len(cands) - 1)]
        self.set_block(b, new)
        self._sync()
        goal = (self.r + 1) // 2
        if p > goal and quick:
            p = self.climb(quick, p, stop_at=goal)
        t, last, tot = 0, 0, None
        swap_only = os.environ.get("C2_QSWAP") == "1"     # no trial(): its exact scan stalls on large line tables
        while p > goal and t - last < trials:
            if swap_only or rng.random() < self.mix:
                q, tot = self.swap2(p, tot)
            else:
                q = self.trial(rng.randrange(self.nblk), p, lateral=lateral)
                tot = None
            if q < p:
                last = t
            p = q
            t += 1
        if p == goal:
            return p, True
        self.restore(base)
        return cur, False

    def grow_scan(self, ntop=96, pick=0):
        # sixteenth holder: env C2_SCANLINES = no exact scan when the line table is larger (the level step then
        # fails and the search retries through its quick path); unset = the parent's behaviour
        if SCANLINES > 0:
            self.prepare()
            if self.tr[0] > SCANLINES:
                return None
        return super().grow_scan(ntop=ntop, pick=pick)

    def scan(self, b, Q, need=None, ntop=64, cap=256, cands=None):
        # C2_SCANLINES also guards every exact scan (trial() and the walks call it): no candidate when the line table
        # is larger, so the move simply finds nothing instead of stalling for hours [M 21:10Z: wired / layout loaders]
        if SCANLINES > 0 and getattr(self, "tr", None) is not None and self.tr[0] > SCANLINES:
            return []
        if SCANCAP > 0:
            ntop, cap = min(ntop, SCANCAP), min(cap, SCANCAP)
            if cands is not None:
                cands = list(cands)[:SCANCAP]
        if os.environ.get("C2_SCANDBG") != "1":
            return super().scan(b, Q, need=need, ntop=ntop, cap=cap, cands=cands)
        import time as _t
        t0 = _t.time()
        n_ = self.tr[0] if getattr(self, "tr", None) is not None else -1
        print(f"  scan b={b} |Q|={len(Q)} need={need} lines={n_} rho={self.tr[1] if n_ >= 0 else -1} ...", flush=True)
        out = super().scan(b, Q, need=need, ntop=ntop, cap=cap, cands=cands)
        print(f"  scan b={b} done: {len(out)} results, {_t.time() - t0:.1f} s", flush=True)
        return out

    def _products(self, pairs, singles):
        an = self.an
        if not getattr(an, "free_top", False):
            return super()._products(pairs, singles)
        points = self.points
        binfo = self.binfo
        nb = an.X0.shape[1]
        nr = an.nr

        def g(u, o):
            return 1 << (u * nb + o)

        def fvec(i, s):
            yt = points[i][5]
            v = 0
            for u in yt:
                for o in bits(s):
                    v |= g(u, o)
            return v
        feats = []
        for n, (i, j, s) in enumerate(pairs):
            feats.append((fvec(i, s), ("p", n, 0)))
            feats.append((fvec(j, s), ("p", n, 1)))
        for n, (i, s) in enumerate(singles):
            feats.append((fvec(i, s), ("s", n)))
        fs = Span()
        for idx, (v, _) in enumerate(feats):
            if not fs.add(v, 1 << idx):
                return None
        coef = [0] * len(feats)
        X = an.X0
        for M in range(nr - 1):
            v = 0
            for o in range(nb):
                x = int(X[M, o])
                for u in bits(x):
                    v |= g(u, o)
            rr, c = fs.reduce(v)
            if rr:
                return None
            for idx in bits(c):
                coef[idx] |= 1 << M
        # top row: true row + don't-care vector in the span of the features
        t = 0
        for o in range(nb):
            for u in bits(int(an.top_true[o])):
                t |= g(u, o)
        rt, ct = fs.reduce(t)
        ds = Span()                          # residuals of the don't-care vectors, tagged by feature combos
        dl = [z for z in an.dlow if z]
        for o in range(nb):
            for z in dl:
                d = 0
                for u in bits(z):
                    d |= g(u, o)
                rd, cd = fs.reduce(d)
                if rd:
                    # carry the feature combination in the high part of the tag
                    r2, c2 = ds.reduce(rd, cd)
                    if r2:
                        ds.piv[r2.bit_length() - 1] = (r2, c2)
        r3, c3 = ds.reduce(rt, ct)
        if r3:
            self.log.append("free top row: no don't-care choice puts the top row in the features' span")
            return None
        for idx in bits(c3):
            coef[idx] |= 1 << (nr - 1)
        A = {key: coef[idx] for idx, (_, key) in enumerate(feats)}

        def phase_of(*pts):
            common = None
            for p in pts:
                ph = set(binfo[points[p][0]]["phases"])
                common = ph if common is None else common & ph
            return min(common)

        prods = []
        for n, (i, j, s) in enumerate(pairs):
            a1, a2 = A[("p", n, 0)], A[("p", n, 1)]
            if a1 and a2:
                prods.append([points[i][4], a1, points[j][4], a2, s, phase_of(i, j)])
            elif a1:
                prods.append([points[i][4], a1, 0, 0, s, phase_of(i)])
            elif a2:
                prods.append([points[j][4], a2, 0, 0, s, phase_of(j)])
        for n, (i, s) in enumerate(singles):
            a1 = A[("s", n)]
            if a1:
                prods.append([points[i][4], a1, 0, 0, s, phase_of(i)])
        return prods


    # ---- degenerate relations (relations of a block's pool that vanish on one of its points)
    def degenerate(self):
        """Per block: a basis of the span of all pool relations that vanish on some point of the block."""
        if getattr(self, "_deg", None) is not None:
            return self._deg
        out = []
        for b, bi in enumerate(self.binfo):
            ts, cell, rel = bi["ts"], bi["cell"], bi["rel"]
            span = Span()
            for ymask in range(1, 1 << len(ts)):
                yt = [ts[i] for i in bits(ymask)]
                # constraints: r(y (x) e_o) = 0 for every output o with all cells (t, o), t in yt
                cons = []
                for o in range(self.nb):
                    if all((t, o) in cell for t in yt):
                        v = 0
                        for t in yt:
                            v |= 1 << cell[(t, o)]
                        cons.append(v)
                # combinations c of rel with (sum c_i rel_i) . v = 0 for all v in cons
                m = len(rel)
                rows = []
                for v in cons:
                    row = 0
                    for i, r in enumerate(rel):
                        if bin(r & v).count("1") & 1:
                            row |= 1 << i
                    rows.append(row)
                from chem_loaders_core import nullspace
                for c in nullspace(rows, m):
                    vec = 0
                    for i in bits(c):
                        vec ^= rel[i]
                    span.add(vec)
            out.append([p[0] for p in span.piv.values()])
        self._deg = out
        return out

    def _complement(self, b, R, cap=13):
        """Parent's candidate directions, but the block's degenerate relations (not yet in span(R)) first."""
        if not getattr(self, "prefer_degenerate", False):
            return super()._complement(b, R, cap)
        sp = Span()
        for v in R:
            sp.add(v)
        first = []
        for v in self.degenerate()[b]:
            if len(first) < cap and sp.add(v):
                first.append(v)
        rest = super()._complement(b, list(R) + first, cap=max(0, cap - len(first))) if len(first) < cap else []
        return first + rest


def install():
    """Make chem_own_engine.search build Engine2 objects."""
    eng.Engine = Engine2
    if int(os.environ.get("C2_LINECAP", "0") or 0) > 0:   # seventeenth holder: capped line table (chem_ld2_e17cap.py)
        import chem_ld2_e17cap
        chem_ld2_e17cap.patch(eng)
    return eng


def search_seeded(an2, seed=0, effort=1, fill=0, verbose=False):
    """chem_own_engine.search from spans that start with every degenerate relation of each block (and `fill`
    random relations more), with degenerate directions preferred in every scan."""
    install()
    E = Engine2(an2, seed)
    E.prefer_degenerate = True
    deg = E.degenerate()
    spans = []
    for b in range(E.nblk):
        R = list(deg[b])[:E.max_e(b)]
        spans.append(R)
    E.restore(spans)
    # fill every block up to the engine's usual starting codimension (+ fill): starting from the degenerate
    # relations alone (5-7 per block) makes the line enumeration explode (OOM at 10 GB, measured 14:12Z)
    evec = list(E.bd.default_e())
    for b in range(E.nblk):
        while len(E.R[b]) < min(E.max_e(b), evec[b] + fill):
            new = E.grow(b)
            if new is None:
                break
            E.set_block(b, new)
    cur, _ = E.price()
    cur = E.climb(100 * effort, cur, stop_at=(E.r + 1) // 2)
    return eng.search(an2, seed, effort, engine=E, verbose=verbose)
