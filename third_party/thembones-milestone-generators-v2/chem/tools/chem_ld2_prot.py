"""chem_ld2_prot: the block engine with PROTECTED planar relations (ours).

Engine3 = chem_ld2_engine.Engine2 whose relation span of block b is P_b + R_b: P_b (the block's planar relations,
Engine2.degenerate) is fixed and always priced, R_b is what the search moves.  Every move of the parent works on
R_b only (opening a direction never drops a planar relation); complements, random relations and the
close-then-open exchange (swap2) are taken modulo P_b + R_b.  Measured reason (ws_chem_ld2.md, 14:55Z): the plain
seeded search keeps 20 of its 25 planar relations on k1.

Use: search_protected(an2, seed, effort, fill) -> (engine, price, history), like chem_own_engine.search.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
if os.environ.get("C2_ENGX") == "1":    # copy with CMAX = 6 (merged blocks), chem_ld2_engx
    import chem_ld2_engx as eng  # noqa: E402
else:
    import chem_own_engine as eng  # noqa: E402
import chem_ld2_engine as e2  # noqa: E402
from chem_loaders_core import Span, bits  # noqa: E402

_Base = eng.Engine            # the original class (chem_ld2_engine.install replaces eng.Engine later)


class Engine3(e2.Engine2):
    prot = None

    def _P(self, b):
        return self.prot[b] if self.prot is not None else []

    def set_protected(self, prot):
        # tenth holder: on short tables (B < 143) a block can hold more planar relations than the engine stores
        # (EMAX = 48: k0 at B 139 crashed in _sync); only then keep the first EMAX - 8 of them (room for the moves).
        # A block with at most EMAX planar relations (every table compiled before) is unchanged.
        cap = int(os.environ.get("C2_PROTCAP", eng.EMAX - 8))      # env C2_PROTCAP: the kept count (test knob)
        self.prot = [list(p) if len(p) <= eng.EMAX else list(p)[:cap] for p in prot]
        self.dirty = [True] * self.nblk

    # ---- spans: protected part first, then the movable part
    def _sync(self):
        for b in range(self.nblk):
            if not self.dirty[b]:
                continue
            bi = self.binfo[b]
            R = self._P(b) + self.R[b]
            self.Rm[b] = 0
            sp = Span()
            loc = self.loc[b]
            for j, v in enumerate(R):
                ok = sp.add(v)
                assert ok, "dependent relations"
                for idx in bits(v):
                    ci, o = loc[idx]
                    self.Rm[b, j, ci] |= 1 << o
            self.ecount[b] = len(R)
            m = np.ones(bi["amb"], np.bool_)
            for p in sp.piv:
                m[p] = False
            self.pivmask[b] = m
            self.dirty[b] = False
        free = np.concatenate(self.pivmask + [np.zeros(1, np.bool_)])
        idx = np.cumsum(free) - 1
        idx[~free] = -1
        self.freeidx = idx.astype(np.int64)
        self.r = int(free.sum())

    @property
    def codim(self):
        return [len(self._P(b)) + len(self.R[b]) for b in range(self.nblk)]

    def max_e(self, b):
        return min(eng.EMAX, len(self.binfo[b]["rel"])) - len(self._P(b))

    def block_dim(self, b):
        return self.binfo[b]["amb"] - len(self._P(b)) - len(self.R[b])

    def random_relation(self, b, sp):
        sp2 = Span()
        for v in self._P(b):
            sp2.add(v)
        for v, _ in list(sp.piv.values()):
            sp2.add(v)
        return _Base.random_relation(self, b, sp2)

    def _complement(self, b, R, cap=13):
        return _Base._complement(self, b, self._P(b) + list(R), cap)

    def swap2(self, cur, tot=None, nopen=3, nexact=6):
        """Parent's close-then-open exchange; directions that touch a protected relation are never opened."""
        rng = self.rng
        if tot is None:
            self.fresh_lines()
            tot = self.total_dim()
        snap = self.snapshot()
        b = rng.randrange(self.nblk)
        if len(self.R[b]) >= self.max_e(b):
            return cur, tot
        self.fresh_lines()
        Q = self._complement(b, self.R[b], cap=12)
        if not Q:
            return cur, tot
        cap = self.caps(b, Q)
        cap[0] = -1
        top = [int(z) for z in np.argsort(-cap, kind="stable")[:3] if cap[z] >= 0]
        if not top:
            return cur, tot
        a = top[rng.randrange(len(top))]
        added = self._vec(Q, a)
        self.set_block(b, self.R[b] + [added])
        self._sync()
        self.swaps2 += 1
        cands = []
        for b2 in range(self.nblk):
            p = len(self._P(b2))
            e2 = len(self.R[b2])
            if e2 < 2 or p + e2 > 22:
                continue
            hist = eng._open_hist(b2, self.pt_blk, self.pt_ymask, self.pt_outs, self.pt_phase, self.ncoord,
                                  self.Rm, self.ecount, self.nb)
            hist[0] = -1
            if p:
                lams = np.arange(len(hist))
                hist[(lams & ((1 << p) - 1)) != 0] = -1
            if b2 == b:
                hist[1 << (p + e2 - 1)] = -1
            for lam in np.argsort(-hist, kind="stable")[:nopen]:
                if hist[lam] < 0:
                    continue
                cands.append((int(hist[lam]), rng.random(), b2, int(lam) >> p))
        cands.sort(reverse=True)
        best = None
        for gain, _, b2, lam in cands[:nexact]:
            R2 = self.R[b2]
            j0 = (lam & -lam).bit_length() - 1
            opened = [R2[i] ^ (R2[j0] if (lam >> i) & 1 else 0) for i in range(len(R2)) if i != j0]
            self.set_block(b2, opened)
            p_, _ = self.price()
            self.exacts += 1
            t2 = self.total_dim()
            if best is None or (p_, -t2) < (best[0], -best[1]):
                best = (p_, t2, b2, opened)
            self.set_block(b2, R2)
            if p_ < cur:
                break
        if best is not None and (best[0] < cur or (best[0] == cur and best[1] > tot)):
            self.set_block(best[2], best[3])
            self._sync()
            if best[0] == cur:
                self.laterals += 1
            return best[0], best[1]
        self.restore(snap)
        return cur, tot


def map_relations(E_from, b_from, E_to, b_to, vecs):
    """Relation vectors of block b_from of engine E_from (masks over its cells) as masks over the cells of block
    b_to of E_to, through the (coordinate, output) of every cell."""
    bi_f, bi_t = E_from.binfo[b_from], E_to.binfo[b_to]
    out = []
    for v in vecs:
        w = 0
        for idx in bits(v):
            ci, o = E_from.loc[b_from][idx]
            t = bi_f["ts"][ci]
            w |= 1 << bi_t["cell"][(t, o)]
        out.append(w)
    return out


def search_protected(an2, seed=0, effort=1, fill=0, verbose=False, prot=None, start_R=None, skip_search=False):
    """chem_own_engine.search from P_b = every planar relation of block b (protected; or the given `prot`) and
    R_b random up to the engine's usual starting codimension (+ fill)."""
    e2.install()
    E = Engine3(an2, seed)
    E.prefer_degenerate = False
    if prot is None:
        deg = E.degenerate()
        prot = [list(deg[b]) for b in range(E.nblk)]
    E.set_protected(prot)
    E.R = [[] for _ in range(E.nblk)]
    if start_R is not None:                 # continue from given movable relations (no fill, no climb)
        E.restore([list(R) for R in start_R])
        cur, _ = E.price()
        if verbose:
            print("start from given spans: codim", E.codim, "price", cur, flush=True)
        if skip_search:                     # e.g. ruin-and-recreate follows: no greedy repairs at the frontier
            return E, cur, [(sum(E.codim), E.r, cur, E.evals, E.scans, 0.0)]
        return eng.search(an2, seed, effort, engine=E, verbose=verbose)
    evec = list(E.bd.default_e())
    for b in range(E.nblk):
        while E.codim[b] < evec[b] + fill:
            new = E.grow(b)
            if new is None:
                break
            E.set_block(b, new)
    cur, _ = E.price()
    cur = E.climb(100 * effort, cur, stop_at=(E.r + 1) // 2)
    return eng.search(an2, seed, effort, engine=E, verbose=verbose)


# ---- double steps (two relations at a time between states that carry one unpaired feature) -----------------
def level_up2(E, cur, trials=150, lateral=0.15, verbose=False):
    """Add the two best single relations (exact scans), then repair with the parent's moves until the price is
    (r + 1) // 2 again. From a state with one unpaired feature this costs one pair for two relations, instead of
    the zero-loss step a single relation needs there. Returns (price, reached?)."""
    rng = E.rng
    p = E.grow_scan(ntop=96)
    if p is None:
        return cur, False
    p = E.grow_scan(ntop=96)
    if p is None:
        return cur, False
    E._sync()
    goal = (E.r + 1) // 2
    if verbose:
        print(f"    double: two relations, price {p} goal {goal} codim {E.codim} lines {E.nlines}", flush=True)
    t = 0
    last = 0
    tot = None
    while p > goal and t - last < trials:
        if rng.random() < E.mix:
            q, tot = E.swap2(p, tot)
        else:
            b = rng.randrange(E.nblk)
            q = E.trial(b, p, lateral=lateral)
            tot = None
        if q < p:
            last = t
            if verbose:
                print(f"    double repair {t}: price {p} -> {q} codim {E.codim}", flush=True)
        p = q
        t += 1
    return p, p == goal


def search_double(E, an2, seed=0, effort=1, verbose=False, max_levels=10 ** 6):
    """From the engine's current state: single steps while the pairing is perfect with an even feature count,
    double steps (level_up2) from states with one unpaired feature. Keeps the best (price, then codimension)."""
    import time as _t
    t0 = _t.time()
    cur, _ = E.price()
    best = (cur, sum(E.codim), E.snapshot())
    hist = [(sum(E.codim), E.r, cur, E.evals, E.scans, 0.0)]
    fails = 0
    max_fail = 1 + effort
    levels = 0
    while fails < max_fail and levels < max_levels:
        base = E.snapshot()
        base_price = cur
        E._sync()
        if E.r % 2 == 0:
            p, ok = E.level_up(cur, quick=60 * effort, trials=150 * effort, pick=3 * fails, walk=0, verbose=verbose)
        else:
            p, ok = level_up2(E, cur, trials=150 * effort, verbose=verbose)
        hist.append((sum(E.codim), E.r, p, E.evals, E.scans, round(_t.time() - t0, 1)))
        if verbose:
            print("level", hist[-1], E.codim, "perfect" if ok else "deficient", flush=True)
        if ok:
            cur = p
            fails = 0
            levels += 1
            if (cur, -sum(E.codim)) < (best[0], -best[1]):
                best = (cur, sum(E.codim), E.snapshot())
        else:
            fails += 1
            E.restore(base)
            cur = base_price
    E.restore(best[2])
    return E, best[0], hist


def search_protected_double(an2, seed=0, effort=1, fill=0, verbose=False, prot=None, start_R=None):
    """search_protected, then double steps from its end state."""
    E, price, hist = search_protected(an2, seed, effort, fill=fill, verbose=verbose, prot=prot, start_R=start_R)
    if verbose:
        print("single-step search ended:", price, E.codim, "; double steps from here", flush=True)
    E2, price2, hist2 = search_double(E, an2, seed, effort, verbose=verbose)
    return E2, price2, hist + hist2


# ---- one-step lookahead (choose, among the relations that keep the pairing at its goal, one whose next step is
# ---- also good; measured reason: from the public k1 state minus two relations the greedy choice reaches 99 and
# ---- then cannot reach 100, which the public state shows to exist) ---------------------------------------------
def good_cands(E, need_nu, cap=256, ntop=96):
    """All (nu', killed, block, vector) over all blocks with nu' >= need_nu after adding that one relation (exact;
    E.prepare() must be current). Zero-loss candidates come from the engine's essential-line filter; otherwise the
    `ntop` candidates per block that kill the fewest lines are priced."""
    n, rho, ET, Rnp, pivrow, npidx = E.tr
    out = []
    for b in range(E.nblk):
        if len(E.R[b]) >= E.max_e(b):
            continue
        Q = E._complement(b, E.R[b])
        if not Q:
            continue
        tb = E.pt_blk
        sel = np.flatnonzero((tb[E.li[:n]] == b) | (tb[E.lj[:n]] == b)).astype(np.int64)
        nq = len(Q)
        Qm = E._qm(b, Q)
        T, Zr, Zc, dinv, h1, h2 = eng._scan_block(rho, ET, Rnp, pivrow, npidx, sel, E.li, E.lj, E.ls, E.lx,
                                                  E.pt_blk, E.pt_ymask, E.ncoord, E.cellid, E.freeidx, E.nb, b, Qm,
                                                  nq, eng.GF_EXP, eng.GF_LOG)
        killed = eng._scan_counts(h1, h2, 1 << nq)
        if need_nu >= rho // 2:
            ess = eng._essential(T, Zr, Zc, dinv)
            obs = np.zeros((len(dinv), (len(dinv) + 61) // 62 + 1), np.int64)
            buf = np.zeros((3, 8192), np.int64)
            cnt = eng._filter_cands(h1, h2, ess, obs, 1 << nq, E.vmax, buf)
            order = sorted(range(cnt), key=lambda q: (int(buf[1, q]), int(buf[2, q])))
            cands = [int(buf[0, q]) for q in order[:cap]]
        else:
            order = np.argsort(killed[1:], kind="stable") + 1
            cands = [int(z) for z in order[:ntop]]
        if not cands:
            continue
        nus = eng._scan_exact(T, Zr, Zc, dinv, h1, h2, np.array(cands, np.int64), rho, 10 ** 9, eng.GF_EXP,
                              eng.GF_LOG)
        E.exacts += len(cands)
        for nu, z in zip(nus, cands):
            if nu >= need_nu:
                out.append((int(nu), int(killed[z]), b, E._vec(Q, z)))
    return out


def best_next(E):
    """Score of a state for the next single relation. If the next step must lose no pair (the state carries one
    unpaired feature), the number of zero-loss relations available (0 = dead end; essential-line filter, cheap);
    otherwise 1 (an ordinary step, not discriminated)."""
    nu0 = E.prepare()
    r1 = E.r - 1
    need = r1 - (r1 + 1) // 2                     # nu needed after one more relation
    if need < nu0:
        return 1
    return len(good_cands(E, need, cap=256))


def level_up_la(E, cur, effort=1, ncand=16, verbose=False, fast_fail=False):
    """One more relation with the goal price, chosen by one-step lookahead; falls back to E.level_up."""
    base = E.snapshot()
    nu0 = E.prepare()
    r1 = E.r - 1
    goal = (r1 + 1) // 2
    need = r1 - goal
    cands = good_cands(E, need)
    if cands:
        cands.sort(key=lambda z: (-z[0], z[1]))
        cands = cands[:ncand]
        scored = []
        for nu, kl, b, v in cands:
            E.set_block(b, E.R[b] + [v])
            s = best_next(E)
            scored.append((s, -kl, E.rng.random(), b, v))
            E.restore(base)
        scored.sort(reverse=True)
        s, _, _, b, v = scored[0]
        E.set_block(b, E.R[b] + [v])
        p, _ = E.price()
        if verbose:
            print(f"    lookahead: {len(cands)} candidates keep the goal; chosen next-step slack {s}; price {p} goal "
                  f"{goal} codim {E.codim}", flush=True)
        if p == goal:
            return p, True
        E.restore(base)
    if fast_fail:
        return cur, False
    return E.level_up(cur, quick=60 * effort, trials=150 * effort, pick=0, walk=0, verbose=verbose)


def search_la(E, an2, seed=0, effort=1, verbose=False, fast_fail=False):
    """Level after level with level_up_la; keeps the best (price, then codimension)."""
    import time as _t
    t0 = _t.time()
    cur, _ = E.price()
    best = (cur, sum(E.codim), E.snapshot())
    hist = [(sum(E.codim), E.r, cur, E.evals, E.scans, 0.0)]
    fails = 0
    while fails < (1 if fast_fail else 1 + effort):
        base = E.snapshot()
        base_price = cur
        p, ok = level_up_la(E, cur, effort, verbose=verbose, fast_fail=fast_fail)
        hist.append((sum(E.codim), E.r, p, E.evals, E.scans, round(_t.time() - t0, 1)))
        if verbose:
            print("level", hist[-1], E.codim, "perfect" if ok else "deficient", flush=True)
        if ok:
            cur = p
            fails = 0
            if (cur, -sum(E.codim)) < (best[0], -best[1]):
                best = (cur, sum(E.codim), E.snapshot())
        else:
            fails += 1
            E.restore(base)
            cur = base_price
    E.restore(best[2])
    return E, best[0], hist


def search_protected_la(an2, seed=0, effort=1, fill=0, verbose=False, prot=None, start_R=None):
    """Protected start (planar relations + random fill to the default codimension, short climb), then lookahead
    levels (search_la) instead of the parent's greedy levels."""
    e2.install()
    E = Engine3(an2, seed)
    E.prefer_degenerate = False
    if prot is None:
        deg = E.degenerate()
        prot = [list(deg[b]) for b in range(E.nblk)]
    E.set_protected(prot)
    E.R = [[] for _ in range(E.nblk)]
    if start_R is not None:
        E.restore([list(R) for R in start_R])
    else:
        evec = list(E.bd.default_e())
        for b in range(E.nblk):
            while E.codim[b] < evec[b] + fill:
                new = E.grow(b)
                if new is None:
                    break
                E.set_block(b, new)
        cur, _ = E.price()
        E.climb(100 * effort, cur, stop_at=(E.r + 1) // 2)
    return search_la(E, an2, seed, effort, verbose=verbose)


def search_lns(E, an2, iters=20, drop=3, effort=1, seed=0, verbose=False, log=print):
    """Ruin and recreate from the engine's current state: drop `drop` random searched (unprotected) relations,
    climb again with lookahead levels (search_la), keep the result if it is cheaper, or as cheap with a higher
    codimension. Returns (engine at the best state, best price, history of (iteration, price, codim))."""
    import random as _r
    rng = _r.Random(7919 * seed + 17)
    cur, _ = E.price()
    best = (cur, sum(E.codim), E.snapshot())
    hist = [(0, cur, sum(E.codim))]
    for it in range(1, iters + 1):
        E.restore(best[2])
        for _ in range(drop):
            bs = [b for b in range(E.nblk) if E.R[b]]
            if not bs:
                break
            b = rng.choice(bs)
            R = list(E.R[b])
            R.pop(rng.randrange(len(R)))
            E.set_block(b, R)
        E2, p, h = search_la(E, an2, seed=seed + it, effort=effort, verbose=False, fast_fail=True)
        c = sum(E2.codim)
        hist.append((it, p, c))
        better = (p, -c) < (best[0], -best[1])
        log(f"lns {it}: price {p} codim {c} {E2.codim} (best {best[0]} at {best[1]}) {h[-1][-1]} s"
            + (" NEW BEST" if better else ""))
        if better:
            best = (p, c, E2.snapshot())
    E.restore(best[2])
    return E, best[0], hist


def search_lean(E, an2, seed=0, quick=60, trials=100, max_fail=1, verbose=False):
    """chem_own_engine.search's level loop with explicit budgets (for portfolios: one failed attempt at the
    frontier instead of 1 + effort). Same moves, same acceptance."""
    import time as _t
    t0 = _t.time()
    cur, _ = E.price()
    cur = E.climb(100, cur, stop_at=(E.r + 1) // 2)
    best = (cur, E.snapshot())
    hist = [(sum(E.codim), E.r, cur, E.evals, E.scans, 0.0)]
    fails = 0
    while fails < max_fail:
        base = E.snapshot()
        base_price = cur
        p, perfect = E.level_up(cur, quick=quick, trials=trials, pick=3 * fails, walk=0, verbose=verbose)
        hist.append((sum(E.codim), E.r, p, E.evals, E.scans, round(_t.time() - t0, 1)))
        if sum(E.codim) == sum(len(R) for R in base) + sum(len(p_) for p_ in E.prot or []):
            break
        if p < best[0] or (p == best[0] and perfect):
            best = (p, E.snapshot())
        if perfect:
            fails = 0
            cur = p
        else:
            fails += 1
            E.restore(base)
            cur = base_price
    E.restore(best[1])
    return E, best[0], hist


def search_protected_lean(an2, seed=0, quick=60, trials=100, max_fail=1, verbose=False):
    e2.install()
    E = Engine3(an2, seed)
    E.prefer_degenerate = False
    deg = E.degenerate()
    E.set_protected([list(deg[b]) for b in range(E.nblk)])
    E.R = [[] for _ in range(E.nblk)]
    evec = list(E.bd.default_e())
    for b in range(E.nblk):
        while E.codim[b] < evec[b]:
            new = E.grow(b)
            if new is None:
                break
            E.set_block(b, new)
    return search_lean(E, an2, seed, quick, trials, max_fail, verbose)
