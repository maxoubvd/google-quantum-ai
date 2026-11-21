#!/usr/bin/env python3
"""r56_inv_39row_levj.py -- R56-INV-39ROW (designer): LEVEL ROWS of a divider tail in the J form.

Loaded by MY copy of the generator (~/qcc/ws/r56_inv_39row/src/invsqrt_fc_g8.py, st_divL, env R56_LEV=J) as module
`r56_levj` (copy this file there under that name). Nothing here touches a builder's file.

A tail row k with taper kk < k - 2 (level l = k - 2 - kk) adds l more bits of d at positions -l .. -1 under the window.
The J form, per level row:
  FORWARD (inside one BBEG / BEND block): the fraction register lives IN PLACE on the digit wires that the model
  reads at those positions (digit wire of row p + 1 at position -p); the low cells' carries (l ANDs) are made, the
  top one is the carry-in of the integer add (its own top carry done in place: n - 2 temporaries), then all are
  erased for free. Every carry that enters a wire non-linearly is also CNOT-copied into a pollution wire (one per
  polluted position; X for the window's lowest wire). At the tail's end the pollution wires clear the digit wires
  and the leftover (CNOTs), the affine steps are replayed (CNOTs), and the pollution wires are X-measured: debts.
  THE EXPLICIT UNDO (the mirror's side): (1) the fraction recursion alone in forward order (the same l ANDs per
  row): every carry gets its Z by the forward outcome bits (debts repaid), the top carries go into the leftover
  (CNOT) and are HELD, internal carries are held as far as the clean wires allow (plan below); (2) the integer
  rows are un-added in reverse order with the held top carry as carry-in; after each, the row's fraction step is
  taken back: a carry that is still alive is erased for free, an erased one is made again (one AND).
Count: forward N = sum of levels; undo N + (internal carries not held)."""
import itertools
import os

import invsqrt_fc_arith as A

ONE = A.ONE
LOG = {}


def _top_maj_inplace(c, a, b, cc, top):
    """top ^= MAJ(a, b, cc); b and cc wires, a a wire / ONE / None. One CCX, no temporary."""
    if A._is_w(a):
        c.cx(cc, a)
        c.cx(cc, b)
        c.ccx(a, b, top)
        c.cx(cc, top)
        c.cx(cc, a)
        c.cx(cc, b)
    elif a == ONE:
        c.ccx(b, cc, top)
        c.cx(b, top)
        c.cx(cc, top)
    else:
        c.ccx(b, cc, top)


def hi_add(c, src, tgt, cin):
    """tgt += src + cin (mod 2^n), cin a WIRE. Gidney carries with the top carry in place: n - 1 Toffoli-class
    gates, n - 2 temporaries. Self-contained (Circ.mirror inverts it)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    assert A._is_w(cin)
    if n == 1:
        A._sum_into(c, src[0], cin, tgt[0])
        return
    carries = [(cin, "in")]
    for i in range(n - 2):
        t, kind = A._maj_compute(c, src[i], tgt[i], carries[-1][0])
        assert kind == "w"
        carries.append((t, kind))
    _top_maj_inplace(c, src[n - 2], tgt[n - 2], carries[-1][0], tgt[n - 1])
    A._sum_into(c, src[n - 1], None, tgt[n - 1])
    for i in range(n - 2, -1, -1):
        A._sum_into(c, src[i], carries[i][0], tgt[i])
        if i >= 1:
            t, kind = carries[i]
            A._maj_erase(c, src[i - 1], tgt[i - 1], carries[i - 1][0], t, kind)


def _peak(recs, hold):
    """extra wires the explicit undo holds at its peak for a hold set of internal carries (k, q)."""
    held, peak = 0, 0
    for r in recs:
        if r["lvl"] == 0:
            continue
        peak = max(peak, held + r["lvl"])
        held += 1 + sum(1 for q in range(r["lvl"] - 1) if (r["k"], q) in hold)
    for r in reversed(recs):
        if r["lvl"] == 0:
            peak = max(peak, held + len(r["tgt"]) - 1)
            continue
        nh = sum(1 for q in range(r["lvl"] - 1) if (r["k"], q) in hold)
        peak = max(peak, held + max(0, len(r["hi_tgt"]) - 2))
        peak = max(peak, held + (r["lvl"] - 1 - nh))
        held -= 1 + nh
    return peak


def _plan(recs, budget):
    internal = [(r["k"], q) for r in recs for q in range(r["lvl"] - 1)]
    env = os.environ.get("R56_HOLD")
    if env is not None:
        hold = set()
        for tok in env.split(","):
            if tok:
                k_, q_ = tok.split(":")
                hold.add((int(k_), int(q_)))
        return hold, _peak(recs, hold), internal
    for size in range(len(internal), -1, -1):
        for hs in itertools.combinations(internal, size):
            p = _peak(recs, set(hs))
            if p <= budget:
                return set(hs), p, internal
    return set(), _peak(recs, set()), internal


def tail(c, K, u, v, Vb, sched, tw, x2, digit):
    """the rotation rows k = x2 + 2 .. len(sched) of st_divL; returns the final u."""
    assert x2 == 0
    nrow = len(sched)
    first = None
    for k in range(2, nrow + 1):
        if (k - 2) - sched[k - 1] > 0:
            first = k
            break
    if first is None:
        first = nrow + 1
    for k in range(2, first):                              # plain rows before the first level row (as today)
        kk = sched[k - 1]
        u = [u[tw - 1]] + u[0:tw - 1]
        src = Vb[kk + 1:]
        A.fadd(c, u[0], src + [None] * (tw - 1 - kk - len(src)), u[kk + 1:tw], cin=Vb[kk])
        c.cx(v[0], u[0])
        digit[-k] = u[0]
    if first > nrow:
        return u
    budget = getattr(c, "cap", 10 ** 9) - c.live
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    t0 = c.toffolis()
    Dw, Ipol, X, recs = {}, {}, None, []
    fpeak = 0
    for k in range(first, nrow + 1):
        kk = sched[k - 1]
        lvl = max(0, (k - 2) - kk)
        u = [u[tw - 1]] + u[0:tw - 1]
        ctl = u[0]
        if lvl == 0:
            src = Vb[kk + 1:]
            tgt = u[kk + 1:tw]
            src = src + [None] * (len(tgt) - len(src))
            A.fadd(c, ctl, src, tgt, cin=Vb[kk], force="gidney")
            recs.append(dict(k=k, lvl=0, ctl=ctl, src=src, tgt=tgt, cin=Vb[kk]))
        else:
            for p in range(lvl, 0, -1):
                if p not in Dw:
                    Dw[p] = u[k - 1 - p]
                assert Dw[p] == u[k - 1 - p]
            Fq = [Dw[p] for p in range(lvl, 0, -1)]              # q = 0 is the lowest position (-lvl)
            a_low = [Vb[kk + 1 + q] for q in range(lvl)]
            cin = Vb[kk]
            assert all(A._is_w(w_) for w_ in a_low) and A._is_w(cin)
            hi_src = Vb[kk + lvl + 1:]
            hi_tgt = u[k - 1:tw]
            live0 = c.live
            A.frame(c, ctl, Fq + hi_tgt)
            car = []
            cc = cin
            for q in range(lvl):
                t, kind = A._maj_compute(c, a_low[q], Fq[q], cc)
                assert kind == "w"
                car.append(t)
                cc = t
            fpeak = max(fpeak, c.live + max(0, len(hi_tgt) - 2))
            hi_add(c, hi_src, hi_tgt, car[-1])
            if X is None:
                X = c.alloc(1)[0]
            c.cx(car[-1], X)
            for q in range(lvl - 1):                              # car[q] enters position -(lvl - q - 1)
                p = lvl - (q + 1)
                if p not in Ipol:
                    Ipol[p] = c.alloc(1)[0]
                c.cx(car[q], Ipol[p])
            fpeak = max(fpeak, c.live)
            for q in range(lvl - 1, -1, -1):
                ccq = cin if q == 0 else car[q - 1]
                A._maj_erase(c, a_low[q], Fq[q], ccq, car[q], "w")
                A._sum_into(c, a_low[q], ccq, Fq[q])
            A.frame(c, ctl, Fq + hi_tgt)
            recs.append(dict(k=k, lvl=lvl, ctl=ctl, Fq=Fq, a_low=a_low, cin=cin, hi_src=hi_src, hi_tgt=hi_tgt))
            del live0
        c.cx(v[0], u[0])
        digit[-k] = u[0]
    lo = u[nrow - 1]                                             # the window's lowest wire = the tail's leftover
    for p, w in Ipol.items():                                    # pollution out of the digit wires
        c.cx(w, Dw[p])
    for r in recs:                                               # the affine steps back: the digit wires are digits
        if r["lvl"] == 0:
            continue
        for q in range(r["lvl"]):
            A._sum_into(c, r["a_low"][q], r["cin"] if q == 0 else None, r["Fq"][q])
    c.cx(X, lo)                                                  # the leftover is affine again
    bX = c.hmr(X)
    c.free_([X])
    bI = {}
    for p, w in Ipol.items():
        bI[p] = c.hmr(w)
        c.free_([w])
    c.ops.append(("BEND", bid))
    n_fwd = c.toffolis() - t0
    # ---------------------------------------------------------------- the explicit undo
    save = c.ops
    c.ops = []
    hold, upeak, internal = _plan(recs, budget)
    alive = {}
    for r in recs:                                               # (1) the fraction recursion, forward order
        if r["lvl"] == 0:
            continue
        k, lvl, ctl, Fq, a_low, cin = r["k"], r["lvl"], r["ctl"], r["Fq"], r["a_low"], r["cin"]
        c.cx(v[0], ctl)
        A.frame(c, ctl, Fq)
        car = []
        cc = cin
        for q in range(lvl):
            t, kind = A._maj_compute(c, a_low[q], Fq[q], cc)
            car.append(t)
            cc = t
        c.z(car[-1], bX)                                         # debts of the forward outcomes
        for q in range(lvl - 1):
            c.z(car[q], bI[lvl - (q + 1)])
        c.cx(car[-1], lo)
        for q in range(lvl - 1, -1, -1):
            ccq = cin if q == 0 else car[q - 1]
            if q == lvl - 1 or (k, q) in hold:
                alive[(k, q)] = car[q]
            else:
                A._maj_erase(c, a_low[q], Fq[q], ccq, car[q], "w")
            A._sum_into(c, a_low[q], ccq, Fq[q])
        A.frame(c, ctl, Fq)
        c.cx(v[0], ctl)
    n_fr = c.toffolis()
    for r in reversed(recs):                                     # (2) the integer rows back, last row first
        ctl = r["ctl"]
        c.cx(v[0], ctl)
        real = c.ops
        c.ops = []
        if r["lvl"] == 0:
            A.fadd(c, ctl, r["src"], r["tgt"], cin=r["cin"], force="gidney")
            sub = c.ops
            c.ops = real
            c.ops.extend(K.Circ.mirror(sub))
            continue                                             # (the control keeps the raw sign: the row is undone)
        k, lvl, Fq, a_low, cin = r["k"], r["lvl"], r["Fq"], r["a_low"], r["cin"]
        ck = alive.pop((k, lvl - 1))
        A.frame(c, ctl, r["hi_tgt"])
        hi_add(c, r["hi_src"], r["hi_tgt"], ck)
        A.frame(c, ctl, r["hi_tgt"])
        sub = c.ops
        c.ops = real
        c.ops.extend(K.Circ.mirror(sub))
        A.frame(c, ctl, Fq)                                      # the row's fraction step back
        car = [None] * lvl
        car[lvl - 1] = ck
        for q in range(lvl):
            ccq = cin if q == 0 else car[q - 1]
            A._sum_into(c, a_low[q], ccq, Fq[q])
            if q == lvl - 1:
                break
            if (k, q) in alive:
                car[q] = alive.pop((k, q))
            else:
                t, kind = A._maj_compute(c, a_low[q], Fq[q], ccq)   # an erased carry made again
                car[q] = t
        for q in range(lvl - 1, -1, -1):
            ccq = cin if q == 0 else car[q - 1]
            A._maj_erase(c, a_low[q], Fq[q], ccq, car[q], "w")
        A.frame(c, ctl, Fq)                                      # (the control keeps the raw sign: the row is undone)
    assert not alive
    n_undo = c.toffolis()
    c.blocks[bid] = c.ops
    c.ops = save
    LOG.update(levels=[(r["k"], r["lvl"]) for r in recs], fwd_block=n_fwd, undo_block=n_undo, undo_fr=n_fr,
               internal=len(internal), held=sorted(hold), undo_peak_extra=upeak, budget=budget,
               fwd_peak_live=fpeak, pollution_wires=1 + len(Ipol))
    LOG["_wires"] = dict(lo=lo, rows=[(r["k"], r["lvl"], r["ctl"], r.get("Fq"), r.get("hi_tgt") or r.get("tgt"),
                                       r.get("a_low"), r.get("cin")) for r in recs])
    return u
