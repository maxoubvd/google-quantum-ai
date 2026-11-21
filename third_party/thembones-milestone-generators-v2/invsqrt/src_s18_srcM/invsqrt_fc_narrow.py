#!/usr/bin/env python3
"""invsqrt_fc_narrow.py -- FC: narrow in-place adders on the kit's Circ (own code).

block(c, src, tgt, cin)  tgt += src + cin (mod 2^n) as a narrow block:
  forward  at most 2 carry wires. Carry c_(i+1) = MAJ(a_i, b_i, c_i) is one AND; it feeds the next carry and the
           sum bit, then it is X-measured (HMR): a measurement debt (-1)^(m_i c_i). The top carry goes into the
           top bit directly (one CCX, no wire).
  undo     generated now, stored in c.blocks[id], spliced in by Circ.mirror_blocks.
           It is X^n . add . X^n: the add re-creates the forward's carries, because the carries of ~T' + s + cin
           equal those of T + s + cin. The first placement of c_i writes sum bit i and repays the debt (a Z
           if m_i).
           Carries follow the checkpoint schedule at C clean wires:
             - blocks of carries; interior carries are removed free, top-down;
             - each block keeps its top carry;
             - after the top, the kept carries are ghosted top-down (X-measured) and each is repaired under
               its ghost bit (re-place the block interior, repay by Z + CZ, remove).
           Expected ANDs: L + 1 + sum over ghosts of (size - 1) / 2.
sched(L, C)  the best checkpoint schedule (block sizes) and its expected AND count.
  invsqrt_fc_narrow.py test        random-value check of every row shape on the kit's simulator
"""
import functools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ONE = "1"


def _w(x):
    return isinstance(x, int)


def _lit(x):
    if _w(x):
        return ("w", x, 0)
    return ("k", 1 if x == ONE else 0)


def _emit(c, ops, cond=None):
    for o in ops:
        if o[0] == "X":
            c.x(o[1], cond)
        else:
            c.cx(o[1], o[2], cond)


def _xor_lit(c, lit, w, cond=None):
    if lit[0] == "k":
        if lit[1]:
            c.x(w, cond)
    else:
        c.cx(lit[1], w, cond)
        if lit[2]:
            c.x(w, cond)


def _xor_src(c, a, w):
    if _w(a):
        c.cx(a, w)
    elif a == ONE:
        c.x(w)


def _form_b(a, bw, cr):
    """b-form (b on bw unchanged): x = a ^ c, y = b ^ c. ('and', xw, yw, ops) | ('const', v) | ('alias', w, p)."""
    if cr[0] == "k":
        kc = cr[1]
        if not _w(a):
            ka = 1 if a == ONE else 0
            if ka == kc:
                return ("const", kc)
            return ("alias", bw, 0)                  # MAJ(k, b, not k) = b
        ops = [("X", a), ("X", bw)] if kc else []
        return ("and", a, bw, ops)
    wc, p = cr[1], cr[2]
    if _w(a):
        ops = [("CX", wc, a), ("CX", wc, bw)]
        if p:
            ops += [("X", a), ("X", bw)]
        return ("and", a, bw, ops)
    ka = 1 if a == ONE else 0
    ops = [("CX", wc, bw)]
    if p:
        ops.append(("X", bw))
    if p ^ ka:
        ops.append(("X", wc))
    return ("and", wc, bw, ops)


def _form_s(a, sw, cr):
    """s-form (sum bit s = a ^ b ^ c already on sw): c' = c ^ x y with x = a ^ c, y = a ^ s."""
    if cr[0] == "k":
        kc = cr[1]
        if not _w(a):
            ka = 1 if a == ONE else 0
            if ka == kc:
                return ("const", kc)
            return ("alias", sw, 1)                  # c' = kc ^ ka ^ s = not s
        ops = [("CX", a, sw)]
        if kc:
            ops.append(("X", a))
        return ("and", a, sw, ops)
    wc, p = cr[1], cr[2]
    if _w(a):
        ops = [("CX", a, sw), ("CX", wc, a)]
        if p:
            ops.append(("X", a))
        return ("and", a, sw, ops)
    ka = 1 if a == ONE else 0
    ops = []
    if ka:
        ops.append(("X", sw))
    if p ^ ka:
        ops.append(("X", wc))
    return ("and", wc, sw, ops)


def fwd(c, src, tgt, cin):
    """tgt += src + cin with at most 2 carry wires. Returns the debts {carry index: HMR bit}."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    debts = {}
    if n == 0:
        return debts
    cur = _lit(cin)
    if n == 1:
        _xor_src(c, src[0], tgt[0])
        _xor_lit(c, cur, tgt[0])
        return debts
    cur_fresh = False
    deferred = None
    for i in range(n - 1):
        a, bw = src[i], tgt[i]
        f = _form_b(a, bw, cur)
        nxt, nfresh = None, False
        if i < n - 2:
            if f[0] == "const":
                nxt = ("k", f[1])
            elif f[0] == "alias":
                nxt = ("w", f[1], f[2])
            else:
                _, xw, yw, ops = f
                _emit(c, ops)
                t = c.alloc(1)[0]
                c.and_compute(xw, yw, t)
                _emit(c, ops[::-1])
                _xor_lit(c, cur, t)
                nxt, nfresh = ("w", t, 0), True
        else:
            top = tgt[n - 1]
            if f[0] == "const":
                if f[1]:
                    c.x(top)
            elif f[0] == "alias":
                c.cx(f[1], top)
                if f[2]:
                    c.x(top)
            else:
                _, xw, yw, ops = f
                _emit(c, ops)
                c.ccx(xw, yw, top)
                _emit(c, ops[::-1])
                _xor_lit(c, cur, top)
            _xor_src(c, src[n - 1], top)
        # sum bit i (deferred while the next carry is an alias of b_i)
        if f[0] == "alias" and i < n - 2:
            new_def = i
        else:
            new_def = None
            _xor_src(c, a, bw)
            _xor_lit(c, cur, bw)
        # retire c_i
        if cur_fresh:
            debts[i] = c.hmr(cur[1])
            c.free_([cur[1]])
        elif deferred is not None:
            c.x(tgt[deferred])                     # s = not b (the alias case: a, c constants, a != c)
        deferred = new_def
        cur, cur_fresh = nxt, nfresh
    assert deferred is None
    return debts


@functools.lru_cache(maxsize=None)
def sched(L, C):
    """checkpoint schedule: block sizes for L game carries at C wires; (sizes, expected ANDs incl. the top)."""
    if L == 0:
        return (), 1.0
    if C <= 0:
        return None, float("inf")
    best = {}

    def go(pos, j):
        key = (pos, j)
        if key in best:
            return best[key]
        res = (None, float("inf"))
        rem = L - pos
        if j + rem <= C:
            res = ((rem,), 0.0)
        for s in range(1, min(C - j, rem - 1) + 1):
            sub, cst = go(pos + s, j + 1)
            if sub is None:
                continue
            tot = cst + 0.5 * (s - 1)
            if tot < res[1] - 1e-12:
                res = ((s,) + sub, tot)
        best[key] = res
        return res

    sizes, extra = go(0, 0)
    if sizes is None:
        return None, float("inf")
    return sizes, L + 1 + extra


def _moves_from_sizes(L, sizes):
    moves = []
    pos = 0
    cps = []
    if L == 0:
        return [("top",)]
    for j, s in enumerate(sizes):
        for r in range(pos + 1, pos + s + 1):
            moves.append(("place", r))
        if j == len(sizes) - 1:
            moves.append(("top",))
            for r in range(pos + s, pos, -1):
                moves.append(("remove", r))
        else:
            for r in range(pos + s - 1, pos, -1):
                moves.append(("remove", r))
            cps.append(pos + s)
        pos += s
    for j in range(len(cps) - 1, -1, -1):
        cp = cps[j]
        prev = cps[j - 1] if j > 0 else 0
        moves.append(("remove", cp) if cp == prev + 1 else ("ghost", cp))
    return moves


_R32 = {}


def _r32_mainline(L, C, solve=False):
    """main line (every ghost returns 0) of R32's optimal adaptive policy (tools/r32_inv-clean_pebble.py, read-only
    import), as moves in game indices. Builds read the cache only (solve=False); the cache builder solves."""
    import importlib.util
    import json
    import re
    cpath = os.path.join(HERE, "invsqrt_fc_narrow_r32.json")
    if "cache" not in _R32:
        _R32["cache"] = json.load(open(cpath)) if os.path.exists(cpath) else {}
    key = "%d,%d" % (L, C)
    if key in _R32["cache"]:
        got = _R32["cache"][key]
        return None if got is None else [tuple(mv) for mv in got]
    if not solve:
        return None
    if "m" not in _R32:
        spec = importlib.util.spec_from_file_location("r32peb", os.path.join(HERE, "r32_inv-clean_pebble.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        _R32["m"] = mod
    M = _R32["m"]
    val, V, tabs = M.solve(L + 2, C, False)
    if not val < 1e6:
        return None
    moves = []
    for name, _w_, _v, _st in M.trace(L + 2, C, False, V, tabs, 0):
        k = re.search(r"c(\d+)", name)
        if name.startswith("place"):
            moves.append(("place", int(k.group(1))))
        elif name.startswith("remove"):
            moves.append(("remove", int(k.group(1))))
        elif name.startswith("GHOST"):
            moves.append(("ghost", int(k.group(1))))
        elif name == "top":
            moves.append(("top",))
        else:
            return None                     # re-place / repay on the main line: not used by the static form
    return moves


@functools.lru_cache(maxsize=None)
def _peb(ell, w):
    """cheapest static repair chain: from the anchor (always on a wire) put node ell (relative) on a wire, repay,
    and clear every node again, with at most w nodes on wires at once. place costs 1 (needs node - 1 on a wire or
    node 1), remove is free (same need). Exact 0-1 BFS over (nodes on wires, repaid). Returns (places, moves) with
    moves ("place", r) / ("remove", r) / ("repay",), or None when w wires cannot reach ell (ell > 2^w - 1)."""
    import collections
    if ell == 0:
        return 0, (("repay",),)
    if w <= 0:
        return None
    start, goal = (frozenset(), 0), (frozenset(), 1)
    dist = {start: 0}
    prev = {start: None}
    dq = collections.deque([start])
    while dq:
        st = dq.popleft()
        S, ph = st
        dcur = dist[st]
        if st == goal:
            break
        nxt = []
        if ph == 0 and ell in S:
            nxt.append(((S, 1), 0, ("repay",)))
        for i in range(1, ell + 1):
            if i == 1 or (i - 1) in S:
                if i in S:
                    nxt.append(((S - {i}, ph), 0, ("remove", i)))
                elif len(S) < w:
                    nxt.append(((S | {i}, ph), 1, ("place", i)))
        for ns, wt, mv in nxt:
            nd = dcur + wt
            if ns not in dist or nd < dist[ns]:
                dist[ns] = nd
                prev[ns] = (st, mv)
                if wt == 0:
                    dq.appendleft(ns)
                else:
                    dq.append(ns)
    if goal not in dist:
        return None
    mv, st = [], goal
    while prev[st] is not None:
        st, m_ = prev[st]
        mv.append(m_)
    return dist[goal], tuple(reversed(mv))


def run_moves(L, C, moves, ops=None):
    """execute a main line with static repairs: a ghost g is repaired under its bit just before its anchor (the
    highest carry below g on a wire) leaves, or at the end (anchor = the carry-in). A repair re-places the chain
    between anchor and ghost; when that chain does not fit the free wires, it is pebbled inside them (`_peb`).
    ops: callbacks (place(k, first, cond), remove(k, cond), ghost(k) -> bit, repay(k, bit), top()) or None to price
    only. Returns the expected AND count or None when a repair does not fit in C wires."""
    present = set()
    pending = {}
    cost = 0.0
    placed = set()

    def anchor(g):
        below = [p for p in present if p < g]
        return max(below) if below else 0

    def flush(leaving):
        nonlocal cost
        for g in sorted(pending, reverse=True):
            a = anchor(g)
            if a != leaving:
                continue
            need = g - a - 1
            if len(present) + need > C:
                pm = _peb(need, C - len(present))
                if pm is None:
                    return False
                pc, pmv = pm
            else:
                pc = need
                pmv = tuple([("place", r) for r in range(1, need + 1)] + [("repay",)] +
                            [("remove", r) for r in range(need, 0, -1)])
            cost += 0.5 * pc
            if ops is not None:
                b = pending[g]
                if "mark" in ops:
                    ops["mark"](("CBEG", b))
                for m_ in pmv:
                    if m_[0] == "place":
                        ops["place"](a + m_[1], False, b)
                    elif m_[0] == "remove":
                        ops["remove"](a + m_[1], b)
                    else:
                        ops["repay"](g, b)
                if "mark" in ops:
                    ops["mark"](("CEND", b))
            del pending[g]
        return True

    for mv in moves:
        if mv[0] == "place":
            k = mv[1]
            assert k - 1 == 0 or (k - 1) in present
            present.add(k)
            if len(present) > C:
                return None
            first = k not in placed
            placed.add(k)
            cost += 1.0
            if ops is not None:
                ops["place"](k, first, None)
        elif mv[0] == "top":
            assert L == 0 or L in present
            cost += 1.0
            if ops is not None:
                ops["top"]()
        else:
            k = mv[1]
            if not flush(k):
                return None
            if mv[0] == "remove":
                assert k - 1 == 0 or (k - 1) in present
                if ops is not None:
                    ops["remove"](k, None)
                present.discard(k)
            else:
                b = ops["ghost"](k) if ops is not None else True
                present.discard(k)
                pending[k] = b
    if not flush(0):
        return None
    assert not present and not pending and len(placed) == L
    return cost


@functools.lru_cache(maxsize=None)
def plan(L, C):
    """best static schedule for L game carries at C wires: (moves, expected ANDs incl. the top)."""
    cands = []
    sizes, _ = sched(L, C)
    if sizes is not None:
        mv = _moves_from_sizes(L, sizes)
        cands.append((run_moves(L, C, mv), 0, mv))
    if not cands and 1 <= L <= 15:
        mv = _r32_mainline(L, C)
        if mv is not None:
            cst = run_moves(L, C, mv)
            if cst is not None:
                cands.append((cst, 1, mv))
    if not cands:
        return None, float("inf")
    cst, _, mv = min(cands)
    return tuple(mv), cst


def game_add(c, src, tgt, cin, debts, C):
    """tgt += src + cin re-creating the carries (s-form) by the checkpoint schedule; repays the debts."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if n == 0:
        return 0.0
    carry = {0: _lit(cin)}
    if n == 1:
        _xor_src(c, src[0], tgt[0])
        _xor_lit(c, carry[0], tgt[0])
        return 0.0

    def write_sum(r):
        _xor_src(c, src[r], tgt[r])
        _xor_lit(c, carry[r], tgt[r])

    write_sum(0)
    F = 0
    while F + 1 <= n - 2:
        f = _form_s(src[F], tgt[F], carry[F])
        if f[0] == "and":
            break
        carry[F + 1] = ("k", f[1]) if f[0] == "const" else ("w", f[1], f[2])
        F += 1
        write_sum(F)
    L = n - 2 - F
    fresh = set()

    def place(r, first, cond=None):
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        assert f[0] == "and"
        _, xw, yw, ops = f
        t = c.alloc(1)[0]
        _emit(c, ops)
        if cond is None:
            c.and_compute(xw, yw, t)
        else:
            c.cand_compute(xw, yw, t, cond)
        _emit(c, ops[::-1])
        _xor_lit(c, carry[r - 1], t, cond)
        carry[r] = ("w", t, 0)
        fresh.add(r)
        if first:
            write_sum(r)
            if r in debts:
                c.z(t, debts[r])

    def remove(r, cond=None):
        t = carry[r][1]
        _xor_lit(c, carry[r - 1], t, cond)
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        _, xw, yw, ops = f
        _emit(c, ops)
        if cond is None:
            c.and_erase(xw, yw, t)
        else:
            c.cand_erase(xw, yw, t, cond)
            c.free_([t])
        _emit(c, ops[::-1])
        del carry[r]
        fresh.discard(r)

    def ghost(r):
        t = carry[r][1]
        g = c.hmr(t)
        c.free_([t])
        del carry[r]
        fresh.discard(r)
        return g

    def repay(r, g):
        lo = carry[r - 1]
        if lo[0] == "k":
            if lo[1]:
                c.neg(g)
        else:
            c.z(lo[1], g)
            if lo[2]:
                c.neg(g)
        f = _form_s(src[r - 1], tgt[r - 1], carry[r - 1])
        _, xw, yw, ops = f
        _emit(c, ops)
        c.cz(xw, yw, g)
        _emit(c, ops[::-1])

    def top():
        tw = tgt[n - 1]
        f = _form_s(src[n - 2], tgt[n - 2], carry[n - 2])
        if f[0] == "const":
            if f[1]:
                c.x(tw)
        elif f[0] == "alias":
            c.cx(f[1], tw)
            if f[2]:
                c.x(tw)
        else:
            _, xw, yw, ops = f
            _emit(c, ops)
            c.ccx(xw, yw, tw)
            _emit(c, ops[::-1])
            _xor_lit(c, carry[n - 2], tw)
        _xor_src(c, src[n - 1], tw)

    moves, cost = plan(L, C)
    if moves is None:
        raise ValueError("narrow undo: no schedule for L=%d at C=%d" % (L, C))
    cb = {"place": lambda k, first, cond: place(F + k, first, cond),
          "remove": lambda k, cond: remove(F + k, cond),
          "ghost": lambda k: ghost(F + k),
          "repay": lambda k, b: repay(F + k, b),
          "top": top,
          "mark": lambda m: c.ops.append(m)}
    got = run_moves(L, C, list(moves), cb)
    assert got is not None and abs(got - cost) < 1e-9
    assert not fresh
    return cost


def _cuccaro_alias(src, n, cin):
    """constant source bit 0 and constant carry-in that differ: c_1 = b_0 (a wire): start the chain at bit 1."""
    return n >= 2 and not _w(cin) and not _w(src[0]) and ((src[0] == ONE) != (cin == ONE))


def cuccaro_need(src, tgt, cin):
    """clean wires the Cuccaro-style undo needs: one per constant source bit below the top, and one for a constant
    carry-in (none for bit 0 and the carry-in when they alias c_1 = b_0)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    ws = [a for a in src if _w(a)]
    if len(set(ws)) != len(ws) or (_w(cin) and cin in ws) or set(ws) & set(tgt) or (_w(cin) and cin in tgt):
        return 10 ** 6                                  # carries live in the source wires: they must be distinct
    if _cuccaro_alias(src, n, cin):
        return sum(1 for a in src[1:n - 1] if not _w(a))
    return sum(1 for a in src[:n - 1] if not _w(a)) + (0 if _w(cin) else 1)


def cuccaro_add(c, src, tgt, cin, debts):
    """tgt += src + cin in place (Cuccaro, Draper, Kutin, Petrie Moulton 2004 form, own code): MAJ sweep with the
    carries held in the source wires, the top bit, UMA sweep back. Constant source bits / carry-in get clean
    wires. Re-creates the forward's carries c_1..c_(n-2) (c_(i+1) sits in a_i after MAJ_i) and repays the debts
    there. 2 CCX per position; returns the CCX count."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    mats = []

    def mat(x):
        if _w(x):
            return x
        w = c.alloc(1)[0]
        if x == ONE:
            c.x(w)
        mats.append((w, x))
        return w

    if n == 1:
        _xor_src(c, src[0], tgt[0])
        _xor_lit(c, _lit(cin), tgt[0])
        return 0
    i0 = 0
    if _cuccaro_alias(src, n, cin):                     # c_1 = b_0 on its own wire: the chain starts at bit 1,
        i0 = 1                                          # and s_0 = NOT b_0 is written after the sweeps
        c0 = tgt[0]
        A_ = [None] + [mat(src[i]) for i in range(1, n - 1)]
    else:
        c0 = mat(cin)
        A_ = [mat(src[i]) for i in range(n - 1)]
    chain = [c0] + A_[:]                                # chain[i] = carry c_i's wire after MAJ_(i-1)
    if i0 == 1:
        chain[1] = tgt[0]
    cnt = 0
    for i in range(i0, n - 1):                          # MAJ(chain[i], tgt[i], A_[i]): A_[i] -> c_(i+1)
        x, y, z = chain[i], tgt[i], A_[i]
        c.cx(z, y)
        c.cx(z, x)
        c.ccx(x, y, z)
        cnt += 1
        if (i + 1) in debts:
            c.z(z, debts[i + 1])
    a_top = src[n - 1]                                  # top: tgt[n-1] ^= a_(n-1) ^ c_(n-1)
    _xor_src(c, a_top, tgt[n - 1])
    c.cx(A_[n - 2] if n - 2 >= i0 else chain[n - 1], tgt[n - 1])
    if n - 2 < i0:                                      # n == 2 with the alias: the top used c_1 = b_0 directly
        pass
    for i in range(n - 2, i0 - 1, -1):                  # UMA(chain[i], tgt[i], A_[i])
        x, y, z = chain[i], tgt[i], A_[i]
        c.ccx(x, y, z)
        c.cx(z, x)
        c.cx(x, y)
        cnt += 1
    if i0 == 1:
        c.x(tgt[0])                                     # s_0 = a_0 xor b_0 xor c_0 = NOT b_0
    for w, x in reversed(mats):
        if x == ONE:
            c.x(w)
        c.free_([w])
    return cnt


def hybrid_add(c, src, tgt, cin, debts, C):
    """tgt += src + cin in place, Cuccaro form where wires are short and clean-wire carries where `C` allows (own
    code). Position i makes carry c_(i+1). A 'clean' position computes it into a fresh wire, c ^ (a ^ c)(b ^ c) (1
    AND), leaves b ^ c in the target, and in the downward pass erases the wire by measurement (CZ repair, free)
    while c_i is still on its wire, then writes the sum. Every other position is Cuccaro's MAJ / UMA with the carry
    in its source wire (2 CCX). Constant source bits are clean first (no materialisation), then wire positions from
    the top while clean wires remain. Re-creates the forward's carries and repays their debts. Returns the count."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if n == 1:
        _xor_src(c, src[0], tgt[0])
        _xor_lit(c, _lit(cin), tgt[0])
        return 0
    alias = _cuccaro_alias(src, n, cin)
    i0 = 1 if alias else 0
    budget = C
    mats = []
    if alias:
        c0 = tgt[0]
    elif _w(cin):
        c0 = cin
    else:
        c0 = c.alloc(1)[0]
        if cin == ONE:
            c.x(c0)
        mats.append((c0, cin))
        budget -= 1
    pos = list(range(i0, n - 1))
    clean = {i for i in pos if not _w(src[i])}
    budget -= len(clean)
    assert budget >= 0, "hybrid_add: constants need more clean wires than C"
    for i in sorted((i for i in pos if _w(src[i])), reverse=True):
        if budget <= 0:
            break
        clean.add(i)
        budget -= 1
    chain = {i0: c0}
    cnt = 0
    for i in pos:                                       # upward: chain[i] holds c_i exactly here
        x, y = chain[i], tgt[i]
        if i in clean:
            a = src[i]
            c.cx(x, y)                                  # y = b ^ c
            if _w(a):
                c.cx(x, a)                              # a' = a ^ c
                t = c.and_compute(a, y)
                c.cx(x, a)
            elif a == ONE:
                c.x(x)                                  # a' = 1 ^ c
                t = c.and_compute(x, y)
                c.x(x)
            else:
                t = c.and_compute(x, y)                 # a' = c
            c.cx(x, t)                                  # t = c ^ a'(b ^ c) = MAJ(a, b, c)
            cnt += 1
            chain[i + 1] = t
        else:
            z = src[i]
            c.cx(z, y)
            c.cx(z, x)
            c.ccx(x, y, z)                              # z = c_(i+1); x = c ^ a; y = b ^ a
            cnt += 1
            chain[i + 1] = z
        if (i + 1) in debts:
            c.z(chain[i + 1], debts[i + 1])
    _xor_src(c, src[n - 1], tgt[n - 1])                 # top: tgt[n-1] ^= a_(n-1) ^ c_(n-1)
    c.cx(chain[n - 1], tgt[n - 1])
    for i in reversed(pos):                             # downward
        x, y, nxt = chain[i], tgt[i], chain[i + 1]
        if i in clean:
            a = src[i]
            c.cx(x, nxt)                                # nxt = a'(b ^ c)
            if _w(a):
                c.cx(x, a)
                c.and_erase(a, y, nxt)
                c.cx(x, a)
                c.cx(a, y)                              # y = b ^ c -> a ^ b ^ c
            elif a == ONE:
                c.x(x)
                c.and_erase(x, y, nxt)
                c.x(x)
                c.x(y)
            else:
                c.and_erase(x, y, nxt)
        else:
            z = nxt
            c.ccx(x, y, z)
            c.cx(z, x)
            c.cx(x, y)
            cnt += 1
    if alias:
        c.x(tgt[0])                                     # s_0 = NOT b_0
    for w, xv in reversed(mats):
        if xv == ONE:
            c.x(w)
        c.free_([w])
    return cnt


def free_prefix(src, n, cin):
    """number of leading carries c_1..c_F that are constants or aliases (no AND, no wire)."""
    src = list(src) + [None] * (n - len(src))
    cur = _lit(cin)
    F = 0
    while F + 1 <= n - 2:
        a = src[F]
        if cur[0] != "k" or _w(a):
            break
        ka = 1 if a == ONE else 0
        cur = ("k", ka) if ka == cur[1] else ("w", -1, 1)
        F += 1
    return F


def feasible(src, tgt, cin, C):
    n = len(tgt)
    if C >= 2 and plan(max(0, n - 2 - free_prefix(src, n, cin)), C)[0] is not None:
        return True
    return C >= 2 and cuccaro_need(src, tgt, cin) <= C


def block(c, src, tgt, cin):
    """narrow block: forward now (between BBEG / BEND markers), undo stored in c.blocks."""
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    debts = fwd(c, src, tgt, cin)
    c.ops.append(("BEND", bid))
    C = c.cap - c.live
    save = c.ops
    c.ops = []
    for w in tgt:
        c.x(w)
    n = len(tgt)
    if plan(max(0, n - 2 - free_prefix(src, n, cin)), C)[0] is not None:
        cost = game_add(c, src, tgt, cin, debts, C)
    elif C > cuccaro_need(src, tgt, cin):               # no schedule, spare clean wires: Cuccaro form with some
        cost = float(hybrid_add(c, src, tgt, cin, debts, C))   # carries on clean wires (1 AND instead of 2 CCX)
    else:                                               # no checkpoint / R32 schedule at C: Cuccaro form
        cost = float(cuccaro_add(c, src, tgt, cin, debts))
    for w in tgt:
        c.x(w)
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((len(tgt), C, cost))


def _test():
    import numpy as np
    import invsqrt_fc_kit as K
    rng = np.random.default_rng(5)
    shapes = []
    nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    for n in range(1, nmax + 1):
        for kind in ("wires", "main", "tail", "constcin"):
            shapes.append((n, kind))
    nbad = 0
    for n, kind in shapes:
        for C in (2, 3, 4, 5, 6, 9, 16):
            c = K.Circ(n_fixed=0)
            c.cap = 10 ** 6
            tgt = c.alloc(n)
            ns = n
            if kind == "wires":
                srcw = c.alloc(ns)
                src = list(srcw)
                cin = c.alloc(1)[0]
            elif kind == "main":
                v = c.alloc(max(0, n - 3))
                src = ([ONE] + list(v) + [ONE] + [None])[:n]
                cin = None
            elif kind == "tail":
                v = c.alloc(max(0, n - 1))
                src = (list(v) + [ONE])[:n]
                cin = c.alloc(1)[0]
            else:
                v = c.alloc(max(0, n - 2))
                src = ([ONE, None] + list(v))[:n]
                cin = ONE
            live0 = c.live
            c.cap = live0 + C
            c.blocks = {}
            if not feasible(src, tgt, cin, C):
                continue
            fwd_start = len(c.ops)
            block(c, src, tgt, cin)
            fwd_ops = c.ops[fwd_start:]
            undo = K.Circ.mirror_blocks(fwd_ops, c.blocks)
            nsh = 4000
            sim = K.Sim(nsh, c.next + 4, seed=n * 31 + C)
            vals = {}
            for w in range(live0):
                vals[w] = rng.integers(0, 2, nsh).astype(bool)
                sim.q[w] = vals[w]
            T = sim.get_word(tgt)
            S = np.zeros(nsh, dtype=np.int64)
            for k, a in enumerate(src):
                if _w(a):
                    S |= sim.q[a].astype(np.int64) << k
                elif a == ONE:
                    S |= np.int64(1) << k
            ci = sim.q[cin].astype(np.int64) if _w(cin) else (1 if cin == ONE else 0)
            sim.run(fwd_ops)
            out = sim.get_word(tgt)
            want = (T + S + ci) % (1 << n)
            ok1 = np.array_equal(out, want)
            peak_ok = c.peak <= c.cap
            sim.run(undo)
            back = sim.get_word(tgt)
            ok2 = np.array_equal(back, T)
            others = [w for w in range(live0, c.next)]
            ok3 = not sim.q[others].any() if others else True
            ok4 = not sim.sign.any()
            for w in range(live0):
                if not np.array_equal(sim.q[w], vals[w]):
                    ok2 = False
            if not (ok1 and ok2 and ok3 and ok4 and peak_ok):
                nbad += 1
                print("BAD", n, kind, C, ok1, ok2, ok3, ok4, peak_ok, c.peak, c.cap)
    print("narrow test: %d shapes x C, bad %d" % (len(shapes), nbad))
    for L in range(0, nmax - 1):
        print("L %2d: " % L + " ".join("%5.1f" % plan(L, C)[1] for C in range(1, 10)))


def _r32cache(lmax=10):
    """precompute R32's main lines where the single-level family has no schedule (L <= lmax, C 2..4)."""
    import json
    out = {}
    for C in (2, 3, 4):
        for L in range(1, lmax + 1):
            if sched(L, C)[0] is not None:
                continue
            _R32.pop("cache", None)
            _R32["cache"] = {}
            mv = _r32_mainline(L, C, solve=True)
            ok = mv is not None and run_moves(L, C, mv) is not None
            out["%d,%d" % (L, C)] = [list(x) for x in mv] if ok else None
            print(L, C, "ok" if ok else "none", run_moves(L, C, mv) if ok else "", flush=True)
    p = os.path.join(HERE, "invsqrt_fc_narrow_r32.json")
    with open(p + ".tmp", "w") as f:
        json.dump(out, f)
    os.replace(p + ".tmp", p)


if __name__ == "__main__":
    if sys.argv[1] == "test":
        _test()
    elif sys.argv[1] == "r32cache":
        _r32cache(int(sys.argv[2]) if len(sys.argv) > 2 else 10)
