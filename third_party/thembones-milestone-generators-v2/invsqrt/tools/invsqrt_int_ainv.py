"""invsqrt_int (G2, 2026-11-12): stage A^-1 (reverse rows of stage A) with fewer clean wires.

The rows module (tools/invsqrt_gen_rowsb3.py, the generator's pinned copy of workstream invsqrt_rows; imported,
never edited) undoes stage A's nine rows in joint blocks (_rows_rev_joint): inside a block every carry stays alive
until the block's top chain has run, so a block needs (row widths + carry-ins + chain length) clean wires. With
the 17 clean wires of today's layout its planner (plan_join) takes rows 0, 1 flat and blocks (2,3), (4,5),
(6,7,8): 71 Toffolis, against 63 with unlimited wires [M: 71 at 17-18, 68 at 19-22, 65 at 26, 63 at 40].

Here a joint block may run its top chain (the increment of the window bits above the block by the block's merged
digit) with p < L carry wires, in blocks, as the module's own bounded reverse rows do (_row_rev_core): inner
carries are erased for free right after their block's write; a block's top carry (boundary) stays for the next
block and is measured at the end and repaired by recomputing its block under the outcome (executed half the
time). Mean extra (L - p) / 2 Toffolis. plan_join_ext chooses the blocks and p. Everything else (rows, merges,
repairs of the forward outcomes, erasure order) is the module's code path, unchanged.
Model of the planner [E, same counting as plan_join]: 71 -> 70 at 17 wires, 71 -> 68.5 at 18, 78 -> 74 at 14.

Call (builder option ainv=1): bind(<name of the pinned rows module>) once, then
emit_A_inv_ext(c, t, win, c2, ctx, mP, pool) exactly where rows2p1.emit_A_inv(c, t, win, c2, ctx, mP, pool) is
called (same context from emit_A, same repairs)."""
from __future__ import annotations

import importlib

RB = None


def bind(name="invsqrt_gen_rowsb3"):
    global RB
    RB = importlib.import_module(name)
    return RB


def _pmin(n):
    """Smallest p for which plan_blocks(n, p) exists."""
    if n <= 0:
        return 0
    for p in range(1, n + 1):
        cap, q = p, 1
        while cap < n and p - q >= 1:
            q += 1
            cap += p - q + 1
        if cap >= n:
            return p
    return n


SHED = 1      # v2 (G2): a joint block may shed the carry-in of its widest row (0 = v1 behaviour)
_DBG = bool(__import__("os").environ.get("INVSQRT_INT_AINV_DEBUG"))


def plan_join_ext(specs, top_bit, pool):
    """As RB.plan_join, returning (first, last, p_chain, shed): p_chain None = the module's own block (chain with
    all its wires) or a flat row; an int = the block's top chain runs with that many carry wires. shed = 1: the
    widest row's carry-in is erased right after that row's writes and re-made, under the outcome of the carry
    above it, for that carry's erasure (mean +0.5, one wire less from the second row on)."""
    r = len(specs)
    a = [sp.m + 1 for sp in specs]
    cin = [1 if sp.has_cin else 0 for sp in specs]
    flat = [c_ + top_bit - 1 for c_ in cin]

    def price(i, j):
        if any(a[k] != a[k + 1] + 1 for k in range(i, j)):
            return None
        L = top_bit - 1 - a[i]
        p0 = sum(a[i:j + 1]) + sum(cin[i:j + 1]) + (j - i) + L
        best_ = None
        for s in range(0, (cin[i] if SHED else 0) + 1):
            kept, peak = 0, 0
            for k in range(i, j + 1):
                peak = max(peak, kept + a[k] + cin[k])
                kept += a[k] + cin[k] - (s if k == i else 0)
            for p in ([L] if L <= 0 else range(_pmin(L), L + 1)):
                if max(peak, kept + p) > pool:
                    continue
                pr = p0 + 0.5 * s + (L - p) / 2
                if best_ is None or pr < best_[0]:
                    best_ = (pr, None if p == L else p, s)
        return best_

    best = [(0, [])] * (r + 1)
    for i in range(r - 1, -1, -1):
        cand = (flat[i] + best[i + 1][0], [(i, i, None, 0)] + best[i + 1][1])
        for j in range(i + 1, r):
            pr = price(i, j)
            if pr is not None and pr[0] + best[j + 1][0] < cand[0]:
                cand = (pr[0] + best[j + 1][0], [(i, j, pr[1], pr[2])] + best[j + 1][1])
        best[i] = cand
    return best[0][1]


def _phase_of(T, mask, cond=None):
    """Apply (-1)^(content mask) (under the classical bit cond, or unconditionally: inside a controlled section
    the gates are conditioned by it). First the AND terms that a combination of live wires carries are cancelled
    with those wires (a Z on each: the phase of a wire's own content), so only the remaining AND terms need their
    factors as functions of the live wires (T.erase expands every term, which fails once a carry-in is gone)."""
    c = T.c
    gens = [(1, ("c", None))]                    # the constant (bit 0 of a mask)
    for w, m in T.f.items():
        gens.append((m, ("w", w)))               # a live wire: Z on it gives the phase of its whole content
    for v, (X, Y) in T.ands.items():
        try:
            A, ca = T._solve(X, None)
            B, cb = T._solve(Y, None)
        except AssertionError:
            continue                             # an AND term whose factors are not functions of live wires
        gens.append((v, ("a", (A, ca, B, cb))))
    basis = {}                                   # leading bit -> (vector, set of generator indices)
    for gi, (m, _) in enumerate(gens):
        vec, comb = m, {gi}
        while vec:
            h = vec.bit_length() - 1
            if h not in basis:
                basis[h] = (vec, comb)
                break
            vec ^= basis[h][0]
            comb = comb ^ basis[h][1]
    vec, used = mask, set()
    while vec:
        h = vec.bit_length() - 1
        assert h in basis, "content is not a function of the live wires (smart)"
        vec ^= basis[h][0]
        used ^= basis[h][1]
    for gi in sorted(used):
        kind, data = gens[gi][1]
        if kind == "c":
            c.neg(cond=cond)
        elif kind == "w":
            c.z(data, cond=cond)
        else:
            A, ca, B, cb = data
            T._phase_prod(A, ca, B, cb, cond)


def _erase_smart(T, w):
    """T.erase with _phase_of's reduction (same measurement, same freeing)."""
    c = T.c
    mask = T.f.pop(w)
    m = c.hmr(w)
    if _DBG:
        print(f"  [ainv] smart erase of wire {w}", flush=True)
    _phase_of(T, mask, m)
    c.bit_free(m)
    c.free(w)


def _chain_blocked(c, T, bits, d0, a, top_bit, mP, zfix_up, p):
    """bits (= reg[a:top_bit], already in the chain's frame) += d0 at bits[0]; carry out of bits[-1] repaired under
    mP (or dropped when mP is None); Z of the forward node outcomes zfix_up[position] on the carry into that
    position. At most p carry wires. Leaves no carry wire."""
    L = len(bits) - 1                            # carries into bits[1..L] (positions a+1 .. top_bit-1)
    sizes = RB.plan_blocks(L, p)
    blocks, pos = [], 1
    for s in sizes:
        blocks.append(list(range(pos, pos + s)))
        pos += s
    cw = [d0] + [None] * L
    for bi, blk in enumerate(blocks):
        last = bi == len(blocks) - 1
        for j in blk:
            w = T.alloc()
            T.ccx(bits[j - 1], cw[j - 1], w)    # carry into bits[j] = bits[j-1] AND carry into bits[j-1]
            for bq in zfix_up.get(a + j, ()):
                c.z(w, cond=bq)
            cw[j] = w
        if last and mP is not None:
            c.cz(bits[L], cw[L], cond=mP)        # carry out of the top window bit: a phase under mP
        top_w = L if last else blk[-1] - 1        # bits whose carry in and carry out are both known now
        for k in range(blk[0] - 1, top_w + 1):
            T.cx(cw[k], bits[k])
        for j in reversed(blk if last else blk[:-1]):
            T.erase(cw[j])                        # free: bits[j-1] before its write = bits[j-1] xor cw[j-1]
            cw[j] = None
    for bi in range(len(blocks) - 2, -1, -1):
        blk = blocks[bi]
        jb, inner = blk[-1], blk[:-1]
        if not inner:                            # its lower neighbour is the previous boundary or d0: alive
            T.erase(cw[jb])
            cw[jb] = None
            continue
        bw = cw[jb]
        m = c.hmr(bw)
        T.f.pop(bw)
        c.free(bw)
        cw[jb] = None
        c.push(m)                                 # under the outcome: rebuild the block below the boundary
        for j in inner:
            w = T.alloc()
            T.ccx(bits[j - 1], cw[j - 1], w)     # bits[j-1] is written (old xor cw[j-1]): (old xor c) c = old c xor c
            T.cx(cw[j - 1], w)                   # so the true carry old AND c needs c added back
            cw[j] = w
        # phase of the boundary: bits[jb-1] (before its write) AND cw[jb-1]; that bit holds it xor cw[jb-1] now
        c.cz(bits[jb - 1], cw[jb - 1])
        c.z(cw[jb - 1])
        for j in reversed(inner):
            T.erase(cw[j])
            cw[j] = None
        c.pop()
        c.bit_free(m)


def rows_rev_joint_ext(c, reg, specs, top_bit, mP, ghosts, qghost, zfix_up, p_chain, shed=0):
    """RB._rows_rev_joint with the top chain run by _chain_blocked when p_chain is an int, and (shed = 1) the
    widest row's carry-in erased after that row's writes (else identical)."""
    T = RB._Trk(c)
    D, cars, chain = [], [], []
    shed_mask = None
    for ri, sp in enumerate(specs):
        a = sp.m + 1
        win = reg[:a]
        for w in win:
            T.x(w)
            T.cx(sp.t, w)
        car = [sp.cin_make(T)]
        for j in range(a):
            RB._enter(T, win[j], sp.opb(j), car[j], False)
            w = T.alloc()
            RB._carry_into(T, win[j], sp.opb(j), car[j], w)
            RB._exit(T, win[j], sp.opb(j), car[j], False)
            if j + 1 < a and (sp.key, j + 1) in ghosts:
                c.z(w, cond=ghosts[(sp.key, j + 1)])
                c.bit_free(ghosts.pop((sp.key, j + 1)))
            car.append(w)
        for j in range(a):
            RB._write(T, win[j], sp.opb(j), car[j])
        for w in win:
            T.cx(sp.t, w)
            T.x(w)
        D.append(car[a])
        cars.append(car[:a] if sp.has_cin else car[1:a])
        if shed and ri == 0:
            # the carry-in is t' AND carrier (+ ca) over wires that never change here: erased for free now
            assert sp.has_cin
            shed_mask = T.f[car[0]]
            if _DBG:
                print(f"  [ainv] shed carry-in wire {car[0]} of row {sp.key}; carries {car[1:]}", flush=True)
            T.erase(car[0])
            cars[0][0] = None
    for i in range(len(specs) - 1, -1, -1):
        sp = specs[i]
        if i + 1 < len(specs):
            low = specs[i + 1]
            assert low.m + 2 == sp.m + 1, "joined rows must have consecutive shifts"
            bit = reg[low.m + 1]
            T.x(bit)
            T.cx(low.t, bit)
            T.ccx(bit, D[i + 1], D[i])
            T.cx(D[i + 1], bit)
            T.cx(low.t, bit)
            T.x(bit)
        if sp.key in qghost:
            c.z(D[i], cond=qghost[sp.key][0])
    sp = specs[0]
    a = sp.m + 1
    bits = reg[a:top_bit]
    if not bits:
        if mP is not None:
            c.z(D[0], cond=mP)
    elif p_chain is not None and p_chain < len(bits) - 1:
        for w in bits:
            T.x(w)
            T.cx(sp.t, w)
        _chain_blocked(c, T, bits, D[0], a, top_bit, mP, zfix_up, p_chain)
        for w in bits:
            T.cx(sp.t, w)
            T.x(w)
    else:
        for w in bits:
            T.x(w)
            T.cx(sp.t, w)
        cw = [D[0]]
        for idx, w in enumerate(bits):
            pos = a + idx + 1
            if pos < top_bit:
                nw = T.alloc()
                T.ccx(w, cw[idx], nw)
                for bq in zfix_up.get(pos, ()):
                    c.z(nw, cond=bq)
                cw.append(nw)
                chain.append(nw)
            elif mP is not None:
                c.cz(w, cw[idx], cond=mP)
        for idx, w in enumerate(bits):
            T.cx(cw[idx], w)
        for w in bits:
            T.cx(sp.t, w)
            T.x(w)
    er = (lambda w_: _erase_smart(T, w_)) if shed else T.erase
    for w in reversed(chain):
        er(w)
    for i in range(len(specs)):
        er(D[i])
        lst = cars[i]
        for idx in range(len(lst) - 1, -1, -1):
            w = lst[idx]
            if w is None:
                continue
            if i == 0 and shed and idx == 1:
                # the carry above the shed carry-in: measure it; under the outcome re-make the carry-in (same
                # function, same tracked content), apply this carry's phase, erase the carry-in again
                mask = T.f.pop(w)
                m = c.hmr(w)
                c.free(w)
                c.push(m)
                w2 = specs[0].cin_make(c)
                T.f[w2] = shed_mask
                _phase_of(T, mask)
                _erase_smart(T, w2)
                c.pop()
                c.bit_free(m)
            else:
                er(w)


def group_rev_ext(c, reg, specs, top_bit, mP, ghosts, qghost, pool, tail=True):
    """RB._group_rev(..., top=None, tree=True, join=True) with plan_join_ext."""
    keys = [sp.key for sp in specs]

    def zfix_for(idx):
        zfix = {}
        for k2 in keys[:idx + 1]:
            if k2 in qghost:
                zfix.setdefault(qghost[k2][1], []).append(qghost[k2][0])
        return zfix

    for i, j, p, s in plan_join_ext(specs, top_bit, pool):
        if i == j:
            RB._row_rev_spec(c, reg[:top_bit], None, specs[i], ghosts, mP, pool, zfix_for(i), tail)
        else:
            rows_rev_joint_ext(c, reg, specs[i:j + 1], top_bit, mP, ghosts, qghost, zfix_for(i - 1), p, s)
    for sp in specs:
        if sp.key in qghost:
            bq, _ = qghost.pop(sp.key)
            c.z(sp.t, cond=bq)
            c.bit_free(bq)


def emit_A_inv_ext(c, t, win, c2, ctx, mP, pool=10, tail=True):
    """Drop-in for rows2p1.emit_A_inv (same context, same repairs) with the extended planner."""
    nt = len(t)
    specs = [RB._spec_A(c2, k, t[nt - 1 - k], ctx["kg"]) for k in range(ctx["nrows"])]
    group_rev_ext(c, win, specs, len(win), mP, ctx["ghosts"], ctx["qghost"], pool, tail)
    assert not ctx["ghosts"] and not ctx["qghost"], "unrepaired forward outcomes"
