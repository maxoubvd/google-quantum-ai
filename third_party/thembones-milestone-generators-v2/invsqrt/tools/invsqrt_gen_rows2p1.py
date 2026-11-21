"""invsqrt_rows2_b: second rows module for stage B (lo -+= u >> k) of our invsqrt generator: "scheduled groups".

Own code of workstream invsqrt_rows2, written against tools/invsqrt_gen_circ.py (builder + exact simulator).
Row primitives (temporary-AND rows with unrepaired forward carries, tree nodes, reverse rows in blocks, tail
variant) are imported read-only from the first rows module (tools/invsqrt_rows_b.py, workstream invsqrt_rows;
method credits there: manwithacat "HW-25", gengsr, maoxin123; Gidney 1709.06648; Litinski 2410.00899).

What is new here (research/ws_invsqrt_rows2.md):
  - the pending wire of a tree group is made only when the group's widest row runs (`lazy`): until then the
    group holds a temporary node and two carries, so the forward pass peaks at (pending wires + 2) clean wires
    instead of (pending wires + 3), in every group;
  - merges with unrepaired carries (`dirty_merge`): two carry wires; the un-merge recomputes the same carries
    and repairs them;
  - un-merge with a wire limit (`pool_un`): every pending wire is measured as soon as its own un-merge is
    done; the carry of the next un-merge into it is then a phase (no Toffoli), and a chain that is short of
    wires runs in blocks;
  - with the forward peak gone, PH2 is what limits the number of pending wires (Q = 56 + pending), and one
    more group than the first module fits the same circuit (partitions for 58 to 61 wires in the log);
  - a group may run as flat rows (`flat`: two carry wires, no node). The order of the groups is NOT free: a
    first variant that ran a wide group after the rows of t4..t0 failed the exact test (order rule in the
    log); `order` is kept as an argument but only the natural order is accepted.
Arithmetic: exactly rows_model of the first module (= the generator's row_virt rows). No change in the table fit.
"""
from __future__ import annotations

import importlib

RB = None                  # the first rows module (row primitives); see bind()
rows_model = rows_model_R21 = rows_model_A = None
PUBLIC_GROUPS = ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11, 12, 13))


def bind(base="invsqrt_rows_b"):
    """Choose the module that provides the row primitives: the first rows module as it is now (it must have
    RowSpec(..., cfree=...), i.e. the file of 12:35Z, sha256 4a063c2c...) or a pinned verbatim copy of it
    (the generator calls bind("<name of its pinned copy>") once after the import)."""
    global RB, rows_model, rows_model_R21, rows_model_A
    RB = importlib.import_module(base)
    rows_model = RB.rows_model
    rows_model_R21 = RB.rows_model_R21
    rows_model_A = RB.rows_model_A
    return RB


bind()


def groups_from_sizes(sizes):
    out, k = [], 0
    for s in sizes:
        out.append(tuple(range(k, k + s)))
        k += s
    return tuple(out)


# ------------------------------------------------------------------ merge of a pending wire (one-bit increment)

def _incr_fwd_dirty(c, bits, x, top, gh):
    """bits += x at bits[0], carry out of bits[-1] XORed into top. Two carry wires; the carry into bits[j]
    (j >= 1) is measured without repair, outcome in gh[j]."""
    car = x
    n = len(bits)
    for j in range(n):
        if j == n - 1:
            c.ccx(bits[j], car, top)
            w = None
        else:
            w = c.alloc()
            c.ccx(bits[j], car, w)
        c.cx(car, bits[j])
        if j >= 1:
            gh[j] = c.hmr(car)
            c.free(car)
        car = w


def _incr_clean(c, bits, x, top, gh=None):
    """Same increment with every carry alive (len(bits) - 1 wires) and erased for free; gh: outcomes of a
    forward _incr_fwd_dirty on the same values, repaired here (Z on the recomputed carry)."""
    car = [x]
    n = len(bits)
    for j in range(n):
        if j == n - 1:
            c.ccx(bits[j], car[j], top)
        else:
            w = c.alloc()
            c.ccx(bits[j], car[j], w)
            if gh is not None and (j + 1) in gh:
                c.z(w, cond=gh[j + 1])
                c.bit_free(gh.pop(j + 1))
            car.append(w)
    for j in range(n - 1, -1, -1):
        if j < n - 1:
            m = c.hmr(car[j + 1])
            c.cz(bits[j], car[j], cond=m)
            c.bit_free(m)
            c.free(car[j + 1])
        c.cx(car[j], bits[j])


def _merge(c, P, tf, bits, top, inverse, gh, dirty):
    """Forward: bits +- (pending digit) in the frame of the group's first digit tf, carry into top; inverse: the
    same adder on the complemented frame (same carries)."""
    c.cx(tf, P)
    for a in bits:
        c.cx(tf, a)
        if inverse:
            c.x(a)
    if inverse:
        _incr_clean(c, bits, P, top, gh)
    elif dirty:
        _incr_fwd_dirty(c, bits, P, top, gh)
    else:
        _incr_clean(c, bits, P, top, None)
    for a in bits:
        if inverse:
            c.x(a)
        c.cx(tf, a)
    c.cx(tf, P)


# ------------------------------------------------------------------ forward tree group, pending wire made late

def _group_fwd_tree(c, reg, specs, top_bit, ghosts, qghost, top_wire=None, make_top=None):
    """Tree group as in the first module (rows narrowest first, row i on reg[:m_i + 1] with its own node on top,
    the node of the row below merged one bit up and measured), with one change: the group's pending wire is
    only made (make_top(), returns the wire) when the row that carries into it is about to run, i.e. for the
    widest row. Until then the group holds a temporary node and two carries, not the pending wire.
    Returns the pending wire (or top_wire)."""
    prev = None
    made = top_wire
    for sp in reversed(specs):
        a = sp.m + 1
        assert a <= top_bit
        if a == top_bit:
            if made is None:
                made = make_top()
            node, temp = made, False
        else:
            node, temp = c.alloc(), True
            c.cx(sp.t, node)
        if prev is not None:
            pw, pa, psp = prev
            RB._merge(c, pw, psp.t, reg[pa:a], node, False)
            qghost[psp.key] = (c.hmr(pw), pa)
            c.free(pw)
        RB._row_fwd(c, reg[:a], node, sp, ghosts)
        prev = (node, a, sp) if temp else None
    if prev is not None:
        pw, pa, psp = prev
        if made is None:
            made = make_top()
        RB._merge(c, pw, psp.t, reg[pa:top_bit], made, False)
        qghost[psp.key] = (c.hmr(pw), pa)
        c.free(pw)
    return made


# ------------------------------------------------------------------ stage B

G59 = groups_from_sizes((2, 4, 4, 4))          # 3 pending wires: 59-wire circuit (pools 7 / 5, pool_un 5)
G60 = groups_from_sizes((2, 3, 3, 2, 4))       # 4 pending wires: 60-wire circuit (pools 8 / 6, pool_un 6)


def emit_B(c, t, u, lo, groups=PUBLIC_GROUPS, order=None, flat=(), hook=None, dirty_merge=True, lazy=True,
           tree=None):
    """Forward rows and merges. groups: partition of the shifts 0..len(t)-1 into consecutive runs, widest first
    (group 0 carries into lo[-1], every other group has a pending wire). order: group indices in the order in
    time in which they run (default 0, 1, 2, ...). flat: group indices that run as flat rows (widest first, two
    carry wires, no node); the others run as tree groups (narrowest first, one node and two carries).
    hook(event, group) at 'group_start', 'group_end', ('merge', None). Returns the context for emit_B_inv;
    ctx['P'] are the pending wires, alive until emit_B_inv measures them. `tree` is accepted and ignored (the
    first module's argument: groups are tree groups unless listed in `flat`)."""
    nt, nu, wlo = len(t), len(u), len(lo)
    assert sorted(k for g in groups for k in g) == list(range(nt))
    lay = RB._layout(groups, nu, wlo)
    order = list(range(len(lay))) if order is None else list(order)
    assert sorted(order) == list(range(len(lay)))
    # Order rule (research/ws_invsqrt_rows2.md, failed variant v0): between the first carry a pending wire
    # receives and its merge no row of a group above it may run. With the merges at the end this leaves the
    # natural order only; another order needs the lower pending wires merged first (not built).
    assert order == sorted(order), "groups must run widest first: see the order rule"
    flat = set(flat)
    ghosts, pend, qghost, specs, mgh = {}, {}, {}, {}, {}
    for gi in order:
        grp, w, has_p = lay[gi]
        if hook:
            hook("group_start", gi)
        specs[gi] = [RB._spec_B(u, k, t[nt - 1 - k]) for k in grp]
        tree = gi not in flat
        if has_p:
            tf = t[nt - 1 - grp[0]]

            def make_top(tf=tf):
                P = c.alloc()
                c.cx(tf, P)
                return P

            if tree and lazy:
                P = _group_fwd_tree(c, lo, specs[gi], w, ghosts, qghost, None, make_top)
            else:
                P = make_top()
                RB._group_fwd(c, lo, specs[gi], P, w, ghosts, qghost, tree)
            pend[gi] = (P, tf, w)
        else:
            RB._group_fwd(c, lo, specs[gi], lo[w], w, ghosts, qghost, tree)
        if hook:
            hook("group_end", gi)
    if hook:
        hook("merge", None)
    for gi in range(len(lay) - 1, 0, -1):
        P, tf, w = pend[gi]
        mgh[gi] = {}
        if gi == 1:
            _merge(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], False, mgh[gi], dirty_merge)
        else:
            Pup, _, wup = pend[gi - 1]
            _merge(c, P, tf, lo[w:wup], Pup, False, mgh[gi], dirty_merge)
    return dict(ghosts=ghosts, pend=pend, lay=lay, qghost=qghost, specs=specs, order=order, flat=flat, mgh=mgh,
                P=[pend[g][0] for g in sorted(pend)])


def _unmerge_pool(c, P, tf, bits, top, m_up, gh, pool, key):
    """Inverse merge with at most `pool` carry wires: the increment chain is run as a reverse row without
    operand whose carry-in is the pending wire (first module's block machinery: a boundary carry is measured
    and repaired under its outcome). top: wire that takes the carry out of bits[-1] back, or None when that
    carry went into a pending wire that is already measured (outcome m_up: repaired by the carry's phase)."""
    sp = RB.RowSpec(key, [], 0, tf, None, None, cfree=P)
    ghosts = {(key, j): b for j, b in gh.items()}
    gh.clear()
    c.cx(tf, P)
    RB._row_rev_core(c, bits, top, sp, ghosts, m_up, pool, None, False, True)
    c.cx(tf, P)
    assert not ghosts


def emit_B_inv(c, t, u, lo, ctx, pool=10, hook=None, tail=True, pool_un=None, join=False):
    """Exact inverse of emit_B: un-merge (clean, repairs the merge carries), pending wires measured, then the
    groups in the mirror of their forward order as flat reverse rows. pool: clean wires a group's rows may hold
    (int, dict {group: n, 'default': n} or function). hook at ('unmerged', None), 'group_start', 'group_end'.
    pool_un: clean wires this module may hold during the un-merge, pending wires included (None: no limit,
    3 pending + 3 carries with the public partition). With a limit every pending wire is measured as soon as
    its own un-merge is done (the carry of the next un-merge into it is then repaired by a phase), and an
    un-merge that is short of wires runs in blocks. `join` is accepted and ignored (first module's argument)."""
    wlo = len(lo)
    lay, pend, ghosts, qghost = ctx["lay"], ctx["pend"], ctx["ghosts"], ctx["qghost"]
    mP = {}
    if pool_un is None:
        for gi in range(1, len(lay)):
            P, tf, w = pend[gi]
            if gi == 1:
                _merge(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], True, ctx["mgh"][gi], False)
            else:
                Pup, _, wup = pend[gi - 1]
                _merge(c, P, tf, lo[w:wup], Pup, True, ctx["mgh"][gi], False)
            assert not ctx["mgh"][gi]
        for gi in range(1, len(lay)):
            P, tf, w = pend[gi]
            mP[gi] = c.hmr(P)
            c.z(tf, cond=mP[gi])
            c.free(P)
    else:
        alive = len(lay) - 1
        for gi in range(1, len(lay)):
            P, tf, w = pend[gi]
            p = pool_un - alive
            assert p >= 1, "pool_un too small"
            if gi == 1:
                _unmerge_pool(c, P, tf, lo[w:wlo - 1], lo[wlo - 1], None, ctx["mgh"][gi], p, ("M", gi))
            else:
                _, _, wup = pend[gi - 1]
                _unmerge_pool(c, P, tf, lo[w:wup], None, mP[gi - 1], ctx["mgh"][gi], p, ("M", gi))
            mP[gi] = c.hmr(P)
            c.z(tf, cond=mP[gi])
            c.free(P)
            alive -= 1
    if hook:
        hook("unmerged", None)
    for gi in reversed(ctx["order"]):
        grp, w, has_p = lay[gi]
        if hook:
            hook("group_start", gi)
        RB._group_rev(c, lo, ctx["specs"][gi], w, None if has_p else lo[w], mP.get(gi), ghosts, qghost,
                      RB._pool_of(pool, gi), gi not in ctx["flat"], tail, False)
        if has_p:
            c.bit_free(mP[gi])
        if hook:
            hook("group_end", gi)
    assert not ghosts and not qghost, "unrepaired forward outcomes"


# ------------------------------------------------------------------ stage A (rows on the low window of u)

def emit_A(c, t, win, pend, c2, nrows=8, kg=7):
    """Stage A rows as a tree group: (win + pending wire as top bit) -+= (c2 >> (k + 1)) + cin, k = 0..nrows-1,
    the first module's arithmetic (rows_model_A) and row circuits. One change: `pend` may be a FUNCTION that
    returns the pending wire; it is called when the rows are done, just before the last node is merged into
    the top window bit. Until then no row touches win[-1] (every row window ends below it), so the caller can
    still load the pending wire's initial content from that wire there (the generator's P0 = 1 - win[-1]).
    The stage then holds 3 clean wires (node + 2 carries) instead of 4. Returns the context; ctx['pend'] is
    the pending wire."""
    nt = len(t)
    assert max(RB._spec_A(c2, k, t[nt - 1 - k], kg).m + 1 for k in range(nrows)) < len(win), \
        "a row window reaches the top window bit: the pending wire cannot be made late"
    specs = [RB._spec_A(c2, k, t[nt - 1 - k], kg) for k in range(nrows)]
    ghosts, qghost = {}, {}
    if callable(pend):
        P = _group_fwd_tree(c, win, specs, len(win), ghosts, qghost, None, pend)
    else:
        P = _group_fwd_tree(c, win, specs, len(win), ghosts, qghost, pend, None)
    return dict(ghosts=ghosts, qghost=qghost, tree=True, pend=P, nrows=nrows, kg=kg)


def emit_A_inv(c, t, win, c2, ctx, mP, pool=10, tail=True, join=True):
    """Exact inverse of the rows of emit_A on the window wires (first module's reverse rows, joined as far as
    `pool` clean wires allow). mP: outcome bit of the pending wire's X measurement; the caller repairs the
    pending wire's initial content. The rows are described with the c2 wires passed HERE (c2 may sit on other
    wires than in emit_A)."""
    nt = len(t)
    specs = [RB._spec_A(c2, k, t[nt - 1 - k], ctx["kg"]) for k in range(ctx["nrows"])]
    RB._group_rev(c, win, specs, len(win), None, mP, ctx["ghosts"], ctx["qghost"], pool, True, tail, join)
    assert not ctx["ghosts"] and not ctx["qghost"], "unrepaired forward outcomes"
