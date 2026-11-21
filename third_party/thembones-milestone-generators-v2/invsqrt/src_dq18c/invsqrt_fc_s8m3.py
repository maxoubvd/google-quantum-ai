#!/usr/bin/env python3
"""invsqrt_fc_s8m3.py -- FC (eighth holder, 2026-11-16): SEVERAL WORDS IN ONE FORWARD PASS (R45-INV-S7B ledger 9).

A run of consecutive narrow adds into the same target R (tgt += src_k + cin_k, k = 1..m, only X / CX ops between
them, which act on T's wires and conjugate R by X^n) is replaced by ONE merged forward block:
  - position by position, the chains are processed in their separate order; chain k's cell runs in chain k's own
    T state (the ops between the blocks, applied and undone around the cell; CNOTs are free), so every literal and
    borrow-in is the one the separate pass used;
  - two NEIGHBOURING chains (in that order) that both have a bare carry at a position (no literal there) are merged:
    MAJ(R[p], x, y) is one AND for both. Exact, and the merged carry is the XOR of the separate passes' carries
    there (two carries out of one position with no literal in between are never both 1); a merged group stays
    merged up to the top;
  - at the position below the top every carry goes into the top bit by CCX (mod 2^n).
The UNDO is unchanged: the separate passes' undo blocks in reverse order, with the mirrored ops between them. A
merged carry's X-measurement bit is repaid by Z on every separate carry of its group (the undo of each absorbed pass
reads the owner's bit there).
HELD CARRY ("hold", the six-bit word): the word's pass runs only up to its last literal and its carry out is HELD on
a wire (word_hold); the word block, the ops after it and the rows' run become ONE block whose forward absorbs the held
carry at its position as the first chain (the word comes first in the separate order; the rows are then generated
in the order 8.1, 8.2, 9 so that the word's chain is next to 8.1's); its undo ends with the word's own undo, which
re-creates the full pass's carries (their debts: the merged forward's bits).
Wires: the carries use wires of the free list at the run (the free list is not changed); the held wire is not
counted live after the word's pass (the undo blocks run where it does not exist) and goes back to the free list.

Hook (src copy of invsqrt_fc2_patch_h1c.py, cfg h1c "m3": 1, "m3rev": 1, h1c.sw "hold": 1): st = begin(c) before the
rows loop; merge(c, st) after it; word_hold in the word block's forward.
"""
LOG = {}            # bid -> (src, tgt, cin)
AUX = {}            # R48-INV-S7C: bid -> (pre ops, post ops, wire): a chain's auxiliary carry-in wire
HOLD = {}           # the held carry of the word block: wire, rpos, bits {R position: bit}, bid
PH = {}             # "c0m": c0's wire, the phantom top id P, the bit m0 of c0's early measurement in the mirror
STATE = {"inst": False}
LIN = ("X", "CX")


# ---- "c0m": c0 measured right after its last read in the mirror (the rule's undo); the H1 undo generators run with
# its wire as one more clean wire and write the top bit into a PHANTOM id, whose writes become phases under m0
def c0m_begin(c, c0):
    PH.clear()
    PH.update(c0=c0, P=10 ** 6, m0=c.bit())


def ph_target(R):
    """an undo's target list with c0 replaced by the phantom."""
    return [PH["P"] if (PH and w == PH["c0"]) else w for w in R]


def ph_real(R):
    """the target wires without c0 (for the X / CX frames outside the blocks: the top is write-only in a pass)."""
    return [w for w in R if not (PH and w == PH["c0"])]


class ph_free:
    """c0's wire is free for the allocator while an undo is generated (the undo runs after c0's measurement)."""

    def __init__(self, c):
        self.c = c

    def __enter__(self):
        if PH:
            c = self.c
            c.free_([PH["c0"]])
            c.free.remove(PH["c0"])
            c.free.insert(0, PH["c0"])
        return self

    def __exit__(self, *a):
        if PH:
            c = self.c
            assert PH["c0"] in c.free, "c0's wire is still held at the end of an undo"
            c.free.remove(PH["c0"])
            c.live += 1
            c.peak = max(c.peak, c.live)
        return False


def ph_fix(ops):
    """writes into the phantom -> their phases under m0 (CCX a b P -> CZ a b, CX x P -> Z x, X P -> NEG)."""
    if not PH:
        return ops
    P, m = PH["P"], PH["m0"]
    out = []
    for o in ops:
        g = o[0]
        if g == "CCX" and o[3] == P:
            assert P not in o[1:3]
            out.append(("CZ", o[1], o[2], m))
        elif g == "CX" and o[2] == P:
            assert o[3] is None, ("a conditioned write into the phantom", o)
            out.append(("Z", o[1], m))
        elif g == "X" and o[1] == P:
            assert o[2] is None, ("a conditioned write into the phantom", o)
            out.append(("NEG", m))
        else:
            assert P not in o[1:], ("the phantom is read or written otherwise", o)
            out.append(o)
    return out


def install():
    import invsqrt_fc_narrow as NW
    if STATE["inst"]:
        return
    prev = NW.block

    def block(c, src, tgt, cin):
        nb0 = len(c.blocks) if hasattr(c, "blocks") else 0
        r = prev(c, src, tgt, cin)
        if hasattr(c, "blocks") and len(c.blocks) == nb0 + 1:
            LOG[nb0] = (list(src), list(tgt), cin)
        return r
    NW.block = block
    STATE["inst"] = True


def begin(c):
    install()
    return {"start": len(c.ops)}


class _LazyTop(list):
    """a target list whose last wire is allocated when it is first read (the held carry, at the top step)."""

    def __init__(self, items, c, box):
        super().__init__(list(items) + [None])
        self._c, self._box = c, box

    def __getitem__(self, i):
        if isinstance(i, int) and (i == len(self) - 1 or i == -1) and self._box.get("w") is None:
            self._box["w"] = self._c.alloc(1)[0]
            list.__setitem__(self, len(self) - 1, self._box["w"])
        return list.__getitem__(self, i)


def word_hold(c, NW, fwd_stream, lits_, Rt, cw, srcw, bid):
    """the word's forward pass up to its last literal; the carry out of that position is HELD (the top step's CCX
    writes it into a fresh wire); returns the debts of the FULL pass (the carries above the last literal get bits
    that the rows' merged forward will write) for the word's undo generator."""
    last = max(i for i, w in enumerate(srcw) if w is not None)
    box = {"w": None}
    tgt = _LazyTop(Rt[:last + 1], c, box)
    lt = list(lits_[:last + 1]) + [((lambda: None), (lambda: None))]
    d = fwd_stream(c, NW, lt, tgt, NW._lit(cw))
    h = box["w"]
    assert h is not None
    c.free_([h])                                  # free for the allocator from here (the undo blocks run where it does
    c.free.remove(h)                              # not exist), at the BOTTOM of the free list so that the forward's
    c.free.insert(0, h)                           # allocations before the rows' pass take other wires
    full = dict(d)
    bits = {}
    for t in range(last + 1, len(Rt) - 1):        # the full pass's carries into target index t (R position t + 1)
        full[t] = bits[t + 1] = c.bit()
    HOLD.clear()
    HOLD.update(wire=h, rpos=last + 2, bits=bits, bid=bid)
    install()
    return full


def _top_blocks(ops, start):
    out = []
    i = start
    while i < len(ops):
        if ops[i][0] == "BBEG":
            bid = ops[i][1]
            k = i
            while not (ops[k][0] == "BEND" and ops[k][1] == bid):
                k += 1
            out.append((i, k, bid))
            i = k + 1
            continue
        i += 1
    return out


def _lin(seg):
    return all(o[0] in LIN and o[-1] is None for o in seg)


def _runs(ops, blks):
    runs, cur = [], blks[:1]
    for a, b in zip(blks[:-1], blks[1:]):
        ok = (a[2] in LOG and b[2] in LOG and _suffix(LOG[a[2]][1], LOG[b[2]][1]) and _lin(ops[a[1] + 1:b[0]]))
        if ok:
            cur.append(b)
        else:
            if len(cur) > 1:
                runs.append(cur)
            cur = [b]
    if len(cur) > 1:
        runs.append(cur)
    return runs


def _suffix(t1, t2):
    """R48-INV-S7C: one target is a suffix of the other (a chain that starts at a higher position)."""
    a, b = (t1, t2) if len(t1) <= len(t2) else (t2, t1)
    return len(a) > 0 and list(b[len(b) - len(a):]) == list(a)


def _inv(seg):
    return [o for o in reversed(seg)]          # X / CX only: self-inverse ops in reverse order


def _remap(ops, mp):
    """conditions (classical bits) renamed by mp in a list of ops."""
    out = []
    for o in ops:
        g = o[0]
        if g in ("X", "CX", "Z", "CZ") and o[-1] in mp:
            o = o[:-1] + (mp[o[-1]],)
        elif g == "NEG" and o[1] in mp:
            o = (g, mp[o[1]])
        elif g == "ZP" and o[3] in mp:
            o = o[:3] + (mp[o[3]],)
        elif g in ("CANDC", "CANDE") and o[4] in mp:
            o = o[:4] + (mp[o[4]],)
        elif g in ("CBEG", "CEND") and o[1] in mp:
            o = (g, mp[o[1]])
        out.append(o)
    return out


def merge(c, st, verbose=True):
    import invsqrt_fc_narrow as NW
    import invsqrt_fc_kit as K
    ops = c.ops
    blks = _top_blocks(ops, st["start"])
    runs = _runs(ops, blks)
    if not runs:
        if verbose:
            print("m3: no run of adds into one target")
        return 0.0
    saved = 0.0
    hold = dict(HOLD) if HOLD else None
    for ri, run in enumerate(reversed(runs)):    # from the back: earlier op indices stay valid
        wb = None
        if hold and ri == len(runs) - 1:          # the first run takes the held carry: the word block before it
            i0w = next(i for i, o in enumerate(ops) if o[0] == "BBEG" and o[1] == hold["bid"])
            i1w = next(i for i in range(i0w, len(ops)) if ops[i][0] == "BEND" and ops[i][1] == hold["bid"])
            assert _lin(ops[i1w + 1:run[0][0]]), "m3: ops between the word block and the rows are not linear"
            wb = (i0w, i1w, hold)
        saved += _merge_run(c, K, NW, run, verbose, wb)
        ops = c.ops
    if hold:
        HOLD.clear()
    return saved


def _merge_run(c, K, NW, run, verbose, wb=None):
    ops = c.ops
    tgt = max((LOG[r_[2]][1] for r_ in run), key=len)   # R48-INV-S7C: the longest target (the others are suffixes)
    if PH:                                         # c0m: the separate passes' targets carry the phantom top; the
        tgt = [PH["c0"] if w == PH["P"] else w for w in tgt]   # merged forward writes the real c0
    n = len(tgt)
    srcs, cins, bits, offs, auxs = [], [], [], [], []
    hw = wb is not None
    if hw:                                         # chain 0: the held carry (no literal, no carry-in)
        srcs.append([None] * n)
        cins.append(None)
        bits.append([])
        offs.append(0)
        auxs.append(None)
    for (i0, i1, bid) in run:
        s_, t_, ci = LOG[bid]
        of_ = n - len(t_)                          # R48-INV-S7C: this chain starts at position of_
        srcs.append([None] * of_ + list(s_) + [None] * (n - of_ - len(s_)))
        cins.append(ci)
        offs.append(of_)
        auxs.append(AUX.get(bid))
        bits.append([o[2] for o in ops[i0 + 1:i1] if o[0] == "HMR"])
    M = len(srcs)
    off = 1 if hw else 0
    gaps = [ops[run[k][1] + 1:run[k + 1][0]] for k in range(len(run) - 1)]
    gT = [[o for o in g if not (o[0] == "X" and o[1] in tgt)] for g in gaps]
    for g in gaps:
        nx = sum(1 for o in g if o[0] == "X" and o[1] in tgt)
        assert nx in (2 * n, 2 * (n - 1)), ("m3: the gap is not a pair of X^n on R", nx)
        assert all(not (o[0] == "CX" and (o[1] in tgt or o[2] in tgt)) for o in g), "m3: a gap op touches R"
    to_state = [[] for _ in range(off)]            # ops that take the first row block's state to chain k's
    acc = []
    for k in range(len(run)):
        to_state.append(list(acc))
        if k < len(run) - 1:
            acc = acc + gT[k]
    # the separate forwards' debts per chain: carry index (R position) -> bit
    debts = []
    if hw:
        debts.append(dict(wb[2]["bits"]))
    for k in range(off, M):
        d = {}
        cur = NW._lit(cins[k])
        idx = 0
        for i in range(offs[k], n - 1):            # R48-INV-S7C: from the chain's first position
            f = NW._form_b(srcs[k][i], tgt[i], cur)
            if cur[0] == "w" and cur[1] == -1:
                d[i] = bits[k][idx]
                idx += 1
            if i < n - 2:
                if f[0] == "and":
                    cur = ("w", -1, 0)
                elif f[0] == "const":
                    cur = ("k", f[1])
                else:
                    cur = ("w", -2, f[2])
        assert idx == len(bits[k]), ("m3: debt count", k, idx, len(bits[k]))
        debts.append(d)
    # ---- the merged forward
    out = []
    _auxw = [a_[2] for a_ in auxs if a_ is not None]              # R48-INV-S7C: auxiliary wires are live too
    pool = [w for w in c.free if not (hw and w == wb[2]["wire"]) and w not in _auxw]
    live0 = c.live + (1 if hw else 0) + len(_auxw)
    if hw:                                         # the forward between the word block and the rows must not use it
        assert all(wb[2]["wire"] not in o[1:] for o in ops[wb[1] + 1:run[0][0]]), "m3: the held wire is touched"
    nlive = [0, 0]

    def alloc():
        w = pool.pop()
        nlive[0] += 1
        nlive[1] = max(nlive[1], nlive[0])
        return w

    def free(w):
        pool.append(w)
        nlive[0] -= 1

    save_ops = c.ops
    c.ops = out

    def in_state(k, fn):
        for o in to_state[k]:
            c.ops.append(o)
        fn()
        for o in _inv(to_state[k]):
            c.ops.append(o)

    carry = [NW._lit(cins[k]) if offs[k] == 0 else ("k", 0) for k in range(M)]  # chain k's carry into the position
    fresh = [False] * M
    cidx = [0] * M
    grp = [{k} for k in range(M)]                 # the chains whose separate carries a live carry stands for
    rename = {}

    def retire(k):
        """X-measure chain k's incoming carry with its own bit; the group's other chains read that bit there."""
        if not fresh[k]:
            return
        b = debts[k][cidx[k]]
        c.ops.append(("HMR", carry[k][1], b))
        for j in grp[k] - {k}:
            rename[debts[j][cidx[k]]] = b
        if carry[k][1] != (wb[2]["wire"] if hw else None):
            free(carry[k][1])
        fresh[k] = False

    for k_ in range(M):                           # R48-INV-S7C: auxiliary wires, in their chain's T state
        if auxs[k_] is not None:
            in_state(k_, lambda k_=k_: c.ops.extend(auxs[k_][0]))
    for p in range(n - 1):
        top = (p == n - 2)
        for k_ in range(M):                       # R48-INV-S7C: a late chain's carry-in enters at its first position
            if offs[k_] > 0 and p == offs[k_]:
                carry[k_] = NW._lit(cins[k_])
        if hw and p == wb[2]["rpos"]:
            carry[0] = ("w", wb[2]["wire"], 0)
            fresh[0] = True
            cidx[0] = p
        order = [k for k in range(M) if srcs[k][p] is not None or carry[k] != ("k", 0)]
        i = 0
        while i < len(order):
            k = order[i]
            bare = srcs[k][p] is None and carry[k][0] == "w"
            k2 = order[i + 1] if i + 1 < len(order) else None
            if bare and k2 is not None and srcs[k2][p] is None and carry[k2][0] == "w":
                # MAJ(R[p], x, y): one AND, one carry for both groups
                x, y, r = carry[k][1], carry[k2][1], tgt[p]
                assert carry[k][2] == 0 and carry[k2][2] == 0
                c.cx(x, r)
                c.cx(x, y)
                t = None
                if top:
                    c.ccx(r, y, tgt[n - 1])
                    c.cx(x, tgt[n - 1])
                else:
                    t = alloc()
                    c.and_compute(r, y, t)
                    c.cx(x, t)                    # t = MAJ(r, x, y)
                c.cx(x, y)
                c.cx(y, r)                        # R[p] = r ^ x ^ y
                retire(k)
                retire(k2)
                if not top:
                    carry[k] = ("w", t, 0)
                    fresh[k] = True
                    cidx[k] = p + 1
                    grp[k] |= grp[k2]
                    grp[k2] = set()
                    carry[k2] = ("k", 0)
                    cidx[k2] = p + 1
                i += 2
                continue
            # own cell (a literal here, a constant carry, or a lone bare carry), in chain k's T state
            assert grp[k], ("m3: an absorbed chain has a literal later", k, p)
            a, bw, cur = srcs[k][p], tgt[p], carry[k]
            res = {}

            def cell(a=a, bw=bw, cur=cur, res=res):
                f = NW._form_b(a, bw, cur)
                if f[0] == "const":
                    if top and f[1]:
                        c.x(tgt[n - 1])
                    res["nxt"] = ("k", f[1])
                elif f[0] == "alias":
                    raise AssertionError("m3: alias carries are not handled")
                else:
                    _, xw, yw, ops_ = f
                    NW._emit(c, ops_)
                    if top:
                        c.ccx(xw, yw, tgt[n - 1])
                    else:
                        res["t"] = alloc()
                        c.and_compute(xw, yw, res["t"])
                    NW._emit(c, ops_[::-1])
                    NW._xor_lit(c, cur, tgt[n - 1] if top else res["t"])
                    res["nxt"] = ("w", res.get("t"), 0)
                NW._xor_src(c, a, bw)             # sum bit p
                NW._xor_lit(c, cur, bw)
            in_state(k, cell)
            retire(k)
            if not top:
                carry[k] = res["nxt"]
                fresh[k] = res["nxt"][0] == "w" and "t" in res
                cidx[k] = p + 1
            i += 1
    for k in range(M):                            # literals of the top position (the separate passes' last step)
        if srcs[k][n - 1] is not None:
            in_state(k, lambda k=k: NW._xor_src(c, srcs[k][n - 1], tgt[n - 1]))
    assert all(not f for f in fresh), "m3: a carry is left"
    for k_ in range(M):                           # R48-INV-S7C: auxiliary wires erased
        if auxs[k_] is not None:
            in_state(k_, lambda k_=k_: c.ops.extend(auxs[k_][1]))
    for o in to_state[M - 1]:                     # end in the last pass's T state (its BEND)
        c.ops.append(o)
    c.ops = save_ops
    written = set(o[2] for o in out if o[0] == "HMR")
    allb = set(b for bl in bits for b in bl) | set(debts[0].values() if hw else [])
    assert written | set(rename) == allb and not (written & set(rename)), (
        "m3: bits", sorted(allb - written - set(rename)), sorted(written & set(rename)))
    # ---- the undo: separate undo blocks in reverse order, mirrored gaps between, renamed bits
    und = []
    for k in range(len(run) - 1, -1, -1):
        ax_ = AUX.get(run[k][2])                  # R48-INV-S7C: the chain's auxiliary wire around its undo block
        if ax_ is not None:
            und.extend(ax_[3] if len(ax_) > 3 else ax_[0])
        und.extend(ph_fix(_remap(c.blocks[run[k][2]], rename)))
        if ax_ is not None:
            und.extend(ax_[4] if len(ax_) > 3 else ax_[1])
        if k > 0:
            und.extend(K.Circ.mirror(gaps[k - 1]))
    i0, i1 = run[0][0], run[-1][1]
    fw = list(out)
    if hw:
        i0w, i1w, _h = wb
        g0 = ops[i1w + 1:i0]
        und.extend(K.Circ.mirror(g0))
        und.extend(ph_fix(_remap(c.blocks[_h["bid"]], rename)))
        fw = list(ops[i0w + 1:i1w]) + list(g0) + fw
        i0 = i0w
    bid = len(c.blocks)
    c.blocks[bid] = und
    before = c.toffolis_exp(ops[i0:i1 + 1])
    c.ops = ops[:i0] + [("BBEG", bid)] + fw + [("BEND", bid)] + ops[i1 + 1:]
    after = c.toffolis_exp(fw)
    c.peak = max(c.peak, live0 + nlive[1])
    if verbose:
        print("m3: merged blocks %s%s (n %d): forward %.1f -> %.1f, %d carries renamed, peak carries %d%s" % (
            ([wb[2]["bid"]] if hw else []), [r[2] for r in run], n, before, after, len(rename), nlive[1],
            " + held" if hw else ""))
    return before - after
