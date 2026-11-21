#!/usr/bin/env python3
"""invsqrt_fc2_patch_frac.py -- FC2 patch "frac" for FC's g8 (own code): DIV1's fraction erase without Toffolis.

At DIV1's transition the fe = 2 fraction wires hold F = (sum of the fractional literal bits and constants) mod 1,
and FC subtracts every term again (one borrow Toffoli per literal at 2^-2: 5 per pass on sc8P). Here, forward:
  F_lo (2^-2 wire) = XOR of the 2^-2 entries               -> cleared by CNOTs from the digit wires;
  F_hi (2^-1 wire) = e2(2^-2 entries) + e1(2^-1 entries)    -> X-measured (HMR, bit b); the phase (-1)^(b F_hi) is
     repaid by Z / CZ / NEG gates conditioned on b over the digit wires (F_hi is a quadratic form in them: every
     entry is digit_i xor digit_t xor pol, or a constant 1).
Undo (FC's explicit-undo blocks, BBEG / BEND, Circ.mirror_blocks): the exact mirror of FC's own erase sequence,
so the undo half is unchanged. The block is self-checked in FC's simulator by the full circuit's exactness test.
This file carries a copy of FC's st_div1 (tools/invsqrt_fc_g8.py, 11:24Z) with only the transition replaced.
"""
import itertools

FC_ONE = "1"


def make_st_div1(F8):
    A = F8.A
    K = F8.K
    ONE = F8.ONE

    def transition(c, cfg, digit, R, a, m, hx, fe):
        fracs = []
        for t in range(a, -1, -1):
            if cfg.get("nocorr"):
                break
            bits, Kq = F8.corr_spec(cfg, t, -fe)
            if cfg.get("rowopt") and a > t:
                lits, const = F8.row_lits(bits, Kq, a, t, t + m + 2 + hx, fe, cfg["_rowchoice"][t][0])
            else:
                lits = [(e, i, 1) for e, i in bits]
                const = -Kq
            for e, i, pol in lits:
                if e < 0:
                    fracs.append((e, t, i, pol))
            cmod = const % (1 << fe)
            for p in range(fe):
                if (cmod >> p) & 1:
                    fracs.append((p - fe, t, None, 1))
        assert fe == 2, "patch frac is written for fe = 2"
        # FC's own erase, captured as the undo block's source (its mirror re-creates F from the digits)
        orig_ops = []
        saved = (list(c.ops), list(c.free), c.next, c.live, c.peak, dict(c.stage_peaks))
        c.ops = []
        for e, t, i, pol in fracs:
            if i is not None:
                c.cx(digit[t], digit[i])
                if pol < 0:
                    c.x(digit[i])
                w = digit[i]
            else:
                w = ONE
            F8._sub_frac(c, dict(R), fe, e, w)
            if i is not None:
                if pol < 0:
                    c.x(digit[i])
                c.cx(digit[t], digit[i])
        orig_ops = c.ops
        c.ops, c.free, c.next, c.live, c.peak, c.stage_peaks = saved
        undo = K.Circ.mirror(orig_ops)
        hi_id = max([q for o in undo for q in o[1:] if isinstance(q, int)] + [-1])
        if hi_id >= c.next:                      # a temporary of the undo above every wire in use: count it
            c.next = hi_id + 1
        # entries of the two columns as affine forms over digit wires: (set of wires, constant)
        col = {-2: [], -1: []}
        for e, t, i, pol in fracs:
            if i is None:
                col[e].append((frozenset(), 1))
            else:
                col[e].append((frozenset({digit[i], digit[t]}) if digit[i] != digit[t] else frozenset(),
                               1 if pol < 0 else 0))
        lo, hi = R[-2], R[-1]
        if cfg.get("frac_undo", "fa") == "fa":
            undo = fa_undo(c, col, lo, hi)
        if not hasattr(c, "blocks") or c.blocks is None:
            c.blocks = {}
        bid = ("fc2frac", len(c.blocks))
        c.blocks[bid] = undo
        c.ops.append(("BBEG", bid))
        for ws, k in col[-2]:                      # F_lo -> 0
            for w in ws:
                c.cx(w, lo)
            if k:
                c.x(lo)
        b = c.hmr(hi)                              # F_hi measured away
        # phase repair: (-1)^(b F_hi), F_hi = e2(col -2) + e1(col -1) over GF(2), as a quadratic form
        quad = {}
        lin = {}
        cst = 0

        def add_prod(f, g):
            nonlocal cst
            (ws1, k1), (ws2, k2) = f, g
            for u in ws1:
                for w in ws2:
                    if u == w:
                        lin[u] = lin.get(u, 0) ^ 1
                    else:
                        key = (min(u, w), max(u, w))
                        quad[key] = quad.get(key, 0) ^ 1
            if k2:
                for u in ws1:
                    lin[u] = lin.get(u, 0) ^ 1
            if k1:
                for w in ws2:
                    lin[w] = lin.get(w, 0) ^ 1
            if k1 and k2:
                cst ^= 1
        for f, g in itertools.combinations(col[-2], 2):
            add_prod(f, g)
        for ws, k in col[-1]:
            for u in ws:
                lin[u] = lin.get(u, 0) ^ 1
            cst ^= k
        for (u, w), q in sorted(quad.items()):
            if q:
                c.cz(u, w, cond=b)
        for u, q in sorted(lin.items()):
            if q:
                c.z(u, cond=b)
        if cst:
            c.neg(cond=b)
        c.ops.append(("BEND", bid))
        return col

    def parity_carry(c, col, target):
        """target ^= bit 2 of N, N = sum(col[-2]) + 2 sum(col[-1]) (entries: (set of digit wires, const));
        the carry of the fraction literals into R_0's last bit. In-place formation of the 2^-2 literals along a
        spanning forest of their graph (copies for cycle edges and for every 2^-1 literal), then full adders
        (carry = 1 AND on a fresh wire, sum in place), the weight-4 carries XORed into target, all undone."""
        undo = []

        def cx(a_, b_):
            c.cx(a_, b_)
            undo.append(("cx", a_, b_))

        def xg(a_):
            c.x(a_)
            undo.append(("x", a_))
        fresh = []

        def new():
            w = c.alloc(1)[0]
            fresh.append(w)
            return w
        w1, w2, k1, k2 = [], [], 0, 0
        for ws, k in col[-1]:                     # copies first (they need the original digit values)
            if not ws:
                k2 += k
                continue
            y = new()
            for w in ws:
                cx(w, y)
            if k:
                xg(y)
            w2.append(y)
        # 2^-2 literals in place along a spanning forest (edges = literal wire pairs)
        adj = {}
        edges = []
        for ws, k in col[-2]:
            if not ws:
                k1 += k
                continue
            u_, v_ = sorted(ws)
            edges.append((u_, v_, k))
        parent = {}
        order_ = []
        seen = set()
        tree_edges = set()
        for u_, v_, k in edges:
            adj.setdefault(u_, []).append(v_)
            adj.setdefault(v_, []).append(u_)
        for r0 in sorted(adj):
            if r0 in seen:
                continue
            seen.add(r0)
            q = [r0]
            while q:
                x0 = q.pop(0)
                order_.append(x0)
                for y0 in adj[x0]:
                    if y0 not in seen:
                        seen.add(y0)
                        parent[y0] = x0
                        tree_edges.add((min(x0, y0), max(x0, y0)))
                        q.append(y0)
        holder = {}
        for x0 in reversed(order_):               # leaves first: the parent is still original
            if x0 in parent:
                cx(parent[x0], x0)
                holder[(min(x0, parent[x0]), max(x0, parent[x0]))] = x0
        used = set()
        for u_, v_, k in edges:
            key = (u_, v_)
            if key in holder and key not in used:
                used.add(key)
                w = holder[key]
            else:                                 # a cycle edge or a repeated literal: copy from originals
                # originals: undo the in-place formation is not possible here; recompute from the tree:
                # the value u_ ^ v_ is the XOR of the in-place wires along the tree path between u_ and v_
                path_u, path_v = [u_], [v_]
                while path_u[-1] in parent:
                    path_u.append(parent[path_u[-1]])
                while path_v[-1] in parent:
                    path_v.append(parent[path_v[-1]])
                su, sv = set(path_u), set(path_v)
                w = new()
                for z in path_u[:-1]:
                    if z not in sv or z == path_u[-1]:
                        cx(z, w)
                for z in path_v[:-1]:
                    if z not in su:
                        cx(z, w)
            if k:
                xg(w)
            w1.append(w)
        # weight 1: reduce to one bit, carries go to weight 2
        ones = {1: k1 % 2, 2: (k1 // 2 + k2) % 2, 4: ((k1 // 2 + k2) // 2) % 2}
        cols = {1: list(w1), 2: list(w2), 4: []}

        def fa(a_, b_, cc, dst):
            t = new()
            cx(cc, a_); cx(cc, b_)
            c.and_compute(a_, b_, t)
            undo.append(("and", a_, b_, t))
            cx(cc, t)
            cx(cc, a_); cx(cc, b_)
            cx(a_, cc); cx(b_, cc)
            cols[dst].append(t)

        def ha(a_, b_, dst):
            t = new()
            c.and_compute(a_, b_, t)
            undo.append(("and", a_, b_, t))
            cx(a_, b_)
            cols[dst].append(t)
        for wt, nxt in ((1, 2), (2, 4)):
            if ones[wt] and wt == 2 and cols[2]:
                # x + 1 = NOT x (weight 2) + x (weight 4): the weight-4 part goes straight into the target
                x0 = cols[2][0]
                c.cx(x0, target)
                xg(x0)
            elif ones[wt]:                        # a constant 1 at weight 1: a wire holding 1
                y = new()
                xg(y)
                cols[wt].append(y)
            while len(cols[wt]) >= 3:
                a_, b_, cc = cols[wt].pop(), cols[wt].pop(), cols[wt].pop()
                fa(a_, b_, cc, nxt)
                cols[wt].append(cc)
            if len(cols[wt]) == 2:
                a_, b_ = cols[wt].pop(), cols[wt].pop()
                ha(a_, b_, nxt)
                cols[wt].append(b_)
        for t in cols[4]:
            c.cx(t, target)
        if ones[4]:
            c.x(target)
        for op in reversed(undo):
            if op[0] == "cx":
                c.cx(op[1], op[2])
            elif op[0] == "x":
                c.x(op[1])
            else:
                c.and_erase(op[1], op[2], op[3], free=False)
        c.free_(fresh)

    def fa_undo(c, col, lo, hi):
        """explicit undo: (lo, hi) = (0, 0) -> (F_lo, F_hi). F_lo by CNOTs; F_hi = e2(col -2) + e1(col -1): copies of
        the -2 entries on free wires, full adders (1 AND each, carry = weight 2), hi ^= every carry, then all undone.
        Built on a scratch op list with the allocator state saved and restored (wires free at this point)."""
        saved = (list(c.ops), list(c.free), c.next, c.live, c.peak, dict(c.stage_peaks))
        c.ops = []
        for ws, k in col[-2]:
            for w in ws:
                c.cx(w, lo)
            if k:
                c.x(lo)
        for ws, k in col[-1]:
            for w in ws:
                c.cx(w, hi)
            if k:
                c.x(hi)
        ys = []
        prog = []                                   # forward ops of the e2 computation, mirrored afterwards
        ops0 = len(c.ops)
        for ws, k in col[-2]:
            y = c.alloc(1)[0]
            for w in ws:
                c.cx(w, y)
            if k:
                c.x(y)
            ys.append(y)
        carries = []
        while len(ys) >= 3:
            a, b, cc = ys.pop(), ys.pop(), ys.pop()
            t = c.alloc(1)[0]
            c.cx(cc, a); c.cx(cc, b)
            c.and_compute(a, b, t)
            c.cx(cc, t)
            c.cx(cc, a); c.cx(cc, b)
            c.cx(a, cc); c.cx(b, cc)                # cc = a ^ b ^ c (weight 1), t = carry (weight 2)
            carries.append(t)
            ys.append(cc)
        if len(ys) == 2:
            t = c.alloc(1)[0]
            c.and_compute(ys[0], ys[1], t)
            carries.append(t)
        prog = c.ops[ops0:]
        for t in carries:
            c.cx(t, hi)
        c.ops.extend(K.Circ.mirror(prog))
        for t in carries:
            c.free_([t])
        undo = c.ops
        c.ops, c.free, c.next, c.live, c.peak, c.stage_peaks = saved
        hi_id = max([q for o in undo for q in o[1:] if isinstance(q, int)] + [-1])
        if hi_id >= c.next:
            c.next = hi_id + 1
        return undo

    def st_div1(c, X, cfg):
        """DIV1-local adder regime: cfg["div1_cap"] = N runs DIV1 with Gidney adders under cap N (narrow off),
        the rest of the circuit keeps FC's cap / narrow; P1's squeeze is expected to take the peak back."""
        lc = cfg.get("div1_cap")
        if lc is None:
            return st_div1_body(c, X, cfg)
        saved = (getattr(c, "cap", None), getattr(c, "narrow", False))
        c.cap, c.narrow = int(lc), False
        try:
            return st_div1_body(c, X, cfg)
        finally:
            c.cap, c.narrow = saved

    def st_div1_body(c, X, cfg):
        m, fe, hx, tk, F1 = cfg["m"], cfg["fe"], cfg["hx"], cfg["tk"], cfg["F1"]
        a = F8.a_top(cfg)
        s1 = F8.sched1_of(cfg)
        v = X[22 - m:22]
        cfg["_v0"] = v[0]
        nd = 22 - m
        top_pos = a + m + 3 + hx
        Frw = c.alloc(fe)
        Z = c.alloc(F1)
        R = {}
        for i, w in enumerate(Frw):
            R[-fe + i] = w
        for i, w in enumerate(Z):
            R[i] = w
        for i in range(nd):
            R[F1 + i] = X[i]
        ext = c.alloc(top_pos - (F1 + nd) + 1)
        for i, w in enumerate(ext):
            R[F1 + nd + i] = w
            c.cx(X[nd - 1], w)
        dsrc = [ONE] + list(v) + [ONE]
        digit = {}
        for t in range(a, -1, -1):
            sig = R[t + m + 3 + hx]
            A.fadd(c, sig, dsrc + [None] * (1 + hx), [R[p] for p in range(t, t + m + 3 + hx)])
            digit[t] = sig
            if cfg.get("nocorr"):
                continue
            bits, Kq = F8.corr_spec(cfg, t, -fe)
            for i in range(t + 1, a + 1):
                c.cx(digit[t], digit[i])
            top = t + m + 2 + hx
            if cfg.get("rowopt") and a > t:
                comp = cfg["_rowchoice"][t][0]
                lits, const = F8.row_lits(bits, Kq, a, t, top, fe, comp)
                F8._add_lits(c, [(e, digit[i], pol) for e, i, pol in lits], const, R, -fe, top)
            else:
                F8._add_bits(c, [(e, digit[i]) for e, i in bits], -Kq, R, -fe, top)
            for i in range(t + 1, a + 1):
                c.cx(digit[t], digit[i])
        assert not cfg.get("odd", 0) and not cfg.get("topc", 0)
        col = transition(c, cfg, digit, R, a, m, hx, fe)
        c.free_(Frw)
        tw = m + 2 + hx
        Vb = list(v) + [ONE]
        sig = R[tw]
        tgt = [R[p] for p in range(0, tw)]
        A.frame(c, sig, tgt)
        c.x(sig)
        A.add(c, Vb + [None] * (tw - len(Vb)), tgt, cin=sig)
        c.x(sig)
        A.frame(c, sig, tgt)
        digit[-1] = sig
        u = list(tgt)
        F8.tailcorr(c, cfg, 1, digit, u)
        for k in range(2, len(s1) + 1):
            kk = s1[k - 1]
            u = [u[tw - 1]] + u[0:tw - 1]
            su = u[0]
            src = Vb[kk + 1:]
            A.fadd(c, su, src + [None] * (tw - 1 - kk - len(src)), u[kk + 1:tw], cin=Vb[kk])
            c.cx(v[0], u[0])
            digit[-k] = u[0]
            if k <= tk:
                F8.tailcorr(c, cfg, k, digit, u)
        digit["fin"] = u[tw - 1]
        for k in range(2, len(s1) + 1):
            c.cx(v[0], digit[-k])
        if cfg.get("fc2_u8"):
            # the window bit that carries R_0's parity = affine(digits, v) xor bit 2 of the fraction sum N [M]:
            # XOR that carry in, so FC's relations see an affine bit and clear it (one held wire less)
            parity_carry(c, col, u[int(cfg.get("fc2_u8_pos", 8))])
        order = list(range(a, -1, -1)) + [-k for k in range(1, len(s1) + 1)] + ["fin"]
        return digit, order, u, ext
    return st_div1


def apply(F8, cfg):
    done = F8.__dict__.setdefault("_fc2_done", set())
    if "frac" in done:
        return
    F8.st_div1 = make_st_div1(F8)
    done.add("frac")
