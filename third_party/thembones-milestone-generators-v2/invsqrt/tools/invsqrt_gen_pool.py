"""invsqrt_gen: pool lookups (own code, workstream invsqrt_gen).

Idea credited to the best public invsqrt circuit as dissected in research/ws_invsqrt_new.md: a table wire is
built as a XOR of products of two GF(2)-linear forms over all wires that currently hold functions of the address
(the "pool"), and a wire whose content is a quadratic form over the remaining pool is erased by HMR plus
conditioned CZs at no Toffoli cost. The synthesiser, the solver and the gate emission below are ours; the tables
come from our own fit.

A function of the 8-bit address is a 256-bit Python int (bit a = value at address a).

REPLACEABLE LOOKUP BLOCKS (interface for an external synthesiser)
-----------------------------------------------------------------
Every lookup of the generator goes through lookup_block(). For each block it writes
    ~/qcc/ws/invsqrt_gen/lut/<block>.spec.json   what the block must compute
    ~/qcc/ws/invsqrt_gen/lut/<block>.own.json    the program our own synthesiser found (same format as below)
and, if a file ~/qcc/ws/invsqrt_gen/lut/<block>.prog.json with the same "spec_sha" exists, it replays that
program instead of running our synthesiser. The replay goes through the same truth-table bookkeeping, so a
program that does not produce the specified outputs stops the build.

spec:   {"block", "spec_sha", "n_addr": 8, "address": text,
         "inputs":  [{"name", "tt"}],   wires that are alive when the block starts, with their content
         "outputs": [{"name", "tt"}],   wires the block must fill; they start at |0>, or, when the entry has
                                         "init" (pre-filled output, since 2026-11-12), with that function of
                                         the address on them (pool until the program finishes them)
         "max_scratch": n}              clean wires the block may use in addition ("s0", "s1", ...)
        "tt" is a 64-digit hex string; bit a of int(tt, 16) is the value at address a. Address bit k (k = 0..6)
        is bit 14 + k of the normalised input xn = x << (21 - e); address bit 7 is e & 1, e = floor(log2 x).
        In the first block the inputs are the eight address wires "a0".."a7"; later blocks list the table wires
        that are alive instead (the address wires are erased after the first block in the address-free schedule).
program: {"block", "spec_sha", "n_scratch": s, "steps": [...]}
        ["and", target, form1, form2]   target ^= (XOR of form1) AND (XOR of form2)      one Toffoli
        ["xor", target, form]           target ^= XOR of form                             free
        A form is a list of wire names; the name "1" is the constant one. The target must not occur in the forms
        of its own step. Forms may use any input, any output (with whatever it holds at that moment) and any
        scratch wire. After the last step every output must hold its "tt". Scratch wires may hold anything: the
        generator keeps them as pool and erases them itself when their content is a quadratic form over the
        other live wires (HMR + conditioned CZ), which is the case for any product of two pool forms.
The cost of a block is its number of "and" steps; the wire budget is len(outputs) + max_scratch.
"""
from __future__ import annotations

import hashlib
import json
import os
import random

N = 256
ONE = (1 << N) - 1
LUT_DIR = os.environ.get("INVSQRT_GEN_LUT") or os.path.expanduser("~/qcc/ws/invsqrt_gen/lut")
# INVSQRT_GEN_LUT: another lookup folder (specs, programs, cache) for experiments, so that lut/ of the work dir,
# which the unattended loops read, is not touched by them


_KERNEL = None


def _kernel():
    """C kernel of the last-multiplier search (tools/invsqrt_gen_kernel.c), compiled on first use into
    ~/qcc/ws/invsqrt_gen/kernel/. Returns the ctypes library or None (no compiler: the Python walk is used; it
    finds the same candidates, only slower). INVSQRT_GEN_NOKERNEL=1 switches the library off."""
    global _KERNEL
    if _KERNEL is None:
        _KERNEL = False
        if os.environ.get("INVSQRT_GEN_NOKERNEL") == "1":
            return None
        try:
            import ctypes
            import subprocess
            src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "invsqrt_gen_kernel.c")
            so = os.path.expanduser("~/qcc/ws/invsqrt_gen/kernel/libisk.so")
            if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
                os.makedirs(os.path.dirname(so), exist_ok=True)
                tmp = so + f".{os.getpid()}.tmp"
                subprocess.run(["gcc", "-O3", "-march=native", "-shared", "-fPIC", "-o", tmp, src], check=True,
                               capture_output=True)
                os.replace(tmp, so)
            lib = ctypes.CDLL(so)
            lib.isk_search.restype = ctypes.c_int64
            lib.isk_search.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                                       ctypes.c_void_p, ctypes.c_int64, ctypes.c_int64]
            _KERNEL = lib
        except Exception:
            _KERNEL = False
    return _KERNEL or None


def tt_bit(table, bit):
    """Truth table of one bit of a 256-entry integer table."""
    v = 0
    for a in range(N):
        if (table[a] >> bit) & 1:
            v |= 1 << a
    return v


class Span:
    """Incremental GF(2) basis of N-bit vectors; every vector carries a tag (XOR of generator bits)."""

    def __init__(self):
        self.lead = {}

    def add(self, vec, tag):
        lead = self.lead
        while vec:
            h = vec.bit_length() - 1
            e = lead.get(h)
            if e is None:
                lead[h] = (vec, tag)
                return True
            vec ^= e[0]
            tag ^= e[1]
        return False

    def reduce(self, vec):
        tag = 0
        lead = self.lead
        while vec:
            h = vec.bit_length() - 1
            e = lead.get(h)
            if e is None:
                return vec, tag
            vec ^= e[0]
            tag ^= e[1]
        return 0, tag


def bits_of(tag):
    out = []
    while tag:
        low = tag & -tag
        out.append(low.bit_length() - 1)
        tag ^= low
    return out


class NotQuadratic(Exception):
    pass


class Pool:
    """Tracks which wires hold which function of the address and emits pool ANDs / free erasures into a Circ."""

    def __init__(self, c, rng, mult_weight=0):
        self.c = c
        self.rng = rng
        self.val = {}
        self.ands = 0
        self.n = N               # dimension of the space of address functions (198 in the zone format, where the
                                 # functions are kept on 256 formal addresses with the unused ones copied)
        self.redundant = {}      # zone format: unused address -> the address whose value it copies
        self.mult_weight = mult_weight
        self.rec = None          # list of recorded steps while a block is being synthesised
        self.names = {}          # wire -> name, for recording and replay
        self.empties = []        # wires waiting for a free fill (solver 2)

    def _form_names(self, A, const):
        return sorted(self.names[w] for w in A) + (["1"] if const else [])

    def run_program(self, steps, wires):
        """Replay an external program (see the module docstring); wires: name -> wire."""
        for st in steps:
            t = wires[st[1]]
            forms = []
            for f in st[2:]:
                forms.append(({wires[n] for n in f if n != "1"}, sum(1 for n in f if n == "1") & 1))
            if st[0] == "and":
                self.and_into(t, forms[0][0], forms[0][1], forms[1][0], forms[1][1])
            elif st[0] == "xor":
                self.xor_into(t, forms[0][0], forms[0][1])
            else:
                raise ValueError(st[0])

    # ---- bookkeeping
    def set(self, w, val):
        self.val[w] = val

    def drop(self, w):
        self.val.pop(w, None)

    def form_val(self, wires, const=0):
        v = ONE if const else 0
        for w in wires:
            v ^= self.val[w]
        return v

    def rand_form(self, wires, weight=None):
        """Random GF(2)-linear form over `wires`: dense (every wire with probability 1/2) by default, because
        sparse forms keep the algebraic degree of the pool low and a low-degree pool cannot reach a table bit;
        `weight` caps the number of wires."""
        weight = weight or self.mult_weight
        while True:
            if weight <= 0:
                A = {w for w in wires if self.rng.random() < 0.5}
            else:
                A = set(self.rng.sample(wires, min(len(wires), self.rng.randint(1, weight))))
            v = self.form_val(A)
            if v not in (0, ONE):
                return A

    # ---- gates
    def xor_into(self, t, A, const=0):
        for w in sorted(A):
            self.c.cx(w, t)
        if const:
            self.c.x(t)
        self.val[t] = self.val.get(t, 0) ^ self.form_val(A, const)
        if self.rec is not None and (A or const):
            self.rec.append(["xor", self.names[t], self._form_names(A, const)])

    def and_into(self, t, A, ca, B, cb):
        """t ^= (xor of A ^ ca) AND (xor of B ^ cb); t must not occur in A or B. One Toffoli."""
        A, B = set(A), set(B)
        assert t not in A and t not in B
        fa, fb = self.form_val(A, ca), self.form_val(B, cb)
        if fa == 0 or fb == 0:
            return
        if A == B:
            if ca == cb:
                self.xor_into(t, A, ca)
            return
        if A - B:
            X, cx_, Y, cy = A, ca, B, cb
        else:
            X, cx_, Y, cy = B, cb, A, ca
        p1, p2 = min(X - Y), min(Y)
        ops = [("cx", w, p2) for w in sorted(Y - {p2})]
        if cy:
            ops.append(("x", p2))
        S, cst = set(X - {p1}), cx_
        if p2 in S:
            S.remove(p2)
            S ^= (Y - {p2})
            ops.append(("cx", p2, p1))
            cst ^= cy
        ops += [("cx", w, p1) for w in sorted(S)]
        if cst:
            ops.append(("x", p1))
        c = self.c
        for o in ops:
            c.cx(o[1], o[2]) if o[0] == "cx" else c.x(o[1])
        c.ccx(p1, p2, t)
        for o in reversed(ops):
            c.cx(o[1], o[2]) if o[0] == "cx" else c.x(o[1])
        self.val[t] = self.val.get(t, 0) ^ (fa & fb)
        self.ands += 1
        if self.rec is not None:
            self.rec.append(["and", self.names[t], self._form_names(A, ca), self._form_names(B, cb)])

    # ---- free erasure
    def quad_solve(self, f, others):
        """Express f as const + linear + quadratic form over `others`; returns list of generators or None."""
        sp, gens = Span(), []

        def add(vec, g):
            gens.append(g)
            sp.add(vec, 1 << (len(gens) - 1))

        add(ONE, ("c",))
        for a in others:
            add(self.val[a], ("l", a))
        rem, tag = sp.reduce(f)
        if rem:
            for i, a in enumerate(others):
                va = self.val[a]
                for b in others[i + 1:]:
                    add(va & self.val[b], ("q", a, b))
                if len(sp.lead) == self.n:
                    break
            rem, tag = sp.reduce(f)
        if rem:
            return None
        return [gens[i] for i in bits_of(tag)]

    def erase(self, w, others=None):
        """Measure w away; its content must be a quadratic form over the other pool wires (0 Toffolis)."""
        others = sorted(k for k in self.val if k != w) if others is None else sorted(others)
        sol = self.quad_solve(self.val[w], others)
        if sol is None:
            raise NotQuadratic(w)
        c = self.c
        m = c.hmr(w)
        for g in sol:
            if g[0] == "c":
                c.neg(cond=m)
            elif g[0] == "l":
                c.z(g[1], cond=m)
            else:
                c.cz(g[1], g[2], cond=m)
        c.bit_free(m)
        del self.val[w]

    # ---- deferred erasure (the public file keeps outcomes "for the unlookup", ws_invsqrt_decode.md 2b.1)
    def phase_basis(self, others):
        """Span of {1, w, w w'} over `others`, built once, for fix_phase()."""
        others = sorted(others)
        sp, gens = Span(), []

        def add(vec, g):
            gens.append(g)
            sp.add(vec, 1 << (len(gens) - 1))

        add(ONE, ("c",))
        for a in others:
            add(self.val[a], ("l", a))
        for i, a in enumerate(others):
            va = self.val[a]
            for b in others[i + 1:]:
                add(va & self.val[b], ("q", a, b))
        return sp, gens

    def fix_phase(self, f, cond, basis):
        """Apply the phase (-1)^(cond AND f(address)) with conditioned NEG / Z / CZ over the wires of `basis`.
        Used for a table wire that was measured earlier with its outcome kept in `cond`."""
        sp, gens = basis
        rem, tag = sp.reduce(f)
        if rem:
            raise NotQuadratic("deferred")
        c = self.c
        for i in bits_of(tag):
            g = gens[i]
            if g[0] == "c":
                c.neg(cond=cond)
            elif g[0] == "l":
                c.z(g[1], cond=cond)
            else:
                c.cz(g[1], g[2], cond=cond)

    # ---- synthesis
    def finalize(self, k, target, max_r=40, tries=8):
        """Add products to wire k until it holds `target`: r random multipliers, co-factors and the linear part
        by one linear solve over the pool. Returns the number of ANDs used."""
        others = sorted(w for w in self.val if w != k)
        d = len(others) + 1
        diff = target ^ self.val.get(k, 0)
        base, gens0 = Span(), [("c",)]
        base.add(ONE, 1)
        for w in others:
            gens0.append(("l", w))
            base.add(self.val[w], 1 << (len(gens0) - 1))
        for r in range(0, max_r + 1):
            # span of {1, w_j, u_i w_j}: (r + 1)(d - 1) + 1 generators minus the relations u_i u_i' (twice)
            # and u_i u_i = u_i
            reach = (r + 1) * (d - 1) + 1 - r * (r + 1) // 2
            if reach < self.n - 6 and r > 0:
                continue                         # success probability 2^-(N - reach) per draw: not worth trying
            for _ in range((6 * tries if reach < self.n + 3 else tries) if r else 1):
                sp = Span()
                sp.lead = dict(base.lead)
                gens = list(gens0)
                us = [self.rand_form(others) for _ in range(r)]
                for i, A in enumerate(us):
                    fa = self.form_val(A)
                    for w in others:
                        gens.append(("p", i, w))
                        sp.add(fa & self.val[w], 1 << (len(gens) - 1))
                rem, tag = sp.reduce(diff)
                if rem:
                    continue
                lin, const, V = set(), 0, [set() for _ in range(r)]
                for i in bits_of(tag):
                    g = gens[i]
                    if g[0] == "c":
                        const ^= 1
                    elif g[0] == "l":
                        lin ^= {g[1]}
                    else:
                        V[g[1]] ^= {g[2]}
                before = self.ands
                for A, Vi in zip(us, V):
                    if Vi:
                        self.and_into(k, A, 0, Vi, 0)
                self.xor_into(k, lin, const)
                assert self.val[k] == target
                return self.ands - before
        raise RuntimeError(f"finalize failed for wire {k}")

    # ---- solver 2: target freedom + search for the last multiplier
    def _search_last(self, others, base_lead, diff, r, budget):
        """r - 1 random multipliers are fixed; the last one is searched by a Gray-code walk. Everything is
        projected onto the quotient by the span of the fixed part, where a candidate multiplier u contributes
        the columns pi(u w_j) = XOR over the wires i of u of pi(w_i w_j). Returns the list of r forms or None."""
        rng = self.rng
        fixed = [self.rand_form(others) for _ in range(r - 1)]
        lead = dict(base_lead)

        def insert(v):
            while v:
                h = v.bit_length() - 1
                e = lead.get(h)
                if e is None:
                    lead[h] = v
                    return
                v ^= e

        val = [self.val[w] for w in others]
        for A in fixed:
            fa = self.form_val(A)
            for v in val:
                insert(fa & v)
        order = sorted(lead, reverse=True)
        basis = [lead[h] for h in order]

        def proj(v):
            for h, bv in zip(order, basis):
                if (v >> h) & 1:
                    v ^= bv
            return v

        b = proj(diff)
        if b == 0:
            return fixed
        d = len(val)
        P = [[0] * d for _ in range(d)]
        for i in range(d):
            for j in range(i, d):
                P[i][j] = P[j][i] = proj(val[i] & val[j])
        alpha = [rng.random() < 0.5 for _ in range(d)]
        perm = list(range(d))
        rng.shuffle(perm)
        budget = min(budget, 1 << d)         # the walk visits every multiplier once in 2^d steps
        lib = _kernel()
        # coordinates of the quotient. Zone format: an unused address is a copy of another one (self.redundant),
        # so it is no coordinate of its own, and where the copy is a pivot the original is zero after proj().
        dead = {s for h, s in self.redundant.items() if h in lead}
        free_pos = [h for h in range(N) if h not in lead and h not in self.redundant and h not in dead]
        if lib is not None and len(free_pos) <= 128 and d <= 64:
            # same walk in C (invsqrt_gen_kernel.c): vectors packed onto the quotient coordinates
            import ctypes

            def pack(v):
                out = 0
                for k, h in enumerate(free_pos):
                    if (v >> h) & 1:
                        out |= 1 << k
                return out

            m64 = (1 << 64) - 1
            flat = (ctypes.c_uint64 * (2 * d * d))()
            for i in range(d):
                for j in range(i, d):
                    pv = pack(P[i][j])
                    for a_, b_ in ((i, j), (j, i)):
                        flat[2 * (a_ * d + b_)] = pv & m64
                        flat[2 * (a_ * d + b_) + 1] = pv >> 64
            pb = pack(b)
            barr = (ctypes.c_uint64 * 2)(pb & m64, pb >> 64)
            aarr = (ctypes.c_uint8 * d)(*[1 if a_ else 0 for a_ in alpha])
            parr = (ctypes.c_int32 * d)(*perm)
            t0 = 0
            while t0 < budget:
                t = lib.isk_search(d, flat, barr, aarr, parr, t0, budget)
                if t < 0:
                    return None
                A = {others[i] for i in range(d) if aarr[i]}
                if self.form_val(A) not in (0, ONE):
                    return fixed + [A]
                t0 = t + 1
            return None
        cols = [0] * d
        for i in range(d):
            if alpha[i]:
                Pi = P[i]
                for j in range(d):
                    cols[j] ^= Pi[j]
        for t in range(budget):
            if t:
                i = perm[(t & -t).bit_length() - 1]
                alpha[i] = not alpha[i]
                Pi = P[i]
                for j in range(d):
                    cols[j] ^= Pi[j]
            led = {}
            for v in cols:
                while v:
                    h = v.bit_length() - 1
                    e = led.get(h)
                    if e is None:
                        led[h] = v
                        break
                    v ^= e
            x = b
            while x:
                e = led.get(x.bit_length() - 1)
                if e is None:
                    break
                x ^= e
            if x == 0:
                A = {others[i] for i in range(d) if alpha[i]}
                if self.form_val(A) not in (0, ONE):
                    return fixed + [A]
        return None

    def finalize2(self, k, target, free, deep):
        """Like finalize(), with two additions (both in the lookup workstream's bookkeeping as "R" and "E"):
        the wire may be finished to target xor any combination of the targets in `free` (wires finished later;
        the caller removes the combination by CNOTs at the end), and when the reach of r products misses by
        delta <= deep dimensions the last multiplier is searched over about 2^(delta + 2) candidates.
        Returns the set of free wires whose target was mixed in."""
        others = sorted(w for w in self.val if w != k and self.val[w])
        d = len(others)
        diff = target ^ self.val.get(k, 0)
        base, gens0 = Span(), [("c",)]
        base.add(ONE, 1)
        for w in others:
            gens0.append(("l", w))
            base.add(self.val[w], 1 << (len(gens0) - 1))
        for j, tj in free:
            gens0.append(("t", j))
            base.add(tj, 1 << (len(gens0) - 1))
        rank0 = len(base.lead)
        base_lead = {h: e[0] for h, e in base.lead.items()}
        for r in range(0, 60):
            reach = rank0 + r * d - r * (r + 1) // 2
            delta = self.n - reach
            if r and delta > deep:
                continue
            us = None
            if r == 0 or delta <= 2:
                cand = [[self.rand_form(others) for _ in range(r)] for _ in range(48 if r else 1)]
            else:
                found = None
                for _ in range(2):
                    found = self._search_last(others, base_lead, diff, r, 1 << min(delta + 2, deep + 2))
                    if found is not None:
                        break
                cand = [found] if found is not None else []
            for us in cand:
                sp = Span()
                sp.lead = dict(base.lead)
                gens = list(gens0)
                for i, A in enumerate(us):
                    fa = self.form_val(A)
                    for w in others:
                        gens.append(("p", i, w))
                        sp.add(fa & self.val[w], 1 << (len(gens) - 1))
                rem, tag = sp.reduce(diff)
                if rem:
                    continue
                lin, const, V, J = set(), 0, [set() for _ in us], set()
                for i in bits_of(tag):
                    g = gens[i]
                    if g[0] == "c":
                        const ^= 1
                    elif g[0] == "l":
                        lin ^= {g[1]}
                    elif g[0] == "t":
                        J ^= {g[1]}
                    else:
                        V[g[1]] ^= {g[2]}
                prods = [(A, Vi) for A, Vi in zip(us, V) if Vi]
                for n, (A, Vi) in enumerate(prods):
                    self.and_into(k, A, 0, Vi, 0)
                    if self.empties and n + 1 < len(prods):
                        # "free fill" (seen in the public lookup by the lut_synth workstream): the partial sum is
                        # a new function; a CNOT copy into a still empty wire makes it a pool element at no cost
                        self.xor_into(self.empties.pop(0), {k})
                self.xor_into(k, lin, const)
                want = target
                for j, tj in free:
                    if j in J:
                        want ^= tj
                assert self.val[k] == want
                return J
        raise RuntimeError(f"finalize2 failed for wire {k}")

    def synth(self, new_wires, targets, deep=0, empty=0, pre=frozenset()):
        """new_wires: clean wires. Phase 1 puts one random product into each (pool building); phase 2 finalises
        the wires listed in `targets` (wire -> function). Wires not in targets are scratch with arbitrary content.
        deep = 0: random multipliers only (solver 1). deep > 0: solver 2 with searches up to 2^(deep + 2).
        empty (solver 2 only): the last `empty` wires of new_wires get no product of their own; they are filled by
        CNOT copies of partial sums while the first outputs are finished.
        pre (solver 2 only; owner B1, 2026-11-12): output wires that are NOT clean but already hold a function of
        the address (pre-filled outputs). They get no product in phase 1, are pool while the other outputs are
        made, and are finished in their turn from what they hold (target xor content)."""
        assert not pre or deep, "pre-filled outputs need solver 2"
        # Table ties (workstream invsqrt_ties make such tables): a target that is an XOR of pool wires and of
        # other targets of this block needs no product. These wires are set aside and made by CNOTs at the end.
        # Without such a target nothing changes (no random number is drawn here), so cached programs stay valid.
        copies = {}
        if deep:
            base, gens = Span(), [None]
            base.add(ONE, 1)
            for w in sorted(w for w in self.val if self.val[w] and w not in pre):
                gens.append(w)
                base.add(self.val[w], 1 << (len(gens) - 1))
            for k in new_wires:
                if k not in targets:
                    continue
                rem, tag = base.reduce(targets[k])
                if rem == 0 and k not in pre:
                    copies[k] = [gens[i] for i in bits_of(tag)]
                else:
                    gens.append(k)
                    base.add(targets[k], 1 << (len(gens) - 1))
            new_wires = [k for k in new_wires if k not in copies]
        self.empties = []
        if deep and empty:
            cand = [k for k in new_wires if k not in pre]
            self.empties = list(cand[len(cand) - empty:])
            for k in self.empties:
                self.val[k] = 0
        for k in new_wires:
            if k in self.empties or k in pre:
                continue
            below = sorted(w for w in self.val if self.val[w])
            A, B = self.rand_form(below), self.rand_form(below)
            while self.form_val(A) & self.form_val(B) in (0, self.form_val(A), self.form_val(B)):
                A, B = self.rand_form(below), self.rand_form(below)
            self.val[k] = 0
            self.and_into(k, A, 0, B, 0)
        order = [k for k in new_wires if k in targets]
        if not deep:
            return [self.finalize(k, targets[k]) for k in order]
        mixed = {}
        for n, k in enumerate(order):
            if k in self.empties:
                self.empties.remove(k)
            mixed[k] = self.finalize2(k, targets[k], [(j, targets[j]) for j in order[n + 1:]], deep)
        for k in reversed(order):                # the wires finished later are exact by now
            if mixed[k]:
                self.xor_into(k, mixed[k])
            assert self.val[k] == targets[k]
        for k, comb in copies.items():           # tied targets: CNOT copies of finished wires
            self.val[k] = 0
            self.xor_into(k, {w for w in comb if w is not None}, sum(1 for w in comb if w is None) & 1)
            assert self.val[k] == targets[k]
        return None


def tt_hex(v):
    return f"{v:064x}"


SOLVER_REV = 2      # bump when the synthesiser changes, so that cached programs of older code are not replayed


def lookup_block(c, pool, block, inputs, outputs, max_scratch, seed=1, log=None, deep=0, empty=0):
    """One replaceable lookup block (interface in the module docstring).
    inputs:  list of (name, wire) for every pool wire alive at the start of the block.
    outputs: list of (name, wire, truth table); wires are clean and allocated by the caller.
    deep:    0 = solver 1 (random multipliers); > 0 = solver 2 (target freedom, searches up to 2^(deep + 2)).
             Programs found by solver 2 are cached in lut/cache/ (the search is deterministic, the cache only
             saves time) and replayed on later builds.
    Returns the scratch wires (left in the pool with whatever they hold; the caller erases or keeps them)."""
    os.makedirs(LUT_DIR, exist_ok=True)
    # pre-filled outputs (owner B1, 2026-11-12): an output wire that is already a pool wire holds a function of the
    # address ("init" in the spec); the program must bring it from that content to its "tt". Without such an
    # output the spec (and its hash) is exactly as before.
    init = {w: pool.val[w] for _, w, _ in outputs if pool.val.get(w)}
    assert {w for _, w in inputs} == set(pool.val) - set(init), "inputs must list exactly the live pool wires"

    def _out(n, w, v):
        d = dict(name=n, tt=tt_hex(v))
        if w in init:
            d["init"] = tt_hex(init[w])
        return d
    spec = dict(block=block, n_addr=8,
                address="bit k (0..6) = xn[14+k], bit 7 = e & 1; tt bit a = value at address a",
                inputs=[dict(name=n, tt=tt_hex(pool.val[w])) for n, w in inputs],
                outputs=[_out(n, w, v) for n, w, v in outputs], max_scratch=max_scratch)
    sha = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:16]
    spec["spec_sha"] = sha
    with open(os.path.join(LUT_DIR, block + ".spec.json"), "w") as fh:
        json.dump(spec, fh)
    names = {w: n for n, w in inputs}
    names.update({w: n for n, w, _ in outputs})
    ext = os.path.join(LUT_DIR, block + ".prog.json")
    prog = None
    if os.path.exists(ext):
        with open(ext) as fh:
            prog = json.load(fh)
        if prog.get("spec_sha") != sha:
            print(f"  [{block}] external program ignored: spec_sha {prog.get('spec_sha')} != {sha}")
            prog = None
    before = pool.ands
    if prog is not None:
        ns = int(prog.get("n_scratch", 0))
        assert ns <= max_scratch, "program uses more scratch wires than the budget"
        scr = c.allocn(ns)
        wires = {n: w for w, n in names.items()}
        wires.update({f"s{i}": w for i, w in enumerate(scr)})
        for _, w, _ in outputs:
            pool.val[w] = init.get(w, 0)
        for w in scr:
            pool.val[w] = 0
        pool.run_program(prog["steps"], wires)
        # a program of our own solver that invsqrt_gen_lutcmp.py installed from the cache says so in "origin"
        own_inst = str(prog.get("origin", "")).startswith("own solver")
        src = (prog["origin"] if own_inst else
               "external " + hashlib.sha256(json.dumps(prog, sort_keys=True).encode()).hexdigest()[:12])
    else:
        scr = c.allocn(max_scratch)
        names.update({w: f"s{i}" for i, w in enumerate(scr)})
        etag = f".e{empty}" if empty else ""
        cache = os.path.join(LUT_DIR, "cache", f"{block}.{sha}.s{seed}.d{deep}{etag}.r{SOLVER_REV}.json")
        if deep and os.path.exists(cache):
            with open(cache) as fh:
                own = json.load(fh)
            wires = {n: w for w, n in names.items()}
            for _, w, _ in outputs:
                pool.val[w] = init.get(w, 0)
            for w in scr:
                pool.val[w] = 0
            pool.run_program(own["steps"], wires)
            src = f"own d{deep} (cached)"
        else:
            pool.names, pool.rec = names, []
            pool.rng = random.Random(f"{seed}:{block}")
            pool.synth([w for _, w, _ in outputs] + scr, {w: v for _, w, v in outputs}, deep=deep, empty=empty,
                       pre=frozenset(init))
            own = dict(block=block, spec_sha=sha, n_scratch=max_scratch, steps=pool.rec)
            pool.rec = None
            if deep:
                os.makedirs(os.path.dirname(cache), exist_ok=True)
                with open(cache, "w") as fh:
                    json.dump(own, fh)
            src = f"own d{deep}" if deep else "own"
        with open(os.path.join(LUT_DIR, block + ".own.json"), "w") as fh:
            json.dump(own, fh)
    for n, w, v in outputs:
        assert pool.val[w] == v, f"{block}: output {n} is wrong"
    if log is not None:
        log[block] = dict(ands=pool.ands - before, outputs=len(outputs), scratch=len(scr), source=src, spec_sha=sha)
    return scr
