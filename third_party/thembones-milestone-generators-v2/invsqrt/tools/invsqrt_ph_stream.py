"""invsqrt_ph: phase streams (own code, workstream invsqrt_ph).

A stream delivers bit k of a source word at rotation position pos0 - k - d, d = stream delay 0..10 held on tier
control wires. One rotation wire r: per step the next source bit is copied onto r (CNOT), r passes the delay
tiers (controlled swaps with cells that are zero low bits of the normalised input), gets its Z_POW and is measured
out (HMR). Nothing is un-streamed: the measured copy at step n is sum over delays [d = d'] src[n - d'], a
polynomial "monomial of tier controls x source bit", repaired after the stream by conditioned Z / CZ and one
Toffoli per product monomial.

Credits: delay cells in the zero low bits of the normalised input, the delay of 8 as a ring of 8 places visited on
alternate half cycles, one rotation wire for both streams, repair by one CCZ per product class: best public
circuit (manwithacat "HW-25" = gengsr + maoxin123), structure read in research/ws_invsqrt_decode.md 3.4, 3.5.
Temporary AND and measurement erasure: Gidney 1709.06648. Symbolic planner, need analysis, flush steps, repair
polynomial and the special-bit synthesiser are our own.

Tier controls (functions of d): a = d & 1, b = (d >> 1) & 1, c14 = [d >= 4], g = [d >= 8], mid = c14 xor g.
Cells: A0 (tier 1), B0 B1 (tier 2), R0..R3 (ring, delay 4), X0..X3 (extra places, with the ring a delay of 8).
"""
from __future__ import annotations

from fractions import Fraction

ONE = "one"
DELAYS = list(range(11))


def ctl_val(name, d):
    return {"a": d & 1, "b": (d >> 1) & 1, "c14": int(d >= 4), "g": int(d >= 8), "mid": int(4 <= d < 8)}[name]


def tier_swaps(n, off):
    ph = (n + off) % 8
    out = [("a", "A0"), ("b", "B%d" % (n % 2))]
    if ph < 4:
        out.append(("c14", "R%d" % ph))
    else:
        out += [("mid", "R%d" % (ph - 4)), ("g", "X%d" % (ph - 4))]
    return out


def plan(symbols, last_rot, off=0, delays=DELAYS, max_steps=64):
    """symbols[n] = source symbol loaded at step n. Symbolic run for every delay. Returns the list of steps
    {load, swaps: [(ctl, cell)], out: {d: symbol or None}}; checks out[d] at step n is symbol n - d for every step
    up to last_rot (the steps that carry a rotation) and that every cell ends empty (flush steps are appended)."""
    cells = {d: {} for d in delays}
    steps = []
    n = 0
    while True:
        r = {d: (symbols[n] if n < len(symbols) else None) for d in delays}
        sw = []
        for ctl, cell in tier_swaps(n, off):
            on = [d for d in delays if ctl_val(ctl, d)]
            if any(r[d] is not None or cells[d].get(cell) is not None for d in on):
                sw.append((ctl, cell))
                for d in on:
                    r[d], cells[d][cell] = cells[d].get(cell), r[d]
        if n <= last_rot:
            for d in delays:
                want = symbols[n - d] if 0 <= n - d < len(symbols) else None
                assert r[d] == want, (n, d, r[d], want)
        steps.append({"load": symbols[n] if n < len(symbols) else None, "swaps": sw, "out": dict(r)})
        n += 1
        if n >= len(symbols) and n > last_rot and all(v is None for d in delays for v in cells[d].values()):
            break
        assert n < max_steps, "stream does not drain"
    return steps


def best_plan(symbols, last_rot, delays=DELAYS):
    best = None
    for off in range(8):
        try:
            st = plan(symbols, last_rot, off, delays)
        except AssertionError:
            continue
        cost = sum(len(s["swaps"]) for s in st)
        if best is None or cost < best[0]:
            best = (cost, off, st)
    return best


# ------------------------------------------------------------------ monomials of the tier controls

MONO_VARS = ("a", "b", "c14", "g")


def mono_val(mk, d):
    return int(all(ctl_val(v, d) for i, v in enumerate(MONO_VARS) if (mk >> i) & 1))


def mono_basis(delays=DELAYS):
    """Monomials (masks over MONO_VARS, low degree first) independent on the valid delays, and for each delay
    the mask (over the chosen list) of monomials whose XOR is the indicator [d = d']."""
    cands = sorted(range(16), key=lambda m: (bin(m).count("1"), m))
    chosen, rows = [], []          # rows: (vector over delays, combination over chosen)
    piv = {}
    for mk in cands:
        vec = sum(mono_val(mk, d) << i for i, d in enumerate(delays))
        v, comb = vec, 1 << len(chosen)
        for b in sorted(piv, reverse=True):
            if (v >> b) & 1:
                v ^= piv[b][0]
                comb ^= piv[b][1]
        if v:
            piv[v.bit_length() - 1] = (v, comb)
            chosen.append(mk)
        if len(chosen) == len(delays):
            break
    assert len(chosen) == len(delays)
    ind = {}
    for i, d in enumerate(delays):
        v, comb = 1 << i, 0
        for b in sorted(piv, reverse=True):
            if (v >> b) & 1:
                v ^= piv[b][0]
                comb ^= piv[b][1]
        assert v == 0
        ind[d] = comb
    return chosen, ind


# ------------------------------------------------------------------ emission

def rotate(c, q, pos, angles=None):
    if angles is not None and pos in angles:
        c.zpow(q, angles[pos])
    else:
        c.zpow(q, Fraction(1, 1 << pos))


class Stream:
    """One stream. src: dict symbol -> wire (ONE needs no wire). ctl: dict a, b, e, g -> wire (e holds c14 and is
    toggled to mid by a CNOT from g around the swaps that need it). cells: dict cell name -> wire."""

    def __init__(self, c, r, symbols, src, ctl, cells, pos0, min_pos, top_pos, angles=None, delays=DELAYS,
                 ind=None, one_idx=None, ext=None):
        """ind (optional): dict delay -> (wires, constant) whose XOR is the indicator [d = delay]. With it the
        repair needs no Toffoli, and a source bit that is the constant 1 (index one_idx, not in `symbols`) is
        XORed onto the rotation wire by CNOTs at its arrival step instead of travelling through the cells."""
        self.c, self.r, self.symbols, self.src, self.ctl, self.cells = c, r, symbols, src, ctl, cells
        self.pos0, self.min_pos, self.top_pos, self.angles, self.delays = pos0, min_pos, top_pos, angles, delays
        self.ind, self.one_idx, self.one_steps = ind, one_idx, {}
        # ext = (first index, le): constant 1 at every source index >= first index ("sign extension" of a word
        # that is delivered negated); le(q) = (wires, constant) whose XOR is [d <= q]. Needs no Toffoli.
        self.ext, self.ext_steps = ext, {}
        last_rot = pos0 - min_pos
        self.cost, self.off, self.steps = best_plan(symbols, last_rot, delays)
        self.bits = []
        self.n_rot = 0

    def _ctl(self, name):
        return self.ctl["e" if name in ("c14", "mid") else name]

    def run(self, inject=None):
        """inject: dict position -> callable(r) -> callable(m) or None; called after the tiers (it may XOR extra
        bits onto r), the returned callable is called with the outcome bit of that step for its own repair."""
        c, r = self.c, self.r
        for n, st in enumerate(self.steps):
            pos = self.pos0 - n
            if st["load"] == ONE:
                c.x(r)
            elif st["load"] is not None:
                c.cx(self.src[st["load"]], r)
            for name, cell in st["swaps"]:
                if name == "mid":
                    c.cx(self.ctl["g"], self.ctl["e"])
                c.cswap(self._ctl(name), r, self.cells[cell])
                if name == "mid":
                    c.cx(self.ctl["g"], self.ctl["e"])
            if self.ind is not None and self.one_idx is not None and (n - self.one_idx) in self.ind:
                ws, k = self.ind[n - self.one_idx]
                for w in ws:
                    c.cx(w, r)
                if k:
                    c.x(r)
                self.one_steps[n] = n - self.one_idx
            if self.ext is not None and n >= self.ext[0]:
                ws, k = self.ext[1](n - self.ext[0])
                for w in ws:
                    c.cx(w, r)
                if k:
                    c.x(r)
                self.ext_steps[n] = (ws, k)
            fix = inject[pos](r) if inject and pos in inject else None
            if self.min_pos <= pos <= self.top_pos:
                rotate(c, r, pos, self.angles)
                self.n_rot += int(pos >= 3)
            m = c.hmr(r)
            if fix is not None:
                fix(m)
            self.bits.append(m)

    def repair(self, after_g=None):
        """Phase repair of all kept outcomes; r is used as the collecting wire, one temporary for a AND b.
        after_g: called once every monomial that contains the control g has been repaired (the caller may then
        release the wire of g, so that the temporary for a AND b does not add to the wire count)."""
        c, r = self.c, self.r
        if self.ind is not None:                              # indicators on wires: conditioned Z / CZ only
            for n, st in enumerate(self.steps):
                m = self.bits[n]
                for d in self.delays:
                    s = st["out"][d]
                    if s is None:
                        continue
                    ws, k = self.ind[d]
                    for w in ws:
                        c.cz(w, self.src[s], cond=m)
                    if k:
                        c.z(self.src[s], cond=m)
                if n in self.one_steps:
                    ws, k = self.ind[self.one_steps[n]]
                    for w in ws:
                        c.z(w, cond=m)
                    if k:
                        c.neg(cond=m)
                if n in self.ext_steps:
                    ws, k = self.ext_steps[n]
                    for w in ws:
                        c.z(w, cond=m)
                    if k:
                        c.neg(cond=m)
            for m in self.bits:
                c.bit_free(m)
            self.bits = []
            return
        monos, ind = mono_basis(self.delays)
        # coef[mono index][n] = set of symbols (XOR)
        coef = [dict() for _ in monos]
        for n, st in enumerate(self.steps):
            for d in self.delays:
                s = st["out"][d]
                if s is None:
                    continue
                for i in range(len(monos)):
                    if (ind[d] >> i) & 1:
                        cur = coef[i].setdefault(n, set())
                        cur ^= {s}
        ab = None
        order = sorted(range(len(monos)), key=lambda i: (0 if monos[i] & 8 else 1, bin(monos[i]).count("1"), monos[i]))
        pending_hook = after_g
        for i in order:
            mk = monos[i]
            if pending_hook is not None and not mk & 8:
                pending_hook()
                pending_hook = None
            vs = [self._ctl(v) for k, v in enumerate(MONO_VARS) if (mk >> k) & 1]
            items = [(n, s) for n in sorted(coef[i]) for s in sorted(coef[i][n], key=str)]
            if not items:
                continue
            if len(vs) == 0:
                for n, s in items:
                    if s == ONE:
                        c.neg(cond=self.bits[n])
                    else:
                        c.z(self.src[s], cond=self.bits[n])
            elif len(vs) == 1:
                for n, s in items:
                    if s == ONE:
                        c.z(vs[0], cond=self.bits[n])
                    else:
                        c.cz(vs[0], self.src[s], cond=self.bits[n])
            else:
                def collect():
                    for n, s in items:
                        if s == ONE:
                            c.x(r, cond=self.bits[n])
                        else:
                            c.cx(self.src[s], r, cond=self.bits[n])
                collect()
                if mk & 3 == 3:                       # contains a AND b
                    if ab is None:
                        ab = c.and_new(self.ctl["a"], self.ctl["b"])
                    if len(vs) == 2:
                        c.cz(ab, r)
                    else:
                        assert len(vs) == 3
                        c.ccz(ab, vs[2], r)
                else:
                    assert len(vs) == 2, mk
                    c.ccz(vs[0], vs[1], r)
                collect()
        if pending_hook is not None:
            pending_hook()
        if ab is not None:
            c.and_del(ab, self.ctl["a"], self.ctl["b"])
        for m in self.bits:
            c.bit_free(m)
        self.bits = []
