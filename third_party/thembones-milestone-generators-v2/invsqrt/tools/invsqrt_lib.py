"""invsqrt workstream library: kmx parsing, exact bitset simulation over many x, reference values.

A circuit is a list of ops: (name, qs, bit, cond, angle) where
  name  : base op name (X, CX, CCX, Z, CZ, CCZ, Z_POW, NEG, SWAP, CSWAP, R, HMR, BIT_INVERT, BIT_STORE0,
          BIT_STORE1, PUSH_CONDITION, POP_CONDITION)
  qs    : tuple of qubit ids
  bit   : target classical bit id or None
  cond  : condition bit id or None
  angle : Fraction (half turns) for Z_POW else None
The simulator keeps one Python int per qubit / bit: bit k of the int is the value in shot k.
"""
from __future__ import annotations

import math
import random
from fractions import Fraction

import numpy as np

MAX_X = 3 * (2 ** 10 - 2) ** 2  # 3,133,452
TOF = ("CCX", "CCZ", "CSWAP")


def parse_kmx(path):
    header, ops = [], []
    with open(path) as fh:
        for line in fh:
            s = line.split("#")[0].strip()
            if not s:
                continue
            tok = s.split()
            name = tok[0]
            if name in ("REGISTER", "APPEND_TO_REGISTER"):
                header.append(s)
                continue
            cond = None
            if "if" in tok:
                i = tok.index("if")
                cond = int(tok[i + 1][1:])
                tok = tok[:i]
            qs, bit, angle = [], None, None
            for t in tok[1:]:
                if t[0] == "q":
                    qs.append(int(t[1:]))
                elif t[0] == "b":
                    bit = int(t[1:])
                else:
                    angle = Fraction(t)
            ops.append((name, tuple(qs), bit, cond, angle))
    return header, ops


def op_text(op):
    name, qs, bit, cond, angle = op
    s = name + "".join(f" q{q}" for q in qs)
    if bit is not None:
        s += f" b{bit}"
    if angle is not None:
        s += " " + frac_text(angle)
    if cond is not None:
        s += f" if b{cond}"
    return s


def frac_text(a: Fraction) -> str:
    a = a % 2
    n, d = a.numerator, a.denominator
    assert d & (d - 1) == 0, "angle must be dyadic"
    k = d.bit_length() - 1
    digits = str(n * 5 ** k).rjust(k, "0") if k else ""
    ip = digits[:-k] if k else str(n)
    fp = digits[-k:] if k else "0"
    return (ip or "0") + "." + fp


def write_kmx(path, header, ops):
    with open(path, "w") as fh:
        for h in header:
            fh.write(h + "\n")
        for op in ops:
            fh.write(op_text(op) + "\n")


def kmx_to_kmb(kmx, kmb):
    import sys
    if "<HOME>/qcc/data/kickmix" not in sys.path:
        sys.path.insert(0, "<HOME>/qcc/data/kickmix")
    import kickmix as km
    c = km.Circuit.from_file(kmx, format="kmx")
    c.to_file(kmb, format="kmb")
    return c


def pack_inputs(xs, width=22):
    """xs: iterable of ints -> list of per-qubit bitsets."""
    xs = np.asarray(xs, dtype=np.uint64)
    out = []
    for k in range(width):
        bits = ((xs >> np.uint64(k)) & np.uint64(1)).astype(np.uint8)
        out.append(int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little"))
    return out


def unpack(v: int, n: int) -> np.ndarray:
    b = np.frombuffer(v.to_bytes((n + 7) // 8, "little"), dtype=np.uint8)
    return np.unpackbits(b, bitorder="little")[:n]


class Sim:
    """Bitset simulator. hmr: 'zero' (all outcomes 0), 'one', or 'rand'."""

    def __init__(self, xs, nq=64, hmr="zero", seed=1):
        self.n = len(xs)
        self.xs = np.asarray(xs, dtype=np.uint64)
        self.full = (1 << self.n) - 1
        self.q = [0] * nq
        for k, v in enumerate(pack_inputs(xs)):
            self.q[k] = v
        self.b = {}
        self.stack = [self.full]
        self.phase = 0
        self.rot = []  # (op index, qubit, angle, bitset)
        self.hmr = hmr
        self.rng = random.Random(seed)
        self.tof_exec = 0  # sum over shots of executed toffolis
        self.i = 0

    def rnd(self):
        if self.hmr == "zero":
            return 0
        if self.hmr == "one":
            return self.full
        return self.rng.getrandbits(self.n)

    def run(self, ops, start=0, stop=None, hook=None):
        q, b = self.q, self.b
        stop = len(ops) if stop is None else stop
        for i in range(start, stop):
            name, qs, bit, cond, angle = ops[i]
            m = self.stack[-1]
            if cond is not None and name != "PUSH_CONDITION":
                m &= b.get(cond, 0)
            if name == "CX":
                q[qs[1]] ^= q[qs[0]] & m
            elif name == "CCX":
                q[qs[2]] ^= q[qs[0]] & q[qs[1]] & m
                self.tof_exec += m.bit_count()
            elif name == "X":
                q[qs[0]] ^= m
            elif name == "CZ":
                self.phase ^= q[qs[0]] & q[qs[1]] & m
            elif name == "Z":
                self.phase ^= q[qs[0]] & m
            elif name == "CCZ":
                self.phase ^= q[qs[0]] & q[qs[1]] & q[qs[2]] & m
                self.tof_exec += m.bit_count()
            elif name == "CSWAP":
                d = (q[qs[1]] ^ q[qs[2]]) & q[qs[0]] & m
                q[qs[1]] ^= d
                q[qs[2]] ^= d
                self.tof_exec += m.bit_count()
            elif name == "SWAP":
                d = (q[qs[0]] ^ q[qs[1]]) & m
                q[qs[0]] ^= d
                q[qs[1]] ^= d
            elif name == "HMR":
                r = self.rnd()
                b[bit] = (b.get(bit, 0) & ~m) | (r & m)
                self.phase ^= q[qs[0]] & r & m
                q[qs[0]] &= ~m
            elif name == "R":
                r = self.rnd() if self.hmr != "zero" else 0
                self.phase ^= q[qs[0]] & r & m
                q[qs[0]] &= ~m
            elif name == "Z_POW":
                self.rot.append((i, qs[0], angle, q[qs[0]] & m))
            elif name == "NEG":
                self.phase ^= m
            elif name == "BIT_INVERT":
                b[bit] = b.get(bit, 0) ^ m
            elif name == "BIT_STORE0":
                b[bit] = b.get(bit, 0) & ~m
            elif name == "BIT_STORE1":
                b[bit] = b.get(bit, 0) | m
            elif name == "PUSH_CONDITION":
                self.stack.append(self.stack[-1] & b.get(cond, 0))
            elif name == "POP_CONDITION":
                if len(self.stack) > 1:
                    self.stack.pop()
            else:
                raise ValueError(name)
            if hook is not None:
                hook(i, self)
        self.i = stop

    def phase_units(self):
        """Total phase per shot in units of 2^-64 turns (uint64, mod 1 turn), as the verifier computes it."""
        tot = np.zeros(self.n, dtype=np.uint64)
        for _, _, ang, bits in self.rot:
            u = half_turns_to_units(ang)
            tot += unpack(bits, self.n).astype(np.uint64) * np.uint64(u)
        tot += unpack(self.phase, self.n).astype(np.uint64) * np.uint64(1 << 63)
        return tot

    def value(self, qubits):
        """Integer value of a little-endian list of qubits per shot (uint64)."""
        v = np.zeros(self.n, dtype=np.uint64)
        for k, qq in enumerate(qubits):
            v |= unpack(self.q[qq], self.n).astype(np.uint64) << np.uint64(k)
        return v


def half_turns_to_units(half_turns) -> int:
    n, d = Fraction(half_turns).as_integer_ratio()
    return (((n << 64) // d + 1) >> 1) & ((1 << 64) - 1)


_REF_CACHE = {}


def reference(xs):
    """target (uint64 units of 2^-64 turn) and tol for each x, exactly as the verifier."""
    tg, tl = [], []
    for x in np.asarray(xs).tolist():
        r = math.isqrt((1 << 130) // x)
        tg.append(((r + 1) >> 1) & ((1 << 64) - 1))
        tl.append(min(r >> 25, 1 << 63))
    return np.array(tg, dtype=np.uint64), np.array(tl, dtype=np.uint64)


def reference_all(cache="<HOME>/qcc/ws/invsqrt/scratch/ref_all.npz"):
    import os
    if os.path.exists(cache):
        d = np.load(cache)
        return d["target"], d["tol"]
    t, l = reference(np.arange(1, MAX_X + 1, dtype=np.uint64))
    np.savez(cache, target=t, tol=l)
    return t, l


def check(sim, target=None, tol=None):
    """Return dict with pass flag, worst ratio, min margin, for a finished Sim."""
    if target is None:
        target, tol = reference(sim.xs)
    th = sim.phase_units()
    err = (th - target).view(np.int64)
    used = np.abs(err).view(np.uint64)
    ratio = used / np.maximum(tol, 1)
    margin = tol.astype(np.int64) - used.astype(np.int64)
    k = int(np.argmax(ratio))
    scratch_dirty = [i for i in range(22, len(sim.q)) if sim.q[i]]
    xr = sim.value(range(22))
    return dict(ok=bool((used <= tol).all()) and not scratch_dirty and bool((xr == sim.xs).all()),
                worst=float(ratio[k]), worst_x=int(sim.xs[k]), min_margin=int(margin.min()),
                min_margin_x=int(sim.xs[int(np.argmin(margin))]), dirty=scratch_dirty,
                x_restored=bool((xr == sim.xs).all()), err=err, tol=tol)
