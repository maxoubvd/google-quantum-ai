"""invsqrt_gen: circuit builder and exact bitset simulator (own code, workstream invsqrt_gen).

Builder: a flat op list (name, qubits, bit, cond, angle) with a smallest-free-id qubit allocator, so that
num_qubits (largest id + 1) equals the peak number of simultaneously allocated wires. Pseudo-ops starting with
"@" are bookkeeping for the simulator and are not written to the circuit file:
  ("@free", (q,))   the wire must be |0> for every input here (checked by the simulator)
  ("@stage", name)  stage boundary for the executed-Toffoli table
Simulator: one Python int per wire, bit k = value for input k (same technique as tools/invsqrt_lib.py of the
invsqrt workstream, rewritten here with stage accounting and wire checks). Semantics follow
research/05_kickmix_and_techniques.md A1: HMR q b: b := random, phase ^= q & b, q := 0.
"""
from __future__ import annotations

import heapq
import random
from fractions import Fraction

import numpy as np

TOF = ("CCX", "CCZ", "CSWAP")


class Circ:
    def __init__(self, n_io=22):
        self.n_io = n_io
        self.ops = []
        self._free = []
        self._next = n_io
        self.next_bit = 0
        self.peak = n_io
        self.live = set(range(n_io))
        self._free_bits = []
        self.meta = {}
        self.lendable = set()
        self.stage_peaks = []            # [stage name, largest number of live wires inside the stage]

    def _count(self):
        n = len(self.live)
        self.peak = max(self.peak, n)
        if self.stage_peaks and n > self.stage_peaks[-1][1]:
            self.stage_peaks[-1][1] = n

    # ---- wires
    def alloc(self):
        q = heapq.heappop(self._free) if self._free else self._bump()
        self.live.add(q)
        self._count()
        return q

    def _bump(self):
        self._next += 1
        return self._next - 1

    def allocn(self, n):
        return [self.alloc() for _ in range(n)]

    def free(self, q):
        assert q in self.live and (q >= self.n_io or q in self.lendable), q
        self.live.discard(q)
        heapq.heappush(self._free, q)
        self.ops.append(("@free", (q,), None, None, None))

    def lend(self, q):
        """An I/O wire that is |0> from here on may be handed out by alloc() until reclaim(q)."""
        assert q in self.live and q < self.n_io, q
        self.lendable.add(q)
        self.free(q)

    def reclaim(self, q):
        """Take a lent I/O wire back (it must be free, hence |0>)."""
        assert q in self._free, q
        self._free.remove(q)
        heapq.heapify(self._free)
        self.live.add(q)
        self._count()
        self.lendable.discard(q)

    def alloc_high(self):
        """Like alloc(), but never returns a lent I/O wire."""
        cand = sorted(w for w in self._free if w >= self.n_io)
        if cand:
            self._free.remove(cand[0])
            heapq.heapify(self._free)
            q = cand[0]
        else:
            q = self._bump()
        self.live.add(q)
        self._count()
        return q

    def is_free(self, q):
        return q in self._free

    def bit(self):
        if self._free_bits:
            return heapq.heappop(self._free_bits)
        self.next_bit += 1
        return self.next_bit - 1

    def bit_free(self, b):
        """The classical bit is not read any more: its id may be reused by a later measurement."""
        heapq.heappush(self._free_bits, b)

    def stage(self, name):
        self.ops.append(("@stage", name, None, None, len(self.live)))
        self.stage_peaks.append([name, len(self.live)])

    # ---- gates
    def _op(self, name, qs, bit=None, cond=None, angle=None):
        assert len(set(qs)) == len(qs), (name, qs)
        self.ops.append((name, tuple(qs), bit, cond, angle))

    def x(self, q, cond=None):
        self._op("X", (q,), cond=cond)

    def cx(self, c, t, cond=None):
        self._op("CX", (c, t), cond=cond)

    def ccx(self, a, b, t, cond=None):
        self._op("CCX", (a, b, t), cond=cond)

    def z(self, q, cond=None):
        self._op("Z", (q,), cond=cond)

    def cz(self, a, b, cond=None):
        self._op("CZ", (a, b), cond=cond)

    def ccz(self, a, b, c, cond=None):
        self._op("CCZ", (a, b, c), cond=cond)

    def swap(self, a, b, cond=None):
        self._op("SWAP", (a, b), cond=cond)

    def cswap(self, c, a, b, cond=None):
        self._op("CSWAP", (c, a, b), cond=cond)

    def zpow(self, q, half_turns: Fraction):
        self._op("Z_POW", (q,), angle=Fraction(half_turns))

    def hmr(self, q):
        """Measure q in the X basis and reset it; returns the outcome bit (fix-up is the caller's job)."""
        b = self.bit()
        self._op("HMR", (q,), bit=b)
        return b

    def neg(self, cond=None):
        self._op("NEG", (), cond=cond)

    def push(self, b):
        """Everything until the matching pop() is executed only in shots where classical bit b is 1."""
        self.ops.append(("PUSH_CONDITION", (), None, b, None))

    def pop(self):
        self.ops.append(("POP_CONDITION", (), None, None, None))

    # ---- composites
    def and_new(self, a, b):
        """Fresh wire holding a AND b (1 Toffoli)."""
        t = self.alloc()
        self.ccx(a, b, t)
        return t

    def and_del(self, t, a, b):
        """Erase t = a AND b by measurement (0 Toffolis) and release the wire."""
        m = self.hmr(t)
        self.cz(a, b, cond=m)
        self.bit_free(m)
        self.free(t)

    # ---- output
    def real_ops(self):
        return [o for o in self.ops if not o[0].startswith("@")]

    def num_qubits(self):
        return 1 + max(max(o[1]) for o in self.real_ops() if o[1])

    def kmx_lines(self):
        out = [f"APPEND_TO_REGISTER q{k} r0" for k in range(self.n_io)]
        out.append('REGISTER r0 "x"')
        out += [op_text(o) for o in self.real_ops()]
        return out

    def write_kmx(self, path):
        with open(path, "w") as fh:
            fh.write("\n".join(self.kmx_lines()) + "\n")

    def static_counts(self):
        c = {}
        for o in self.real_ops():
            c[o[0]] = c.get(o[0], 0) + 1
        return c


def frac_text(a: Fraction) -> str:
    a = a % 2
    n, d = a.numerator, a.denominator
    assert d & (d - 1) == 0, "angle must be dyadic"
    k = d.bit_length() - 1
    if k == 0:
        return f"{n}.0"
    digits = str(n * 5 ** k).rjust(k + 1, "0")
    return digits[:-k] + "." + digits[-k:]


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


def kmx_to_kmb(kmx, kmb):
    import sys
    if "<HOME>/qcc/data/kickmix" not in sys.path:
        sys.path.insert(0, "<HOME>/qcc/data/kickmix")
    import kickmix as km
    c = km.Circuit.from_file(kmx, format="kmx")
    c.to_file(kmb, format="kmb")
    return c


# ------------------------------------------------------------------ simulator

def pack_bits(vals, k):
    bits = ((vals >> np.uint64(k)) & np.uint64(1)).astype(np.uint8)
    return int.from_bytes(np.packbits(bits, bitorder="little").tobytes(), "little")


def unpack(v: int, n: int) -> np.ndarray:
    b = np.frombuffer(v.to_bytes((n + 7) // 8, "little"), dtype=np.uint8)
    return np.unpackbits(b, bitorder="little")[:n]


def half_turns_to_units(half_turns) -> int:
    n, d = Fraction(half_turns).as_integer_ratio()
    return (((n << 64) // d + 1) >> 1) & ((1 << 64) - 1)


class Sim:
    """hmr: 'zero' (all outcomes 0), 'one', or 'rand'."""

    def __init__(self, xs, nq, hmr="rand", seed=1, width=22):
        self.n = len(xs)
        self.xs = np.asarray(xs, dtype=np.uint64)
        self.full = (1 << self.n) - 1
        self.q = [0] * nq
        for k in range(width):
            self.q[k] = pack_bits(self.xs, k)
        self.b = {}
        self.stack = [self.full]
        self.phase = 0
        self.rot = []            # (angle, bitset)
        self.hmr = hmr
        self.rng = random.Random(seed)
        self.tof = 0
        self.stages = []         # (name, tof at entry, live wires at entry)
        self.bad_free = []

    def rnd(self):
        if self.hmr == "zero":
            return 0
        if self.hmr == "one":
            return self.full
        return self.rng.getrandbits(self.n)

    def run(self, ops, start=0, stop=None):
        q, b = self.q, self.b
        stop = len(ops) if stop is None else stop
        for i in range(start, stop):
            name, qs, bit, cond, angle = ops[i]
            if name[0] == "@":
                if name == "@free":
                    if q[qs[0]]:
                        self.bad_free.append((i, qs[0]))
                elif name == "@stage":
                    self.stages.append((qs, self.tof, angle))
                continue
            m = self.stack[-1]
            if cond is not None and name != "PUSH_CONDITION":
                m &= b.get(cond, 0)
            if name == "CX":
                q[qs[1]] ^= q[qs[0]] & m
            elif name == "CCX":
                q[qs[2]] ^= q[qs[0]] & q[qs[1]] & m
                self.tof += m.bit_count()
            elif name == "X":
                q[qs[0]] ^= m
            elif name == "CZ":
                self.phase ^= q[qs[0]] & q[qs[1]] & m
            elif name == "Z":
                self.phase ^= q[qs[0]] & m
            elif name == "CCZ":
                self.phase ^= q[qs[0]] & q[qs[1]] & q[qs[2]] & m
                self.tof += m.bit_count()
            elif name == "CSWAP":
                d = (q[qs[1]] ^ q[qs[2]]) & q[qs[0]] & m
                q[qs[1]] ^= d
                q[qs[2]] ^= d
                self.tof += m.bit_count()
            elif name == "SWAP":
                d = (q[qs[0]] ^ q[qs[1]]) & m
                q[qs[0]] ^= d
                q[qs[1]] ^= d
            elif name == "HMR":
                r = self.rnd()
                b[bit] = (b.get(bit, 0) & ~m) | (r & m)
                self.phase ^= q[qs[0]] & r & m
                q[qs[0]] &= ~m
            elif name == "Z_POW":
                self.rot.append((angle, q[qs[0]] & m))
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

    def value(self, wires):
        v = np.zeros(self.n, dtype=np.uint64)
        for k, w in enumerate(wires):
            v |= unpack(self.q[w], self.n).astype(np.uint64) << np.uint64(k)
        return v

    def phase_units(self):
        tot = np.zeros(self.n, dtype=np.uint64)
        for ang, bits in self.rot:
            if bits:
                tot += unpack(bits, self.n).astype(np.uint64) * np.uint64(half_turns_to_units(ang))
        tot += unpack(self.phase, self.n).astype(np.uint64) * np.uint64(1 << 63)
        return tot

    def stage_table(self):
        out = []
        for k, (name, tof, live) in enumerate(self.stages):
            nxt = self.stages[k + 1][1] if k + 1 < len(self.stages) else self.tof
            out.append((name, (nxt - tof) / self.n, live))
        return out
