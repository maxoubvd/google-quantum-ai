"""F3 (invsqrt_fw, 2026-11-13): exact bit-sliced simulator of a .kmb circuit over ALL inputs at once.

Our own code (reads the circuit file as data with tools/kmb.py; nothing of the public source is imported).
Every wire is a packed bit vector over the inputs x = 1..N (N = 3,133,452 by default), so one numpy op per
gate. Tracks: wire values, classical bits, the condition stack, the +-1 sign and the Z_POW angle sum (units of
2^-64 turn; the low word of every angle must be 0, checked), and, per op, which wires are non-zero.

HMR outcomes: mode "zero" (m = 0 on every shot: no phase kick, conditioned repairs inactive) or "rand" (seeded).

    python invsqrt_fw_sim.py <file.kmb> [N]     runs the check: x restored, scratch 0, phase within tolerance
"""
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import kmb as K  # noqa: E402

MAXX = 3133452
MASK64 = (1 << 64) - 1


def pack(bits_bool):
    """bool array (N,) -> packed uint64 words, little-endian bit order inside each word."""
    n = len(bits_bool)
    w = (n + 63) // 64
    pad = np.zeros(w * 64, dtype=np.uint8)
    pad[:n] = bits_bool
    return np.packbits(pad.reshape(-1, 8), axis=1, bitorder="little").reshape(-1).view(np.uint64).copy()


def unpack(words, n):
    return np.unpackbits(words.view(np.uint8), bitorder="little")[:n].astype(bool)


class Circ:
    """Decoded op list: per op (base opcode, cond bit or -1, qubits tuple, out bit or -1, angle units or None)."""

    def __init__(self, path):
        k = K.read(path, mmap=False)
        self.k = k
        self.nq, self.nb = k.num_qubits, k.num_bits
        self.regs = {name: ids.astype(int).tolist() for name, ids, isq in k.registers}
        ops = k.ops
        n = len(ops)
        nqv = ((ops & 0x30) >> 4).astype(int)
        i3 = np.cumsum(nqv == 3) - 1
        i2 = np.cumsum(nqv == 2) - 1
        i1 = np.cumsum(nqv == 1) - 1
        hasb = (ops & 0x40) != 0
        ib = np.cumsum(hasb) - 1
        hasc = (ops & 0x80) != 0
        ic = np.cumsum(hasc) - 1
        isz = (ops & 0x7F) == 0x13
        ia = np.cumsum(isz) - 1
        self.ops = []
        tag_i = 0
        for i in range(n):
            o = int(ops[i])
            base = 0x81 if o == 0x81 else o & 0x7F
            c = int(k.cond[ic[i]]) if hasc[i] else -1
            m = nqv[i]
            if m == 3:
                q = (int(k.q3[0][i3[i]]), int(k.q3[1][i3[i]]), int(k.q3[2][i3[i]]))
            elif m == 2:
                q = (int(k.q2[0][i2[i]]), int(k.q2[1][i2[i]]))
            elif m == 1:
                q = (int(k.q1[i1[i]]),)
            else:
                q = ()
            b = int(k.bit_t[ib[i]]) if hasb[i] else -1
            a = None
            if base == 0x13:
                w0, w1 = int(k.angles[ia[i]][0]), int(k.angles[ia[i]][1])
                a = (w0, w1)
            t = None
            if base == 0x06:
                t = k.tags[tag_i] if tag_i < len(k.tags) else None
                tag_i += 1
            self.ops.append((base, c, q, b, a, t))

    def name(self, i):
        base, c, q, b, a, t = self.ops[i]
        s = "PUSH_CONDITION" if base == 0x81 else K.NAMES.get(base, hex(base))
        return s


class Sim:
    def __init__(self, circ, n=MAXX, hmr="zero", seed=1, phase=True, track=True):
        self.c = circ
        self.n = n
        self.W = (n + 63) // 64
        self.V = np.zeros((circ.nq, self.W), dtype=np.uint64)
        self.B = np.zeros((max(circ.nb, 1), self.W), dtype=np.uint64)
        self.full = pack(np.ones(n, dtype=bool))
        self.stack = []
        self.hmr = hmr
        self.rng = np.random.default_rng(seed)
        self.phase = phase
        self.S = np.zeros(self.W, dtype=np.uint64)  # sign bit per input
        self.P = np.zeros(n, dtype=np.uint64) if phase else None  # angle sum, units 2^-64 turn
        self.track = track
        self.nz = np.zeros(circ.nq, dtype=bool)
        self.pos = 0
        self.tof = 0.0  # mean Toffoli-like executed per shot
        self.rot = []  # (op index, packed active bits) of every Z_POW executed
        self.live = []  # per op: number of non-zero wires after the op
        self.events = []  # (op, wire, now_nonzero)
        xs = np.arange(1, n + 1, dtype=np.int64)
        xr = circ.regs.get("x")
        assert xr is not None, circ.regs.keys()
        for j, w in enumerate(xr):
            self.V[w] = pack(((xs >> j) & 1).astype(bool))
            self.nz[w] = True
        self.xs = xs

    def mask(self, c):
        m = self.full
        for s in self.stack:
            m = m & s
        if c >= 0:
            m = m & self.B[c]
        return m

    def _upd(self, i, ws):
        if not self.track:
            return
        for w in ws:
            z = bool(self.V[w].any())
            if z != self.nz[w]:
                self.nz[w] = z
                self.events.append((i, w, z))

    def run_to(self, stop):
        V, B = self.V, self.B
        ops = self.c.ops
        while self.pos < stop:
            i = self.pos
            base, c, q, b, a, t = ops[i]
            if base == 0x06 or base == 0x07 or base == 0x17 or base == 0x47:
                pass
            elif base == 0x81:
                self.stack.append(B[c].copy())
            elif base == 0x01:
                self.stack.pop()
            else:
                m = self.mask(c)
                if base == 0x11:  # X
                    V[q[0]] ^= m
                    self._upd(i, q)
                elif base == 0x21:  # CX
                    V[q[1]] ^= V[q[0]] & m
                    self._upd(i, q[1:])
                elif base == 0x31:  # CCX
                    act = V[q[0]] & V[q[1]] & m
                    V[q[2]] ^= act
                    self.tof += self._mean(m)
                    self._upd(i, q[2:])
                elif base == 0x22:  # SWAP
                    d = (V[q[0]] ^ V[q[1]]) & m
                    V[q[0]] ^= d
                    V[q[1]] ^= d
                    self._upd(i, q)
                elif base == 0x32:  # CSWAP
                    d = (V[q[1]] ^ V[q[2]]) & V[q[0]] & m
                    V[q[1]] ^= d
                    V[q[2]] ^= d
                    self.tof += self._mean(m)
                    self._upd(i, q[1:])
                elif base == 0x10:  # Z
                    self.S ^= V[q[0]] & m
                elif base == 0x20:  # CZ
                    self.S ^= V[q[0]] & V[q[1]] & m
                elif base == 0x30:  # CCZ
                    self.S ^= V[q[0]] & V[q[1]] & V[q[2]] & m
                    self.tof += self._mean(m)
                elif base == 0x00:  # NEG
                    self.S ^= m
                elif base == 0x13:  # Z_POW
                    w0, w1 = a
                    act = V[q[0]] & m
                    self.rot.append((i, act.copy()))
                    if self.phase:
                        assert w0 == 0, ("angle below 2^-64 turn", i, w0)
                        ab = unpack(act, self.n)
                        self.P[ab] += np.uint64(w1)
                elif base == 0x12:  # R
                    V[q[0]] &= ~m
                    self._upd(i, q)
                elif base == 0x50:  # HMR
                    if self.hmr == "zero":
                        mb = np.zeros(self.W, dtype=np.uint64)
                    else:
                        mb = pack(self.rng.integers(0, 2, self.n, dtype=np.uint8).astype(bool))
                    mb &= m
                    self.S ^= mb & V[q[0]]
                    B[b] = (B[b] & ~m) | mb
                    V[q[0]] &= ~m
                    self._upd(i, q)
                elif base == 0x40:
                    B[b] &= ~m
                elif base == 0x41:
                    B[b] |= m
                elif base == 0x42:
                    B[b] ^= m
                else:
                    raise ValueError(f"op {i}: unhandled opcode {base:#x}")
            if self.track:
                self.live.append(int(self.nz.sum()))
            self.pos += 1
        return self

    def _mean(self, m):
        if m is self.full:
            return 1.0
        return float(np.unpackbits(m.view(np.uint8)).sum()) / self.n

    def phase_units(self):
        s = unpack(self.S, self.n).astype(np.uint64) << np.uint64(63)
        return self.P + s


def reference(n=MAXX):
    """target mod 1 and tol(x), uint64 units of 2^-64 turn, exactly as the official verifier (beta = 1, 24 bits)."""
    import math
    num = 1 << 130
    shift = 24 + 1
    rs = [math.isqrt(num // x) for x in range(1, n + 1)]
    target = np.array([((r + 1) >> 1) & MASK64 for r in rs], dtype=np.uint64)
    tol = np.array([min(r >> shift, 1 << 63) for r in rs], dtype=np.uint64)
    return target, tol


def check(path, n=MAXX, hmr="zero", seed=1):
    t0 = time.time()
    c = Circ(path)
    print(f"{path}: Q {c.nq}, bits {c.nb}, ops {len(c.ops)}, registers {[ (k, len(v)) for k, v in c.regs.items()]}"
          f" (decoded in {time.time() - t0:.1f} s), HMR outcomes {hmr}", flush=True)
    s = Sim(c, n=n, hmr=hmr, seed=seed)
    s.run_to(len(c.ops))
    xr = c.regs["x"]
    ok_x = all(np.array_equal(s.V[w], Sim(c, n=n, phase=False, track=False).V[w]) for w in xr)
    scratch = [w for w in range(c.nq) if w not in xr]
    dirty = [w for w in scratch if s.V[w].any()]
    target, tol = reference(n)
    err = (s.phase_units() - target).view(np.int64)
    used = np.abs(err).astype(np.uint64)
    margin = tol.astype(np.int64) - used.astype(np.int64)
    print(f"simulated in {time.time() - t0:.1f} s: x restored {ok_x}, dirty scratch {dirty}, mean Toffolis "
          f"{s.tof:.4f}, Z_POW executed {len(s.rot)}, max live {max(s.live)}; failing {(margin < 0).sum()}, "
          f"min margin {int(margin.min())} at x = {int(np.argmin(margin)) + 1}", flush=True)
    return c, s, err, tol


if __name__ == "__main__":
    check(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else MAXX,
          sys.argv[3] if len(sys.argv) > 3 else "zero", int(sys.argv[4]) if len(sys.argv) > 4 else 1)
