#!/usr/bin/env python3
"""invsqrt_fc_kit.py -- FC gate kit (own code): circuit builder with LIFO wire allocation and peak tracking,
a bit-sliced exact simulator over a batch of inputs (every wire = numpy bool array; HMR outcomes random per
shot; phase sign and rotation angle tracked per shot), and a kickmix .kmx text writer.

Gates: X, CX, CCX, SWAP, CSWAP, Z, CZ, CCZ, ZP (Z_POW, angle as an integer numerator over 2^ANG_BITS turns),
HMR (returns a classical bit), and conditioned versions of X / CX / Z / CZ / ZP (cond = classical bit id).
Counted Toffoli-like: CCX, CCZ, CSWAP (also when conditioned; this kit never conditions them).
Temporary AND (Gidney): and_compute(a, b) -> t  (1 CCX);  and_erase(a, b, t) (HMR + CZ(a, b) if m: 0 CCX).
"""
import numpy as np

ANG_BITS = 64          # rotation angles are integers over 2^64 turns (uint64, wraps)


class Circ:
    def __init__(self, n_fixed=0):
        self.ops = []
        self.n_alloc = n_fixed
        self.free = []
        self.next = n_fixed
        self.live = n_fixed
        self.peak = n_fixed
        self.nbits = 0
        self.tags = []          # (op index, name) stage marks
        self.stage_peaks = {}   # stage name -> peak live wires inside it

    # ---- wires
    def alloc(self, k=1):
        out = []
        for _ in range(k):
            if self.free:
                q = self.free.pop()
            else:
                q = self.next
                self.next += 1
            out.append(q)
            self.live += 1
        self.peak = max(self.peak, self.live)
        if self.tags:
            nm = self.tags[-1][1]
            self.stage_peaks[nm] = max(self.stage_peaks.get(nm, 0), self.live)
        return out

    def free_(self, qs):
        for q in reversed(list(qs)):
            self.free.append(q)
            self.live -= 1

    def bit(self):
        b = self.nbits
        self.nbits += 1
        return b

    def tag(self, name):
        self.tags.append((len(self.ops), name, self.live))
        self.stage_peaks[name] = max(self.stage_peaks.get(name, 0), self.live)

    # ---- gates
    def x(self, q, cond=None):
        self.ops.append(("X", q, cond))

    def cx(self, c, t, cond=None):
        assert c != t
        self.ops.append(("CX", c, t, cond))

    def ccx(self, a, b, t):
        assert len({a, b, t}) == 3
        self.ops.append(("CCX", a, b, t))

    def swap(self, a, b):
        assert a != b
        self.ops.append(("SWAP", a, b))

    def cswap(self, c, a, b):
        assert len({c, a, b}) == 3
        self.ops.append(("CSWAP", c, a, b))

    def z(self, q, cond=None):
        self.ops.append(("Z", q, cond))

    def cz(self, a, b, cond=None):
        assert a != b
        self.ops.append(("CZ", a, b, cond))

    def ccz(self, a, b, c):
        assert len({a, b, c}) == 3
        self.ops.append(("CCZ", a, b, c))

    def zp(self, q, num, cond=None):
        """Z_POW with angle num / 2^ANG_BITS turns."""
        self.ops.append(("ZP", q, int(num) % (1 << ANG_BITS), cond))
        return len(self.ops) - 1

    def neg(self, cond=None):
        self.ops.append(("NEG", cond))

    def hmr(self, q):
        b = self.bit()
        self.ops.append(("HMR", q, b))
        return b

    # ---- temporary AND (logical ops: ANDC = CCX onto a zero target; ANDE = HMR + CZ(a, b) if m; the
    # mirror of a stage swaps them, so every self-contained stage has an automatic inverse)
    def and_compute(self, a, b, t=None):
        if t is None:
            t = self.alloc(1)[0]
        assert len({a, b, t}) == 3
        self.ops.append(("ANDC", a, b, t))
        return t

    def and_erase(self, a, b, t, free=True):
        assert len({a, b, t}) == 3
        self.ops.append(("ANDE", a, b, t))
        if free:
            self.free_([t])

    # ---- conditioned temporary AND (repairs of the narrow undo): CANDC = CCX if bit g onto a zero target;
    # CANDE = erase of t = g a b (HMR, then CZ(a, b) if g and the outcome)
    def cand_compute(self, a, b, t, g):
        assert len({a, b, t}) == 3
        self.ops.append(("CANDC", a, b, t, g))

    def cand_erase(self, a, b, t, g):
        assert len({a, b, t}) == 3
        self.ops.append(("CANDE", a, b, t, g))

    def toffolis(self, ops=None):
        """static count (a conditioned AND counts 1)."""
        return sum(1 for o in (self.ops if ops is None else ops) if o[0] in ("CCX", "CCZ", "CSWAP", "ANDC", "CANDC"))

    def toffolis_exp(self, ops=None):
        """expected count per shot: a CANDC fires with probability 1/2 (its bit is an X-measurement outcome)."""
        ops = self.ops if ops is None else ops
        return (sum(1 for o in ops if o[0] in ("CCX", "CCZ", "CSWAP", "ANDC"))
                + 0.5 * sum(1 for o in ops if o[0] == "CANDC"))

    @staticmethod
    def mirror_blocks(ops, blocks):
        """inverse of a forward op list in which narrow adders sit between ("BBEG", id) and ("BEND", id): every
        such block is replaced by its explicit undo blocks[id]; the rest is mirrored."""
        out = []
        i = len(ops) - 1
        seg = []
        while i >= 0:
            o = ops[i]
            if o[0] == "BEND":
                out.extend(Circ.mirror(list(reversed(seg))))
                seg = []
                bid = o[1]
                out.extend(blocks[bid])
                while ops[i][0] != "BBEG" or ops[i][1] != bid:
                    i -= 1
                i -= 1
                continue
            seg.append(o)
            i -= 1
        out.extend(Circ.mirror(list(reversed(seg))))
        return out

    @staticmethod
    def mirror(ops):
        """inverse of a self-contained op list: reversed, ANDC <-> ANDE (no HMR / ZP / NEG allowed inside)."""
        out = []
        for o in reversed(ops):
            g = o[0]
            if g in ("BBEG", "BEND"):
                continue
            if g == "ANDC":
                out.append(("ANDE",) + o[1:])
            elif g == "ANDE":
                out.append(("ANDC",) + o[1:])
            elif g in ("X", "CX", "CCX", "SWAP", "CSWAP", "Z", "CZ", "CCZ"):
                if g in ("X", "CX", "Z", "CZ") and o[-1] is not None:
                    raise ValueError("conditioned op in a mirrored stage")
                out.append(o)
            else:
                raise ValueError("op %s cannot be mirrored" % g)
        return out

    def stage_counts(self):
        """Toffoli-like count between tags."""
        out = []
        marks = self.tags + [(len(self.ops), "END", self.live)]
        for (i0, nm, _), (i1, _, _) in zip(marks[:-1], marks[1:]):
            out.append((nm, sum(1 for o in self.ops[i0:i1] if o[0] in ("CCX", "CCZ", "CSWAP", "ANDC"))))
        return out


class Sim:
    """bit-sliced simulator: n shots, wires as bool arrays."""

    def __init__(self, n, nq, seed=1):
        self.n = n
        self.q = np.zeros((nq, n), dtype=bool)
        self.rng = np.random.default_rng(seed)
        self.sign = np.zeros(n, dtype=bool)
        self.ang = np.zeros(n, dtype=np.uint64)  # numerators over 2^64 turns (uint64 wraps)
        self.bits = {}
        self.check = True
        self.fire = None          # set to [] to record the firing pattern of every ZP in order

    def set_word(self, wires, vals):
        vals = np.asarray(vals, dtype=np.int64)
        for k, w in enumerate(wires):
            self.q[w] = ((vals >> k) & 1).astype(bool)

    def get_word(self, wires, signed=False):
        v = np.zeros(self.n, dtype=np.int64)
        for k, w in enumerate(wires):
            v |= self.q[w].astype(np.int64) << k
        if signed and len(wires):
            top = 1 << (len(wires) - 1)
            v = np.where(v >= top, v - 2 * top, v)
        return v

    def run(self, ops, start=0, stop=None):
        q = self.q
        for o in ops[start:stop]:
            g = o[0]
            if g == "X":
                _, a, c = o
                if c is None:
                    q[a] ^= True
                else:
                    q[a] ^= self.bits[c]
            elif g == "CX":
                _, a, t, c = o
                if c is None:
                    q[t] ^= q[a]
                else:
                    q[t] ^= q[a] & self.bits[c]
            elif g == "CCX":
                _, a, b, t = o
                q[t] ^= q[a] & q[b]
            elif g == "SWAP":
                _, a, b = o
                tmp = q[a].copy()
                q[a] = q[b]
                q[b] = tmp
            elif g == "CSWAP":
                _, c, a, b = o
                d = (q[a] ^ q[b]) & q[c]
                q[a] ^= d
                q[b] ^= d
            elif g == "Z":
                _, a, c = o
                self.sign ^= q[a] if c is None else (q[a] & self.bits[c])
            elif g == "CZ":
                _, a, b, c = o
                v = q[a] & q[b]
                self.sign ^= v if c is None else (v & self.bits[c])
            elif g == "CCZ":
                _, a, b, cc = o
                self.sign ^= q[a] & q[b] & q[cc]
            elif g == "ZP":
                _, a, num, c = o
                on = q[a] if c is None else (q[a] & self.bits[c])
                self.ang += np.where(on, np.uint64(num), np.uint64(0))
                if self.fire is not None:
                    self.fire.append(on.copy())
            elif g == "NEG":
                _, c = o
                self.sign ^= True if c is None else self.bits[c]
            elif g == "HMR":
                _, a, b = o
                m = self.rng.integers(0, 2, self.n).astype(bool)
                self.bits[b] = m
                self.sign ^= q[a] & m
                q[a] = False
            elif g == "ANDC":
                _, a, b, t = o
                if self.check and q[t].any():
                    raise AssertionError("ANDC onto a nonzero target q%d" % t)
                q[t] ^= q[a] & q[b]
            elif g == "ANDE":
                _, a, b, t = o
                if self.check and (q[t] != (q[a] & q[b])).any():
                    raise AssertionError("ANDE: target q%d is not the AND of q%d q%d" % (t, a, b))
                q[t] = False          # HMR + CZ(a, b) if m: the phase cancels exactly when t = a b
            elif g == "CANDC":
                _, a, b, t, gb = o
                if self.check and q[t].any():
                    raise AssertionError("CANDC onto a nonzero target q%d" % t)
                q[t] ^= q[a] & q[b] & self.bits[gb]
            elif g == "CANDE":
                _, a, b, t, gb = o
                if self.check and (q[t] != (q[a] & q[b] & self.bits[gb])).any():
                    raise AssertionError("CANDE: target q%d is not g AND q%d q%d" % (t, a, b))
                q[t] = False
            elif g in ("BBEG", "BEND", "CBEG", "CEND"):
                pass
            else:
                raise ValueError(g)


KMX_CBLOCKS = False   # write each marked repair span as one PUSH_CONDITION block (P1's tools see it)
KMX_CSWAP_SPLIT = False   # FC3 (18:45Z): write CSWAP as CX CCX CX (set by g8's build_full from cfg "kmx_cswap_split")


def kmx_text(circ, x_wires, qmax=None):
    """kickmix text: register r0 "x" = x_wires (LSB first); ops; Z_POW in half-turns as exact decimals."""
    from fractions import Fraction
    lines = []
    for w in x_wires:
        lines.append(f"APPEND_TO_REGISTER q{w} r0")
    lines.append('REGISTER r0 "x"')

    def fmt_half_turns(num):
        fr = Fraction(2 * num, 1 << ANG_BITS) % 2
        # exact decimal of a dyadic fraction
        den = fr.denominator
        k = den.bit_length() - 1
        assert den == 1 << k
        intpart = fr.numerator // den
        rem = fr.numerator % den
        digits = (rem * 10 ** k) // den if k else 0
        s = f"{intpart}." + (str(digits).rjust(k, "0") if k else "0")
        return s

    bitc = [circ.nbits]
    blk = [None]                     # KMX_CBLOCKS: the bit of the open repair block (its ops are written plain)

    def cnd(c):
        if c is None or (blk[0] is not None and c == blk[0]):
            return ""
        return f" if b{c}"
    for o in circ.ops:
        g = o[0]
        if g == "CBEG":
            if KMX_CBLOCKS:
                lines.append(f"PUSH_CONDITION if b{o[1]}")
                blk[0] = o[1]
            continue
        if g == "CEND":
            if KMX_CBLOCKS:
                lines.append("POP_CONDITION")
                blk[0] = None
            continue
        if g == "X":
            _, a, c = o
            lines.append(f"X q{a}" + cnd(c))
        elif g == "CX":
            _, a, t, c = o
            lines.append(f"CX q{a} q{t}" + cnd(c))
        elif g == "CCX":
            lines.append(f"CCX q{o[1]} q{o[2]} q{o[3]}")
        elif g == "SWAP":
            lines.append(f"SWAP q{o[1]} q{o[2]}")
        elif g == "CSWAP":
            if KMX_CSWAP_SPLIT:          # FC3: CX(b, a) CCX(c, a, b) CX(b, a): same count; P1's pass can rewrite the CCX
                lines += [f"CX q{o[3]} q{o[2]}", f"CCX q{o[1]} q{o[2]} q{o[3]}", f"CX q{o[3]} q{o[2]}"]
            else:
                lines.append(f"CSWAP q{o[1]} q{o[2]} q{o[3]}")
        elif g == "Z":
            _, a, c = o
            lines.append(f"Z q{a}" + cnd(c))
        elif g == "CZ":
            _, a, b, c = o
            lines.append(f"CZ q{a} q{b}" + cnd(c))
        elif g == "CCZ":
            lines.append(f"CCZ q{o[1]} q{o[2]} q{o[3]}")
        elif g == "ZP":
            _, a, num, c = o
            lines.append(f"Z_POW q{a} {fmt_half_turns(num)}" + (f" if b{c}" if c is not None else ""))
        elif g == "NEG":
            _, c = o
            lines.append("NEG" + cnd(c))
        elif g == "HMR":
            lines.append(f"HMR q{o[1]} b{o[2]}")
        elif g == "ANDC":
            # temporary AND onto a zero target: one CCX
            lines.append(f"CCX q{o[1]} q{o[2]} q{o[3]}")
        elif g == "ANDE":
            # measurement-based erase of t = a b: X-basis measure, phase fixed by CZ(a, b) on outcome 1
            b = bitc[0]
            bitc[0] += 1
            lines.append(f"HMR q{o[3]} b{b}")
            lines.append(f"CZ q{o[1]} q{o[2]} if b{b}")
        elif g == "CANDC":
            lines.append(f"CCX q{o[1]} q{o[2]} q{o[3]}" + cnd(o[4]))
        elif g == "CANDE":
            # t = g a b: the X-measurement kicks (-1)^(m g a b); repaired by CZ(a, b) when g and m
            b = bitc[0]
            bitc[0] += 1
            lines.append(f"HMR q{o[3]} b{b}")
            if blk[0] is not None and o[4] == blk[0]:
                lines.append(f"CZ q{o[1]} q{o[2]} if b{b}")
            else:
                lines.append(f"PUSH_CONDITION if b{o[4]}")
                lines.append(f"CZ q{o[1]} q{o[2]} if b{b}")
                lines.append("POP_CONDITION")
        elif g in ("BBEG", "BEND"):
            pass
        else:
            raise ValueError("kmx_text: unhandled op %s" % g)
    return "\n".join(lines) + "\n"
