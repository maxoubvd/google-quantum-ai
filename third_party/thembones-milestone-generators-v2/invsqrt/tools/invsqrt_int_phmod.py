"""invsqrt_int_phmod (workstream invsqrt_int, builder B1, 2026-11-12): copy of tools/invsqrt_zpol_phmod.py as pinned
in tools/invsqrt_gen_pin/ph6/ (sha256 2f912832dc9852d7), gate for gate the same at its defaults, with ONE addition
for DESIGN3 item I2 (zone flag as a lookup column, S3-INV-E 1.3):
  unnorm(..., zext=<wire>, hook2=<callable>): the backward zone rotation takes the zone flag z from the caller's
  wire `zext` (looked up by the caller, e.g. as one more output of its address lookup) instead of making it with
  5 ANDs; after the rotation `hook2()` is called, where the caller measures and repairs z (its helper chain is
  still alive). Forward (norm) is unchanged. Without zext nothing changes.
  unnorm(..., zt4=<wire>, hook2=<callable>): the caller holds t4 = NOT a5 AND NOT a4 AND NOT (a3 AND a2) AND NOT p
  (the flag's own AND tree, made by the caller as the first helpers of its unlookup chain); the backward flag is
  then z = a6 AND NOT t4, 1 AND instead of 5; hook2() after the rotation lets the caller erase t4 and its parts.

Header of the copy follows.
invsqrt_zpol_phmod: copy of the phase-delivery module v4d (tools/invsqrt_ph_mod.py, sha256 13b6ebd7f14bc348,
workstream invsqrt_ph; copied 2026-11-11 23:50Z by workstream invsqrt_zpol, the original is not edited) with two
options for the zone format, everything else gate for gate as v4d:
  zb  = "z9c" (default): zone flag z = a6 AND (p OR a5 OR a4 OR (a3 AND a2)), 5 ANDs each way (v4d);
        "z9d": zone 1 = a7 >= 80 for even exponents, a7 >= 64 for odd ones: z = a6 AND (p OR a5 OR a4), 3 ANDs.
  zsb = "pos" (default): zone special bits as in v4d, [e < 8] AND xq13 AND B_q at 2 Toffolis per position;
        "gate": the gate product t13 = [e < 8] AND xq13 is made once at the first zone position and kept to the
        last; each bit is t13 AND B_q (1 Toffoli, 0 if B_q = 1) and is erased against t13 at its own step:
        1 + (positions) Toffolis, the same pool targets as v4d, one wire held across those steps;
        "chain": the generator owner's carry chain c_0 = t13 AND B_top, c_k = c_(k-1) AND N_(k-1) (N = hi wires):
        1 + (positions) Toffolis and ONE pool target, but c_0..c_(k-1) stay alive until the last position.
Interface otherwise unchanged (PhaseDelivery(c, xq, neg=True, half=True, zone=True, zb=..., zsb=...)). The
imports of invsqrt_ph_net / _stream / _sx are the v4d files, unchanged. Log: research/ws_invsqrt_zpol.md.

Original header of v4d follows.
invsqrt_ph: phase-delivery module of the invsqrt generator (own code, workstream invsqrt_ph).

Emitters written against the generator's circuit builder (tools/invsqrt_gen_circ.py, imported by the caller and
never edited): normaliser NORM / UNNORM, the exponent code on the low wires of the normalised input, the hi
stream PH1 with the small-x path, the lo stream PH2. Interface note: research/ws_invsqrt_ph.md.

Notation: e = floor(log2 x), s = 21 - e, j = e >> 1, d = 10 - j (stream delay), p = e & 1, xn = x << s.
Tier controls a = d & 1, b = (d >> 1) & 1, E14 = [e < 14], E6 = [e < 6]; F = [e < 20], G = [e < 18].
Exponent code between NORM and UNNORM: wire F, and XORed onto the wires of xn bits 0..3: b, E14, a AND G, E6
(each is 0 whenever its wire can hold a data bit). Parity p sits on wire 21 after NORM (xn bit 21 is cleared).

Public ideas used, credited: thermometer of the leading one, reflection swaps with CNOT-loaded controls, code
bits parked on zero wires of xn plus one flag wire, half thermometer for the inverse, small-x path through the
hi stream: best public circuit, manwithacat "HW-25" (gengsr + maoxin123), structure read in
research/ws_invsqrt_decode.md part 3 (no gate list, table or constant of that file is used here). Temporary AND
and measurement-based erasure: Gidney 1709.06648, 1905.07682. Own: the swap networks (invsqrt_ph_netsearch.py),
the code layout and its decode, free flag erasure by quadratic forms over kept wires, the stream planner and
repair (invsqrt_ph_stream.py), the special-bit synthesiser (invsqrt_ph_sx.py).
"""
from __future__ import annotations

import invsqrt_ph_net as NET
import invsqrt_ph_stream as ST
import invsqrt_ph_sx as SX


def d_of(e):
    return 10 - (e >> 1)


E_FUN = {
    "F": lambda e: int(e < 20),
    "G": lambda e: int(e < 18),
    "P": lambda e: e & 1,
    "B": lambda e: (d_of(e) >> 1) & 1,
    "E14": lambda e: int(e < 14),
    "C2": lambda e: (d_of(e) & 1) & int(e < 18),
    "C3": lambda e: int(e < 6),
    "BE": lambda e: ((d_of(e) >> 1) & 1) & int(e < 14),
}
KEPT = ("F", "P", "B", "E14", "C2", "C3", "G")
CELLS = {"A0": 0, "B0": 1, "B1": 2, "R0": 4, "R1": 5, "R2": 6, "R3": 7, "X0": 8, "X1": 9, "X2": 10, "X3": 11}


class PhaseDelivery:
    def __init__(self, c, xq, angles_hi=None, angles_lo=None, top_pos=35, lo_min_pos=10, lo_sign=1, neg=False,
                 half=False, zone=False, zb="z9c", zsb="pos", sx_tries=60):
        """angles_hi / angles_lo: dict position -> Fraction of a half turn for the rotations of the hi / lo stream
        (default 2^-position, or -2^-position for the lo stream if lo_sign = -1: phase word = hi - lo).
        lo_min_pos: highest-weight position that has a rotation in the lo stream (10; 11 if the tables keep the
        top lo bit clear for every x < 256, which saves one rotation). Lo bits of x < 64 above it go through
        the small-x path of the hi stream."""
        from fractions import Fraction
        # neg: "negated word" form (40 rotations). The arithmetic computes X = 2^28 - (phase word) in the
        # generator's reg units: hi word N = 7 wires (0..127), lo as before; the module delivers -X: every angle
        # negative, constant 1 at the positions above the hi word (sign extension, free), no rotation at position
        # 10 in the lo stream (X is smallest at t = 0, so the top lo wire is 0 for every x in 64..255).
        # Only with the hi stream inside unnorm (the constants need the half thermometer).
        self.neg = neg
        # half: UNNORM on the half thermometer (network NET.BWD, one control per layer = CNOT sum of half flags
        # and parity, leading one moved by swaps) instead of the full thermometer and NET.FWD in reverse.
        self.half = half
        # zone: table format Z9c of the table-format builder (research/ws_invsqrt_tab.md). Zone flag
        # z = a6 AND (p OR a5 OR a4 OR (a3 AND a2)) on the address wires xq[16..20] and the parity. When z, norm
        # leaves the wires xq[0..14] rotated by one position (xq[i] holds xn bit i + 1, xq[14] holds 0) and
        # unnorm takes the rotation back. The exponent code then has another layout (own; the layout without
        # zones reads two code bits through gates that do not imply a zero wire in zone 1): F on its wire and,
        # on xq[0..4], G = [e < 18] (gate F), E14 and a AND G (gate G), E6 and b AND E14 (gate E14). The ring
        # cell on xq[7] moves to xq[3].
        self.zone = zone
        assert zb in ("z9c", "z9d") and zsb in ("pos", "gate", "chain"), (zb, zsb)
        self.zb, self.zsb = zb, zsb      # invsqrt_zpol options, see the module header
        self.sx_tries = sx_tries         # invsqrt_zpol: random multipliers per product count in the small-x synthesis
                                         # (60 = v4d); more tries find fewer products where they exist
        self._H = None                   # hi wires of the late hi stream (needed by zsb = "chain")
        self._gate8 = None
        self._own = None         # (wire, function of e) pairs offered to the small-x path (late hi stream: flags)
        if neg:
            lo_sign, lo_min_pos = -1, 11
            if angles_hi is None:
                angles_hi = {q: Fraction(-1, 1 << q) for q in range(0, 18)}
        self.c, self.xq, self.top_pos = c, xq, top_pos
        self.lo_min_pos, self.lo_sign = lo_min_pos, lo_sign
        self.angles_hi = angles_hi
        if angles_lo is None and lo_sign < 0:
            angles_lo = {q: Fraction(-1, 1 << q) for q in range(lo_min_pos, top_pos + 1)}
        self.angles_lo = angles_lo
        self.F = None
        self.ctl = None          # tier control wires while the code is off the t wires
        self._gpark = False      # control E6 parked on the wire of xn bit 6 (park_g)
        self.n_rot = 0
        self.info = {}

    # ------------------------------------------------------------ helpers
    def _or_new(self, x, y):
        c = self.c
        w = c.alloc()
        c.x(x)
        c.x(y)
        c.ccx(x, y, w)
        c.x(x)
        c.x(y)
        c.x(w)
        return w

    def _or_del(self, w, x, y):
        c = self.c
        c.x(w)
        m = c.hmr(w)
        c.x(x)
        c.x(y)
        c.cz(x, y, cond=m)
        c.x(x)
        c.x(y)
        c.bit_free(m)
        c.free(w)

    def _lin(self, f, g, w):
        """XOR g(e) onto w by CNOTs from the thermometer flags f[k] = [e < k], k = 1..21 (a function of e >> 1
        needs the even flags only)."""
        c = self.c
        if g(21):
            c.x(w)
        for k in range(1, 22):
            if g(k - 1) != g(k):
                c.cx(f[k], w)

    # ------------------------------------------------------------ NORM
    def norm(self, keep_tree=False):
        """keep_tree (invsqrt_int, zone format): the forward zone flag's AND tree t1..t4 and z stay alive and are
        handed to the caller in self.ztree (wires holding a3 a2, NOT a5 NOT a4, t2 NOT t1, t3 NOT p, a6 NOT t4)."""
        self.keep_tree = keep_tree
        c, xq = self.c, self.xq
        f = {}
        w = c.alloc()
        c.cx(xq[21], w)
        c.x(w)
        f[21] = w
        w = c.alloc()
        c.cx(xq[21], w)
        c.cx(xq[20], w)
        c.x(w)
        f[20] = w
        for k in range(19, 0, -1):
            w = c.alloc()
            c.x(xq[k])
            c.ccx(f[k + 1], xq[k], w)
            c.x(xq[k])
            f[k] = w
        for k in range(22):                                   # clear the leading one
            self._lin(f, lambda e, k=k: int(e == k), xq[k])
        cw = c.alloc()
        cur = frozenset()
        n_sw = 0
        for layer in NET.FWD["layers"]:
            for a, b, on, _off in layer["pairs"]:
                on = frozenset(on)
                self._lin(f, lambda e: int((e in on) != (e in cur)), cw)
                cur = on
                c.cswap(cw, xq[a], xq[b])
                n_sw += 1
        self._lin(f, lambda e: int(e in cur), cw)
        c.free(cw)
        self.info["norm_swaps"] = n_sw
        # Flags are released by measurement without a Toffoli: a flag that is a quadratic form (XOR of wires and
        # of products of two wires) of wires that stay is repaired by conditioned Z / CZ. Two rounds:
        # 1. parity onto wire 21; the odd flags are quadratic in parity and the even flags;
        # 2. clean code wires by CNOTs from the even flags; the even flags are quadratic in the code wires.
        self._lin(f, E_FUN["P"], xq[21])
        if self.zone:
            self._zone_rot(xq[21], forward=True)

        def release(flags, kept):
            names = list(kept)
            vec = {n: sum(kept[n][1](e) << e for e in range(22)) for n in names}
            cols, what = [(1 << 22) - 1], [()]
            for n in names:
                cols.append(vec[n])
                what.append((n,))
            for i, n in enumerate(names):
                for n2 in names[i + 1:]:
                    cols.append(vec[n] & vec[n2])
                    what.append((n, n2))
            for k in flags:
                sol = SX.solve(cols, sum(int(e < k) << e for e in range(22)))
                assert sol is not None, ("flag not in the span", k)
                m = c.hmr(f[k])
                for i in sol:
                    t = what[i]
                    if len(t) == 0:
                        c.neg(cond=m)
                    elif len(t) == 1:
                        c.z(kept[t[0]][0], cond=m)
                    else:
                        c.cz(kept[t[0]][0], kept[t[1]][0], cond=m)
                c.bit_free(m)
                c.free(f[k])
                del f[k]

        kept = {"P": (xq[21], E_FUN["P"])}
        for m in range(1, 11):
            kept["h%d" % m] = (f[2 * m], lambda e, m=m: int(e < 2 * m))
        release(range(1, 22, 2), kept)
        K = {"F": f[20]}
        for name in (("BE", "E14", "C2", "C3", "G") if self.zone else ("B", "E14", "C2", "C3", "G")):
            K[name] = c.alloc()
            self._lin(f, E_FUN[name], K[name])
        release(range(2, 20, 2), {n: (K[n], E_FUN[n]) for n in K})
        if self.zone:                                         # zone layout of the code, see __init__
            c.cx(K["G"], xq[0])
            c.cx(K["E14"], xq[1])
            c.cx(K["C2"], xq[2])
            c.cx(K["C3"], xq[3])
            c.cx(K["BE"], xq[4])
            c.and_del(K["C3"], xq[3], K["E14"])
            c.and_del(K["BE"], xq[4], K["E14"])
            c.and_del(K["C2"], xq[2], K["G"])
            c.and_del(K["E14"], xq[1], K["G"])
            c.and_del(K["G"], xq[0], K["F"])
            self.F = K["F"]
            return xq[21]
        # code onto the wires of xn bits 0..3, clean copies erased for free
        c.cx(K["B"], xq[0])
        c.cx(K["E14"], xq[1])
        c.cx(K["C2"], xq[2])
        c.cx(K["C3"], xq[3])
        c.and_del(K["C2"], xq[2], K["G"])
        c.and_del(K["C3"], xq[3], K["G"])
        self._or_del(K["G"], K["B"], K["E14"])
        c.and_del(K["B"], xq[0], K["F"])
        c.and_del(K["E14"], xq[1], K["F"])
        self.F = K["F"]
        return xq[21]

    # ------------------------------------------------------------ zone flag and rotation (format Z9c)
    def _zone_rot(self, pw, forward, zext=None, zt4=None):
        """z = a6 AND (p OR a5 OR a4 OR (a3 AND a2)) (5 ANDs), the wires xq[0..14] moved by one position under z
        (two reflections, sums 14 and 13: 14 swaps), z and its intermediates erased against the address wires.
        zext (invsqrt_int): z is on this wire already (the caller made it and erases it): swaps only.
        zt4 (invsqrt_int): the caller holds t4 = NOT a5 AND NOT a4 AND NOT (a3 AND a2) AND NOT p on this wire
        (e.g. as helpers of its unlookup chain): z = a6 AND NOT t4 is 1 AND, erased here; t4 stays the caller's."""
        c, xq = self.c, self.xq
        if zext is not None:
            if callable(zext):                                # the caller makes z in its hook, before this
                zext = zext()
            for sm in ((14, 13) if forward else (13, 14)):
                for q in range(7):
                    c.cswap(zext, xq[q], xq[sm - q])
            return
        if zt4 is not None:
            if callable(zt4):                                 # the caller makes t4 in its hook, before this
                zt4 = zt4()
            a6 = xq[20]
            c.x(zt4)
            z = c.and_new(a6, zt4)
            for sm in ((14, 13) if forward else (13, 14)):
                for q in range(7):
                    c.cswap(z, xq[q], xq[sm - q])
            c.and_del(z, a6, zt4)
            c.x(zt4)
            return
        a2, a3, a4, a5, a6 = (xq[16 + i] for i in range(5))
        if self.zb == "z9d":
            # invsqrt_zpol: z = a6 AND (p OR a5 OR a4) = a6 AND NOT (NOT a5 AND NOT a4 AND NOT p): 3 ANDs
            c.x(a5)
            c.x(a4)
            t2 = c.and_new(a5, a4)                            # NOT a5 AND NOT a4
            c.x(pw)
            t4 = c.and_new(t2, pw)                            # ... AND NOT p: the bracket is false
            c.x(t4)
            z = c.and_new(a6, t4)
            for sm in ((14, 13) if forward else (13, 14)):
                for q in range(7):
                    c.cswap(z, xq[q], xq[sm - q])
            c.and_del(z, a6, t4)
            c.x(t4)
            c.and_del(t4, t2, pw)
            c.x(pw)
            c.and_del(t2, a5, a4)
            c.x(a4)
            c.x(a5)
            return
        t1 = c.and_new(a3, a2)
        c.x(a5)
        c.x(a4)
        t2 = c.and_new(a5, a4)                                # NOT a5 AND NOT a4
        c.x(t1)
        t3 = c.and_new(t2, t1)                                # ... AND NOT (a3 a2)
        c.x(pw)
        t4 = c.and_new(t3, pw)                                # ... AND NOT p: the bracket is false
        c.x(t4)
        z = c.and_new(a6, t4)
        for sm in ((14, 13) if forward else (13, 14)):
            for q in range(7):
                c.cswap(z, xq[q], xq[sm - q])
        if forward and self.keep_tree:
            # invsqrt_int: the tree is handed to the caller as five functions of the address (pool of its first
            # lookup); the caller erases them (each is a quadratic form over its own two factors)
            c.x(t4)
            c.x(pw)
            c.x(t1)
            c.x(a4)
            c.x(a5)
            self.ztree = dict(t1=t1, t2=t2, t3=t3, t4=t4, z=z)
            return
        c.and_del(z, a6, t4)
        c.x(t4)
        c.and_del(t4, t3, pw)
        c.x(pw)
        c.and_del(t3, t2, t1)
        c.x(t1)
        c.and_del(t2, a5, a4)
        c.x(a4)
        c.x(a5)
        c.and_del(t1, a3, a2)

    # ------------------------------------------------------------ exponent code off / on the t wires
    def code_off(self, lean=True, park=False):
        """Tier controls onto 4 fresh wires (5 Toffolis); wires of xn bits 0..3 hold plain data again.
        lean: the flag wire F is released as well (4 wires in all); code_on then costs 2 Toffolis instead of 0.
        park: the control E6 is parked on the wire of xn bit 6 at once (3 wires in all), see park_g."""
        c, xq = self.c, self.xq
        assert self.ctl is None
        self._gpark = False
        if self.zone:
            G = c.and_new(xq[0], self.F)
            c.cx(G, xq[0])
            e14 = c.and_new(xq[1], G)
            c.cx(e14, xq[1])
            a = c.and_new(xq[2], G)
            c.cx(a, xq[2])
            c.cx(self.F, a)
            c.cx(G, a)                                        # a = (a AND G) xor F xor G
            self._lean = lean
            if lean:
                self._or_del(self.F, G, a)                    # at most 5 wires at a time
                self.F = None
            g = c.and_new(xq[3], e14)
            c.cx(g, xq[3])
            b = c.and_new(xq[4], e14)
            c.cx(b, xq[4])
            c.cx(G, b)
            c.cx(e14, b)                                      # b = G xor E14 xor (b AND E14)
            self.ctl = {"a": a, "b": b, "e": e14, "g": g}
            if lean:
                self._or_del(G, b, e14)
                self._G = None
            else:
                self._G = G
            if park:
                self.park_g()
            return self.ctl
        b = c.and_new(xq[0], self.F)
        c.cx(b, xq[0])
        e14 = c.and_new(xq[1], self.F)
        c.cx(e14, xq[1])
        G = self._or_new(b, e14)
        a = c.and_new(xq[2], G)
        c.cx(a, xq[2])
        c.cx(self.F, a)
        c.cx(G, a)
        self._lean = lean
        if lean:
            self._or_del(self.F, G, a)                        # F = G or a; at most 5 wires at a time
            self.F = None
        g = c.and_new(xq[3], G)
        c.cx(g, xq[3])
        self.ctl = {"a": a, "b": b, "e": e14, "g": g}
        if lean:
            self._or_del(G, b, e14)
            self._G = None
        else:
            self._G = G
        if park:
            self.park_g()
        return self.ctl

    def park_g(self):
        """Release the wire of the control E6 = [e < 6] while the code is off: E6 is XORed onto the wire of xn
        bit 6 (zero whenever E6 = 1) and its own wire is erased for free, because E6 = (that wire) AND E14.
        Allowed while nothing reads xn bit 6 as data and no stream runs: the rows of t4..t0 and the merges of B
        and B^-1 (bit 6 is no first digit of a group). 3 control wires stay. Taking it back: unpark_g, 1 Toffoli."""
        c, xq = self.c, self.xq
        assert self.ctl is not None and not self._gpark
        g = self.ctl["g"]
        c.cx(g, xq[6])
        c.and_del(g, xq[6], self.ctl["e"])
        self.ctl["g"] = None
        self._gpark = True

    def unpark_g(self):
        c, xq = self.c, self.xq
        assert self._gpark
        g = c.and_new(xq[6], self.ctl["e"])
        c.cx(g, xq[6])
        self.ctl["g"] = g
        self._gpark = False

    def code_on(self):
        c, xq = self.c, self.xq
        if self._gpark:
            self.unpark_g()
        a, b, e14, g = (self.ctl[k] for k in ("a", "b", "e", "g"))
        G = self._or_new(b, e14) if self._lean else self._G
        if self.zone:
            c.cx(G, b)
            c.cx(e14, b)                                      # wire holds b AND E14
            c.cx(g, xq[3])
            c.and_del(g, xq[3], e14)
            c.cx(b, xq[4])
            c.and_del(b, xq[4], e14)
            if self._lean:
                self.F = self._or_new(G, a)                   # at most 5 wires at a time
            c.cx(G, a)
            c.cx(self.F, a)                                   # wire holds a AND G
            c.cx(a, xq[2])
            c.and_del(a, xq[2], G)
            c.cx(e14, xq[1])
            c.and_del(e14, xq[1], G)
            c.cx(G, xq[0])
            c.and_del(G, xq[0], self.F)
            self.ctl = None
            return
        c.cx(g, xq[3])
        c.and_del(g, xq[3], G)
        if self._lean:
            self.F = self._or_new(G, a)
        c.cx(G, a)
        c.cx(self.F, a)
        c.cx(a, xq[2])
        c.and_del(a, xq[2], G)
        self._or_del(G, b, e14)
        c.cx(b, xq[0])
        c.and_del(b, xq[0], self.F)
        c.cx(e14, xq[1])
        c.and_del(e14, xq[1], self.F)
        self.ctl = None

    def _cells(self):
        cells = dict(CELLS, R3=3) if self.zone else CELLS     # zone 1: xq[7] holds a live bit at e = 13
        return {k: self.xq[v] for k, v in cells.items()}

    # ------------------------------------------------------------ PH1: hi stream and small-x path
    def ph1(self, H, small=None):
        """H: 8 wires, LSB first; H[7] may be None (bit 7 is the constant 1). Bit h goes to position j + 7 - h.
        small: None, or dict with H[x], lo[x] (integers for x = 1..63, lo as in ph2) and pool = list of
        (wire, truth table over x = 1..63 as an int with bit x): wires that stay as they are during PH1."""
        c = self.c
        assert not self.neg, "negated form: deliver the hi word through unnorm(pw, H, small, hook)"
        own = self.ctl is None
        if own:
            self.code_off()
        r = c.alloc()
        symbols = list(range(7)) + ([ST.ONE] if H[7] is None else [7])
        src = {k: H[k] for k in range(8) if H[k] is not None}
        st = ST.Stream(c, r, symbols, src, self.ctl, self._cells(), 17, 0, 17, self.angles_hi)
        inject = self._special(small) if small is not None else None
        st.run(inject)
        st.repair()
        c.free(r)
        self.n_rot += st.n_rot
        self.info["ph1_swaps"] = st.cost
        if own:
            self.code_on()

    def _special(self, small):
        c = self.c
        dom = sum(1 << x for x in range(1, 64))
        mp = self.lo_min_pos                                  # positions 0..mp-1 are served here for x < 64
        D = [0] * mp
        for x in range(1, 64):
            j = (x.bit_length() - 1) >> 1
            if self.neg:                                      # wire bits of the negated word: N, the constant
                have = (small["H"][x] << (mp - 8 - j)) + sum(1 << (mp - 1 - q) for q in range(j))    # ones above
                want = have + (small["lo"][x] >> (27 + j - mp))
            else:
                have = small["H"][x] << (mp - 8 - j)
                want = have + self.lo_sign * (small["lo"][x] >> (27 + j - mp))
            diff = (want ^ have) & ((1 << mp) - 1)
            for k in range(mp):
                if (diff >> k) & 1:
                    D[mp - 1 - k] |= 1 << x
        wires = [w for w, _ in small["pool"]]
        tts = [t & dom for _, t in small["pool"]]
        own = self._own                                       # own wires: functions of the exponent
        if own is None:
            own = [(self.ctl[n], lambda e, n=n: ST.ctl_val(n, d_of(e))) for n in ("a", "b")]
        for w_, fn in own:
            wires.append(w_)
            tts.append(sum(fn(x.bit_length() - 1) << x for x in range(1, 64)))
        inject = {}
        n_and = 0
        for q in range(mp):
            if D[q] == 0:
                continue
            L, prods = SX.synth(D[q], tts, dom, seed=1000 + q, tries=self.sx_tries)
            n_and += len(prods) + 1
            inject[q] = (lambda r, L=L, prods=prods: SX.emit_special(c, wires, L, prods, self.ctl["g"], r))
        self.info["special_positions"] = sorted(inject)
        if self.zone and self._gate8 is not None and any(x >= 64 for x in small["H"]):
            # Zone format: for x = 64..255 (j = 3, digits all zero in zone 0) the top lo wire can be set where
            # the digit word is not zero: the odd x of zone 1 at e = 7, whose bit a0 sits on xq[13]. The bits to
            # add at the positions above lo_min_pos are synthesised as  [e < 8] AND xq[13] AND B  with B an XOR
            # of pool wires (a linear solve on these inputs): 2 Toffolis per position.
            # (Only with the hi stream inside unnorm, where the flag [e < 8] is a wire.)
            xs2 = sorted(x for x in small["H"] if x >= 64)
            D2 = [0] * mp
            for x in xs2:
                j = (x.bit_length() - 1) >> 1
                if self.neg:
                    have = (small["H"][x] << (mp - 8 - j)) + sum(1 << (mp - 1 - q) for q in range(j))
                    want = have + (small["lo"][x] >> (27 + j - mp))
                else:
                    have = small["H"][x] << (mp - 8 - j)
                    want = have + self.lo_sign * (small["lo"][x] >> (27 + j - mp))
                diff = (want ^ have) & ((1 << mp) - 1)
                for k in range(mp):
                    if (diff >> k) & 1:
                        D2[mp - 1 - k] |= 1 << x
            full = [t for _, t in small["pool"]]
            i_a = [w for w, _ in small["pool"]].index(self.xq[13])
            t_a = full[i_a] & sum(1 << x for x in xs2)
            n2 = 0
            if self.zsb != "pos":
                n2 = self._zone_bits_shared(inject, wires, full, t_a, D2, mp, small)
                n_and += n2
                self.info["zone_special_positions"] = [q for q in range(mp) if D2[q]]
                self.info["zone_special_ands"] = n2
                self.info["special_ands"] = n_and
                return inject
            for q in range(mp):
                if D2[q] == 0:
                    continue
                assert D2[q] & ~t_a == 0, ("zone special bit outside the odd inputs of zone 1", q)
                sol = SX.solve([t_a] + [t & t_a for t in full], D2[q])
                assert sol is not None, ("zone special bit not in the span of the pool", q)
                form_b = SX.Form([i - 1 for i in sol if i > 0], const=int(0 in sol))
                prev = inject.get(q)

                def both(r, prev=prev, form_b=form_b):
                    f1 = prev(r) if prev is not None else None
                    if form_b.idx:
                        f2 = SX.emit_special(c, wires, SX.Form([]), [(SX.Form([i_a]), form_b)], self._gate8, r)
                    else:                                     # B is the constant 1: the bit is [e < 8] AND xq[13]
                        f2 = SX.emit_special(c, wires, SX.Form([i_a]), [], self._gate8, r)

                    def fix(m):
                        if f1 is not None:
                            f1(m)
                        f2(m)
                    return fix

                inject[q] = both
                n2 += 2 if form_b.idx else 1
            n_and += n2
            self.info["zone_special_positions"] = [q for q in range(mp) if D2[q]]
            self.info["zone_special_ands"] = n2
        self.info["special_ands"] = n_and
        return inject

    # ------------------------------------------------------------ zone special bits, shared gate (invsqrt_zpol)
    def _zone_bits_shared(self, inject, wires, full, t_a, D2, mp, small):
        """Zone special bits with ONE gate product t13 = [e < 8] AND xq13 for all positions (zsb = "gate" or
        "chain"); returns the Toffoli count. The stream visits positions from deep to shallow (step n = pos0 - n),
        so the first zone position is the largest q and the last the smallest.
        gate : bit_q = t13 AND B_q, B_q = an XOR of pool wires equal to D2[q] on the odd inputs of zone 1 at e = 7
               (the same targets as v4d); 1 Toffoli per position (0 if B_q = 1), erased at its step against t13.
        chain: bit at the top position = t13 AND B_top; the next one = that AND N_k (hi wire whose carry it is),
               and so on: needs D2[q] = D2[q + 1] AND (hi bit mp - 2 - q) on those inputs (asserted; else "gate").
        t13 (and in "chain" the whole chain) is erased after the fix of the last position, in reverse order."""
        c, xq, gate = self.c, self.xq, self._gate8
        zpos = sorted((q for q in range(mp) if D2[q]), reverse=True)
        if not zpos:
            return 0
        cols = [t_a] + [t & t_a for t in full]

        def form_of(target):
            sol = SX.solve(cols, target)
            assert sol is not None, ("zone special bit not in the span of the pool", target)
            return SX.Form([i - 1 for i in sol if i > 0], const=int(0 in sol))

        def collect(form):
            """the XOR form onto one of its pool wires by CNOTs (in place); returns (wire, undo)."""
            ops = []
            rb = min(form.idx)
            for i in sorted(form.idx - {rb}):
                c.cx(wires[i], wires[rb])
                ops.append(("cx", wires[i], wires[rb]))
            if form.const:
                c.x(wires[rb])
                ops.append(("x", wires[rb]))

            def undo():
                for op in reversed(ops):
                    if op[0] == "cx":
                        c.cx(op[1], op[2])
                    else:
                        c.x(op[1])
            return wires[rb], undo

        mode = self.zsb
        hi_tt = {}
        if mode == "chain":
            H = self._H or []
            pw_ = [w for w, _ in small["pool"]]
            for k, w in enumerate(H):
                if w is not None and w in pw_:
                    hi_tt[k] = full[pw_.index(w)] & t_a
            ok = zpos == list(range(zpos[0], zpos[-1] - 1, -1))
            for q in zpos[1:]:
                k = mp - 2 - q
                ok = ok and k in hi_tt and D2[q] == (D2[q + 1] & hi_tt[k])
            if not ok:
                mode = "gate"
        self.info["zone_special_mode"] = mode
        first, last = zpos[0], zpos[-1]
        st = {"t13": None, "stack": [], "prev": None}       # stack: (wire, a, b, form or None) to erase in reverse
        n_tof = [0]

        def make_bit(q, form_b, prev):
            def zbit(r):
                f1 = prev(r) if prev is not None else None
                if st["t13"] is None:
                    st["t13"] = c.and_new(gate, xq[13])
                    n_tof[0] += 1
                t13 = st["t13"]
                v, keep = None, None
                if mode == "gate" or q == first:
                    if form_b.idx:
                        bw, undo = collect(form_b)
                        v = c.and_new(t13, bw)
                        undo()
                        n_tof[0] += 1
                    else:                                    # B = 1 on these inputs: the bit is t13 itself
                        assert form_b.const
                    w_out = v if v is not None else t13
                else:                                        # chain: c_k = c_(k-1) AND N_(mp - 2 - q)
                    hk = self._H[mp - 2 - q]
                    v = c.and_new(st["prev"], hk)
                    n_tof[0] += 1
                    keep = (v, st["prev"], hk, None)
                    w_out = v
                c.cx(w_out, r)

                def fix(m):
                    if f1 is not None:
                        f1(m)
                    c.z(w_out, cond=m)
                    if mode == "gate":
                        if v is not None:
                            bw, undo = collect(form_b)
                            c.and_del(v, t13, bw)
                            undo()
                    else:
                        if q == first:
                            if v is not None:
                                st["stack"].append((v, t13, None, form_b))
                        else:
                            st["stack"].append(keep)
                        st["prev"] = w_out
                    if q == last:
                        for w_, a_, b_, fb in reversed(st["stack"]):
                            if fb is not None:
                                bw, undo = collect(fb)
                                c.and_del(w_, a_, bw)
                                undo()
                            else:
                                c.and_del(w_, a_, b_)
                        st["stack"] = []
                        c.and_del(t13, gate, xq[13])
                        st["t13"] = None
                        st["prev"] = None
                return fix
            return zbit

        for q in zpos:
            form_b = form_of(D2[q]) if (mode == "gate" or q == first) else None
            inject[q] = make_bit(q, form_b, inject.get(q))
        # Toffolis are counted at emission; the static count is known in advance:
        if mode == "gate":
            return 1 + sum(1 for q in zpos if form_of(D2[q]).idx)
        return 1 + int(bool(form_of(D2[first]).idx)) + (len(zpos) - 1)

    def _late_ph1(self, hf, H, small):
        c = self.c

        def lin_m(h, w):                                      # XOR h(e >> 1) onto w, hf[k] = [e >> 1 < k]
            if h(10):
                c.x(w)
            for k in range(1, 11):
                if h(k - 1) != h(k):
                    c.cx(hf[k], w)

        a2, b2 = c.alloc(), c.alloc()
        lin_m(lambda m: (10 - m) & 1, a2)
        lin_m(lambda m: ((10 - m) >> 1) & 1, b2)
        self.ctl = {"a": a2, "b": b2, "e": hf[7], "g": hf[3]}
        self._gate8 = hf[4]                                   # [e < 8], gate of the zone special bits
        self._H = list(H)                                     # invsqrt_zpol: hi wires, for zsb = "chain"
        self._own = [(hf[1], lambda e: int(e < 2)), (hf[2], lambda e: int(e < 4))]
        f_a, f_b = (lambda m: (10 - m) & 1), (lambda m: ((10 - m) >> 1) & 1)

        def lend(fn):
            """A special bit holds 2 to 4 temporaries from its step to the outcome of that step; no tier swap
            sits in between, so the two control copies are taken off for that time (CNOTs, free) and their
            wires are lent: the hi stream then needs no wire for its special bits beyond these two."""
            def inj(r):
                lin_m(f_a, self.ctl["a"])
                lin_m(f_b, self.ctl["b"])
                c.free(self.ctl["a"])
                c.free(self.ctl["b"])
                fix0 = fn(r)

                def fix(m):
                    fix0(m)
                    self.ctl["a"], self.ctl["b"] = c.alloc(), c.alloc()
                    lin_m(f_a, self.ctl["a"])
                    lin_m(f_b, self.ctl["b"])
                return fix
            return inj

        def T(k):                                             # [d >= k] as (wires, constant)
            return ([], 1) if k <= 0 else (([], 0) if k >= 11 else ([hf[11 - k]], 0))

        ind = {d: (T(d)[0] + T(d + 1)[0], T(d)[1] ^ T(d + 1)[1]) for d in ST.DELAYS}
        r = c.alloc()
        const7 = H[7] is None
        symbols = list(range(7)) + ([] if const7 else [7])
        src = {k: H[k] for k in range(8) if H[k] is not None}
        ext = None
        if self.neg:                                          # no bit at index 7; ones at every index >= 8
            assert const7, "negated form: the hi word has 7 wires"
            ext = (8, lambda q: ([hf[10 - q]], 1) if q <= 9 else ([], 1))      # [d <= q] = 1 xor [d >= q + 1]
        st = ST.Stream(c, r, symbols, src, self.ctl, self._cells(), 17, 0, 17, self.angles_hi,
                       ind=ind, one_idx=7 if (const7 and not self.neg) else None, ext=ext)
        inject = self._special(small) if small is not None else None
        if inject:
            inject = {q: lend(fn) for q, fn in inject.items()}
        st.run(inject)
        st.repair()
        c.free(r)
        self.n_rot += st.n_rot
        self.info["ph1_swaps"] = st.cost
        lin_m(f_a, self.ctl["a"])
        lin_m(f_b, self.ctl["b"])
        c.free(self.ctl["a"])
        c.free(self.ctl["b"])
        self.ctl = None
        self._own = None

    # ------------------------------------------------------------ PH2: lo stream
    def ph2(self, lo, top_one=False):
        """lo: wires LSB first; wire i goes to position 26 + j - i when that position is in lo_min_pos..top_pos.
        top_one (test of the one-bank design only): a constant 1 above the last wire; with 26 wires, top_one and
        lo_min_pos = 0 the whole phase word goes through this one stream (33 rotations, T, S, Z)."""
        c = self.c
        own = self.ctl is None
        if own:
            self.code_off()
        parked = self._gpark                                  # E6 parked by the caller: taken back for the stream
        if parked:                                            # and parked again as soon as the repair allows
            self.unpark_g()
        r = c.alloc()
        src = {k: w for k, w in enumerate(lo)}
        symbols = list(range(len(lo))) + ([ST.ONE] if top_one else [])
        st = ST.Stream(c, r, symbols, src, self.ctl, self._cells(), 36, self.lo_min_pos, self.top_pos,
                       self.angles_lo)
        st.run()
        st.repair(after_g=self.park_g if parked else None)
        c.free(r)
        self.n_rot += st.n_rot
        self.info["ph2_swaps"] = st.cost
        if own:
            self.code_on()

    # ------------------------------------------------------------ UNNORM
    def unnorm(self, pw=None, H=None, small=None, hook=None, zext=None, hook2=None, zt4=None):
        """pw: wire holding the parity (default: wire 21). It is erased here; wire 21 must otherwise be |0>.
        H (optional): the hi word (as in ph1) is delivered here instead of by ph1, after the half thermometer is
        back: tier controls and delay indicators are then CNOT sums of flags, so this "late PH1" needs neither
        the code moves nor Toffolis for the repair, and the constant bit H[7] = None costs nothing.
        The full thermometer is rebuilt from the exponent code and the parity (5 + 4 + 10 Toffolis), so the swaps
        are those of NORM in reverse order with the same free controls, the leading one is set by CNOTs and the
        flags are erased for free once x is back."""
        c, xq = self.c, self.xq
        if pw is None or pw == xq[21]:
            pw = c.alloc()
            c.cx(xq[21], pw)
            c.cx(pw, xq[21])
        self.code_off(lean=False)
        a, b, e14, g = (self.ctl[k] for k in ("a", "b", "e", "g"))
        G, F = self._G, self.F
        self.ctl = None
        # half thermometer hf[m] = [e < 2m] = [d >= 11 - m], 4 Toffolis
        c.cx(G, b)
        c.cx(e14, b)                                          # b AND E14 = b_mid xor b_g
        bg = c.and_new(g, b)                                  # [d = 10]
        c.cx(bg, b)                                           # b_mid = [d in 6, 7]
        c.cx(e14, G)
        t3 = c.and_new(a, G)                                  # [d = 3]
        c.cx(e14, G)
        c.cx(F, a)
        c.cx(G, a)
        c.cx(t3, a)                                           # a AND E14 = a_mid xor a_g
        ag = c.and_new(g, a)                                  # [d = 9]
        c.cx(ag, a)                                           # a_mid = [d in 5, 7]
        t7 = c.and_new(a, b)                                  # [d = 7]
        c.cx(b, a)
        c.cx(t7, a)
        c.cx(g, a)                                            # [d >= 5]
        c.cx(g, b)                                            # [d >= 6]
        c.cx(g, t7)                                           # [d >= 7]
        c.cx(e14, t3)                                         # [d >= 3]
        c.cx(bg, ag)                                          # [d >= 9]
        hf = {10: F, 9: G, 8: t3, 7: e14, 6: a, 5: b, 4: t7, 3: g, 2: ag, 1: bg}
        if H is not None:
            self._late_ph1(hf, H, small)
            if hook is not None:                              # caller may release its hi wires here
                hook()
        if self.zone:
            self._zone_rot(pw, forward=False, zext=zext, zt4=zt4)
        if hook2 is not None:                                 # invsqrt_int: caller erases its zone flag here
            hook2()
        if self.half:
            self._unnorm_half(hf, pw)
            return
        # odd flags: [e < 2m + 1] = hf[m + 1] xor p AND [e >> 1 = m], 10 Toffolis
        f = {2 * m: hf[m] for m in range(1, 11)}
        for m in range(10):
            if m >= 1:
                c.cx(hf[m], hf[m + 1])                        # wire holds [e >> 1 = m]
            w = c.and_new(pw, hf[m + 1])
            if m >= 1:
                c.cx(hf[m], hf[m + 1])
            f[2 * m + 1] = w
        for m in range(10):
            c.cx(f[2 * m + 1], pw)                            # pw = p AND [e >= 20] = [e = 21]
        c.x(pw)
        f[21] = pw
        for m in range(10):
            c.cx(hf[m + 1], f[2 * m + 1])
        cw = c.alloc()
        cur = frozenset()
        n_sw = 0
        for layer in reversed(NET.FWD["layers"]):
            for qa, qb, on, _off in reversed(layer["pairs"]):
                on = frozenset(on)
                self._lin(f, lambda e: int((e in on) != (e in cur)), cw)
                cur = on
                c.cswap(cw, xq[qa], xq[qb])
                n_sw += 1
        self._lin(f, lambda e: int(e in cur), cw)
        c.free(cw)
        self.info["unnorm_swaps"] = n_sw
        for k in range(22):                                   # set the leading one
            self._lin(f, lambda e, k=k: int(e == k), xq[k])
        for k in range(1, 20):                                # f[k] = f[k+1] AND NOT x_k
            m = c.hmr(f[k])
            c.cz(f[k + 1], xq[k], cond=m)
            c.z(f[k + 1], cond=m)
            c.bit_free(m)
            c.free(f[k])
        c.x(f[20])
        c.cx(xq[20], f[20])
        c.cx(xq[21], f[20])
        c.free(f[20])
        c.x(f[21])
        c.cx(xq[21], f[21])
        c.free(f[21])
        self.F = None

    def _unnorm_half(self, hf, pw):
        """Swaps of NET.BWD under layer controls alpha(e >> 1) xor beta p, then parity and half flags erased
        against x: the parity for free, the flags [e < 2m], m = 1..9, by one CCZ each under the outcome."""
        c, xq = self.c, self.xq

        def lin_m(h, w):                                      # XOR h(e >> 1) onto w, hf[k] = [e >> 1 < k]
            if h[10]:
                c.x(w)
            for k in range(1, 11):
                if h[k - 1] != h[k]:
                    c.cx(hf[k], w)

        # The leading one. Without further data in NET.BWD it starts on wire 21 and is moved by the swaps.
        # "lead_in": per e >> 1 a [layer index, wire] where it enters instead, by CNOTs from the flags (free).
        # "lead_and": the e >> 1 for which no leading one enters at all; it is set by one AND after the swaps.
        lin = NET.BWD.get("lead_in")
        la = NET.BWD.get("lead_and", [])

        def flag_k(k, w):                                     # XOR [e >> 1 = k] onto w
            if k == 10:
                c.x(w)
                c.cx(hf[10], w)
            else:
                c.cx(hf[k + 1], w)
                if k >= 1:
                    c.cx(hf[k], w)

        def lead_in(idx):
            for k in range(11):
                if k not in la and lin[k][0] == idx:
                    flag_k(k, xq[lin[k][1]])

        if lin is None:
            c.x(xq[21])
            for k in la:
                flag_k(k, xq[21])
        cw = c.alloc()
        cur_a, cur_b, n_sw = [0] * 11, 0, 0
        for li, layer in enumerate(NET.BWD["layers"]):
            if lin is not None:
                lead_in(li)
            al, be = layer["alpha"], layer["beta"]
            lin_m([x ^ y for x, y in zip(al, cur_a)], cw)
            if be != cur_b:
                c.cx(pw, cw)
            cur_a, cur_b = al, be
            for qa, qb in layer["pairs"]:
                c.cswap(cw, xq[qa], xq[qb])
                n_sw += 1
        lin_m(cur_a, cw)
        if cur_b:
            c.cx(pw, cw)
        c.free(cw)
        if lin is not None:
            lead_in(len(NET.BWD["layers"]))
        self.info["unnorm_swaps"] = n_sw
        for k in la:                                          # x_(2k+1) ^= [e = 2k + 1], x_(2k) ^= [e = 2k]
            if k == 10:
                c.x(hf[10])                                   # wire holds [e >> 1 = 10]
                iw = hf[10]
            else:
                if k >= 1:
                    c.cx(hf[k], hf[k + 1])                    # wire holds [e >> 1 = k]
                iw = hf[k + 1]
            t = c.and_new(pw, iw)
            c.cx(t, xq[2 * k + 1])
            c.cx(iw, xq[2 * k])
            c.cx(t, xq[2 * k])
            c.and_del(t, pw, iw)
            if k == 10:
                c.x(hf[10])
            elif k >= 1:
                c.cx(hf[k], hf[k + 1])
        self.info["unnorm_lead_ands"] = len(la)
        mb_ = c.hmr(pw)                                       # p = sum over m of [e >> 1 = m] x_(2m+1)
        for m in range(11):
            if m == 10:
                c.z(xq[21], cond=mb_)
            else:
                c.cz(hf[m + 1], xq[2 * m + 1], cond=mb_)
            if m >= 1:
                c.cz(hf[m], xq[2 * m + 1], cond=mb_)
        c.bit_free(mb_)
        c.free(pw)
        for m in range(1, 10):                                # hf[m] = hf[m+1] AND NOT x_(2m+1) AND NOT x_(2m)
            mm = c.hmr(hf[m])
            c.x(xq[2 * m + 1])
            c.x(xq[2 * m])
            c.ccz(hf[m + 1], xq[2 * m + 1], xq[2 * m], cond=mm)
            c.x(xq[2 * m + 1])
            c.x(xq[2 * m])
            c.bit_free(mm)
            c.free(hf[m])
        mm = c.hmr(hf[10])                                    # [e < 20] = NOT x21 AND NOT x20
        c.neg(cond=mm)
        c.z(xq[21], cond=mm)
        c.z(xq[20], cond=mm)
        c.cz(xq[21], xq[20], cond=mm)
        c.bit_free(mm)
        c.free(hf[10])
        self.F = None
