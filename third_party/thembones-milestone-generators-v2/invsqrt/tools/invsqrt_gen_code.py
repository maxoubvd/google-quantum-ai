"""invsqrt_gen: the exponent between NORM and UNNORM on ONE flag wire plus five t wires (own code).

Idea credited to the best public invsqrt circuit (read in research/ws_invsqrt.md and ws_invsqrt_decode.md 3.3: it
keeps its exponent code XORed onto t wires that are zero for the exponents where the code bit is set, plus one
flag). The code below, its read-out and the free erasures are our own construction for our own tier controls
(binary delay d = 10 - (e >> 1), bits b0..b3).

Facts used [P]: the normalised input xn = x << s has s = 2 d + 1 - p >= 2 d zero low bits, so
    d >= 1  =>  t0, t1 are zero        d >= 2  =>  t0..t3 are zero        d >= 4  =>  t0..t7 are zero.
Parked form (state "parked", one wire F = [d >= 1] besides the t wires):
    t0 ^= G   (G = [d >= 2])      t1 ^= b0      t2 ^= b1      t3 ^= B2   (B2 = [d >= 4] = b2 xor b3)      t4 ^= b3
Every parked bit is non-zero only where its t wire is zero, and it can be read back with one AND each:
    G = F AND t0     b0 = F AND t1     b1 = G AND t2     B2 = G AND t3     b3 = B2 AND t4          (5 Toffolis)
Putting the code back costs nothing: after the CNOTs each of the five wires is the AND of two live wires again
and is erased by HMR + conditioned CZ.
While the code is off the t wires, F and G are quadratic forms of the control wires (b2 b3 = 0 for d <= 10):
    G = b1 xor b2 xor b3 xor b1 b2 xor b1 b3            F = b0 xor G xor b0 G
so they can be measured away for free (state "ctl": only the four control wires) and cost one AND each to
raise again (G = b1 OR (b2 xor b3), F = b0 OR G).

States:  "off"    F, G, b0..b3 on six wires, t0..t4 clean (as after NORM)
         "ctl"    b0..b3 on four wires, t0..t4 clean
         "parked" F on one wire, code on t0..t4 (rows of t13..t5 and every stage that does not read t0..t4 or the
                  tier controls can run in this state)
Hooks for a stage emitter: code.unpark(flags=True|False), code.drop_flags(), code.park(), code.b (the four
control wires, valid in "off" and "ctl"), code.F.
"""
from __future__ import annotations


class ExpCode:
    def __init__(self, c, xq, F, G, b):
        self.c = c
        self.t = [xq[k] for k in range(5)]
        self.F, self.G, self.b = F, G, list(b)
        self.state = "off"
        self.moves = dict(unpark=0, park=0, drop=0, raise_=0)

    def drop_flags(self):
        assert self.state == "off"
        c, b = self.c, self.b
        m = c.hmr(self.F)                      # F = b0 xor G xor b0 G
        c.z(b[0], cond=m)
        c.z(self.G, cond=m)
        c.cz(b[0], self.G, cond=m)
        c.bit_free(m)
        c.free(self.F)
        m = c.hmr(self.G)                      # G = b1 xor b2 xor b3 xor b1 b2 xor b1 b3
        for w in b[1:]:
            c.z(w, cond=m)
        c.cz(b[1], b[2], cond=m)
        c.cz(b[1], b[3], cond=m)
        c.bit_free(m)
        c.free(self.G)
        self.F = self.G = None
        self.state = "ctl"
        self.moves["drop"] += 1

    def raise_flags(self):
        assert self.state == "ctl"
        c, b = self.c, self.b
        c.cx(b[3], b[2])                       # b[2] holds B2 = b2 xor b3 for a moment
        G = c.and_new(b[1], b[2])
        c.cx(b[1], G)
        c.cx(b[2], G)                          # G = b1 OR B2
        c.cx(b[3], b[2])
        F = c.and_new(b[0], G)
        c.cx(b[0], F)
        c.cx(G, F)                             # F = b0 OR G
        self.F, self.G = F, G
        self.state = "off"
        self.moves["raise_"] += 1

    def park(self):
        if self.state == "ctl":
            self.raise_flags()
        assert self.state == "off"
        c, t, b, F, G = self.c, self.t, self.b, self.F, self.G
        c.cx(b[3], b[2])                       # B2
        for src, dst in ((G, t[0]), (b[0], t[1]), (b[1], t[2]), (b[2], t[3]), (b[3], t[4])):
            c.cx(src, dst)
        c.and_del(b[3], b[2], t[4])
        c.and_del(b[2], G, t[3])
        c.and_del(b[1], G, t[2])
        c.and_del(b[0], F, t[1])
        c.and_del(G, F, t[0])
        self.G, self.b = None, [None] * 4
        self.state = "parked"
        self.moves["park"] += 1

    def unpark(self, flags=True):
        assert self.state == "parked"
        c, t, F = self.c, self.t, self.F
        G = c.and_new(F, t[0])
        b0 = c.and_new(F, t[1])
        b1 = c.and_new(G, t[2])
        B2 = c.and_new(G, t[3])
        b3 = c.and_new(B2, t[4])
        for src, dst in ((b3, t[4]), (B2, t[3]), (b1, t[2]), (b0, t[1]), (G, t[0])):
            c.cx(src, dst)
        c.cx(b3, B2)                           # b2 = B2 xor b3
        self.G, self.b = G, [b0, b1, B2, b3]
        self.state = "off"
        self.moves["unpark"] += 1
        if not flags:
            self.drop_flags()
