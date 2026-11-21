"""dqi_int_j6blk.py -- workstream dqi_int (B5): the prefix block at J = 6 on the kit's primitives.

Runs schedule B5-1 of tools/dqi_int_j6sched.py (the same code that the phase model checks: `sched` and the negative
control) through a machine adapter whose methods are kit calls:
  alloc / free -> b.alloc_qubits(12) / b.free;  cx -> b.cx or the namespace's cx (constant multiply-add);
  scale; sq -> square12; mulx -> ar.mul_xor0 (42 ANDs); hmr -> b.alloc_bits(12) + b.hmr + b.free;
  cz -> b.cz; pp -> ar.phase_product; czc / ppc -> scale by the constant, cz / phase_product, scale back;
  loop -> the polynomial loop (patch P) for U' of degree 5: 33 products, plus the zero node's product.
Registers: heads a1..b5 = words[2..11], N_g / L_g = words[g.word + 0, 1] / [+2, +3], the zero node's X / L = words[68] /
words[69]. The six U' words are the clean points PTS[0..4] and U'(0) in place in b1.
Installed by tools/dqi_int_pfdp.py under DQI_INT_J6=1 (it also makes head 5 N-only in the chain); the zero node's
chain skips head 5 with the overlay copy of dqi_pub_zt.py (tools/dqi_int_zt6.py) under the same switch.
"""
import sys

import dqi_int_j6sched as S

STATS = dict(forward=0, reverse=0, products=[], marks=0)


class Kit:
    def __init__(self, b, words, G):
        self.b = b
        self.ar, self.cxk, self.scalek, self.square12 = G['ar'], G['cx'], G['scale'], G['square12']
        self.mul, self.inv = G['mul'], G['inv']
        self.G = G
        self.words = words
        groups = G['groups']
        self.regs = {}
        for i in range(1, 6):
            self.regs['a%d' % i] = words[2 * i]
            self.regs['b%d' % i] = words[2 * i + 1]
        for g, _ in S.NODES:
            base = groups[g]['word']
            for h in (0, 1):
                self.regs[('N', g, h)] = words[base + h]
                self.regs[('L', g, h)] = words[base + 2 + h]
        self.regs[('N', 'z', 0)] = words[68]
        self.regs[('L', 'z', 0)] = words[69]
        self.prods = 0

    def w(self, n):
        return self.regs[n]

    def alloc(self, n):
        assert n not in self.regs, n
        self.regs[n] = self.b.alloc_qubits(12)

    def free(self, n):
        self.b.free(self.regs.pop(n))

    def cx(self, s, d, c=1):
        if c == 1:
            self.b.cx(self.w(s), self.w(d))
        elif c:
            self.cxk(self.b, self.w(s), self.w(d), c)

    def scale(self, n, c):
        assert c
        if c != 1:
            self.scalek(self.b, self.w(n), c)

    def sq(self, n, undo=False):
        self.square12(self.b, self.w(n), undo)

    def mulx(self, x, y, z):
        self.ar.mul_xor0(self.b, self.w(x), self.w(y), self.w(z))
        self.prods += 1

    def hmr(self, n):
        m = self.b.alloc_bits(12)
        self.b.hmr(self.w(n), m)
        self.free(n)
        return m

    def cz(self, n, m):
        self.b.cz(self.w(n), m)

    def pp(self, x, y, m):
        self.ar.phase_product(self.b, self.w(x), self.w(y), m)

    def czc(self, n, m, c):
        if c:
            self.scale(n, c)
            self.cz(n, m)
            self.scale(n, self.inv(c))

    def ppc(self, x, y, m, c):
        if c:
            self.scale(y, c)
            self.pp(x, y, m)
            self.scale(y, self.inv(c))

    def mfree(self, m):
        self.b.free(m)

    def form(self, a, bb, beta, body):
        self.cx(bb, a, beta)
        self.sq(a)
        body(a)
        self.sq(a, True)
        self.cx(bb, a, beta)

    def cz_lin(self, terms, m):
        terms = [(r, c) for r, c in terms if c]
        if not terms:
            return
        (h, c0), rest = terms[0], terms[1:]
        self.scale(h, c0)
        for r, c in rest:
            self.cx(r, h, c)
        self.cz(h, m)
        for r, c in rest:
            self.cx(r, h, c)
        self.scale(h, self.inv(c0))

    def loop(self):
        kit_loop(self)


def _econst(c):
    """The model's E constant (int lo | hi << 12) as the kit's (lo, hi)."""
    return c & S.MASK, c >> S.MB


def kit_loop(kit):
    """Patch P for U' of degree 5: all 14 N_g hold X_g^2; for each of 33 points t (F elements 0..31 and infinity) U'(t)
    is exposed on one clean U' word, XX(t) = sum_g Tr(X_g^2 lam_g(t)) on one N word, one 42-AND product into one L word
    inside a free spread that adds (1 / sc_g) W_g,t p_t to every L_g (sc_g = Q0^2 Q1^2 Q2^2). Then the zero node:
    zL ^= U'(0) (x_0 + a_0)^2 / (k0^2 k1^2 k2^2)."""
    b, G = kit.b, kit.G
    mul, inv, K, combine, square24, groups = G['mul'], G['inv'], G['K'], G['combine'], G['external_square24'], G['groups']
    cx, scale, ar, square12 = kit.cxk, kit.scalek, kit.ar, kit.square12
    words = kit.words
    a0, h0 = words[0], words[1]

    def emul(u, v):
        t = mul(u[1], v[1])
        return (mul(u[0], v[0]) ^ mul(K, t), mul(u[0], v[1]) ^ mul(u[1], v[0]) ^ t)

    def einv(u):
        a, c = u
        n = mul(a, a) ^ mul(a, c) ^ mul(K, mul(c, c))
        ni = inv(n)
        return (mul(a ^ c, ni), mul(c, ni))

    def eadd(u, v):
        return (u[0] ^ v[0], u[1] ^ v[1])

    def escale(u, s):
        return (mul(u[0], s), mul(u[1], s))

    def cxc(src, dst, c):
        if c == 1:
            b.cx(src, dst)
        elif c:
            cx(b, src, dst, c)

    hom = [g for g, _ in S.NODES]
    W = {gi: (groups[gi]['c'] ^ groups[gi]['kap'], 1) for gi in hom}
    allnodes = [W[gi] for gi in hom] + [(W[gi][0] ^ W[gi][1], W[gi][1]) for gi in hom]
    assert len(set(allnodes)) == len(allnodes)
    pts = list(range(len(allnodes) + 4))           # deg U' XX = 5 + 27 = 32: 32 finite points + infinity
    upts = list(S.UPTS)
    Uw = {p: kit.w(('U', p)) for p in upts[:5]}
    Uw[0] = kit.w('b1')                             # U'(0) in place
    host = upts[0]

    def lam_u(t):                                   # weights of the six words in U'(t); t None: leading coefficient
        out = {}
        for p in upts:
            num, den = 1, 1
            for r in upts:
                if r != p:
                    if t is not None:
                        num = mul(num, t ^ r)
                    den = mul(den, p ^ r)
            out[p] = mul(num, inv(den))
        return out

    def lam_x(node, t):
        num, den = (1, 0), (1, 0)
        for nu in allnodes:
            if nu != node:
                if t is not None:
                    num = emul(num, eadd((t, 0), nu))
                den = emul(den, eadd(node, nu))
        return emul(num, einv(den))

    def spread_c(t, w):
        if t is None:
            out = (1, 0)
            for s in pts:
                out = emul(out, eadd(w, (s, 0)))
            return out
        num, den = (1, 0), 1
        for s in pts:
            if s != t:
                num = emul(num, eadd(w, (s, 0)))
                den = mul(den, t ^ s)
        return escale(num, inv(den))

    def sq(v):
        return mul(v, v)

    k0, k1, k2 = S.K0, S.K1, S.K2
    sc = {gi: mul(mul(sq(k0 ^ groups[gi]['kap']), sq(k1 ^ groups[gi]['kap'])), sq(k2 ^ groups[gi]['kap'])) for gi in hom}
    base = {gi: groups[gi]['word'] for gi in hom}
    Nw = {gi: (words[base[gi]], words[base[gi] + 1]) for gi in hom}
    Lw = {gi: (words[base[gi] + 2], words[base[gi] + 3]) for gi in hom}

    def n_to_x2(gi, undo):
        N = combine(words[base[gi]:base[gi] + 2])
        ct = groups[gi]['c']
        if undo:
            square24(b, N, True)
        cx(b, h0, a0, ct)
        b.cx(combine([a0, h0]), N)
        cx(b, h0, a0, ct)
        if not undo:
            square24(b, N)

    for gi in hom:
        n_to_x2(gi, False)
    for t in pts + [None]:
        lu = None
        if t is not None and t in Uw:
            ux = Uw[t]
        else:
            lu = lam_u(t)
            ux = Uw[host]
            scale(b, ux, lu[host])
            for p in upts:
                if p != host:
                    cxc(Uw[p], ux, lu[p])
        cf = []
        for gi in hom:
            lam = lam_x(W[gi], t)
            cf.append((Nw[gi][0], lam[1]))
            cf.append((Nw[gi][1], lam[0] ^ lam[1]))
        yi = next(i for i, (_, c) in enumerate(cf) if c)
        yw, yc = cf[yi]
        scale(b, yw, yc)
        for i, (w_, c) in enumerate(cf):
            if i != yi and c:
                cx(b, w_, yw, c)
        cc = []
        for gi in hom:
            c = escale(spread_c(t, W[gi]), inv(sc[gi]))
            cc.append((Lw[gi][0], c[0]))
            cc.append((Lw[gi][1], c[1]))
        ti = next(i for i, (_, c) in enumerate(cc) if c)
        tw, tc = cc[ti]
        scale(b, tw, inv(tc))
        for i, (w_, c) in enumerate(cc):
            if i != ti and c:
                cx(b, tw, w_, c)
        ar.mul_xor0(b, ux, yw, tw)
        kit.prods += 1
        for i, (w_, c) in enumerate(cc):
            if i != ti and c:
                cx(b, tw, w_, c)
        scale(b, tw, tc)
        for i, (w_, c) in enumerate(cf):
            if i != yi and c:
                cx(b, w_, yw, c)
        scale(b, yw, inv(yc))
        if lu is not None:
            for p in upts:
                if p != host:
                    cxc(Uw[p], ux, lu[p])
            scale(b, ux, inv(lu[host]))
    for gi in hom:
        n_to_x2(gi, True)
    zX, zL = words[68], words[69]                   # the zero node (w = 0, kappa 0)
    zc = mul(mul(sq(k0), sq(k1)), sq(k2))
    b.cx(a0, zX)
    square12(b, zX)
    scale(b, zL, zc)
    ar.mul_xor0(b, Uw[0], zX, zL)
    kit.prods += 1
    scale(b, zL, inv(zc))
    square12(b, zX, True)
    b.cx(a0, zX)


def check_field(G):
    """The model's constants are valid for this kit namespace only if the field, the nodes and K agree."""
    mul, groups = G['mul'], G['groups']
    for x, y in ((3, 5), (4095, 2047), (1234, 3210), (17, 4000), (2048, 2048)):
        assert mul(x, y) == S.cmul(x, y), ('field differs', x, y)
    assert G['K'] == S.ch.K, 'K differs'
    for g in range(3):
        assert groups[g]['kap'] == S.CKAP[g], ('kappa differs', g)
    for g, _ in S.NODES:
        assert groups[g]['kap'] == S.CKAP[g] and groups[g]['c'] == S.CCS[g], ('node differs', g)
    assert G['ks'][0] == S.K0 and G['ks'][2] == S.K1 and G['ks'][4] == S.K2, 'ks differs'


def make_j6_block(G, direction, mask_stack):
    """block(b, words, homogeneous) for pfd's install wrapper: direction(b) is True on the forward call; the forward's
    masks go on mask_stack[id(b)] and the return pops them."""
    check_field(G)

    def block(b, words, homogeneous):
        assert sorted(homogeneous) == sorted(g for g, _ in S.NODES), 'homogeneous groups differ from the model'
        kit = Kit(b, words, G)
        fwd = direction(b)
        print(f'dqb_prefix6 (J = 6, B5-1) block entry ({"forward" if fwd else "return"}): allocated '
              f'{b.num_allocated_qubits}', file=sys.stderr, flush=True)
        with b.mark('dqb_prefix5'):
            if fwd:
                mask_stack.setdefault(id(b), []).append(S.block_forward(kit))
                STATS['forward'] += 1
            else:
                S.block_return(kit, mask_stack[id(b)].pop())
                STATS['reverse'] += 1
        STATS['products'].append(kit.prods)
    return block
