"""r46_dqi_syn_t4.py -- seat R46-DQI-SYN: the synthesis tail's correction with its masked operand in an 11-WIRE
register that walks three hyperplanes of forms (mode t4 of tools/r46_dqi_syn_tail.py imports this file).

The correction (kit `a07/decode_correction_and.term_correction_with_inverse`) adds v ^= F p / q and u ^= F v / q,
F = [d = 0] AND NOT pole, qi = 1 / q live.  The kit's product is 42 rank-one terms, each (form m of x) AND (the same
form m of y) XORed into the word val_m.  A register of 11 wires holding F * lambda_k(qi) for a basis lambda of a
hyperplane H = {m : par(m & h) = 0} serves every term with m in H; ONE Toffoli CCX(F, g(qi), wire) moves it to another
hyperplane.  No two hyperplanes hold all 42 forms [M tools/r46_dqi_syn_pair.py: 0 covering pairs]; three do (plan()).
Per call: two products, the first walks H1 > H2 > H3, the second H3 > H2 > H1: 4 toggles.
  forward  F (1); Z X-measured (mask deferred, as today); register (11); 84 + 4; register and F X-measured with
           DEFERRED masks (the return pays them by classically conditioned CZ(F, form(qi)) and CZ(F): free)     100
  reverse  Z = [d = 0] (11); F (1); the forward's masks (free); Z X-measured; register (11); 84 + 4; register
           erased (free: F and qi live); F X-measured; Z tested again (11; d is untouched); F again (1) to pay
           F's mask, erased free                                                                              123
Live beside qi 12, Zq, pole during the products: F + 11 = base + 26 = 902 in both directions [P].
Line (J's tailctl): 169 + 180 = 349.  Expected: -69 forward, -57 reverse = -126 per decode at Q 902.

Self-test (classical bits; phases are not modelled):  python r46_dqi_syn_t4.py test [kit decoder dir]

RESULTS [M box 6, Nov 15; research/panel/R46-DQI-SYN.md]: best_q23j5's recipe with tailctl out and this module in
(R46_TAIL_MODE=t4): raw 688684a4825b8169, Q 902 (every stage at or under 902), 350,102.3 against 350,228.2 (-125.9);
quick 0 / 0; official verifier 0 / 0 on both salts.  Negative control (mode t4m: the return does not pay the 11
deferred register masks): 2,580 phase errors in 10,240 shots, so F = 1 on about every second shot.  After P1's dwx
pass: 53f4a5c8cbc1a960, 348,670; paired against the line's file on 6 salts: child-only failures 0 of 614,400 shots;
re-roll r46t4a (986599ed01bc1170): 0 / 0 on both official salts, 3.37936e11.
"""
import os
import sys

LOG = {}


def _bits(x):
    return [i for i in range(x.bit_length()) if x >> i & 1]


def _par(x):
    return bin(x).count('1') & 1


_PLAN = {}


def plan(terms):
    """Three hyperplanes (h1, h2, h3) of GF(2)^12 such that every operand form lies in at least one of them."""
    masks = tuple(m for m, _ in terms)
    if masks in _PLAN:
        return _PLAN[masks]
    out = {h: frozenset(i for i, m in enumerate(masks) if _par(m & h)) for h in range(1, 4096)}
    h1 = min(out, key=lambda h: (len(out[h]), h))
    best = None
    for h2 in sorted(out, key=lambda h: (len(out[h1] & out[h]), h)):
        if h2 == h1:
            continue
        rest = out[h1] & out[h2]
        for h3 in range(1, 4096):
            if h3 != h1 and h3 != h2 and not rest & out[h3]:
                best = (h1, h2, h3)
                break
        if best:
            break
    assert best, 'no three hyperplanes hold all the forms'
    _PLAN[masks] = best
    return best


class Reg:
    """11 wires holding F * rows[k](qi); rows = a basis of the hyperplane {m : par(m & h) = 0}."""

    def __init__(self, b, F, qi, h):
        self.b, self.F, self.qi, self.h = b, F, qi, h
        e = _bits(h)[0]
        J = [j for j in range(12) if j != e]
        self.w = b.alloc_qubits(11)
        for j in J:
            if h >> j & 1:
                b.cx(qi[e], qi[j])
        for k, j in enumerate(J):
            b.init_and(F, qi[j], self.w[k])
        for j in reversed(J):
            if h >> j & 1:
                b.cx(qi[e], qi[j])
        self.rows = [(1 << j) | ((h >> j & 1) << e) for j in J]

    def solve(self, m):
        basis = {}
        for k, r in enumerate(self.rows):
            v, c = r, 1 << k
            while v:
                p = v.bit_length() - 1
                if p not in basis:
                    basis[p] = (v, c)
                    break
                v ^= basis[p][0]
                c ^= basis[p][1]
            else:
                raise AssertionError('register rows dependent')
        v, c = m, 0
        while v:
            p = v.bit_length() - 1
            if p not in basis:
                return None
            v ^= basis[p][0]
            c ^= basis[p][1]
        return _bits(c)

    def mul(self, x, z, terms, done):
        """z ^= x * (F qi): the terms whose form lies in the current hyperplane and are not in `done`."""
        b = self.b
        for n, (mask, val) in enumerate(terms):
            if n in done or _par(mask & self.h):
                continue
            S = self.solve(mask)
            assert S, (mask, self.h)
            ix, iz = _bits(mask), _bits(val)
            p, t, r = ix[0], iz[0], S[0]
            for a in ix[1:]:
                b.cx(x[a], x[p])
            for s in S[1:]:
                b.cx(self.w[s], self.w[r])
            for j in iz[1:]:
                b.cx(z[t], z[j])
            b.ccx(x[p], self.w[r], z[t])
            for j in reversed(iz[1:]):
                b.cx(z[t], z[j])
            for s in reversed(S[1:]):
                b.cx(self.w[s], self.w[r])
            for a in reversed(ix[1:]):
                b.cx(x[a], x[p])
            done.add(n)

    def toggle(self, hnew):
        """One Toffoli: the register's span moves from hyperplane h to hyperplane hnew."""
        b = self.b
        odd = [k for k in range(11) if _par(self.rows[k] & hnew)]
        w0 = odd[0]
        for k in odd[1:]:
            b.cx(self.w[w0], self.w[k])
            self.rows[k] ^= self.rows[w0]
        both = self.h & hnew
        g = (both & -both) if both else ((self.h & -self.h) | (hnew & -hnew))
        ix = _bits(g)
        for a in ix[1:]:
            b.cx(self.qi[a], self.qi[ix[0]])
        b.ccx(self.F, self.qi[ix[0]], self.w[w0])
        for a in reversed(ix[1:]):
            b.cx(self.qi[a], self.qi[ix[0]])
        self.rows[w0] ^= g
        self.h = hnew
        assert not any(_par(r & hnew) for r in self.rows)

    def product(self, x, z, terms, walk):
        done = set()
        assert self.h == walk[0]
        self.mul(x, z, terms, done)
        for h in walk[1:]:
            self.toggle(h)
            self.mul(x, z, terms, done)
        assert len(done) == len(terms)

    def erase(self):
        b = self.b
        for k in range(11):
            ix = _bits(self.rows[k])
            for a in ix[1:]:
                b.cx(self.qi[a], self.qi[ix[0]])
            b.del_and(self.F, self.qi[ix[0]], self.w[k])
            for a in reversed(ix[1:]):
                b.cx(self.qi[a], self.qi[ix[0]])
        b.free(self.w)


def _make_f(b, Z, pole):
    F = b.alloc_qubits(1)[0]
    b.x(pole)
    b.init_and(Z, pole, F)
    b.x(pole)
    return F


def _drop_f(b, Z, pole, F):
    b.x(pole)
    b.del_and(Z, pole, F)
    b.x(pole)
    b.free(F)


def t4_forward(b, ar, p, u, v, pole, qi, Z):
    terms = list(ar.tensor(12))
    h1, h2, h3 = plan(terms)
    F = _make_f(b, Z, pole)
    zm = b.alloc_bits(1)
    b.hmr(Z, zm[0])
    b.free(Z)
    R = Reg(b, F, qi, h1)
    R.product(p, v, terms, (h1, h2, h3))
    R.product(v, u, terms, (h3, h2, h1))
    fm = b.alloc_bits(11)
    b.hmr(R.w, fm)
    b.free(R.w)
    ff = b.alloc_bits(1)
    b.hmr(F, ff[0])
    b.free(F)
    LOG['t4_forward'] = LOG.get('t4_forward', 0) + 1
    LOG['plan'] = (h1, h2, h3)
    return (zm, fm, ff, tuple(R.rows))


def t4_reverse(b, ar, guarded, p, d, u, v, pole, qi, token, mutate=False):
    zm, fm, ff, rows = token
    terms = list(ar.tensor(12))
    h1, h2, h3 = plan(terms)
    Z = guarded.fresh_zero_flag(b, d)
    b.cz(Z, zm[0])
    b.free(zm)
    F = _make_f(b, Z, pole)
    b.cz(F, ff[0])
    b.free(ff)
    for k, r in enumerate(rows):
        ix = _bits(r)
        for a in ix[1:]:
            b.cx(qi[a], qi[ix[0]])
        if not mutate:
            b.ccz(F, qi[ix[0]], fm[k])
        for a in reversed(ix[1:]):
            b.cx(qi[a], qi[ix[0]])
    b.free(fm)
    m2 = b.alloc_bits(1)
    b.hmr(Z, m2[0])
    b.free(Z)
    R = Reg(b, F, qi, h1)
    R.product(v, u, terms, (h1, h2, h3))
    R.product(p, v, terms, (h3, h2, h1))
    R.erase()
    f2 = b.alloc_bits(1)
    b.hmr(F, f2[0])
    b.free(F)
    Z = guarded.fresh_zero_flag(b, d)
    b.cz(Z, m2[0])
    b.free(m2)
    F = _make_f(b, Z, pole)
    b.cz(F, f2[0])
    b.free(f2)
    _drop_f(b, Z, pole, F)
    LOG['t4_reverse'] = LOG.get('t4_reverse', 0) + 1
    return Z


# ------------------------------------------------------------------------------------------------- self-test
class _Mock:
    """Classical bits only: checks the VALUES (not the phases) and counts Toffolis and the widest live set."""

    def __init__(self):
        self.q, self.live, self.peak, self.tof, self.n = {}, 0, 0, 0, 0

    def alloc_qubits(self, n):
        out = []
        for _ in range(n):
            self.n += 1
            self.q[self.n] = 0
            out.append(self.n)
        self.live += n
        self.peak = max(self.peak, self.live)
        return out

    def alloc_bits(self, n):
        return [None] * n

    def free(self, x):
        if isinstance(x, int):
            x = [x]
        for a in x:
            if a is None:
                continue
            assert self.q[a] == 0, 'freed a wire that is not 0'
            del self.q[a]
            self.live -= 1

    def x(self, a):
        self.q[a] ^= 1

    def cx(self, a, t):
        self.q[t] ^= self.q[a]

    def ccx(self, a, c, t):
        self.tof += 1
        self.q[t] ^= self.q[a] & self.q[c]

    def init_and(self, a, c, t):
        assert self.q[t] == 0
        self.tof += 1
        self.q[t] = self.q[a] & self.q[c]

    def del_and(self, a, c, t):
        assert self.q[t] == self.q[a] & self.q[c], 'del_and on a wire that is not the AND'
        self.q[t] = 0

    def hmr(self, a, m):
        for w in (a if isinstance(a, list) else [a]):
            self.q[w] = 0

    def cz(self, *a):
        pass

    def ccz(self, *a):
        pass

    def word(self, value):
        w = self.alloc_qubits(12)
        for i, a in enumerate(w):
            self.q[a] = value >> i & 1
        return w

    def read(self, w):
        return sum(self.q[a] << i for i, a in enumerate(w))


class _Guard:
    @staticmethod
    def fresh_zero_flag(b, d):
        Z = b.alloc_qubits(1)[0]
        ladder = b.alloc_qubits(11)          # the kit's ladder: 11 ANDs, 11 wires beside the flag
        b.tof += 11
        b.q[Z] = int(b.read(d) == 0)
        b.free(ladder)
        return Z


def _test(kit):
    import random
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    sys.path.insert(0, kit)
    from exp.arithmetic import field as ar
    terms = list(ar.tensor(12))

    def mul(x, y):
        r = 0
        for m, val in terms:
            if _par(x & m) & _par(y & m):
                r ^= val
        return r
    print('plan (h1, h2, h3):', plan(terms), ' forms inside each:',
          [sum(1 for m, _ in terms if not _par(m & h)) for h in plan(terms)])
    rnd = random.Random(46)
    stats, fired = set(), 0
    for trial in range(400):
        b = _Mock()
        pv, qv, uv, vv = (rnd.randrange(4096) for _ in range(4))
        dv = 0 if trial % 2 == 0 else rnd.randrange(1, 4096)
        polev = 1 if trial % 8 == 6 else 0
        p, qi, u, v, d = b.word(pv), b.word(qv), b.word(uv), b.word(vv), b.word(dv)
        pole = b.alloc_qubits(1)[0]
        b.q[pole] = polev
        Fv = int(dv == 0 and not polev)
        fired += Fv
        base = b.live
        b.peak = base
        Z = _Guard.fresh_zero_flag(b, d)
        b.tof = 0
        tok = t4_forward(b, ar, p, u, v, pole, qi, Z)
        fwd, fpeak = b.tof, b.peak - base
        v1 = vv ^ (mul(pv, qv) if Fv else 0)
        u1 = uv ^ (mul(v1, qv) if Fv else 0)
        assert (b.read(v), b.read(u), b.read(p), b.read(qi), b.read(d)) == (v1, u1, pv, qv, dv), trial
        assert b.live == base
        b.tof, b.peak = 0, base
        Z = t4_reverse(b, ar, _Guard, p, d, u, v, pole, qi, tok)
        assert (b.read(v), b.read(u), b.read(p), b.read(qi), b.read(d)) == (vv, uv, pv, qv, dv), trial
        assert b.q[Z] == int(dv == 0) and b.live == base + 1
        stats.add((fwd, fpeak, b.tof, b.peak - base))
    print('400 random cases (F = 1 in %d): values exact in both directions' % fired)
    print('(forward Toffolis, forward live beside qi and pole, reverse Toffolis, reverse live):', sorted(stats))


if __name__ == '__main__':
    assert sys.argv[1:2] == ['test'], __doc__
    _test(sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
        '~/qcc/ws/dqi_int/src/source_cap902/playbook/kit/decoder'))
