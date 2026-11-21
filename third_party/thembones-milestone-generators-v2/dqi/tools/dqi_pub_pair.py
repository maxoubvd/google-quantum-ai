"""dqi_pub_pair.py -- builder D1 (holder 4): the payload pair program of R33-DQI-FIELD as a kickmix kernel.

One 12-bit word g times two 12-bit scalars d1, d2, both products XORed into two targets:
    t1 ^= g * d1,   t2 ^= g * d2          (GF(2^12), modulus 0x1053, the kit's polynomial basis)
through the 7 GF(8) points of the kit's 42-term scheme, 10 mixed-form ANDs per point (a SAT program of
tools/r33_dqi-field_quad.py, rule yfree1 / yfree2): 70 CCX per word instead of 84. Each AND is ONE CCX on in-place
parities of (g, d1, d2) with its output spread over (t1, t2) by CNOT conjugation (the kit's mul_xor0 pattern):
no ancilla, no phase. The program leaves a fixed quadratic function of the scalars alone (the garbage); it is the
same for every word the scalars multiply, so it is paid ONCE per (d1, d2) pair: one CCX per basis form of its
span (mod linear functions), its target spread over the targets of ALL the words, plus CNOTs.

  prog = Program(path)                       expand a program file (k3_n2_L10_yfree*.json)
  pair_level(b, prog, words, d1, d2, T1, T2)  T1[i] ^= words[i] * d1, T2[i] ^= words[i] * d2 for all i
  prog.counts(n)                             CCX for n words (70 n + cache)
Registers are lists of 12 qubits (bit i = coefficient of x^i). Inputs and targets must be pairwise disjoint.
"""
import collections
import importlib.util
import json
import os

TOOLS = os.path.dirname(os.path.abspath(__file__))


def _designer():
    """The designer's checker (read-only): its 7-point GF(8) scheme build_scheme() -> (E, J)."""
    spec = importlib.util.spec_from_file_location('r33_dqi_field_prog', os.path.join(TOOLS, 'r33_dqi-field_prog.py'))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _bits(x):
    i = 0
    while x:
        if x & 1:
            yield i
        x >>= 1
        i += 1


def _gauss_basis(vecs):
    """Indices of a basis subset of vecs (ints over GF(2)) and, for each vec, its combination of basis indices."""
    piv = {}            # pivot bit -> (vector, combination mask over basis positions)
    basis = []
    combos = []
    for idx, v in enumerate(vecs):
        c = 0
        w = v
        for pb in sorted(piv, reverse=True):
            if (w >> pb) & 1:
                w ^= piv[pb][0]
                c ^= piv[pb][1]
        if w:
            pos = len(basis)
            basis.append(idx)
            piv[w.bit_length() - 1] = (w, c ^ (1 << pos))
            combos.append(1 << pos)
        else:
            combos.append(c)
    return basis, combos


def _express(forms, q):
    """Indices S with XOR of forms[S] == q (forms independent; q in their span)."""
    piv = {}
    for idx, v in enumerate(forms):
        c = 1 << idx
        for pb in sorted(piv, reverse=True):
            if (v >> pb) & 1:
                v ^= piv[pb][0]
                c ^= piv[pb][1]
        assert v
        piv[v.bit_length() - 1] = (v, c)
    w, c = q, 0
    for pb in sorted(piv, reverse=True):
        if (w >> pb) & 1:
            w ^= piv[pb][0]
            c ^= piv[pb][1]
    assert w == 0, 'garbage outside the cache span'
    return list(_bits(c))


class Program:
    """A per-point program expanded to gates on the 36 input bits (g 0-11, d1 12-23, d2 24-35) and the 24 target
    bits (t1 0-11, t2 12-23)."""

    def __init__(self, path):
        P = json.load(open(path))
        assert P['k'] == 3 and P['n'] == 2 and P['result'] == 'SAT', path
        self.path = path
        L = P['L']
        nv = 9
        l = [[int(ch) for ch in row] for row in P['program']['l']]
        lp = [[int(ch) for ch in row] for row in P['program']['lp']]
        W = [[int(ch) for ch in row] for row in P['program']['W']]
        lin = [[0] * nv for _ in range(6)]
        quad = [dict() for _ in range(6)]
        for c in range(6):
            for i in range(L):
                if not W[c][i]:
                    continue
                for a in range(nv):
                    if l[i][a] & lp[i][a]:
                        lin[c][a] ^= 1
                    for bb in range(a + 1, nv):
                        if (l[i][a] & lp[i][bb]) ^ (l[i][bb] & lp[i][a]):
                            quad[c][(a, bb)] = quad[c].get((a, bb), 0) ^ 1
        garbage = [sorted(p for p, v in quad[c].items() if v and p[0] >= 3) for c in range(6)]
        assert not [p for c in range(6) for p, v in quad[c].items() if v and p[1] < 3], 'word garbage'
        for c in range(6):
            for (a, bb) in garbage[c]:
                assert (a - 3) // 3 == (bb - 3) // 3, 'cross-scalar garbage'
        D = _designer()
        E, J = D.build_scheme()
        self.E, self.J = E, J

        def vmask(p, q):                       # 36-bit input mask of point bit q (0..8) at point p
            if q < 3:
                return E[p][q]
            s, k = (q - 3) // 3, (q - 3) % 3
            return E[p][k] << (12 * (s + 1))

        def outvec(p, c):                      # 24-bit target vector of output bit c (r = c // 3, j = c % 3)
            r, j = c // 3, c % 3
            return J[p][j] << (12 * r)
        self.ands = []                         # (formA36, formB36, out24)
        lin_m = [0] * 24                       # target bit e ^= parity(inputs & lin_m[e])
        gmono = {}                             # (s, formA12, formB12) -> out24 (garbage monomials, per word)
        for p in range(7):
            for i in range(L):
                fa = 0
                fb = 0
                for q in range(nv):
                    if l[i][q]:
                        fa ^= vmask(p, q)
                    if lp[i][q]:
                        fb ^= vmask(p, q)
                out = 0
                for c in range(6):
                    if W[c][i]:
                        out ^= outvec(p, c)
                self.ands.append((fa, fb, out))
            for c in range(6):
                ov = outvec(p, c)
                for q in range(nv):
                    if lin[c][q]:
                        m = vmask(p, q)
                        for e in _bits(ov):
                            lin_m[e] ^= m
                for (a, bb) in garbage[c]:
                    s = (a - 3) // 3
                    key = (s, E[p][(a - 3) % 3], E[p][(bb - 3) % 3], p)
                    gmono[key] = gmono.get(key, 0) ^ ov
        self._finish(lin_m, gmono)

    def _finish(self, lin_m, gmono):
        """Shared tail: simplify the word gates (self.ands), then plan the scalar cache from the garbage monomials
        gmono[(s, formA12, formB12)] = out24."""
        # simplify the word gates: drop zero ones, fold linear ones (A == B) into lin_m
        kept = []
        for fa, fb, out in self.ands:
            if not out or not fa or not fb:
                continue
            if fa == fb:
                for e in _bits(out):
                    lin_m[e] ^= fa
                continue
            kept.append((fa, fb, out))
        self.ands = kept
        self.lin = lin_m
        # cache plan, per scalar s: the garbage each target bit e needs is a quadratic function G_e(d_s) = Q_e + L_e.
        # Its quadratic parts span a space N_s; a basis of N_s made of rank-2 forms u ^ v costs one AND each
        # (u(d) v(d) has quadratic part u ^ v and linear part u & v); the rest is CNOTs.
        self.cache = []                        # (s, formA12, formB12, out24 per word)
        self.cache_lin = [0] * 24              # per word: target bit e ^= parity(d-bits (24) & cache_lin[e])
        self.cache_stats = {}
        mode = os.environ.get('DQI_PUB_PAIR_CACHE', 'rank2')
        # keys (s, formA12, formB12[, place]): with a place tag the plan is made per (scalar, place), whose local
        # garbage space is small and spanned by rank-2 forms more often than the global one
        for s in (0, 1):
            groups = collections.OrderedDict()
            for k, v in gmono.items():
                if k[0] == s and v:
                    groups.setdefault(k[3] if len(k) > 3 else None, []).append((k[:3], v))
            n_mono = n_and = 0
            qall = [0] * 24
            for monos in groups.values():
                qe = [0] * 24
                le = [0] * 24
                for (_, fa, fb), ov in monos:
                    qm, lm = self._qpart(fa, fb), fa & fb
                    for e in _bits(ov):
                        qe[e] ^= qm
                        le[e] ^= lm
                for e in range(24):
                    qall[e] ^= qe[e]
                plan = self._plan_rank2(qe) if mode == 'rank2' else None
                if plan is None:               # fall back: the monomials themselves (one AND each)
                    plan = [(fa, fb) for (_, fa, fb), _ in monos]
                forms = [self._qpart(u, v) for u, v in plan]
                basis, _ = _gauss_basis(forms)
                plan = [plan[i] for i in basis]
                forms = [forms[i] for i in basis]
                outs = [0] * len(plan)
                for e in range(24):
                    if not qe[e] and not le[e]:
                        continue
                    S = _express(forms, qe[e])
                    lin_e = le[e]
                    for k in S:
                        outs[k] |= 1 << e
                        lin_e ^= plan[k][0] & plan[k][1]
                    self.cache_lin[e] ^= lin_e << (12 * s)
                for k, (u, v) in enumerate(plan):
                    if outs[k]:
                        self.cache.append((s, u, v, outs[k]))
                        n_and += 1
                n_mono += len(monos)
            self.cache_stats['scalar%d' % s] = dict(monomials=n_mono, needed_rank=len(_gauss_basis(qall)[0]),
                                                    ands=n_and, groups=len(groups))

    @staticmethod
    def _rank_alt(q):
        """GF(2) rank and two independent rows of the alternating 12 x 12 matrix of a quadratic part."""
        rows = [0] * 12
        for i in range(12):
            for j in range(i + 1, 12):
                if (q >> (i * 12 + j)) & 1:
                    rows[i] |= 1 << j
                    rows[j] |= 1 << i
        basis, _ = _gauss_basis(rows)
        return len(basis), [rows[i] for i in basis[:2]]

    def _plan_rank2(self, qe):
        """A basis of span(qe) made of rank-2 forms, as (u, v) pairs; None if the rank-2 elements do not span it."""
        gens = [qe[i] for i in _gauss_basis(qe)[0]]
        r = len(gens)
        if r == 0 or r > 16:
            return None
        cand = []
        for mask in range(1, 1 << r):
            q = 0
            for i in _bits(mask):
                q ^= gens[i]
            rk, rows = self._rank_alt(q)
            if rk == 2:
                u, v = rows
                assert self._qpart(u, v) == q
                cand.append((u, v))
        forms = [self._qpart(u, v) for u, v in cand]
        basis, _ = _gauss_basis(forms)
        if len(basis) < r:
            return None
        return [cand[i] for i in basis]

    @staticmethod
    def _qpart(fa, fb):
        """Quadratic part of parity(x & fa) * parity(x & fb) on 12 bits: bit (i, j), i < j, as one int."""
        v = 0
        for i in range(12):
            for j in range(i + 1, 12):
                if (((fa >> i) & (fb >> j)) ^ ((fa >> j) & (fb >> i))) & 1:
                    v |= 1 << (i * 12 + j)
        return v

    @staticmethod
    def _lpart(fa, fb):
        return fa & fb

    def counts(self, n):
        return dict(per_word=len(self.ands), cache=len(self.cache), total=n * len(self.ands) + len(self.cache),
                    bilinear=84 * n)


class MixedProgram(Program):
    """One word times two scalars through GF(2^12) = GF(4)[t]/(f), deg f = 6: the product of two degree-5
    polynomials over GF(4) is fixed by 11 conditions, taken at the 5 rational places (t = 0, 1, w, w^2, infinity:
    GF(4) products) and at 3 places of degree 2 (GF(16) values). Per word: 5 GF(4) pair programs (one word x two
    cached scalars, 5 gates each; rule yfree1) + 3 GF(16) places, each two copies of a two-sided 5-gate product
    (rule xyfree) plus the ANDs that cancel the word's garbage, shared by both copies: 5 x 5 + 3 x (10 + G).
    Same output format as Program (ands / lin / cache / cache_lin)."""

    def __init__(self, prog4_path, prog16_path):
        P4 = json.load(open(prog4_path))
        P16 = json.load(open(prog16_path))
        assert P4['k'] == 2 and P4['n'] == 2 and P4['rule'] == 'yfree1' and P4['result'] == 'SAT'
        # D1 holder 5: a GF(16) program may also be a DIRECT pair program (n = 2, rule yfree / yfree1: one word x two
        # cached scalars), used at the GF(16) places like the GF(4) program at the rational places.
        direct16 = P16['k'] == 4 and P16['n'] == 2 and P16['rule'] in ('yfree', 'yfree1') and P16['result'] == 'SAT'
        assert direct16 or (P16['k'] == 4 and P16['n'] == 1 and P16['rule'] == 'xyfree' and P16['result'] == 'SAT')
        self.path = '%s+%s' % (os.path.basename(prog4_path), os.path.basename(prog16_path))
        D = _designer()
        mul = D.mul

        def pw(a, e):
            r = 1
            for _ in range(e):
                r = mul(r, a)
            return r
        om = next(v for v in range(2, 4096) if mul(v, v) ^ v ^ 1 == 0)            # GF(4) = {0, 1, w, w^2}
        beta = next(v for v in range(2, 4096) if pw(v, 4) ^ v ^ 1 == 0)          # GF(16) = F2[beta], beta^4 + beta + 1
        g4 = [0, 1, om, mul(om, om)]
        c4 = {g4[i]: i for i in range(4)}                                        # 2-bit coords in (1, w)
        g16 = {}
        for i in range(16):
            v = 0
            for j in range(4):
                if (i >> j) & 1:
                    v ^= pw(beta, j)
            g16[v] = i                                                           # 4-bit coords in (1, b, b^2, b^3)
        theta = None
        for cand in range(2, 4096):
            basis = [mul(g4[1 + j], pw(cand, d)) for d in range(6) for j in range(2)]
            if len(_gauss_basis(basis)[0]) == 12:
                theta = cand
                break
        basis = [mul(g4[1 + j], pw(theta, d)) for d in range(6) for j in range(2)]  # index 2 d + j
        piv = {}
        for i, v in enumerate(basis):
            c = 1 << i
            for pb in sorted(piv, reverse=True):
                if (v >> pb) & 1:
                    v ^= piv[pb][0]
                    c ^= piv[pb][1]
            piv[v.bit_length() - 1] = (v, c)

        def digits(x):                       # x -> 6 GF(4) digits (as field elements)
            w, c = x, 0
            for pb in sorted(piv, reverse=True):
                if (w >> pb) & 1:
                    w ^= piv[pb][0]
                    c ^= piv[pb][1]
            assert w == 0
            return [g4[(c >> (2 * d)) & 3] for d in range(6)]
        alphas = []
        seen_cls = set()
        for v in sorted(g16):
            if v in c4:
                continue
            cls = frozenset((v, pw(v, 4)))
            if cls in seen_cls:
                continue
            seen_cls.add(cls)
            alphas.append(v)
            if len(alphas) == 3:
                break
        places = [('r', t) for t in g4] + [('inf', None)] + [('q', a) for a in alphas]

        def local(x, place):                 # local value bits of x at a place
            dg = digits(x)
            kind, t = place
            if kind == 'inf':
                return c4[dg[5]]
            val, tp = 0, 1
            for d in range(6):
                val ^= mul(dg[d], tp)
                tp = mul(tp, t)
            return c4[val] if kind == 'r' else g16[val]
        nb = [2 if k in ('r', 'inf') else 4 for k, _ in places]
        Emask = []                           # per place: list of 12-bit masks, one per local bit
        for pi, pl in enumerate(places):
            ms = [0] * nb[pi]
            for e in range(12):
                v = local(1 << e, pl)
                for j in range(nb[pi]):
                    if (v >> j) & 1:
                        ms[j] |= 1 << e
            Emask.append(ms)
        # local products in coordinates (must match the SAT programs' fields: GF(4) mod x^2+x+1, GF(16) mod x^4+x+1)
        from_c16 = {i: v for v, i in g16.items()}

        def lprod(place, xv, yv):
            kind = place[0]
            if kind in ('r', 'inf'):
                return c4[mul(g4[xv], g4[yv])]
            return g16[mul(from_c16[xv], from_c16[yv])]
        # reconstruction: 12-bit product = L . (local product bits, 22 of them), solved and verified
        import random
        rng = random.Random(5)
        offs = [sum(nb[:i]) for i in range(len(places))]
        NL = sum(nb)
        rows = []
        for _ in range(400):
            x, y = rng.randrange(4096), rng.randrange(4096)
            v = 0
            for pi, pl in enumerate(places):
                v |= lprod(pl, local(x, pl), local(y, pl)) << offs[pi]
            rows.append((v, mul(x, y)))
        pv = {}
        for v, pay in rows:
            for pb in sorted(pv, reverse=True):
                if (v >> pb) & 1:
                    v ^= pv[pb][0]
                    pay ^= pv[pb][1]
            if v:
                pv[v.bit_length() - 1] = (v, pay)
            else:
                assert pay == 0, 'local products do not determine the product'
        assert len(pv) == NL, (len(pv), NL)
        # back-substitute: Lcol[b] = contribution of local bit b
        Lcol = [0] * NL
        for pb in sorted(pv):
            v, pay = pv[pb]
            for b2 in _bits(v ^ (1 << pb)):
                pay ^= Lcol[b2]
            Lcol[pb] = pay
        for _ in range(2000):
            x, y = rng.randrange(4096), rng.randrange(4096)
            acc = 0
            for pi, pl in enumerate(places):
                lp_ = lprod(pl, local(x, pl), local(y, pl))
                for j in range(nb[pi]):
                    if (lp_ >> j) & 1:
                        acc ^= Lcol[offs[pi] + j]
            assert acc == mul(x, y), 'reconstruction wrong'
        self.ands = []
        lin_m = [0] * 24
        gmono = {}
        wg_stats = []

        def parse(P):
            L = P['L']
            l = [[int(ch) for ch in row] for row in P['program']['l']]
            lp = [[int(ch) for ch in row] for row in P['program']['lp']]
            W = [[int(ch) for ch in row] for row in P['program']['W']]
            return L, l, lp, W

        def local_terms(L, l, lp, W, nout, nv):
            lin = [[0] * nv for _ in range(nout)]
            quad = [dict() for _ in range(nout)]
            for c in range(nout):
                for i in range(L):
                    if not W[c][i]:
                        continue
                    for a in range(nv):
                        if l[i][a] & lp[i][a]:
                            lin[c][a] ^= 1
                        for bb in range(a + 1, nv):
                            if (l[i][a] & lp[i][bb]) ^ (l[i][bb] & lp[i][a]):
                                quad[c][(a, bb)] = quad[c].get((a, bb), 0) ^ 1
            return lin, quad
        L4, l4, lp4, W4 = parse(P4)
        lin4, quad4 = local_terms(L4, l4, lp4, W4, 4, 6)
        L16, l16, lp16, W16 = parse(P16)
        if direct16:
            lin16, quad16 = local_terms(L16, l16, lp16, W16, 8, 12)
        else:
            lin16, quad16 = local_terms(L16, l16, lp16, W16, 4, 8)
        for pi, pl in enumerate(places):
            em = Emask[pi]
            if nb[pi] == 4 and direct16:
                # direct GF(16) pair program: vars X 0-3 (word), Y1 4-7, Y2 8-11; outputs c = 4 r + bit
                def gmask(q, em=em):
                    if q < 4:
                        return em[q]
                    s, kk = (q - 4) // 4, (q - 4) % 4
                    return em[kk] << (12 * (s + 1))

                def ovec(c, pi=pi):
                    r, bit = c // 4, c % 4
                    return Lcol[offs[pi] + bit] << (12 * r)
                for i in range(L16):
                    fa = fb = 0
                    for q in range(12):
                        if l16[i][q]:
                            fa ^= gmask(q)
                        if lp16[i][q]:
                            fb ^= gmask(q)
                    out = 0
                    for c in range(8):
                        if W16[c][i]:
                            out ^= ovec(c)
                    self.ands.append((fa, fb, out))
                for c in range(8):
                    ov = ovec(c)
                    for q in range(12):
                        if lin16[c][q]:
                            for e in _bits(ov):
                                lin_m[e] ^= gmask(q)
                    for (a, bb), val in quad16[c].items():
                        if not val or (a < 4 <= bb):            # word x scalar: the product itself
                            continue
                        assert a >= 4, 'GF(16) pair program: word garbage left'
                        sa, sb = (a - 4) // 4, (bb - 4) // 4
                        assert sa == sb or P16['rule'] == 'yfree', 'GF(16) pair program: cross-scalar garbage'
                        if sa != sb:
                            raise NotImplementedError('cross-scalar cache (rule yfree) at GF(16) places')
                        key = (sa, em[(a - 4) % 4], em[(bb - 4) % 4], pi)
                        gmono[key] = gmono.get(key, 0) ^ ov
                wg_stats.append(0)
                continue
            if nb[pi] == 2:
                # vars: X 0-1 (word), Y1 2-3, Y2 4-5; outputs c = 2 r + bit
                def gmask(q):
                    if q < 2:
                        return em[q]
                    s, k = (q - 2) // 2, (q - 2) % 2
                    return em[k] << (12 * (s + 1))

                def ovec(c):
                    r, bit = c // 2, c % 2
                    return Lcol[offs[pi] + bit] << (12 * r)
                for i in range(L4):
                    fa = fb = 0
                    for q in range(6):
                        if l4[i][q]:
                            fa ^= gmask(q)
                        if lp4[i][q]:
                            fb ^= gmask(q)
                    out = 0
                    for c in range(4):
                        if W4[c][i]:
                            out ^= ovec(c)
                    self.ands.append((fa, fb, out))
                for c in range(4):
                    ov = ovec(c)
                    for q in range(6):
                        if lin4[c][q]:
                            for e in _bits(ov):
                                lin_m[e] ^= gmask(q)
                    for (a, bb), val in quad4[c].items():
                        if not val or (a < 2 <= bb):            # word x scalar: the product itself
                            continue
                        assert a >= 2 and (a - 2) // 2 == (bb - 2) // 2, 'GF(4) program: unexpected garbage'
                        s = (a - 2) // 2
                        key = (s, em[(a - 2) % 2], em[(bb - 2) % 2], pi)
                        gmono[key] = gmono.get(key, 0) ^ ov
            else:
                # two copies of the two-sided GF(16) product: vars x 0-3 = the word X, y 4-7 = scalar r
                wgarb = {}                   # (formA12, formB12) on the word -> out24 (both copies)
                for r in (0, 1):
                    def gmask(q, r=r):
                        if q < 4:
                            return em[q]
                        return em[q - 4] << (12 * (r + 1))

                    def ovec(c, r=r):
                        return Lcol[offs[pi] + c] << (12 * r)
                    for i in range(L16):
                        fa = fb = 0
                        for q in range(8):
                            if l16[i][q]:
                                fa ^= gmask(q)
                            if lp16[i][q]:
                                fb ^= gmask(q)
                        out = 0
                        for c in range(4):
                            if W16[c][i]:
                                out ^= ovec(c)
                        self.ands.append((fa, fb, out))
                    for c in range(4):
                        ov = ovec(c)
                        for q in range(8):
                            if lin16[c][q]:
                                for e in _bits(ov):
                                    lin_m[e] ^= gmask(q)
                        for (a, bb), val in quad16[c].items():
                            if not val:
                                continue
                            if a < 4 and bb < 4:            # the word's garbage: cancelled per word below
                                key = (em[a], em[bb])
                                wgarb[key] = wgarb.get(key, 0) ^ ov
                            elif a >= 4 and bb >= 4:        # the scalar's garbage: cached
                                key = (r, em[a - 4], em[bb - 4], pi)
                                gmono[key] = gmono.get(key, 0) ^ ov
                            # a < 4 <= bb: word x scalar, the product itself
                # cancel the word garbage: needed functions per target bit, basis of rank-2 forms (one AND each)
                qe = [0] * 24
                le = [0] * 24
                for (fa, fb), ov in wgarb.items():
                    qm = self._qpart(fa, fb)
                    for e in _bits(ov):
                        qe[e] ^= qm
                        le[e] ^= fa & fb
                plan = self._plan_rank2(qe)
                if plan is None:
                    plan = list(wgarb.keys())
                forms = [self._qpart(u, v) for u, v in plan]
                bsel, _ = _gauss_basis(forms)
                plan = [plan[i] for i in bsel]
                forms = [forms[i] for i in bsel]
                outs = [0] * len(plan)
                for e in range(24):
                    if not qe[e] and not le[e]:
                        continue
                    S = _express(forms, qe[e])
                    lin_e = le[e]
                    for k in S:
                        outs[k] |= 1 << e
                        lin_e ^= plan[k][0] & plan[k][1]
                    lin_m[e] ^= lin_e                       # word bits (0-11) of the 36-bit input
                for k, (u, v) in enumerate(plan):
                    if outs[k]:
                        self.ands.append((u, v, outs[k]))
                wg_stats.append(len(plan))
        self.word_garbage_ands = wg_stats
        self._finish(lin_m, gmono)

def _apply_pair_ccx(b, inq, fa, fb, targets):
    """CCX(parity(inq & fa), parity(inq & fb)) XORed into every qubit of `targets`, zero workspace."""
    n = len(inq)
    content = [1 << k for k in range(n)]           # what qubit k holds, as a mask over the original bits
    done = []                                      # CNOTs applied to the inputs (control, target)

    def cnot(c, t):
        b.cx(inq[c], inq[t])
        content[t] ^= content[c]
        done.append((c, t))
    sa = list(_bits(fa))
    pa = sa[0]
    for k in sa[1:]:
        cnot(k, pa)
    # coordinates of fb in the current basis
    S = _coords(content, fb, n)
    cand = [k for k in S if k != pa]
    assert cand, 'dependent forms'
    pb = cand[0]
    for k in S:
        if k != pb:
            cnot(k, pb)
    t0 = targets[0]
    for t in targets[1:]:
        b.cx(t0, t)
    b.ccx(inq[pa], inq[pb], t0)
    for t in reversed(targets[1:]):
        b.cx(t0, t)
    for c, t in reversed(done):
        b.cx(inq[c], inq[t])


def _coords(content, f, n):
    """Set S of qubit indices whose contents XOR to f (contents form a basis)."""
    rows = [(content[k], 1 << k) for k in range(n)]
    piv = {}
    for v, c in rows:
        for pb in sorted(piv, reverse=True):
            if (v >> pb) & 1:
                v ^= piv[pb][0]
                c ^= piv[pb][1]
        if v:
            piv[v.bit_length() - 1] = (v, c)
    w, c = f, 0
    for pb in sorted(piv, reverse=True):
        if (w >> pb) & 1:
            w ^= piv[pb][0]
            c ^= piv[pb][1]
    assert w == 0
    return list(_bits(c))


def _xorset(lists):
    """Qubits that occur an odd number of times in the given lists (XOR fan-out sets), in first-seen order."""
    cnt = {}
    first = {}
    for lst in lists:
        for q in lst:
            k = str(q)
            cnt[k] = cnt.get(k, 0) ^ 1
            first.setdefault(k, q)
    return [first[k] for k in first if cnt[k]]


def emit_word_map(b, prog, w, d1, d2, tset):
    """Word part for one instance; tset[e] (e = 0..23) = the physical qubits that receive target bit e (bits 0-11:
    the product with d1, 12-23: with d2), any XOR-linear spreading of the two products."""
    inq = list(w) + list(d1) + list(d2)
    for fa, fb, out in prog.ands:
        targets = _xorset([tset[e] for e in _bits(out)])
        if targets:
            _apply_pair_ccx(b, inq, fa, fb, targets)
    for e, m in enumerate(prog.lin):
        if m and tset[e]:
            for k in _bits(m):
                for q in tset[e]:
                    b.cx(inq[k], q)


def emit_cache_map(b, prog, d1, d2, tsets):
    """The scalar garbage of ALL instances at once (one CCX per cache form, fanned out over every instance)."""
    dq = list(d1) + list(d2)
    for s, fa, fb, out in prog.cache:
        targets = _xorset([ts[e] for ts in tsets for e in _bits(out)])
        if targets:
            _apply_pair_ccx(b, dq, fa << (12 * s), fb << (12 * s), targets)
    for e, m in enumerate(prog.cache_lin):
        if m:
            qs = _xorset([ts[e] for ts in tsets])
            for k in _bits(m):
                for q in qs:
                    b.cx(dq[k], q)


def _plain_tset(t1, t2):
    return [[t1[e]] for e in range(12)] + [[t2[e]] for e in range(12)]


def emit_word(b, prog, g, d1, d2, t1, t2):
    """t1 ^= g * d1 + G1(d1, d2), t2 ^= g * d2 + G2(d1, d2): the word part (len(prog.ands) CCX + CNOTs)."""
    emit_word_map(b, prog, g, d1, d2, _plain_tset(t1, t2))


def emit_cache(b, prog, d1, d2, T1, T2):
    """The scalar garbage for ALL words at once: T1[i], T2[i] (lists of target words) get G1, G2."""
    emit_cache_map(b, prog, d1, d2, [_plain_tset(a, c) for a, c in zip(T1, T2)])


# ---- the shifted payload product of the converter, both payloads at once (switch DQI_PUB_PAIR_SHIFT=1) -----------
# The kit's dk.shifted_product(b, g, x, out, F), called once per payload: the head x is split by the gap flag F into
# (x, z) = ((1-F) h, F h) and out ^= g * (x + w z) by evaluation at the points t = 0..N-1 (N = len(g)), so the scalar
# at point t is h * c_F(t) with c_F(t) = V[t][0] (F = 0) or V[t][1] (F = 1). Moving c_F(t) to the WORD side (12 ANDs
# per point: z_t = c_F(t) g(t), none where V[t][0] == V[t][1]) leaves both heads fixed over the whole level, so the
# pair program multiplies z_t by (h1, h2) with one cache per level: 82 n + 28 - 12 (t = 1) against 84 n + 24.

def shifted_pair(b, prog, K, g, x1, out1, x2, out2, F):
    """out1 ^= g * c_F * x1 and out2 ^= g * c_F * x2 as the two kit calls dk.shifted_product would (K: the kit's
    exp.convert.kernels module: vandermonde, expose, unexpose, scale_const, mul, inv)."""
    N = len(g)
    V, W = K.vandermonde(N)
    gf_mul = K.mul
    pts = [t for t in range(N) if any(W[j][t] for j in range(max(len(out1), len(out2))))]

    def tset_for(t):
        ts = []
        for out in (out1, out2):
            for e in range(12):
                qs = []
                for j in range(len(out)):
                    c = W[j][t] if j < len(out) else 0
                    if c:
                        v = gf_mul(c, 1 << e)
                        qs.extend(out[j][i] for i in range(12) if (v >> i) & 1)
                ts.append(_xorset([qs]))
        return ts
    tsets = {t: tset_for(t) for t in pts}
    # Allocation parity with the two kit calls (kickmix reuses freed wires last-in first-out, and the decode's
    # final routing wants exact wires): the first payload's call took and gave back 12 qubits and 12 bits; replay
    # that cycle without gates, then do the whole pair inside the second cycle (one z, one mask for every point).
    z0 = b.alloc_qubits(12)
    m0 = b.alloc_bits(12)
    b.free(z0)
    b.free(m0)
    emit_cache_map(b, prog, x1, x2, [tsets[t] for t in pts])
    z = b.alloc_qubits(12)
    mask = b.alloc_bits(12)
    for t in pts:
        a, c = V[t][0], (V[t][1] if N > 1 else 0)
        xg, p = K.expose(b, g, V[t][:N])
        s = a ^ c
        if s == 0:                                  # c_F(t) = a for both values of F
            if a != 1:
                K.scale_const(b, xg, a)
            emit_word_map(b, prog, xg, x1, x2, tsets[t])
            if a != 1:
                K.scale_const(b, xg, K.inv(a))
        else:                                       # z = a g(t) + F s g(t)
            K.scale_const(b, xg, s)
            for j in range(12):
                b.init_and(F, xg[j], z[j])
            K.scale_const(b, xg, K.inv(s))
            if a:
                if a != 1:
                    K.scale_const(b, xg, a)
                b.cx(xg, z)
                if a != 1:
                    K.scale_const(b, xg, K.inv(a))
            emit_word_map(b, prog, z, x1, x2, tsets[t])
            b.hmr(z, mask)                          # z back to 0; repay (-1)^(mask . z)
            if a:
                if a != 1:
                    K.scale_const(b, xg, a)
                b.cz(xg, mask)
                if a != 1:
                    K.scale_const(b, xg, K.inv(a))
            K.scale_const(b, xg, s)
            for j in range(12):
                b.ccz(F, xg[j], mask[j])
            K.scale_const(b, xg, K.inv(s))
        K.unexpose(b, g, V[t][:N], p)
    b.free(z)
    b.free(mask)


PAIRSP = dict(pending=None, pairs=0, prog=None, K=None)


def _shifted_dispatch(b, g, x, out, F):
    """Replaces dk.shifted_product: the first payload's call is held, the second emits both (they are consecutive
    at every call site: `for v in payloads: dk.shifted_product(...)` and root_action's consecutive 'u' ops)."""
    pend = PAIRSP['pending']
    if pend is None:
        PAIRSP['pending'] = (b, list(g), x, out, F)
        return
    b0, g0, x0, out0, F0 = pend
    PAIRSP['pending'] = None
    same = (b0 is b and len(g0) == len(g) and all(str(p) == str(q) for gp, gq in zip(g0, g) for p, q in zip(gp, gq))
            and str(F0) == str(F))
    if not same:
        raise RuntimeError('pub_pair: unpaired shifted_product calls')
    shifted_pair(b, PAIRSP['prog'], PAIRSP['K'], g, x0, out0, x, out, F)
    PAIRSP['pairs'] += 1


def install_shift(prog, decode=None):
    """Point every direct_kernels.shifted_product function object (module attributes, clones in other scopes and
    captured references: its __code__ is replaced) at the dispatcher. Returns the number of objects patched."""
    import sys
    import types
    cands = []
    for m in list(sys.modules.values()):
        d = getattr(m, '__dict__', None)
        if isinstance(d, dict) and isinstance(d.get('shifted_product'), types.FunctionType):
            cands.append(d['shifted_product'])
    walk = sys.modules.get('walk_bound')
    if walk is not None and decode is not None:
        cands.extend(f for f in walk.functions(decode) if getattr(f, '__name__', '') == 'shifted_product')
    K = None                     # the kit's exp/convert/kernels.py (every copy has the same code)
    for m in list(sys.modules.values()):
        d = getattr(m, '__dict__', None)
        if isinstance(d, dict) and all(callable(d.get(nm)) for nm in
                                       ('vandermonde', 'expose', 'unexpose', 'scale_const', 'rank_update', 'mul', 'inv')):
            if getattr(m, '__file__', '') and m.__file__.endswith(os.path.join('convert', 'kernels.py')):
                K = m
                break
    assert K is not None, 'kit kernels module not found'
    PAIRSP['K'] = K
    done = 0
    seen = set()
    for f in cands:
        if id(f) in seen or not f.__code__.co_filename.endswith('direct_kernels.py'):
            continue
        if f.__code__.co_name != 'shifted_product' or f.__closure__:
            continue
        seen.add(id(f))
        G = f.__globals__
        G['_D1_SHIFT_DISPATCH'] = _shifted_dispatch
        stub = {}
        exec('def shifted_product(b,g,x,out,F):\n return _D1_SHIFT_DISPATCH(b,g,x,out,F)\n', G, stub)
        f.__code__ = stub['shifted_product'].__code__
        done += 1
    PAIRSP['prog'] = prog
    import atexit
    atexit.register(lambda: print('dqi_pub_pair shift: pairs emitted', PAIRSP['pairs'], 'pending',
                                  PAIRSP['pending'] is not None, flush=True))
    return done


def pair_level(b, prog, words, d1, d2, T1, T2):
    """T1[i] ^= words[i] * d1 and T2[i] ^= words[i] * d2 for every i (exact; self-inverse)."""
    assert len(words) == len(T1) == len(T2)
    emit_cache(b, prog, d1, d2, T1, T2)
    for g, a, c in zip(words, T1, T2):
        emit_word(b, prog, g, d1, d2, a, c)


DEFAULT_PROG = os.path.expanduser(os.environ.get('DQI_PUB_PAIR_PROG',
                                                 '~/qcc/ws/dqi_pub/pair/progs/k3_n2_L10_yfree2_0_1_2.json'))


def load_program(spec=None):
    """A program file (GF(8) points, Program) or 'mixed:<gf4 yfree1 file>,<gf16 xyfree file>' (MixedProgram)."""
    spec = spec or DEFAULT_PROG
    if spec.startswith('mixed:'):
        a, b = spec[len('mixed:'):].split(',')
        return MixedProgram(os.path.expanduser(a), os.path.expanduser(b))
    return Program(os.path.expanduser(spec))

# ---- install: the O pair of the converter core (patch token pub_pair, after a3 / var_sync3 / pub_zt / var_ztp) ----
# In my a3 core (dqi_pub_a3.CORE_SRC_V5) the forward runs `w[i] ^= top * O[i]` (17 products) and then _bez_monic's
# first reduce `E[i] ^= E[-1] * O[i]` (17); the reverse ends with _bez_monic's last reduce and then the same top loop.
# O, top and E[-1] do not change between the two loops, so both become one pair_level on the 17 words of O.
_BEZ_SRC = '''def _bez_monic(b, O, E, payloads, reverse, history_scratch12, phase):
    def reduce():
        for i in range(len(O)):
            ar.mul_xor0(b, E[-1], O[i], E[i])
    reduce()
    look_root_cached(b, O, E[:-1], payloads, reverse, phase)
    reduce()
'''


LOG = {}

_FWD_OLD = '        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])\n        bez(False)\n'
_FWD_NEW = ('        if phase[\'_short_monic\']:\n'
            '            _D1_PAIR_LEVEL(b,O,top,E[-1],w,E[:17])\n'
            '            _D1_BEZ(b,O,E,[u,w],False,history_scratch12,phase)\n'
            '        else:\n'
            '            for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])\n'
            '            bez(False)\n')
_REV_OLD = '        bez(True)\n        for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])\n'
_REV_NEW = ('        if phase[\'_short_monic\']:\n'
            '            _D1_BEZ(b,O,E,[u,w],True,history_scratch12,phase)\n'
            '            _D1_PAIR_LEVEL(b,O,top,E[-1],w,E[:17])\n'
            '        else:\n'
            '            bez(True)\n'
            '            for i in range(17):old.ar.mul_xor0(b,top,O[i],w[i])\n')


_DELTA_OLD = '        for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],delta)\n        b.cx(delta,top)\n'
_DELTA_NEW = ('        if phase[\'_short_monic\'] and _D1_DELTA is not None:_D1_DELTA(b,O,E,A,delta)\n'
              '        else:\n'
              '            for j in _D1_J:old.ar.mul_xor0(b,B[j+1],A[j],delta)\n'
              '        b.cx(delta,top)\n')


def core_pair_source(a3):
    """The LOADED a3 module's v5 core with the two loops replaced (each must occur exactly once)."""
    src = a3.CORE_SRC_V5
    assert src.count(_FWD_OLD) == 1 and src.count(_REV_OLD) == 1, 'a3 v5 core text changed'
    assert src.count(_DELTA_OLD) == 1, 'a3 v5 delta loop changed'
    return src.replace(_FWD_OLD, _FWD_NEW).replace(_REV_OLD, _REV_NEW).replace(_DELTA_OLD, _DELTA_NEW)


def make_delta(prog_y1):
    """Forward only: _bez_monic's second reduce E[i] ^= E[-1] O[i] (deferred to here) and the delta loop
    delta ^= O[i] A[2i] (after the swaps) as one program per word O[i]: d1 = E[-1] fixed (garbage cached once),
    d2 = A[2i] per word (its garbage cancels inside the program: rule y1free)."""
    assert all(not ((m >> 12) & 0xFFF) for m in prog_y1.cache_lin), 'y1free program carries scalar-2 garbage'
    assert all(c[0] == 0 for c in prog_y1.cache), 'y1free program caches scalar 2'

    def delta_pair(b, O, E, A, delta):
        tsets = [[[E[i][e]] for e in range(12)] + [[delta[e]] for e in range(12)] for i in range(17)]
        emit_cache_map(b, prog_y1, E[-1], A[0], tsets)
        for i in range(17):
            emit_word_map(b, prog_y1, O[i], E[-1], A[2 * i], tsets[i])
        PAIRSP['delta'] = PAIRSP.get('delta', 0) + 1
    return delta_pair


def _norm(src):
    return '\n'.join(line.rstrip() for line in src.strip().splitlines())


def _make_bez(scope_bez, a3, delta=False):
    """_bez_monic without the reduce next to the pair: the FIRST one forward, the LAST one in reverse (with
    delta=True the forward's second reduce is skipped too: make_delta does it beside the delta loop). Keeps a3's S2
    free-stack sync (its wrapper syncs at the reverse entry). The wrapped function must be the kit's
    final_decoder/conversion.py _bez_monic (source checked)."""
    import inspect
    wrapped = getattr(scope_bez, '_d1_wrapped', False)
    inner = scope_bez
    if wrapped:
        inner = dict(zip(scope_bez.__code__.co_freevars, scope_bez.__closure__))['fn'].cell_contents
    assert _norm(inspect.getsource(inner)) == _norm(_BEZ_SRC), 'unexpected _bez_monic'
    G = inner.__globals__

    def _d1_bez(b, O, E, payloads, reverse, history_scratch12, phase):
        if wrapped and reverse:
            a3._sync(b, 's2')
        ar = G['ar']

        def reduce():
            for i in range(len(O)):
                ar.mul_xor0(b, E[-1], O[i], E[i])
        if reverse:
            reduce()
        G['look_root_cached'](b, O, E[:-1], payloads, reverse, phase)
        if not reverse and not delta:
            reduce()
    return _d1_bez


def install(km, ns):
    import sys
    if os.environ.get('DQI_PUB_PAIR_OFF', '0') == '1':        # control build: no change at all
        return {'off': True}
    a3 = sys.modules['dqi_pub_a3']
    assert a3.core_source()[0] == 'v5', 'pub_pair needs DQI_PUB_A3_CORE=v5'
    prog = load_program()
    cells = a3.converter_cells(ns.decode)
    assert cells, 'no converter core'

    def level(b, O, top, Etop, w, Elow):
        pair_level(b, prog, list(O), top, Etop, list(w), list(Elow))
    done = {}
    mine = set()            # functions created here (converter_cells also lists the wrappers' own 'core' cells)
    swapped = 0
    for cell in cells:
        parent, f = None, cell.cell_contents
        while id(f) not in mine and getattr(f, '__wrapped__', None) is not None:
            parent, f = f, f.__wrapped__
        if id(f) in mine:
            continue
        assert f.__name__ == 'dot_cache_core' and '_D1_FIN' in f.__code__.co_names, 'innermost core is not a3 v5'
        if id(f) not in done:
            scope = f.__globals__
            scope['_D1_PAIR_LEVEL'] = level
            dprog = os.environ.get('DQI_PUB_PAIR_DELTA_PROG')
            scope['_D1_DELTA'] = make_delta(Program(os.path.expanduser(dprog))) if dprog else None
            scope['_D1_BEZ'] = _make_bez(scope['_bez_monic'], a3, delta=bool(dprog))
            ns_exec = dict()
            exec(compile(core_pair_source(a3), __file__, 'exec'), scope, ns_exec)
            core_new = ns_exec['dot_cache_core']

            def make_checked(core_new):
                def dot_cache_core(*args, **kw):
                    out = core_new(*args, **kw)
                    if PAIRSP['pending'] is not None or A14['pending'] is not None:
                        raise RuntimeError('pub_pair: a shifted_product call left unpaired in the converter core')
                    return out
                return dot_cache_core
            new = make_checked(core_new)
            mine.add(id(core_new))
            mine.add(id(new))
            done[id(f)] = (f, new)
        new = done[id(f)][1]
        if parent is None:
            if cell.cell_contents is f:
                cell.cell_contents = new
                swapped += 1
        else:
            pc = dict(zip(parent.__code__.co_freevars, parent.__closure__))['core']
            if pc.cell_contents is f:
                pc.cell_contents = new
                swapped += 1
    LOG.update(program=os.path.basename(prog.path), counts=prog.counts(17), cores=len(done), cells=len(cells),
               swapped=swapped)
    if os.environ.get('DQI_PUB_PAIR_SHIFT', '0') == '1':
        LOG['shift_patched'] = install_shift(prog, ns.decode)
    if os.environ.get('DQI_PUB_PAIR_PLAIN', '0') in ('1', 'fwd'):
        # 'fwd': only the forward ops lists (runops); a12's reverse plain is left to B5's plainpair (its pivot word)
        LOG['plain_patched'] = install_plain(prog, ns.decode,
                                             reverse=os.environ.get('DQI_PUB_PAIR_PLAIN') == '1')
    if os.environ.get('DQI_PUB_PAIR_A14', '0') == '1':
        LOG['a14_patched'] = install_a14(ns.decode)
    return dict(LOG)


# ---- the plain levels (n >= 3): forward ops lists run by runops, reverse a12's plain (switch DQI_PUB_PAIR_PLAIN=1) --
_RUNOPS_SRC = '''def runops(b,ops,a,ai,undo=False):
    for op in (reversed(ops)if undo else ops):
        if op[0]=='s':scale(b,op[1],a,ai,op[2]^undo)
        elif op[0]=='m':ar.mul_xor0(b,op[1],op[2],op[3])
        elif op[0]=='c':convolution_xor(b,op[1],op[2],op[3],undo)
        else:raise ValueError(op[0])
'''
_RUNOPS_NEW = '''def runops(b,ops,a,ai,undo=False):
    ops=_D1_PAIR_OPS(ops)
    for op in (reversed(ops)if undo else ops):
        if op[0]=='s':scale(b,op[1],a,ai,op[2]^undo)
        elif op[0]=='m':ar.mul_xor0(b,op[1],op[2],op[3])
        elif op[0]=='c':convolution_xor(b,op[1],op[2],op[3],undo)
        elif op[0]=='P':_D1_PAIR_RUN(b,op)
        else:raise ValueError(op[0])
'''
_PLAIN_LOOP = ('        for v in payloads:\n'
               '            for i in range(len(g)-1):ar.mul_xor0(b,g[i],v[-1],v[i])\n')
_PLAIN_NEW = '        _D1_PAIR_LOOP(b,g,payloads)\n'
PLAIN = dict(prog=None, groups=0, words=0)


def _key(word):
    return tuple(str(q) for q in word)


def _pair_ops(ops):
    """Group ('m', W_k, h1, T1_k) k < n followed by ('m', W_k, h2, T2_k) k < n (n >= 3, same words, inputs and
    targets disjoint) into one ('P', W, h1, h2, T1, T2) op; everything else unchanged."""
    out = []
    i = 0
    L = len(ops)
    while i < L:
        op = ops[i]
        if op[0] == 'm':
            h1 = _key(op[2])
            j = i
            while j < L and ops[j][0] == 'm' and _key(ops[j][2]) == h1:
                j += 1
            n = j - i
            if n >= 3 and j + n <= L and all(ops[j + k][0] == 'm' for k in range(n)):
                h2 = _key(ops[j][2])
                same_w = all(_key(ops[j + k][1]) == _key(ops[i + k][1]) and _key(ops[j + k][2]) == h2
                             for k in range(n))
                ins = set()
                for k in range(n):
                    ins.update(_key(ops[i + k][1]))
                ins.update(h1)
                ins.update(h2)
                outs = set()
                for k in range(n):
                    outs.update(_key(ops[i + k][3]))
                    outs.update(_key(ops[j + k][3]))
                if same_w and h1 != h2 and not (ins & outs) and not (set(h1) & set(h2)):
                    out.append(('P', [ops[i + k][1] for k in range(n)], ops[i][2], ops[j][2],
                                [ops[i + k][3] for k in range(n)], [ops[j + k][3] for k in range(n)]))
                    PLAIN['groups'] += 1
                    PLAIN['words'] += n
                    i = j + n
                    continue
        out.append(op)
        i += 1
    return out


def _pair_run(b, op):
    _, W, h1, h2, T1, T2 = op
    pair_level(b, PLAIN['prog'], list(W), h1, h2, list(T1), list(T2))


def _pair_loop(b, g, payloads):
    """a12 plain's double loop: payload words v[i] ^= g[i] * v[-1] for i < len(g) - 1, both payloads."""
    n = len(g) - 1
    if len(payloads) == 2 and n >= 3:
        v1, v2 = payloads
        pair_level(b, PLAIN['prog'], list(g[:n]), v1[-1], v2[-1], list(v1[:n]), list(v2[:n]))
        PLAIN['groups'] += 1
        PLAIN['words'] += n
        return
    G = _pair_loop.ar
    for v in payloads:
        for i in range(n):
            G.mul_xor0(b, g[i], v[-1], v[i])


def _swap_code(f, src_old_check, new_src, freevars=()):
    """Replace f.__code__ by new_src compiled with the same free variables (checked)."""
    import inspect
    import textwrap
    raw = inspect.getsource(f)
    assert src_old_check(raw), 'unexpected source for ' + f.__name__
    body = textwrap.dedent(new_src(raw))
    if freevars:
        outer = 'def _d1_outer():\n    ' + ' = '.join(freevars) + ' = None\n' + \
                textwrap.indent(body, '    ') + '    return ' + f.__name__ + '\n'
        ns_ = {}
        exec(compile(outer, f.__code__.co_filename, 'exec'), {}, ns_)
        code = [c for c in ns_['_d1_outer'].__code__.co_consts if hasattr(c, 'co_name') and c.co_name == f.__name__][0]
    else:
        ns_ = {}
        exec(compile(body, f.__code__.co_filename, 'exec'), {}, ns_)
        code = ns_[f.__name__].__code__
    assert code.co_freevars == f.__code__.co_freevars, (code.co_freevars, f.__code__.co_freevars)
    f.__code__ = code


# ---- a14's restored reverse frames: one flag-selected coefficient for both payloads (DQI_PUB_PAIR_A14=1) ----------
# a14 (research/agents/a14/reverse_fusion.py) runs, per payload and per word i >= cutoff,
#   with dk.selected(b, g[i], g[i-1], F) as c: mul_xor0(b, c, x, out[i])
# The selected word is the same for both payloads: select once per word (12 ANDs instead of 24). Allocation
# parity: the first payload's selection cycles (12 qubits + 12 bits per word) are replayed without gates first.
A14 = dict(pending=None, frames_paired=0, words=0)
_A14_STUB = ('def shifted_product(b,g,x,out,F):\n'
             '    return _D1_A14_DISPATCH(b,g,x,out,F,frames,frame_to_top,original_product,scope)\n')


def _a14_dispatch(b, g, x, out, F, frames, frame_to_top, original_product, scope):
    pend = A14['pending']
    if pend is None:
        A14['pending'] = (b, list(g), x, out, F)
        return
    b0, g0, x0, out0, F0 = pend
    A14['pending'] = None
    same = (b0 is b and len(g0) == len(g) and all(str(p) == str(q) for gp, gq in zip(g0, g) for p, q in zip(gp, gq))
            and str(F0) == str(F))
    if not same:
        raise RuntimeError('pub_pair: unpaired a14 shifted_product calls')
    restored = False
    if frames:
        frame = frames[-1]
        m = frame['m']
        restored = frame['reverse'] and frame['restored'] and m in frame_to_top
    if not restored:
        original_product(b, g, x0, out0, F)
        original_product(b, g, x, out, F)
        return
    cutoff = frame_to_top[m] - 3
    assert len(out) == m - 1 and len(out0) == m - 1
    ar = A14['ar']
    idx = list(range(cutoff, len(out)))
    for _ in idx:                                    # the first payload's selection cycles, without gates
        z0 = b.alloc_qubits(12)
        m0 = b.alloc_bits(12)
        b.free(z0)
        b.free(m0)
    for i in idx:
        with scope['dk'].selected(b, g[i], g[i - 1], F) as coefficient:
            ar.mul_xor0(b, coefficient, x0, out0[i])
            ar.mul_xor0(b, coefficient, x, out[i])
        A14['words'] += 1
    A14['frames_paired'] += 1


def install_a14(decode):
    import sys
    walk = sys.modules['walk_bound']
    n = 0
    for f in walk.functions(decode):
        if getattr(f, '__name__', '') != 'shifted_product' or not f.__code__.co_filename.endswith('reverse_fusion.py'):
            continue
        if f.__closure__ is None or '_D1_A14_DISPATCH' in f.__code__.co_names:
            continue
        need = ('frames', 'frame_to_top', 'original_product', 'scope')
        assert set(f.__code__.co_freevars) == set(need), f.__code__.co_freevars
        f.__globals__['_D1_A14_DISPATCH'] = _a14_dispatch
        A14['ar'] = f.__globals__['ar']
        _swap_code(f, lambda s: 'scope[\'dk\'].selected' in s and 'mul_xor0' in s, lambda s: _A14_STUB,
                   f.__code__.co_freevars)
        n += 1
    import atexit
    atexit.register(lambda: print('dqi_pub_pair a14: frames', A14['frames_paired'], 'words', A14['words'],
                                  'pending', A14['pending'] is not None, flush=True))
    return n


def install_plain(prog, decode, reverse=True):
    import sys
    PLAIN['prog'] = prog
    walk = sys.modules['walk_bound']
    fns = list(walk.functions(decode))
    n_runops = n_plain = 0
    seen = set()
    for f in fns:
        if id(f) in seen:
            continue
        seen.add(id(f))
        name = getattr(f, '__name__', '')
        if name == 'runops' and not f.__closure__:
            try:
                ok = _norm(__import__('inspect').getsource(f)) == _norm(_RUNOPS_SRC)
            except (OSError, TypeError):
                ok = False
            if not ok:
                continue
            f.__globals__['_D1_PAIR_OPS'] = _pair_ops
            f.__globals__['_D1_PAIR_RUN'] = _pair_run
            _swap_code(f, lambda s: _norm(s) == _norm(_RUNOPS_SRC), lambda s: _RUNOPS_NEW)
            n_runops += 1
        elif reverse and name == 'plain' and f.__code__.co_filename.endswith('known_zero_inverse.py'):
            f.__globals__['_D1_PAIR_LOOP'] = _pair_loop
            _pair_loop.ar = f.__globals__['ar']
            _swap_code(f, lambda s: s.count(_PLAIN_LOOP) == 1, lambda s: s.replace(_PLAIN_LOOP, _PLAIN_NEW),
                       f.__code__.co_freevars)
            n_plain += 1
    import atexit
    atexit.register(lambda: print('dqi_pub_pair plain: groups', PLAIN['groups'], 'words', PLAIN['words'], flush=True))
    return dict(runops=n_runops, plain=n_plain)
