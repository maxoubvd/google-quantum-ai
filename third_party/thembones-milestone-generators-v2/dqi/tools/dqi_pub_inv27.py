"""dqi_pub_inv27.py -- builder D1 (holder 4): a one-operand inverse phase (-1)^(m . inv0(z)) [+ zm [z = 0]] for a
CLASSICAL 12-bit mask m and a quantum 12-bit word z (kit basis, modulus 0x1053), on the cubic tower.

Mechanism credit: the public DQI decoder's 27-Toffoli kernel (Syed Asad Ali, Apache-2.0), read only for its idea;
everything here is derived and searched by us (tools/dqi_pub_inv27model.py, dqi_pub_qmap.py, dqi_pub_xag.py).
  m . w = Tr(alpha w) with alpha = the trace-dual of m (free, classical); z' = z / alpha (free classical scaling);
  then the FIXED phase (-1)^Tr(1 / z'). With z' = a + b theta + c theta^2 (a, b, c in GF(16) = F2[beta],
  beta^4 = beta + 1; theta in GF(8), theta^3 = theta + 1):
      B = b c + a^2 + b^2 + c^2 = e2,   C = (b^3 + b c^2 + c^3) + a B = N(z'),   Tr(1 / z') = Tr16(B / C).
  ANDs: joint map (b, c) -> (b c, b^3 + b c^2 + c^3) [program file], a B into C (9), the leaf (-1)^Tr16(B / C)
  = prod_j (-1)^(B_j h_j(C)) with h_j = Tr16(beta^j / C) from 5 nested ANDs on C; zero term [C = 0] 2 ANDs.
  B and C are X-measured and repaid by classically conditioned (free) phases.
  phase(b, z, mask, zero_mask=None) ; zero_mask: None, a classical bit, or a qubit (then 1 CCZ more).
"""
import json
import os

TOOLS = os.path.dirname(os.path.abspath(__file__))
MOD = 0x1053
PROGS = os.path.expanduser(os.environ.get('DQI_PUB_INV27_DIR', '~/qcc/ws/dqi_pub/inv27'))


def mul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 4096:
            a ^= MOD
    return r


def pw(a, e):
    r = 1
    while e:
        if e & 1:
            r = mul(r, a)
        a = mul(a, a)
        e >>= 1
    return r


def _bits(x):
    i = 0
    while x:
        if x & 1:
            yield i
        x >>= 1
        i += 1


def _solve_columns(cols, n):
    """cols[i] = image of unit i (n-bit ints); returns the inverse map as columns (image of unit i under inverse)."""
    rows = [(cols[i], 1 << i) for i in range(n)]
    piv = {}
    for v, c in rows:
        for pb in sorted(piv, reverse=True):
            if (v >> pb) & 1:
                v ^= piv[pb][0]
                c ^= piv[pb][1]
        assert v, 'singular'
        piv[v.bit_length() - 1] = (v, c)
    inv = []
    for e in range(n):
        w, c = 1 << e, 0
        for pb in sorted(piv, reverse=True):
            if (w >> pb) & 1:
                w ^= piv[pb][0]
                c ^= piv[pb][1]
        assert w == 0
        inv.append(c)
    return inv


def _apply(cols, v):
    r = 0
    for i in _bits(v):
        r ^= cols[i]
    return r


class Tables:
    def __init__(self, joint_path=None):
        beta = next(v for v in range(2, 4096) if pw(v, 4) ^ v ^ 1 == 0)
        theta = next(v for v in range(2, 4096) if pw(v, 3) ^ v ^ 1 == 0)
        self.basis = [mul(pw(theta, i), pw(beta, j)) for i in range(3) for j in range(4)]   # bit 4 i + j
        # cubic coordinates: to_cubic maps a polynomial-basis word to its (a | b << 4 | c << 8) bits
        self.from_cubic = list(self.basis)                       # columns: cubic unit -> polynomial word
        self.to_cubic = _solve_columns(self.from_cubic, 12)
        # trace-dual: alpha = sum_k m_k delta_k with Tr(delta_k x^j) = [j = k]
        def tr(a):
            s, v = 0, a
            for _ in range(12):
                s ^= v
                v = mul(v, v)
            return s & 1
        trcols = [sum(tr(mul(1 << i, 1 << j)) << j for j in range(12)) for i in range(12)]   # Tr(x^i x^j)
        self.dual = _solve_columns(trcols, 12)                   # mask bits -> alpha (columns)
        sq16 = [0] * 4                                           # squaring in GF(16), beta basis
        for j in range(4):
            sq16[j] = self._g16mul(1 << j, 1 << j)
        self.sq16 = sq16
        self.sq12 = [mul(1 << i, 1 << i) for i in range(12)]
        self.xmul = [[mul(1 << i, 1 << j) for j in range(12)] for i in range(12)]          # x^i x^j
        P = json.load(open(os.path.join(PROGS, 'k4_n1_L9_bil.json')))['program']
        self.prod16 = [(int(P['l'][i][:4][::-1], 2), int(P['lp'][i][4:][::-1], 2),
                        sum(int(P['W'][c][i]) << c for c in range(4))) for i in range(9)]
        self._check_prod16()
        X = json.load(open(os.path.join(PROGS, 'xag_k5.json')))['program']
        self.leaf = X
        self.joint = self._joint(joint_path)

    @staticmethod
    def _g16mul(a, b):
        r = 0
        for i in range(4):
            if (b >> i) & 1:
                r ^= a << i
        for i in (6, 5, 4):
            if (r >> i) & 1:
                r ^= 0b10011 << (i - 4)
        return r

    def _check_prod16(self):
        for x in range(16):
            for y in range(16):
                acc = 0
                for la, lb, out in self.prod16:
                    if bin(la & x).count('1') & bin(lb & y).count('1') & 1:
                        acc ^= out
                assert acc == self._g16mul(x, y), 'GF(16) program'

    def _joint(self, path):
        """(gates [(formA8, formB8, out8)], linear [8 masks over the 8 bc bits]) for (b, c) -> (b c | K(b, c) << 4)."""
        def target(v):
            bb, cc = v & 15, v >> 4
            K = self._g16mul(self._g16mul(bb, bb), bb) ^ self._g16mul(bb, self._g16mul(cc, cc)) ^ \
                self._g16mul(self._g16mul(cc, cc), cc)
            return self._g16mul(bb, cc) | (K << 4)
        gates = []
        if path:
            P = json.load(open(path))['program']
            L = len(P['l'])
            for i in range(L):
                fa = sum(int(ch) << u for u, ch in enumerate(P['l'][i]))
                fb = sum(int(ch) << u for u, ch in enumerate(P['lp'][i]))
                out = sum(int(P['W'][c][i]) << c for c in range(8))
                gates.append((fa, fb, out))
        else:                       # naive: b c, b * b^2, b * c^2, c * c^2 by the 9-term product (exact as functions)
            sqb = [sum(((self.sq16[j] >> k) & 1) << j for j in range(4)) for k in range(4)]   # rows of Sq
            def sqform(f):          # form on y = Sq(x): f . Sq(x) = (Sq^T f) . x
                r = 0
                for k in range(4):
                    if (f >> k) & 1:
                        r ^= sqb[k]
                return r
            for la, lb, out in self.prod16:
                gates.append((la, lb << 4, out))                    # b * c -> B
                gates.append((la, sqform(lb), out << 4))            # b * b^2 -> C
                gates.append((la, sqform(lb) << 4, out << 4))       # b * c^2 -> C
                gates.append((la << 4, sqform(lb) << 4, out << 4))  # c * c^2 -> C
        # canonical gates: (fa, fb) unordered, outputs merged; fa = fb or a zero form is linear (left to the correction)
        merged = {}
        for fa, fb, out in gates:
            if fa == 0 or fb == 0 or fa == fb:
                continue
            key = (min(fa, fb), max(fa, fb))
            merged[key] = merged.get(key, 0) ^ out
        gates = [(fa, fb, out) for (fa, fb), out in sorted(merged.items()) if out]
        # linear correction: target minus the gates' values (difference must be linear: check on all 256 inputs)
        def prog(v):
            acc = 0
            for fa, fb, out in gates:
                if bin(fa & v).count('1') & bin(fb & v).count('1') & 1:
                    acc ^= out
            return acc
        diff = [target(v) ^ prog(v) for v in range(256)]
        assert diff[0] == 0
        lin = [0] * 8
        for u in range(8):
            d = diff[1 << u]
            for c in range(8):
                if (d >> c) & 1:
                    lin[c] |= 1 << u
        for v in range(256):
            acc = 0
            for c in range(8):
                if bin(lin[c] & v).count('1') & 1:
                    acc |= 1 << c
            assert acc == diff[v], 'joint program is not quadratic-exact'
        return gates, lin


T = None


def tables():
    global T
    if T is None:
        jp = os.environ.get('DQI_PUB_INV27_JOINT')
        T = Tables(os.path.expanduser(jp) if jp else None)
    return T


# ---- small emission helpers --------------------------------------------------------------------------------------

def _coords(content, f):
    n = len(content)
    piv = {}
    for k in range(n):
        v, c = content[k], 1 << k
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


class _View:
    """In-place CNOT conjugation of a list of qubits to expose parities (undo restores)."""

    def __init__(self, b, qs):
        self.b, self.qs = b, qs
        self.content = [1 << k for k in range(len(qs))]
        self.done = []

    def cnot(self, c, t):
        self.b.cx(self.qs[c], self.qs[t])
        self.content[t] ^= self.content[c]
        self.done.append((c, t))

    def expose(self, f, avoid=()):
        S = _coords(self.content, f)
        cand = [k for k in S if k not in avoid]
        assert cand, 'cannot expose'
        p = cand[0]
        for k in S:
            if k != p:
                self.cnot(k, p)
        return p

    def undo(self):
        for c, t in reversed(self.done):
            self.b.cx(self.qs[c], self.qs[t])
        self.done = []
        self.content = [1 << k for k in range(len(self.qs))]


def _and_into(b, qs, fa, fb, targets, kind='ccx'):
    """XOR (parity(qs & fa) AND parity(qs & fb)) into every qubit of targets (one CCX), or init_and into one."""
    V = _View(b, qs)
    pa = V.expose(fa)
    pb = V.expose(fb, avoid=(pa,))
    if kind == 'init':
        b.init_and(qs[pa], qs[pb], targets[0])
    elif kind == 'del':
        b.del_and(qs[pa], qs[pb], targets[0])
    elif kind == 'cz':                         # phase gate CZ(pa, pb), conditioned by the caller
        b.cz(qs[pa], qs[pb])
    else:
        t0 = targets[0]
        for t in targets[1:]:
            b.cx(t0, t)
        b.ccx(qs[pa], qs[pb], t0)
        for t in reversed(targets[1:]):
            b.cx(t0, t)
    V.undo()


def _parity_bit(b, bits_list, mask):
    """A fresh classical bit = parity of the selected classical bits."""
    p = b.alloc_bits(1)[0]
    b.bit_store(p, False)
    for i in _bits(mask):
        b.cx(bits_list[i], p)
    return p


def _classical_linear(b, src, cols, n_out):
    """Fresh classical word = linear map of a classical word (columns: image of unit i)."""
    out = b.alloc_bits(n_out)
    b.bit_store(out, False)
    for i, col in enumerate(cols):
        for k in _bits(col):
            b.cx(src[i], out[k])
    return out


def _classical_mul(b, u, v, t):
    """Fresh classical word = u * v (GF(4096)), classical bits only."""
    out = b.alloc_bits(12)
    b.bit_store(out, False)
    for i in range(12):
        b.push_condition(u[i])
        for j in range(12):
            for k in _bits(t.xmul[i][j]):
                b.cx(v[j], out[k])
        b.pop_condition()
    return out


def _classical_inverse(b, a, t):
    """Fresh classical word a^-1 = a^4094 (a = 0 -> 0) and the list of words to free."""
    owned = []
    cur = None
    base = a
    e = 4094
    sq = base
    first = True
    while e:
        if e & 1:
            if first:
                cur = _classical_linear(b, sq, [1 << i for i in range(12)], 12)
                first = False
            else:
                nxt = _classical_mul(b, cur, sq, t)
                owned.append(cur)
                cur = nxt
        e >>= 1
        if e:
            nsq = _classical_linear(b, sq, t.sq12, 12)
            if sq is not base:
                owned.append(sq)
            sq = nsq
    if sq is not base:
        owned.append(sq)
    return cur, owned


def _scale(b, z, m, minv, t):
    """In place z <- z * m for classical m (minv = m^-1): free, 12 temporary qubits."""
    tmp = b.alloc_qubits(12)
    for j in range(12):
        b.push_condition(m[j])
        for i in range(12):
            for k in _bits(t.xmul[i][j]):
                b.cx(z[i], tmp[k])
        b.pop_condition()
    for i in range(12):
        b.swap(z[i], tmp[i])
    for j in range(12):
        b.push_condition(minv[j])
        for i in range(12):
            for k in _bits(t.xmul[i][j]):
                b.cx(z[i], tmp[k])
        b.pop_condition()
    b.free(tmp)


def _linear_inplace(b, x, cols, km_linear):
    km_linear(b, x, cols)


def phase(b, z, mask, zero_mask=None, linear_inplace=None):
    """(-1)^(mask . inv0(z)) * (-1)^(zero_mask [z = 0]); z quantum (12), mask classical (12)."""
    t = tables()
    assert len(z) == len(mask) == 12
    alpha = _classical_linear(b, mask, t.dual, 12)
    ainv, owned = _classical_inverse(b, alpha, t)
    nz = b.alloc_bits(1)[0]                                    # mask != 0
    b.bit_store(nz, False)
    for bit in mask:
        b.bit_store(nz, True, control=bit)
    b.push_condition(nz)
    _scale(b, z, ainv, alpha, t)
    b.pop_condition()
    core(b, z, nz, zero_mask, linear_inplace)
    b.push_condition(nz)
    _scale(b, z, alpha, ainv, t)
    b.pop_condition()
    b.free(nz)
    for w in owned:
        b.free(w)
    b.free(ainv)
    b.free(alpha)


def core(b, z, nz, zero_mask=None, linear_inplace=None):
    """(-1)^(nz Tr(1 / z)) (-1)^(zero_mask [z = 0]) for the quantum word z (nz a classical condition bit)."""
    t = tables()
    lin = linear_inplace
    lin(b, z, t.to_cubic)                                      # z now holds a | b << 4 | c << 8
    a, bw, cw = list(z[:4]), list(z[4:8]), list(z[8:])
    B = b.alloc_qubits(4)
    C = b.alloc_qubits(4)
    BC = list(B) + list(C)
    bc = bw + cw
    gates, jlin = t.joint
    for fa, fb, out in gates:
        _and_into(b, bc, fa, fb, [BC[e] for e in _bits(out)])
    for c in range(8):
        for u in _bits(jlin[c]):
            b.cx(bc[u], BC[c])
    for w in (a, bw, cw):                                      # B += a^2 + b^2 + c^2
        for i in range(4):
            for k in _bits(t.sq16[i]):
                b.cx(w[i], B[k])
    aB = a + list(B)
    for la, lb, out in t.prod16:                               # C += a B
        _and_into(b, aB, la, lb << 4, [C[e] for e in _bits(out)])
    # the leaf: signals = C bits + 5 ANDs; phase (-1)^(B_j h_j) conditioned on nz
    X = t.leaf
    anc = b.alloc_qubits(5)
    sig = list(C) + list(anc)
    for j in range(5):
        fa = sum(v << i for i, v in enumerate(X['A'][j]))
        fb = sum(v << i for i, v in enumerate(X['B'][j]))
        _and_into(b, sig[:4 + j], fa, fb, [anc[j]], kind='init')
    b.push_condition(nz)
    for jj in range(4):
        f = sum(v << i for i, v in enumerate(X['C'][jj]))
        V = _View(b, sig)
        p = V.expose(f)
        b.cz(B[jj], sig[p])
        V.undo()
    b.pop_condition()
    if zero_mask is not None:                                  # (-1)^(zm [C = 0]): [C = 0] = (1+C0)(1+C1)(1+C2)(1+C3)
        zt = b.alloc_qubits(2)
        quantum_zm = isinstance(zero_mask, type(z[0]))
        if not quantum_zm:                                     # a classical zm: the 2 ANDs execute only when zm = 1
            b.push_condition(zero_mask)
        b.x(C)
        b.init_and(C[0], C[1], zt[0])
        b.init_and(C[2], C[3], zt[1])
        if quantum_zm:
            b.ccz(zt[0], zt[1], zero_mask)
        else:
            b.cz(zt[0], zt[1])
        b.del_and(C[2], C[3], zt[1])
        b.del_and(C[0], C[1], zt[0])
        b.x(C)
        if not quantum_zm:
            b.pop_condition()
        b.free(zt)
    for j in reversed(range(5)):
        fa = sum(v << i for i, v in enumerate(X['A'][j]))
        fb = sum(v << i for i, v in enumerate(X['B'][j]))
        _and_into(b, sig[:4 + j], fa, fb, [anc[j]], kind='del')
    b.free(anc)
    # erase C (it holds K(b, c) + a B + joint linear terms): X-measure, repay (-1)^(mC . C) by conditioned phases
    mC = b.alloc_bits(4)
    b.hmr(C, mC)
    b.free(C)
    _repay(b, bc, aB, mC, gates, jlin, t, part='C')
    mB = b.alloc_bits(4)
    b.hmr(B, mB)
    b.free(B)
    _repay(b, bc, None, mB, gates, jlin, t, part='B', words=(a, bw, cw))
    b.free(mC)
    b.free(mB)
    lin(b, z, t.from_cubic)


def _repay(b, bc, aB, m, gates, jlin, t, part, words=None):
    """Phase (-1)^(m . value) of the value the measured register held, by classically conditioned phases."""
    sh = 4 if part == 'C' else 0
    for fa, fb, out in gates:                                  # joint products
        o = (out >> sh) & 15
        if o:
            p = _parity_bit(b, m, o)
            b.push_condition(p)
            _and_into(b, bc, fa, fb, None, kind='cz')
            b.pop_condition()
            b.free(p)
    for c in range(4):                                         # joint linear terms
        for u in _bits(jlin[c + sh]):
            b.push_condition(m[c])
            b.z(bc[u])
            b.pop_condition()
    if part == 'C':                                            # a B products into C
        for la, lb, out in t.prod16:
            p = _parity_bit(b, m, out)
            b.push_condition(p)
            _and_into(b, aB, la, lb << 4, None, kind='cz')
            b.pop_condition()
            b.free(p)
    else:                                                      # squares into B
        for w in words:
            for i in range(4):
                for k in _bits(t.sq16[i]):
                    b.push_condition(m[k])
                    b.z(w[i])
                    b.pop_condition()


# ---- allocation parity with the kit kernel, and the install (patch token pub_inv27) ----------------------------
# The kit's free stacks are LIFO and later routing depends on them. The kernel runs on an arena: every scratch qubit
# comes from one pool of POOL qubits (the kit's 15), bit frees are deferred; at the end everything is freed in reverse
# allocation order with reversed lists (net zero on both stacks, measured), and then the kit kernel's own alloc / free
# sequence (recorded once per zero-mask mode on a scratch builder) is replayed without gates.

POOL = 15


def _items(x):
    try:
        return list(x), False
    except TypeError:
        return [x], True


class _Recorder:
    def __init__(self, b):
        self._b = b
        self.allocs, self.ev, self.where = [], [], {}

    def __getattr__(self, k):
        return getattr(self._b, k)

    def _alloc(self, kind, n):
        r = self._b.alloc_qubits(n) if kind == 'q' else self._b.alloc_bits(n)
        i = len(self.allocs)
        self.allocs.append((kind, n))
        for p, e in enumerate(r):
            self.where[str(e)] = (i, p)
        self.ev.append(('a', i))
        return r

    def alloc_qubits(self, n):
        return self._alloc('q', n)

    def alloc_bits(self, n):
        return self._alloc('b', n)

    def free(self, x):
        items, single = _items(x)
        self.ev.append(('f', [self.where[str(e)] for e in items], single))
        return self._b.free(x)


def _replay(b, allocs, ev):
    got = {}
    for e in ev:
        if e[0] == 'a':
            kind, n = allocs[e[1]]
            got[e[1]] = b.alloc_qubits(n) if kind == 'q' else b.alloc_bits(n)
        else:
            items = [got[i][p] for i, p in e[1]]
            b.free(items[0] if e[2] else items)


class _Arena:
    def __init__(self, b, nq):
        self._b = b
        self._pool_alloc = b.alloc_qubits(nq)
        self._arr = type(self._pool_alloc)
        self._qtype = type(self._pool_alloc[0])
        self._pool = list(self._pool_alloc)[::-1]
        self._bits = []
        self.peak = 0

    def __getattr__(self, k):
        return getattr(self._b, k)

    def alloc_qubits(self, n):
        assert len(self._pool) >= n, 'inv27: scratch pool exhausted'
        out = [self._pool.pop() for _ in range(n)]
        self.peak = max(self.peak, len(self._pool_alloc) - len(self._pool))
        return self._arr(out)

    def alloc_bits(self, n):
        r = self._b.alloc_bits(n)
        self._bits.append(r)
        return r

    def free(self, x):
        items, _single = _items(x)
        if items and isinstance(items[0], self._qtype):
            self._pool.extend(reversed(items))

    def close(self):
        assert len(self._pool) == len(self._pool_alloc), 'inv27: scratch not returned'
        for r in reversed(self._bits):
            self._b.free(list(reversed(list(r))))
        self._b.free(list(reversed(list(self._pool_alloc))))


def _new_builder(b):
    import sys
    for nm in ('kickmix', 'circ'):
        m = sys.modules.get(nm)
        if m is not None and hasattr(m, 'CircuitBuilder'):
            return m.CircuitBuilder()
    return type(b)()


def _mode(x, zero_mask):
    if zero_mask is None:
        return 'none'
    return 'qubit' if isinstance(zero_mask, type(x[0])) else 'bit'


def make_nonzero(orig, lin, counter=None, modes=None):
    """The replacement of the kit's phase_inverse12_virtual_nonzero(b, x, mask, zero_mask=None): same phase, same
    allocation effect as orig (its event sequence recorded per mode on a scratch builder). modes: the zero-mask
    modes this kernel takes ('none', 'bit', 'qubit'; default all); the others run the kit's own kernel."""
    rec = {}
    modes = set(modes or ('none', 'bit', 'qubit'))

    def nonzero(b, x, mask, zero_mask=None):
        mode = _mode(x, zero_mask)
        if mode not in modes:
            if counter is not None:
                counter['kit_' + mode] += 1
            return orig(b, x, mask, zero_mask)
        if mode not in rec:
            sb = _new_builder(b)
            sx = sb.alloc_qubits(12)
            sm = sb.alloc_bits(12)
            szm = None if mode == 'none' else (sb.alloc_qubits(1)[0] if mode == 'qubit' else sb.alloc_bits(1)[0])
            R = _Recorder(sb)
            orig(R, sx, sm, szm)
            rec[mode] = (R.allocs, R.ev)
        A = _Arena(b, POOL)
        phase(A, x, mask, zero_mask, linear_inplace=lin)
        A.close()
        _replay(b, *rec[mode])
        if counter is not None:
            counter[mode] += 1
    return nonzero


LOG = {}


def install(km, ns):
    """Patch token pub_inv27: every kit phase_inverse12_virtual_nonzero (exact_subs.py; module attributes, private
    rebind copies, captured references) gets this kernel through a __code__ swap (stub calls _D1_INV27)."""
    import atexit
    import collections
    import sys
    import types
    if os.environ.get('DQI_PUB_INV27_OFF', '0') == '1':
        return {'off': True}
    name = 'phase_inverse12_virtual_nonzero'
    cands = []
    for m in list(sys.modules.values()):
        d = getattr(m, '__dict__', None)
        if isinstance(d, dict) and isinstance(d.get(name), types.FunctionType):
            cands.append(d[name])
    walk = sys.modules.get('walk_bound')
    if walk is not None and ns is not None and getattr(ns, 'decode', None) is not None:
        cands.extend(f for f in walk.functions(ns.decode) if getattr(f, '__name__', '') == name)
        for f in walk.functions(ns.decode):                    # copies held in scopes (private_rebind)
            g = getattr(f, '__globals__', None)
            if isinstance(g, dict) and isinstance(g.get(name), types.FunctionType):
                cands.append(g[name])
    cands = [f for f in cands if f.__code__.co_filename.endswith('exact_subs.py') and f.__code__.co_name == name
             and not f.__closure__]
    assert cands, 'pub_inv27: no kit phase_inverse12_virtual_nonzero found'
    tables()                                                   # load the programs now (fail at install, not mid-decode)
    f0 = cands[0]
    code0 = f0.__code__
    assert all(f.__code__ is code0 or '_D1_INV27' in f.__code__.co_names for f in cands), 'pub_inv27: code differs'
    orig = types.FunctionType(code0, f0.__globals__, name, f0.__defaults__)
    counter = collections.Counter()
    modes = [m for m in os.environ.get('DQI_PUB_INV27_MODES', 'none,bit,qubit').split(',') if m]
    done, seen = 0, set()
    for f in cands:
        if id(f) in seen or '_D1_INV27' in f.__code__.co_names:
            continue
        seen.add(id(f))
        G = f.__globals__
        G['_D1_INV27'] = make_nonzero(orig, G['linear_inplace'], counter, modes)
        stub = {}
        exec('def %s(b, x, mask, zero_mask=None):\n return _D1_INV27(b, x, mask, zero_mask)\n' % name, G, stub)
        f.__code__ = stub[name].__code__
        done += 1
    jp = os.environ.get('DQI_PUB_INV27_JOINT', '')
    LOG.update(patched=done, joint=os.path.basename(jp) or 'naive', joint_gates=len(tables().joint[0]), modes=modes)
    atexit.register(lambda: print('dqi_pub_inv27 calls by zero-mask mode', dict(counter), LOG, flush=True))
    return dict(LOG)
