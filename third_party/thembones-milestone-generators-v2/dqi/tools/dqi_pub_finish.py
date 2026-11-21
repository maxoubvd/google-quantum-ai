"""dqi_pub_finish.py -- builder D1: the projective finish as a kernel of our converter (transplant a3).

Mechanism credit: the projective finish of the public DQI decoder by Syed Asad Ali (Apache-2.0,
kaggle_public/syed3000_1113/src_dqi, module projective_finish.py), re-implemented here on our own kit's kernels.
No code of that source is imported or copied.

What it computes (all words GF(2^12); the converter core works in its square-root domain, so do we):
  B'(z) = sum_{k<35} B'_k z^k + z^35 (monic), X = (X_0..X_34) = Laurent coefficients, S = sum_j X_j z^-(j+1),
  H(P)X = polynomial part of P * S.  Our assembly computes C' = H(B')X (35 words).  The finish computes instead
      K' = B' H(M')X + M' H(B')X,   deg K' <= 34,   K' = M' C' mod B',
  with M' = prod_{i<17} q_i, q_i = z^2 + z + kappa_i the minimal polynomial of alpha_i = y + t_i (y^2 + y = KX,
  the kit's GF(2^24) generator; t_i = 2i).  At every root b of B': K'(b) = M'(b) C'(b); the scan folds 1/M
  into its classical mask constant (node_polynomial P -> P * M, M = coefficient-wise square of M').
  In residues: K' mod q_i = B'(alpha_i) * Q(alpha_i) with Q = H(M')X free in X, and K'_34 = sum_j M'_j X_j.
Circuit (in place on X; B' read-only, viewed and restored):
  forward  = free map L (X -> residues of Q, K'_34); per node: view B'(alpha_i) in B'[0:2],
             multiply_forward99 (99 T, 24 scratch, old residue's mask deferred), unview; free CRT^-1.
  reverse  = free CRT; per node (reversed): view, divide_lowwidth369 (<= 369 executed T, 24 scratch),
             CZ(mask, residue), unview; free L^-1.
"""
from functools import lru_cache

NODES = 17
TSTEP = 2


class Field:
    """GF(2^12) of the kit (kickmix GF2Field(12)); multiplication through the given object."""

    def __init__(self, gf):
        self.gf = gf

    def mul(self, a, b):
        return self.gf.mul(a, b)

    def inv(self, a):
        return self.gf.invert(a)


def pw_table(F, KX, t, count):
    """Coordinates (a, b) of alpha^k, alpha = y + t, y^2 + y = KX, for k < count."""
    out = [(1, 0)]
    for _ in range(1, count):
        a, b = out[-1]
        out.append((F.mul(t, a) ^ F.mul(KX, b), a ^ F.mul(t ^ 1, b)))
    return out


def pmul(F, a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] ^= F.mul(x, y)
    return out


@lru_cache(None)
def geometry(gf_id, KX):
    F = _FIELDS[gf_id]
    M = [1]
    pws = []
    for i in range(NODES):
        t = TSTEP * i
        kappa = KX ^ F.mul(t, t) ^ t
        M = pmul(F, M, [kappa, 1, 1])
        pws.append(pw_table(F, KX, t, 36))
    assert len(M) == 35 and M[-1] == 1
    # L: rows 2i+c = component c of Q(alpha_i), Q_k = sum_j M_{k+j+1} X_j (k < 34); row 34 = M (dot)
    L = []
    for i in range(NODES):
        for c in range(2):
            row = []
            for j in range(35):
                v = 0
                for k in range(0, 34 - j):
                    v ^= F.mul(pws[i][k][c], M[k + j + 1])
                row.append(v)
            L.append(row)
    L.append(list(M))
    CRT = [[pws[i][k][c] for k in range(35)] for i in range(NODES) for c in range(2)]
    CRT.append([0] * 34 + [1])
    VIEWS = []
    for i in range(NODES):
        V = [[int(r == k) for k in range(35)] for r in range(35)]
        V[0] = [pws[i][k][0] for k in range(35)]
        V[1] = [pws[i][k][1] for k in range(35)]
        VIEWS.append((tuple(map(tuple, V)), pws[i][35]))
    return dict(M=tuple(M), L=tuple(map(tuple, L)), CRT=tuple(map(tuple, CRT)), VIEWS=tuple(VIEWS))


_FIELDS = {}


def register_field(gf):
    _FIELDS[id(gf)] = Field(gf)
    return id(gf)


def _const_xor(b, word, value):
    for bit in range(12):
        if value >> bit & 1:
            b.x(word[bit])


def _join(km, lo, hi):
    return km.array([*lo, *hi])


class Finish:
    """Bind the kit's primitives: linear_words(b, words, rows, reverse), multiply_forward99, divide_lowwidth369."""

    def __init__(self, km, gf, KX, linear_words, multiply_forward99, divide, mul_into=None):
        self.km, self.KX = km, KX
        self.gid = register_field(gf)
        self.geo = geometry(self.gid, KX)
        self.linear_words = linear_words
        self.mul99 = multiply_forward99
        self.divide = divide
        self.mul_into = mul_into          # z ^= x * y (the kit's mul24_xor99), for the held node

    def scan_factor(self, gf):
        """M = coefficient-wise square of M' (normal domain), ascending."""
        return [gf.mul(m, m) for m in self.geo['M']]

    def _view(self, b, Bw, i, undo=False):
        V, c35 = self.geo['VIEWS'][i]
        if not undo:
            self.linear_words(b, Bw, V)
        _const_xor(b, Bw[0], c35[0])
        _const_xor(b, Bw[1], c35[1])
        if undo:
            self.linear_words(b, Bw, V, True)

    def forward(self, b, Bw, Xw, hold=False):
        """X (35 words) -> K' in place; Bw (35 words, monic implicit) restored. Returns the 17 masks (node 0's is
        None with hold=True: its old residue is kept in a fresh 24-qubit word, returned as `held`)."""
        assert len(Bw) == 35 and len(Xw) == 35
        g = self.geo
        self.linear_words(b, Xw, g['L'])
        masks = [None] * NODES
        held = None
        # with hold, node 0 goes last: its kept word is then allocated when every other scratch is free (width)
        for i in (list(range(1, NODES)) + [0]) if hold else range(NODES):
            self._view(b, Bw, i)
            x = _join(self.km, Xw[2 * i], Xw[2 * i + 1])
            y = _join(self.km, Bw[0], Bw[1])
            if hold and i == 0:
                held = b.alloc_qubits(24)
                self.mul_into(b, x, y, held, self.KX)
                b.swap(x, held)
            else:
                masks[i] = self.mul99(b, x, y, self.KX)
            self._view(b, Bw, i, True)
        self.linear_words(b, Xw, g['CRT'], True)
        return (masks, held) if hold else masks

    def reverse_hold(self, b, Bw, Kw, masks, held):
        """Reverse of forward(hold=True): node 0 first (un-multiply from the held word, then free it), then the
        divisions of nodes 16..1."""
        g = self.geo
        self.linear_words(b, Kw, g['CRT'])
        self._view(b, Bw, 0)
        x = _join(self.km, Kw[0], Kw[1])
        y = _join(self.km, Bw[0], Bw[1])
        self.mul_into(b, held, y, x, self.KX)
        b.swap(x, held)
        b.free(held)
        self._view(b, Bw, 0, True)
        for i in reversed(range(1, NODES)):
            self._view(b, Bw, i)
            x = _join(self.km, Kw[2 * i], Kw[2 * i + 1])
            y = _join(self.km, Bw[0], Bw[1])
            self.divide(b, x, y, self.KX)
            b.cz(masks[i], x)
            b.free(masks[i])
            self._view(b, Bw, i, True)
        self.linear_words(b, Kw, g['L'], True)

    def reverse(self, b, Bw, Kw, masks):
        """K' -> X in place, paying the forward's masks."""
        g = self.geo
        self.linear_words(b, Kw, g['CRT'])
        for i in reversed(range(NODES)):
            V, c35 = g['VIEWS'][i]
            self.linear_words(b, Bw, V)
            _const_xor(b, Bw[0], c35[0])
            _const_xor(b, Bw[1], c35[1])
            x = _join(self.km, Kw[2 * i], Kw[2 * i + 1])
            y = _join(self.km, Bw[0], Bw[1])
            self.divide(b, x, y, self.KX)
            b.cz(masks[i], x)
            b.free(masks[i])
            _const_xor(b, Bw[1], c35[1])
            _const_xor(b, Bw[0], c35[0])
            self.linear_words(b, Bw, V, True)
        self.linear_words(b, Kw, g['L'], True)


def model_K(F, M, Bp, X):
    """Classical K' = B' H(M')X + M' H(B')X for one shot (lists of ints); B' monic of degree 35."""
    Bfull = list(Bp) + [1]

    def H(P, n_out):
        return [_dot(F, [P[k + j + 1] if k + j + 1 < len(P) else 0 for j in range(35)], X) for k in range(n_out)]

    Q = H(list(M), 34)
    C = H(Bfull, 35)
    t1 = pmul(F, Bfull, Q)
    t2 = pmul(F, list(M), C)
    out = [(t1[k] if k < len(t1) else 0) ^ (t2[k] if k < len(t2) else 0) for k in range(max(len(t1), len(t2)))]
    assert all(v == 0 for v in out[35:]), 'degree above 34'
    return out[:35], C


def _dot(F, a, b):
    v = 0
    for x, y in zip(a, b):
        v ^= F.mul(x, y)
    return v
