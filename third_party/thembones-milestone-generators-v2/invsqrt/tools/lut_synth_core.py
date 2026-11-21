"""lut_synth core: GF(2) algebra on 256-entry truth tables (Python ints), a wire machine with AND / CX steps,
the linear solver for one finishing step (one-sided "A" products and two-sided "B" products over a closed
algebra), and a numba Gray-code search over one fixed factor.

Slots: 0 = constant 1, 1 + k = address bit k (k = 0..7, bit 7 = parity), 9.. = wires.
A form is a bit mask over slots. Program steps: ["and", umask, vmask, target_slot, tag], ["cx", srcmask, target_slot].
"""
import numpy as np
from numba import njit

N = 256
FULL = (1 << N) - 1
M64 = (1 << 64) - 1


def var_tt(k):
    v = 0
    for a in range(N):
        if (a >> k) & 1:
            v |= 1 << a
    return v


BASE = [FULL] + [var_tt(k) for k in range(8)]
NB = 9


class Span:
    """Echelon basis keyed by leading bit, with a combination mask carried along."""
    __slots__ = ("piv",)

    def __init__(self):
        self.piv = {}

    def reduce(self, v, c=0):
        piv = self.piv
        while v:
            p = piv.get(v.bit_length() - 1)
            if p is None:
                return v, c
            v ^= p[0]
            c ^= p[1]
        return 0, c

    def add(self, v, c=0):
        v, c = self.reduce(v, c)
        if v == 0:
            return False, c
        self.piv[v.bit_length() - 1] = (v, c)
        return True, c

    def rank(self):
        return len(self.piv)

    def contains(self, v):
        return self.reduce(v)[0] == 0


class Mach:
    def __init__(self, nw=39):
        self.nw = nw
        self.ns = NB + nw
        self.val = list(BASE) + [0] * nw
        self.prog = []
        self.ands = 0

    def form(self, mask):
        v = 0
        s = 0
        val = self.val
        while mask:
            if mask & 1:
                v ^= val[s]
            mask >>= 1
            s += 1
        return v

    def AND(self, um, vm, tgt, tag=""):
        assert tgt >= NB and not ((um | vm) >> tgt) & 1, "target wire inside a control form"
        p = self.form(um) & self.form(vm)
        self.val[tgt] ^= p
        self.prog.append(["and", um, vm, tgt, tag])
        self.ands += 1
        return p

    def CX(self, sm, tgt):
        assert tgt >= NB and not (sm >> tgt) & 1
        self.val[tgt] ^= self.form(sm)
        self.prog.append(["cx", sm, tgt])

    def express(self, f, slots):
        sp = Span()
        for s in slots:
            sp.add(self.val[s], 1 << s)
        r, c = sp.reduce(f)
        return c if r == 0 else None


def replay(prog, nw=39):
    M = Mach(nw)
    for st in prog:
        if st[0] == "and":
            M.AND(st[1], st[2], st[3], st[4] if len(st) > 4 else "")
        else:
            M.CX(st[1], st[2])
    return M


def t_rem(vals, T):
    """Targets that are independent modulo span(vals) (greedy); len = number of unfinished target dimensions."""
    sp = Span()
    for v in vals:
        sp.add(v)
    rem = []
    for t in T:
        if sp.add(t)[0]:
            rem.append(t)
    return rem


def solve_cols(pool, cols, J, T, allow_zero=False):
    """pool, cols: lists of (vector, id mask). Find tau != 0 and a combination with
    sum_{j in tau} T[j] = J + sum(combination). Returns (tau mask, combination mask) or None.
    allow_zero: tau = 0 is acceptable (used when J itself is a target)."""
    sp = Span()
    for v, i in pool:
        sp.add(v, i)
    for v, i in cols:
        sp.add(v, i)
    rows = sorted(sp.piv.items(), reverse=True)

    def full(v, c=0):
        # canonical residual modulo span(pool, cols): every pivot bit cleared (Span.reduce stops too early for that)
        for h, (pv, pc) in rows:
            if (v >> h) & 1:
                v ^= pv
                c ^= pc
        return v, c

    rj, cj = full(J)
    nt = len(T)
    tm = (1 << nt) - 1
    s2 = Span()
    dep = None
    for j, t in enumerate(T):
        rt, ct = full(t)
        ok, comb = s2.add(rt, (1 << j) | (ct << nt))
        if not ok and dep is None:
            dep = comb
    if rj:
        r, comb = s2.reduce(rj, cj << nt)
        if r or not ((comb & tm) or allow_zero):
            return None
        return comb & tm, comb >> nt
    if allow_zero:
        return 0, cj
    if dep is None:
        return None
    return dep & tm, (dep >> nt) ^ cj


class Proj:
    """Quotient of GF(2)^256 by span(vecs): coordinates = the non-pivot bits after full reduction."""

    def __init__(self, vecs):
        sp = Span()
        for v in vecs:
            sp.add(v)
        self.rows = sorted(((h, p[0]) for h, p in sp.piv.items()), reverse=True)
        pivs = set(sp.piv)
        self.nonpiv = [b for b in range(N) if b not in pivs]
        self.c = len(self.nonpiv)
        self.wd = max(1, (self.c + 63) // 64)

    def red(self, v):
        for h, row in self.rows:
            if (v >> h) & 1:
                v ^= row
        return v

    def pack(self, vs):
        out = np.zeros((len(vs), self.wd), dtype=np.uint64)
        nonpiv = self.nonpiv
        for i, v in enumerate(vs):
            r = self.red(v)
            x = 0
            if r:
                for j, b in enumerate(nonpiv):
                    if (r >> b) & 1:
                        x |= 1 << j
            for w in range(self.wd):
                out[i, w] = (x >> (64 * w)) & M64
        return out


@njit(cache=True)
def _search(G, B, g0, b0, start, count):
    """Gray-code walk over x in GF(2)^K: generators g_m(x) = g0_m + sum_k x_k G[k,m], rhs(x) = b0 + sum_k x_k B[k].
    Returns the first Gray index >= start (fewer than start + count) with rhs(x) in span(g(x)), else -1."""
    K, M, W = G.shape
    g = g0.copy()
    r = b0.copy()
    gray = start ^ (start >> 1)
    for k in range(K):
        if (gray >> k) & 1:
            for m in range(M):
                for w in range(W):
                    g[m, w] ^= G[k, m, w]
            for w in range(W):
                r[w] ^= B[k, w]
    wk = np.empty((M, W), dtype=np.uint64)
    t = np.empty(W, dtype=np.uint64)
    one = np.uint64(1)
    for it in range(count):
        for m in range(M):
            for w in range(W):
                wk[m, w] = g[m, w]
        for w in range(W):
            t[w] = r[w]
        for i in range(M):
            pw = -1
            for w in range(W):
                if wk[i, w] != 0:
                    pw = w
                    break
            if pw < 0:
                continue
            a = wk[i, pw]
            piv = a & (~a + one)
            if t[pw] & piv:
                for w in range(pw, W):
                    t[w] ^= wk[i, w]
            for j in range(i + 1, M):
                if wk[j, pw] & piv:
                    for w in range(pw, W):
                        wk[j, w] ^= wk[i, w]
        z = True
        for w in range(W):
            if t[w] != 0:
                z = False
        if z:
            return start + it
        nxt = start + it + 1
        k = 0
        while ((nxt >> k) & 1) == 0:
            k += 1
        if k >= K:
            return -1
        for m in range(M):
            for w in range(W):
                g[m, w] ^= G[k, m, w]
        for w in range(W):
            r[w] ^= B[k, w]
    return -1


@njit(cache=True)
def _search2(D, Q, bbase, y_start, n_outer, x_start, x_count):
    """Two-level Gray search for a two-sided product (u0 + u')(v0 + v') whose fixed parts are both enumerated:
    v0 = sum_j y_j p_j (outer), u0 = sum_k x_k p_k (inner). D[k,c] = proj(delta_c p_k), Q[k,j] = proj(p_k p_j),
    bbase = proj(J). Condition: bbase + proj(u0 v0) in span(delta_c v0) + span(delta_c u0).
    Returns (outer Gray index, inner Gray index) or (-1, -1). x_start applies to the first outer index only."""
    K, H, W = D.shape
    one = np.uint64(1)
    f = np.zeros((H, W), dtype=np.uint64)       # delta_c v0
    Bk = np.zeros((K, W), dtype=np.uint64)      # proj(p_k v0)
    gy = y_start ^ (y_start >> 1)
    for j in range(K):
        if (gy >> j) & 1:
            for c in range(H):
                for w in range(W):
                    f[c, w] ^= D[j, c, w]
            for k in range(K):
                for w in range(W):
                    Bk[k, w] ^= Q[k, j, w]
    E = np.empty((H, W), dtype=np.uint64)
    pw_ = np.empty(H, dtype=np.int64)
    pb_ = np.empty(H, dtype=np.uint64)
    Dr = np.empty((K, H, W), dtype=np.uint64)
    Br = np.empty((K, W), dtype=np.uint64)
    br = np.empty(W, dtype=np.uint64)
    g = np.empty((H, W), dtype=np.uint64)
    r = np.empty(W, dtype=np.uint64)
    wk = np.empty((H, W), dtype=np.uint64)
    t = np.empty(W, dtype=np.uint64)
    for oi in range(n_outer):
        # echelon of f
        ne = 0
        for c in range(H):
            for w in range(W):
                t[w] = f[c, w]
            for i in range(ne):
                if t[pw_[i]] & pb_[i]:
                    for w in range(W):
                        t[w] ^= E[i, w]
            p = -1
            for w in range(W):
                if t[w] != 0:
                    p = w
                    break
            if p >= 0:
                for w in range(W):
                    E[ne, w] = t[w]
                pw_[ne] = p
                pb_[ne] = t[p] & (~t[p] + one)
                ne += 1
        # reduce tables
        for k in range(K):
            for c in range(H):
                for w in range(W):
                    t[w] = D[k, c, w]
                for i in range(ne):
                    if t[pw_[i]] & pb_[i]:
                        for w in range(W):
                            t[w] ^= E[i, w]
                for w in range(W):
                    Dr[k, c, w] = t[w]
            for w in range(W):
                t[w] = Bk[k, w]
            for i in range(ne):
                if t[pw_[i]] & pb_[i]:
                    for w in range(W):
                        t[w] ^= E[i, w]
            for w in range(W):
                Br[k, w] = t[w]
        for w in range(W):
            t[w] = bbase[w]
        for i in range(ne):
            if t[pw_[i]] & pb_[i]:
                for w in range(W):
                    t[w] ^= E[i, w]
        for w in range(W):
            br[w] = t[w]
        # inner Gray walk
        xs = x_start if oi == 0 else 1
        for c in range(H):
            for w in range(W):
                g[c, w] = 0
        for w in range(W):
            r[w] = br[w]
        gx = xs ^ (xs >> 1)
        for k in range(K):
            if (gx >> k) & 1:
                for c in range(H):
                    for w in range(W):
                        g[c, w] ^= Dr[k, c, w]
                for w in range(W):
                    r[w] ^= Br[k, w]
        xi = xs
        while xi < x_count:
            for c in range(H):
                for w in range(W):
                    wk[c, w] = g[c, w]
            pre = False
            for w in range(W):
                t[w] = r[w]
                if r[w] != 0:
                    pre = True
            z = False
            if pre:     # a right-hand side that is already zero is a trivial hit (no target involved): skip it
                for i in range(H):
                    p = -1
                    for w in range(W):
                        if wk[i, w] != 0:
                            p = w
                            break
                    if p < 0:
                        continue
                    a = wk[i, p]
                    piv = a & (~a + one)
                    if t[p] & piv:
                        for w in range(p, W):
                            t[w] ^= wk[i, w]
                    for j in range(i + 1, H):
                        if wk[j, p] & piv:
                            for w in range(p, W):
                                wk[j, w] ^= wk[i, w]
                z = True
                for w in range(W):
                    if t[w] != 0:
                        z = False
            if z:
                return y_start + oi, xi
            xi += 1
            k = 0
            while ((xi >> k) & 1) == 0:
                k += 1
            if k >= K:
                break
            for c in range(H):
                for w in range(W):
                    g[c, w] ^= Dr[k, c, w]
            for w in range(W):
                r[w] ^= Br[k, w]
        # next outer
        nxt = y_start + oi + 1
        j = 0
        while ((nxt >> j) & 1) == 0:
            j += 1
        if j >= K:
            return -1, -1
        for c in range(H):
            for w in range(W):
                f[c, w] ^= D[j, c, w]
        for k in range(K):
            for w in range(W):
                Bk[k, w] ^= Q[k, j, w]
    return -1, -1


def gray_search(G, B, g0, b0, limit_log2, verify, chunk_log2=22, stats=None):
    """Runs _search in chunks; verify(gray_bits) -> result or None (keep searching)."""
    K = G.shape[0]
    total = 1 << min(K, limit_log2)
    start = 1
    while start < total:
        cnt = min(1 << chunk_log2, total - start)
        idx = int(_search(G, B, g0, b0, start, cnt))
        if idx < 0:
            if stats is not None:
                stats["tries"] = stats.get("tries", 0) + cnt
            start += cnt
            continue
        if stats is not None:
            stats["tries"] = stats.get("tries", 0) + idx - start + 1
        res = verify(idx ^ (idx >> 1))
        if res is not None:
            return res
        start = idx + 1
    return None


def bits(mask):
    out = []
    s = 0
    while mask:
        if mask & 1:
            out.append(s)
        mask >>= 1
        s += 1
    return out
