"""r48_dqi_ret2_lateflag.py -- seat R48-DQI-RET2: the zero flag of the 15 flagged un-synthesis divisions made AFTER
the 24-bit product, from the 6-bit subnorm (5 ANDs) instead of from the 24-bit divisor before it (23 ANDs).

What the line does (kit `norm_deferred_flag_return.divide_from_norms` + r6n + latem + div6f, called by B5's a1 sweep):
    F = [W = 0] (24-bit test, 23) -> completion: [U, T] += F lead co (12) -> premap(U) -> W made safe (W + F) ->
    U <- U conj(W) (99) -> N, halves, M, 1 / M, parts ... -> body -> F retired through bit 0 of M t.
F is needed before the product only so that an identity division (W + F = 1) keeps the completed dividend.
On every valid input a fold with W = 0 is a TERMINATED fold: its head is a ghost (0), and the raw dividend is
U_raw = lead x (c0, c1), the padding digit of the frame's leader word ((c0, c1) = the first two of the kit's
constant_qadic constants). So the product by conj(0) leaves the right head (0) without any flag, the measured old
dividend is a linear function of a live word (its mask is paid by classically conditioned CZ gates between F and the
leader), and the flag can be read from M = Norm6(Norm(W)) once the kernel has made it: M = 0 <=> W = 0.

Here (flagged divisions only; the 18 unflagged ones are untouched):
    premap_inf(U_raw) -> U <- U conj(W) (99; 0 when W = 0) -> N = Norm(W) (42) -> halves (84) -> M (15) ->
    F = [M = 0] (5) -> N, M made safe -> 1 / M, parts, rebuild (as built) -> masks paid ->
    completion: T' += F lead corr (12; T' = the target words without the leader, corr = co_T + premap_inf(c0, c1)),
    phase (-1)^(F <oldmask, lead (c0, c1)>) (CZ) -> N, M restored -> body_inf -> F retired as built.
premap_inf / body_inf = D1's CRT body with its last point at infinity (tools/r48_dqi_ret2_forms.py): the premap
then never writes the leader word, so the raw dividend of a terminated fold cannot rescale it; same product count.
Per flagged division: -23 + 5 = -18; head 33's completion has no nonzero constant left (-12): -282 per decode [P].

Exactness: identical gates on W != 0 except the body's point set (exact by the floor-map check). On W = 0 it is
exact when the head is 0 and U_raw = lead (c0, c1): the terminated fold of the kit's zero-padding synthesis. A fold
with W = 0 that is NOT terminated cannot occur on a decodable input: K_(i+1)(gamma_i) = 0 gives q_i | K_i, K_(i+1),
hence q_i | K_1, K_2, and the locator has no irreducible quadratic factor (design file, section 5).

Install as a DQI_INT_EXT module AFTER a1 / div6f (token `ext` is behind them in the line's DQI_INT_P) and before
r46_dqi_conv_zpark. R48_LF=0 switches it off (the module then changes nothing).

R48_LATEBODY=1 (second piece, L8; a pure reorder, no Toffoli added or removed): in all 33 divisions the caller's body
`after()` runs AFTER the kernel has retired F and erased t, M and N, instead of before. As built, those 24 to 25
helper wires stay live through every body (594 products, most of un-synthesis's operations: P1 measured 1 to 2 free
wires there); the erasures read only y, N, M and t, none of which the body changes. With the reorder the bodies run
at 26 free wires: holder capacity for P1's links. Default 0 (L1 alone).
R48_KPOINT=1 (L1c): the body of a 3-word target (head 32) takes the finite points {0, k}: both are roots of the padding
polynomial, so its completion has no nonzero constant left (-12).
R48_N3FLAG=1 (L1d): the flag is read from the 3-bit norm inside div6's GF(64) inverse (2 ANDs instead of 5: -45).
R48_BETA0=1 (P2', R37-DQI-DIV's P2 reopened): the sweep's fresh leader beta is the constant 1 until frame 0's body.
With the infinity forms the premap never writes it, so beta is not allocated before frame 0's division (unflagged):
that division enters at 864 live and takes builder B6's dr form (tools/dqi_kern_div, read only: 310 at 36 helpers)
instead of late-M (334.94 at 26), with the body after the kernel's erasures; beta is allocated and set to 1 at the
start of the body and moved onto the line's recorded wires after it: -24.94.
"""
import atexit
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import r48_dqi_ret2_forms as FORMS                               # noqa: E402

ON = os.environ.get('R48_LF', '1') != '0'
LATEBODY = os.environ.get('R48_LATEBODY', '0') == '1'
KPOINT = os.environ.get('R48_KPOINT', '0') == '1'
N3FLAG = os.environ.get('R48_N3FLAG', '0') == '1'
BETA0 = os.environ.get('R48_BETA0', '0') == '1'
DIV6ALL = os.environ.get('R48_DIV6ALL', '0') == '1'   # fifth instance: div6's 12-live inverse in the UNFLAGGED calls too
STATS = dict(lateflag=0, passed=0, pre=0, bodies=0, completions=0, empty_completions=0, sweeps=0, latebody=0,
             n3flag=0, beta0=0)

SWEEP_SRC = '''
def sweep(b, C, B, folds, divide):
    previous.qadic(b, C, B, k=previous.K)
    beta = b.alloc_qubits(12)
    if not _R48_BETA0:
        b.x(beta[0])
    leaders = [C[-2], C[-1], beta]
    H = b.alloc_qubits(24)
    h0 = [H[:12], H[12:]]
    b.free(H)
    beta_home = beta
    if _R48_BETA0:
        b.free(beta)                    # the line's wires for beta, recorded; free through frame 0's division
        beta = leaders[2] = None
    A = h0 + list(C[:-2])
    heads = previous.prior.rr.tail.heads_from_banks(A, B)
    frames, _, _ = previous.frames_for(A, B)
    for j, (a, d, k) in enumerate(frames):
        n = len(a)
        a1, a2, beta = leaders
        U = previous.combine(a[:2])
        W = previous.combine(d[:2])
        ft = list(a[2:]) + [a1, a2]
        fs = list(d) + [beta]
        st = list(d[2:]) + [beta]
        ss = list(a[2:]) + [a1, a2]
        kfirst, ksecond, _ = previous.rank.bodies(b, a, d, leaders, k)
        if j:
            previous.qframe(b, a, d, leaders, frames[j-1][2], True)
            previous.qframe(b, a, d, leaders, k)
            nf = id(folds[2*j][1]) in _B5_NOFLAG
            crt = nf or _B5_MODE != 'noflagonly'
            assert crt and _B5_MODE == 'kitco', _B5_MODE
            if nf:
                def first():
                    _D1A1_BODY(b, a[:2], fs, ft, k)
                _B5_LOG['pre_outside'] += 1
                _D1A1_PRE(b, a[:2], ft, k)
                divide_from_norms(b, U, W, folds[2*j][1], _R48_NF, first, divide, K=k)
            else:
                def first_inf():
                    _R48_BODY(b, a[:2], fs, ft, k)
                _R48_PRE(b, a[:2], ft, k)
                divide_from_norms(b, U, W, folds[2*j][1], _R48_LF(b, a2, ft, k, n + 1, n + 2), first_inf, divide, K=k)
            _B5_LOG['crt_bodies'] += 1
        source = previous.combine(a[2:4]) if n > 3 else previous.combine([a[-1], a1])
        nf2 = id(folds[2*j+1][1]) in _B5_NOFLAG
        if j == 0 and _R48_BETA0:
            assert nf2 and beta is None
            def second_b0():
                bt = b.alloc_qubits(12)
                b.x(bt[0])
                leaders[2] = bt
                _R48_BODY(b, d[:2], ss, list(d[2:]) + [bt], k)
            _R48_PRE(b, d[:2], list(d[2:]) + [None], k)
            _R48_DR(b, W, source, folds[2*j+1][1], second_b0, divide, k)
            leaders[2] = _R48_HOME(b, leaders[2], beta_home)      # back on the line's wires (the tail reclaims h0's)
            _R48_STATS['beta0'] += 1
        elif nf2:
            def second():
                _D1A1_BODY(b, d[:2], ss, st, k)
            _B5_LOG['pre_outside'] += 1
            _D1A1_PRE(b, d[:2], st, k)
            divide_from_norms(b, W, source, folds[2*j+1][1], _R48_NF, second, divide, K=k)
        else:
            def second_inf():
                _R48_BODY(b, d[:2], ss, st, k)
            _R48_PRE(b, d[:2], st, k)
            divide_from_norms(b, W, source, folds[2*j+1][1], _R48_LF(b, beta, st, k, n, n + 1), second_inf, divide, K=k)
        _B5_LOG['crt_bodies'] += 1
    return A, heads, leaders
'''


def _bits(v):
    i = 0
    while v:
        if v & 1:
            yield i
        v >>= 1
        i += 1


class LateFlag:
    """The caller's `before` of a late-flag division: called as before(F, oldmask) AFTER the kernel's arithmetic."""

    def __init__(self, b, lead, T, c0, c1, corr, cl_xor, mul):
        self.b, self.lead, self.T, self.c0, self.c1, self.corr = b, lead, T, c0, c1, corr
        self.cl_xor, self.mul = cl_xor, mul

    def __call__(self, F, oldmask):
        b, lead, mul = self.b, self.lead, self.mul
        assert self.corr[-1] == 0 and len(self.corr) == len(self.T)
        if any(self.corr[:-1]):
            self.cl_xor(b, F, lead, self.T[:-1], self.corr[:-1])           # T' += F lead corr: 12 CCX
            STATS['completions'] += 1
        else:
            STATS['empty_completions'] += 1
        # old dividend on a terminated fold = (c0 lead, c1 lead): phase (-1)^(F <oldmask, that>), CZ only
        mk = b.alloc_bits(12)
        b.bit_store(mk, False)
        for i in range(12):
            for j in _bits(mul(self.c0, 1 << i)):
                b.cx(oldmask[j], mk[i])
            for j in _bits(mul(self.c1, 1 << i)):
                b.cx(oldmask[12 + j], mk[i])
        for i in range(12):
            b.ccz(F, lead[i], mk[i])
        b.free(mk)


class NoFlag:
    """Marker passed as `before` by the sweep for an unflagged division (the unflagged kernel never calls it)."""


NF_MARK = NoFlag()
BODY = "    after()\n"


def _edits(r6n, latem):
    return (
        (r6n.HEAD, "    early = b.alloc_bits(1)\n    b.bit_store(early, False)\n"),
        (r6n.PRODUCT, r6n.LATE_NORM),
        (r6n.SECOND, ""),
        (latem.MBLOCK, ""),
        (latem.LOOP_END, latem.LOOP_END + latem.MBLOCK + ("" if N3FLAG else r6n.SECOND)),
        (r6n.UNPATCH, "    b.cx(F, N[0])\n    b.cx(F, M[0])\n" + ("" if LATEBODY else BODY)),
        ("    b.free(oldmask)\n", "    before(F, oldmask)\n    b.free(oldmask)\n"),
    ) + (((PRODUCE, PRODUCE3),) if N3FLAG else ())


PRODUCE = "        inverse_mask, kept_M = g['produce_retained'](b, M)\n"
PRODUCE3 = "        inverse_mask, kept_M, F = g['produce_retained'](b, M)\n        b.cx(F, N[0])\n"


def _make_n3(G, div6, zero_flag):
    """div6's in-place GF(64) inverse (tools/dqi_junc_r8p_div6.py, read only) with the divisor's zero flag read from
    its own 3-bit norm: Norm3(M) = 0 <=> M = 0 <=> y = 0 (2 ANDs). M and the norm are made safe (the element 1) at
    once; the caller makes N safe from the returned flag. Returns (produce, rebuild); produce -> (mask, kept, F)."""
    ar = G['ar']
    norm = G['norm_xor_noscratch']
    mul = G['mul_xor_noscratch']
    inv3 = sys.modules['exp.noscratch.kernels'].inverse3_noscratch
    TO6, FRM6, K6 = ar.tower(6)
    assert FRM6[0] == 1, 'tower coordinate 0 is not the element 1'
    _, rebuild = div6._make(G)

    def produce(b, M):
        STATS['n3flag'] += 1
        ar.linear_inplace(b, M, TO6)
        N3 = b.alloc_qubits(3)
        norm(b, M, N3, K6)                         # Norm of M over GF(8): 0 exactly when M = 0
        F = zero_flag(b, N3)                       # [y = 0], 2 ANDs
        b.cx(F, M[0])                              # M safe (1 when it was 0)
        b.cx(F, N3[0])                             # its norm (1)
        P = inv3(b, N3)
        ar.erase_inverse(b, P, N3)
        b.free(N3)
        b.cx(M[3:], M[:3])
        masks = []
        for half in (M[:3], M[3:]):
            t = b.alloc_qubits(3)
            mul(b, half, P, t)
            m = b.alloc_bits(3)
            b.hmr(half, m)
            masks.append(m)
            b.swap(half, t)
            b.free(t)
        ar._erase_norm(b, M, P, K6)
        b.free(P)
        ar.linear_inplace(b, M, FRM6)
        zero = b.alloc_bits(6)
        b.bit_store(zero, False)
        return zero, masks, F
    return produce, rebuild


def _edits_noflag_latebody(r6n, latem):
    """r6n's UNFLAGGED + latem's late M, with the body taken out of its place (the caller appends it at the end)."""
    assert sum(1 for b_, _ in r6n.UNFLAGGED if b_ == r6n.UNPATCH) == 1
    out = [(before, "" if before == r6n.UNPATCH else after) for before, after in r6n.UNFLAGGED]
    return tuple(out) + tuple(latem.LATE_M)


def install(km, module):
    if not ON:
        return dict(on=False)
    r6n = sys.modules.get('dqi_build_patch_r6n')
    latem = sys.modules.get('dqi_junc_r8p_latem')
    assert r6n is not None and latem is not None, 'install after r6n and latem'
    div6 = __import__('dqi_junc_r8p_div6')
    M = sys.modules.get('dqi_pub_a1model')
    assert M is not None, 'install after a1 (dqi_pub_a1model not loaded)'
    walk = sys.modules['walk_bound']
    mod = sys.modules['norm_deferred_flag_return']
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    assert src.count(BODY) == 1 and src.endswith("    b.free(N)\n"), 'unexpected divide_from_norms text'
    tail = BODY if LATEBODY else ""
    lf = r6n._variant(mod, src + tail, _edits(r6n, latem), 'divide_flag_lf')
    nf_lb = r6n._variant(mod, src + BODY, _edits_noflag_latebody(r6n, latem), 'divide_noflag_lb') if LATEBODY else None
    sn = next(m for m in list(sys.modules.values())
              if 'rank_affine' in getattr(m, '__dict__', {}) and 'projected_product' in m.__dict__)
    rank_affine = sn.rank_affine
    common = next(m for m in list(sys.modules.values())
                  if 'linear_words' in getattr(m, '__dict__', {}) and 'matrix_ops' in m.__dict__)
    cx = common.cx
    zero_flag = sys.modules['exp.arithmetic.division'].zero_flag
    cache = {}
    KD = __import__('dqi_kern_div') if BETA0 else None

    def dr(b, x, y, fm, after, owner, K):
        """B6's dr division (tools/dqi_kern_div.divide_dr_dec, unflagged path, read only: norm 42, product 99, cached
        inverse 58, halves 84, inverse phase 27) with the caller's body AFTER u and N are erased: the body allocates
        beta, and the kit's phase_inverse36_pivot takes its virtual form (27 ANDs with pub_inv27) only when live + 15
        <= 889, a literal of an older ceiling: B6's form with the body first measured 320.2 here [M lf5], not 310.
        The virtual form is called directly (876 live + 15 scratch). The body reads y and does not change it."""
        from exp.jet24.pro_jet import cached_inverse, PLAN12
        from exact_subs import phase_inverse12_virtual  # the form the kit's selector takes below its width line
        g = owner.__globals__
        ar = g['ar']
        N = b.alloc_qubits(12)
        g['norm_xor_noscratch'](b, y, N, K)            # N = Norm(y)
        b.cx(y[12:], y[:12])
        oldmask = g['multiply_forward99'](b, x, y, K)   # x = x conj(y); old x = new x * y
        b.cx(y[12:], y[:12])
        u, Fz = cached_inverse(b, N, PLAN12, True)
        b.free(Fz)                                      # y != 0 on an unflagged fold (policy W)
        KD._halves_times(b, g, x, u, N)                 # x = x conj(y) / N
        g['phase_product24'](b, x, y, oldmask, K)
        b.free(oldmask)
        um = b.alloc_bits(12)
        b.hmr(u, um)
        b.free(u)
        phase_inverse12_virtual(b, N, um)               # u = 1 / N erased against N (27 with pub_inv27; 15 scratch)
        b.free(um)
        ar._erase_norm(b, y, N, K)
        b.free(N)
        after()
        b.free(fm)
    kmod = getattr(common, 'km', None) or sys.modules['kickmix']
    reclaim = common.reclaim if hasattr(common, 'reclaim') else sys.modules['exp.decoder.common'].reclaim

    def home(b, reg, wires):
        """Move the 12 qubits of `reg` onto the recorded wires `wires` (swaps only); returns the register there."""
        on = {str(q) for q in reg}
        target = {str(q) for q in wires}
        need = [q for q in wires if str(q) not in on]              # recorded wires that do not hold the word: free
        move = [q for q in reg if str(q) not in target]            # qubits of the word outside the recorded wires
        assert len(need) == len(move) and len(reg) == len(wires) == 12
        if not need:
            return reg
        need_r, move_r = kmod.array(need), kmod.array(move)
        reclaim(b, [need_r])
        b.swap(move_r, need_r)
        b.free(move_r)
        it = iter(need)
        return kmod.array([q if str(q) in target else next(it) for q in reg])

    def forms(m, n, k):
        key = (m, n, k)
        if key not in cache:
            cache[key] = FORMS.crt_forms_inf(M, m, n, k, FORMS.kpoints(n, k) if KPOINT else None)
        return cache[key]

    def premap(b, Uw, T, k):
        pre, _ = forms(len(T) + 1, len(T), k)
        assert pre[-1] == [0, 0]
        for i, (c0, c1) in enumerate(pre):
            if c0:
                cx(b, Uw[0], T[i], c0)
            if c1:
                cx(b, Uw[1], T[i], c1)
        STATS['pre'] += 1

    def body(b, head, source, target, k):
        _, ups = forms(len(source), len(target), k)
        for ac, bc, tc in ups:
            rank_affine(b, head, ac, source, bc, target, tc)
        STATS['bodies'] += 1

    made_div = {}

    def wrap(prev):
        if id(prev) in made_div:
            return made_div[id(prev)]
        cells = {}
        cells_u = {}

        def divide(b, x, y, fm, before, after, owner, K=None):
            if before is NF_MARK:
                if DIV6ALL:
                    # width only: the unflagged kernel's retained-19 inverse holds 14 wires beside M and N (26
                    # helpers, tight 902); div6's in-place GF(64) inverse holds 12 (+5 Toffolis per call)
                    Gu = owner.__globals__
                    assert 'produce_retained' in Gu and 'rebuild_retained' in Gu, 'divider without the retained hooks'
                    if id(Gu) not in cells_u:
                        cells_u[id(Gu)] = div6._make(Gu)
                    saved_u = Gu['produce_retained'], Gu['rebuild_retained']
                    Gu['produce_retained'], Gu['rebuild_retained'] = cells_u[id(Gu)]
                    STATS['div6all'] = STATS.get('div6all', 0) + 1
                    try:
                        if nf_lb is None:
                            STATS['passed'] += 1
                            return prev(b, x, y, fm, None, after, owner, K)
                        STATS['latebody'] += 1
                        return nf_lb(b, x, y, fm, None, after, owner, K)
                    finally:
                        Gu['produce_retained'], Gu['rebuild_retained'] = saved_u
                if nf_lb is None:
                    STATS['passed'] += 1
                    return prev(b, x, y, fm, None, after, owner, K)
                STATS['latebody'] += 1
                return nf_lb(b, x, y, fm, None, after, owner, K)
            if not isinstance(before, LateFlag):
                STATS['passed'] += 1
                return prev(b, x, y, fm, before, after, owner, K)
            STATS['lateflag'] += 1
            G = owner.__globals__
            assert 'produce_retained' in G and 'rebuild_retained' in G, 'divider without the retained hooks'
            if id(G) not in cells:
                cells[id(G)] = _make_n3(G, div6, zero_flag) if N3FLAG else div6._make(G)
            saved = G['produce_retained'], G['rebuild_retained']
            G['produce_retained'], G['rebuild_retained'] = cells[id(G)]
            try:
                return lf(b, x, y, fm, before, after, owner, K)
            finally:
                G['produce_retained'], G['rebuild_retained'] = saved
        divide._r48_lateflag = True
        made_div[id(prev)] = divide
        return divide

    made = {}
    done = 0
    for fn in walk.functions(module.decode):
        g = fn.__globals__
        s = g.get('sweep')
        if fn.__name__ != 'consume_synthesis' or s is None or not hasattr(s, '__code__'):
            continue
        if not s.__code__.co_filename.endswith('dqi_int_a1.py'):
            continue
        if id(s) not in made:
            env = dict(s.__globals__)
            assert env.get('_B5_MODE') == 'kitco', env.get('_B5_MODE')
            prev = env['divide_from_norms']
            assert getattr(prev, '_dqi_int_div6f', False), 'expected div6f\'s dispatcher under the a1 sweep'
            cl_xor = env['previous'].controlled_linear_xor

            def make_lf(b, lead, T, k, degree, ntarget, _cl=cl_xor):
                c0, c1, corr = FORMS.completion(M, k, degree, ntarget, len(T),
                                                FORMS.kpoints(len(T), k) if KPOINT else None)
                return LateFlag(b, lead, T, c0, c1, corr, _cl, M.mul)
            env.update(_R48_PRE=premap, _R48_BODY=body, _R48_LF=make_lf, _R48_NF=NF_MARK, divide_from_norms=wrap(prev),
                       _R48_BETA0=BETA0, _R48_DR=dr, _R48_STATS=STATS, _R48_HOME=home)
            exec(compile(SWEEP_SRC, __file__, 'exec'), env)
            made[id(s)] = env['sweep']
        g['sweep'] = made[id(s)]
        done += 1
    assert done, 'no live consume_synthesis with B5\'s a1 sweep found'
    STATS['sweeps'] = len(made)
    atexit.register(lambda: print('r48_dqi_ret2_lateflag', STATS, file=sys.stderr, flush=True))
    return dict(on=True, installed=done, sweeps=len(made), latebody=LATEBODY, kpoint=KPOINT, n3flag=N3FLAG,
                beta0=BETA0)
