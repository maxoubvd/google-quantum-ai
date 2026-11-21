"""dqi_kern_div.py -- 24-bit divisions x -> x / y (y kept) for un-synthesis (workstream dqi_kern, builder B6, D3).

All variants are built from the kit's own arithmetic (the fold divider's globals g: ar, norm_xor_noscratch,
mul_xor_noscratch, multiply_forward99, phase_product24, phase_inverse6_classical6) and are tested by
tools/dqi_kern_divtest.py. Field: GF(2^24) = GF(2^12)[w] with the kit's extension constant K; 1/y = conj(y) / N,
N = Norm(y) in GF(2^12); GF(2^12) = GF(64)[.] (ar.tower(12)), 1/N = conj6(N) / M with M = Norm6(N) in GF(64).

tw2  (tower, M kept): product x conj(y) (99) -> N (42, after the product) -> halves times conj6(N) (84, masks
     paid at the end against N) -> M (15) -> t = 1/M out of place (21, a GF(8)-tower inverse with 6 scratch;
     M kept) -> the four 6-bit parts times t (60, masks paid at once against M) -> t erased by the kit's
     classical 6-bit inverse phase (M live) -> M, N erased free. No rebuild of M. Width: 24 at the product,
     N 12 + target 12 at the halves, N 12 + M 6 + t 6 + target 6 = 30 at the parts.
"""
import sys

VARIANTS = {}


def variant(name, flagged=False):
    def reg(fn):
        fn.flagged = flagged
        VARIANTS[name] = fn
        return fn
    return reg


def inv6_oop(b, g, M):
    """t = 1/M (6 bits, M != 0) out of place, M kept: GF(64) = GF(8)^2 tower, Norm3 6 + inverse3 3 + two GF(8)
    products 12 = 21 Toffolis, at most 6 scratch qubits. Returns t (6 fresh qubits)."""
    ar = g['ar']
    norm = g['norm_xor_noscratch']
    mul = g['mul_xor_noscratch']
    inv3 = sys.modules['exp.noscratch.kernels'].inverse3_noscratch
    to6, frm6, k3 = ar.tower(6)
    ar.linear_inplace(b, M, to6)
    N3 = b.alloc_qubits(3)
    norm(b, M, N3, k3)                              # Norm3(M)
    P = inv3(b, N3)                                 # P = 1/Norm3(M)
    ar.erase_inverse(b, P, N3)                      # N3 = 1/P: free
    b.free(N3)
    t = b.alloc_qubits(6)
    b.cx(M[3:], M[:3])                              # conj3(M)
    mul(b, M[:3], P, t[:3])
    mul(b, M[3:], P, t[3:])                         # t = conj3(M) / Norm3(M) = 1/M (tower coordinates)
    b.cx(M[3:], M[:3])
    ar._erase_norm(b, t, P, k3)                     # P = Norm3(t): free
    b.free(P)
    ar.linear_inplace(b, M, frm6)
    ar.linear_inplace(b, t, frm6)
    return t


def inv6_iso20(b, g, M):
    """t = 1/M (6 bits) out of place, M kept: the kit's width-preserving 20-AND inverse
    (exp/inverse_search/inverse20_iso.inverse6_xor20, fresh output: 3 scratch). Returns t."""
    iso = sys.modules.get('exp.inverse_search.inverse20_iso')
    if iso is None:
        iso = __import__('exp.inverse_search.inverse20_iso', fromlist=['inverse6_xor20'])
    return iso.inverse6_xor20(b, M)


@variant('tw3')
def divide_tw3(b, x, y, K, g):
    """tw2 with the 20-AND 6-bit inverse (inv6_iso20) in place of inv6_oop (21)."""
    divide_tw2(b, x, y, K, g, inv6=inv6_iso20)


@variant('tw2')
def divide_tw2(b, x, y, K, g, inv6=None):
    ar = g['ar']
    norm = g['norm_xor_noscratch']
    mul = g['mul_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    b.cx(y[12:], y[:12])                            # conj(y)
    oldmask = g['multiply_forward99'](b, x, y, K)   # x = x conj(y); mask on old x = x'' y
    b.cx(y[12:], y[:12])
    N = b.alloc_qubits(12)
    norm(b, y, N, K)                                # N = Norm(y)
    conj = [ar.field(12).frobenius(1 << j, 6) for j in range(12)]
    ar.linear_inplace(b, N, conj)                   # conj6(N)
    masks = []
    for off in (0, 12):
        half = x[off:off + 12]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, half, N, out)
        m = b.alloc_bits(12)
        b.hmr(half, m)                              # old half = final half * N: paid at the end
        b.swap(half, out)
        b.free(out)
        masks.append(m)
    ar.linear_inplace(b, N, conj)                   # N again
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)                               # M = Norm6(N)
    ar.linear_inplace(b, N, frm)
    t = (inv6 or inv6_oop)(b, g, M)
    for off in (0, 12):
        half = x[off:off + 12]
        ar.linear_inplace(b, half, to)
        for j in (0, 6):
            part = half[j:j + 6]
            out = b.alloc_qubits(6)
            mul(b, part, t, out)
            ar.erase_product(b, out, M, part)       # old part = new part * M, M live: free
            b.swap(part, out)
            b.free(out)
        ar.linear_inplace(b, half, frm)
    tm = b.alloc_bits(6)
    b.hmr(t, tm)
    b.free(t)
    g['phase_inverse6_classical6'](b, M, tm, zero_mask=tm[0])
    b.free(tm)
    ar.linear_inplace(b, N, to)
    ar._erase_norm(b, N, M, k6)                     # M = Norm6(N): free
    b.free(M)
    ar.linear_inplace(b, N, frm)
    for off, m in zip((0, 12), masks):
        ar.phase_product(b, x[off:off + 12], N, m)
        b.free(m)
    ar._erase_norm(b, y, N, K)                      # N = Norm(y): free
    b.free(N)
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)


def divide_dr_dec(b, x, y, fm, before, after, owner, K=None, flagged=True):
    """dr in the decoder's interface (that of norm_deferred_flag_return.divide_from_norms: x -> x / y, y kept;
    flagged: y = 0 allowed with x = 0, the flag F = [y = 0] is handed to `before`, retired into the caller's
    mask bit fm as the kit does). Norm FIRST (the flag is a 12-bit zero test, 11 ANDs, instead of r6n's 24-bit
    one), y made safe (y = 1 when y = 0), x = x conj(y) (99), u = 1/N by the cached inverse (58, N kept), halves
    times u (84, masks paid at once against N), F retired FREE: with N restored to Norm(y) and u turned into
    inverse0(N), bit 0 of N u is [y != 0]; then u erased by the kit's inverse phase (N live).
    Width: N 12 + u 12 + a 12-qubit target (36) + F (37 flagged)."""
    from exp.jet24.pro_jet import cached_inverse, PLAN12
    from exp.mul99.phase36 import phase_inverse36_pivot
    from exp.arithmetic.division import zero_flag
    g = owner if isinstance(owner, dict) else owner.__globals__
    ar = g['ar']
    N = b.alloc_qubits(12)
    g['norm_xor_noscratch'](b, y, N, K)            # N = Norm(y)
    if flagged:
        F = zero_flag(b, N)                         # [y = 0]
        before(F)
        b.cx(F, y[0])                               # y safe
        b.cx(F, N[0])                               # N = Norm(y safe) (Norm(1) = 1)
    b.cx(y[12:], y[:12])
    oldmask = g['multiply_forward99'](b, x, y, K)   # x = x conj(y); old x = new x * y
    b.cx(y[12:], y[:12])
    u, Fz = cached_inverse(b, N, PLAN12, True)
    b.free(Fz)                                      # N != 0 here
    _halves_times(b, g, x, u, N)                    # x = x conj(y) / N
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)
    if flagged:
        b.cx(F, y[0])
        b.cx(F, N[0])                               # N = Norm(y)
        after()
        b.cx(F, u[0])                               # u = inverse0(N) (u was 1 when y = 0)
        fmask = b.alloc_bits(1)
        b.hmr(F, fmask[0])
        b.free(F)
        b.cx(fm[0], fmask[0])
        b.free(fm)
        b.z(fmask[0])                               # (-1)^fmask ...
        pm = b.alloc_bits(12)
        b.bit_store(pm, False)
        b.cx(fmask[0], pm[0])
        ar.phase_product(b, N, u, pm)               # ... times (-1)^(fmask [y != 0]): bit 0 of N u
        b.free(pm)
        b.free(fmask)
    else:                                           # r6n's unflagged contract: after(), the caller's fm consumed
        after()                                     # (its phase fm [y = 0] is 1: y != 0 on these folds)
        b.free(fm)
    um = b.alloc_bits(12)
    b.hmr(u, um)
    b.free(u)
    phase_inverse36_pivot(b, N, um)                 # u = inverse0(N)
    b.free(um)
    ar._erase_norm(b, y, N, K)
    b.free(N)


def divide_tw3_dec(b, x, y, fm, before, after, owner, K=None, flagged=True):
    """tw3 in the decoder's interface (that of divide_dr_dec): x -> x / y, y kept; flagged: y = 0 allowed, and
    then x is kept (x ^ the caller's correction), as in the decoder's own flagged form. Flagged: F = [y = 0] from
    the 24-bit divisor BEFORE the product (23 ANDs; the product must see the safe y so that x is kept), before(F),
    y made safe. Product (24-qubit target), N = Norm(y) (safe), halves times conj6(N) (masks paid against N),
    M = Norm6(N), t = 1/M (20-AND inverse), the four 6-bit parts times t (masks paid at once against M), old x paid
    against the safe y; y, N, M restored; F retired free through bit 0 of M t (t stays 1/M_safe; the kit's zero
    term in the 6-bit inverse phase accounts for it); t erased; M, N free. Width 30 (+F: 31 flagged)."""
    from exp.arithmetic.division import zero_flag
    g = owner if isinstance(owner, dict) else owner.__globals__
    ar = g['ar']
    norm = g['norm_xor_noscratch']
    mul = g['mul_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    if flagged:
        F = zero_flag(b, y)                         # [y = 0]
        before(F)
        b.cx(F, y[0])                               # y safe
    b.cx(y[12:], y[:12])
    oldmask = g['multiply_forward99'](b, x, y, K)   # x = x conj(y); old x = x_final * y (safe y)
    b.cx(y[12:], y[:12])
    N = b.alloc_qubits(12)
    norm(b, y, N, K)                                # N = Norm(y) (1 when y was 0)
    conj = [ar.field(12).frobenius(1 << j, 6) for j in range(12)]
    ar.linear_inplace(b, N, conj)
    masks = []
    for off in (0, 12):
        half = x[off:off + 12]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, half, N, out)
        m = b.alloc_bits(12)
        b.hmr(half, m)                              # old half = final half * N
        b.swap(half, out)
        b.free(out)
        masks.append(m)
    ar.linear_inplace(b, N, conj)
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)                               # M = Norm6(N)
    ar.linear_inplace(b, N, frm)
    t = inv6_iso20(b, g, M)                         # t = 1/M
    for off in (0, 12):
        half = x[off:off + 12]
        ar.linear_inplace(b, half, to)
        for j in (0, 6):
            part = half[j:j + 6]
            out = b.alloc_qubits(6)
            mul(b, part, t, out)
            ar.erase_product(b, out, M, part)       # old part = new part * M: free
            b.swap(part, out)
            b.free(out)
        ar.linear_inplace(b, half, frm)
    for off, m in zip((0, 12), masks):
        ar.phase_product(b, x[off:off + 12], N, m)
        b.free(m)
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)
    if flagged:
        b.cx(F, y[0])                               # y restored
        b.cx(F, N[0])                               # N = Norm(y)
        b.cx(F, M[0])                               # M = Norm6(N)
        after()
        fmask = b.alloc_bits(1)
        b.hmr(F, fmask[0])
        b.free(F)
        b.cx(fm[0], fmask[0])
        b.free(fm)
        b.z(fmask[0])                               # (-1)^fmask ...
        pm = b.alloc_bits(6)
        b.bit_store(pm, False)
        b.cx(fmask[0], pm[0])
        ar.phase_product(b, M, t, pm)               # ... times (-1)^(fmask [y != 0]): bit 0 of M t
        b.free(pm)
        b.free(fmask)
    else:
        after()
        b.free(fm)
    tm = b.alloc_bits(6)
    b.hmr(t, tm)
    b.free(t)
    g['phase_inverse6_classical6'](b, M, tm, zero_mask=tm[0])   # t = 1/M_safe: 1 when M = 0
    b.free(tm)
    ar.linear_inplace(b, N, to)
    ar._erase_norm(b, N, M, k6)                     # M = Norm6(N): free
    b.free(M)
    ar.linear_inplace(b, N, frm)
    ar._erase_norm(b, y, N, K)                      # N = Norm(y): free
    b.free(N)


def test_before(b, x):
    """The caller's before(F) in tests: XOR the constant DQI_KERN_BEFORE into x when F = 1 (the decoder's
    correct_first / correct_second add a constant to words that include x); 0 = no-op."""
    import os
    c = int(os.environ.get('DQI_KERN_BEFORE', '0'))

    def before(F):
        for j in range(24):
            if c >> j & 1:
                b.cx(F, x[j])
    return before


def _dec_test(flagged, fn=None):
    def run(b, x, y, K, g):
        fm = b.alloc_bits(1)
        b.bit_store(fm, False)
        (fn or divide_dr_dec)(b, x, y, fm, test_before(b, x), lambda: None, g, K, flagged)
    return run


VARIANTS['twf'] = _dec_test(True, divide_tw3_dec)
VARIANTS['twf'].flagged = True
VARIANTS['twu'] = _dec_test(False, divide_tw3_dec)
VARIANTS['twu'].flagged = False


VARIANTS['drf'] = _dec_test(True)
VARIANTS['drf'].flagged = True
VARIANTS['dru'] = _dec_test(False)
VARIANTS['dru'].flagged = False


def _product_and_norm(b, x, y, K, g):
    """x = x conj(y) (oldmask on the old x = x'' y, paid at the end), N = Norm(y) after the product."""
    b.cx(y[12:], y[:12])
    oldmask = g['multiply_forward99'](b, x, y, K)
    b.cx(y[12:], y[:12])
    N = b.alloc_qubits(12)
    g['norm_xor_noscratch'](b, y, N, K)
    return oldmask, N


def _halves_times(b, g, x, u, N):
    """x halves times u (12-bit target each); old half = new half * N with N live: paid at once."""
    ar = g['ar']
    for off in (0, 12):
        half = x[off:off + 12]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, half, u, out)
        ar.erase_product(b, out, N, half)
        b.swap(half, out)
        b.free(out)


def _finish(b, x, y, K, g, N, oldmask):
    ar = g['ar']
    ar._erase_norm(b, y, N, K)
    b.free(N)
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)


@variant('dr')
def divide_direct(b, x, y, K, g):
    """Direct form: u = 1/N out of place by the cached inverse (exp/jet24 PLAN12, N kept), halves times u with
    their masks paid at once, u erased by the kit's classical 12-bit inverse phase (N live). Width 36."""
    from exp.jet24.pro_jet import cached_inverse, PLAN12
    from exp.mul99.phase36 import phase_inverse36_pivot
    oldmask, N = _product_and_norm(b, x, y, K, g)
    u, F = cached_inverse(b, N, PLAN12, True)
    b.free(F)                                       # N != 0: the flag is identically 0
    _halves_times(b, g, x, u, N)
    um = b.alloc_bits(12)
    b.hmr(u, um)
    b.free(u)
    phase_inverse36_pivot(b, N, um)
    b.free(um)
    _finish(b, x, y, K, g, N, oldmask)


@variant('dt42')
def divide_direct_tower42(b, x, y, K, g):
    """dt with t erased right after u is formed (t = Norm6(u): free) and recomputed (15) to erase u at the end:
    width 42 (N 12 + M 6 + u 12 + target 12)."""
    ar = g['ar']
    mul = g['mul_xor_noscratch']
    norm = g['norm_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    oldmask, N = _product_and_norm(b, x, y, K, g)
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)
    ar.linear_inplace(b, N, frm)
    t = inv6_oop(b, g, M)
    ar.linear_inplace(b, N, to)
    b.cx(N[6:], N[:6])
    u = b.alloc_qubits(12)
    mul(b, N[:6], t, u[:6])
    mul(b, N[6:], t, u[6:])                         # u = conj6(N) t, tower coordinates
    ar._erase_norm(b, u, t, k6)                     # t = Norm6(u): free
    b.free(t)
    b.cx(N[6:], N[:6])
    ar.linear_inplace(b, N, frm)
    ar.linear_inplace(b, u, frm)
    _halves_times(b, g, x, u, N)
    ar.linear_inplace(b, u, to)
    t = b.alloc_qubits(6)
    norm(b, u, t, k6)                               # t = Norm6(u) again (15)
    ar.linear_inplace(b, N, to)
    b.cx(N[6:], N[:6])
    ar.erase_product(b, N[:6], t, u[:6])
    ar.erase_product(b, N[6:], t, u[6:])
    b.free(u)
    b.cx(N[6:], N[:6])
    ar.linear_inplace(b, N, frm)
    tm = b.alloc_bits(6)
    b.hmr(t, tm)
    b.free(t)
    g['phase_inverse6_classical6'](b, M, tm, zero_mask=tm[0])
    b.free(tm)
    ar.linear_inplace(b, N, to)
    ar._erase_norm(b, N, M, k6)
    b.free(M)
    ar.linear_inplace(b, N, frm)
    _finish(b, x, y, K, g, N, oldmask)


@variant('dt')
def divide_direct_tower(b, x, y, K, g):
    """Direct form with 1/N from the GF(64) tower and kept factors: M = Norm6(N), t = 1/M (out of place, M kept),
    u = conj6(N) t; halves times u (masks paid at once with N); u erased FREE (bilinear in conj6(N), t); t erased
    by the classical 6-bit inverse phase (M live). Width 48 (N 12 + M 6 + t 6 + u 12 + target 12)."""
    ar = g['ar']
    mul = g['mul_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    oldmask, N = _product_and_norm(b, x, y, K, g)
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    g['norm_xor_noscratch'](b, N, M, k6)            # M = Norm6(N)
    ar.linear_inplace(b, N, frm)
    t = inv6_oop(b, g, M)
    ar.linear_inplace(b, N, to)
    b.cx(N[6:], N[:6])                              # conj6(N) in tower coordinates
    u = b.alloc_qubits(12)
    mul(b, N[:6], t, u[:6])
    mul(b, N[6:], t, u[6:])                         # u = conj6(N) t = 1/N (tower coordinates)
    b.cx(N[6:], N[:6])
    ar.linear_inplace(b, N, frm)
    ar.linear_inplace(b, u, frm)
    _halves_times(b, g, x, u, N)
    ar.linear_inplace(b, u, to)
    ar.linear_inplace(b, N, to)
    b.cx(N[6:], N[:6])
    ar.erase_product(b, N[:6], t, u[:6])            # free: conj6(N) and t live
    ar.erase_product(b, N[6:], t, u[6:])
    b.free(u)
    b.cx(N[6:], N[:6])
    ar.linear_inplace(b, N, frm)
    tm = b.alloc_bits(6)
    b.hmr(t, tm)
    b.free(t)
    g['phase_inverse6_classical6'](b, M, tm, zero_mask=tm[0])
    b.free(tm)
    ar.linear_inplace(b, N, to)
    ar._erase_norm(b, N, M, k6)
    b.free(M)
    ar.linear_inplace(b, N, frm)
    _finish(b, x, y, K, g, N, oldmask)
