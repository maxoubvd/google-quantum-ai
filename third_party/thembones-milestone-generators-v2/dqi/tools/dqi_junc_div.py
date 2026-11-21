"""dqi_junc_div.py -- in-place 12-bit division in a 12-qubit window (workstream dqi_junc, synthesis part of D5).

  inv6_inplace(b, N) -> masks          N -> 1/N in place in GF(64) (0 -> 0): 21 Toffolis, 6 private qubits.
        Through GF(8) coordinates: N' = Norm(N) 6, P = 1/N' 3, N' erased free (quadratic in P), conj(N)'s two
        coordinates multiplied by P in place through a 3-qubit target (2 x 6), P = Norm(new N) erased free.
        The X-measured old coordinates leave a phase LINEAR in the old N: repay_inv6(b, N_old_register, masks).
  divide_narrow(b, x, y)               x -> x / y in place, y kept (promise y != 0 or x = 0).
        1. T = x conj(y) into a 12-qubit target (42); x X-measured (paid at the end: x_old = x_new y, bilinear);
        2. N = Norm(y) (15), N -> M = 1/N in place (21);
        3. T's two GF(64) coordinates times M in place through a 6-qubit target (2 x 15): T M = x / y;
        4. M X-measured; N = Norm(y) rebuilt (15) to pay every mask: M by the kit's GF(64) inverse phase
           (phase_inverse6_classical6, <= 6 executed), the inversion's masks (linear in N), the coordinates'
           masks (T_i = x_i N, bilinear), then N erased (quadratic in y).
        Toffolis 42 + 15 + 21 + 30 + 15 + 5.9 = 128.9 executed; window 12.
  divide_zs(b, x, y, dy=None)          zero-safe: a zero divisor acts as 1 (G = [y = 0] added to y's bit 0);
        with dy (a register holding y's derivative), first y ^= F dy with F = [y = 0] (the derivative fallback
        of even levels). G and F are released by X-measurement and the zero phase of y under the measured bit.
        Window 12 + G (+ F). This reproduces x12c's divisor rules (+1 where both nodes are double zeros, +1 where
        parked, derivative where D_t) without holding the fold's flags.
"""
from dqi_junc_kit import gf, ar, norm_xor_noscratch, mul_xor_noscratch, inverse3_noscratch, measured  # noqa: F401
from exp.phase24_tradeoff.phase6small import phase_inverse6_classical6
from exp.arithmetic.division import zero_flag, phase_zero

TO12, FRM12, K12 = ar.tower(12)
TO6, FRM6, K6 = ar.tower(6)
CONJ12 = tuple(gf.pow(1 << j, 64) for j in range(12))          # y -> y^64 (an involution on GF(2^12))


def inv6_inplace(b, N):
    ar.linear_inplace(b, N, TO6)
    N3 = b.alloc_qubits(3)
    norm_xor_noscratch(b, N, N3, K6)          # 6
    P = inverse3_noscratch(b, N3)             # 3
    ar.erase_inverse(b, P, N3)                # N3 = 1/P: quadratic in P, free
    b.free(N3)
    b.cx(N[3:], N[:3])                        # conj(N) = (a + c, c)
    masks = []
    for half in (N[:3], N[3:]):
        t = b.alloc_qubits(3)
        mul_xor_noscratch(b, half, P, t)      # 6 each
        masks.append(measured(b, half))       # phase linear in conj(N_old)
        b.swap(half, t)
        b.free(t)
    ar._erase_norm(b, N, P, K6)               # P = Norm(conj(N_old) P): erased free
    b.free(P)
    ar.linear_inplace(b, N, FRM6)
    return masks


def repay_inv6(b, N, masks):
    """N holds the OLD value again (polynomial basis): pay the phase left by inv6_inplace (free)."""
    ar.linear_inplace(b, N, TO6)
    b.cx(N[3:], N[:3])
    b.cz(N[:3], masks[0])
    b.cz(N[3:], masks[1])
    b.cx(N[3:], N[:3])
    ar.linear_inplace(b, N, FRM6)
    for m in masks:
        b.free(m)


def divide_narrow(b, x, y):
    # 1. T = x conj(y)
    ar.linear_inplace(b, y, CONJ12)
    T = b.alloc_qubits(12)
    ar.mul_xor0(b, x, y, T)                   # 42
    m1 = measured(b, x)                       # x_old = x_new y at the end
    b.swap(x, T)
    b.free(T)
    ar.linear_inplace(b, y, CONJ12)
    # 2. N = Norm(y), M = 1/N in place
    ar.linear_inplace(b, y, TO12)
    N = b.alloc_qubits(6)
    norm_xor_noscratch(b, y, N, K12)          # 15
    ar.linear_inplace(b, y, FRM12)
    m2 = inv6_inplace(b, N)                   # 21; N holds M
    # 3. T M, coordinate by coordinate
    ar.linear_inplace(b, x, TO12)
    m3 = []
    for half in (x[:6], x[6:]):
        t = b.alloc_qubits(6)
        mul_xor_noscratch(b, half, N, t)      # 15 each
        m3.append(measured(b, half))          # T_i = x_new_i N_old
        b.swap(half, t)
        b.free(t)
    ar.linear_inplace(b, x, FRM12)
    # 4. M out, then every mask paid with N rebuilt
    m4 = measured(b, N)                       # (-1)^(m4 . inv6(Norm(y)))
    b.free(N)
    ar.linear_inplace(b, y, TO12)
    N = b.alloc_qubits(6)
    norm_xor_noscratch(b, y, N, K12)          # 15
    phase_inverse6_classical6(b, N, m4)       # <= 6 executed
    b.free(m4)
    repay_inv6(b, N, m2)
    ar.linear_inplace(b, x, TO12)
    ar.phase_product(b, x[:6], N, m3[0])
    ar.phase_product(b, x[6:], N, m3[1])
    ar.linear_inplace(b, x, FRM12)
    for m in m3:
        b.free(m)
    ar._erase_norm(b, y, N, K12)              # N = Norm(y): erased free
    b.free(N)
    ar.linear_inplace(b, y, FRM12)
    ar.phase_product(b, x, y, m1)             # x_old = x_new y
    b.free(m1)


def release_zero(b, F, y):
    """F = [y == 0] with y live: X-measure F, repay with the zero phase of y under the measured bit."""
    m = b.alloc_bits(1)
    b.hmr(F, m[0])
    b.free(F)
    b.push_condition(m[0])
    phase_zero(b, y, m[0])
    b.pop_condition()
    b.free(m)


def divide_zs(b, x, y, dy=None):
    F = None
    if dy is not None:
        F = zero_flag(b, y)                   # 11
        for d, v in zip(dy, y):
            b.ccx(F, d, v)                    # y <- y + F dy (= dy where y = 0)
    G = zero_flag(b, y)                       # 11
    b.cx(G, y[0])                             # a zero divisor acts as 1
    divide_narrow(b, x, y)
    b.cx(G, y[0])
    release_zero(b, G, y)                     # 5.5 executed
    if F is not None:
        for d, v in zip(dy, y):
            b.ccx(F, d, v)
        release_zero(b, F, y)
