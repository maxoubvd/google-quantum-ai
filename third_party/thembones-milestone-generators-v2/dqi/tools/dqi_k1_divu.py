"""dqi_k1_divu.py -- K1 (round 3, door 2 fix rounds): division forms x -> x / y (GF(2^24), y kept) that keep the
divisor's 12-bit inverse u = 1/N as the halves' multiplier without holding N beside it (the width door at 26).
Built from the kit's own arithmetic exactly as tools/dqi_kern_div.py does (g = the fold divider's globals); run
through tools/dqi_k1_divprof.py (forms named here are looked up in VARIANTS).

u26: product x conj(y) (99) -> N = Norm(y) (42) -> M = Norm6(N) (15) -> t = 1/M out of place (20, M then erased free
     against N) -> N's register turned into u = conj6(N) t IN PLACE (2 x 15; the old parts are conj6(N)'s, paid at the
     end against the recomputed N) -> t erased free (t = Norm6(u)) -> halves x u (84; old half = new half * N, paid
     at the end) -> N recomputed (42) -> masks paid (free) -> u erased by the kit's 12-bit inverse phase (N live)
     -> N erased free (Norm(y)) -> old x paid (free). Width: 27 at the 6-bit inverse, 24 at the halves.
     Counted 99 + 42 + 15 + 20 + 30 + 84 + 42 + 32 = 364 (K1 log, round 3: 'N needed twice').
"""
import sys

VARIANTS = {}


def variant(name, flagged=False):
    def reg(fn):
        fn.flagged = flagged
        VARIANTS[name] = fn
        return fn
    return reg


@variant('tw3g')
def divide_tw3g(b, x, y, K, g, Gw):
    """tw3 (B6's tw2 with the 20-AND inverse) with M PARKED ON BORROWED WIRES of G during the four part scalings:
    G6 ^= M (CNOT), M's clean wires erased free against N; each part's erasure phase is taken against G6 = g ^ M and
    its g-term cancelled after G6 is restored; M recomputed clean (15) to restore G6 and to erase t. Counted: tw3 +
    15 = 338.9; width: the inverse still needs 27 (N + M + t + 3), the parts 24 instead of 30."""
    import dqi_kern_div as D
    ar = g['ar']
    norm = g['norm_xor_noscratch']
    mul = g['mul_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    G6 = Gw[:6]
    b.cx(y[12:], y[:12])
    oldmask = g['multiply_forward99'](b, x, y, K)
    b.cx(y[12:], y[:12])
    N = b.alloc_qubits(12)
    norm(b, y, N, K)
    conj = [ar.field(12).frobenius(1 << j, 6) for j in range(12)]
    ar.linear_inplace(b, N, conj)
    masks = []
    for off in (0, 12):
        half = x[off:off + 12]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, half, N, out)
        m = b.alloc_bits(12)
        b.hmr(half, m)
        b.swap(half, out)
        b.free(out)
        masks.append(m)
    ar.linear_inplace(b, N, conj)
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)
    t = D.inv6_iso20(b, g, M)
    b.cx(M, G6)                                         # park: G6 = g ^ M
    ar._erase_norm(b, N, M, k6)                         # M's clean wires freed (M = Norm6(N), N live)
    b.free(M)
    ar.linear_inplace(b, N, frm)
    pmasks = []
    for off in (0, 12):
        half = x[off:off + 12]
        ar.linear_inplace(b, half, to)
        for j in (0, 6):
            part = half[j:j + 6]
            out = b.alloc_qubits(6)
            mul(b, part, t, out)
            pm = b.alloc_bits(6)
            b.hmr(part, pm)                             # old part = new part * M
            b.swap(part, out)
            b.free(out)
            ar.phase_product(b, part, G6, pm)           # against g ^ M; the g-term is cancelled below
            pmasks.append(pm)
        ar.linear_inplace(b, half, frm)
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)                                   # M again (15)
    b.cx(M, G6)                                         # G6 = g restored
    for off in (0, 12):
        half = x[off:off + 12]
        ar.linear_inplace(b, half, to)
        for j in (0, 6):
            pm = pmasks.pop(0)
            ar.phase_product(b, half[j:j + 6], G6, pm)  # cancel the g-term
            b.free(pm)
        ar.linear_inplace(b, half, frm)
    tm = b.alloc_bits(6)
    b.hmr(t, tm)
    b.free(t)
    g['phase_inverse6_classical6'](b, M, tm, zero_mask=tm[0])
    b.free(tm)
    ar._erase_norm(b, N, M, k6)
    b.free(M)
    ar.linear_inplace(b, N, frm)
    for off, m in zip((0, 12), masks):
        ar.phase_product(b, x[off:off + 12], N, m)
        b.free(m)
    ar._erase_norm(b, y, N, K)
    b.free(N)
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)


divide_tw3g.borrowed = 24


@variant('u26')
def divide_u26(b, x, y, K, g):
    from exp.mul99.phase36 import phase_inverse36_pivot
    import dqi_kern_div as D
    ar = g['ar']
    norm = g['norm_xor_noscratch']
    mul = g['mul_xor_noscratch']
    to, frm, k6 = ar.tower(12)
    # 1. x <- x conj(y); the old x = x_final * y is paid at the end (bilinear)
    b.cx(y[12:], y[:12])
    oldmask = g['multiply_forward99'](b, x, y, K)
    b.cx(y[12:], y[:12])
    # 2. N = Norm(y)
    N = b.alloc_qubits(12)
    norm(b, y, N, K)
    # 3. M = Norm6(N); 4. t = 1/M out of place; M erased free against N
    ar.linear_inplace(b, N, to)
    M = b.alloc_qubits(6)
    norm(b, N, M, k6)
    ar.linear_inplace(b, N, frm)
    t = D.inv6_iso20(b, g, M)
    ar.linear_inplace(b, N, to)
    ar._erase_norm(b, N, M, k6)
    b.free(M)
    # 5. N's register -> u = conj6(N) t in place (tower coordinates), old parts X-measured (deferred)
    b.cx(N[6:], N[:6])                                  # conj6(N) in tower coordinates
    nmasks = []
    for j in (0, 6):
        part = N[j:j + 6]
        out = b.alloc_qubits(6)
        mul(b, part, t, out)
        m = b.alloc_bits(6)
        b.hmr(part, m)
        b.swap(part, out)
        b.free(out)
        nmasks.append(m)
    # 6. t = Norm6(u): erased free
    ar._erase_norm(b, N, t, k6)
    b.free(t)
    ar.linear_inplace(b, N, frm)
    u = N                                               # u = 1/N, standard coordinates
    # 7. halves x u; old half = new half * N (deferred)
    hmasks = []
    for off in (0, 12):
        half = x[off:off + 12]
        out = b.alloc_qubits(12)
        ar.mul_xor0(b, half, u, out)
        m = b.alloc_bits(12)
        b.hmr(half, m)
        b.swap(half, out)
        b.free(out)
        hmasks.append(m)
    # 8. N again; the deferred masks; u erased by the inverse phase; N erased; old x paid
    N2 = b.alloc_qubits(12)
    norm(b, y, N2, K)
    for off, m in zip((0, 12), hmasks):
        ar.phase_product(b, x[off:off + 12], N2, m)
        b.free(m)
    ar.linear_inplace(b, N2, to)
    b.cx(N2[6:], N2[:6])
    for j, m in zip((0, 6), nmasks):
        b.cz(m, N2[j:j + 6])
        b.free(m)
    b.cx(N2[6:], N2[:6])
    ar.linear_inplace(b, N2, frm)
    um = b.alloc_bits(12)
    b.hmr(u, um)
    b.free(u)
    phase_inverse36_pivot(b, N2, um)
    b.free(um)
    ar._erase_norm(b, y, N2, K)
    b.free(N2)
    g['phase_product24'](b, x, y, oldmask, K)
    b.free(oldmask)
