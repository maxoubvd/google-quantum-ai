#!/usr/bin/env python3
"""dqi_pp_gvn.py -- workstream dqi_pp (P1): EXACT value numbering of an emitted kickmix .kmb (no sampling).

Every wire and classical bit holds an affine form over GF(2): a constant bit plus an XOR of atoms. Atoms are the
inputs (syndrome qubits q0..q839, mask bits b0..b49139), each HMR outcome, and each nonlinear product. A form is
stored as a 128-bit XOR hash of its atoms' labels plus the constant. A product of two forms is canonicalised on the
2-D span of their linear parts: A*B, A*(A+B) = A + A*B and B*(A+B) = B + A*B share one atom labelled by the two
smallest of {A, B, A+B}; dependent parts (A = 0, B = 0, A = B) give an affine result with no atom. Sound by
construction: two equal hashes are equal Boolean functions up to a 2^-128 collision, so a near-identity (a wire
that differs from a product on rare inputs) can never pass. Incomplete: identities that need deeper algebra are
missed, which only loses sites.

Classes per Toffoli-class op (CCX / CCZ / CSWAP, free or inside PUSH_CONDITION blocks; weight 0.5 per level):
  CCX  1 DEG   product affine in a, b          -> CX / X network (pass 6)
       2 HELD  another wire w holds the product -> CX(w, t) (+ X)  (pass 2)
       3 UNC0  target holds the product (t -> 0)   -> HMR(t); CZ(a, b) if m (pass 3, generalised)
       4 UNC1  target holds product + 1 (t -> 1)   -> HMR(t); CZ(a, b) if m; NEG if m; X(t)
  CCZ  1 DEG   linear parts of a, b, c dependent  -> CZ / Z / NEG
  CSWAP 1 DEG  d = c*(x+y) affine                 -> SWAP / CX / X
        2 HELD another wire w holds d             -> CX(w, x); CX(w, y)
Usage: python dqi_pp_gvn.py <circuit.kmb> <out.npz> [seed=1]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

U = np.uint64
NIN = 840
NMASK = 49140


@njit(cache=True, inline='always')
def smix(x):
    x = x + U(0x9E3779B97F4A7C15)
    x = (x ^ (x >> U(30))) * U(0xBF58476D1CE4E5B9)
    x = (x ^ (x >> U(27))) * U(0x94D049BB133111EB)
    return x ^ (x >> U(31))


@njit(cache=True, inline='always')
def lab(k0, k1, k2, k3, salt):
    h = smix(salt ^ k0)
    h = smix(h ^ k1)
    h = smix(h ^ k2)
    h = smix(h ^ k3)
    g = smix(h ^ U(0xD6E8FEB86659FD93))
    g = smix(g ^ k3 ^ salt)
    return h, g


@njit(cache=True, inline='always')
def less(a0, a1, b0, b1):
    return a1 < b1 or (a1 == b1 and a0 < b0)


@njit(cache=True)
def prod(a0, a1, ac, b0, b1, bc, salt):
    """(A+ac)(B+bc) -> (r0, r1, rc, deg): deg 1 when affine (no atom)."""
    za = a0 == 0 and a1 == 0
    zb = b0 == 0 and b1 == 0
    if za and zb:
        return U(0), U(0), ac & bc, 1
    if za:
        if ac:
            return b0, b1, bc, 1
        return U(0), U(0), U(0), 1
    if zb:
        if bc:
            return a0, a1, ac, 1
        return U(0), U(0), U(0), 1
    if a0 == b0 and a1 == b1:
        if ((U(1) ^ ac ^ bc) & U(1)) != U(0):
            return a0, a1, ac & bc, 1
        return U(0), U(0), ac & bc, 1
    s0 = a0 ^ b0
    s1 = a1 ^ b1
    # exclude the largest of A, B, S; the other two (ordered) key the atom
    if less(a0, a1, b0, b1):
        if less(b0, b1, s0, s1):        # A < B < S: exclude S
            u0, u1, v0, v1, ex = a0, a1, b0, b1, 2
        elif less(a0, a1, s0, s1):      # A < S < B: exclude B
            u0, u1, v0, v1, ex = a0, a1, s0, s1, 1
        else:                           # S < A < B: exclude B
            u0, u1, v0, v1, ex = s0, s1, a0, a1, 1
    else:
        if less(a0, a1, s0, s1):        # B < A < S: exclude S
            u0, u1, v0, v1, ex = b0, b1, a0, a1, 2
        elif less(b0, b1, s0, s1):      # B < S < A: exclude A
            u0, u1, v0, v1, ex = b0, b1, s0, s1, 0
        else:                           # S < B < A: exclude A
            u0, u1, v0, v1, ex = s0, s1, b0, b1, 0
    l0, l1 = lab(u0, u1, v0, v1, salt)
    r0, r1 = l0, l1
    if ex == 1:                         # {A, A+B}: UV = A + AB
        r0 ^= a0
        r1 ^= a1
    elif ex == 0:                       # {B, A+B}: UV = B + AB
        r0 ^= b0
        r1 ^= b1
    if bc:
        r0 ^= a0
        r1 ^= a1
    if ac:
        r0 ^= b0
        r1 ^= b1
    return r0, r1, ac & bc, 0


@njit(cache=True)
def atom_label(a0, a1, b0, b1, salt):
    if (a0 == 0 and a1 == 0) or (b0 == 0 and b1 == 0) or (a0 == b0 and a1 == b1):
        return U(0), U(0)
    s0 = a0 ^ b0
    s1 = a1 ^ b1
    if less(a0, a1, b0, b1):
        if less(b0, b1, s0, s1):
            u0, u1, v0, v1 = a0, a1, b0, b1
        elif less(a0, a1, s0, s1):
            u0, u1, v0, v1 = a0, a1, s0, s1
        else:
            u0, u1, v0, v1 = s0, s1, a0, a1
    else:
        if less(a0, a1, s0, s1):
            u0, u1, v0, v1 = b0, b1, a0, a1
        elif less(b0, b1, s0, s1):
            u0, u1, v0, v1 = b0, b1, s0, s1
        else:
            u0, u1, v0, v1 = s0, s1, b0, b1
    return lab(u0, u1, v0, v1, salt)


@njit(cache=True)
def run(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, nq, nbits, seed, ntoff, zf=None, kz=None):
    salt = smix(U(seed) * U(0x2545F4914F6CDD1D))
    H0 = np.zeros(nq, np.uint64)
    H1 = np.zeros(nq, np.uint64)
    HC = np.zeros(nq, np.uint64)
    B0 = np.zeros(nbits, np.uint64)
    B1 = np.zeros(nbits, np.uint64)
    BC = np.zeros(nbits, np.uint64)
    for i in range(min(NIN, nq)):
        H0[i], H1[i] = lab(U(i), U(1), U(0), U(0), salt)
    for i in range(min(NMASK, nbits)):
        B0[i], B1[i] = lab(U(i), U(2), U(0), U(0), salt)
    # condition stack: effective condition form per level
    S0 = np.zeros(64, np.uint64)
    S1 = np.zeros(64, np.uint64)
    SC = np.zeros(64, np.uint64)
    SC[0] = U(1)
    sp = 0
    # outputs per Toffoli-class op
    t_op = np.zeros(ntoff, np.int64)
    t_kind = np.zeros(ntoff, np.int8)
    t_w = np.full(ntoff, -1, np.int32)
    t_cd = np.zeros(ntoff, np.int8)      # HELD: constant differs (add X)
    t_lvl = np.zeros(ntoff, np.int8)
    t_rel = np.zeros(ntoff, np.int16)    # DEG detail: bit flags and constants
    t_p0 = np.zeros(ntoff, np.uint64)    # CCX: product form (hash and constant)
    t_p1 = np.zeros(ntoff, np.uint64)
    t_pc = np.zeros(ntoff, np.uint8)
    t_c0 = np.zeros(ntoff, np.uint64)    # CCZ: canonical key of the cubic part (span of the three linear parts)
    t_c1 = np.zeros(ntoff, np.uint64)
    t_k0 = np.zeros(ntoff, np.uint64)    # effective condition form of the op
    t_k1 = np.zeros(ntoff, np.uint64)
    t_kc = np.zeros(ntoff, np.uint8)
    t_l0 = np.zeros(ntoff, np.uint64)    # canonical atom label of the product (affine-equivalence class)
    t_l1 = np.zeros(ntoff, np.uint64)
    t_z = np.zeros((ntoff, 3, 2), np.uint64)   # CCZ: forms (hash only) of the pair products ab, ac, bc
    e7 = np.zeros((7, 2), np.uint64)
    k = 0
    i3 = 0
    i2 = 0
    i1 = 0
    ib = 0
    ic = 0
    for oi in range(ops.shape[0]):
        op = ops[oi]
        if op == 0x81:
            c = cd[ic]
            ic += 1
            r0, r1, rc, dg = prod(S0[sp], S1[sp], SC[sp], B0[c], B1[c], BC[c], salt)
            sp += 1
            S0[sp] = r0
            S1[sp] = r1
            SC[sp] = rc
            continue
        base = op & 0x7F
        if base == 0x01:
            if sp > 0:
                sp -= 1
            continue
        nqk = (op & 0x30) >> 4
        a = 0
        b = 0
        t = 0
        if nqk == 3:
            a = q3a[i3]
            b = q3b[i3]
            t = q3c[i3]
            i3 += 1
        elif nqk == 2:
            a = q2a[i2]
            b = q2b[i2]
            i2 += 1
        elif nqk == 1:
            a = q1[i1]
            i1 += 1
        bit = 0
        if (op & 0x40) != 0 and base != 0x06 and base != 0x07:
            bit = bt[ib]
            ib += 1
        # effective condition K
        k0 = S0[sp]
        k1 = S1[sp]
        kc = SC[sp]
        lvl = sp
        if (op & 0x80) != 0:
            c = cd[ic]
            ic += 1
            k0, k1, kc, dg = prod(k0, k1, kc, B0[c], B1[c], BC[c], salt)
            lvl += 1
        kzero = k0 == 0 and k1 == 0 and kc == 0
        kone = k0 == 0 and k1 == 0 and kc == 1
        if base == 0x31 or base == 0x30 or base == 0x32:
            t_op[k] = oi
            t_lvl[k] = lvl
            t_k0[k] = k0
            t_k1[k] = k1
            t_kc[k] = np.uint8(kc)
        if kzero:
            if kz is not None:
                kz[oi] = np.uint8(1)                  # never executes: touches nothing
            if base == 0x31 or base == 0x30 or base == 0x32:
                k += 1
            continue
        if base == 0x11:                          # X
            H0[a] ^= k0
            H1[a] ^= k1
            HC[a] ^= kc
        elif base == 0x12 or base == 0x50:        # R / HMR: q = (1+K) q
            if kone:
                H0[a] = U(0)
                H1[a] = U(0)
                HC[a] = U(0)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc ^ U(1), H0[a], H1[a], HC[a], salt)
                H0[a] = r0
                H1[a] = r1
                HC[a] = rc
            if base == 0x50:                      # outcome: fresh atom (or old value where K = 0)
                l0, l1 = lab(U(oi), U(3), U(0), U(0), salt)
                if kone:
                    B0[bit] = l0
                    B1[bit] = l1
                    BC[bit] = U(0)
                else:
                    d0, d1, dc, dg = prod(k0, k1, kc, l0 ^ B0[bit], l1 ^ B1[bit], BC[bit], salt)
                    B0[bit] ^= d0
                    B1[bit] ^= d1
                    BC[bit] ^= dc
        elif base == 0x21:                        # CX a -> b
            r0, r1, rc, dg = prod(k0, k1, kc, H0[a], H1[a], HC[a], salt)
            H0[b] ^= r0
            H1[b] ^= r1
            HC[b] ^= rc
        elif base == 0x22:                        # SWAP
            r0, r1, rc, dg = prod(k0, k1, kc, H0[a] ^ H0[b], H1[a] ^ H1[b], HC[a] ^ HC[b], salt)
            H0[a] ^= r0
            H1[a] ^= r1
            HC[a] ^= rc
            H0[b] ^= r0
            H1[b] ^= r1
            HC[b] ^= rc
        elif base == 0x31:                        # CCX a b -> t
            p0, p1, pc, dg = prod(H0[a], H1[a], HC[a], H0[b], H1[b], HC[b], salt)
            t_l0[k], t_l1[k] = atom_label(H0[a], H1[a], H0[b], H1[b], salt)
            t_p0[k] = p0
            t_p1[k] = p1
            t_pc[k] = np.uint8(pc)
            kind = 0
            if dg == 1:
                kind = 1
                # relation detail: bit0 A=0, bit1 B=0, bit2 A=B; bit3 ac, bit4 bc
                rel = 0
                if H0[a] == 0 and H1[a] == 0:
                    rel |= 1
                if H0[b] == 0 and H1[b] == 0:
                    rel |= 2
                if H0[a] == H0[b] and H1[a] == H1[b]:
                    rel |= 4
                rel |= int(HC[a]) << 3
                rel |= int(HC[b]) << 4
                t_rel[k] = rel
            elif H0[t] == p0 and H1[t] == p1:
                kind = 3 if HC[t] == pc else 4
            else:
                for w in range(nq):
                    if w != t and H0[w] == p0 and H1[w] == p1:
                        kind = 2
                        t_w[k] = w
                        t_cd[k] = 0 if HC[w] == pc else 1
                        break
            t_kind[k] = kind
            k += 1
            r0, r1, rc, dg = prod(k0, k1, kc, p0, p1, pc, salt)
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x32:                        # CSWAP a: x=b, y=t
            x0 = H0[b] ^ H0[t]
            x1 = H1[b] ^ H1[t]
            xc = HC[b] ^ HC[t]
            p0, p1, pc, dg = prod(H0[a], H1[a], HC[a], x0, x1, xc, salt)
            t_l0[k], t_l1[k] = atom_label(H0[a], H1[a], x0, x1, salt)
            t_p0[k] = p0
            t_p1[k] = p1
            t_pc[k] = np.uint8(pc)
            kind = 0
            if dg == 1:
                kind = 1
                rel = 0
                if H0[a] == 0 and H1[a] == 0:
                    rel |= 1
                if x0 == 0 and x1 == 0:
                    rel |= 2
                if H0[a] == x0 and H1[a] == x1:
                    rel |= 4
                rel |= int(HC[a]) << 3
                rel |= int(xc) << 4
                t_rel[k] = rel
            else:
                for w in range(nq):
                    if w != b and w != t and H0[w] == p0 and H1[w] == p1:
                        kind = 2
                        t_w[k] = w
                        t_cd[k] = 0 if HC[w] == pc else 1
                        break
            t_kind[k] = kind
            k += 1
            r0, r1, rc, dg = prod(k0, k1, kc, p0, p1, pc, salt)
            H0[b] ^= r0
            H1[b] ^= r1
            HC[b] ^= rc
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x30:                        # CCZ: phase only; degenerate if A, B, C dependent
            zq0, zq1, zqc, dg = prod(H0[a], H1[a], HC[a], H0[b], H1[b], HC[b], salt)
            t_z[k, 0, 0] = zq0
            t_z[k, 0, 1] = zq1
            zq0, zq1, zqc, dg = prod(H0[a], H1[a], HC[a], H0[t], H1[t], HC[t], salt)
            t_z[k, 1, 0] = zq0
            t_z[k, 1, 1] = zq1
            zq0, zq1, zqc, dg = prod(H0[b], H1[b], HC[b], H0[t], H1[t], HC[t], salt)
            t_z[k, 2, 0] = zq0
            t_z[k, 2, 1] = zq1
            a0, a1 = H0[a], H1[a]
            b0, b1 = H0[b], H1[b]
            c0, c1 = H0[t], H1[t]
            rel = 0
            if a0 == 0 and a1 == 0:
                rel |= 1
            if b0 == 0 and b1 == 0:
                rel |= 2
            if c0 == 0 and c1 == 0:
                rel |= 4
            if a0 == b0 and a1 == b1:
                rel |= 8
            if a0 == c0 and a1 == c1:
                rel |= 16
            if b0 == c0 and b1 == c1:
                rel |= 32
            if (a0 ^ b0) == c0 and (a1 ^ b1) == c1:
                rel |= 64
            if rel != 0:
                t_kind[k] = 1
                t_rel[k] = rel | (int(HC[a]) << 7) | (int(HC[b]) << 8) | (int(HC[t]) << 9)
            else:
                e7[0, 0] = a0
                e7[0, 1] = a1
                e7[1, 0] = b0
                e7[1, 1] = b1
                e7[2, 0] = c0
                e7[2, 1] = c1
                e7[3, 0] = a0 ^ b0
                e7[3, 1] = a1 ^ b1
                e7[4, 0] = a0 ^ c0
                e7[4, 1] = a1 ^ c1
                e7[5, 0] = b0 ^ c0
                e7[5, 1] = b1 ^ c1
                e7[6, 0] = a0 ^ b0 ^ c0
                e7[6, 1] = a1 ^ b1 ^ c1
                for x in range(7):                    # insertion sort on (hi, lo)
                    y = x
                    while y > 0 and less(e7[y, 0], e7[y, 1], e7[y - 1, 0], e7[y - 1, 1]):
                        tmp0 = e7[y, 0]
                        tmp1 = e7[y, 1]
                        e7[y, 0] = e7[y - 1, 0]
                        e7[y, 1] = e7[y - 1, 1]
                        e7[y - 1, 0] = tmp0
                        e7[y - 1, 1] = tmp1
                        y -= 1
                h0 = salt ^ U(0x1234567)
                h1 = salt ^ U(0x7654321)
                for x in range(7):
                    h0 = smix(h0 ^ e7[x, 0])
                    h1 = smix(h1 ^ e7[x, 1] ^ h0)
                t_c0[k] = h0
                t_c1[k] = h1
            k += 1
        elif base == 0x40:                        # BIT_STORE0: b = (1+K) b
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(0)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit], salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x41:                        # BIT_STORE1: b = b + K (1 + b)
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(1)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit] ^ U(1), salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x42:                        # BIT_INVERT: b ^= K
            B0[bit] ^= k0
            B1[bit] ^= k1
            BC[bit] ^= kc
        # phase ops (Z, CZ, NEG, Z_POW), INFO_TAG, DEBUG: no value change
        if zf is not None and nqk > 0:
            zz = np.uint8(0)
            if H0[a] == 0 and H1[a] == 0 and HC[a] == 0:
                zz |= np.uint8(1)
            if nqk >= 2 and H0[b] == 0 and H1[b] == 0 and HC[b] == 0:
                zz |= np.uint8(2)
            if nqk == 3 and H0[t] == 0 and H1[t] == 0 and HC[t] == 0:
                zz |= np.uint8(4)
            # constant ONE after the op: bits 4..6
            if H0[a] == 0 and H1[a] == 0 and HC[a] == 1:
                zz |= np.uint8(16)
            if nqk >= 2 and H0[b] == 0 and H1[b] == 0 and HC[b] == 1:
                zz |= np.uint8(32)
            if nqk == 3 and H0[t] == 0 and H1[t] == 0 and HC[t] == 1:
                zz |= np.uint8(64)
            zf[oi] = zz
    return (t_op[:k], t_kind[:k], t_w[:k], t_cd[:k], t_lvl[:k], t_rel[:k], t_p0[:k], t_p1[:k], t_pc[:k], k,
            t_c0[:k], t_c1[:k], t_k0[:k], t_k1[:k], t_kc[:k], t_l0[:k], t_l1[:k], t_z[:k])


def load(path):
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    return c, ops, arr


def stages(c):
    """Top-level stage boundaries from '@enter:' / '@exit:' tags: list of (name, start op, end op)."""
    ops = np.asarray(c.ops)
    tp = np.flatnonzero(ops == 0x06)
    out = []
    names = ('mask_projection', 'initial', 'extract', 'synthesis', 'convert', 'scan', 'unconvert',
             'padded_unsynthesis', 'unextract')
    open_ = {}
    for k, p in enumerate(tp.tolist()):
        t = c.tags[k]
        if t.startswith('@enter:') and t[7:] in names:
            open_[t[7:]] = p
        elif t.startswith('@exit:') and t[6:] in open_:
            out.append((t[6:], open_.pop(t[6:]), p))
    return out


def main():
    t0 = time.time()
    path, out = sys.argv[1], sys.argv[2]
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    c, ops, arr = load(path)
    ntoff = int((((ops & 0x7C) == 0x30) & (ops != 0x81)).sum())
    (t_op, t_kind, t_w, t_cd, t_lvl, t_rel, t_p0, t_p1, t_pc, k, t_c0, t_c1, t_k0, t_k1,
     t_kc, t_l0, t_l1, t_z) = run(ops, *arr, c.num_qubits, c.num_bits, seed, ntoff)
    assert k == ntoff, (k, ntoff)
    base = (ops[t_op] & 0x7F).astype(np.int64)
    wt = np.power(0.5, t_lvl.astype(np.float64))
    print('ops %d, Toffoli-class %d, %.1fs' % (len(ops), k, time.time() - t0))
    st = stages(c)
    names = {0x31: 'CCX', 0x30: 'CCZ', 0x32: 'CSWAP'}
    kinds = {1: 'DEG', 2: 'HELD', 3: 'UNC0', 4: 'UNC1'}
    for code in (0x31, 0x30, 0x32):
        for kd in (1, 2, 3, 4):
            m = (base == code) & (t_kind == kd)
            if not m.any():
                continue
            per = []
            for nm, a, b in st:
                mm = m & (t_op >= a) & (t_op < b)
                if mm.any():
                    per.append('%s %d/%.1f' % (nm[:8], int(mm.sum()), float(wt[mm].sum())))
            print('%-5s %-4s sites %7d  expected %9.1f  free %6d | %s' % (
                names[code], kinds[kd], int(m.sum()), float(wt[m].sum()), int((m & (t_lvl == 0)).sum()),
                '; '.join(per)))
    np.savez(out, t_op=t_op, t_kind=t_kind, t_w=t_w, t_cd=t_cd, t_lvl=t_lvl, t_rel=t_rel, t_p0=t_p0, t_p1=t_p1,
             t_pc=t_pc, t_c0=t_c0, t_c1=t_c1, t_k0=t_k0, t_k1=t_k1, t_kc=t_kc, t_l0=t_l0, t_l1=t_l1, t_z=t_z)
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
