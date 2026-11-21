#!/usr/bin/env python3
"""r55_dqi_gpeak_zero.py -- R55-DQI-GPEAK (designer): EXACT idle-zero intervals with PUSH blocks followed as blocks.

P1's tools/dqi_pp_zero.py (value numbering of tools/dqi_pp_gvn.py; both read only, their prod / lab are imported)
models an op under a classical condition K as  wire ^= K x form  with K x form ONE product atom. That is sound but
it loses K a + K b = K (a + b) and K K = K: a wire that a PUSH block writes and erases again (by gates or by HMR)
keeps a nonzero formal value until some LATER unconditional HMR, and counts as live in every "tight" census.

Here the outermost PUSH block is run in a SHADOW state "as if the block executes" (the same value numbering, with
only the inner conditions as products) and merged at its POP:
    value after the block = old + K (shadow + old)      for every wire and classical bit the block wrote,
so a wire whose shadow form equals its old form is unchanged (in particular: zero stays zero), exactly.
Inside the block a wire is zero iff its shadow form AND its pre-block form are both zero (the two branches).
Sound for the same reason P1's is: equal 128-bit forms are equal Boolean functions; the merge formula is the
definition of a conditional block. Nothing sampled. Ops whose own _IF condition stands outside any block keep P1's
treatment. Output format = P1's zero.npz (wire, start, end, structural): a drop-in for every census / pass.
Check (optional 3rd argument): my empirical file (tools/r55_dqi_gpeak_emp.py): every interval here must start at an
event after which the wire was zero on every sample; violations are printed and must be 0.
Usage: python r55_dqi_gpeak_zero.py <circuit.kmb> <out.npz> [empirical.npz] [P1 zero.npz to compare]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_gvn as g  # noqa: E402
import dqi_pp_prove as P  # noqa: E402
from dqi_pp_gvn import prod, lab, smix, U, NIN, NMASK  # noqa: E402

SH = 27


@njit(cache=True)
def run_blocks(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, nq, nbits, seed, zf, kz, stats):
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
    # inner condition stack (relative to the outermost block): S[0] = 1 outside, S[1] = 1 inside the block
    S0 = np.zeros(66, np.uint64)
    S1 = np.zeros(66, np.uint64)
    SC = np.zeros(66, np.uint64)
    SC[0] = U(1)
    sp = 0
    # the outermost block: its condition form and the saved pre-block forms of what it wrote
    K0 = U(0)
    K1 = U(0)
    KC = U(1)
    serial = 0
    skip = 0                                  # > 0: inside a block that never executes (depth counter)
    wst = np.zeros(nq, np.int64)
    o0 = np.zeros(nq, np.uint64)
    o1 = np.zeros(nq, np.uint64)
    oc = np.zeros(nq, np.uint64)
    wl = np.zeros(nq, np.int64)
    nw = 0
    bst = np.zeros(nbits, np.int64)
    p0 = np.zeros(nbits, np.uint64)
    p1 = np.zeros(nbits, np.uint64)
    pc = np.zeros(nbits, np.uint64)
    bl = np.zeros(nbits, np.int64)
    nb = 0
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
            if skip > 0:
                skip += 1
                continue
            if sp == 0:
                # outermost block: condition form at PUSH time
                if B0[c] == 0 and B1[c] == 0 and BC[c] == 0:
                    skip = 1
                    continue
                K0, K1, KC = B0[c], B1[c], BC[c]
                serial += 1
                nw = 0
                nb = 0
                sp = 1
                S0[1], S1[1], SC[1] = U(0), U(0), U(1)
                stats[0] += 1
            else:
                r0, r1, rc, dg = prod(S0[sp], S1[sp], SC[sp], B0[c], B1[c], BC[c], salt)
                sp += 1
                S0[sp], S1[sp], SC[sp] = r0, r1, rc
            continue
        base = op & 0x7F
        if base == 0x01:
            if skip > 0:
                skip -= 1
                continue
            if sp > 1:
                sp -= 1
            elif sp == 1:
                sp = 0
                kone = K0 == 0 and K1 == 0 and KC == 1
                for j in range(nw):
                    w = wl[j]
                    if H0[w] == o0[w] and H1[w] == o1[w] and HC[w] == oc[w]:
                        stats[1] += 1
                        continue
                    stats[2] += 1
                    if not kone:
                        r0, r1, rc, dg = prod(K0, K1, KC, H0[w] ^ o0[w], H1[w] ^ o1[w], HC[w] ^ oc[w], salt)
                        H0[w] = o0[w] ^ r0
                        H1[w] = o1[w] ^ r1
                        HC[w] = oc[w] ^ rc
                for j in range(nb):
                    x = bl[j]
                    if B0[x] == p0[x] and B1[x] == p1[x] and BC[x] == pc[x]:
                        continue
                    if not kone:
                        r0, r1, rc, dg = prod(K0, K1, KC, B0[x] ^ p0[x], B1[x] ^ p1[x], BC[x] ^ pc[x], salt)
                        B0[x] = p0[x] ^ r0
                        B1[x] = p1[x] ^ r1
                        BC[x] = pc[x] ^ rc
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
        k0 = S0[sp]
        k1 = S1[sp]
        kc = SC[sp]
        if (op & 0x80) != 0:
            c = cd[ic]
            ic += 1
            if skip == 0:
                k0, k1, kc, dg = prod(k0, k1, kc, B0[c], B1[c], BC[c], salt)
        if skip > 0:
            kz[oi] = np.uint8(1)
            continue
        kzero = k0 == 0 and k1 == 0 and kc == 0
        kone = k0 == 0 and k1 == 0 and kc == 1
        if kzero:
            kz[oi] = np.uint8(1)
            continue
        inblk = sp >= 1
        # save pre-block forms of what this op writes
        if inblk:
            wa = -1
            wb = -1
            if base == 0x11 or base == 0x12 or base == 0x50:
                wa = a
            elif base == 0x21:
                wa = b
            elif base == 0x22:
                wa = a
                wb = b
            elif base == 0x31:
                wa = t
            elif base == 0x32:
                wa = b
                wb = t
            if wa >= 0 and wst[wa] != serial:
                wst[wa] = serial
                o0[wa], o1[wa], oc[wa] = H0[wa], H1[wa], HC[wa]
                wl[nw] = wa
                nw += 1
            if wb >= 0 and wst[wb] != serial:
                wst[wb] = serial
                o0[wb], o1[wb], oc[wb] = H0[wb], H1[wb], HC[wb]
                wl[nw] = wb
                nw += 1
            if base == 0x50 or base == 0x40 or base == 0x41 or base == 0x42:
                if bst[bit] != serial:
                    bst[bit] = serial
                    p0[bit], p1[bit], pc[bit] = B0[bit], B1[bit], BC[bit]
                    bl[nb] = bit
                    nb += 1
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
            if base == 0x50:
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
            q0, q1_, qc, dg = prod(H0[a], H1[a], HC[a], H0[b], H1[b], HC[b], salt)
            r0, r1, rc, dg = prod(k0, k1, kc, q0, q1_, qc, salt)
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x32:                        # CSWAP a: x = b, y = t
            x0 = H0[b] ^ H0[t]
            x1 = H1[b] ^ H1[t]
            xc = HC[b] ^ HC[t]
            q0, q1_, qc, dg = prod(H0[a], H1[a], HC[a], x0, x1, xc, salt)
            r0, r1, rc, dg = prod(k0, k1, kc, q0, q1_, qc, salt)
            H0[b] ^= r0
            H1[b] ^= r1
            HC[b] ^= rc
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x40:                        # BIT_STORE0
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(0)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit], salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x41:                        # BIT_STORE1
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(1)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit] ^ U(1), salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x42:                        # BIT_INVERT
            B0[bit] ^= k0
            B1[bit] ^= k1
            BC[bit] ^= kc
        if nqk > 0:
            zz = np.uint8(0)
            za = H0[a] == 0 and H1[a] == 0 and HC[a] == 0
            if inblk and za and wst[a] == serial:
                za = o0[a] == 0 and o1[a] == 0 and oc[a] == 0
            if za:
                zz |= np.uint8(1)
            if nqk >= 2:
                zb = H0[b] == 0 and H1[b] == 0 and HC[b] == 0
                if inblk and zb and wst[b] == serial:
                    zb = o0[b] == 0 and o1[b] == 0 and oc[b] == 0
                if zb:
                    zz |= np.uint8(2)
            if nqk == 3:
                zt = H0[t] == 0 and H1[t] == 0 and HC[t] == 0
                if inblk and zt and wst[t] == serial:
                    zt = o0[t] == 0 and o1[t] == 0 and oc[t] == 0
                if zt:
                    zz |= np.uint8(4)
            zf[oi] = zz
    return 0


@njit(cache=True)
def _active(fserial, sp, st):
    for d in range(1, sp + 1):
        if fserial[d] == st:
            return True
    return False


@njit(cache=True)
def run_frames(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, nq, nbits, seed, zf, kz, stats):
    """As run_blocks, with EVERY PUSH level its own shadow frame (a nested block is merged into its parent's
    shadow at its POP, the outermost into the real state). An op's own _IF bit keeps P1's product treatment."""
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
    MAXD = 66
    FK0 = np.zeros(MAXD, np.uint64)
    FK1 = np.zeros(MAXD, np.uint64)
    FKC = np.zeros(MAXD, np.uint64)
    fserial = np.zeros(MAXD, np.int64)
    fjs = np.zeros(MAXD, np.int64)
    fjb = np.zeros(MAXD, np.int64)
    sp = 0
    serial = 0
    skip = 0
    # journal of saved pre-frame forms (a stack; the innermost frame's entries are on top)
    JN = 4 * nq + 1024
    jw = np.zeros(JN, np.int64)
    j0 = np.zeros(JN, np.uint64)
    j1 = np.zeros(JN, np.uint64)
    jc = np.zeros(JN, np.uint64)
    jp = np.zeros(JN, np.int64)          # the wire's stamp before this frame saved it
    jz = np.zeros(JN, np.uint8)          # the wire's zok flag before this frame saved it
    nj = 0
    JB = 4 * nbits + 1024
    kb = np.zeros(JB, np.int64)
    k0b = np.zeros(JB, np.uint64)
    k1b = np.zeros(JB, np.uint64)
    kcb = np.zeros(JB, np.uint64)
    kpb = np.zeros(JB, np.int64)
    nbj = 0
    wst = np.zeros(nq, np.int64)         # serial of the innermost frame that saved the wire
    zok = np.ones(nq, np.uint8)          # every saved pre-frame form of the wire (in the active frames) is zero
    bst = np.zeros(nbits, np.int64)
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
            if skip > 0:
                skip += 1
                continue
            if B0[c] == 0 and B1[c] == 0 and BC[c] == 0:
                skip = 1
                continue
            sp += 1
            serial += 1
            FK0[sp], FK1[sp], FKC[sp] = B0[c], B1[c], BC[c]
            fserial[sp] = serial
            fjs[sp] = nj
            fjb[sp] = nbj
            stats[0] += 1
            continue
        base = op & 0x7F
        if base == 0x01:
            if skip > 0:
                skip -= 1
                continue
            if sp == 0:
                continue
            K0, K1, KC = FK0[sp], FK1[sp], FKC[sp]
            kone = K0 == 0 and K1 == 0 and KC == 1
            par = fserial[sp - 1] if sp > 1 else -1
            wi = fjs[sp]
            for j in range(fjs[sp], nj):
                w = jw[j]
                if H0[w] == j0[j] and H1[w] == j1[j] and HC[w] == jc[j]:
                    stats[1] += 1
                else:
                    stats[2] += 1
                    if not kone:
                        r0, r1, rc, dg = prod(K0, K1, KC, H0[w] ^ j0[j], H1[w] ^ j1[j], HC[w] ^ jc[j], salt)
                        H0[w] = j0[j] ^ r0
                        H1[w] = j1[j] ^ r1
                        HC[w] = jc[j] ^ rc
                if sp > 1:
                    if jp[j] == par:                 # the parent frame had saved it already
                        wst[w] = par
                        zok[w] = jz[j]
                    else:                            # the parent now holds this entry
                        pz = np.uint8(1)
                        if _active(fserial, sp - 1, jp[j]):
                            pz = jz[j]
                        oldz = j0[j] == 0 and j1[j] == 0 and jc[j] == 0
                        o0_, o1_, oc_, op_, oz_ = j0[j], j1[j], jc[j], jp[j], jz[j]
                        jw[wi] = w
                        j0[wi] = o0_
                        j1[wi] = o1_
                        jc[wi] = oc_
                        jp[wi] = op_
                        jz[wi] = oz_
                        wi += 1
                        wst[w] = par
                        zok[w] = np.uint8(1) if (oldz and pz == 1) else np.uint8(0)
                else:
                    wst[w] = 0
                    zok[w] = np.uint8(1)
            nj = wi
            wi = fjb[sp]
            for j in range(fjb[sp], nbj):
                x = kb[j]
                if not (B0[x] == k0b[j] and B1[x] == k1b[j] and BC[x] == kcb[j]):
                    if not kone:
                        r0, r1, rc, dg = prod(K0, K1, KC, B0[x] ^ k0b[j], B1[x] ^ k1b[j], BC[x] ^ kcb[j], salt)
                        B0[x] = k0b[j] ^ r0
                        B1[x] = k1b[j] ^ r1
                        BC[x] = kcb[j] ^ rc
                if sp > 1:
                    if kpb[j] == par:
                        bst[x] = par
                    else:
                        o0_, o1_, oc_, op_ = k0b[j], k1b[j], kcb[j], kpb[j]
                        kb[wi] = x
                        k0b[wi] = o0_
                        k1b[wi] = o1_
                        kcb[wi] = oc_
                        kpb[wi] = op_
                        wi += 1
                        bst[x] = par
                else:
                    bst[x] = 0
            nbj = wi
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
        k0 = U(0)
        k1 = U(0)
        kc = U(1)
        if (op & 0x80) != 0:
            c = cd[ic]
            ic += 1
            if skip == 0:
                k0, k1, kc = B0[c], B1[c], BC[c]
        if skip > 0:
            kz[oi] = np.uint8(1)
            continue
        kzero = k0 == 0 and k1 == 0 and kc == 0
        kone = k0 == 0 and k1 == 0 and kc == 1
        if kzero:
            kz[oi] = np.uint8(1)
            continue
        if sp >= 1:
            cur = fserial[sp]
            wa = -1
            wb = -1
            if base == 0x11 or base == 0x12 or base == 0x50:
                wa = a
            elif base == 0x21:
                wa = b
            elif base == 0x22:
                wa = a
                wb = b
            elif base == 0x31:
                wa = t
            elif base == 0x32:
                wa = b
                wb = t
            for rep in range(2):
                w = wa if rep == 0 else wb
                if w >= 0 and wst[w] != cur:
                    pz = np.uint8(1)
                    if _active(fserial, sp - 1, wst[w]):
                        pz = zok[w]
                    jw[nj] = w
                    j0[nj], j1[nj], jc[nj] = H0[w], H1[w], HC[w]
                    jp[nj] = wst[w]
                    jz[nj] = zok[w]
                    nj += 1
                    oldz = H0[w] == 0 and H1[w] == 0 and HC[w] == 0
                    wst[w] = cur
                    zok[w] = np.uint8(1) if (oldz and pz == 1) else np.uint8(0)
            if base == 0x50 or base == 0x40 or base == 0x41 or base == 0x42:
                if bst[bit] != cur:
                    kb[nbj] = bit
                    k0b[nbj], k1b[nbj], kcb[nbj] = B0[bit], B1[bit], BC[bit]
                    kpb[nbj] = bst[bit]
                    nbj += 1
                    bst[bit] = cur
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
            if base == 0x50:
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
            q0, q1_, qc, dg = prod(H0[a], H1[a], HC[a], H0[b], H1[b], HC[b], salt)
            r0, r1, rc, dg = prod(k0, k1, kc, q0, q1_, qc, salt)
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x32:                        # CSWAP a: x = b, y = t
            x0 = H0[b] ^ H0[t]
            x1 = H1[b] ^ H1[t]
            xc = HC[b] ^ HC[t]
            q0, q1_, qc, dg = prod(H0[a], H1[a], HC[a], x0, x1, xc, salt)
            r0, r1, rc, dg = prod(k0, k1, kc, q0, q1_, qc, salt)
            H0[b] ^= r0
            H1[b] ^= r1
            HC[b] ^= rc
            H0[t] ^= r0
            H1[t] ^= r1
            HC[t] ^= rc
        elif base == 0x40:                        # BIT_STORE0
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(0)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit], salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x41:                        # BIT_STORE1
            if kone:
                B0[bit] = U(0)
                B1[bit] = U(0)
                BC[bit] = U(1)
            else:
                r0, r1, rc, dg = prod(k0, k1, kc, B0[bit], B1[bit], BC[bit] ^ U(1), salt)
                B0[bit] ^= r0
                B1[bit] ^= r1
                BC[bit] ^= rc
        elif base == 0x42:                        # BIT_INVERT
            B0[bit] ^= k0
            B1[bit] ^= k1
            BC[bit] ^= kc
        if nqk > 0:
            zz = np.uint8(0)
            za = H0[a] == 0 and H1[a] == 0 and HC[a] == 0
            if za and sp >= 1 and _active(fserial, sp, wst[a]):
                za = zok[a] == 1
            if za:
                zz |= np.uint8(1)
            if nqk >= 2:
                zb = H0[b] == 0 and H1[b] == 0 and HC[b] == 0
                if zb and sp >= 1 and _active(fserial, sp, wst[b]):
                    zb = zok[b] == 1
                if zb:
                    zz |= np.uint8(2)
            if nqk == 3:
                zt = H0[t] == 0 and H1[t] == 0 and HC[t] == 0
                if zt and sp >= 1 and _active(fserial, sp, wst[t]):
                    zt = zok[t] == 1
                if zt:
                    zz |= np.uint8(4)
            zf[oi] = zz
    return 0


def intervals(c, ops, zf, kz):
    """As tools/dqi_pp_zero.py main (P1), from the zero flags."""
    _, W0, W1, W2, BT, CB, DEP = P.decode(c)
    n = len(ops)
    base = ops & 0x7F
    Wm = np.stack([W0, W1, W2], axis=1)
    ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
    ev_w = Wm.reshape(-1).astype(np.int64)
    ev_slot = np.tile(np.arange(3), n)
    kk = (ev_w >= 0) & (kz[ev_op] == 0)
    ev_op, ev_w, ev_slot = ev_op[kk], ev_w[kk], ev_slot[kk]
    zero_after = ((zf[ev_op] >> ev_slot.astype(np.uint8)) & 1).astype(bool)
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, zero_after = ev_op[o], ev_w[o], zero_after[o]
    ne = len(ev_op)
    nxt = np.full(ne, n, np.int64)
    same = ev_w[1:] == ev_w[:-1]
    nxt[:-1][same] = ev_op[1:][same]
    sel = np.flatnonzero(zero_after)
    wire, start, end = ev_w[sel], ev_op[sel], nxt[sel]
    structural = np.isin(base[start], (0x50, 0x12)) & (CB[start] < 0) & (DEP[start] == 0)
    first = np.ones(ne, bool)
    first[1:] = ~same
    fw, fo = ev_w[first], ev_op[first]
    m = fw >= 840
    touched = np.zeros(c.num_qubits, bool)
    touched[fw] = True
    never = np.flatnonzero(~touched & (np.arange(c.num_qubits) >= 840))
    wire = np.concatenate([wire, fw[m], never])
    start = np.concatenate([start, np.full(int(m.sum()), -1), np.full(len(never), -1)])
    end = np.concatenate([end, fo[m], np.full(len(never), n)])
    structural = np.concatenate([structural, np.ones(int(m.sum()) + len(never), bool)])
    keep = end > start + 1
    return wire[keep], start[keep], end[keep], structural[keep]


def live_of(nq, n, start, end):
    d = np.zeros(n + 2, np.int32)
    s = np.clip(start, 0, n)
    e = np.clip(end, 0, n)
    ok = e > s
    np.add.at(d, s[ok], 1)
    np.add.at(d, e[ok], -1)
    return (nq - np.cumsum(d[:n])).astype(np.int32)


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    empp = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != '-' else None
    refp = sys.argv[4] if len(sys.argv) > 4 else None
    gc, ops, arr = g.load(sys.argv[1])
    n = len(ops)
    zf = np.zeros(n, np.uint8)
    kz = np.zeros(n, np.uint8)
    stats = np.zeros(4, np.int64)
    if os.environ.get('R55_ZERO_RULE', 'frames') == 'blocks':     # the first form: outermost blocks only
        run_blocks(ops, *arr, c.num_qubits, c.num_bits, 1, zf, kz, stats)
    else:                                                        # every PUSH level its own frame
        run_frames(ops, *arr, c.num_qubits, c.num_bits, 1, zf, kz, stats)
    print('PUSH blocks that can execute: %d; wires written by a block and UNCHANGED at its POP: %d; '
          'changed: %d; ops that never execute: %d (%.1fs)' % (stats[0], stats[1], stats[2], int(kz.sum()),
                                                               time.time() - t0), flush=True)
    wire, start, end, structural = intervals(c, ops, zf, kz)
    print('exact idle-zero intervals %d (structural %d), wire-ops %d' % (
        len(wire), int(structural.sum()), int((end - start).sum())))
    np.savez(out, wire=wire, start=start, end=end, structural=structural)
    L = live_of(c.num_qubits, n, start, end)
    tof = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    levels = (864, 877, 889, 896, 899, 900)
    Lr = None
    if refp:
        R = np.load(refp)
        Lr = live_of(c.num_qubits, n, R['start'].astype(np.int64), R['end'].astype(np.int64))
        # superset check: every P1 interval start must be zero-after here too
        key_new = set(((wire.astype(np.int64) << SH) | np.clip(start, 0, n)).tolist())
        rs = np.clip(R['start'].astype(np.int64), 0, n)
        rk = ((R['wire'].astype(np.int64) << SH) | rs)
        miss = sum(1 for x in rk.tolist() if x not in key_new)
        print('P1 intervals %d; of them not starting an interval here: %d (same events are compared by wire and start)' % (
            len(rk), miss))
        print('wire-ops free: P1 %d, here %d' % (int((np.clip(R['end'], 0, n).astype(np.int64) - rs).sum()),
                                                 int((np.clip(end, 0, n) - np.clip(start, 0, n)).sum())))
    print('tight peak of the file: %s here %d' % (('P1 %d,' % int(Lr.max())) if Lr is not None else '', int(L.max())))
    print('%-20s %-10s %6s  %s' % ('stage', 'zero file', 'peak', 'Toffoli-class ops at live > ' + ' / '.join(map(str, levels))))
    for nm, a, b in g.stages(c):
        T = tof[a:b]
        for lab_, LL in ((('P1', Lr),) if Lr is not None else ()) + (('blocks', L),):
            X = LL[a:b]
            print('%-20s %-10s %6d  %s' % (nm, lab_, int(X.max()), ' / '.join('%d' % int(T[X > lv].sum()) for lv in levels)))
    if empp:
        E = np.load(empp)
        ek = set(((E['wire'].astype(np.int64) << SH) | np.clip(E['start'].astype(np.int64), 0, n)).tolist())
        sel = start >= 0
        kk = ((wire[sel].astype(np.int64) << SH) | start[sel]).tolist()
        bad = [x for x in kk if x not in ek]
        print('SOUNDNESS against %d samples: intervals checked %d, starting at an event after which the wire was NOT '
              'zero on every sample: %d' % (int(E['samples']), len(kk), len(bad)))
        for x in bad[:20]:
            print('   wire %d op %d' % (x >> SH, x & ((1 << SH) - 1)))
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
