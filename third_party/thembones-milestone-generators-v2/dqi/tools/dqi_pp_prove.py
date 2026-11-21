#!/usr/bin/env python3
"""dqi_pp_prove.py -- workstream dqi_pp (P1): EXACT local prover for candidate rewrites (no sampling).

For a claim F(values at op tau) == 0 (Z: product never 1; EA / EB: one control implies the other; U: target equals
the product; H: wire w equals the product), walk BACKWARD from tau and keep the set S of wires and classical bits
whose values at that point still influence F (liveness: a write into S adds the operands it reads; an unconditioned
reset (R / HMR / BIT_STORE) removes its target; an HMR outcome in S becomes a fresh free variable). Choose the
earliest cut (outside every PUSH_CONDITION block) at which |S| + fresh <= K. Then enumerate ALL 2^n assignments of
the cut variables (bit-parallel) and simulate the window exactly (values, blocks, classical bits). If F is zero on
every assignment the claim is PROVED for every reachable state, since the reachable states at the cut are a subset
of all assignments. A near-identity is false on some reachable state, hence on some assignment: it cannot pass.

Usage: python dqi_pp_prove.py <circuit.kmb> <fp.npz> <out.npz> [K=16] [L=400000] [classes=Z,EA,EB,U,H]
Writes per candidate: op, class, w, result (1 proved, 0 not proved, -1 cone too large), nvars.
"""
import os
import sys
import time

import numpy as np
from numba import njit, prange

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)
CL = {'Z': 1, 'EA': 2, 'EB': 3, 'U': 4, 'H': 5, 'W0': 6, 'PAIR': 7}


def decode(c):
    """Per-op arrays: W0, W1, W2 (wires, -1), BT (bit written, -1), CB (own condition bit, -1), DEP (block depth)."""
    ops = np.asarray(c.ops)
    n = len(ops)
    nq = (ops & 0x30) >> 4
    W = np.full((3, n), -1, dtype=np.int32)
    for k, arrs in ((3, c.q3), (2, c.q2), (1, (c.q1,))):
        sel = np.flatnonzero(nq == k)
        for j in range(k):
            W[j, sel] = np.asarray(arrs[j], dtype=np.int32)
    base = ops & 0x7F
    hb = ((ops & 0x40) != 0) & (base != 0x06) & (base != 0x07)
    BT = np.full(n, -1, dtype=np.int32)
    BT[hb] = np.asarray(c.bit_t, dtype=np.int32)
    hc = (ops & 0x80) != 0
    CB = np.full(n, -1, dtype=np.int32)
    CB[hc] = np.asarray(c.cond, dtype=np.int32)
    push = ops == 0x81
    pop = ops == 0x01
    DEP = (np.cumsum(push.astype(np.int64) - pop.astype(np.int64)) - push).astype(np.int8)
    return np.ascontiguousarray(ops), W[0].copy(), W[1].copy(), W[2].copy(), BT, CB, DEP


@njit(cache=True)
def claim_wires(cls, a, b, t, w):
    """Wires the claim reads at tau."""
    if cls == 4:
        return a, b, t
    if cls == 5:
        return a, b, w
    return a, b, t                                   # Z on CCZ / CSWAP uses three; CCX Z uses a, b (t harmless)


@njit(cache=True)
def walk(ops, W0, W1, W2, BT, CB, DEP, tau, cls, w, K, L, nq, nbits, inS, inB, mode, cut_req,
         outw, outb, fresh_ops, uw, ub):
    """Backward liveness walk from tau. mode 0: find the cut (returns it); mode 1: walk to cut_req and fill the S
    sets (outw / outb: membership at the cut), fresh_ops (HMR ops whose outcome is a variable) and the unions uw / ub
    of every wire / bit that was ever in S (the only ones simulated); returns the fresh count."""
    inS[:] = False
    inB[:] = False
    a = np.int64(W0[tau])
    b = np.int64(W1[tau])
    t = np.int64(W2[tau])
    w = np.int64(w)
    base_t = ops[tau] & 0x7F
    ns = 0
    nb = 0
    # claim operands
    fi = np.int64(-1)
    if cls == 7:                                     # product at op w (= fi) equals product at op tau
        fi = w
        for x in (a, b):
            if not inS[x]:
                inS[x] = True
                ns += 1
    elif cls == 6:                                   # wire w is zero before op tau
        inS[w] = True
        ns = 1
    elif base_t == 0x31:
        if cls == 1 or cls == 2 or cls == 3:
            for x in (a, b):
                if not inS[x]:
                    inS[x] = True
                    ns += 1
        else:
            x3 = t if cls == 4 else w
            for x in (a, b, x3):
                if not inS[x]:
                    inS[x] = True
                    ns += 1
    else:
        for x in (a, b, t):
            if not inS[x]:
                inS[x] = True
                ns += 1
    if CB[tau] >= 0 and ops[tau] != 0x81 and cls != 6:
        if not inB[CB[tau]]:
            inB[CB[tau]] = True
            nb += 1
    if mode == 1:
        uw[a] = True
        uw[b] = True
        if t >= 0:
            uw[t] = True
        if w >= 0 and cls != 7:
            uw[w] = True
        if CB[tau] >= 0:
            ub[CB[tau]] = True
    # enclosing blocks of tau are relevant: their PUSH bits must be tracked
    rel = np.zeros(64, np.bool_)
    rd = DEP[tau] if cls != 6 else 0
    for lv in range(rd):
        rel[lv] = True
    fresh = 0
    best = -1
    stop = max(0, tau - L)
    j = tau - 1
    nfo = 0
    while j >= stop:
        if mode == 1 and j < cut_req:
            break
        op = ops[j]
        base = op & 0x7F
        d = DEP[j]
        if mode == 1:                                # union: operands needed just after op j
            for x in (W0[j], W1[j], W2[j]):
                if x >= 0 and inS[x]:
                    uw[x] = True
            if BT[j] >= 0 and inB[BT[j]]:
                ub[BT[j]] = True
        if op == 0x81:                               # PUSH (backward: closes the block opened at its POP)
            lv = d                                   # this block's level
            if rel[lv]:
                bb = CB[j]
                if not inB[bb]:
                    inB[bb] = True
                    nb += 1
                if lv > 0:
                    rel[lv - 1] = True
            rel[lv] = False
        elif base == 0x01:                           # POP: block of level d - 1 ... nothing to do
            pass
        else:
            touched = False
            ccb = CB[j]
            # quantum targets
            if base == 0x11 or base == 0x12 or base == 0x50:          # X, R, HMR on W0
                x = W0[j]
                if inS[x]:
                    touched = True
                    if (base == 0x12 or base == 0x50) and ccb < 0 and d == 0:
                        inS[x] = False
                        ns -= 1
            elif base == 0x21:                                         # CX W0 -> W1
                x = W1[j]
                if inS[x]:
                    touched = True
                    y = W0[j]
                    if not inS[y]:
                        inS[y] = True
                        ns += 1
            elif base == 0x22:                                         # SWAP
                x = W0[j]
                y = W1[j]
                if inS[x] or inS[y]:
                    touched = True
                    if not (inS[x] and inS[y]):
                        cond_any = ccb >= 0 or d > 0
                        if cond_any:
                            if not inS[x]:
                                inS[x] = True
                                ns += 1
                            if not inS[y]:
                                inS[y] = True
                                ns += 1
                        else:
                            tmp = inS[x]
                            inS[x] = inS[y]
                            inS[y] = tmp
            elif base == 0x31:                                         # CCX
                x = W2[j]
                if inS[x]:
                    touched = True
                    for y in (W0[j], W1[j]):
                        if not inS[y]:
                            inS[y] = True
                            ns += 1
            elif base == 0x32:                                         # CSWAP
                x = W1[j]
                y = W2[j]
                if inS[x] or inS[y]:
                    touched = True
                    for z in (W0[j], x, y):
                        if not inS[z]:
                            inS[z] = True
                            ns += 1
            # classical targets
            bw = BT[j]
            if bw >= 0 and inB[bw]:
                if base == 0x50:
                    touched = True
                    if ccb < 0 and d == 0:
                        inB[bw] = False
                        nb -= 1
                    fresh += 1
                    if mode == 1:
                        fresh_ops[nfo] = j
                        nfo += 1
                elif base == 0x40 or base == 0x41:
                    touched = True
                    if ccb < 0 and d == 0:
                        inB[bw] = False
                        nb -= 1
                elif base == 0x42:
                    touched = True
            if touched:
                if ccb >= 0:
                    if not inB[ccb]:
                        inB[ccb] = True
                        nb += 1
                if d > 0:
                    rel[d - 1] = True
        if cls == 7 and j == fi:                     # the claim reads the controls of the first CCX here
            for x in (W0[fi], W1[fi]):
                if not inS[x]:
                    inS[x] = True
                    ns += 1
        if mode == 1:                                # union: operands needed just before op j
            for x in (W0[j], W1[j], W2[j]):
                if x >= 0 and inS[x]:
                    uw[x] = True
            if BT[j] >= 0 and inB[BT[j]]:
                ub[BT[j]] = True
            if CB[j] >= 0 and inB[CB[j]]:
                ub[CB[j]] = True
        # a cut is allowed before op j when op j is outside every block
        if mode == 0:
            if DEP[j] == 0 and op != 0x81 and base != 0x01 and ns + nb + fresh <= K and (cls != 7 or j <= fi):
                best = j
            if ns + nb + fresh > K + 14:
                break
        j -= 1
    if mode == 0:
        return best
    for x in range(nq):
        outw[x] = inS[x]
    for x in range(nbits):
        outb[x] = inB[x]
    return nfo


@njit(cache=True)
def pattern(var, word):
    if var < 6:
        if var == 0:
            return np.uint64(0xAAAAAAAAAAAAAAAA)
        if var == 1:
            return np.uint64(0xCCCCCCCCCCCCCCCC)
        if var == 2:
            return np.uint64(0xF0F0F0F0F0F0F0F0)
        if var == 3:
            return np.uint64(0xFF00FF00FF00FF00)
        if var == 4:
            return np.uint64(0xFFFF0000FFFF0000)
        return np.uint64(0xFFFFFFFF00000000)
    if (word >> (var - 6)) & 1:
        return ALL
    return np.uint64(0)


@njit(cache=True)
def prove_one(ops, W0, W1, W2, BT, CB, DEP, tau, cls, w, K, L, nq, nbits):
    inS = np.zeros(nq, np.bool_)
    inB = np.zeros(nbits, np.bool_)
    outw = np.zeros(nq, np.bool_)
    outb = np.zeros(nbits, np.bool_)
    uw = np.zeros(nq, np.bool_)
    ub = np.zeros(nbits, np.bool_)
    fresh_ops = np.zeros(4096, np.int64)
    cut = walk(ops, W0, W1, W2, BT, CB, DEP, tau, cls, w, K, L, nq, nbits, inS, inB, 0, 0, outw, outb, fresh_ops,
               uw, ub)
    if cut < 0:
        return -1, 99
    nfo = walk(ops, W0, W1, W2, BT, CB, DEP, tau, cls, w, K, L, nq, nbits, inS, inB, 1, cut, outw, outb,
               fresh_ops, uw, ub)
    nv = nfo
    for x in range(nq):
        if outw[x]:
            nv += 1
    for x in range(nbits):
        if outb[x]:
            nv += 1
    if nv > K:
        return -1, nv
    nwd = 1
    if nv > 6:
        nwd = 1 << (nv - 6)
    # slots only for the union of live wires / bits; row 0 is a zero row for everything else
    wslot = np.zeros(nq, np.int64)
    nws = 1
    for x in range(nq):
        if uw[x] or outw[x]:
            wslot[x] = nws
            nws += 1
    bslot = np.zeros(nbits, np.int64)
    nbs = 1
    for x in range(nbits):
        if ub[x] or outb[x]:
            bslot[x] = nbs
            nbs += 1
    V = np.zeros((nws, nwd), np.uint64)
    B = np.zeros((nbs, nwd), np.uint64)
    var = 0
    for x in range(nq):
        if outw[x]:
            s = wslot[x]
            for q in range(nwd):
                V[s, q] = pattern(var, q)
            var += 1
    for x in range(nbits):
        if outb[x]:
            s = bslot[x]
            for q in range(nwd):
                B[s, q] = pattern(var, q)
            var += 1
    fvar = np.full(tau + 1 - cut, -1, np.int64)
    for i in range(nfo):
        fvar[fresh_ops[i] - cut] = var
        var += 1
    stack = np.zeros((64, nwd), np.uint64)
    top = np.full(nwd, ALL, np.uint64)
    sp = 0
    m = np.zeros(nwd, np.uint64)
    Pi = np.zeros(nwd, np.uint64)
    for j in range(cut, tau):
        if cls == 7 and j == w:
            for q in range(nwd):
                Pi[q] = V[wslot[W0[j]], q] & V[wslot[W1[j]], q]
        op = ops[j]
        base = op & 0x7F
        if op == 0x81:
            s = bslot[CB[j]]
            for q in range(nwd):
                stack[sp, q] = top[q]
                top[q] = top[q] & B[s, q]
            sp += 1
            continue
        if base == 0x01:
            if sp > 0:
                sp -= 1
                for q in range(nwd):
                    top[q] = stack[sp, q]
            continue
        if not (base == 0x11 or base == 0x12 or base == 0x21 or base == 0x22 or base == 0x31 or base == 0x32
                or base == 0x50 or base == 0x40 or base == 0x41 or base == 0x42):
            continue
        # skip ops that write nothing simulated (targets outside the union)
        x = W0[j]
        y = W1[j]
        z = W2[j]
        if base == 0x11 or base == 0x12:
            if wslot[x] == 0:
                continue
        elif base == 0x21:
            if wslot[y] == 0:
                continue
        elif base == 0x22:
            if wslot[x] == 0 and wslot[y] == 0:
                continue
        elif base == 0x31:
            if wslot[z] == 0:
                continue
        elif base == 0x32:
            if wslot[y] == 0 and wslot[z] == 0:
                continue
        elif base == 0x50:
            if wslot[x] == 0 and bslot[BT[j]] == 0:
                continue
        else:
            if bslot[BT[j]] == 0:
                continue
        cb = CB[j]
        for q in range(nwd):
            m[q] = top[q]
        if cb >= 0:
            s = bslot[cb]
            for q in range(nwd):
                m[q] &= B[s, q]
        sx = wslot[x] if x >= 0 else 0
        sy = wslot[y] if y >= 0 else 0
        sz = wslot[z] if z >= 0 else 0
        if base == 0x11:
            for q in range(nwd):
                V[sx, q] ^= m[q]
        elif base == 0x12:
            for q in range(nwd):
                V[sx, q] &= ~m[q]
        elif base == 0x21:
            for q in range(nwd):
                V[sy, q] ^= V[sx, q] & m[q]
        elif base == 0x22:
            for q in range(nwd):
                dd = (V[sx, q] ^ V[sy, q]) & m[q]
                V[sx, q] ^= dd
                V[sy, q] ^= dd
        elif base == 0x31:
            for q in range(nwd):
                V[sz, q] ^= V[sx, q] & V[sy, q] & m[q]
        elif base == 0x32:
            for q in range(nwd):
                dd = (V[sy, q] ^ V[sz, q]) & V[sx, q] & m[q]
                V[sy, q] ^= dd
                V[sz, q] ^= dd
        elif base == 0x50:
            s = bslot[BT[j]]
            fv = fvar[j - cut]
            for q in range(nwd):
                o = np.uint64(0)
                if fv >= 0:
                    o = pattern(fv, q)
                B[s, q] = (B[s, q] & ~m[q]) | (o & m[q])
                V[sx, q] &= ~m[q]
        elif base == 0x40:
            s = bslot[BT[j]]
            for q in range(nwd):
                B[s, q] &= ~m[q]
        elif base == 0x41:
            s = bslot[BT[j]]
            for q in range(nwd):
                B[s, q] |= m[q]
        else:
            s = bslot[BT[j]]
            for q in range(nwd):
                B[s, q] ^= m[q]
        # row 0 must stay zero (it stands for every unsimulated wire / bit)
        for q in range(nwd):
            V[0, q] = np.uint64(0)
            B[0, q] = np.uint64(0)
    a = wslot[W0[tau]]
    b = wslot[W1[tau]]
    t = wslot[W2[tau]] if W2[tau] >= 0 else 0
    sw = wslot[w] if (w >= 0 and cls != 7) else 0
    base_t = ops[tau] & 0x7F
    for q in range(nwd):
        mm = top[q]
        if CB[tau] >= 0:
            mm &= B[bslot[CB[tau]], q]
        if cls == 6:
            f = V[sw, q]
            mm = ALL
        elif cls == 7:
            f = Pi[q] ^ (V[a, q] & V[b, q])
            mm = ALL
        elif base_t == 0x31:
            p = V[a, q] & V[b, q]
            if cls == 1:
                f = p
            elif cls == 2:
                f = p ^ V[a, q]
            elif cls == 3:
                f = p ^ V[b, q]
            elif cls == 4:
                f = p ^ V[t, q]
            else:
                f = p ^ V[sw, q]
        elif base_t == 0x32:
            f = (V[b, q] ^ V[t, q]) & V[a, q]
        else:
            f = V[a, q] & V[b, q] & V[t, q]
        if (f & mm) != 0:
            return 0, nv
    return 1, nv


@njit(cache=True, parallel=True)
def prove_all(ops, W0, W1, W2, BT, CB, DEP, taus, clss, ws, K, L, nq, nbits):
    n = len(taus)
    res = np.zeros(n, np.int8)
    nvs = np.zeros(n, np.int16)
    for i in prange(n):
        r, nv = prove_one(ops, W0, W1, W2, BT, CB, DEP, taus[i], clss[i], ws[i], K, L, nq, nbits)
        res[i] = r
        nvs[i] = nv
    return res, nvs


def candidates(fp, classes):
    t_op = fp['t_op']
    out_op, out_cls, out_w = [], [], []
    fz, fea, feb, fu, fh = fp['fz'], fp['fea'], fp['feb'], fp['fu'], fp['fh']
    sets = {'Z': fz, 'EA': fea & ~fz, 'EB': feb & ~fz & ~fea, 'U': fu & ~fz, 'H': (fh >= 0) & ~fz & ~fu}
    for name in classes:
        m = sets[name]
        out_op.append(t_op[m])
        out_cls.append(np.full(int(m.sum()), CL[name], np.int8))
        out_w.append(np.where(fh[m] >= 0, fh[m], -1))
    return np.concatenate(out_op), np.concatenate(out_cls), np.concatenate(out_w).astype(np.int64)


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    fp = np.load(sys.argv[2])
    out = sys.argv[3]
    K = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    L = int(sys.argv[5]) if len(sys.argv) > 5 else 400000
    classes = (sys.argv[6] if len(sys.argv) > 6 else 'Z,EA,EB,U,H').split(',')
    ops, W0, W1, W2, BT, CB, DEP = decode(c)
    taus, clss, ws = candidates(fp, classes)
    base = ops[taus] & 0x7F
    keep = (base == 0x31) | (clss == 1)            # CCZ / CSWAP only for Z
    taus, clss, ws = taus[keep], clss[keep], ws[keep]
    print('decoded %.1fs; candidates %d' % (time.time() - t0, len(taus)), flush=True)
    res, nvs = prove_all(ops, W0, W1, W2, BT, CB, DEP, taus.astype(np.int64), clss.astype(np.int64), ws, K, L,
                         c.num_qubits, c.num_bits)
    np.savez(out, op=taus, cls=clss, w=ws, res=res, nv=nvs, K=K, L=L)
    inv = {v: k for k, v in CL.items()}
    for code, nm in ((0x31, 'CCX'), (0x30, 'CCZ'), (0x32, 'CSWAP')):
        for cl in sorted(set(clss.tolist())):
            m = (ops[taus] & 0x7F == code) & (clss == cl)
            if m.any():
                print('%-5s %-2s candidates %6d  proved %6d  not proved %6d  cone too large %6d' % (
                    nm, inv[cl], int(m.sum()), int((res[m] == 1).sum()), int((res[m] == 0).sum()),
                    int((res[m] == -1).sum())), flush=True)
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
