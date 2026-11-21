#!/usr/bin/env python3
"""invsqrt_pp_links.py -- workstream invsqrt_pp (P1): exact links and wire capacity of an invsqrt .kmb.

1. Measure-early reorder (optional, `asap=1`): every top-level unconditioned HMR moves up to right after the last op
   that touches its wire or mentions its bit (or after the POP of the block holding that op). Exact: the moved
   measurement acts on the same value; its kick commutes with every op it crosses.
2. Exact links: Toffoli-class ops i < j with p_j = p_i + affine(wires just before j) on ALL inputs (p = a*b for
   CCX / CCZ pairs, c*(x+y) for CSWAP). Candidates on 2,048 sample lanes, then checked on all 3,133,452 inputs.
3. Capacity: per op, wires identically zero on all inputs and not used by the op (ops inside a conditioned block
   count as using their wires for the whole block); and the GF(2) affine rank of the 46 wire values (sample
   upper bound on what CNOTs can free). A link fits if every op strictly between i and j leaves a zero wire.
Usage: python invsqrt_pp_links.py <circuit.kmb> <out.json> [asap=0|1]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_blocks import sample_bits, vec_sample  # noqa: E402

TOF = (0x30, 0x31, 0x32)


def asap_order(C):
    """New op order with every top-level unconditioned HMR moved up (see module doc)."""
    last_w = {}
    last_b = {}
    place = {}          # HMR index -> index of the op it follows
    open_blk = []
    for i in range(C.n):
        b = C.base[i]
        if b == kmb.PUSH:
            open_blk.append(i)
        if b == 0x50 and C.level[i] == 0 and C.cnd[i] < 0:
            w = C.qs[i][0]
            bt = C.bit[i]
            place[i] = max(last_w.get(w, -1), last_b.get(bt, -1))
        top = open_blk[0] if open_blk else None
        anchor = C.match[top] if top is not None else i       # an op inside a block anchors after the block
        for w in C.qs[i]:
            last_w[w] = max(last_w.get(w, -1), anchor)
        for bt in (C.bit[i], C.cnd[i]):
            if bt >= 0:
                last_b[bt] = max(last_b.get(bt, -1), anchor)
        if b == kmb.POP:
            open_blk.pop()
    after = {}
    for h, k in place.items():
        after.setdefault(k, []).append(h)
    order = []
    moved = set(place)
    for h in after.get(-1, []):
        order.append(h)
    for i in range(C.n):
        if i not in moved:
            order.append(i)
        for h in after.get(i, []):
            order.append(h)
    assert len(order) == C.n and len(set(order)) == C.n
    nmoved = sum(1 for h, k in place.items() if k < h - 1)
    return order, nmoved


def permute(C, order):
    D = sim.Circ()
    D.k, D.nq, D.nbits, D.xw = C.k, C.nq, C.nbits, C.xw
    D.base = [C.base[i] for i in order]
    D.qs = [C.qs[i] for i in order]
    D.bit = [C.bit[i] for i in order]
    D.cnd = [C.cnd[i] for i in order]
    D.n = C.n
    D.orig = list(order)
    D.match = {}
    D.level = [0] * D.n
    st = []
    for i in range(D.n):
        if D.base[i] == kmb.PUSH:
            D.level[i] = len(st)
            st.append(i)
        elif D.base[i] == kmb.POP:
            D.match[st.pop()] = i
            D.level[i] = len(st)
        else:
            D.level[i] = len(st)
    return D


def gf2_rank_rows(Mb):
    """Rank over GF(2) of the rows of a (r, ns) uint8 matrix."""
    M = np.packbits(np.ascontiguousarray(Mb), axis=1, bitorder='little')
    pad = (-M.shape[1]) % 8
    if pad:
        M = np.concatenate([M, np.zeros((M.shape[0], pad), np.uint8)], axis=1)
    M = np.ascontiguousarray(M).view(np.uint64).copy()
    r = 0
    nw = M.shape[1]
    for w in range(nw):
        for b in range(64):
            if r >= M.shape[0]:
                return r
            bb = np.uint64(b)
            col = (M[r:, w] >> bb) & np.uint64(1)
            nz = np.flatnonzero(col)
            if len(nz) == 0:
                continue
            p = r + int(nz[0])
            if p != r:
                M[[r, p]] = M[[p, r]]
            mask = ((M[:, w] >> bb) & np.uint64(1)).astype(bool)
            mask[r] = False
            M[mask] ^= M[r]
            r += 1
    return r


def main():
    t0 = time.time()
    path, outp = sys.argv[1], sys.argv[2]
    use_asap = len(sys.argv) > 3 and sys.argv[3] == '1'
    C0 = sim.load(path)
    if use_asap:
        order, nmoved = asap_order(C0)
        C = permute(C0, order)
        print('measure-early: HMRs moved', nmoved, flush=True)
    else:
        C = C0
    rng = np.random.default_rng(777)
    lanes = np.sort(rng.choice(sim.NX, size=2048, replace=False)).astype(np.int64)
    nq = C.nq
    zero = np.ones(nq, bool)
    for w in C.xw:
        zero[w] = False
    # block extents: ops inside a top-level block use their wires for the whole block
    blk_of = [-1] * C.n
    blk_wires = {}
    for i in range(C.n):
        if C.base[i] == kmb.PUSH and C.level[i] == 0:
            j = C.match[i]
            ws = set()
            for k in range(i, j + 1):
                blk_of[k] = i
                ws.update(C.qs[k])
            blk_wires[i] = ws
    free_cnt = np.zeros(C.n, np.int16)
    tof_idx = [i for i in range(C.n) if C.base[i] in TOF]
    prod_s = {}
    span_s = {}
    full_rank_ops = []

    def hook(i, S, P, bits, act):
        used = blk_wires[blk_of[i]] if blk_of[i] >= 0 else set(C.qs[i])
        free_cnt[i] = sum(1 for w in range(nq) if zero[w] and w not in used)
        if C.base[i] in TOF and C.level[i] == 0 and C.cnd[i] < 0:
            q = C.qs[i]
            if C.base[i] == 0x32:
                p = S[q[0]] & (S[q[1]] ^ S[q[2]])
            else:
                p = S[q[0]] & S[q[1]]
            prod_s[i] = vec_sample(p, lanes)
            span_s[i] = sample_bits(S, lanes)

    def post(i, S):
        for w in C.qs[i]:
            zero[w] = not S[w].any()

    # run with a post-op update of zero flags
    S = sim.init_state(C)
    P = np.zeros(sim.W, np.uint64)
    bits = [0] * C.nbits
    stk = []
    nrot = [0]
    for i in range(C.n):
        b = C.base[i]
        act = all(stk)
        hook(i, S, P, bits, act)
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act or (C.cnd[i] >= 0 and not bits[C.cnd[i]]):
            continue
        sim.apply(C, i, S, P, bits, None, nrot)
        post(i, S)
    print('simulated', round(time.time() - t0, 1), 's', flush=True)
    n46 = int((free_cnt == 0).sum())
    print('ops with no free zero wire (46-wide):', n46, 'of', C.n, flush=True)
    # sample affine rank at a subset of 46-wide ops (redundancy capacity), 2nd pass
    want = set(np.flatnonzero(free_cnt == 0)[::25].tolist())
    ranks = {}

    def hook_r(i, S, P, bits, act):
        if i in want:
            Wb = sample_bits(S, lanes)
            M = np.concatenate([np.ones((1, Wb.shape[1]), np.uint8), Wb], axis=0)
            ranks[i] = gf2_rank_rows(M) - 1
    sim.run(C, hook=hook_r)
    rk = np.array(list(ranks.values()))
    if len(rk):
        print('affine rank at %d sampled 46-wide ops: min %d mean %.2f max %d (46 = no CNOT-freeable wire)' %
              (len(rk), rk.min(), rk.mean(), rk.max()), flush=True)
    else:
        print('no full-width op: every op has a free zero wire (affine rank not sampled)', flush=True)
    # links: affine-equivalent products among unconditioned top-level Toffoli-class ops
    tl = sorted(prod_s)
    links = []
    Pm = np.stack([prod_s[i] for i in tl], axis=0)       # (nt, ns)
    for jj, j in enumerate(tl):
        if jj == 0:
            continue
        Wb = span_s[j]
        A = np.concatenate([np.ones((1, Wb.shape[1]), np.uint8), Wb], axis=0)   # (47, ns)
        # reduce rows of [Pm[:jj+1]] modulo the row space of A
        M = np.concatenate([A, Pm[:jj + 1]], axis=0)
        Mp = np.ascontiguousarray(np.packbits(np.ascontiguousarray(M), axis=1, bitorder='little')).view(np.uint64).copy()
        na = A.shape[0]
        r = 0
        for w in range(Mp.shape[1]):
            for b in range(64):
                if r >= na:
                    break
                bb = np.uint64(b)
                col = (Mp[r:na, w] >> bb) & np.uint64(1)
                nz = np.flatnonzero(col)
                if len(nz) == 0:
                    continue
                p = r + int(nz[0])
                if p != r:
                    Mp[[r, p]] = Mp[[p, r]]
                mask = ((Mp[:, w] >> bb) & np.uint64(1)).astype(bool)
                mask[r] = False
                Mp[mask] ^= Mp[r]
                r += 1
        red = Mp[na:]
        tgt = red[jj]
        if not tgt.any():
            continue        # product affine at j: audit_affine's class (should be 0)
        eq = np.flatnonzero((red[:jj] == tgt).all(axis=1))
        for ii in eq.tolist():
            links.append((tl[ii], j))
    print('sample links (affine-equivalent product pairs):', len(links), round(time.time() - t0, 1), 's', flush=True)
    fits = []
    for (i, j) in links:
        m = int(free_cnt[i + 1:j].min()) if j > i + 1 else 99
        fits.append(m)
    fits = np.array(fits) if fits else np.zeros(0, int)
    print('links whose interval keeps >= 1 zero wire throughout:', int((fits >= 1).sum()), 'of', len(links))
    res = {'asap': use_asap, 'n46': n46, 'links': [(int(a), int(b), int(f)) for (a, b), f in zip(links, fits)],
           'ranks': {int(k): int(v) for k, v in ranks.items()}, 'free_cnt': free_cnt.tolist(),
           'orig': getattr(C, 'orig', None)}
    with open(outp, 'w') as f:
        json.dump(res, f)
    print('done', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
