#!/usr/bin/env python3
"""invsqrt_pp_blocks.py -- workstream invsqrt_pp (P1): which conditioned correction blocks of an invsqrt .kmb
could be replaced by Clifford phase gates (0 Toffolis), now or deferred to a later point of the circuit.

A top-level block B (PUSH m ... POP) that contains Toffoli-class ops, executed, must restore every wire and apply
a phase (-1)^g(x). If g equals a quadratic form q of the wire values at a top-level position t after the HMR that
writes m (q = c0 + sum c_i w_i + sum c_ij w_i w_j on ALL inputs), then B can be replaced by
"PUSH m; Z / CZ / NEG network of q; POP" at t: Clifford, 0 Toffolis. The test is first run on 2,048 sample lanes
(GF(2) elimination), then every hit is checked exactly on all 3,133,452 inputs (pass 2).

Usage: python invsqrt_pp_blocks.py <circuit.kmb> <out.json> [step=8] [nsamp=2048]
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402

TOF = (0x30, 0x31, 0x32)


def sample_bits(S, lanes):
    """(nq, ns) uint8 bits of every wire at the sample lanes."""
    w = lanes >> 6
    b = (lanes & 63).astype(np.uint64)
    return ((S[:, w] >> b) & np.uint64(1)).astype(np.uint8)


def vec_sample(v, lanes):
    w = lanes >> 6
    b = (lanes & 63).astype(np.uint64)
    return ((v[w] >> b) & np.uint64(1)).astype(np.uint8)


def quad_basis(Wb):
    """Wb: (ns, nq) uint8 -> (ns, 1 + nq + nq(nq-1)/2) uint8, with the index list of every column."""
    ns, nq = Wb.shape
    iu, ju = np.triu_indices(nq, 1)
    prods = Wb[:, iu] & Wb[:, ju]
    A = np.concatenate([np.ones((ns, 1), np.uint8), Wb, prods], axis=1)
    cols = [('1',)] + [('w', i) for i in range(nq)] + [('ww', int(i), int(j)) for i, j in zip(iu, ju)]
    return A, cols


def gf2_solve(A, G):
    """Solve A c = g for every column g of G over GF(2). Returns (solvable bool[m], solutions uint8[m, k])."""
    n, k = A.shape
    m = G.shape[1]
    M = np.concatenate([A, G], axis=1)
    tot = k + m
    padc = (-tot) % 64
    if padc:
        M = np.concatenate([M, np.zeros((n, padc), np.uint8)], axis=1)
    Mp = np.packbits(M, axis=1, bitorder='little').view(np.uint64).copy()
    row = 0
    piv = []
    for col in range(k):
        if row >= n:
            break
        w, b = col >> 6, np.uint64(col & 63)
        colbits = (Mp[row:, w] >> b) & np.uint64(1)
        nz = np.flatnonzero(colbits)
        if len(nz) == 0:
            continue
        p = row + int(nz[0])
        if p != row:
            Mp[[row, p]] = Mp[[p, row]]
        mask = ((Mp[:, w] >> b) & np.uint64(1)).astype(bool)
        mask[row] = False
        Mp[mask] ^= Mp[row]
        piv.append(col)
        row += 1
    full = np.unpackbits(Mp.view(np.uint8), axis=1, bitorder='little')[:, :tot]
    Gp = full[:, k:]
    solvable = ~(Gp[row:].any(axis=0))
    sol = np.zeros((m, k), np.uint8)
    for r, col in enumerate(piv):
        sol[:, col] = Gp[r]
    return solvable, sol


def find_blocks(C):
    """Top-level blocks holding Toffoli-class ops, with the HMR that writes their condition bit."""
    out = []
    last_writer = {}
    for i in range(C.n):
        b = C.base[i]
        if C.level[i] == 0 and b == kmb.PUSH:
            j = C.match[i]
            nt = sum(1 for k in range(i + 1, j) if C.base[k] in TOF)
            if nt:
                texp = sum(0.5 ** (C.level[k] + (1 if C.cnd[k] >= 0 else 0)) for k in range(i + 1, j)
                           if C.base[k] in TOF)
                wr = last_writer.get(C.cnd[i], -1)
                out.append({'push': i, 'pop': j, 'cond': C.cnd[i], 'writer': wr,
                            'writer_kind': C.base[wr] if wr >= 0 else -1,
                            'mwire': C.qs[wr][0] if wr >= 0 and C.base[wr] == 0x50 else -1,
                            'tof': nt, 'tof_exp': texp})
        if C.bit[i] >= 0 and C.base[i] in (0x50, 0x40, 0x41, 0x42):
            last_writer[C.bit[i]] = i
    return out


def main():
    t0 = time.time()
    path, outp = sys.argv[1], sys.argv[2]
    step = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    ns = int(sys.argv[4]) if len(sys.argv) > 4 else 2048
    C = sim.load(path)
    blocks = find_blocks(C)
    print('blocks with Toffolis', len(blocks), 'expected Toffolis', sum(b['tof_exp'] for b in blocks), flush=True)
    rng = np.random.default_rng(12345)
    lanes = np.sort(rng.choice(sim.NX, size=ns, replace=False)).astype(np.int64)
    by_writer = {}
    for k, b in enumerate(blocks):
        by_writer.setdefault(b['writer'], []).append(k)
    by_push = {b['push']: k for k, b in enumerate(blocks)}
    V = {}
    G = {}
    pos = []          # (op index, sample matrix (ns, nq))
    info = {}

    def hook(i, S, P, bits, act):
        if i in by_writer:
            for k in by_writer[i]:
                V[k] = S[blocks[k]['mwire']].copy() if blocks[k]['mwire'] >= 0 else None
        if i in by_push:
            k = by_push[i]
            r0, g0, nr0 = sim.block_effect(C, S, bits, i, outcome=lambda j: 0)
            r1, g1, nr1 = sim.block_effect(C, S, bits, i, outcome=lambda j: 1)
            G[k] = g0
            v = V.get(k)
            info[k] = {'restored0': r0, 'restored1': r1, 'g_indep': bool(np.array_equal(g0, g1)),
                       'g_eq_v': bool(v is not None and np.array_equal(g0, v)), 'rot_in_block': nr0 + nr1,
                       'g_weight': int(np.unpackbits(g0.view(np.uint8)).sum())}
        if C.level[i] == 0 and C.base[i] != kmb.POP and (i % step == 0 or i in by_push):
            pos.append((i, sample_bits(S, lanes).T.copy()))

    sim.run(C, hook=hook)
    print('pass 1 done', len(pos), 'positions', round(time.time() - t0, 1), 's', flush=True)
    order = sorted(range(len(blocks)), key=lambda k: blocks[k]['writer'])
    Gs = np.stack([vec_sample(G[k], lanes) for k in range(len(blocks))], axis=1)   # (ns, nb)
    hits = {k: [] for k in range(len(blocks))}
    sols = {}
    for pi, (i, Wb) in enumerate(pos):
        act = [k for k in order if blocks[k]['writer'] < i]
        if not act:
            continue
        A, cols = quad_basis(Wb)
        ok, sol = gf2_solve(A, Gs[:, act])
        for t, k in enumerate(act):
            if ok[t]:
                hits[k].append(i)
                if k not in sols:
                    sols[k] = (i, sol[t].copy())
        if pi % 200 == 0:
            print('  pos', pi, i, 'hits so far', sum(1 for k in hits if hits[k]), round(time.time() - t0, 1), 's',
                  flush=True)
    # pass 2: exact check of the first hit of every block on all inputs
    want = {}
    for k, (i, s) in sols.items():
        want.setdefault(i, []).append(k)
    exact = {}

    def hook2(i, S, P, bits, act):
        if i in want:
            for k in want[i]:
                _, s = sols[k]
                nq = C.nq
                iu, ju = np.triu_indices(nq, 1)
                q = np.zeros(sim.W, np.uint64)
                if s[0]:
                    q ^= sim.ONES
                for w in range(nq):
                    if s[1 + w]:
                        q ^= S[w]
                for t in np.flatnonzero(s[1 + nq:]).tolist():
                    q ^= S[iu[t]] & S[ju[t]]
                d = q ^ G[k]
                # padded lanes are copies of x = 1: compare only real lanes
                nbad = int(np.unpackbits(d.view(np.uint8), bitorder='little')[:sim.NX].sum())
                exact[k] = nbad

    sim.run(C, hook=hook2)
    res = []
    for k, b in enumerate(blocks):
        r = dict(b)
        r.update(info.get(k, {}))
        r['hits'] = len(hits[k])
        r['first_hit'] = hits[k][0] if hits[k] else -1
        r['hit_at_push'] = b['push'] in hits[k]
        r['exact_bad'] = exact.get(k, -1)
        if k in sols:
            s = sols[k][1]
            r['form_terms'] = int(s.sum())
        res.append(r)
    with open(outp, 'w') as f:
        json.dump(res, f, indent=0)
    n_ok = [r for r in res if r['first_hit'] >= 0 and r['exact_bad'] == 0]
    print('blocks', len(res), 'restored all', all(r.get('restored0') and r.get('restored1') for r in res),
          'g independent of inner outcomes', sum(r.get('g_indep', False) for r in res),
          'g == measured value', sum(r.get('g_eq_v', False) for r in res))
    print('quadratic somewhere (exact on all inputs):', len(n_ok), 'expected Toffolis',
          sum(r['tof_exp'] for r in n_ok), '; at the block itself:', sum(1 for r in n_ok if r['hit_at_push']))
    print('sample hits refuted exactly:', sum(1 for r in res if r['first_hit'] >= 0 and r['exact_bad'] != 0))
    print('done', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
