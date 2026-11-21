#!/usr/bin/env python3
"""invsqrt_pp_squeeze.py -- workstream invsqrt_pp (P1): ONE WIRE LESS, by dropping idle values at the peak, then a
re-allocation of wires. Exact on all 3,133,452 inputs.

1. Occupancy (as tools/invsqrt_pp_width.py): at every op the wires that are nonzero on some input, touched by the op,
   or touched anywhere in the open top-level block. The true peak is its maximum; "relief ops" are the ops above the
   target width (default: peak - 1).
2. Drop candidates: a wire w untouched between two top-level events e1 < e2 (an event = a top-level op or a whole
   top-level block) while it holds a nonzero value v. If v, just before e2, is an EXACT affine form of the other wires
   (GF(2) elimination on all inputs, w excluded), w can be X-measured after e1 (HMR -> bit m, debt (-1)^(m v)) and
   rebuilt by X / CNOT gates just before e2, the debt paid by "Z(w) if m": 0 Toffolis, w free strictly between.
3. Greedy cover of the relief ops by candidates (largest new coverage first), refusing a candidate whose rebuild
   uses a wire that is dropped at its rebuild time, or whose drop span contains another chosen rebuild that uses it.
4. Child with the drops; its occupancy must be <= target; then wires are re-allocated by interval colouring of the
   occupied stretches (greedy by start, optimal for intervals), the x register first on its own ids.
5. Check: final x register, zero ancillas and Z-phase equal to the parent for the zero and 3 random HMR transcripts,
   and every Z_POW sees the parent's value. The output is written only if the check passes.
Usage: python invsqrt_pp_squeeze.py <circuit.kmb> <out.kmb> [target width]
"""
import dataclasses
import heapq
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_asap import zpow_values  # noqa: E402
import invsqrt_pp_smallfn as sf  # noqa: E402

X_, CX, Z, HMR, PUSH, POP = 0x11, 0x21, 0x10, 0x50, 0x81, 0x01
LASTMASK = np.uint64((1 << (sim.NX - 64 * (sim.W - 1))) - 1)


def events(C):
    """Top-level events: (start, end, wires) for each top-level op (with qubits) or whole top-level block."""
    ev = []
    i = 0
    while i < C.n:
        if C.base[i] == PUSH and C.level[i] == 0:
            j = C.match[i]
            ws = set()
            for k in range(i, j + 1):
                ws.update(C.qs[k])
            ev.append((i, j, ws))
            i = j + 1
            continue
        if C.qs[i]:
            ev.append((i, i, set(C.qs[i])))
        i += 1
    return ev


def occupancy(C):
    """occ[i] = bitmask (python int) of occupied wires at op i (width semantics), nonzero status tracked per op."""
    blk = {}
    for i in range(C.n):
        if C.base[i] == PUSH and C.level[i] == 0:
            ws = set()
            for j in range(i, C.match[i] + 1):
                ws.update(C.qs[j])
            blk[i] = (C.match[i], sum(1 << w for w in ws))
    occ = [0] * C.n
    cur = [None]
    st = dict(nzm=0, prev=None)

    def hook(i, S, P, bits, act):
        if st['prev'] is None:
            m0 = 0
            for w in np.flatnonzero(S.any(axis=1)).tolist():
                m0 |= 1 << w
            st['nzm'] = m0
        else:
            for w in C.qs[st['prev']]:
                if S[w].any():
                    st['nzm'] |= 1 << w
                else:
                    st['nzm'] &= ~(1 << w)
        st['prev'] = i
        if C.level[i] == 0 and C.base[i] == PUSH:
            cur[0] = blk[i]
        m = st['nzm']
        for w in C.qs[i]:
            m |= 1 << w
        if cur[0] is not None:
            m |= cur[0][1]
            if i >= cur[0][0]:
                cur[0] = None
        occ[i] = m
    sim.run(C, hook=hook)
    return occ


def solve_excl(S, p, excl):
    """p = c0 + sum_{u in R} S[u] exactly (u not in excl)? -> (c0, R) or None."""
    nq = S.shape[0]
    rows = [(np.full(sim.W, sim.ONES, np.uint64), 1)]
    for u in range(nq):
        if u not in excl and S[u].any():
            rows.append((S[u].copy(), 1 << (u + 1)))
    piv = []
    for v, tg in rows:
        v[-1] &= LASTMASK
        for (w, bb, pv, pt) in piv:
            if v[w] & bb:
                v ^= pv
                tg ^= pt
        nz = np.flatnonzero(v)
        if len(nz) == 0:
            continue
        w = int(nz[0])
        piv.append((w, np.uint64(int(v[w]) & -int(v[w])), v, tg))
    t = p.copy()
    t[-1] &= LASTMASK
    comb = 0
    for (w, bb, pv, pt) in piv:
        if t[w] & bb:
            t ^= pv
            comb ^= pt
    if t.any():
        return None
    return comb & 1, [u for u in range(nq) if (comb >> (u + 1)) & 1]


LANES = np.sort(np.random.default_rng(7).choice(sim.NX, size=2048, replace=False)).astype(np.int64)
LW, LB = LANES // 64, (LANES % 64).astype(np.uint64)


def sampled(S, rows):
    """(len(rows), 32) packed sample bits of the given wires."""
    bitsm = ((S[rows][:, LW] >> LB[None, :]) & np.uint64(1)).astype(np.uint8)
    return np.ascontiguousarray(np.packbits(bitsm, axis=1, bitorder='little')).view(np.uint64)


def solve_p1(S, w, excl):
    """v = S[w] = c0 + sum R + S[u1] S[u2] exactly, u's not in excl? -> (c0, R, (u1, u2)) or None.
    Pairs screened on 2,048 sample lanes against the echelon basis of the other wires, then checked exactly."""
    others = [u for u in range(S.shape[0]) if u not in excl and S[u].any()]
    if len(others) < 2:
        return None
    M = sampled(S, others)
    t = sampled(S, [w])[0]
    piv = []
    rows = [np.full(32, np.uint64(0xFFFFFFFFFFFFFFFF))] + [M[k].copy() for k in range(len(others))]
    for v in rows:
        v = v.copy()
        for (wi, bb, pv) in piv:
            if v[wi] & bb:
                v ^= pv
        nz = np.flatnonzero(v)
        if len(nz) == 0:
            continue
        wi = int(nz[0])
        piv.append((wi, np.uint64(int(v[wi]) & -int(v[wi])), v))
    iu, ju = np.triu_indices(len(others), 1)
    X = (M[iu] & M[ju]) ^ t[None, :]
    for (wi, bb, pv) in piv:
        hit = (X[:, wi] & bb) != 0
        X[hit] ^= pv[None, :]
    good = np.flatnonzero(~X.any(axis=1))
    for g in good[:64].tolist():
        u1, u2 = others[iu[g]], others[ju[g]]
        r = solve_excl(S, S[w] ^ (S[u1] & S[u2]), excl)
        if r is not None:
            return r[0], r[1], (u1, u2)
    return None



def solve_p2(S, w, excl, kmax=6):
    """v = S[w] = c0 + sum R + S[a]S[b] + S[c]S[d] exactly? Pairs drawn from a small support of v (greedy on the
    sample lanes: the fewest wires of which v is a function), screened on samples, checked exactly."""
    others = [u for u in range(S.shape[0]) if u not in excl and S[u].any()]
    if len(others) < 2:
        return None
    bitsm = ((S[others][:, LW] >> LB[None, :]) & np.uint64(1)).astype(np.int64)
    vs = ((S[w][LW] >> LB) & np.uint64(1)).astype(np.int64)
    sup = []
    for _ in range(kmax):
        best, bc = None, None
        for t in range(len(others)):
            if t in sup:
                continue
            key = np.zeros(len(vs), np.int64)
            for u in sup + [t]:
                key = key * 2 + bitsm[u]
            k0 = set(key[vs == 0].tolist())
            k1 = set(key[vs == 1].tolist())
            c = len(k0 & k1)
            if bc is None or c < bc:
                best, bc = t, c
        sup.append(best)
        if bc == 0:
            break
    if bc != 0:
        return None
    U = [others[t] for t in sup]
    pairs = [(U[i], U[j]) for i in range(len(U)) for j in range(i + 1, len(U))]
    M = sampled(S, others)
    t = sampled(S, [w])[0]
    piv = []
    rows = [np.full(32, np.uint64(0xFFFFFFFFFFFFFFFF))] + [M[k].copy() for k in range(len(others))]
    for v in rows:
        v = v.copy()
        for (wi, bb, pv) in piv:
            if v[wi] & bb:
                v ^= pv
        nz = np.flatnonzero(v)
        if len(nz) == 0:
            continue
        wi = int(nz[0])
        piv.append((wi, np.uint64(int(v[wi]) & -int(v[wi])), v))
    pos = {u: k for k, u in enumerate(others)}
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            (a, b), (c, d) = pairs[i], pairs[j]
            x = t ^ (M[pos[a]] & M[pos[b]]) ^ (M[pos[c]] & M[pos[d]])
            for (wi, bb, pv) in piv:
                if x[wi] & bb:
                    x = x ^ pv
            if x.any():
                continue
            r = solve_excl(S, S[w] ^ (S[a] & S[b]) ^ (S[c] & S[d]), excl)
            if r is not None:
                return r[0], r[1], [(a, b), (c, d)]
    return None



def solve_small(S, w, excl):
    """('SF', support wires, program) with at most two ANDs, or None (invsqrt_pp_smallfn)."""
    cands = [u for u in range(S.shape[0]) if u not in excl and S[u].any()]
    r = sf.support_tt(S, w, cands)
    if r is None:
        return None
    U, tt, care = r
    prog = sf.synth(tt, care, len(U), allow_nest=False)
    if prog is None:
        return None
    return ('SF', U, prog)



def solve_p1x(S, w, excl, kmax=16):
    """v = S[w] = c0 + sum R + A B exactly, A and B XORs of at most two live wires (a carry MAJ(x,y,c) =
    (x+c)(y+c)+c is one such AND). Factors drawn from the greedy support of v; screened on samples, checked exactly."""
    others = [u for u in range(S.shape[0]) if u not in excl and S[u].any()]
    if len(others) < 2:
        return None
    bitsm = ((S[others][:, LW] >> LB[None, :]) & np.uint64(1)).astype(np.int64)
    vs = ((S[w][LW] >> LB) & np.uint64(1)).astype(np.int64)
    sup = []
    for _ in range(kmax):
        best, bc = None, None
        for t in range(len(others)):
            if t in sup:
                continue
            key = np.zeros(len(vs), np.int64)
            for u in sup + [t]:
                key = key * 2 + bitsm[u]
            c = len(set(key[vs == 0].tolist()) & set(key[vs == 1].tolist()))
            if bc is None or c < bc:
                best, bc = t, c
        sup.append(best)
        if bc == 0:
            break
    U = [others[t] for t in sup]
    M = sampled(S, others)
    pos = {u: k for k, u in enumerate(others)}
    forms = [(u,) for u in U] + [(u, u2) for i, u in enumerate(U) for u2 in U[i + 1:]]
    FV = np.array([M[pos[f[0]]] ^ (M[pos[f[1]]] if len(f) > 1 else 0) for f in forms], dtype=np.uint64)
    t = sampled(S, [w])[0]
    piv = []
    rows = [np.full(32, np.uint64(0xFFFFFFFFFFFFFFFF))] + [M[k].copy() for k in range(len(others))]
    for vv in rows:
        vv = vv.copy()
        for (wi, bb, pv) in piv:
            if vv[wi] & bb:
                vv ^= pv
        nz = np.flatnonzero(vv)
        if len(nz) == 0:
            continue
        wi = int(nz[0])
        piv.append((wi, np.uint64(int(vv[wi]) & -int(vv[wi])), vv))
    iu, ju = np.triu_indices(len(forms), 1)
    X = (FV[iu] & FV[ju]) ^ t[None, :]
    for (wi, bb, pv) in piv:
        hit = (X[:, wi] & bb) != 0
        X[hit] ^= pv[None, :]
    good = np.flatnonzero(~X.any(axis=1))
    for g in good[:64].tolist():
        fa, fb = forms[iu[g]], forms[ju[g]]
        if set(fa) & set(fb):
            continue
        A = S[fa[0]] ^ (S[fa[1]] if len(fa) > 1 else 0)
        B = S[fb[0]] ^ (S[fb[1]] if len(fb) > 1 else 0)
        r = solve_excl(S, S[w] ^ (A & B), excl)
        if r is not None:
            return ('PX', r[0], r[1], fa, fb)
    return None


def popcount(m):
    return bin(m).count('1')


def reallocate(D, k, target):
    """Interval colouring of D's occupied stretches into `target` wires; returns a new Kmb or None."""
    occ = occupancy(D)
    n = D.n
    stretches = []                       # (start, end, wire)
    for w in range(D.nq):
        s = None
        for i in range(n):
            o = (occ[i] >> w) & 1
            if o and s is None:
                s = i
            elif not o and s is not None:
                stretches.append((s, i - 1, w))
                s = None
        if s is not None:
            stretches.append((s, n - 1, w))
    xw = list(D.xw)
    stretches.sort(key=lambda t: (t[0], 0 if t[2] in xw else 1, t[2]))
    reserved = {}                        # x wire id -> its own stretches (pinned to its own id)
    for w in xw:
        if w >= target:
            return None, 'x wire %d >= target' % w
        sw = sorted(t for t in stretches if t[2] == w)
        if not sw or sw[0][0] != 0 or sw[-1][1] != n - 1:
            return None, 'x wire %d is not occupied at the start and the end' % w
        # only the first (input) and the last (output) stretch must sit on the register's own id
        reserved[w] = [(sw[0][0], sw[0][1]), (sw[-1][0], sw[-1][1])]
    busy_until = {p: -1 for p in range(target)}   # last op of the stretch currently on each id
    phys_of = {}

    def clash(p, s, e):
        for (a2, b2) in reserved.get(p, ()):
            if not (e < a2 or b2 < s):
                return True
        return False
    for s, e, w in stretches:
        if w in xw and (s == 0 or e == n - 1):
            p = w
            if busy_until[p] >= s:
                return None, 'x wire %d re-enters at op %d while its id is busy' % (w, s)
        else:
            cands = [p for p in range(target) if busy_until[p] < s and (p not in reserved or not clash(p, s, e))]
            if not cands:
                return None, 'out of wires at op %d' % s

            def slack(q):
                nxt = [a2 for (a2, b2) in reserved.get(q, ()) if a2 > e]
                return (min(nxt) - e) if nxt else 10 ** 9
            # an x wire's free gap that fits goes first (tightest gap), plain wires after
            p = min(cands, key=lambda q: (0 if q in reserved else 1, slack(q), q))
        phys_of[(w, s)] = (p, e)
        busy_until[p] = e
    # per wire, sorted stretch starts for lookup
    starts = {}
    for (w, s), (p, e) in phys_of.items():
        starts.setdefault(w, []).append((s, e, p))
    for w in starts:
        starts[w].sort()

    def m(w, i):
        for s, e, p in starts.get(w, []):
            if s <= i <= e:
                return p
        raise RuntimeError('op %d touches wire %d outside its stretches' % (i, w))
    q3 = [[], [], []]
    q2 = [[], []]
    q1 = []
    for i in range(n):
        q = D.qs[i]
        if len(q) == 3:
            for j in range(3):
                q3[j].append(m(q[j], i))
        elif len(q) == 2:
            for j in range(2):
                q2[j].append(m(q[j], i))
        elif len(q) == 1:
            q1.append(m(q[0], i))
    regs = []
    for name, ids, isq in k.registers:
        ids = np.asarray(ids).copy()
        for t in range(len(ids)):
            if isq[t]:
                ids[t] = m(int(ids[t]), 0) if int(ids[t]) in xw else min(int(ids[t]), target - 1)
        regs.append((name, ids, isq))
    k2 = dataclasses.replace(k, registers=regs, num_qubits=target,
                             q3=tuple(np.asarray(a, dtype=np.uint32) for a in q3),
                             q2=tuple(np.asarray(a, dtype=np.uint32) for a in q2),
                             q1=np.asarray(q1, dtype=np.uint32))
    return k2, 'ok'


def main():
    t0 = time.time()
    src, dst = sys.argv[1], sys.argv[2]
    C = sim.load(src)
    k0 = kmb.read(src)
    occ = occupancy(C)
    cnt = [popcount(o) for o in occ]
    peak = max(cnt)
    target = int(sys.argv[3]) if len(sys.argv) > 3 else peak - 1
    relief = [i for i in range(C.n) if cnt[i] > target]
    print('declared %d, true peak %d, target %d, relief ops %d (%.1fs)' % (C.nq, peak, target, len(relief),
                                                                         time.time() - t0), flush=True)
    print('registers:', [(r[0], len(r[1]), int(np.asarray(r[2]).sum())) for r in k0.registers])
    ev = events(C)
    # per wire: indices of the events touching it
    touch = {w: [] for w in range(C.nq)}
    for e, (a, b, ws) in enumerate(ev):
        for w in ws:
            touch[w].append(e)
    relief_set = set(relief)
    rel_sorted = np.array(relief, np.int64)
    cands = []                           # (w, drop after op, rebuild before op, covered relief ops)
    for w in range(C.nq):
        if w in C.xw and os.environ.get('SQ_X', '1') != '1':
            continue
        tl = touch[w]
        for e1, e2 in zip(tl[:-1], tl[1:]):
            a = ev[e1][1]
            b = ev[e2][0]
            if b - a < 2:
                continue
            lo = np.searchsorted(rel_sorted, a + 1)
            hi = np.searchsorted(rel_sorted, b)
            if hi <= lo:
                continue
            if not ((occ[a + 1] >> w) & 1):
                continue                 # already free there
            cands.append((w, a, b, set(rel_sorted[lo:hi].tolist())))
    print('candidate idle stretches covering relief ops: %d' % len(cands), flush=True)
    # exact affine relation of each candidate's value at its rebuild point
    want = {}
    for c in cands:
        want.setdefault(c[2], []).append(c)
    rel = {}

    P1MODE = os.environ.get('SQ_P1', '0') == '1'

    def hook(i, S, P, bits, act):
        if i in want:
            for c in want[i]:
                w = c[0]
                r = solve_excl(S, S[w], {w})
                if r is not None:
                    rel[(w, c[1], c[2])] = (r[0], r[1], None)
                elif P1MODE:
                    r = solve_p1(S, w, {w})
                    if r is None and os.environ.get('SQ_P2', '0') == '1':
                        r = solve_p2(S, w, {w})
                    if r is None and os.environ.get('SQ_SMALL', '0') == '1':
                        r = solve_small(S, w, {w})
                    if r is None and os.environ.get('SQ_P1X', '0') == '1':
                        r = solve_p1x(S, w, {w})
                    if r is not None:
                        rel[(w, c[1], c[2])] = r
    sim.run(C, hook=hook)
    good = [c for c in cands if rel.get((c[0], c[1], c[2])) is not None]
    n1 = sum(1 for c in good if rel[(c[0], c[1], c[2])][0] in ('SF', 'PX') or rel[(c[0], c[1], c[2])][2] is not None)
    print('of them with an exact rebuild: %d (affine %d, one CCX %d) (%.1fs)' % (
        len(good), len(good) - n1, n1, time.time() - t0), flush=True)

    def prs(r):
        if r[0] == 'PX':
            return [None]
        if r[0] == 'SF':
            return [None] * {'aff': 0, 'one': 1, 'two': 2}[r[2][0]]
        if r[2] is None:
            return []
        return list(r[2]) if isinstance(r[2], list) else [r[2]]

    def supp(r):
        if r[0] == 'PX':
            return list(r[2]) + list(r[3]) + list(r[4])
        if r[0] == 'SF':
            return list(r[1])
        return list(r[1]) + [u for pq in prs(r) for u in pq]

    def greedy(good):
        chosen = []
        uncovered = set(relief_set)
        for price in (0, 1, 2):
            pool = [c for c in good if len(prs(rel[(c[0], c[1], c[2])])) == price]
            while uncovered:
                best, bestn = None, 0
                for c in pool:
                    if c in chosen:
                        continue
                    nn = len(c[3] & uncovered)
                    if nn <= bestn:
                        continue
                    w, a, b = c[0], c[1], c[2]
                    R = supp(rel[(w, a, b)])
                    ok = True
                    for d in chosen:
                        dw, da, db = d[0], d[1], d[2]
                        if dw == w and not (b < da or db < a):
                            ok = False
                            break
                        if dw in R and da < b <= db:
                            ok = False
                            break
                        if w in supp(rel[(dw, da, db)]) and a < db <= b:
                            ok = False
                            break
                    if ok:
                        best, bestn = c, nn
                if best is None:
                    break
                chosen.append(best)
                uncovered -= best[3]
        return chosen, uncovered
    chosen, uncovered = greedy(good)
    if uncovered and os.environ.get('SQ_EARLY', '0') == '1':
        # shorter stretches: rebuild at an earlier top-level event start, after the uncovered relief ops covered
        starts = np.array([e[0] for e in ev], np.int64)
        extra = []
        for (w, a, b, cov) in cands:
            U = [r for r in cov if r in uncovered]
            if not U:
                continue
            lo = np.searchsorted(starts, max(U) + 1)
            hi = np.searchsorted(starts, b)
            for c in starts[lo:hi][:int(os.environ.get('SQ_EARLYN', '12'))].tolist():
                if (w, a, c) in rel:
                    continue
                cv = set(r for r in cov if r < c)
                extra.append((w, a, c, cv))
        want2 = {}
        for c in extra:
            want2.setdefault(c[2], []).append(c)

        def hook2(i, S, P, bits, act):
            if i in want2:
                for c in want2[i]:
                    w = c[0]
                    r = solve_excl(S, S[w], {w})
                    if r is not None:
                        rel[(w, c[1], c[2])] = (r[0], r[1], None)
                    elif P1MODE:
                        r = solve_p1(S, w, {w})
                        if r is None and os.environ.get('SQ_P2', '0') == '1':
                            r = solve_p2(S, w, {w})
                        if r is None and os.environ.get('SQ_SMALL', '0') == '1':
                            r = solve_small(S, w, {w})
                        if r is None and os.environ.get('SQ_P1X', '0') == '1':
                            r = solve_p1x(S, w, {w})
                        if r is not None:
                            rel[(w, c[1], c[2])] = r
        sim.run(C, hook=hook2)
        g2 = [c for c in extra if rel.get((c[0], c[1], c[2])) is not None]
        print('early rebuild points tried %d, exact %d (%.1fs)' % (len(extra), len(g2), time.time() - t0), flush=True)
        chosen, uncovered = greedy(good + g2)
    nccx = sum(len(prs(rel[(c[0], c[1], c[2])])) for c in chosen)
    print('chosen drops: %d (with one CCX: %d = +%d Toffolis); relief ops left uncovered: %d' % (
        len(chosen), nccx, nccx, len(uncovered)), flush=True)
    for c in chosen:
        print('  drop wire %d after op %d, rebuild before op %d (%d ops), relation %s' % (
            c[0], c[1], c[2], c[2] - c[1], rel[(c[0], c[1], c[2])]))
    CAP = os.environ.get('SQ_CAP', '0') == '1'
    if CAP:
        chosen = [c for c in chosen if rel[(c[0], c[1], c[2])][2] is None]
        nccx = 0
        print('CAPACITY MODE: %d affine drops kept, width unchanged' % len(chosen), flush=True)
    if uncovered and not CAP:
        u = sorted(uncovered)
        print('uncovered relief ops (first 20):', u[:20])
        print('NO CHILD: the target width needs drops with nonlinear rebuilds')
        return
    maxccx = int(os.environ.get('SQ_MAXCCX', '30'))
    if nccx > maxccx:
        print('NO CHILD: %d rebuild CCX exceed SQ_MAXCCX=%d (a wire is worth about Toffoli-eq / (Q - 1))' % (
            nccx, maxccx))
        return
    # child with the drops
    nb = k0.num_bits
    inserts = {}
    for c in chosen:
        w, a, b = c[0], c[1], c[2]
        rr = rel[(w, a, b)]
        m = nb
        nb += 1
        inserts.setdefault(a, []).append((HMR, (w,), m, 0))
        if rr[0] == 'PX':
            _, c0x, Rx, fa, fb = rr
            gx = ([(X_, (w,), 0, 0)] if c0x else []) + [(CX, (u, w), 0, 0) for u in Rx]
            pre = []
            if len(fa) > 1:
                pre.append((CX, (fa[1], fa[0]), 0, 0))
            if len(fb) > 1:
                pre.append((CX, (fb[1], fb[0]), 0, 0))
            gx += pre + [(0x31, (fa[0], fb[0], w), 0, 0)] + list(reversed(pre))
            inserts.setdefault(b - 1, []).extend(gx + [(PUSH, (), 0, m), (Z, (w,), 0, 0), (POP, (), 0, 0)])
            continue
        if rr[0] == 'SF':
            gates, _ = sf.emit(rr[2], rr[1], w)
            inserts.setdefault(b - 1, []).extend(gates + [(PUSH, (), 0, m), (Z, (w,), 0, 0), (POP, (), 0, 0)])
            continue
        c0, R, pr = rr
        seq = ([(X_, (w,), 0, 0)] if c0 else []) + [(CX, (u, w), 0, 0) for u in R] + \
              [(0x31, (q0, q1_, w), 0, 0) for (q0, q1_) in prs((c0, R, pr))] + \
              [(PUSH, (), 0, m), (Z, (w,), 0, 0), (POP, (), 0, 0)]
        inserts.setdefault(b - 1, []).extend(seq)
    kD = dataclasses.replace(kw.edit(k0, [], inserts), num_bits=nb)
    tmp = dst + '.drops.kmb'
    kw.write(tmp, kD)
    D = sim.load(tmp)
    if CAP:
        S0, P0, _, _ = sim.run(C)
        ok = True
        for seed in (0, 1, 2, 3):
            oc = None
            if seed:
                dr = np.random.default_rng(seed).integers(0, 2, size=D.n)
                oc = (lambda j, d=dr: int(d[j]))
            S1, P1, _, _ = sim.run(D, outcome=oc)
            e = bool(np.array_equal(S1, S0)) and bool(np.array_equal(P1, P0))
            print('  seed', seed, 'final state and Z-phase equal to the parent:', e, flush=True)
            ok &= e
        zp = zpow_values(C) == zpow_values(D)
        print('  every Z_POW sees the parent value:', zp)
        print('CHECK', 'PASS' if ok and zp else 'FAIL', round(time.time() - t0, 1), 's')
        if ok and zp:
            os.replace(tmp, dst)
        else:
            os.remove(tmp)
        return
    occD = occupancy(D)
    pD = max(popcount(o) for o in occD)
    print('child with drops: true peak %d (%.1fs)' % (pD, time.time() - t0), flush=True)
    if pD > target:
        print('NO CHILD: the drops do not reach the target')
        return
    k2, why = reallocate(D, kmb.read(tmp), target)
    if k2 is None:
        print('NO CHILD: re-allocation failed:', why)
        return
    tmp2 = dst + '.tmp'
    kw.write(tmp2, k2)
    E = sim.load(tmp2)
    print('re-allocated: %d qubits' % E.nq, flush=True)
    # check against the parent
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=E.n)
            oc = (lambda j, d=dr: int(d[j]))
        S1, P1, _, _ = sim.run(E, outcome=oc)
        e = bool(np.array_equal(S1[E.xw], S0[C.xw])) and bool(np.array_equal(P1, P0))
        e &= not S1[[w for w in range(E.nq) if w not in E.xw]].any()
        e &= not S0[[w for w in range(C.nq) if w not in C.xw]].any()
        print('  seed', seed, 'final x, zero ancillas and Z-phase equal to the parent:', e, flush=True)
        ok &= e
    zp = zpow_values(C) == zpow_values(E)
    print('  every Z_POW sees the parent value:', zp)
    print('CHECK', 'PASS' if ok and zp else 'FAIL', round(time.time() - t0, 1), 's')
    if ok and zp:
        os.replace(tmp2, dst)
        os.remove(tmp)
    else:
        os.remove(tmp2)


if __name__ == '__main__':
    main()
