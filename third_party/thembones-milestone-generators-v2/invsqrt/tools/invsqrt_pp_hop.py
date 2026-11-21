#!/usr/bin/env python3
"""invsqrt_pp_hop.py -- workstream invsqrt_pp (P1): retained products with HOLDERS THAT HOP between free wires.

Same chains and exact relations as tools/invsqrt_pp_hold.py (links file of invsqrt_pp_links.py, asap=0), but a
chain's holder need not stay on one wire: before every op where the parent needs the holder's wire, the held value
moves to another free wire by two CNOTs (CX h h'; CX h' h). A wire is free at op t when it is zero on all inputs
before t and not touched by op t (a top-level block counts as touching all its wires for its whole length).
Admission: chains by earliest end, admitted if at every op of the chain the free wires outnumber the holders already
admitted there (with free moves this count is also sufficient for an assignment). Emission as in hold.py:
  first op   CCX a b t -> CCX a b h ; CX h t      (CSWAP c x y -> CX x y ; CCX c y h ; CX x y ; CX h x ; CX h y)
  later ops  lambda into h (CX / X), then CX h t (CSWAP: CX h x ; CX h y): 0 Toffolis
  last op    as later, then HMR h -> m ; CZ a b if m (CSWAP: CZ c x if m ; CZ c y if m)
Check: all-input simulation of the child (zero and 3 random HMR transcripts) gives the parent's final state and
Z-phase, and every Z_POW sees the parent's value. The output is written only if the check passes.
Usage: python invsqrt_pp_hop.py <parent.kmb> <links.json> <out.kmb> [maxchains=0 (all)]
"""
import dataclasses
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_blocks import gf2_solve, sample_bits, vec_sample  # noqa: E402
from invsqrt_pp_hold import chains_from, prod_of  # noqa: E402
from invsqrt_pp_squeeze import occupancy  # noqa: E402
from invsqrt_pp_asap import zpow_values  # noqa: E402

CX, CCX, X, HMR = 0x21, 0x31, 0x11, 0x50
CZ_IF = 0x20 | 0x80


def relations(C, chains):
    """Exact affine lambda of every non-first chain op (None if not exact)."""
    member = {}
    for ci, ch in enumerate(chains):
        for m, i in enumerate(ch):
            member[i] = (ci, m)
    rng = np.random.default_rng(99)
    lanes = np.sort(rng.choice(sim.NX, size=2048, replace=False)).astype(np.int64)
    D, rel = {}, {}

    def hook(i, S, P, bits, act):
        if i not in member or not act:
            return
        ci, m = member[i]
        d = prod_of(C, i, S)
        if m > 0:
            tgt = d ^ D[chains[ci][m - 1]]
            Wb = sample_bits(S, lanes).T
            A = np.concatenate([np.ones((Wb.shape[0], 1), np.uint8), Wb], axis=1)
            ok, sol = gf2_solve(A, vec_sample(tgt, lanes)[:, None])
            lam = None
            if ok[0]:
                s = sol[0]
                q = np.zeros(sim.W, np.uint64)
                if s[0]:
                    q ^= sim.ONES
                ws = []
                for w in range(C.nq):
                    if s[1 + w] and S[w].any():
                        q ^= S[w]
                        ws.append(w)
                dd = q ^ tgt
                if int(np.unpackbits(dd.view(np.uint8), bitorder='little')[:sim.NX].sum()) == 0:
                    lam = (int(s[0]), ws)
            rel[i] = lam
        D[i] = d
    sim.run(C, hook=hook)
    return rel


def main():
    t0 = time.time()
    ppath, lpath, outp = sys.argv[1], sys.argv[2], sys.argv[3]
    maxch = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    L = json.load(open(lpath))
    assert not L['asap'], 'links must be computed on the parent order (asap=0)'
    C = sim.load(ppath)
    links = [tuple(x) for x in L['links']]
    chains = chains_from(links, set())
    rel = relations(C, chains)
    good = []
    for ch in chains:
        if any(rel.get(i) is None for i in ch[1:]):
            k = next(m for m in range(1, len(ch)) if rel.get(ch[m]) is None)
            ch = ch[:k]
        if len(ch) >= 2 and all(C.level[i] == 0 and C.cnd[i] < 0 for i in ch):
            good.append(ch)
    print('chains %d, exact and top-level %d (%.1fs)' % (len(chains), len(good), time.time() - t0), flush=True)
    occ = occupancy(C)
    xm = sum(1 << w for w in C.xw)
    full = (1 << C.nq) - 1
    freem = [full & ~(o | xm) for o in occ]
    cap = np.array([bin(f).count('1') for f in freem], np.int64)
    used = np.zeros(C.n, np.int64)
    good.sort(key=lambda ch: (ch[-1], ch[0]))
    adm = []
    for ch in good:
        lo, hi = ch[0], ch[-1]
        if np.all(cap[lo:hi + 1] - used[lo:hi + 1] >= 1):
            adm.append(ch)
            used[lo:hi + 1] += 1
            if maxch and len(adm) >= maxch:
                break
    print('admitted %d chains (%d Toffolis) of %d; capacity: ops with no free wire %d of %d' % (
        len(adm), sum(len(c) - 1 for c in adm), len(good), int((cap == 0).sum()), C.n), flush=True)
    # assignment with hops: sweep ops in order
    starts = {}
    for ci, ch in enumerate(adm):
        starts.setdefault(ch[0], []).append(ci)
    endop = {ci: ch[-1] for ci, ch in enumerate(adm)}
    where = {}                    # chain -> current wire
    holder_at = {}                # (chain, op) -> wire used at that chain op
    moves = {}                    # op t -> [(chain, from, to)] inserted after op t-1 (i.e. before op t)
    member = {}
    for ci, ch in enumerate(adm):
        for i in ch:
            member.setdefault(i, []).append(ci)
    nmoves = 0
    for t in range(C.n):
        fr = freem[t]
        held = set(where.values())
        # holders whose wire is not free at t must move
        for ci in sorted(where):
            w = where[ci]
            if not (fr >> w) & 1:
                cand = [u for u in range(C.nq) if (fr >> u) & 1 and u not in held]
                if not cand:
                    raise RuntimeError('no free wire for chain %d at op %d' % (ci, t))
                # prefer the wire that stays free longest
                best, bl = None, -1
                for u in cand:
                    k = t
                    while k < C.n and k <= endop[ci] and (freem[k] >> u) & 1:
                        k += 1
                    if k > bl:
                        best, bl = u, k
                moves.setdefault(t, []).append((ci, w, best))
                held.discard(w)
                held.add(best)
                where[ci] = best
                nmoves += 1
        for ci in starts.get(t, []):
            cand = [u for u in range(C.nq) if (fr >> u) & 1 and u not in held]
            if not cand:
                raise RuntimeError('no free wire to start chain %d at op %d' % (ci, t))
            best, bl = None, -1
            for u in cand:
                k = t
                while k < C.n and k <= endop[ci] and (freem[k] >> u) & 1:
                    k += 1
                if k > bl:
                    best, bl = u, k
            where[ci] = best
            held.add(best)
        for ci in member.get(t, []):
            holder_at[(ci, t)] = where[ci]
        for ci in [c for c in where if endop[c] == t]:
            del where[ci]
    print('holder moves: %d (%.1fs)' % (nmoves, time.time() - t0), flush=True)
    # emit
    k = kmb.read(ppath)
    newbit = k.num_bits
    delete = []
    inserts = {}
    saved = 0
    for t, mv in moves.items():
        seq = []
        for ci, a, b in mv:
            seq += [(CX, (a, b), 0, 0), (CX, (b, a), 0, 0)]
        inserts.setdefault(t - 1, []).extend(seq)
    for ci, ch in enumerate(adm):
        for m, i in enumerate(ch):
            h = holder_at[(ci, i)]
            q = C.qs[i]
            seq = []
            if m == 0:
                if C.base[i] == 0x32:
                    c, x, y = q
                    seq = [(CX, (x, y), 0, 0), (CCX, (c, y, h), 0, 0), (CX, (x, y), 0, 0), (CX, (h, x), 0, 0),
                           (CX, (h, y), 0, 0)]
                else:
                    a, bq, tt = q
                    seq = [(CCX, (a, bq, h), 0, 0), (CX, (h, tt), 0, 0)]
            else:
                const, ws = rel[i]
                for w in ws:
                    seq.append((CX, (w, h), 0, 0))
                if const:
                    seq.append((X, (h,), 0, 0))
                if C.base[i] == 0x32:
                    c, x, y = q
                    seq += [(CX, (h, x), 0, 0), (CX, (h, y), 0, 0)]
                else:
                    a, bq, tt = q
                    seq += [(CX, (h, tt), 0, 0)]
                saved += 1
                if m == len(ch) - 1:
                    seq.append((HMR, (h,), newbit, 0))
                    if C.base[i] == 0x32:
                        c, x, y = q
                        seq += [(CZ_IF, (c, x), 0, newbit), (CZ_IF, (c, y), 0, newbit)]
                    else:
                        a, bq, tt = q
                        seq += [(CZ_IF, (a, bq), 0, newbit)]
                    newbit += 1
            delete.append(i)
            # moves scheduled before op i are inserted after op i-1, which runs before this op's rewrite
            inserts[i] = seq + inserts.get(i, [])
    k2 = dataclasses.replace(kw.edit(k, delete, inserts), num_bits=newbit)
    tmp = outp + '.tmp'
    kw.write(tmp, k2)
    print('written; Toffoli-class ops removed %d, new bits %d' % (saved, newbit - k.num_bits), flush=True)
    Cc = sim.load(tmp)
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=Cc.n)
            oc = (lambda j, d=dr: int(d[j]))
        S2, P2, _, _ = sim.run(Cc, outcome=oc)
        e = bool(np.array_equal(S2, S0)) and bool(np.array_equal(P2, P0))
        print('  check seed', seed, 'final state and Z-phase equal:', e, flush=True)
        ok &= e
    zp = zpow_values(C) == zpow_values(Cc)
    print('  every Z_POW sees the parent value:', zp)
    print('CHECK', 'PASS' if ok and zp else 'FAIL', 'saved Toffolis', saved, round(time.time() - t0, 1), 's')
    if ok and zp:
        os.replace(tmp, outp)
    else:
        os.remove(tmp)


if __name__ == '__main__':
    main()
