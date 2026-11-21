#!/usr/bin/env python3
"""invsqrt_pp_pubrw.py -- workstream P1 (sixth holder): census of three PUBLIC exact rewrites on a small .kmb.

Exact over every input lane (invsqrt: all 3,133,452 x; chemistry angle circuit: all 2,048 selection values).
Wire values are simulated with every measurement outcome = 1 (every conditioned block runs, so values inside
blocks are known); a second run with outcome = 0 checks that the top-level values do not depend on outcomes.

  c3  DEAD WRITE: CCX a b t whose target t is next touched only by XOR-writes (X t, CX . t, CCX . . t) and then by
      an HMR t at the same or an outer block context (unconditioned op).  Rewrite: delete the CCX, add CZ a b
      conditioned on the HMR's bit after the HMR (inside PUSH of the CCX's block bits when the CCX was deeper).
      Exact iff a*b at the HMR equals a*b at the CCX on every lane (checked) and the block bits are unchanged.
  c2  MEASURE-OUT: CCX a b t ... CCX a' b' t (same block context and condition), equal products on every lane,
      t touched in between, and a wire z that is 0 on every lane at the first CCX and untouched until after the
      second.  Rewrite: CCX a b z; CX z t ... CX z t; HMR z; CZ a' b' if m.  -1 Toffoli.  (Pairs with t untouched
      in between are reported as P2: both Toffolis cancel.)
  c4  FEWEST WIRES: the largest number of wires non-zero on some lane at one op boundary (interval colouring of
      the busy segments; interval graphs are perfect, so this is the exact minimum for relabelled segments).
Executed weight of an op: 0.5 per enclosing block plus 0.5 for its own condition (invsqrt outcomes are fair coins).

Usage: python invsqrt_pp_pubrw.py <circuit.kmb> [json_out]
"""
import hashlib
import json
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402

ONES = np.uint64(0xFFFFFFFFFFFFFFFF)


def lanes(C):
    """Initial state over all lanes: invsqrt x register or chemistry selection register."""
    regs = {r[0]: np.asarray(r[1]).tolist() for r in C.k.registers}
    if 'selection' in regs:
        # chemistry, the verifier's chain: loader k runs forward on words = batch k-1 (0 for k = 0) and again in
        # the reverse pass on words = batch k (it must restore batch k-1). Lanes 0..ns-1 = forward,
        # ns..2ns-1 = reverse, over every selection value (a superset of the verified addresses).
        sel = regs['selection']
        ns = 1 << len(sel)
        Wd = max(1, (ns + 63) // 64)
        a = np.arange(Wd * 64, dtype=np.int64) % ns
        S = np.zeros((C.nq, Wd), np.uint64)
        for j, w in enumerate(sel):
            S[w] = np.packbits(((a >> j) & 1).astype(np.uint8), bitorder='little').view(np.uint64)
        prev = chain_words(C)  # {(register, bit): vector} = batch k-1
        for (r, j), v in prev.items():
            S[regs[r][j]] = v
        F = run_values(C, S.copy())
        S2 = np.concatenate([S, S], axis=1)
        for r in regs:
            if r.startswith('word_'):
                for w in regs[r]:
                    S2[w, Wd:] = F[w]
        return S2, 'selection', regs
    return sim.init_state(C), 'x', regs


def run_values(C, S):
    """Run every op on the lane state S (outcome 0; values of the loaders do not depend on outcomes)."""
    bits, stk = [0] * C.nbits, []
    for i in range(C.n):
        b = C.base[i]
        act = all(stk)
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act or (C.cnd[i] >= 0 and not bits[C.cnd[i]]):
            continue
        if b == 0x50:
            S[C.qs[i][0]][:] = 0
            bits[C.bit[i]] = 0
        else:
            sim.apply(C, i, S, np.zeros(S.shape[1], np.uint64), bits, None, [0])
    return S


def chain_words(C):
    """batch k-1 as {(register, bit): lane vector} for C = angle_k_<k>.kmb (empty for k = 0)."""
    import re
    m = re.search(r'angle_k_(\d+)\.kmb$', C.path)
    if not m or int(m.group(1)) == 0:
        return {}
    d = os.path.dirname(C.path)
    words = {}
    for k in range(int(m.group(1))):
        Ck = sim.load(os.path.join(d, 'angle_k_%d.kmb' % k))
        Ck.path = ''
        rg = {r[0]: np.asarray(r[1]).tolist() for r in Ck.k.registers}
        ns = 1 << len(rg['selection'])
        Wd = max(1, (ns + 63) // 64)
        a = np.arange(Wd * 64, dtype=np.int64) % ns
        S = np.zeros((Ck.nq, Wd), np.uint64)
        for j, w in enumerate(rg['selection']):
            S[w] = np.packbits(((a >> j) & 1).astype(np.uint8), bitorder='little').view(np.uint64)
        for (r, j), v in words.items():
            S[rg[r][j]] = v
        S = run_values(Ck, S)
        words = {(r, j): S[w].copy() for r in rg if r.startswith('word_') for j, w in enumerate(rg[r])}
    return words


def hsh(v):
    return hashlib.blake2b(v.tobytes(), digest_size=12).digest()


def contexts(C):
    """Block stack (tuple of PUSH indices) of every op."""
    ctx = [()] * C.n
    st = []
    for i in range(C.n):
        if C.base[i] == kmb.PUSH:
            ctx[i] = tuple(st)
            st.append(i)
        elif C.base[i] == kmb.POP:
            st.pop()
            ctx[i] = tuple(st)
        else:
            ctx[i] = tuple(st)
    return ctx


def roles(C, i):
    """(wire, role) touched by op i: role 'x' XOR-written, 'h' HMR, 'o' anything else (read / other write)."""
    b, q = C.base[i], C.qs[i]
    if b == 0x11:
        return [(q[0], 'x')]
    if b == 0x21:
        return [(q[0], 'o'), (q[1], 'x')]
    if b == 0x31:
        return [(q[0], 'o'), (q[1], 'o'), (q[2], 'x')]
    if b == 0x50:
        return [(q[0], 'h')]
    return [(w, 'o') for w in q]


def simulate(C, outcome_one, want_prod_at, want_zero=True):
    """Run all ops; returns per-op product hashes for CCX, requested product hashes, zero flags after each op."""
    S, kind, regs = lanes(C)
    bits = [0] * C.nbits
    stk = []
    prod = {}
    extra = {}
    zafter = {}
    nz = np.array([bool(S[w].any()) for w in range(C.nq)])
    peak, peak_at = int(nz.sum()), 0
    for i in range(C.n):
        b = C.base[i]
        act = all(stk)
        if i in want_prod_at:
            for (a, bb) in want_prod_at[i]:
                extra[(i, a, bb)] = hsh(S[a] & S[bb])
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07) or not act:
            continue
        if C.cnd[i] >= 0 and not bits[C.cnd[i]]:
            continue
        q = C.qs[i]
        if b == 0x31:
            p = S[q[0]] & S[q[1]]
            prod[i] = (hsh(p), not p.any())
        if b == 0x50:
            m = 1 if outcome_one else 0
            S[q[0]][:] = 0
            bits[C.bit[i]] = m
        else:
            sim.apply(C, i, S, np.zeros(S.shape[1], np.uint64), bits, None, [0])
        for w in q:
            z = not S[w].any()
            zafter[(i, w)] = z
            nz[w] = not z
        c = int(nz.sum())
        if c > peak:
            peak, peak_at = c, i
    return prod, extra, zafter, peak, peak_at, kind, regs


def assign_c2(rows, taken_ops=()):
    """Disjoint c2 sites: each op once, each free wire holds one pair at a time (earliest second op first).
    rows: (i, k, zs, wt, nbetween).  Returns [(i, k, z, wt)]."""
    used = set(taken_ops)
    occ = {}
    out = []
    for (i, k, zs, w, nb) in sorted(rows, key=lambda r: (r[1], r[0])):
        if i in used or k in used:
            continue
        for z in zs:
            if all(k2 < i or i2 > k for (i2, k2) in occ.get(z, [])):
                occ.setdefault(z, []).append((i, k))
                used.update((i, k))
                out.append((i, k, z, w))
                break
    return out


def main(path=None, quiet=False):
    path = path or sys.argv[1]
    C = sim.load(path)
    C.path = path
    ctx = contexts(C)
    wt = [0.5 ** (len(ctx[i]) + (1 if C.cnd[i] >= 0 else 0)) for i in range(C.n)]
    # touch lists per wire
    touch = [[] for _ in range(C.nq)]
    for i in range(C.n):
        for (w, r) in roles(C, i):
            touch[w].append((i, r))
    pos = {}
    for w in range(C.nq):
        for j, (i, r) in enumerate(touch[w]):
            pos[(i, w)] = j
    ccx = [i for i in range(C.n) if C.base[i] == 0x31]
    # ---- c3 candidates (structural)
    c3 = []
    want = {}
    for i in ccx:
        a, b, t = C.qs[i]
        j = pos[(i, t)] + 1
        while j < len(touch[t]) and touch[t][j][1] == 'x':
            j += 1
        if j >= len(touch[t]) or touch[t][j][1] != 'h':
            continue
        h = touch[t][j][0]
        if C.cnd[h] >= 0 or ctx[h] != ctx[i][:len(ctx[h])]:
            continue
        c3.append((i, h))
        want.setdefault(h, []).append((a, b))
    # ---- c2 candidates (structural): same target, same context and condition
    c2 = []
    last = {}
    for i in ccx:
        a, b, t = C.qs[i]
        key = (t, ctx[i], C.cnd[i])
        for i0 in last.get(key, [])[-8:]:
            c2.append((i0, i))
        last.setdefault(key, []).append(i)
    # ---- simulate (outcomes 1 and 0)
    prod1, extra1, z1, peak1, pat1, kind, regs = simulate(C, True, want)
    prod0, extra0, z0, peak0, pat0, _, _ = simulate(C, False, want)
    top_mismatch = sum(1 for i in prod0 if len(ctx[i]) == 0 and C.cnd[i] < 0 and prod0[i] != prod1.get(i))
    zf = {k: z1[k] and z0.get(k, True) for k in z1}  # zero after op on every lane in both runs

    def zero_at(w, i):
        """wire w is 0 just before op i (every lane) and the touch list index of its next touch at/after i."""
        lst = touch[w]
        lo, hi = 0, len(lst)
        while lo < hi:
            m = (lo + hi) // 2
            if lst[m][0] < i:
                lo = m + 1
            else:
                hi = m
        nxt = lst[lo][0] if lo < len(lst) else C.n
        if lo == 0:
            z = not (kind == 'x' and w in regs.get('x', [])) and not (
                kind == 'selection' and any(w in regs[r] for r in regs if r == 'selection' or r.startswith('word_')))
        else:
            z = zf.get((lst[lo - 1][0], w), False)
        return z, nxt

    res = {'circuit': path, 'nq': C.nq, 'ops': C.n, 'top_level_outcome_dependence': top_mismatch}
    # c3
    rows3 = []
    for (i, h) in c3:
        if i not in prod1:
            continue
        a, b, t = C.qs[i]
        ok = extra1.get((h, a, b)) == prod1[i][0] and extra0.get((h, a, b)) == prod0.get(i, (None,))[0]
        rows3.append((i, h, ok, wt[i], prod1[i][1], ctx[i] == ctx[h] and C.cnd[i] < 0))
    res['c3_candidates'] = len(rows3)
    res['c3_sites_exact'] = sum(1 for r in rows3 if r[2])
    res['c3_saved_per_shot'] = sum(r[3] for r in rows3 if r[2])
    res['c3_sites_exact_same_block'] = sum(1 for r in rows3 if r[2] and r[5])
    res['c3_list'] = [[r[0], r[1], r[2], r[3], r[5]] for r in rows3]
    # c2
    rows2, p2 = [], []
    for (i, k) in c2:
        if i not in prod1 or k not in prod1 or prod1[i][0] != prod1[k][0] or prod0.get(i) != prod0.get(k):
            continue
        t = C.qs[i][2]
        between = touch[t][pos[(i, t)] + 1:pos[(k, t)]]
        if not between:
            p2.append((i, k, wt[i]))
            continue
        zs = []
        for z in range(C.nq):
            if z in C.qs[i] or z in C.qs[k]:
                continue
            zz, nxt = zero_at(z, i)
            if zz and nxt > k:
                zs.append(z)
        if zs:
            # the condition bit (if any) must not be rewritten between the two
            cb = C.cnd[i]
            if cb >= 0 and any(C.bit[j] == cb for j in range(i + 1, k)):
                continue
            rows2.append((i, k, zs, wt[i], len(between)))
    fit2 = assign_c2(rows2)
    res['c2_pairs_with_free_wire'] = len(rows2)
    res['c2_sites'] = len(fit2)
    res['c2_saved_per_shot'] = sum(r[3] for r in fit2)
    res['c2_list'] = [[r[0], r[1], r[2], r[3], r[4]] for r in rows2]
    p2d, usedp = [], set()
    for (i, k, w) in p2:
        if i not in usedp and k not in usedp:
            usedp.update((i, k))
            p2d.append((i, k, w))
    res['P2_pairs'] = len(p2d)
    res['P2_saved_per_shot'] = sum(2 * r[2] for r in p2d)
    res['P2_list'] = [list(r) for r in p2d]
    # c2 without a free wire: how many equal-product pairs exist at all (to say why absent)
    res['c2_equal_product_pairs_touched'] = sum(
        1 for (i, k) in c2 if i in prod1 and k in prod1 and prod1[i][0] == prod1[k][0]
        and touch[C.qs[i][2]][pos[(i, C.qs[i][2])] + 1:pos[(k, C.qs[i][2])]])
    res['c4_peak_busy_outcome1'] = peak1
    res['c4_peak_busy_outcome0'] = peak0
    res['c4_peak_op'] = pat1
    res['width'] = C.nq
    res['ccx_static'] = len(ccx)
    res['ccx_zero_product'] = sum(1 for i in prod1 if prod1[i][1] and prod0.get(i, (0, True))[1])
    out = {k: v for k, v in res.items() if not k.endswith('_list')}
    if not quiet:
        print(json.dumps(out, indent=1))
    if not quiet and len(sys.argv) > 2:
        with open(sys.argv[2], 'w') as f:
            json.dump(res, f, indent=1)
    return res


if __name__ == '__main__':
    main()
