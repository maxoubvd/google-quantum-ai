#!/usr/bin/env python3
"""invsqrt_pp_alap.py -- workstream invsqrt_pp (P1): compute-late + measure-early reorder of an invsqrt .kmb (exact).

Compute-late: every top-level unconditioned X / CX / CCX / CSWAP that writes a wire which is ZERO on all inputs just
before it (an initialisation) moves down to just before the next op that touches any of its wires (an op inside a
conditioned block counts as the whole block). Disjoint gates commute and these gates make no phase, so the circuit is
unchanged; the initialised wire stays zero longer. Then measure-early (tools/invsqrt_pp_links.asap_order). Writes the
kmb and checks it against the input on all 3,133,452 inputs (4 transcripts, every Z_POW value).
Usage: python invsqrt_pp_alap.py <in.kmb> <out.kmb>
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402
from invsqrt_pp_asap import zpow_values  # noqa: E402

WRITERS = (0x11, 0x21, 0x31, 0x32)


def main():
    t0 = time.time()
    src, dst = sys.argv[1], sys.argv[2]
    C = sim.load(src)
    # zero flags before each op (zero transcript; top-level values do not depend on the transcript)
    zero_before = {}
    z = np.ones(C.nq, bool)
    for w in C.xw:
        z[w] = False

    def hook(i, S, P, bits, act):
        if C.base[i] in WRITERS and C.level[i] == 0 and C.cnd[i] < 0:
            tgt = C.qs[i][-1] if C.base[i] != 0x32 else None
            if tgt is not None:
                zero_before[i] = not S[tgt].any()
    sim.run(C, hook=hook)
    # anchors: the next op touching each wire, block ops count as their whole block (anchor = the block's PUSH)
    blk_start = [-1] * C.n
    blk_wires = {}
    i = 0
    while i < C.n:
        if C.base[i] == kmb.PUSH and C.level[i] == 0:
            j = C.match[i]
            ws = set()
            for k in range(i, j + 1):
                ws.update(C.qs[k])
                blk_start[k] = i
            blk_wires[i] = ws
            i = j + 1
            continue
        i += 1
    # simple forward scan instead: for each candidate, walk to the first later op (top-level or block) touching its wires
    moved = {}
    for i, zb in zero_before.items():
        if not zb:
            continue
        ws = set(C.qs[i])
        k = i + 1
        while k < C.n:
            if blk_start[k] >= 0:
                bs = blk_start[k]
                if ws & blk_wires[bs]:
                    break
                k = C.match[bs] + 1
                continue
            if ws & set(C.qs[k]) or C.base[k] in (0x06, 0x07):
                break
            k += 1
        if k > i + 1:
            moved[i] = k            # place op i right before op k
    # build the edit: delete moved ops, insert each before its anchor (= after op anchor-1)
    k0 = kmb.read(src)
    ops = np.asarray(k0.ops)
    delete = sorted(moved)
    inserts = {}
    for i in sorted(moved):
        a = moved[i] - 1
        while a in moved:            # the op before the anchor is itself moved away: attach to an earlier stay op
            a -= 1
        inserts.setdefault(a, []).append((int(ops[i]), C.qs[i], C.bit[i], 0))
    k2 = kw.edit(k0, delete, inserts)
    tmp = dst + '.cl.kmb'
    kw.write(tmp, k2)
    # measure-early on top
    D = sim.load(tmp)
    from invsqrt_pp_links import asap_order
    order, nm = asap_order(D)
    pos = {op: q for q, op in enumerate(order)}
    mv = [h for h in range(D.n) if D.base[h] == 0x50 and D.level[h] == 0 and D.cnd[h] < 0
          and pos[h] > 0 and order[pos[h] - 1] != h - 1]
    kk = kmb.read(tmp)
    ops2 = np.asarray(kk.ops)
    ins2 = {}
    for h in mv:
        ins2.setdefault(order[pos[h] - 1], []).append((int(ops2[h]), D.qs[h], D.bit[h], 0))
    kw.write(dst, kw.edit(kk, mv, ins2))
    os.remove(tmp)
    E = sim.load(dst)
    print('initialisations moved later: %d; HMRs moved earlier: %d' % (len(moved), len(mv)), flush=True)
    S0, P0, _, _ = sim.run(C)
    ok = True
    for seed in (0, 1, 2, 3):
        oc = None
        if seed:
            dr = np.random.default_rng(seed).integers(0, 2, size=E.n)
            oc = (lambda j, d=dr: int(d[j]))
        S1, P1, _, _ = sim.run(E, outcome=oc)
        ok &= bool(np.array_equal(S1, S0)) and bool(np.array_equal(P1, P0))
    zp = zpow_values(C) == zpow_values(E)
    print('CHECK', 'PASS' if ok and zp else 'FAIL', round(time.time() - t0, 1), 's')


if __name__ == '__main__':
    main()
