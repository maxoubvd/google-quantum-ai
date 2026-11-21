#!/usr/bin/env python3
"""invsqrt_pp_width.py -- workstream invsqrt_pp (P1): the TRUE peak width of an invsqrt .kmb.

At every op, the wires that are nonzero on some input or touched by the op (ops of a conditioned block count for the
whole block) are occupied. The maximum over the circuit is the width a re-allocation of wires (interval colouring
of the occupied stretches, greedy by start, which is optimal for intervals) reaches without changing any gate.
Also prints the widest stage cuts (INFO_TAG names) and how many ops sit at the peak.
Usage: python invsqrt_pp_width.py <circuit.kmb>
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import invsqrt_pp_sim as sim  # noqa: E402


def main():
    C = sim.load(sys.argv[1])
    k = kmb.read(sys.argv[1])
    tags = list(k.tags)
    occ = np.zeros(C.n, np.int64)
    nz = np.zeros(C.nq, bool)
    tagname = {}
    ti = 0
    for i in range(C.n):
        if C.base[i] == 0x06:
            tagname[i] = tags[ti] if ti < len(tags) else ''
            ti += 1
    blk_touch = {}
    for i in range(C.n):
        if C.base[i] == kmb.PUSH and C.level[i] == 0:
            ws = set()
            for j in range(i, C.match[i] + 1):
                ws.update(C.qs[j])
            blk_touch[i] = (C.match[i], ws)
    cur_blk = [None]

    def hook(i, S, P, bits, act):
        if C.level[i] == 0 and C.base[i] == kmb.PUSH:
            cur_blk[0] = blk_touch[i]
        live = set(w for w in range(C.nq) if S[w].any())
        live.update(C.qs[i])
        if cur_blk[0] is not None:
            live.update(cur_blk[0][1])
            if i >= cur_blk[0][0]:
                cur_blk[0] = None
        occ[i] = len(live)
    sim.run(C, hook=hook)
    pk = int(occ.max())
    at = np.flatnonzero(occ == pk)
    print('declared qubits', C.nq, '; true peak width', pk, '; ops at the peak', len(at), 'of', C.n,
          '; first / last at the peak', int(at[0]), int(at[-1]))
    for w in range(pk, max(pk - 6, 0), -1):
        print('  ops with width >= %d: %d' % (w, int((occ >= w).sum())))
    # stage of the peak ops by the last tag before them
    tpos = sorted(tagname)
    seen = {}
    for i in at[:2000]:
        j = np.searchsorted(tpos, i) - 1
        nm = tagname[tpos[j]] if j >= 0 else '-'
        seen[nm] = seen.get(nm, 0) + 1
    print('  peak ops by last tag:', sorted(seen.items(), key=lambda x: -x[1])[:8])


if __name__ == '__main__':
    main()
