#!/usr/bin/env python3
"""r86_chem_g35_dump.py: seat R86-CHEM-G35. Runs builder C4's chem_c4_build.main() UNCHANGED (imported, never edited)
with one wrapper around C4's dead-row hook (chem_c4_dead.solve_dead_row, installed by C4_DEAD_ROW): when the engine
reaches generic group R86_GI (default 3 = group 4), the kill system the engine hands the hook (pre-kill residual,
particular residual R0, shared kernel N, pool rows, solutions) and EVERY usable dead cell of the group's outputs with
its effect compensated through the pool on the killed column (C4's own rule, chem_c4_dead L465-489, applied to every
output instead of one) are written to R86_DUMP (npz + json). C4's hook then runs as before.

  R86_DUMP=out.npz CHEM_BROT=7 CHEM_HB=77 C4_DROP=1 C4_DEAD_ROW=3:2:72 python -I r86_chem_g35_dump.py \
      table.npz --k 9 --wirebasis --plan 6,6,6:1:0:0,6:1:0:1 [--emit DIR]
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_loaders_core as core  # noqa: E402
import chem_c4_dead as deadm  # noqa: E402
import chem_c4_drop as dr  # noqa: E402

bits = core.bits
OUT = os.environ["R86_DUMP"]
GI = int(os.environ.get("R86_GI", "3"))


def comp_cells(dec, outs, live, killed, P, usable, nr, cells):
    """every usable dead cell of the outputs, compensated on the killed column through the pool (C4's rule)."""
    an = dec.an
    t = killed[0]
    liveu = np.uint32(live)
    piv = {}
    for pi in range(len(usable)):
        v = 0
        for r in range(nr):
            if (int(P[pi][r]) >> t) & 1:
                v |= 1 << r
        tag = 1 << pi
        while v:
            hh = v.bit_length() - 1
            if hh not in piv:
                piv[hh] = (v, tag)
                break
            v ^= piv[hh][0]
            tag ^= piv[hh][1]

    def express(c):
        tag = 0
        while c:
            hh = c.bit_length() - 1
            if hh not in piv:
                return None
            c ^= piv[hh][0]
            tag ^= piv[hh][1]
        return tag
    out = []
    for o in outs:
        pend = None
        for c in cells:
            if c[0] != o:
                continue
            E = deadm._effect(an, nr, c[1], c[2], c[4])
            cv = 0
            for r in range(nr):
                if (int(E[r]) >> t) & 1:
                    cv |= 1 << r
            tag = express(cv)
            group = [c]
            if tag is None:
                if pend is None:
                    pend = (E, cv, c)
                    continue
                E = E ^ pend[0]
                tag = express(cv ^ pend[1])
                group = [pend[2], c]
                pend = None
                if tag is None:
                    continue
            Ec = E.copy()
            for pi in bits(tag):
                Ec ^= P[pi]
            assert not ((Ec >> np.uint32(t)) & 1).any()
            out.append(dict(o=int(o), E=[int(x) for x in E], tag=str(tag),
                            cells=[[int(c_[0]), int(c_[1]), int(c_[2]), [int(c_[3][0]), int(c_[3][1])], bool(c_[4]),
                                    int(c_[5])] for c_ in group],
                            Ec=[int(x) for x in (Ec & liveu)]))
    return out


def dump(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, base):
    an = dec.an
    cells = [c for c in deadm.HOOK_ROW["cells"] if c[0] in outs]
    comp = comp_cells(dec, outs, live, killed, P, usable, nr, cells)
    lc = bits(live)
    ylocal = []
    for ym in an.ymasks:
        m = sum(1 << ci for ci, tt in enumerate(lc) if ym >> tt & 1)
        if m and m not in ylocal:
            ylocal.append(m)
    np.savez(OUT, X=dec.X[:, outs].T.copy(), R0=np.stack(R0), N=np.stack(N) if N else np.zeros((0, nr), np.uint32),
             P=P, usable=np.array(usable), outs=np.array(outs), live=live, killed=np.array(killed), nr=nr, b=b,
             base=base, K=K)
    meta = dict(k=int(an.k), gi=GI, outs=[int(o) for o in outs], live=int(live), killed=[int(x) for x in killed],
                nr=int(nr), b=int(b), base=int(base), K=int(K), npool=len(usable),
                sols=[str(c) for c in sols], kern=[str(c) for c in kern], ylocal=ylocal,
                ymasks=[int(m) for m in an.ymasks], ncells_raw=len(cells), comp=comp,
                prods=[[int(p[0]), hex(int(p[1])), int(p[2]), hex(int(p[3])), int(p[4])] for p in
                       [dec.prods[i] for i in usable]])
    json.dump(meta, open(OUT.replace(".npz", ".json"), "w"))
    cnt = {}
    for c in comp:
        cnt[c["o"]] = cnt.get(c["o"], 0) + 1
    print(f"R86 DUMP gi {GI}: outs {list(outs)} pool {len(usable)} K {K} base {base} nr {nr} live {live:#x} killed "
          f"{list(killed)}; dead cells raw {len(cells)}, compensated unknowns per output {cnt}; ylocal {ylocal} -> {OUT}",
          flush=True)


_orig = deadm.install_row


def install_row(cells, m, target_gi, conj, rk=None, specs=None):
    _orig(cells, m, target_gi, conj, rk, specs)
    hook = dr.DEAD_ROW_HOOK

    def wrapped(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base):
        if getattr(dec, "_gi", -1) == GI and not os.path.exists(OUT):
            dump(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, base)
        return hook(dec, outs, live, killed, R0, N, K, P, usable, sols, kern, nr, b, rows_of, base)
    dr.DEAD_ROW_HOOK = wrapped


deadm.install_row = install_row

if __name__ == "__main__":
    import chem_c4_build
    chem_c4_build.main()
