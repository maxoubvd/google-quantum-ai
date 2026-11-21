#!/usr/bin/env python3
"""r46_dqi_scan_held.py -- R46-DQI-SCAN: HELD IDLE WIRES of a DQI decoder, static (no samples).

One pass over the ops. Per wire: the writes since its last reset (R / HMR), the op that made it (first write after
the reset: kind and operands) and every IDLE GAP = two consecutive touches of the wire at least G ops apart while it
may be nonzero. SWAP exchanges the two wires' records (the value moves). A gap of a wire made by ONE CCX and nothing
else is a single AND bit held idle: it can be measured out after the touch that opens the gap (HMR, CZ on its two
operands) and remade before the touch that closes it for 1 Toffoli, if its operands are then on wires (checked on
samples by tools/r46_dqi_scan_heldchk.py). With P1's exact zero intervals the tool gives, stage by stage, the ops at
tight live 902 / 901 / 900 before and after parking each class.
Usage: python r46_dqi_scan_held.py scan <dqi.kmb> <gaps.npz> [G=500]
       python r46_dqi_scan_held.py rep <dqi.kmb> <gaps.npz> <zero.npz> [min gap=500] [list=40]
"""
import os
import sys
import time

import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

DATA = 1000000


@njit(cache=True)
def scan(ops, q3a, q3b, q3c, q2a, q2b, q1, nq, G, out):
    last = np.full(nq, -1, np.int64)
    nwr = np.zeros(nq, np.int64)
    bop = np.full(nq, -1, np.int64)
    bk = np.zeros(nq, np.int64)
    ba = np.full(nq, -1, np.int64)
    bb = np.full(nq, -1, np.int64)
    bna = np.zeros(nq, np.int64)
    bnb = np.zeros(nq, np.int64)
    for w in range(840):
        nwr[w] = DATA
    i3 = 0
    i2 = 0
    i1 = 0
    ng = 0
    cap = out.shape[0]
    tw = np.zeros(3, np.int64)
    for oi in range(ops.shape[0]):
        op = ops[oi]
        if op == 0x81:
            continue
        base = op & 0x7F
        if base == 0x01:
            continue
        nqk = (op & 0x30) >> 4
        a = -1
        b = -1
        t = -1
        nt = 0
        if nqk == 3:
            a = q3a[i3]
            b = q3b[i3]
            t = q3c[i3]
            i3 += 1
            tw[0] = a
            tw[1] = b
            tw[2] = t
            nt = 3
        elif nqk == 2:
            a = q2a[i2]
            b = q2b[i2]
            i2 += 1
            tw[0] = a
            tw[1] = b
            nt = 2
        elif nqk == 1:
            a = q1[i1]
            i1 += 1
            tw[0] = a
            nt = 1
        else:
            continue
        for x in range(nt):
            w = tw[x]
            if last[w] >= 0 and oi - last[w] >= G and nwr[w] >= 1 and nwr[w] <= 64:
                if ng < cap:
                    out[ng, 0] = w
                    out[ng, 1] = last[w]
                    out[ng, 2] = oi
                    out[ng, 3] = nwr[w]
                    out[ng, 4] = bop[w]
                    out[ng, 5] = bk[w]
                    out[ng, 6] = ba[w]
                    out[ng, 7] = bb[w]
                    out[ng, 8] = bna[w]
                    out[ng, 9] = bnb[w]
                    out[ng, 10] = base
                    out[ng, 11] = x
                ng += 1
            last[w] = oi
        # writes
        if base == 0x12 or base == 0x50:
            nwr[a] = 0
            bop[a] = -1
        elif base == 0x11:
            if nwr[a] == 0:
                bop[a] = oi
                bk[a] = base
                ba[a] = -1
                bb[a] = -1
            nwr[a] += 1
        elif base == 0x21:
            if nwr[b] == 0:
                bop[b] = oi
                bk[b] = base
                ba[b] = a
                bb[b] = -1
                bna[b] = nwr[a]
                bnb[b] = 0
            nwr[b] += 1
        elif base == 0x22:
            x1 = nwr[a]
            nwr[a] = nwr[b]
            nwr[b] = x1
            x1 = bop[a]
            bop[a] = bop[b]
            bop[b] = x1
            x1 = bk[a]
            bk[a] = bk[b]
            bk[b] = x1
            x1 = ba[a]
            ba[a] = ba[b]
            ba[b] = x1
            x1 = bb[a]
            bb[a] = bb[b]
            bb[b] = x1
            x1 = bna[a]
            bna[a] = bna[b]
            bna[b] = x1
            x1 = bnb[a]
            bnb[a] = bnb[b]
            bnb[b] = x1
        elif base == 0x31:
            if nwr[t] == 0:
                bop[t] = oi
                bk[t] = base
                ba[t] = a
                bb[t] = b
                bna[t] = nwr[a]
                bnb[t] = nwr[b]
            nwr[t] += 1
        elif base == 0x32:
            nwr[b] += DATA
            nwr[t] += DATA
    return ng


def stage_ranges(c):
    out = []
    depth = 0
    start = name = None
    for p, t in zip(c.tag_positions().tolist(), c.tags):
        if t.startswith("@enter:"):
            if depth == 0:
                start, name = p, t[7:]
            depth += 1
        elif t.startswith("@exit:"):
            depth -= 1
            if depth == 0:
                out.append((name, start, p))
    return out


def main():
    t0 = time.time()
    mode = sys.argv[1]
    c = kmb.read(sys.argv[2])
    n = len(c.ops)
    nq = c.num_qubits
    if mode == "scan":
        G = int(sys.argv[4]) if len(sys.argv) > 4 else 500
        ops = np.ascontiguousarray(np.asarray(c.ops))
        arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64) & kmb.ID_MASK) for x in (c.q3[0], c.q3[1], c.q3[2],
                                                                                         c.q2[0], c.q2[1], c.q1)]
        out = np.zeros((6000000, 12), np.int64)
        ng = scan(ops, *arr, nq, G, out)
        print("gaps >= %d ops: %d (kept %d) (%.0fs)" % (G, ng, min(ng, len(out)), time.time() - t0))
        np.savez(sys.argv[3], gaps=out[:min(ng, len(out))], G=G)
        return
    Z = np.load(sys.argv[4])
    g = np.load(sys.argv[3])["gaps"]
    ming = int(sys.argv[5]) if len(sys.argv) > 5 else 500
    nlist = int(sys.argv[6]) if len(sys.argv) > 6 else 40
    zw = Z["wire"].astype(np.int64)
    s = np.clip(Z["start"].astype(np.int64), 0, n)
    e = np.clip(Z["end"].astype(np.int64), 0, n)
    ok = e > s
    zw, s, e = zw[ok], s[ok], e[ok]
    d = np.zeros(n + 2, np.int32)
    np.add.at(d, s, 1)
    np.add.at(d, e, -1)
    live = nq - np.cumsum(d[:n])
    # a gap is HELD when the wire is not exactly zero at its middle
    order = np.lexsort((s, zw))
    zw, s, e = zw[order], s[order], e[order]
    first = np.searchsorted(zw, np.arange(nq + 1))
    w, gs, ge, nwr, bop, bk, ba, bb, bna, bnb, ebase, erole = [g[:, k] for k in range(12)]
    ln = ge - gs
    keep = ln >= ming
    mid = (gs + ge) // 2
    held = np.zeros(len(g), bool)
    for i in np.flatnonzero(keep):
        a, b = first[w[i]], first[w[i] + 1]
        k = np.searchsorted(s[a:b], mid[i], side="right") - 1
        held[i] = not (k >= 0 and e[a + k] > mid[i])
    keep &= held
    one_and = keep & (nwr == 1) & (bk == 0x31)
    and_data = one_and & (bna >= DATA) & (bnb >= DATA)
    and_other = one_and & ~and_data
    few = keep & (nwr >= 1) & (nwr <= 4) & ~one_and
    copy1 = keep & (nwr == 1) & (bk == 0x21)
    print("gaps >= %d ops with the wire not exactly zero: %d; made by ONE CCX: %d (both operands data words: %d, "
          "other: %d); other wires with <= 4 writes since reset: %d (one CX: %d)" % (
              ming, int(keep.sum()), int(one_and.sum()), int(and_data.sum()), int(and_other.sum()), int(few.sum()),
              int(copy1.sum())))
    classes = (("one AND of two data wires", and_data), ("one AND, an operand not data (ladder end / flag)", and_other),
               ("<= 4 writes, not one AND", few))
    cover = {}
    for nm, m in classes:
        dd = np.zeros(n + 2, np.int32)
        np.add.at(dd, gs[m] + 1, 1)
        np.add.at(dd, ge[m], -1)
        cover[nm] = np.cumsum(dd[:n])
    print("STAGES: ops at tight live 902 / 901 / 900 today | after parking the idle single ANDs of data | after parking "
          "every idle single AND | idle wire-ops of each class in the stage")
    for snm, a, b in stage_ranges(c):
        if snm in ("mask_projection", "restore_masks", "initial", "uninitial"):
            continue
        L = live[a:b]
        c1 = cover[classes[0][0]][a:b]
        c2 = c1 + cover[classes[1][0]][a:b]
        row = []
        for LL in (L, L - c1, L - c2):
            row.append("max %d: %d / %d / %d" % (int(LL.max()), int((LL == 902).sum()), int((LL == 901).sum()),
                                                 int((LL == 900).sum())))
        print("%-18s %s | %s | %s | wire-ops %d / %d / %d" % (
            snm, row[0], row[1], row[2], int(c1.sum()), int(cover[classes[1][0]][a:b].sum()),
            int(cover[classes[2][0]][a:b].sum())))
    names = {0x31: "CCX", 0x21: "CX", 0x11: "X", 0x22: "SWAP", 0x32: "CSWAP", 0x50: "HMR", 0x12: "R", 0x20: "CZ",
             0x30: "CCZ", 0x10: "Z"}
    st = stage_ranges(c)

    def stage_of(x):
        for nm, a, b in st:
            if a <= x < b:
                return nm
        return "-"
    for nm, m in classes[:2]:
        idx = np.flatnonzero(m)
        pk = np.array([int((live[gs[i] + 1:ge[i]] >= 902).sum()) for i in idx])
        pk1 = np.array([int((live[gs[i] + 1:ge[i]] >= 901).sum()) for i in idx])
        print("CLASS %s: %d gaps, %d wire-ops; gaps with ops at 902 inside: %d (902-ops covered, with multiplicity: %d); "
              "at >= 901: %d" % (nm, len(idx), int(ln[idx].sum()), int((pk > 0).sum()), int(pk.sum()),
                                 int((pk1 > 0).sum())))
        o = np.argsort(-pk * 1000000 - ln[idx])[:nlist]
        for j in o:
            i = idx[j]
            print("   wire %3d made at op %9d by CCX %d %d (writes since reset of the operands %s, %s); idle %9d..%9d "
                  "(%7d ops, %7d at 902, %7d at >= 901) %s; closed by %s (operand %d)" % (
                      w[i], bop[i], ba[i], bb[i], "data" if bna[i] >= DATA else str(bna[i]),
                      "data" if bnb[i] >= DATA else str(bnb[i]), gs[i], ge[i], ln[i], pk[j], pk1[j],
                      stage_of(int(gs[i])), names.get(int(ebase[i]), hex(int(ebase[i]))), erole[i]))


if __name__ == "__main__":
    main()
