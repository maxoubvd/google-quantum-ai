#!/usr/bin/env python3
"""dqi_pp_zero.py -- workstream dqi_pp (P1): EXACT idle-zero intervals of every wire (no sampling).

A wire is free between two consecutive events on it when its exact form after the first event is 0 (value numbering
of tools/dqi_pp_gvn.py: constant 0, no atom): this covers resets and every gate-level uncompute that GVN can follow.
Untouched non-input wires are free from the start. Output in the format of tools/dqi_int_ppzero.py (wire, start, end,
structural), read by tools/dqi_int_ppret2.py: free at after-op times start .. end-1.
Usage: python dqi_pp_zero.py <circuit.kmb> <out.npz>
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_pp_gvn as g  # noqa: E402
import dqi_pp_prove as P  # noqa: E402


def main():
    t0 = time.time()
    c = kmb.read(sys.argv[1])
    out = sys.argv[2]
    gc, ops, arr = g.load(sys.argv[1])
    ntoff = int((((ops & 0x7C) == 0x30) & (ops != 0x81)).sum())
    zf = np.zeros(len(ops), np.uint8)
    kz = np.zeros(len(ops), np.uint8)
    g.run(ops, *arr, c.num_qubits, c.num_bits, 1, ntoff, zf, kz)
    print('ops that never execute (condition identically 0): %d' % int(kz.sum()), flush=True)
    print('value numbering with zero flags %.1fs' % (time.time() - t0), flush=True)
    _, W0, W1, W2, BT, CB, DEP = P.decode(c)
    n = len(ops)
    base = ops & 0x7F
    Wm = np.stack([W0, W1, W2], axis=1)
    ev_op = np.repeat(np.arange(n, dtype=np.int64), 3)
    ev_w = Wm.reshape(-1).astype(np.int64)
    ev_slot = np.tile(np.arange(3), n)
    kk = (ev_w >= 0) & (kz[ev_op] == 0)                # events of never-executed ops touch nothing
    ev_op, ev_w, ev_slot = ev_op[kk], ev_w[kk], ev_slot[kk]
    zero_after = ((zf[ev_op] >> ev_slot.astype(np.uint8)) & 1).astype(bool)
    o = np.lexsort((ev_op, ev_w))
    ev_op, ev_w, zero_after = ev_op[o], ev_w[o], zero_after[o]
    ne = len(ev_op)
    nxt = np.full(ne, n, np.int64)
    same = ev_w[1:] == ev_w[:-1]
    nxt[:-1][same] = ev_op[1:][same]
    sel = np.flatnonzero(zero_after)
    wire, start, end = ev_w[sel], ev_op[sel], nxt[sel]
    structural = np.isin(base[start], (0x50, 0x12)) & (CB[start] < 0) & (DEP[start] == 0)
    # untouched non-input wires from the start
    first = np.ones(ne, bool)
    first[1:] = ~same
    fw, fo = ev_w[first], ev_op[first]
    m = fw >= 840
    touched = np.zeros(c.num_qubits, bool)
    touched[fw] = True
    never = np.flatnonzero(~touched & (np.arange(c.num_qubits) >= 840))
    wire = np.concatenate([wire, fw[m], never])
    start = np.concatenate([start, np.full(int(m.sum()), -1), np.full(len(never), -1)])
    end = np.concatenate([end, fo[m], np.full(len(never), n)])
    structural = np.concatenate([structural, np.ones(int(m.sum()) + len(never), bool)])
    keep = end > start + 1
    wire, start, end, structural = wire[keep], start[keep], end[keep], structural[keep]
    print('exact idle-zero intervals %d (structural %d), wire-ops %d (structural %d)' % (
        len(wire), int(structural.sum()), int((end - start).sum()), int((end - start)[structural].sum())))
    np.savez(out, wire=wire, start=start, end=end, structural=structural)
    print('done %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
