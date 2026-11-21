#!/usr/bin/env python3
"""dqi_int_kmbw.py -- workstream dqi_int (B5): an independent .kmb writer for tools/kmb.py's reader, and an op-level
editor (delete ops, insert ops after given ops), per research/05_kickmix_and_techniques.md section A3.

  write(path, k)            k: a kmb.Kmb (as read by tools/kmb.py, or edited); header counts recomputed.
  edit(k, delete, inserts)  delete: op indices to drop; inserts: {op index: [(opcode, qubits tuple, bit_t, cond)]}
                            new ops placed right AFTER the given original op (after it is deleted, if it is).
Self-test: python dqi_int_kmbw.py <in.kmb> <out.kmb>  (read, write, compare bytes).
"""
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402


def _counts(k):
    ops = np.asarray(k.ops)
    nq = (ops & 0x30) >> 4
    n3, n2, n1 = int((nq == 3).sum()), int((nq == 2).sum()), int((nq == 1).sum())
    nbit = len(k.bit_t)
    ncond = len(k.cond)
    nang = len(k.angles)
    ntag = len(k.tags)
    nchar = sum(len(t.encode('utf-8')) for t in k.tags)
    return n3, n2, n1, nbit, ncond, nang, ntag, nchar


def write(path, k):
    ops = np.ascontiguousarray(np.asarray(k.ops, dtype=np.uint8))
    n3, n2, n1, nbit, ncond, nang, ntag, nchar = _counts(k)
    q3 = [np.asarray(x, dtype='<u4') for x in k.q3]
    q2 = [np.asarray(x, dtype='<u4') for x in k.q2]
    q1 = np.asarray(k.q1, dtype='<u4')
    assert len(q3[0]) == n3 and len(q2[0]) == n2 and len(q1) == n1, 'operand arrays do not match the opcodes'
    nqub = 0
    for arr in q3 + q2 + [q1]:
        if len(arr):
            nqub = max(nqub, int(arr.max()) + 1)
    for name, ids, isq in k.registers:
        if len(ids) and isq.any():
            nqub = max(nqub, int(np.asarray(ids)[isq].max()) + 1)
    nbits = k.num_bits                                   # bits are never removed by our edits
    with open(path, 'wb') as f:
        f.write(kmb.MAGIC)
        f.write(struct.pack('<IIII', 2, nqub, nbits, len(k.registers)))
        f.write(struct.pack('<9Q', len(ops), n3, n2, n1, nbit, ncond, nang, ntag, nchar))

        def packet(t, payload):
            f.write(struct.pack('<QQ', t, len(payload)))
            f.write(payload)
        for name, ids, isq in k.registers:
            nb = name.encode('utf-8')
            raw = (np.asarray(ids, dtype=np.uint64) | (np.asarray(isq, dtype=np.uint64) << np.uint64(31)))
            packet(0, struct.pack('<II', len(nb), len(raw)) + nb + raw.astype('<u4').tobytes())
        packet(1, ops.tobytes())
        packet(2, q3[0].tobytes())
        packet(3, q3[1].tobytes())
        packet(4, q3[2].tobytes())
        packet(5, q2[0].tobytes())
        packet(6, q2[1].tobytes())
        packet(7, q1.tobytes())
        packet(8, np.asarray(k.bit_t, dtype='<u4').tobytes())
        packet(9, np.asarray(k.cond, dtype='<u4').tobytes())
        packet(10, np.asarray(k.angles, dtype='<u8').reshape(-1).tobytes())
        data = [t.encode('utf-8') for t in k.tags]
        ends = np.cumsum([len(x) for x in data], dtype=np.uint64) if data else np.zeros(0, np.uint64)
        packet(11, struct.pack('<Q', len(data)) + ends.astype('<u8').tobytes() + b''.join(data))
    return nqub


def _class_index(ops):
    """For each op, the index of its element in each array class (valid only where the op uses that class)."""
    nq = (ops & 0x30) >> 4
    idx = {}
    for k in (1, 2, 3):
        idx[k] = np.cumsum(nq == k) - 1
    hb = ((ops & 0x40) != 0) & ((ops & 0x7F) != 0x06) & ((ops & 0x7F) != 0x07)
    idx['b'] = np.cumsum(hb) - 1
    hc = (ops & 0x80) != 0
    idx['c'] = np.cumsum(hc) - 1
    return nq, hb, hc, idx


def edit(k, delete, inserts):
    """Returns a new kmb.Kmb with the ops in `delete` removed and `inserts` added (see module docstring)."""
    ops = np.asarray(k.ops, dtype=np.uint8)
    n = len(ops)
    nq, hb, hc, idx = _class_index(ops)
    keep = np.ones(n, dtype=bool)
    keep[np.asarray(sorted(delete), dtype=np.int64)] = False
    # per-op operand tuples only where needed: build new arrays by concatenating kept slices and inserted elements
    cols3 = [np.asarray(x) for x in k.q3]
    cols2 = [np.asarray(x) for x in k.q2]
    col1 = np.asarray(k.q1)
    colb = np.asarray(k.bit_t)
    colc = np.asarray(k.cond)
    pts = sorted(inserts)
    assert not pts or pts[-1] < n
    new_ops, n3, n2, n1, nb, nc = [], [[], [], []], [[], []], [], [], []
    prev = 0
    for p in pts + ([n - 1] if not pts or pts[-1] != n - 1 else []):
        hi = p + 1
        seg = slice(prev, hi)
        kseg = keep[seg]
        o = ops[seg][kseg]
        new_ops.append(o)
        sel = np.arange(prev, hi)[kseg]
        m3 = nq[sel] == 3
        m2 = nq[sel] == 2
        m1 = nq[sel] == 1
        for j in range(3):
            n3[j].append(cols3[j][idx[3][sel[m3]]])
        for j in range(2):
            n2[j].append(cols2[j][idx[2][sel[m2]]])
        n1.append(col1[idx[1][sel[m1]]])
        nb.append(colb[idx['b'][sel[hb[sel]]]])
        nc.append(colc[idx['c'][sel[hc[sel]]]])
        prev = hi
        if p in inserts and p < n:
            for (code, qs, bt, cd) in inserts[p]:
                new_ops.append(np.array([code], dtype=np.uint8))
                kq = (code & 0x30) >> 4
                assert len(qs) == kq, (code, qs)
                if kq == 3:
                    for j in range(3):
                        n3[j].append(np.array([qs[j]], dtype=np.uint32))
                elif kq == 2:
                    for j in range(2):
                        n2[j].append(np.array([qs[j]], dtype=np.uint32))
                elif kq == 1:
                    n1.append(np.array([qs[0]], dtype=np.uint32))
                if (code & 0x40) and (code & 0x7F) not in (0x06, 0x07):
                    nb.append(np.array([bt], dtype=np.uint32))
                if code & 0x80:
                    nc.append(np.array([cd], dtype=np.uint32))
    cat = lambda L, dt: np.concatenate(L).astype(dt) if L else np.zeros(0, dt)  # noqa: E731
    out = kmb.Kmb(k.num_qubits, k.num_bits, k.registers, cat(new_ops, np.uint8),
                  tuple(cat(n3[j], np.uint32) for j in range(3)), tuple(cat(n2[j], np.uint32) for j in range(2)),
                  cat(n1, np.uint32), cat(nb, np.uint32), cat(nc, np.uint32), np.asarray(k.angles), list(k.tags))
    return out


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    k = kmb.read(src)
    write(dst, k)
    a = open(src, 'rb').read()
    b = open(dst, 'rb').read()
    print('round trip identical:', a == b, len(a), len(b))
