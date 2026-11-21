#!/usr/bin/env python3
"""dqi_int_runwho.py -- B5 (16th holder, Nov 16): for named runs of a decoder file, the IDLE wires (no touch inside the
run, live through it) and how each could be remade at its next touch. Read only.

Per idle wire w of run [ra, rb): lt = last touch before ra, nt = next touch at or after rb, birth = the op that ends
w's last exact zero interval before ra (P1's zero.npz), and the WRITES to w between birth and lt (CCX / CSWAP / CX /
SWAP / X / R targets). Classes:
  copy  : one write, CX(a -> w), a not written in (birth, nt)           -> park cost 0 (CX at the open and the close)
  and1  : one write, CCX(a, b -> w), a, b not written in (birth, nt)     -> park cost 1 (HMR, CCX, Z_IF)
  andx  : one write, CCX with a control written before nt
  flag  : birth key holds zero_toggle / flag (a zero test; see r48_dqi_ret2_rp_flagchk.py)
  other : anything else (count of writes printed)
Usage: python dqi_int_runwho.py <dqi.kmb> <zero.npz> <segments_real.json|-> <ra:rb,ra:rb,...> [maxlist=40]
"""
import collections
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

SH = 27


def dd(names):
    out = []
    for x in names:
        if not out or out[-1] != x:
            out.append(x)
    return out


def main():
    t0 = time.time()
    path, zp, segp, runs = sys.argv[1:5]
    maxlist = int(sys.argv[5]) if len(sys.argv) > 5 else 40
    c = kmb.read(path)
    ops = np.asarray(c.ops)
    n = len(ops)
    nq = c.num_qubits
    M = kmb.ID_MASK
    base = ops & 0x7F
    k = ((ops & 0x30) >> 4).astype(np.int8)
    i1 = np.cumsum(k == 1, dtype=np.int64) - 1
    i2 = np.cumsum(k == 2, dtype=np.int64) - 1
    i3 = np.cumsum(k == 3, dtype=np.int64) - 1
    q1 = np.asarray(c.q1).astype(np.int64) & M
    q2 = [np.asarray(x).astype(np.int64) & M for x in c.q2]
    q3 = [np.asarray(x).astype(np.int64) & M for x in c.q3]
    Z = np.load(zp)
    zw = Z['wire'].astype(np.int64)
    zs = np.clip(Z['start'].astype(np.int64), 0, n)
    ze = np.clip(Z['end'].astype(np.int64), 0, n)
    ok = ze > zs
    zw, zs, ze = zw[ok], zs[ok], ze[ok]
    o = np.lexsort((ze, zw))
    zw, zs, ze = zw[o], zs[o], ze[o]
    zke = (zw << SH) | ze
    real = (ops != 0x81) & (base != 0x17) & (base != 0x07)
    # touches and writes
    tparts, wparts = [], []
    for cls, arrs in ((1, (q1,)), (2, q2), (3, q3)):
        idx = np.flatnonzero((k == cls) & real)
        ci = {1: i1, 2: i2, 3: i3}[cls][idx]
        for slot, a in enumerate(arrs):
            w = a[ci]
            tparts.append((w << SH) | idx)
            b = base[idx]
            if cls == 3:
                wr = ((b == 0x31) & (slot == 2)) | ((b == 0x32) & (slot >= 1))
            elif cls == 2:
                wr = ((b == 0x21) & (slot == 1)) | (b == 0x22)
            else:
                wr = (b == 0x11) | (b == 0x12) | (b == 0x50)
            wparts.append((w[wr] << SH) | idx[wr])
    tk = np.sort(np.concatenate(tparts))
    wk = np.sort(np.concatenate(wparts))
    del tparts, wparts
    keyof = lambda i: '-'  # noqa: E731
    if segp != '-':
        S = json.load(open(segp))
        stacks, segs = S['stacks'], S['segments']
        mark = np.flatnonzero(base == 0x07)

        def keyof(i, depth=3):
            j = int(np.searchsorted(mark, i, 'right')) - 1
            if j < 0 or j >= len(segs):
                return '-'
            return ' > '.join(dd([x.split('@')[0] for x in stacks[segs[j][0]]])[-depth:])
    print('%s: Q %d, ops %d (%.0fs)' % (path, nq, n, time.time() - t0), flush=True)

    def writes_between(w, a, b):
        lo = np.searchsorted(wk, (w << SH) | (a + 1), 'left')
        hi = np.searchsorted(wk, (w << SH) | b, 'left')
        return (wk[lo:hi] & ((1 << SH) - 1))

    W = np.arange(nq, dtype=np.int64)
    for rr in runs.split(','):
        ra, rb = (int(x) for x in rr.split(':'))
        p = np.searchsorted(tk, (W << SH) | ra, 'left')
        pc = np.minimum(p, len(tk) - 1)
        same = (p < len(tk)) & ((tk[pc] >> SH) == W)
        nxt = np.where(same, tk[pc] & ((1 << SH) - 1), n)
        idle = nxt >= rb
        q = np.searchsorted(zke, (W << SH) | ra, 'right')
        qc = np.minimum(q, len(zke) - 1)
        hit = (q < len(zke)) & (zw[qc] == W) & (zs[qc] < rb)
        idle &= ~hit
        iw = np.flatnonzero(idle)
        cnt = collections.Counter()
        rows = []
        for w in iw.tolist():
            pl = p[w] - 1
            lt = int(tk[pl] & ((1 << SH) - 1)) if pl >= 0 and (tk[pl] >> SH) == w else -1
            nt = int(nxt[w])
            qb = np.searchsorted(zke, (w << SH) | ra, 'right') - 1
            birth = int(ze[qb]) if qb >= 0 and zw[qb] == w else -1
            wr = writes_between(w, birth - 1, lt + 1) if birth >= 0 else np.zeros(0, np.int64)
            bk = keyof(birth) if birth >= 0 else 'input'
            cls = 'other'
            det = '%d writes' % len(wr)
            if 'zero_toggle' in bk or 'flag' in bk:
                cls = 'flag'
            if len(wr) == 1:
                j = int(wr[0])
                if base[j] == 0x21 and k[j] == 2:
                    a = int(q2[0][i2[j]])
                    aw = writes_between(a, j, nt)
                    cls = 'copy' if len(aw) == 0 else 'copyx'
                    det = 'CX %d -> %d at %d; source writes before close %d' % (a, w, j, len(aw))
                elif base[j] == 0x31 and k[j] == 3:
                    a, b_ = int(q3[0][i3[j]]), int(q3[1][i3[j]])
                    na, nb = len(writes_between(a, j, nt)), len(writes_between(b_, j, nt))
                    cls = 'and1' if na + nb == 0 else 'andx'
                    det = 'CCX %d,%d -> %d at %d [%s]; control writes before close %d / %d' % (a, b_, w, j, keyof(j), na, nb)
            cnt[cls] += 1
            rows.append((cls, w, birth, lt, nt, bk, det))
        print('\nRUN %d..%d (%d ops) [%s]: idle %d; classes %s' % (ra, rb, rb - ra, keyof((ra + rb) // 2), len(iw),
                                                                   dict(cnt)), flush=True)
        order = {'copy': 0, 'and1': 1, 'flag': 2, 'copyx': 3, 'andx': 4, 'other': 5}
        rows.sort(key=lambda r: (order[r[0]], r[1]))
        for (cls, w, birth, lt, nt, bk, det) in rows[:maxlist]:
            print('  %-5s w %4d birth %d [%s] lt %d nt %d [%s] | %s' % (cls, w, birth, bk, lt, nt, keyof(nt), det))
        oth = collections.Counter(r[5] for r in rows if r[0] == 'other')
        print('  other by birth key: %s' % '; '.join('%s x%d' % kv for kv in oth.most_common(8)))
    print('runtime %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
