"""dqi_int_idle901.py -- B5 (17th holder, Nov 16): a COPY of the RET2 re-pricer's r48_dqi_ret2_rp_idle.py (credit:
R48-DQI-RET2 re-pricer; the original is unchanged) with ONE switch: B5_IDLE_EXCLUDE=<comma list of wire ids> drops
those wires from the flag candidates (a cover flag whose re-test at the close is not exact), so that the greedy cover
picks another idle flag for its runs.  Everything else is the original.
"""
#!/usr/bin/env python3
"""r48_dqi_ret2_rp_idle.py -- R48-DQI-RET2 re-pricer: IDLE HOLDERS of the runs at tight live >= LEVEL. Read only.

For one decoder file and its exact zero intervals (tools/dqi_pp_zero.py; live(t) = Q - wires exactly zero after op t,
as tools/r48_dqi_ret2_occ.py), in the named stages: every maximal run of ops at live >= LEVEL. Per run, the wires that
are live through the whole run and that no op inside it touches (marker ops DEBUG_PRINT_Q excluded): its IDLE
HOLDERS. Per idle wire: its gap (last touch before the run .. next touch after it), its BIRTH (the op that ends its
last zero interval before the run; 'input' if it was never zero before) with the marker key there, and its DEATH (the
op that starts its next zero interval). A wire is a FLAG candidate when its birth key holds 'zero_toggle' or
'flag' and no other idle wire of the run has the same birth op within +-64 ops (a one-wire value).
Cover: candidate parks = (flag wire, its gap); a park covers every run that lies inside its gap; greedy (most runs
first) until every run with an idle flag is covered. A flag park costs one zero test of its word at the close (11
ANDs for a 12-bit word) IF that word is intact there: NOT checked here (an estimate for the builders).
Usage: python r48_dqi_ret2_rp_idle.py <dqi.kmb> <zero.npz> <stages comma> [LEVEL=902] [segments_real.json|-] [maxlist=60]
"""
import collections
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_gvn import stages  # noqa: E402

SH = 27


def dd(names):
    out = []
    for x in names:
        if not out or out[-1] != x:
            out.append(x)
    return out


_B5_EXCL = set(int(x) for x in os.environ.get('B5_IDLE_EXCLUDE', '').split(',') if x.strip())

def main():
    t0 = time.time()
    path, zp, stl = sys.argv[1:4]
    LEVEL = int(sys.argv[4]) if len(sys.argv) > 4 else 902
    segp = sys.argv[5] if len(sys.argv) > 5 else '-'
    maxlist = int(sys.argv[6]) if len(sys.argv) > 6 else 60
    want = stl.split(',')
    c = kmb.read(path)
    ops = np.asarray(c.ops)
    n = len(ops)
    nq = c.num_qubits
    assert n < (1 << SH)
    M = kmb.ID_MASK
    # live count, as r48_dqi_ret2_occ.py
    Z = np.load(zp)
    zw = Z['wire'].astype(np.int64)
    zs = np.clip(Z['start'].astype(np.int64), 0, n)
    ze = np.clip(Z['end'].astype(np.int64), 0, n)
    ok = ze > zs
    zw, zs, ze = zw[ok], zs[ok], ze[ok]
    d = np.zeros(n + 2, np.int32)
    np.add.at(d, zs, 1)
    np.add.at(d, ze, -1)
    live = nq - np.cumsum(d[:n])
    o = np.lexsort((ze, zw))          # by wire, then end (intervals of one wire may overlap before P1's merge)
    zw, zs, ze = zw[o], zs[o], ze[o]
    zke = (zw << SH) | ze
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    # touch index: key = wire << SH | op, every qubit operand of every op except the marker DEBUG_PRINT_Q (0x17)
    k = ((ops & 0x30) >> 4).astype(np.int8)
    base = ops & 0x7F
    real = (ops != 0x81) & (base != 0x17)
    parts = []
    for cls, arrs in ((1, (c.q1,)), (2, c.q2), (3, c.q3)):
        sel = k == cls
        idx = np.flatnonzero(sel)
        keep = real[idx]
        for a in arrs:
            w = (np.asarray(a).astype(np.int64) & M)[:len(idx)]
            parts.append((w[keep] << SH) | idx[keep])
    tk = np.sort(np.concatenate(parts))
    del parts
    print('%s: Q %d, ops %d, touches %d (%.0fs)' % (path, nq, n, len(tk), time.time() - t0), flush=True)
    # marker keys
    keyof = None
    if segp != '-':
        S = json.load(open(segp))
        stacks, segs = S['stacks'], S['segments']
        mark = np.flatnonzero(base == 0x07)

        def keyof(i, depth=3):
            j = int(np.searchsorted(mark, i, 'right')) - 1
            if j < 0 or j >= len(segs):
                return '-'
            return ' > '.join(dd([x.split('@')[0] for x in stacks[segs[j][0]]])[-depth:])
    W = np.arange(nq, dtype=np.int64)
    names = {0x31: 'CCX', 0x32: 'CSWAP', 0x30: 'CCZ', 0x21: 'CX', 0x22: 'SWAP', 0x20: 'CZ', 0x50: 'HMR', 0x12: 'R',
             0x11: 'X'}
    allruns = []      # (stage, ra, rb)
    cand = []         # (wire, lt, nt, birth, birthkey)
    per_run = []
    for nm, a, b in stages(c):
        if nm not in want:
            continue
        hi = (live[a:b] >= LEVEL).astype(np.int8)
        dh = np.diff(np.concatenate(([0], hi, [0])))
        starts = np.flatnonzero(dh == 1) + a
        ends = np.flatnonzero(dh == -1) + a
        nflag_runs = 0
        for ra, rb in zip(starts.tolist(), ends.tolist()):
            ri = len(allruns)
            allruns.append((nm, ra, rb))
            # untouched in [ra, rb)
            p = np.searchsorted(tk, (W << SH) | ra, 'left')
            pc = np.minimum(p, len(tk) - 1)
            same = (p < len(tk)) & ((tk[pc] >> SH) == W)
            nxt = np.where(same, tk[pc] & ((1 << SH) - 1), n)
            idle = nxt >= rb
            # live through [ra, rb): at LEVEL = Q every wire is nonzero at every op of the run; below Q an
            # approximate test (first interval of w ending after ra must start at or after rb)
            if LEVEL < nq:
                q = np.searchsorted(zke, (W << SH) | ra, 'right')
                qc = np.minimum(q, len(zke) - 1)
                hit = (q < len(zke)) & ((zw[qc]) == W) & (zs[qc] < rb)
                idle &= ~hit
            iw = np.flatnonzero(idle)
            # gap ends and birth / death
            pl = p[iw] - 1
            plc = np.maximum(pl, 0)
            lt = np.where((pl >= 0) & ((tk[plc] >> SH) == iw), tk[plc] & ((1 << SH) - 1), -1)
            nt = nxt[iw]
            qb = np.searchsorted(zke, (iw << SH) | ra, 'right') - 1   # last interval with end <= ra
            qbc = np.maximum(qb, 0)
            birth = np.where((qb >= 0) & (zw[qbc] == iw), ze[qbc], -1)
            # flag candidates: birth key with zero_toggle / flag, one wire per birth (no other idle wire born +-64)
            bo = np.sort(birth[birth >= 0])
            nf = 0
            fl = []
            cls = collections.Counter()
            for j, w in enumerate(iw.tolist()):
                bj = int(birth[j])
                bk = 'input' if bj < 0 else (keyof(bj) if keyof else names.get(int(base[bj]), hex(int(base[bj]))))
                cls[bk] += 1
                if bj < 0:
                    continue
                lo = np.searchsorted(bo, bj - 64, 'left')
                hi2 = np.searchsorted(bo, bj + 64, 'right')
                single = (hi2 - lo) == 1
                if single and ('zero_toggle' in bk or 'flag' in bk) and w not in _B5_EXCL:
                    nf += 1
                    fl.append((w, bk))
                    cand.append((w, int(lt[j]), int(nt[j]), bj, bk, ri))
            nflag_runs += nf > 0
            midkey = keyof((ra + rb) // 2) if keyof else '-'
            per_run.append((nm, ra, rb, midkey, len(iw), nf, fl, cls, int(toff[ra:rb].sum())))
        print('stage %-20s runs at >= %d: %d; runs with an idle one-wire flag: %d' % (nm, LEVEL, len(starts),
                                                                                    nflag_runs), flush=True)
    # greedy cover of the runs by flag gaps
    R = len(allruns)
    covered = np.zeros(R, bool)
    ra_arr = np.array([r[1] for r in allruns], np.int64)
    rb_arr = np.array([r[2] for r in allruns], np.int64)
    gaps = {}
    for (w, lt, nt, bj, bk, ri) in cand:
        gaps.setdefault((w, lt, nt), (bk, bj))
    cov = {g: np.flatnonzero((ra_arr > g[1]) & (rb_arr <= g[2])) for g in gaps}
    chosen = []
    while True:
        best, bc = None, 0
        for g, rr in cov.items():
            m = int((~covered[rr]).sum())
            if m > bc:
                best, bc = g, m
        if best is None:
            break
        chosen.append((best, bc))
        covered[cov[best]] = True
    print('\nCOVER (one-wire flag parks, greedy): %d parks cover %d of %d runs; %d runs have no idle one-wire flag'
          % (len(chosen), int(covered.sum()), R, int((~covered).sum())))
    for (w, lt, nt), m in chosen:
        bk, bj = gaps[(w, lt, nt)]
        print('  park wire %4d  gap %d..%d (%d ops)  covers %d runs  born %d in [%s]  closes at [%s]' % (
            w, lt, nt, nt - lt, m, bj, bk, keyof(nt) if keyof else names.get(int(base[min(nt, n - 1)]), '?')))
    tc = sum(r[8] for i, r in enumerate(per_run) if covered[i])
    tu = sum(r[8] for i, r in enumerate(per_run) if not covered[i])
    print('Toffoli-class ops at >= %d: in covered runs %d, in runs without an idle flag %d' % (LEVEL, tc, tu))
    byst = collections.Counter((r[0], bool(covered[i])) for i, r in enumerate(per_run))
    print('runs by stage (covered / not): %s' % ', '.join('%s %d / %d' % (s, byst[(s, True)], byst[(s, False)])
                                                         for s in want))
    print('\nRUNS (stage, ops, C = covered, key at the middle, Toffoli ops, idle wires, idle one-wire flags, '
          'top birth classes of idle wires)')
    for i, (nm, ra, rb, mk, ni, nf, fl, cls, nt_) in enumerate(per_run[:maxlist]):
        top = '; '.join('%s x%d' % kv for kv in cls.most_common(4))
        print('%-10s %9d..%-9d %6d ops %s [%s] T %d idle %d flags %d %s | %s' % (
            nm, ra, rb, rb - ra, 'C' if covered[i] else '-', mk, nt_, ni, nf,
            ','.join('%d' % w for w, _ in fl[:6]), top))
    print('runtime %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
