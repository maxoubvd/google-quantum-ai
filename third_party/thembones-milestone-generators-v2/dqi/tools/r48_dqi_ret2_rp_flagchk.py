#!/usr/bin/env python3
"""r48_dqi_ret2_rp_flagchk.py -- R48-DQI-RET2 re-pricer: can each candidate flag park of tools/r48_dqi_ret2_rp_idle.py
be closed by ONE zero test of the word the flag was made from? Read only.

Per park (flag wire f, gap lt..nt, birth op, from the census's "park wire" lines): the tested word W = the controls of
the CCX chain in the birth's marker segment (up to the birth op) minus the chain's targets and f (the kit's zero_toggle:
11 CCX over the 12 negated bits). STATIC: the ops that WRITE a wire of W between the birth and the close (CCX target,
CSWAP / SWAP pair, CX target, X, R, HMR). SAMPLES (SCAN's simulate_snapop, 64 lanes per run; the verifier's law and
uniform weights alternately, as tools/r46_dqi_scan_heldchk.py): f against [W == 0] in both polarities just after the
open (lt) and just before the close (nt), and the number of lanes with f = 1 at the close.
ROOM (with a zero.npz): the re-test's chain needs 11 ancillas + the flag at once, so it must sit at an op t after the
last run at live >= LEVEL inside the gap and before the close, with live(t) <= LEVEL - 12 (the flag counted in live);
t* = the latest such op; f against [W == 0] is also sampled at t*.
Usage: python r48_dqi_ret2_rp_flagchk.py <dqi.kmb> <segments_real.json> <idle census txt> [R=16] [seed=7]
       [zero.npz] [LEVEL=902]
"""
import json
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
from dqi_pp_act import pack, v, km  # noqa: E402
from r46_dqi_scan_rank import simulate_snapop  # noqa: E402

ALL = np.uint64(0xFFFFFFFFFFFFFFFF)


def pc(x):
    return bin(int(x)).count('1')


def main():
    t0 = time.time()
    path, segp, cens = sys.argv[1:4]
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 16
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 7
    parks = []
    for line in open(cens):
        m = re.match(r'\s+park wire\s+(\d+)\s+gap (\d+)\.\.(\d+) .* born (\d+) in \[(.*?)\]', line)
        if m:
            parks.append(tuple(int(m.group(k)) for k in range(1, 5)) + (m.group(5),))
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
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
    S = json.load(open(segp))
    segs = S['segments']
    mark = np.flatnonzero(base == 0x07)
    print('%s: %d parks from %s (%.0fs)' % (path, len(parks), cens, time.time() - t0), flush=True)

    def word_of(f, birth):
        j = int(np.searchsorted(mark, birth, 'right')) - 1
        for jj in (j, j - 1):
            lo = int(mark[jj]) if jj >= 0 else 0
            sl = np.arange(lo, birth + 1)
            cc = sl[(base[sl] == 0x31) & (k[sl] == 3)]
            ctl = set(q3[0][i3[cc]].tolist()) | set(q3[1][i3[cc]].tolist())
            tg = set(q3[2][i3[cc]].tolist())
            W = sorted(ctl - tg - {f})
            if len(W) == 12:
                return W, jj
        return W, j

    def writes(W, a, b):
        sl = np.arange(a + 1, b)
        Ws = np.zeros(nq, bool)
        Ws[W] = True
        cnt = 0
        s3 = sl[k[sl] == 3]
        b3 = base[s3]
        x = s3[b3 == 0x31]
        cnt += int(Ws[q3[2][i3[x]]].sum())
        x = s3[b3 == 0x32]
        cnt += int(Ws[q3[1][i3[x]]].sum() + Ws[q3[2][i3[x]]].sum())
        s2 = sl[k[sl] == 2]
        b2 = base[s2]
        x = s2[b2 == 0x21]
        cnt += int(Ws[q2[1][i2[x]]].sum())
        x = s2[b2 == 0x22]
        cnt += int(Ws[q2[0][i2[x]]].sum() + Ws[q2[1][i2[x]]].sum())
        s1 = sl[k[sl] == 1]
        b1 = base[s1]
        x = s1[(b1 == 0x11) | (b1 == 0x12) | (b1 == 0x50)]
        cnt += int(Ws[q1[i1[x]]].sum())
        return cnt

    live = None
    LEVEL = int(sys.argv[7]) if len(sys.argv) > 7 else 902
    if len(sys.argv) > 6:
        Z = np.load(sys.argv[6])
        zs = np.clip(Z['start'].astype(np.int64), 0, n)
        ze = np.clip(Z['end'].astype(np.int64), 0, n)
        okz = ze > zs
        d = np.zeros(n + 2, np.int32)
        np.add.at(d, zs[okz], 1)
        np.add.at(d, ze[okz], -1)
        live = nq - np.cumsum(d[:n])
    info = []
    for (f, lt, nt, birth, bk) in parks:
        W, jj = word_of(f, birth)
        wb = writes(W, birth, nt) if len(W) == 12 else -1
        ts, room = -1, ''
        if live is not None:
            seg = live[lt + 1:nt]
            hi = np.flatnonzero(seg >= LEVEL)
            last_end = lt + 1 + (int(hi[-1]) + 1 if len(hi) else 0)
            cand = np.flatnonzero(live[last_end:nt] <= LEVEL - 12)
            if len(cand):
                ts = last_end + int(cand[-1])
            room = 'last run at >= %d ends %d; latest op with live <= %d before the close: %s (%s ops before it)' % (
                LEVEL, last_end, LEVEL - 12, ts if ts >= 0 else 'NONE', nt - ts if ts >= 0 else '-')
        info.append((f, lt, nt, birth, bk, W, wb, ts))
        print('park f %4d gap %d..%d born %d [%s]: word %s (%d wires); writes to the word between birth and close: %d; %s'
              % (f, lt, nt, birth, bk, ','.join(map(str, W[:12])), len(W), wb, room), flush=True)
    sites = np.unique(np.array([x[1] + 1 for x in info] + [x[2] for x in info] + [x[7] for x in info if x[7] >= 0],
                               np.int64))
    ns = len(sites)
    snap = np.zeros((ns, nq), np.uint64)
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    law = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    uni = np.full(spec.max_errors + 1, 1.0 / (spec.max_errors + 1))
    rng = np.random.default_rng(seed)
    P = len(info)
    stat = np.zeros((P, 8), np.int64)    # open eq / neg mismatch, close eq / neg, ones at close, ones at open, t* eq / neg
    for r in range(R):
        vinit = np.zeros(nq, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, law if r % 2 == 0 else uni)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        js = simulate_snapop(ops, *arr, vinit, binit, seed * 1000 + r, sites, snap)
        assert js == ns
        for p, (f, lt, nt, birth, bk, W, wb, ts) in enumerate(info):
            if len(W) != 12:
                continue
            for col, site in ((0, lt + 1), (2, nt), (6, ts)):
                if site < 0:
                    continue
                row = snap[int(np.searchsorted(sites, site))]
                z = (~np.bitwise_or.reduce(row[W])) & ALL
                fv = row[f]
                stat[p, col] += pc(fv ^ z)
                stat[p, col + 1] += pc(fv ^ (~z & ALL))
                if col != 6:
                    stat[p, 4 if col == 2 else 5] += pc(fv)
        print('run %d done (%.0fs)' % (r, time.time() - t0), flush=True)
    print('\nRESULT on %d samples (f = flag at the open / close; z = [W == 0])' % (64 * R))
    for p, (f, lt, nt, birth, bk, W, wb, ts) in enumerate(info):
        oe, on, ce, cn, ones_c, ones_o, te, tn = stat[p].tolist()
        pol = 'f = z' if oe == 0 else ('f = not z' if on == 0 else 'neither')
        ok = (oe == 0 and ce == 0) or (on == 0 and cn == 0)
        okt = ts >= 0 and ((oe == 0 and te == 0) or (on == 0 and tn == 0))
        print('f %4d gap %9d..%-9d writes %3d | open: %s (mismatch %d / %d), close mismatch %d / %d, f=1 at close on %d '
              'lanes | RE-TEST AT THE CLOSE %s | re-test op t* %s: %s' % (
                  f, lt, nt, wb, pol, oe, on, ce, cn, ones_c,
                  'EXACT ON SAMPLES' if ok and len(W) == 12 else 'NOT VALID', ts,
                  'EXACT ON SAMPLES (live %d)' % int(live[ts]) if okt else ('no room' if ts < 0 else 'NOT VALID')))
    print('runtime %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
