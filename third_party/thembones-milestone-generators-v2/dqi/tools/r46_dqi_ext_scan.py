#!/usr/bin/env python3
"""r46_dqi_ext_scan.py -- R46-DQI-EXT: gate-level scans of a DQI decoder built with the r7prof markers.

A marker file (DEBUG_PRINT ops, one per profile segment; segments_real.json beside it) lets every Toffoli-class op
be charged to the call stack that emitted it. Three modes:

  stat <dqi.kmb> <out.npz> [runs=64] [seed=5]
      P1's bit-parallel value model (tools/dqi_pp_act.py: valid verifier samples, random HMR outcomes), with, for
      every static Toffoli-class op: act (samples on which it executes), ones (its product is 1), pa / pb (first /
      second control is 1; CSWAP: control, x + y), pt (target wire is 1 before the op), and seg (index of the
      marker that opened its segment).
  table <stat.npz> <segments_real.json> <stage> [depth=3] [top=60]
      executed Toffolis per decode of one stage by the innermost <depth> frames of the stack, with the operand
      statistics of each row (min / max of P(a), P(b), P(product) over its ops; ops with a control never 1).
  open <dqi.kmb> <segments_real.json> <segment index> [count=1] [detail=0] [skip Toffolis=0] [take Toffolis=0]
      symbolic scan of the ops of <count> consecutive segments: every wire is a GF(2) linear form of the wires at
      entry plus one fresh variable per AND / measurement; for every CCX the two control forms and the set of wires
      that receive the product at the end; prints the counts, the rank of the control forms, the supports by
      12-wire word, and the rank of the products as quadratic polynomials modulo affine functions (a deficit = a
      Toffoli that is a sum of others and of live wires).
  rankall <dqi.kmb> <segments_real.json> <stage> [frames=3] [top=40]
      the same rank test on every segment of a stage, chunk by chunk (see the function's note).
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402


def stat(path, out, R, seed):
    from numba import njit
    from dqi_pp_act import pack, v, km, _rng, popc64
    ALL = np.uint64(0xFFFFFFFFFFFFFFFF)

    @njit
    def sim(ops, q3a, q3b, q3c, q2a, q2b, q1, bt, cd, vinit, binit, seed, act, ones, pa, pb, pt):
        vv = vinit.copy()
        bits = binit.copy()
        s = np.uint64(seed) * np.uint64(0x9E3779B97F4A7C15) + np.uint64(1)
        stack = np.zeros(256, np.uint64)
        sp = 0
        active = ALL
        i3 = 0
        i2 = 0
        i1 = 0
        ib = 0
        ic = 0
        k = 0
        for oi in range(ops.shape[0]):
            op = ops[oi]
            if op == 0x81:
                c = cd[ic]
                ic += 1
                stack[sp] = active
                sp += 1
                active = active & bits[c]
                continue
            base = op & 0x7F
            if base == 0x01:
                if sp > 0:
                    sp -= 1
                    active = stack[sp]
                continue
            nqk = (op & 0x30) >> 4
            a = 0
            b = 0
            t = 0
            if nqk == 3:
                a = q3a[i3]
                b = q3b[i3]
                t = q3c[i3]
                i3 += 1
            elif nqk == 2:
                a = q2a[i2]
                b = q2b[i2]
                i2 += 1
            elif nqk == 1:
                a = q1[i1]
                i1 += 1
            bit = 0
            if (op & 0x40) != 0 and base != 0x06 and base != 0x07:
                bit = bt[ib]
                ib += 1
            m = active
            if (op & 0x80) != 0:
                m = m & bits[cd[ic]]
                ic += 1
            if base == 0x11:
                vv[a] ^= m
            elif base == 0x12:
                vv[a] &= ~m
            elif base == 0x21:
                vv[b] ^= vv[a] & m
            elif base == 0x22:
                d = (vv[a] ^ vv[b]) & m
                vv[a] ^= d
                vv[b] ^= d
            elif base == 0x31:
                d = vv[a] & vv[b] & m
                act[k] += popc64(m)
                ones[k] += popc64(d)
                pa[k] += popc64(vv[a] & m)
                pb[k] += popc64(vv[b] & m)
                pt[k] += popc64(vv[t] & m)
                vv[t] ^= d
                k += 1
            elif base == 0x32:
                x = vv[b] ^ vv[t]
                d = x & vv[a] & m
                act[k] += popc64(m)
                ones[k] += popc64(d)
                pa[k] += popc64(vv[a] & m)
                pb[k] += popc64(x & m)
                pt[k] += popc64(vv[t] & m)
                vv[b] ^= d
                vv[t] ^= d
                k += 1
            elif base == 0x30:
                d = vv[a] & vv[b] & vv[t] & m
                act[k] += popc64(m)
                ones[k] += popc64(d)
                pa[k] += popc64(vv[a] & m)
                pb[k] += popc64(vv[b] & m)
                pt[k] += popc64(vv[t] & m)
                k += 1
            elif base == 0x50:
                s = _rng(s)
                bits[bit] = (bits[bit] & ~m) | (s & m)
                vv[a] &= ~m
            elif base == 0x40:
                bits[bit] &= ~m
            elif base == 0x41:
                bits[bit] |= m
            elif base == 0x42:
                bits[bit] ^= m
        return k

    t0 = time.time()
    c = kmb.read(path)
    ops = np.ascontiguousarray(np.asarray(c.ops))
    toff = ((ops & 0x7C) == 0x30) & (ops != 0x81)
    t_op = np.flatnonzero(toff)
    nt = len(t_op)
    mark = np.flatnonzero(ops == 0x07)
    seg = np.searchsorted(mark, t_op, side='right') - 1
    arr = [np.ascontiguousarray(np.asarray(x, dtype=np.int64)) for x in (c.q3[0], c.q3[1], c.q3[2], c.q2[0], c.q2[1],
                                                                        c.q1, c.bit_t, c.cond)]
    spec = v.ProblemSpecDQI.benchmark_instance(num_shots=v.DEFAULT_NUM_SHOTS)
    tables = v.FieldTables(km.GF2Field(spec.degree))
    wprob = v.optimal_dqi_weights(spec.m, spec.max_errors, spec.degree, spec.r) ** 2
    if os.environ.get('R46_WLAW'):
        # second session: a UNIFORM law on the weights lo..hi instead of the verifier's (R46_WLAW=18-35)
        lo, hi = (int(x) for x in os.environ['R46_WLAW'].split('-'))
        wprob = np.zeros_like(wprob)
        wprob[lo:hi + 1] = 1.0 / (hi + 1 - lo)
        print('weight law: uniform on %d..%d (%d entries)' % (lo, hi, len(wprob)), flush=True)
    rng = np.random.default_rng(seed)
    cnt = [np.zeros(nt, np.int64) for _ in range(5)]
    print('read %.0fs: %d ops, %d Toffoli-class static, %d markers' % (time.time() - t0, len(ops), nt, len(mark)),
          flush=True)
    for r in range(R):
        vinit = np.zeros(c.num_qubits, np.uint64)
        binit = np.zeros(c.num_bits, np.uint64)
        syn, mbuc, _ = v.sample_shots(spec, tables, rng, 64, wprob)
        vinit[:840] = pack(v.to_bits(syn, spec.degree))
        mb = pack(v.to_bits(mbuc, spec.degree))
        binit[:len(mb)] = mb
        k = sim(ops, *arr, vinit, binit, seed * 1000 + r, *cnt)
        assert k == nt
        if (r + 1) % 16 == 0 or r == R - 1:
            n = 64 * (r + 1)
            print('run %d: %d samples, mean executed %.1f (%.0fs)' % (r + 1, n, cnt[0].sum() / n, time.time() - t0),
                  flush=True)
            np.savez(out + '.tmp.npz', t_op=t_op, base=(ops[t_op] & 0x7F), act=cnt[0], ones=cnt[1], pa=cnt[2],
                     pb=cnt[3], pt=cnt[4], seg=seg, samples=n)
            os.replace(out + '.tmp.npz', out)
    print('done %.0fs' % (time.time() - t0))


def table(npz, segf, stage, depth, top):
    d = np.load(npz)
    S = json.load(open(segf))
    marks, st, segs = S['marks'], S['stacks'], S['segments']
    n = int(d['samples'])
    act, ones, pa, pb, seg = d['act'], d['ones'], d['pa'], d['pb'], d['seg']
    assert seg.max() < len(segs), (int(seg.max()), len(segs))
    rows = {}
    tot = 0.0
    for k in range(len(act)):
        s, mk, _c = segs[seg[k]]
        if not marks[mk] or marks[mk][0] != stage:
            continue
        names = [x.split('@')[0] for x in st[s]]
        key = ' > '.join(names[-depth:])
        e = act[k] / n
        tot += e
        r = rows.setdefault(key, [0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0, 0])
        r[0] += 1
        r[1] += e
        if act[k]:
            fa, fb, fo = pa[k] / act[k], pb[k] / act[k], ones[k] / act[k]
            r[2] = min(r[2], fa)
            r[3] = max(r[3], fa)
            r[4] = min(r[4], fb)
            r[5] = max(r[5], fb)
            r[6] = min(r[6], fo)
            r[7] = max(r[7], fo)
            if pa[k] == 0 or pb[k] == 0:
                r[8] += 1
            if ones[k] == 0:
                r[9] += 1
    print('stage %s: executed per decode %.1f (%d samples)' % (stage, tot, n))
    print('%-78s %6s %9s %6s  %-11s %-11s %-11s %5s %5s' % ('innermost %d frames' % depth, 'static', 'executed',
                                                              'share', 'P(a) range', 'P(b) range', 'P(ab) range',
                                                              'ctl0', 'prod0'))
    for key in sorted(rows, key=lambda q: -rows[q][1])[:top]:
        r = rows[key]
        print('%-78s %6d %9.1f %5.1f%%  %.3f-%.3f %.3f-%.3f %.3f-%.3f %5d %5d'
              % (key[-78:], r[0], r[1], 100 * r[1] / tot, r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]))


def _rank(vecs):
    piv = {}
    rk = 0
    for x in vecs:
        while x:
            h = x.bit_length() - 1
            if h in piv:
                x ^= piv[h]
            else:
                piv[h] = x
                rk += 1
                break
    return rk


def _words(mask, nq):
    out = {}
    i = 0
    while mask:
        if mask & 1:
            w = 'w%d' % (i // 12) if i < 840 else 'anc'
            out[w] = out.get(w, 0) + 1
        mask >>= 1
        i += 1
    return out


def open_(path, segf, s0, count, detail, skip_t=0, take_t=0):
    c = kmb.read(path)
    S = json.load(open(segf))
    ops = np.asarray(c.ops)
    mark = np.flatnonzero(ops == 0x07)
    p0 = int(mark[s0])
    p1 = int(mark[s0 + count]) if s0 + count < len(mark) else len(ops)
    if skip_t or take_t:                        # a slice of the range by Toffoli-class count (one call of a kernel)
        sl = np.asarray(ops[p0:p1])
        tp = np.flatnonzero(((sl & 0x7C) == 0x30) & (sl != 0x81))
        q0 = p0 + (int(tp[skip_t - 1]) + 1 if skip_t else 0)
        q1 = p0 + int(tp[skip_t + take_t]) if take_t and skip_t + take_t < len(tp) else p1
        print('slice: Toffolis %d..%d of %d in the range' % (skip_t, skip_t + take_t, len(tp)))
        p0, p1 = q0, q1
    nq = c.num_qubits
    i3 = int(np.count_nonzero(((ops[:p0] & 0x30) >> 4) == 3))
    i2 = int(np.count_nonzero(((ops[:p0] & 0x30) >> 4) == 2))
    i1 = int(np.count_nonzero(((ops[:p0] & 0x30) >> 4) == 1))
    q3 = [np.asarray(x[i3:i3 + (p1 - p0)]) for x in c.q3]
    q2 = [np.asarray(x[i2:i2 + (p1 - p0)]) for x in c.q2]
    q1 = np.asarray(c.q1[i1:i1 + (p1 - p0)])
    seg = np.asarray(ops[p0:p1])
    for j in range(count):
        s, mk, stc = S['segments'][s0 + j]
        print('segment %d: static %d, marks %s, stack ... %s' % (s0 + j, stc, '/'.join(S['marks'][mk]),
                                                                 ' > '.join(x.split('@')[0] for x in S['stacks'][s][-5:])))
    ONE = 1 << nq
    form = [1 << i for i in range(nq)]
    nvar = nq + 1
    j3 = j2 = j1 = 0
    terms = []          # (kind, form a, form b, target wire, variable)
    kinds = {}
    cond_ops = 0
    hmr = 0
    depth = 0
    for op in seg.tolist():
        if op == 0x81:
            depth += 1
            continue
        base = op & 0x7F
        if base == 0x01:
            depth = max(0, depth - 1)
            continue
        nqk = (op & 0x30) >> 4
        a = b = t = 0
        if nqk == 3:
            a, b, t = int(q3[0][j3]), int(q3[1][j3]), int(q3[2][j3])
            j3 += 1
        elif nqk == 2:
            a, b = int(q2[0][j2]), int(q2[1][j2])
            j2 += 1
        elif nqk == 1:
            a = int(q1[j1])
            j1 += 1
        cond = depth > 0 or (op & 0x80) != 0
        kinds[base] = kinds.get(base, 0) + 1
        if base in (0x11, 0x21, 0x22) and cond:
            cond_ops += 1                       # classically conditioned Clifford: an unknown bit times a known form
            v_ = 1 << nvar
            nvar += 1
            if base == 0x11:
                form[a] ^= v_
            elif base == 0x21:
                form[b] ^= v_
            else:
                form[a] ^= v_
                form[b] ^= v_
            continue
        if base == 0x11:
            form[a] ^= ONE
        elif base == 0x12:
            form[a] = 0
        elif base == 0x21:
            form[b] ^= form[a]
        elif base == 0x22:
            form[a], form[b] = form[b], form[a]
        elif base == 0x31:
            v_ = 1 << nvar
            nvar += 1
            terms.append(('ccx', form[a], form[b], t, v_))
            form[t] ^= v_
        elif base == 0x32:
            v_ = 1 << nvar
            nvar += 1
            terms.append(('cswap', form[a], form[b] ^ form[t], t, v_))
            form[b] ^= v_
            form[t] ^= v_
        elif base == 0x30:
            terms.append(('ccz', form[a], form[b], t, 0))
        elif base == 0x50:
            hmr += 1
            form[a] = 0
    names = {0x11: 'X', 0x12: 'R', 0x21: 'CX', 0x22: 'SWAP', 0x31: 'CCX', 0x32: 'CSWAP', 0x30: 'CCZ', 0x50: 'HMR',
             0x20: 'CZ', 0x10: 'Z', 0x07: 'MARK', 0x40: 'BIT0', 0x41: 'BIT1', 0x42: 'BITINV', 0x06: 'TAG'}
    print('ops %d..%d (%d): %s; conditioned Cliffords %d' % (p0, p1, p1 - p0, ', '.join(
        '%s %d' % (names.get(k, hex(k)), n) for k, n in sorted(kinds.items(), key=lambda q: -q[1])), cond_ops))
    lin = (1 << (nq + 1)) - 1
    tt = [x for x in terms if x[0] != 'ccz']
    if not tt:
        return
    A = [x[1] for x in tt]
    B = [x[2] for x in tt]
    nl_a = sum(1 for x in A if x >> (nq + 1))
    nl_b = sum(1 for x in B if x >> (nq + 1))
    ua = ub = 0
    for x in A:
        ua |= x & (ONE - 1)
    for x in B:
        ub |= x & (ONE - 1)
    print('Toffolis with a product: %d (CCX %d, CSWAP %d); CCZ %d; HMR %d' % (len(tt), sum(1 for x in tt if x[0] == 'ccx'),
                                                                          sum(1 for x in tt if x[0] == 'cswap'),
                                                                          len(terms) - len(tt), hmr))
    print('control a: rank %d, entry-wire support %d %s, forms with an earlier product in them %d'
          % (_rank(A), bin(ua).count('1'), _words(ua, nq), nl_a))
    print('control b: rank %d, entry-wire support %d %s, forms with an earlier product in them %d'
          % (_rank(B), bin(ub).count('1'), _words(ub, nq), nl_b))
    # the products as quadratic polynomials modulo affine functions
    mono = {}
    vecs = []
    for kind, fa, fb, t, v_ in tt:
        ia = [i for i in range(fa.bit_length()) if fa >> i & 1]
        ib_ = [i for i in range(fb.bit_length()) if fb >> i & 1]
        acc = {}
        for i in ia:
            for j in ib_:
                if i == j or i == nq or j == nq:
                    continue                      # x x = x, x 1 = x: affine
                key = (i, j) if i < j else (j, i)
                acc[key] = acc.get(key, 0) ^ 1
        x = 0
        for key, on in acc.items():
            if on:
                if key not in mono:
                    mono[key] = len(mono)
                x |= 1 << mono[key]
        vecs.append(x)
    rk = _rank(vecs)
    print('products as quadratic polynomials modulo affine: %d products, rank %d (deficit %d), %d monomials, zero %d'
          % (len(vecs), rk, len(vecs) - rk, len(mono), sum(1 for x in vecs if x == 0)))
    # where each product ends
    outs = {}
    for w in range(nq):
        f = form[w] >> (nq + 1)
        i = 0
        while f:
            if f & 1:
                outs.setdefault(nq + 1 + i, []).append(w)
            f >>= 1
            i += 1
    dest = {}
    lost = 0
    for kind, fa, fb, t, v_ in tt:
        ws = outs.get(v_.bit_length() - 1, [])
        if not ws:
            lost += 1
        key = tuple(sorted({('w%d' % (w // 12) if w < 840 else 'anc') for w in ws}))
        dest[key] = dest.get(key, 0) + 1
    print('products by the words that hold them at the end of the range (none = measured out or cancelled: %d):' % lost)
    for key, n in sorted(dest.items(), key=lambda q: -q[1])[:24]:
        print('   %4d  %s' % (n, ' '.join(key) if key else '(none)'))
    chg = [w for w in range(nq) if form[w] != (1 << w)]
    print('wires whose value differs from entry at the end: %d %s' % (len(chg), _words(sum(1 << w for w in chg), nq)))
    if detail:
        for n_, (kind, fa, fb, t, v_) in enumerate(tt[:detail]):
            print('  %3d %s |a| %2d |b| %2d target %d outs %s' % (n_, kind, bin(fa).count('1'), bin(fb).count('1'), t,
                                                              outs.get(v_.bit_length() - 1, [])))


def rankall(path, segf, stage, frames, top):
    """Every segment of one stage, cut into chunks of products whose controls are linear in the wires at the chunk's
    start (a product that reads an AND output of its own chunk opens a new chunk; the wire values at a chunk's start
    are taken as independent variables: every dependency reported holds as an identity, some across chunks are not
    seen). Per kernel (innermost frames): products, chunks, and the deficit = products minus the rank of the
    products as quadratic polynomials modulo affine functions (a product that is a sum of others and of wires)."""
    c = kmb.read(path)
    S = json.load(open(segf))
    marks, st, segs = S['marks'], S['stacks'], S['segments']
    ops = np.asarray(c.ops)
    mark = np.flatnonzero(ops == 0x07)
    idx = [i for i, (s_, mk, n) in enumerate(segs) if marks[mk] and marks[mk][0] == stage]
    s0, s1 = idx[0], idx[-1] + 1
    assert idx == list(range(s0, s1)), 'stage segments not contiguous'
    p0 = int(mark[s0])
    p1 = int(mark[s1]) if s1 < len(mark) else len(ops)
    nq = c.num_qubits
    nk = (ops[:p0] & 0x30) >> 4
    i3, i2, i1 = (int(np.count_nonzero(nk == k)) for k in (3, 2, 1))
    del nk
    n = p1 - p0
    q3 = [np.asarray(x[i3:i3 + n]).tolist() for x in c.q3]
    q2 = [np.asarray(x[i2:i2 + n]).tolist() for x in c.q2]
    q1 = np.asarray(c.q1[i1:i1 + n]).tolist()
    seg_ops = np.asarray(ops[p0:p1]).tolist()
    ident = [1 << i for i in range(nq)]
    ONE = 1 << nq
    rows = {}
    state = dict(form=list(ident), nvar=nq + 1, prods=[], key=None)

    def flush():
        prods = state['prods']
        if prods:
            mono = {}
            vecs = []
            zero = 0
            for fa, fb in prods:
                ia = [i for i in range(fa.bit_length()) if fa >> i & 1]
                ib_ = [i for i in range(fb.bit_length()) if fb >> i & 1]
                acc = {}
                for i in ia:
                    for j in ib_:
                        if i == j or i == nq or j == nq:
                            continue
                        key = (i, j) if i < j else (j, i)
                        acc[key] = acc.get(key, 0) ^ 1
                x = 0
                for key, on in acc.items():
                    if on:
                        if key not in mono:
                            mono[key] = len(mono)
                        x |= 1 << mono[key]
                zero += x == 0
                vecs.append(x)
            rk = _rank(vecs)
            r = rows.setdefault(state['key'], [0, 0, 0, 0])
            r[0] += len(prods)
            r[1] += 1
            r[2] += len(prods) - rk
            r[3] += zero
        state['prods'] = []
        state['form'] = list(ident)
        state['nvar'] = nq + 1

    j3 = j2 = j1 = 0
    depth = 0
    si = s0 - 1
    for op in seg_ops:
        if op == 0x07:
            flush()
            si += 1
            s_, mk, _n = segs[si]
            state['key'] = ' > '.join(x.split('@')[0] for x in st[s_][-frames:])
            continue
        if op == 0x81:
            depth += 1
            continue
        base = op & 0x7F
        if base == 0x01:
            depth = max(0, depth - 1)
            continue
        nqk = (op & 0x30) >> 4
        a = b = t = 0
        if nqk == 3:
            a, b, t = q3[0][j3], q3[1][j3], q3[2][j3]
            j3 += 1
        elif nqk == 2:
            a, b = q2[0][j2], q2[1][j2]
            j2 += 1
        elif nqk == 1:
            a = q1[j1]
            j1 += 1
        form = state['form']
        cond = depth > 0 or (op & 0x80) != 0
        if base in (0x11, 0x21, 0x22) and cond:
            v_ = 1 << state['nvar']
            state['nvar'] += 1
            if base == 0x11:
                form[a] ^= v_
            elif base == 0x21:
                form[b] ^= v_
            else:
                form[a] ^= v_
                form[b] ^= v_
            continue
        if base == 0x11:
            form[a] ^= ONE
        elif base == 0x12:
            form[a] = 0
        elif base == 0x21:
            form[b] ^= form[a]
        elif base == 0x22:
            form[a], form[b] = form[b], form[a]
        elif base == 0x31 or base == 0x32:
            fa = form[a]
            fb = form[b] if base == 0x31 else form[b] ^ form[t]
            if (fa | fb) >> (nq + 1):               # reads a product (or an unknown) of this chunk: new chunk
                flush()
                form = state['form']
                fa = form[a]
                fb = form[b] if base == 0x31 else form[b] ^ form[t]
            state['prods'].append((fa, fb))
            v_ = 1 << state['nvar']
            state['nvar'] += 1
            form[t] ^= v_
            if base == 0x32:
                form[b] ^= v_
        elif base == 0x50:
            v_ = 1 << state['nvar']                  # measured and reset: 0 afterwards (kept exact)
            state['nvar'] += 1
            form[a] = 0
    flush()
    tp = sum(r[0] for r in rows.values())
    td = sum(r[2] for r in rows.values())
    print('stage %s: %d products (CCX + CSWAP, static), %d chunks, deficit %d, products that are affine %d'
          % (stage, tp, sum(r[1] for r in rows.values()), td, sum(r[3] for r in rows.values())))
    print('%-70s %8s %7s %7s %7s' % ('kernel (innermost %d frames)' % frames, 'products', 'chunks', 'deficit', 'affine'))
    for key in sorted(rows, key=lambda k: (-rows[k][2], -rows[k][0]))[:top]:
        r = rows[key]
        print('%-70s %8d %7d %7d %7d' % (key[-70:], r[0], r[1], r[2], r[3]))



def main():
    mode = sys.argv[1]
    if mode == 'stat':
        stat(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 64,
             int(sys.argv[5]) if len(sys.argv) > 5 else 5)
    elif mode == 'table':
        table(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 3,
              int(sys.argv[6]) if len(sys.argv) > 6 else 60)
    elif mode == 'rankall':
        rankall(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]) if len(sys.argv) > 5 else 3,
                int(sys.argv[6]) if len(sys.argv) > 6 else 40)
    elif mode == 'open':
        open_(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]) if len(sys.argv) > 5 else 1,
              int(sys.argv[6]) if len(sys.argv) > 6 else 0, int(sys.argv[7]) if len(sys.argv) > 7 else 0,
              int(sys.argv[8]) if len(sys.argv) > 8 else 0)


if __name__ == '__main__':
    main()
