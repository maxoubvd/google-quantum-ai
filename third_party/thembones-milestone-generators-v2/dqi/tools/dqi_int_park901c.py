#!/usr/bin/env python3
"""dqi_int_park901c.py -- B5 (17th holder, Nov 16): = dqi_int_park901b.py (unchanged there) + B5_PARK_MOVEZ=1: a 'trans'
line whose ops are all DIAGONAL (CZ / Z / CCZ, with or without a classical condition) on a flag that is remade later
gets no transient: its ops are moved to just after the flag's remake (exact: a diagonal gate commutes with every gate
that leaves the computational values of its wires alone; checked: no op between writes its other wires, the parked
flag is constant; each condition bit is copied into a fresh bit at the original place (BIT_STORE0, then BIT_STORE1
conditioned on it) and the moved op is conditioned on the copy). Saves the transient's 11 Toffolis.
Original docstring of park901b: wire parks for the 901 file as a .kmb post-pass, from a SPEC file.
Spec lines (op indices of the INPUT file; '#' comments):
  park  <f> <lt> <nt> <W12>      HMR f right after op lt (bit b); before op nt: the kit's 11-CCX zero test of W into f
                                  (f = [W == 0]), its ancillas erased by HMR + conditioned CZ, then Z_IF(b) on f.
  open  <op> <f> <name>          HMR f right after op (bit saved under <name>).
  trans <op1> <op2> <f> <W12>    the ops op1..op2 (adjacent) that read f get a TRANSIENT z = [W == 0] instead of f:
                                  zero test into a fresh z, the ops with f replaced by z, z and the chain erased by
                                  measurement (11 Toffolis). Exact where z equals f on the inputs those ops act on.
  remake <op> <f> <W12> <name>   before op: zero test of W into f (f must be |0>), then Z_IF(bit <name>) on f.
  hmr   <op> <w>                 HMR w right after op (a wire that is |0> on every input; bit unused).
Fresh ancillas: wires num_qubits .. +10 (a re-assignment pass places them). Every inserted HMR takes a new bit.
Checks: no op in (open, remake) of a parked wire reads it except the trans ops; inserts outside PUSH_CONDITION.
Usage: python dqi_int_park901c.py <in.kmb> <spec> <out.kmb>   [B5_PARK_MOVEZ=1]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

HMR, X, CCX, CZ_IF, Z_IF = 0x50, 0x11, 0x31, 0xA0, 0x90


def main():
    t0 = time.time()
    src, spec, dst = sys.argv[1:4]
    k = kmb.read(src)
    ops = np.asarray(k.ops)
    n = len(ops)
    dep = np.cumsum((ops == 0x81).astype(np.int64) - (ops == 0x01).astype(np.int64))
    kk = ((ops & 0x30) >> 4).astype(np.int8)
    i1 = np.cumsum(kk == 1) - 1
    i2 = np.cumsum(kk == 2) - 1
    i3 = np.cumsum(kk == 3) - 1
    hc = (ops & 0x80) != 0
    ci = np.cumsum(hc) - 1
    hb = ((ops & 0x40) != 0) & ((ops & 0x7F) != 0x06) & ((ops & 0x7F) != 0x07)
    bi = np.cumsum(hb) - 1
    M = kmb.ID_MASK
    q1 = np.asarray(k.q1)
    q2 = [np.asarray(x) for x in k.q2]
    q3 = [np.asarray(x) for x in k.q3]
    cdv = np.asarray(k.cond)
    btv = np.asarray(k.bit_t)

    def operands(t):
        c = int(kk[t])
        if c == 3:
            return [int(q3[j][i3[t]]) & M for j in range(3)]
        if c == 2:
            return [int(q2[j][i2[t]]) & M for j in range(2)]
        if c == 1:
            return [int(q1[i1[t]]) & M]
        return []

    nb = [k.num_bits]
    anc = list(range(k.num_qubits, k.num_qubits + 11))     # 10 chain ancillas + z
    inserts, delete, bits = {}, [], {}
    toff = [0]

    def newbit():
        b = nb[0]
        nb[0] += 1
        return b

    def chain(W, target):
        """X on W, 11 CCX into target (= [W == 0]); returns (compute ops, erase ops for the 10 ancillas + X)."""
        a = anc[:10]
        comp = [(X, (w,), 0, 0) for w in W]
        comp.append((CCX, (W[0], W[1], a[0]), 0, 0))
        for i in range(1, 10):
            comp.append((CCX, (a[i - 1], W[i + 1], a[i]), 0, 0))
        comp.append((CCX, (a[9], W[11], target), 0, 0))
        toff[0] += 11
        return comp, a

    def erase_chain(W, a):
        out = []
        for i in range(9, -1, -1):
            cb = newbit()
            out.append((HMR, (a[i],), cb, 0))
            out.append((CZ_IF, (W[0], W[1]) if i == 0 else (a[i - 1], W[i + 1]), 0, cb))
        out += [(X, (w,), 0, 0) for w in W]
        return out

    def after(op):                       # an insertion point at or after op with block depth 0
        o = op
        while dep[o] != 0:
            o += 1
        return o

    def before(op):                      # insertion point (insert after this index) right before op, depth 0
        c = op - 1
        while dep[c] != 0:
            c -= 1
        return c

    MOVEZ = os.environ.get('B5_PARK_MOVEZ', '0') == '1'
    moved = []
    DIAG = (0x10, 0x20, 0x30)
    WRITES = {0x11: [0], 0x12: [0], 0x50: [0], 0x21: [1], 0x31: [2], 0x22: [0, 1], 0x32: [1, 2]}

    def movable(o1, o2, f, lines):
        rs = [int(x[1]) for x in lines if x and x[0] == 'remake' and int(x[2]) == f and int(x[1]) > o2]
        if not rs:
            print('movez %d..%d: no later remake' % (o1, o2))
            return False
        ro = rs[0]
        others, bitsr = set(), set()
        for t in range(o1, o2 + 1):
            if int(ops[t] & 0x7F) not in DIAG or dep[t] != 0:
                print('movez %d..%d: op %d not diagonal or in a block (%s, dep %d)' % (o1, o2, t, hex(int(ops[t])), int(dep[t])))
                return False
            qs = operands(t)
            if f not in qs:
                return False
            others.update(q for q in qs if q != f)
            if hc[t]:
                bitsr.add(int(cdv[ci[t]]))
        for t in range(o2 + 1, ro):
            code = int(ops[t] & 0x7F)
            if code in WRITES:
                qs = operands(t)
                if any(qs[i] in others for i in WRITES[code] if i < len(qs)):
                    print('movez %d..%d: op %d writes %s' % (o1, o2, t, qs))
                    return False
        return True

    parked = []                          # (wire, from op, to op, allowed reader ops)
    trans_readers = {}
    lines = [ln.split('#')[0].split() for ln in open(spec)]
    for p in [x for x in lines if x]:
        kind = p[0]
        if kind == 'park':
            f, lt, nt = int(p[1]), int(p[2]), int(p[3])
            W = [int(x) for x in p[4].split(',')]
            b = newbit()
            inserts.setdefault(after(lt), []).append((HMR, (f,), b, 0))
            comp, a = chain(W, f)
            inserts.setdefault(before(nt), []).extend(comp + erase_chain(W, a) + [(Z_IF, (f,), 0, b)])
            parked.append((f, lt, nt, set()))
        elif kind == 'open':
            op, f, name = int(p[1]), int(p[2]), p[3]
            b = newbit()
            bits[name] = (b, f, op)
            inserts.setdefault(after(op), []).append((HMR, (f,), b, 0))
        elif kind == 'trans' and MOVEZ and movable(int(p[1]), int(p[2]), int(p[3]), lines):
            o1, o2, f = int(p[1]), int(p[2]), int(p[3])
            ro = [int(x[1]) for x in lines if x and x[0] == 'remake' and int(x[2]) == f and int(x[1]) > o2][0]
            mv, cp = [], []
            for t in range(o1, o2 + 1):
                cd = int(cdv[ci[t]]) if hc[t] else 0
                bt = int(btv[bi[t]]) if hb[t] else 0
                if hc[t]:                    # the condition bit may be reused later: copy it now (classical, free)
                    nbit = newbit()
                    cp += [(0x40, (), nbit, 0), (0xC1, (), nbit, cd)]
                    cd = nbit
                mv.append((int(ops[t]), tuple(operands(t)), bt, cd))
                delete.append(t)
            inserts.setdefault(o2, []).extend(cp)
            moved.append((ro, mv))
            trans_readers.setdefault(f, set()).update(range(o1, o2 + 1))
            print('movez: ops %d..%d (%d diagonal) moved after the remake at %d' % (o1, o2, len(mv), ro), flush=True)
        elif kind == 'trans':
            o1, o2, f = int(p[1]), int(p[2]), int(p[3])
            W = [int(x) for x in p[4].split(',')]
            z = anc[10]
            comp, a = chain(W, z)
            mods = []
            for t in range(o1, o2 + 1):
                qs = operands(t)
                assert f in qs, ('trans op does not read the wire', t, qs)
                assert not (int(kk[t]) == 3 and int(ops[t] & 0x7F) == CCX and qs[2] == f), 'trans op writes the wire'
                qs2 = tuple(z if q == f else q for q in qs)
                cd = int(cdv[ci[t]]) if hc[t] else 0
                bt = int(btv[bi[t]]) if hb[t] else 0
                mods.append((int(ops[t]), qs2, bt, cd))
                delete.append(t)
            m = newbit()
            ez = [(HMR, (z,), m, 0), (CZ_IF, (a[9], W[11]), 0, m)]
            assert dep[o1 - 1] == 0 and dep[o2] == 0, 'trans inside a block'
            inserts.setdefault(o2, []).extend(comp + mods + ez + erase_chain(W, a))
            trans_readers.setdefault(f, set()).update(range(o1, o2 + 1))
        elif kind == 'remake':
            op, f, name = int(p[1]), int(p[2]), p[4]
            W = [int(x) for x in p[3].split(',')]
            b, f0, o0 = bits[name]
            assert f0 == f
            comp, a = chain(W, f)
            inserts.setdefault(before(op), []).extend(comp + erase_chain(W, a) + [(Z_IF, (f,), 0, b)])
            parked.append((f, o0, op, set()))
        elif kind == 'hmr':
            op, w = int(p[1]), int(p[2])
            inserts.setdefault(after(op), []).append((HMR, (w,), newbit(), 0))
        else:
            raise SystemExit('unknown spec line %s' % p)
    for ro, mv in moved:
        inserts.setdefault(before(ro), []).extend(mv)
    # readers of parked wires strictly inside their parked spans
    for (f, a0, b0, allowed) in parked:
        allowed = allowed | trans_readers.get(f, set())
        lo, hi = a0 + 1, b0
        bad = [t for t in range(lo, hi) if t not in allowed and f in operands(t)] if hi - lo < 3_000_000 else None
        if bad is None:                  # long span: vectorised check
            sl = np.arange(lo, hi)
            hit = np.zeros(len(sl), bool)
            for cls, arrs, ix in ((3, q3, i3), (2, q2, i2), (1, [q1], i1)):
                m_ = kk[sl] == cls
                for arr_ in arrs:
                    hit[m_] |= (np.asarray(arr_)[ix[sl[m_]]] & M) == f
            bad = [int(t) for t in sl[hit] if int(t) not in allowed]
        print('parked wire %d over %d..%d: other readers %s' % (f, a0, b0, bad[:5]), flush=True)
        assert not bad, 'a parked wire is read inside its span'
    out = kw.edit(k, delete, inserts)
    out.num_bits = nb[0]
    nq = kw.write(dst, out)
    o2 = np.asarray(out.ops)
    tof = int((((o2 & 0x7C) == 0x30) & (o2 != 0x81)).sum())
    t_in = int((((ops & 0x7C) == 0x30) & (ops != 0x81)).sum())
    print('written %s: Q %d, bits %d, ops %d, static Toffoli-class %d (%+d; chains %d) (%.0fs)' % (
        dst, nq, nb[0], len(o2), tof, tof - t_in, toff[0], time.time() - t0))


if __name__ == '__main__':
    main()
