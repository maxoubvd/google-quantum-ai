#!/usr/bin/env python3
"""invsqrt_pp_sim.py -- workstream invsqrt_pp (P1): exact all-input bit-sliced simulator for an invsqrt .kmb.

Every wire is a vector of 3,133,452 bits (input x = 1..3,133,452, one lane per input, padded with x = 1).
Z-type phases (Z, CZ, CCZ, NEG, HMR kick) are tracked as one XOR bit-vector P; rotations (Z_POW) are counted,
not tracked. Classical bits are scalars (an HMR outcome is the same on every lane: it does not depend on data).
Outcome policy: a function op index -> 0/1 (default: the zero transcript).

API
  C = load(path)                 parsed circuit (lists; C.base, C.qs, C.bit, C.cnd, C.n)
  x_lanes()                      (22, W) uint64: bits of x per lane
  run(C, hook=None, outcome=None, start_state=None, lo=0, hi=None)
       hook(i, S, P, bits, active) is called BEFORE op i (S: (nq, W) state, P phase, bits list, active: bool)
  block_effect(C, S, lo, hi, outcome)  run ops lo..hi (a block body incl. its PUSH/POP) as EXECUTED on a copy of
       S; returns (restored: bool, g phase vector, touched wires)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402

NX = 3133452
W = (NX + 63) // 64
ONES = np.uint64(0xFFFFFFFFFFFFFFFF)


class Circ:
    pass


def load(path):
    c = kmb.read(path)
    C = Circ()
    C.k = c
    C.nq = c.num_qubits
    C.nbits = c.num_bits
    ops = np.asarray(c.ops)
    q3 = [np.asarray(x).tolist() for x in c.q3]
    q2 = [np.asarray(x).tolist() for x in c.q2]
    q1 = np.asarray(c.q1).tolist()
    bt = np.asarray(c.bit_t).tolist()
    cd = np.asarray(c.cond).tolist()
    i3 = i2 = i1 = ib = ic = 0
    base, qs, bit, cnd = [], [], [], []
    for oi in range(len(ops)):
        op = int(ops[oi])
        if op == kmb.PUSH:
            base.append(kmb.PUSH)
            qs.append(())
            bit.append(-1)
            cnd.append(cd[ic])
            ic += 1
            continue
        b = op & 0x7F
        n = (op & 0x30) >> 4
        if n == 3:
            q = (q3[0][i3], q3[1][i3], q3[2][i3])
            i3 += 1
        elif n == 2:
            q = (q2[0][i2], q2[1][i2])
            i2 += 1
        elif n == 1:
            q = (q1[i1],)
            i1 += 1
        else:
            q = ()
        bb = -1
        if (op & 0x40) != 0 and b not in (0x06, 0x07):
            bb = bt[ib]
            ib += 1
        cc = -1
        if op & 0x80:
            cc = cd[ic]
            ic += 1
        base.append(b)
        qs.append(q)
        bit.append(bb)
        cnd.append(cc)
    C.base, C.qs, C.bit, C.cnd, C.n = base, qs, bit, cnd, len(base)
    xr = [r for r in c.registers if r[0] == 'x']
    C.xw = np.asarray(xr[0][1]).tolist() if xr else list(range(22))
    # block structure: matching PUSH/POP, level of every op (number of open blocks)
    C.match = {}
    C.level = [0] * C.n
    st = []
    for i in range(C.n):
        if base[i] == kmb.PUSH:
            C.level[i] = len(st)
            st.append(i)
        elif base[i] == kmb.POP:
            j = st.pop()
            C.match[j] = i
            C.level[i] = len(st)
        else:
            C.level[i] = len(st)
    return C


def x_lanes():
    x = np.arange(1, NX + 1, dtype=np.int64)
    x = np.concatenate([x, np.ones(W * 64 - NX, dtype=np.int64)])
    out = np.zeros((22, W), dtype=np.uint64)
    for j in range(22):
        bits = ((x >> j) & 1).astype(np.uint8)
        out[j] = np.packbits(bits, bitorder='little').view(np.uint64)
    return out


def init_state(C):
    S = np.zeros((C.nq, W), dtype=np.uint64)
    xl = x_lanes()
    for j, w in enumerate(C.xw):
        S[w] = xl[j]
    return S


def apply(C, i, S, P, bits, outcome, nrot):
    """Apply op i (unconditionally: the caller decided it is active). Returns updated P (array, in place)."""
    b = C.base[i]
    q = C.qs[i]
    if b == 0x31:
        S[q[2]] ^= S[q[0]] & S[q[1]]
    elif b == 0x21:
        S[q[1]] ^= S[q[0]]
    elif b == 0x11:
        S[q[0]] ^= ONES
    elif b == 0x20:
        P ^= S[q[0]] & S[q[1]]
    elif b == 0x10:
        P ^= S[q[0]]
    elif b == 0x50:
        m = outcome(i) if outcome else 0
        if m:
            P ^= S[q[0]]
        S[q[0]][:] = 0
        bits[C.bit[i]] = m
    elif b == 0x30:
        P ^= S[q[0]] & S[q[1]] & S[q[2]]
    elif b == 0x32:
        d = S[q[0]] & (S[q[1]] ^ S[q[2]])
        S[q[1]] ^= d
        S[q[2]] ^= d
    elif b == 0x22:
        t = S[q[0]].copy()
        S[q[0]][:] = S[q[1]]
        S[q[1]][:] = t
    elif b == 0x12:
        S[q[0]][:] = 0
    elif b == 0x00:
        P ^= ONES
    elif b == 0x40:
        bits[C.bit[i]] = 0
    elif b == 0x41:
        bits[C.bit[i]] = 1
    elif b == 0x42:
        bits[C.bit[i]] ^= 1
    elif b == 0x13:
        nrot[0] += 1
    return P


def run(C, hook=None, outcome=None, S=None, P=None, bits=None, lo=0, hi=None, active_stack=None):
    if S is None:
        S = init_state(C)
    if P is None:
        P = np.zeros(W, dtype=np.uint64)
    if bits is None:
        bits = [0] * C.nbits
    hi = C.n if hi is None else hi
    stk = [] if active_stack is None else list(active_stack)
    nrot = [0]
    for i in range(lo, hi):
        b = C.base[i]
        act = all(stk)
        if hook is not None:
            hook(i, S, P, bits, act)
        if b == kmb.PUSH:
            stk.append(bool(bits[C.cnd[i]]) and act)
            continue
        if b == kmb.POP:
            stk.pop()
            continue
        if b in (0x06, 0x07):
            continue
        if not act:
            continue
        if C.cnd[i] >= 0 and not bits[C.cnd[i]]:
            continue
        apply(C, i, S, P, bits, outcome, nrot)
    return S, P, bits, nrot[0]


def block_effect(C, S, bits, lo, outcome=None):
    """Block opened by the PUSH at lo, executed (its own condition forced true). Returns (restored, g, nrot)."""
    hi = C.match[lo]
    S2 = S.copy()
    P2 = np.zeros(W, dtype=np.uint64)
    b2 = list(bits)
    _, P2, _, nrot = run(C, outcome=outcome, S=S2, P=P2, bits=b2, lo=lo + 1, hi=hi, active_stack=[True])
    restored = bool(np.array_equal(S2, S))
    return restored, P2, nrot


if __name__ == '__main__':
    import time
    t0 = time.time()
    C = load(sys.argv[1])
    print('loaded', C.n, 'ops', C.nq, 'qubits', 'x wires', C.xw, round(time.time() - t0, 1), 's')
    S, P, bits, nrot = run(C)
    xl = x_lanes()
    ok = all(np.array_equal(S[w], xl[j]) for j, w in enumerate(C.xw))
    other = [w for w in range(C.nq) if w not in C.xw and S[w].any()]
    print('zero transcript: x restored', ok, 'dirty wires', other, 'rotations', nrot,
          'Z-phase nonzero lanes', int(np.unpackbits(P.view(np.uint8)).sum()), round(time.time() - t0, 1), 's')
