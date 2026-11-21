#!/usr/bin/env python3
"""dqi_dq3_lite.py -- DQ3 (Nov 18): a light bit-parallel VALUE builder for exhaustive kernel tests (no op records,
no phases): every wire is a bool vector over all test inputs.  The subset of the kickmix builder API the kit's
inverse kernels call; counts Toffoli-class ops with quantum controls (broadcast lengths), tracks live / peak
qubits, and asserts that freed qubits are clean.  Used by tools/dqi_dq3_ff24test.py.
"""
from contextlib import contextmanager

import numpy as np


class W:
    __slots__ = ('v', 'kind')

    def __init__(self, n, kind):
        self.v = np.zeros(n, dtype=bool)
        self.kind = kind


class LB:
    def __init__(self, n, live0=0):
        self.n = n
        self.live_q = 0
        self.peak = 0
        self.toff = 0
        self.base = live0          # what num_allocated_qubits reports beyond the test's own wires
        self.cond = []

    @property
    def num_allocated_qubits(self):
        return self.base + self.live_q

    def alloc_qubits(self, k):
        self.live_q += k
        self.peak = max(self.peak, self.live_q)
        return [W(self.n, 'q') for _ in range(k)]

    def alloc_bits(self, k):
        return [W(self.n, 'c') for _ in range(k)]

    def free(self, ws):
        ws = ws if isinstance(ws, (list, tuple)) else [ws]
        for w in ws:
            if w.kind == 'q':
                assert not w.v.any(), 'freeing a dirty qubit'
                self.live_q -= 1

    def bit_store(self, ws, val, control=None):
        for w in (ws if isinstance(ws, (list, tuple)) else [ws]):
            if control is None:
                w.v[:] = bool(val)
            else:
                w.v[control.v] = bool(val)

    @contextmanager
    def mark(self, *a, **k):
        yield

    def _m(self):
        m = None
        for c in self.cond:
            m = c.v if m is None else (m & c.v)
        return m

    def push_condition(self, c):
        self.cond.append(c)

    def pop_condition(self):
        self.cond.pop()

    @staticmethod
    def _pairs(*args):
        lists = [a if isinstance(a, (list, tuple)) else None for a in args]
        n = max((len(l) for l in lists if l is not None), default=1)
        return [[(a[i] if isinstance(a, (list, tuple)) else a) for a in args] for i in range(n)]

    def _apply(self, t, val):
        m = self._m()
        if m is None:
            t.v ^= val
        else:
            t.v ^= val & m

    def x(self, t):
        for (q,) in self._pairs(t):
            self._apply(q, np.ones(self.n, dtype=bool))

    def cx(self, a, t):
        for p, q in self._pairs(a, t):
            self._apply(q, p.v if isinstance(p, W) else np.full(self.n, bool(p)))

    def swap(self, a, t):
        assert not self.cond
        for p, q in self._pairs(a, t):
            p.v, q.v = q.v.copy(), p.v.copy()

    def _isq(self, w):
        return isinstance(w, W) and w.kind == 'q'

    def ccx(self, a, b, t):
        for p, q, r in self._pairs(a, b, t):
            if self._isq(p) and self._isq(q):
                self.toff += 1
            pv = p.v if isinstance(p, W) else np.full(self.n, bool(p))
            qv = q.v if isinstance(q, W) else np.full(self.n, bool(q))
            self._apply(r, pv & qv)

    def init_and(self, a, b, t):
        for p, q, r in self._pairs(a, b, t):
            assert not r.v.any()
            if self._isq(p) and self._isq(q):
                self.toff += 1
            self._apply(r, p.v & q.v)

    def del_and(self, a, b, t):
        for p, q, r in self._pairs(a, b, t):
            self._apply(r, p.v & q.v)
            assert not r.v.any(), 'del_and on a wrong value'

    def cswap(self, c, a, t):
        for p, q, r in self._pairs(c, a, t):
            if self._isq(p):
                self.toff += 1
            m = p.v.copy()
            if self._m() is not None:
                m &= self._m()
            d = (q.v ^ r.v) & m
            q.v ^= d
            r.v ^= d

    def ccz(self, a, b, c):
        for p, q, r in self._pairs(a, b, c):
            if self._isq(p) and self._isq(q) and self._isq(r):
                self.toff += 1

    def z(self, *a):
        pass

    def s(self, *a):
        pass

    def cz(self, *a):
        pass

    def hmr(self, q, m):
        for p, r in self._pairs(q, m):
            p.v[:] = False
            r.v[:] = False          # values only: the outcome is irrelevant to the values that follow
