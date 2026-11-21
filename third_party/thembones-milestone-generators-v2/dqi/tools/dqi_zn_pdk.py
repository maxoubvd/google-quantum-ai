#!/usr/bin/env python3
"""dqi_zn_pdk.py -- builder Z1 (Nov 14): the paired payload digit kernels in kickmix (R33-DQI-CONV, step 2).

Lower level L of a paired block (U, L), inside L's pivot inverse context (a = g[-1] clipped, ai = 1/a), one payload:
  words: t = X (the upper level's digit, already t / a_U), s = the payload's word at L, low = n = L - 1 words,
  g = the divisor (L words, rotated), c = clip(q1') (q1' + [q1' == 0]), Z0 = [not B] q0'.
pd_down (descent):  D1 = X c (temp);  s ^= D1 g[-2] (= P);  s *= ai (= D0);  low ^= (D0 + D1 w) g[0..L-2] mod w^n
                    (n point products: the kit's evaluation product with N = n points aliases P r(w) into low,
                    r = w^n mod prod(w - p); that is cancelled for free by fanning s out with coefficients r
                    before and after P is added);  D1 erased by measurement;  s ^= Z0 X (= Y).
pd_up (return, the exact transpose): t ^= Z0 s;  [s, D1'] ^= transposed point products (S0 into s, S1 into a temp);
                    s *= ai;  D1' ^= (s + <r, low>) g[-2];  t ^= D1' c;  D1' erased by measurement (free phases).
Prices: pd_down = pd_up = (n + 4) x 42 executed Toffolis per payload (n + 1 of them are the kernel proper:
42 (delta + 1) with delta = n), two temporary words at most (D1 and the kit's scaling temp).
Usage: python dqi_zn_pdk.py test [shots=256] [seed=1]
"""
import os
import sys

import numpy as np

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)
SRC = os.environ.get('DQI_SW_SOURCE') or os.path.expanduser('~/qcc/ws/dqi_sw/inputs/source_cap902')  # line's kit
KIT = SRC + '/playbook/kit/decoder'
sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
import kickmix as km  # noqa: E402
sys.modules['circ'] = km
for p in (KIT, KIT + '/exp/convert', SRC + '/exp/r-main/narrow240'):
    if p not in sys.path:
        sys.path.insert(0, p)
import kernels as K  # noqa: E402  (the kit's exp/convert/kernels.py: cx, expose, rank_update, vandermonde)
import direct_kernels as DK  # noqa: E402  (conv with N points, zero-workspace rank updates)
from hankel_check import mul as fmul, inv as finv  # noqa: E402
from exp.safe_scalar.converter import scale as kit_scale, ar  # noqa: E402

_R = {}


def r_poly(n):
    """r(w) = w^n mod prod_{t<n} (w + t) over the kit's field (points t = 0..n-1 of K.vandermonde): the low
    coefficients of prod (w + t)."""
    if n not in _R:
        p = [1]
        for t in range(n):
            q = [0] * (len(p) + 1)
            for i, c in enumerate(p):
                q[i + 1] ^= c
                q[i] ^= fmul(c, t)
            p = q
        _R[n] = p[:n]
    return _R[n]


def fan_out(b, s, low, r):
    for j, rj in enumerate(r):
        K.cx(b, s, low[j], rj)


def fan_in(b, low, s, r):
    for j, rj in enumerate(r):
        K.cx(b, low[j], s, rj)


def pd_down(b, t, s, low, g, a, ai, c, Z0):
    n = len(low)
    assert len(g) == n + 1
    r = r_poly(n)
    fan_out(b, s, low, r)
    D1 = b.alloc_qubits(12)
    ar.mul_xor0(b, t, c, D1)
    ar.mul_xor0(b, D1, g[-2], s)
    fan_out(b, s, low, r)
    kit_scale(b, s, a, ai, True)
    DK.conv(b, g[:-1], [s, D1], low, n)
    ar.erase_product(b, t, c, D1)
    b.free(D1)
    ar.mul_xor0(b, Z0, t, s)


def conv_transpose(b, g, low, outs, N):
    """Transpose of DK.conv(b, g, [q0, q1], low, N): outs[0] ^= sum_t g(p_t) y_t, outs[1] ^= sum_t p_t g(p_t) y_t,
    y_t = <W[:, t], low> (exposed in place), one product per point with two targets."""
    V, W = K.vandermonde(N)
    for t in range(N):
        yc = [W[j][t] for j in range(len(low))]
        if any(yc):
            DK.rank_update(b, g, V[t][:len(g)], low, yc, outs, V[t][:len(outs)])


def erase_conv_transpose(b, g, low, acc, N, extra=None):
    """Erase the accumulator acc = sum_t p_t g(p_t) y_t (+ extra x * y) by measurement: free phase products,
    the classical mask carried through the public scaling by p_t with classical XORs."""
    V, W = K.vandermonde(N)
    f = ar.field(12)
    m = b.alloc_bits(12)
    b.hmr(acc, m)
    b.free(acc)
    for t in range(N):
        yc = [W[j][t] for j in range(len(low))]
        if not any(yc):
            continue
        pt = V[t][1] if len(V[t]) > 1 else 0
        if pt == 0:
            continue
        mk = b.alloc_bits(12)
        b.bit_store(mk, False)
        for j in range(12):                       # mk = transpose(multiply by p_t) applied to m
            v = f.mul(pt, 1 << j)
            for i in ar.bits(v):
                b.cx(m[i], mk[j])
        x, px = K.expose(b, g, V[t][:len(g)])
        y, py = K.expose(b, low, yc)
        ar.phase_product(b, x, y, mk)
        K.unexpose(b, low, yc, py)
        K.unexpose(b, g, V[t][:len(g)], px)
        b.free(mk)
    if extra is not None:
        ar.phase_product(b, extra[0], extra[1], m)
    b.free(m)


def pd_up(b, t, s, low, g, a, ai, c, Z0):
    n = len(low)
    assert len(g) == n + 1
    r = r_poly(n)
    ar.mul_xor0(b, Z0, s, t)
    D1 = b.alloc_qubits(12)
    conv_transpose(b, g[:-1], low, [s, D1], n)
    kit_scale(b, s, a, ai, True)
    fan_in(b, low, s, r)
    ar.mul_xor0(b, s, g[-2], D1)
    fan_in(b, low, s, r)
    ar.mul_xor0(b, D1, c, t)
    # D1 = S1 + (s + <r, low>) g[-2]: erase by measurement; the second term's phase with s shifted by <r, low>
    fan_in(b, low, s, r)
    erase_conv_transpose(b, g[:-1], low, D1, n, extra=(s, g[-2]))
    fan_in(b, low, s, r)


# ------------------------------------------------------------------------------------------------ test
def model_down(t, s, low, g, a, c, Z0):
    from dqi_x12c_gf import mul, inv0
    n = low.shape[1]
    ai = inv0(np.where(a == 0, 1, a))
    D1 = mul(t, c)
    s1 = s ^ mul(D1, g[:, -2])
    D0 = mul(s1, ai)
    low2 = low.copy()
    for i in range(n):
        acc = mul(D0, g[:, i])
        if i >= 1:
            acc ^= mul(D1, g[:, i - 1])
        low2[:, i] ^= acc
    return t, D0 ^ mul(Z0, t), low2


def model_up(t, s, low, g, a, c, Z0):
    from dqi_x12c_gf import mul, inv0
    n = low.shape[1]
    ai = inv0(np.where(a == 0, 1, a))
    t1 = t ^ mul(Z0, s)
    S0 = np.zeros_like(t)
    S1 = np.zeros_like(t)
    for i in range(n):
        S0 ^= mul(g[:, i], low[:, i])
        if i >= 1:
            S1 ^= mul(g[:, i - 1], low[:, i])
    s2 = mul(s ^ S0, ai)
    acc = S1 ^ mul(s2, g[:, -2])
    return t1 ^ mul(acc, c), s2, low


def build(kind, n):
    b = km.CircuitBuilder()
    W = lambda: b.alloc_qubits(12)  # noqa: E731
    regs = dict(t=W(), s=W(), c=W(), Z0=W(), a=W(), ai=W(), low=[W() for _ in range(n)], g=[W() for _ in range(n)])
    g = regs['g'] + [regs['a']]
    (pd_down if kind == 'down' else pd_up)(b, regs['t'], regs['s'], regs['low'], g, regs['a'], regs['ai'],
                                            regs['c'], regs['Z0'])
    return b.finish_circuit(), regs


def run(c, regs, vals, seed):
    import dqi_int_convtest as CT
    nshot = vals['t'].shape[0]
    sim = km.Simulator(batch_size=nshot, debug_print=False, counting_operations=True, seed=seed)
    sim.use_same_registers_as(c)
    flat = {}
    for k in ('t', 's', 'c', 'Z0', 'a', 'ai'):
        flat[k] = km.array(list(regs[k]))
        sim.write_across_shots(flat[k], CT.to_bits(vals[k][:, None]))
    for k in ('low', 'g'):
        flat[k] = km.array([q for w in regs[k] for q in w])
        sim.write_across_shots(flat[k], CT.to_bits(vals[k]))
    sim.do(c)
    out = {k: CT.from_bits(sim.read_across_shots(flat[k]), 1)[:, 0] for k in ('t', 's')}
    out['low'] = CT.from_bits(sim.read_across_shots(flat['low']), vals['low'].shape[1])
    nreg = 12 * (6 + 2 * vals['low'].shape[1])
    anc = km.array([km.q(i) for i in range(nreg, c.num_qubits)]) if c.num_qubits > nreg else None
    dirty = sim.read_across_shots(anc).any(axis=0) if anc is not None else np.zeros(nshot, dtype=bool)
    phase = np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool)
    return out, dirty, phase, sim.operation_statistics().num_toffoli_like / nshot


def test(shots, seed):
    from dqi_x12c_gf import inv0
    rng = np.random.default_rng(seed)
    for n in (15, 13, 11, 9, 3):
        for kind in ('down', 'up'):
            c, regs = build(kind, n)
            vals = {k: rng.integers(0, 4096, size=shots, dtype=np.int64) for k in ('t', 's', 'c', 'Z0', 'a')}
            vals['a'][: shots // 8] = 0                       # clipped pivots: the register holds 1 in context
            vals['a'] = np.where(vals['a'] == 0, 1, vals['a'])
            vals['ai'] = inv0(vals['a'])
            vals['c'] = np.where(vals['c'] == 0, 1, vals['c'])
            vals['low'] = rng.integers(0, 4096, size=(shots, n), dtype=np.int64)
            vals['g'] = rng.integers(0, 4096, size=(shots, n), dtype=np.int64)
            gfull = np.concatenate([vals['g'], vals['a'][:, None]], axis=1)
            out, dirty, phase, tof = run(c, regs, vals, seed)
            ref = (model_down if kind == 'down' else model_up)(vals['t'], vals['s'], vals['low'], gfull, vals['a'],
                                                                vals['c'], vals['Z0'])
            ok = (out['t'] == ref[0]) & (out['s'] == ref[1]) & (out['low'] == ref[2]).all(axis=1)
            print('n %2d %-4s exact %d/%d, dirty %d, phase %d, qubits %d (registers %d), executed %.1f (%.1f per 42: '
                  'asked n + 4 = %d)' % (n, kind, ok.sum(), shots, dirty.sum(), phase.sum(), c.num_qubits,
                                         12 * (6 + 2 * n), tof, tof / 42, n + 4), flush=True)


if __name__ == '__main__':
    if sys.argv[1] == 'test':
        test(int(sys.argv[2]) if len(sys.argv) > 2 else 256, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
