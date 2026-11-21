#!/usr/bin/env python3
"""dqi_dq3_om36.py -- DQ3 (Nov 18): ROW 'om36': the extraction junction's in-place 24-bit inverse at space 36 in the
cheapest form of its family that fits 36 wires.

At space 36 (j_out with its pole flag live: entry 865, 7 calls a decode) inv24 gives the flag-free cached form
(outer omits 4 raw rows, narrow inner) = 183.  The family's minimum at 36 [M, scratchpad ghw12.py: every split of
PLAN24's and PLAN12's raw rows] is 181: outer omits nothing, the INNER cached 12-bit inverse omits its 4 quotient raw
rows during its GF(64) leaf and rebuilds them after it (9 terms, PLAN12's generalized Hamming weight d4 = 9), narrow
leaf: peak max(21 + 6 + 9, 36) = 36.  Values are the cached form's (only the order of the raw terms changes).
Install: token dqi_dq3_om36 AFTER dqi_dq3_inv24 in DQI_INT_EXT: it wraps the invert24_inplace that inv24 installed and
takes over only the calls with CEIL - live == 36 (keep=False).  Switches: DQI_DQ3_OM36_CEIL (901), DQI_DQ3_OM36_OFF=1,
DQI_DQ3_OM36_OUT=<json>.  Test: `python dqi_dq3_om36.py test [n]` (values against the full cached form on n inputs,
default all 2^24; Toffolis and peak by the light builder).
"""
import atexit
import importlib
import json
import os
import sys

CEIL = int(os.environ.get('DQI_DQ3_OM36_CEIL', '901'))
COUNT = 3            # the outer's omitted raw rows (rebuild 7)
OMIT = 1 << 9        # the inner's omitted raw wire (rebuild 3): peak 36 [M test], 182 on the line
STATS = dict(om36=0, passed=0, installed=0)
CALLS = []
sys.modules.setdefault('dqi_dq3_om36', sys.modules[__name__])


def make_inner(jet, ref, leaf, omit=0b1111000000):
    """cached_inverse(x, PLAN12, with_zero=False) with the raw coordinates in `omit` (a bitmask over PLAN12's 10 raw
    wires, never D's six) computed AFTER the GF(64) leaf: +1 CCX for every term that touches them.
    omit = 1 << 9: raw wire 9 is touched by 3 of the 15 terms (PLAN12's d1 = 3) [M ghw12]."""
    bits, rankone_xor = jet.bits, jet.rankone_xor
    assert not omit & 0b111111

    def inner(b, x, plan, with_zero):
        assert plan is jet.PLAN12 and not with_zero
        n, h = 12, 6
        size = plan['size']
        to, frm, K = jet.tower(n)
        jet.linear_inplace(b, x, to)
        terms = plan['terms']
        now = [j for j in range(size) if not omit >> j & 1]
        later = [j for j in range(size) if omit >> j & 1]
        wires = {}
        for j, w in zip(now, b.alloc_qubits(len(now))):
            wires[j] = w

        def part(vec, cols):
            return sum(1 << i for i, j in enumerate(cols) if vec >> j & 1)
        for left, right, vec in terms:                  # every term's part on the raw wires held now
            v = part(vec, now)
            if v:
                with ref.expose_pair(b, x, left << 1, right << 1) as (p, q):
                    rankone_xor(b, p, q, [wires[j] for j in now], v)
        D = [wires[j] for j in range(h)]
        jet.linear_xor(b, x, D, plan['affine'])
        k, F = leaf(b, D, with_zero)
        jet.linear_xor(b, x, D, plan['affine'])
        if later:
            for j, w in zip(later, b.alloc_qubits(len(later))):
                wires[j] = w
            for left, right, vec in terms:              # the omitted coordinates (the rebuilt part)
                v = part(vec, later)
                if v:
                    with ref.expose_pair(b, x, left << 1, right << 1) as (p, q):
                        rankone_xor(b, p, q, [wires[j] for j in later], v)
        raw = [wires[j] for j in range(size)]
        jet.linear_inplace(b, raw, plan['reduce'])
        kept = plan['size'] - len(plan['retire'])
        for q_, row in zip(raw[kept:], plan['retire']):
            mask = b.alloc_bits(1)[0]
            b.hmr(q_, mask)
            for j in bits(row):
                left, right, _ = terms[j]
                with ref.expose_pair(b, x, left << 1, right << 1) as (p, r):
                    b.ccz(p, r, mask)
            b.free(mask)
            b.free(q_)
        out = list(raw[:kept]) + list(b.alloc_qubits(n - kept))
        jet.linear_inplace(b, out, plan['out'])
        source = list(x) + list(k)
        for forms, vectors, lf in plan['leaves']:
            outputs, gates, low, offset = lf

            def remap(form, forms=forms):
                result = form & 1
                for j in bits(form >> 1):
                    result ^= forms[j] << 1
                return result
            for e, v in zip(outputs, vectors):
                if e & 1:
                    for j in bits(v):
                        b.x(out[j])
                for a_ in bits((e & (1 << low) - 1) >> 1):
                    for i in bits(forms[a_]):
                        for j in bits(v):
                            b.cx(source[i], out[j])
            for i, (left, right) in enumerate(gates):
                vec = 0
                for e, v in zip(outputs, vectors):
                    if e >> offset + i & 1:
                        vec ^= v
                with ref.expose_pair(b, source, remap(left), remap(right)) as (p, q):
                    rankone_xor(b, p, q, out, vec)
        jet.ar._erase_norm(b, out, k, K)
        b.free(k)
        jet.linear_inplace(b, x, frm)
        jet.linear_inplace(b, out, frm)
        return out, F
    return inner


def make_outer(sheared, xs, jet, inner):
    """sheared.cached_inverse_partial with the inner above (count 0 = no outer omission)."""
    return xs.clone(sheared.cached_inverse_partial, dict(sheared.cached_inverse_partial.__globals__,
                                                       cached_inverse=inner))


_K = {}


def _objects():
    if not _K:
        inv24 = sys.modules['dqi_dq3_inv24']
        ff, jet = inv24._kit()                       # the decode's own flag-free objects (narrow leaf = the line's)
        c24n = ff.CACHED24[1]
        inner_n = c24n.__globals__['cached_inverse']
        leaf = inner_n.__globals__['inverse6_fixed18']        # leaf_narrow as bound in the line's clone
        xs = sys.modules['decode_flagfree_cached'].xs
        inner = make_inner(jet, jet.ref, leaf, OMIT)
        _K['outer'] = make_outer(ff.sheared, xs, jet, inner)
        _K['jet'] = jet
        _K['inv24'] = inv24
    return _K


def invert24_om36(k, b, n):
    O = _objects()
    jet, inv24 = O['jet'], O['inv24']
    K = getattr(k, 'K', None)
    columns, back, _ = jet.extension_columns(K)
    jet.linear_inplace(b, n, columns)
    out, flag = O['outer'](b, n, jet.PLAN24, COUNT, False)
    assert flag is None
    t = b.alloc_bits(24)
    b.hmr(n, t)
    b.swap(n, out)
    b.free(out)
    jet.linear_inplace(b, n, back)
    key = ('M', tuple(columns))
    if key not in inv24._K:
        inv24._K[key] = inv24._columns_matrix(jet, columns)
    M = inv24._K[key]
    m0 = b.alloc_bits(12)
    m1 = b.alloc_bits(12)
    b.bit_store(m0, False)
    b.bit_store(m1, False)
    for j in range(24):
        col = M[j]
        for i in range(24):
            if col >> i & 1:
                if j < 12:
                    b.cx(t[i], m0[j])
                    b.cx(t[i], m1[j])
                else:
                    b.cx(t[i], m1[j - 12])
    b.free(t)
    return [m0, m1]


def _dump():
    path = os.environ.get('DQI_DQ3_OM36_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_om36: %s\n' % json.dumps(dict(STATS, ceil=CEIL)))


def install(km, ns):
    if os.environ.get('DQI_DQ3_OM36_OFF', '0') == '1':
        return dict(off=True)
    ip = sys.modules.get('dqi_junc_r8p_ipjunc')
    if ip is None or 'dqi_dq3_inv24' not in sys.modules:
        return dict(installed=0, why='needs ipjunc and dqi_dq3_inv24 (install om36 after inv24)')
    prev = ip.invert24_inplace

    def invert24_inplace(k, b, n, keep=False):
        live = b.num_allocated_qubits
        if not keep and CEIL - live == 36:
            STATS['om36'] += 1
            CALLS.append(live)
            return invert24_om36(k, b, n)
        STATS['passed'] += 1
        return prev(k, b, n, keep) if keep else prev(k, b, n)
    ip.invert24_inplace = invert24_inplace
    STATS['installed'] = 1
    atexit.register(_dump)
    return dict(STATS, ceil=CEIL)


def test(nin):
    import numpy as np
    root = os.path.expanduser('~/qcc/ws/dqi_dq3/gen/zc2p/home/qcc/ws/dqi_zn/src/kv_unfused5')
    tools = os.path.dirname(os.path.abspath(__file__))
    sys.argv = [sys.argv[0], root, str(nin)]
    sys.path.insert(0, tools)
    import dqi_dq3_ff24test as T                    # the kit objects as in the ff24 test (kit's narrow leaf)
    import dqi_dq3_lite as L
    inner = make_inner(T.jet, T.jet.ref, T.leaf_narrow, OMIT)
    outer = make_outer(T.sheared, T.xs, T.jet, inner)

    def form(b, x):
        columns, back, _ = T.jet.extension_columns(None)
        T.jet.linear_inplace(b, x, columns)
        out, flag = outer(b, x, T.jet.PLAN24, COUNT, False)
        t = b.alloc_bits(24)
        b.hmr(x, t)
        b.swap(x, out)
        b.free(out)
        T.jet.linear_inplace(b, x, back)

    def run(fn):
        b = L.LB(nin)
        x = b.alloc_qubits(24)
        idx = np.arange(nin, dtype=np.int64) if nin == 1 << 24 else np.random.default_rng(7).integers(0, 1 << 24, nin)
        for i in range(24):
            x[i].v[:] = (idx >> i) & 1
        fn(b, x)
        val = np.zeros(nin, dtype=np.int64)
        for i in range(24):
            val |= x[i].v.astype(np.int64) << i
        return val, b.toff, b.peak - 24
    a = run(lambda b, x: T.ff24(b, x, 0, 0))
    c = run(form)
    print('full cached (0, fast): Toffolis %d, peak extra %d' % (a[1], a[2]))
    print('om36 (outer omits %d, inner omit mask %s, kit narrow leaf): Toffolis %d, peak extra %d' % (COUNT, bin(OMIT), c[1], c[2]))
    print('values equal on %d inputs: %s (mismatches %d)' % (nin, bool((a[0] == c[0]).all()), int((a[0] != c[0]).sum())))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    test(int(sys.argv[2]) if len(sys.argv) > 2 else 1 << 24)
