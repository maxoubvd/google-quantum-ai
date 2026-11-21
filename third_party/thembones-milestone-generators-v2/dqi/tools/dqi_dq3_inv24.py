#!/usr/bin/env python3
"""dqi_dq3_inv24.py -- DQ3 (Nov 18): ROW 'inv24' of the census: the extraction junction's in-place 24-bit inverse in
the kit's flag-free CACHED form where the width allows.

tools/dqi_junc_r8p_ipjunc.invert24_inplace (read only) inverts the carried scale n in place: iota = 1/Norm(n) (norm 42
+ cached 12-bit inverse 58, peak +34), conj, two K-products (84): 184 Toffolis.  The line runs it 14 times a decode in
j_out_nopole at entry 864 / 865 [M census cen3].  The kit's flag-free cached 24-bit inverse (exp/frontier/root/
decode_flagfree_cached: sheared.cached_inverse_partial, choose_outer) costs 170 + its rebuild at the space it gets:
179 at space 37 (entry 864 under 901), 183 at space 36 [M census + tools/dqi_dq3_ff24test.py].  This module gives
invert24_inplace(k, b, n, keep=False) that form whenever choose_outer(CEIL - live) finds an option cheaper than 184,
and keeps the original otherwise (and always for keep=True).
Contract kept: n -> 1/n in place in the same extension basis; returns [m0, m1], two 12-bit classical masks with the
owed phase (-1)^(m0 . (n0 + n1) + m1 . n1) on the OLD halves (what pay_conj / HalfPay pay).  The flag-free form leaves
a 24-bit transcript t on columns(x_old) (owed (-1)^(t . M x_old)); w = M^T t is computed by classical CX on bits, then
m0 = w[:12], m1 = w[12:] ^ w[:12] (w . x = w_lo . n0 + w_hi . n1 = m0 . (n0 + n1) + m1 . n1).  Free.
Install: token dqi_dq3_inv24 LAST in DQI_INT_EXT.  Switches: DQI_DQ3_INV24_CEIL (901), DQI_DQ3_INV24_OFF=1,
DQI_DQ3_INV24_OUT=<json> (per call: live at entry, (count, narrow) or 'orig').
"""
import atexit
import importlib
import json
import os
import sys

CEIL = int(os.environ.get('DQI_DQ3_INV24_CEIL', '901'))
STATS = dict(cached=0, orig=0, installed=0)
CALLS = []
_K = {}
sys.modules.setdefault('dqi_dq3_inv24', sys.modules[__name__])


class _NS:
    pass


def _leaf_name(c24):
    """name of the narrow leaf the cached 24-bit form finally calls (the line rebinds it in cloned scopes)."""
    try:
        inner = c24[1].__globals__['cached_inverse']
        leaf = inner.__globals__['inverse6_fixed18']
        return leaf.__globals__['inverse6_initialized_nine'].__name__
    except Exception:
        return None


def _kit():
    """the decode's own flag-free objects (choose_outer, CACHED24), preferring the clone whose narrow leaf is the
    line's rebound one (measured: leaf_narrow -> chosen -> inverse6_xor20, 20 Toffolis); module objects otherwise."""
    if not _K:
        import gc
        import types
        ff = sys.modules.get('decode_flagfree_cached') or importlib.import_module('decode_flagfree_cached')
        best = None
        for obj in gc.get_objects():
            if (isinstance(obj, types.FunctionType) and obj.__code__.co_name == 'inverse24_forward_flagfree'
                    and obj.__code__.co_filename.endswith('decode_flagfree_cached.py')):
                g = obj.__globals__
                if 'CACHED24' in g and _leaf_name(g['CACHED24']) not in (None, 'inverse6_initialized_nine'):
                    best = g
                    break
        ns = _NS()
        g = best if best is not None else ff.__dict__
        ns.CACHED24 = g.get('CACHED24', ff.CACHED24)
        ns.choose_outer = g.get('choose_outer', ff.choose_outer)       # pure functions of the plan: the module's
        ns.sheared = g.get('sheared', ff.sheared)
        _K['ff'] = ns
        _K['jet'] = ff.jet
        _K['leaf'] = _leaf_name(g['CACHED24'])
    return _K['ff'], _K['jet']


def _columns_matrix(jet, cols):
    """M e_j as 24-bit masks: apply jet.linear_inplace to unit vectors with a tiny value builder."""
    out = []
    for j in range(24):
        v = [1 if i == j else 0 for i in range(24)]
        out.append(sum(bit << i for i, bit in enumerate(_lin(jet, v, cols))))
    return out


def _lin(jet, v, cols):
    class Q:
        __slots__ = ('x',)

        def __init__(self, x):
            self.x = x

    class B:
        def cx(self, a, t):
            if isinstance(a, (list, tuple)):
                for p, q in zip(a, t):
                    q.x ^= p.x
            else:
                t.x ^= a.x

        def swap(self, a, t):
            if isinstance(a, (list, tuple)):
                for p, q in zip(a, t):
                    p.x, q.x = q.x, p.x
            else:
                a.x, t.x = t.x, a.x
    w = [Q(x) for x in v]
    jet.linear_inplace(B(), w, cols)
    return [q.x for q in w]


def invert24_cached(k, b, n, orig):
    ff, jet = _kit()
    live = b.num_allocated_qubits
    sel = ff.choose_outer(CEIL - live)
    if sel is None:
        STATS['orig'] += 1
        CALLS.append((live, 'orig'))
        return orig(k, b, n)
    count, narrow = sel
    # price of the chosen form against the original's 184: 170 + rebuild(count) + 2 a narrow inner
    columns_l, _ = ff.sheared.partial_layout(count)
    rebuild = sum(bool(jet.apply(columns_l, v) >> (21 - count)) for _, _, v in jet.PLAN24['terms'])
    if 170 + rebuild + 2 * narrow >= 184:
        STATS['orig'] += 1
        CALLS.append((live, 'orig'))
        return orig(k, b, n)
    STATS['cached'] += 1
    CALLS.append((live, [count, narrow]))
    K = getattr(k, 'K', None)
    columns, back, _ = jet.extension_columns(K)
    jet.linear_inplace(b, n, columns)
    out, flag = ff.CACHED24[narrow](b, n, jet.PLAN24, count, False)
    assert flag is None
    t = b.alloc_bits(24)
    b.hmr(n, t)
    b.swap(n, out)
    b.free(out)
    jet.linear_inplace(b, n, back)
    # owed phase (-1)^(t . M x_old), M = the 'columns' map: masks in invert24_inplace's convention
    key = ('M', tuple(columns))
    if key not in _K:
        _K[key] = _columns_matrix(jet, columns)
    M = _K[key]
    m0 = b.alloc_bits(12)
    m1 = b.alloc_bits(12)
    b.bit_store(m0, False)
    b.bit_store(m1, False)
    for j in range(24):
        col = M[j]                       # bit i of col = M_ij;  w_j = parity(t & col)
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
    path = os.environ.get('DQI_DQ3_INV24_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, leaf=_K.get('leaf'), calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_dq3_inv24: %s\n' % json.dumps(dict(STATS, ceil=CEIL)))


def install(km, ns):
    if os.environ.get('DQI_DQ3_INV24_OFF', '0') == '1':
        return dict(off=True)
    ip = sys.modules.get('dqi_junc_r8p_ipjunc')
    if ip is None:
        return dict(installed=0, why='dqi_junc_r8p_ipjunc not loaded')
    orig = ip.invert24_inplace
    if getattr(orig, '_dq3_inv24', False):
        return dict(again=True)

    def invert24_inplace(k, b, n, keep=False):
        if keep:
            STATS['orig'] += 1
            CALLS.append((b.num_allocated_qubits, 'keep'))
            return orig(k, b, n, keep=True)
        return invert24_cached(k, b, n, orig)
    invert24_inplace._dq3_inv24 = True
    ip.invert24_inplace = invert24_inplace       # j_out_nopole / ret_ab read the module global at call time
    STATS['installed'] = 1
    atexit.register(_dump)
    return dict(STATS, ceil=CEIL)
