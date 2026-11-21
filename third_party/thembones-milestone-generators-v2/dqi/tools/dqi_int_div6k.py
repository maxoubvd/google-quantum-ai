#!/usr/bin/env python3
"""dqi_int_div6k.py -- B5 (23rd holder, Nov 18): the un-synthesis division's 6-bit subnorm REBUILT AT 14 ANDs where one
wire is idle (the 17 unflagged late-M divisions sit at slack 1 [M DQ3 census cen3: divide_noflag_lb extra 24, slack 1]).

The line's div6 (tools/dqi_junc_r8p_div6.py, read only) inverts M = Norm6(N) in place (produce, 21 Toffolis) and later
rebuilds M = Norm6(N) into a fresh register (rebuild: the norm's linear part + the GF(64) product N_lo * N_hi, i.e. 15
symmetric rank-one terms p_i = (m_i . N_lo)(m_i . N_hi) XORed into the coordinates v_i), only to pay phases with it.
Here, where the width allows:
  produce: one CX first copies c = M[J] onto a fresh wire, J a coordinate where the output vector v_R of term R has a 1;
  rebuild: the norm WITHOUT term R (14 ANDs) gives M' = M + v_R p_R; c ^= M'[J] leaves c = p_R; M' ^= v_R c gives M
  exactly; c (= p_R, a product of two parities of N) is cleared by kickmix's del_and on the two exposed parities (HMR,
  CZ if the outcome is 1): 0 Toffolis.  -1 Toffoli a division; the same M and the same masks, so every later phase
  payment is unchanged.  Exactness test: `python dqi_int_div6k.py test [shots] [seed]` (K1's division bench, both forms).
Install: token dqi_int_div6k LAST in DQI_INT_EXT (after the line's div6).  Gate at produce entry: live + EXTRA <= CEIL
(DQI_INT_D6K_CEIL 901; DQI_INT_D6K_EXTRA 7 = the span's own extra 6 over produce's entry + the kept wire).  Switches:
DQI_INT_D6K_OFF=1, DQI_INT_D6K_R (term index, default 0), DQI_INT_D6K_OUT=<json> (per call: live at produce entry, form).
"""
import atexit
import json
import os
import sys
import types

CEIL = int(os.environ.get('DQI_INT_D6K_CEIL', '901'))
EXTRA = int(os.environ.get('DQI_INT_D6K_EXTRA', '7'))
RTERM = int(os.environ.get('DQI_INT_D6K_R', '0'))
STATS = dict(kept=0, orig=0, rebuild=0, rebuild_orig=0, scopes=0)
CALLS = []


class _Kept(list):
    """div6's two classical mask arrays, plus the kept wire (attribute b5c)."""
    b5c = None


def _make(G, produce0, rebuild0):
    norm = G['norm_xor_noscratch']
    ng = norm.__globals__
    linear_xor, square_columns, field = ng['linear_xor'], ng['square_columns'], ng['field']
    mul = ng['mul_xor_noscratch']
    mg = mul.__globals__
    tensor, bits, tensor_xor = mg['tensor'], mg['bits'], mg['tensor_xor_noscratch']
    ar = G['ar']
    TO6, FRM6, _K6 = ar.tower(6)
    terms = list(tensor(6))
    assert len(terms) == 15, len(terms)
    mask_r, vec_r = terms[RTERM]
    J = min(bits(vec_r))
    rest = terms[:RTERM] + terms[RTERM + 1:]
    ixr = list(bits(mask_r))

    def produce(b, M):
        live = b.num_allocated_qubits
        if live + EXTRA > CEIL:
            STATS['orig'] += 1
            CALLS.append((live, 0))
            return produce0(b, M)
        c = b.alloc_qubits(1)
        b.cx(M[J], c[0])                                   # c = M[J] (the input of the in-place inverse)
        zero, masks = produce0(b, M)
        kept = _Kept(masks)
        kept.b5c = c
        STATS['kept'] += 1
        CALLS.append((live, 1))
        return zero, kept

    def rebuild(b, N, kept, k6):
        c = getattr(kept, 'b5c', None)
        if c is None:
            STATS['rebuild_orig'] += 1
            return rebuild0(b, N, kept, k6)
        STATS['rebuild'] += 1
        M = b.alloc_qubits(6)
        f = field(6)                                       # norm_xor_noscratch(b, N, M, k6) without term R
        linear_xor(b, N[:6], M, square_columns(6))
        linear_xor(b, N[6:], M, [f.mul(k6, v) for v in square_columns(6)])
        tensor_xor(b, N[:6], N[6:], M, rest)               # M = Norm6(N) + v_R p_R
        b.cx(M[J], c[0])                                   # c = p_R
        for i in bits(vec_r):
            b.cx(c[0], M[i])                               # M = Norm6(N)
        p = ixr[0]                                         # clear c = (m_R . N_lo)(m_R . N_hi): 0 Toffolis
        for j in ixr[1:]:
            b.cx(N[j], N[p])
            b.cx(N[6 + j], N[6 + p])
        b.del_and(N[p], N[6 + p], c[0])
        for j in reversed(ixr[1:]):
            b.cx(N[6 + j], N[6 + p])
            b.cx(N[j], N[p])
        b.free(c)
        ar.linear_inplace(b, M, TO6)                       # div6's rebuild from here on, unchanged
        b.cx(M[3:], M[:3])
        b.cz(M[:3], kept[0])
        b.cz(M[3:], kept[1])
        b.cx(M[3:], M[:3])
        ar.linear_inplace(b, M, FRM6)
        for m in kept:
            b.free(m)
        return M
    produce._b5_d6k = True
    rebuild._b5_d6k = True
    return produce, rebuild


def _dump():
    path = os.environ.get('DQI_INT_D6K_OUT')
    if path:
        try:
            json.dump(dict(STATS, ceil=CEIL, extra=EXTRA, term=RTERM, calls=CALLS), open(path, 'w'))
        except Exception:
            pass
    sys.stderr.write('dqi_int_div6k: %s\n' % json.dumps(dict(STATS, ceil=CEIL, extra=EXTRA, term=RTERM)))


def install(km, ns):
    if os.environ.get('DQI_INT_D6K_OFF', '0') == '1':
        return dict(STATS, off=True)
    # (1) the line: R48's lateflag (R48_DIV6ALL) builds the unflagged calls' produce / rebuild per divider scope from
    # dqi_junc_r8p_div6._make(G) at call time and restores retained19's hooks after the call, so wrap _make itself
    # (DQ5's c12, installed before this module, wraps the same function: its produce is wrapped here in turn).
    import importlib
    div6 = sys.modules.get('dqi_junc_r8p_div6') or importlib.import_module('dqi_junc_r8p_div6')
    mk0 = div6._make
    if not getattr(mk0, '_b5_d6k', False):
        def _make_d6k(G):
            p0, r0 = mk0(G)
            STATS['made'] = STATS.get('made', 0) + 1
            return _make(G, p0, r0)
        _make_d6k._b5_d6k = True
        div6._make = _make_d6k
        STATS['make_wrapped'] = 1
    # (2) scopes whose hooks already hold div6's produce / rebuild (div6.install without R48: the bench)
    walk = sys.modules['walk_bound']
    done = set()
    for fn in walk.functions(ns.decode):
        if not isinstance(fn, types.FunctionType):
            continue
        G = fn.__globals__
        if id(G) in done or not fn.__code__.co_filename.endswith('decode_inverse_retained.py'):
            continue
        p0, r0 = G.get('produce_retained'), G.get('rebuild_retained')
        if p0 is None or r0 is None or getattr(p0, '_b5_d6k', False):
            continue
        if not (getattr(p0, '__code__', None) and p0.__code__.co_filename.endswith('dqi_junc_r8p_div6.py')):
            continue                                       # only on div6's produce / rebuild
        done.add(id(G))
        G['produce_retained'], G['rebuild_retained'] = _make(G, p0, r0)
        STATS['scopes'] += 1
    assert STATS.get('make_wrapped') or STATS['scopes'], 'div6 not found (install after the line patches)'
    atexit.register(_dump)
    return dict(scopes=STATS['scopes'], make_wrapped=STATS.get('make_wrapped', 0), ceil=CEIL, extra=EXTRA, term=RTERM)


def _test(argv):
    """Run K1's division bench (tools/dqi_k1_divprof.py) with this module installed first: lateu and latef."""
    import runpy
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    sys.path.insert(0, here)
    import kickmix as km
    sys.modules['circ'] = km
    import dqi_junc_r8_entry as entry
    print('install:', install(km, entry), flush=True)
    sys.argv = ['dqi_k1_divprof.py', argv[0] if argv else 'lateu,latef'] + list(argv[1:])
    runpy.run_path(os.path.join(here, 'dqi_k1_divprof.py'), run_name='__main__')


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'test':
    sys.modules.setdefault('dqi_int_div6k', sys.modules[__name__])
    _test(sys.argv[2:])
