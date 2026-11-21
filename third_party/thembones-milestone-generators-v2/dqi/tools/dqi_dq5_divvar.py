"""dqi_dq5_divvar.py -- DQ5 (Nov 18): variant forms of the line's un-synthesis division (x -> x / y in GF(2^24)),
for tools/dqi_dq5_divprof.py.

VARIANTS[name] = dict(edits=fn(flagged) -> tuple of (before, after) source replacements applied after r6n + latem,
                      hooks=fn(G) -> (produce_retained, rebuild_retained) swapped into the divider's globals for the call,
                      or None to keep the scope's own hooks).
  div6 : J's div6 inverse (tools/dqi_junc_r8p_div6.py), as R48's late-flag wrapper puts it into the line's 17 unflagged
         divisions (R48_DIV6ALL=1): the line's unflagged form.
  c12  : div6._make wrapped by tools/dqi_dq5_c12.py (the a14 consuming plan, 20 products, same width).
helpers(G, ar, K) -> extra globals visible to the edited source.
"""
import sys


def _div6_hooks(G):
    import dqi_junc_r8p_div6 as div6
    return div6._make(G)


def _c12_hooks(G):
    import dqi_junc_r8p_div6 as div6
    import dqi_dq5_c12 as c12
    c12._wrap_make(div6)
    return div6._make(G)


def _raw(G):
    """G with the profiler's part wrappers removed (div6k reads the kernels' own globals)."""
    G2 = dict(G)
    for k in ('norm_xor_noscratch', 'mul_xor_noscratch'):
        G2[k] = getattr(G[k], '__wrapped__', G[k])
    return G2


def _c12k_hooks(G):
    """c12's produce under B5's div6k (tools/dqi_int_div6k.py, read only): kept wire + rebuild at 14."""
    import dqi_int_div6k as d6k
    p0, r0 = _c12_hooks(G)
    return d6k._make(_raw(G), p0, r0)


def _d6k_hooks(G):
    import dqi_int_div6k as d6k
    p0, r0 = _div6_hooks(G)
    return d6k._make(_raw(G), p0, r0)


VARIANTS = {
    'div6': dict(edits=lambda flagged: (), hooks=_div6_hooks),
    'c12': dict(edits=lambda flagged: (), hooks=_c12_hooks),
    'd6k': dict(edits=lambda flagged: (), hooks=_d6k_hooks),
    'c12k': dict(edits=lambda flagged: (), hooks=_c12k_hooks),
}


def helpers(G, ar, K):
    return {}
