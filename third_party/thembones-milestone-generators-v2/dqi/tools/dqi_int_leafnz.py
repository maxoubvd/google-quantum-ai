#!/usr/bin/env python3
"""dqi_int_leafnz.py -- B5 (22nd holder, Nov 18): the kit's 18-AND GF(64) leaf (pro_jet.inverse6_joint18_12extra) with
its one redundant AND (the rank-one product ('p', 200, 474, 64): a OR b = 1 on every nonzero input, DQ2's census /
R46-DQI-EXT's INV17) written as the affine 1 + a + b (X, CX, CX: no Toffoli), built at generation time: 17 ANDs, exact
on every NONZERO input; its zero-flag output is then identically 0 (F = NOR(a, b) = 1 + a + b + ab) and is freed.
For call sites whose input is proved nonzero (INV17's list A: the junction norms, R46-DQI-EXT).
  make(pj) -> leaf_nz(b, x, with_zero=False) returning (t, None), and cached12_nz(b, x, plan, with_zero=False)
  python dqi_int_leafnz.py test      (tracer: F == 0 on all 64 inputs, t == 1/x on the 63 nonzero ones, 17 ANDs)
"""
import inspect
import sys
import types

OLD = "with ref.expose_pair(b, list(x) + raw, a, c) as (p, q): rankone_xor(b, p, q, out, d)"
NEW = ("with ref.expose_pair(b, list(x) + raw, a, c) as (p, q): "
       "(_b5_affine_xor if (a, c, d) == (200, 474, 64) else rankone_xor)(b, p, q, out, d)")


def _affine_xor(b, p, q, z, vector):
    j = 0
    while vector:
        if vector & 1:
            b.x(z[j])
            b.cx(p, z[j])
            b.cx(q, z[j])
        vector >>= 1
        j += 1


def make(pj):
    fn = pj.inverse6_joint18_12extra
    src = inspect.getsource(fn)
    assert src.count(OLD) == 1, 'leafnz: the rank-one line of inverse6_joint18_12extra not found once'
    gl = dict(fn.__globals__)
    gl['_b5_affine_xor'] = _affine_xor
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + src.replace(OLD, NEW), fn.__code__.co_filename, 'exec'), gl, ns)
    joint_nz = ns['inverse6_joint18_12extra']

    def leaf_nz(b, x, with_zero=False):
        assert not with_zero
        t, F = joint_nz(b, x)
        b.free(F)                                  # identically |0> (see the docstring)
        return t, None
    ci = pj.cached_inverse
    cached12_nz = types.FunctionType(ci.__code__, dict(ci.__globals__, inverse6_fixed18=leaf_nz), 'cached_inverse_nz',
                                     ci.__defaults__, ci.__closure__)
    return leaf_nz, cached12_nz, joint_nz


if __name__ == '__main__' and sys.argv[1:] == ['test']:
    import numpy as np
    import dqi_int_inv6trace as T
    from exp.jet24 import pro_jet as pj
    import kickmix as km
    leaf_nz, cached12_nz, joint_nz = make(pj)
    gf = km.GF2Field(6)
    b = T.TB()
    x = b.alloc_qubits(6)
    for i in range(6):
        x[i].v[:] = (np.arange(64) >> i) & 1
    t, F = joint_nz(b, x)
    val = sum(t[i].v.astype(np.int64) << i for i in range(6))
    want = np.array([0] + [gf.pow(v, 62) for v in range(1, 64)])
    nt = sum(1 for o in b.ops if o[4])
    print('joint_nz: ANDs %d; F identically 0: %s; inverse exact on nonzero inputs: %s (at 0 it reads %d)' % (
        nt, not F.v.any(), bool((val == want)[1:].all()), int(val[0])))
