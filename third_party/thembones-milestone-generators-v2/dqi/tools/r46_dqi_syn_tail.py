"""r46_dqi_syn_tail.py -- seat R46-DQI-SYN: the synthesis tail's correction WITH its factor register at Q 902.

The correction (kit `exp/frontier/a07/decode_correction_and.term_correction_with_inverse`) adds, under
F = [d = 0] AND NOT pole,  v ^= F p / q  and  u ^= F v / q  with qi = 1 / q live.  The kit forms factor = F AND qi
(12 ANDs) and runs two plain rank-42 products: 1 + 12 + 84 = 97.  J's `dqi_junc_r8p_tailctl.py` (in every recipe of
the line since Nov 12) controls each of the 84 CCX terms by F through a temporary AND instead: 1 + 168 = 169 per
call, because factor (12) beside qi (12), Zq, pole, Z and F is 2 wires over the stage's width
[M: best_q23j5's recipe without tailctl = Q 904, 350,084.4 (-143.8), my build `notc`].

This module keeps the factor and removes the two wires (Z and F) from the products' lifetime, exactly:
  forward  F = Z AND NOT pole; factor = F AND qi; F erased at once (free, Z and pole live); Z X-measured at once (its
           mask is deferred to the return today as well); the two products; factor X-measured with 12 DEFERRED mask
           bits (the return pays them with classically conditioned CZ(F, qi_j) while F and qi are live).   97
  reverse  Z = [d = 0] (11); the forward's Z mask (CZ); F; the forward's factor masks (conditioned CZ, free);
           factor = F AND qi; F erased (free); Z X-measured (mask m2); the two products; Z tested AGAIN (11; d is
           untouched by the products), m2 paid by CZ; F again (1); factor erased (free, F and qi live); F erased;
           Z returned to the caller as before.                                              11 + 1 + 12 + 84 + 11 + 1 = 120
Line today (tailctl): forward 169, reverse 11 + 169 = 180.  Expected: -72 forward, -60 reverse = -132.0 per decode.
During the products the live set is factor 12 + qi 12 + Zq + pole: the same 26 as the tail's masked updates
(`masked_update`: its 12-wire register beside qi, Zq, pole), which fit 902 today.

Install as a DQI_INT_EXT module (token `ext`) on a recipe WITHOUT tailctl (DQI_JUNC_R8P without it); it swaps the code
of every live clone of the a07 function, as tailctl does.  R46_TAIL_MODE=fr (default: both directions), f (forward
only; the reverse then runs tailctl's per-term control and pays the deferred factor masks).

RESULTS [M box 6, Nov 15]: mode fr = 350,097.0 (-131.2), quick 0 / 0, but Q 914: the reverse's SECOND zero test runs
its 10-wire AND ladder beside the factor, and the forward still touches 904 at the instant the factor is formed (Z and
F both live).  Mode t3 = the narrowest order of the full factor: forward F, Z measured, factor, F measured (mask
deferred: the return pays CZ(F)); reverse Z kept, F erased during the products and re-made (1): expected base + 27 in
both directions, i.e. Q 903: ONE wire over the line.  See research/panel/R46-DQI-SYN.md.
Measured after that note: t3 = 350,086.1 at Q 904 (synthesis 903, un-synthesis 904).

Mode t4 (THE ONE THAT FITS 902; code in tools/r46_dqi_syn_t4.py, which must lie beside this file): an 11-wire masked
register walking three hyperplanes of forms: 350,102.3 at Q 902 (-125.9), official 0 / 0 on both salts.  Mode t4m is
t4's negative control (the return leaves the deferred register masks unpaid: 2,580 phase errors in 10,240 shots).
"""
import os
import sys

LOG = dict(swapped=0, mode='')

SRC = '''
def term_correction_with_inverse(b, p, q, d, u, v, pole, qi,
                                 reverse=False, token=None, Z=None):
    class _Ctl:
        def __init__(self, inner, F):
            self.inner, self.F = inner, F

        def __getattr__(self, name):
            return getattr(self.inner, name)

        def ccx(self, a, c, t):
            anc = self.inner.alloc_qubits(1)[0]
            self.inner.init_and(self.F, a, anc)
            self.inner.ccx(anc, c, t)
            self.inner.del_and(self.F, a, anc)
            self.inner.free(anc)

    def make_f(Zq):
        F = b.alloc_qubits(1)[0]
        b.x(pole)
        b.init_and(Zq, pole, F)
        b.x(pole)
        return F

    def drop_f(Zq, F):
        b.x(pole)
        b.del_and(Zq, pole, F)
        b.x(pole)
        b.free(F)

    if _R46_MODE in ('t4', 't4m'):
        if not reverse:
            return _r46.t4_forward(b, ar, p, u, v, pole, qi, Z)
        return _r46.t4_reverse(b, ar, guarded, p, d, u, v, pole, qi, token, _R46_MODE == 't4m')
    if not reverse and _R46_MODE == 't3':
        F = make_f(Z)
        zm = b.alloc_bits(1)
        b.hmr(Z, zm[0])
        b.free(Z)
        factor = b.alloc_qubits(12)
        b.init_and(F, qi, factor)
        ff = b.alloc_bits(1)
        b.hmr(F, ff[0])
        b.free(F)
        ar.mul_xor0(b, p, factor, v)
        ar.mul_xor0(b, v, factor, u)
        fm = b.alloc_bits(12)
        b.hmr(factor, fm)
        b.free(factor)
        return (zm, fm, ff)
    if not reverse:
        F = make_f(Z)
        factor = b.alloc_qubits(12)
        b.init_and(F, qi, factor)
        drop_f(Z, F)
        zm = b.alloc_bits(1)
        b.hmr(Z, zm[0])
        b.free(Z)
        ar.mul_xor0(b, p, factor, v)
        ar.mul_xor0(b, v, factor, u)
        fm = b.alloc_bits(12)
        b.hmr(factor, fm)
        b.free(factor)
        return (zm, fm)
    zm, fm = token[0], token[1]
    Z = guarded.fresh_zero_flag(b, d)
    b.cz(Z, zm[0])
    b.free(zm)
    F = make_f(Z)
    if len(token) == 3:
        b.cz(F, token[2][0])
        b.free(token[2])
    for c, bit in zip(qi, fm):
        b.ccz(F, c, bit)
    b.free(fm)
    if _R46_MODE == 't3':
        factor = b.alloc_qubits(12)
        b.init_and(F, qi, factor)
        drop_f(Z, F)
        ar.mul_xor0(b, v, factor, u)
        ar.mul_xor0(b, p, factor, v)
        F = make_f(Z)
        b.del_and(F, qi, factor)
        b.free(factor)
        drop_f(Z, F)
        return Z
    if _R46_MODE == 'f':
        cb = _Ctl(b, F)
        ar.mul_xor0(cb, v, qi, u)
        ar.mul_xor0(cb, p, qi, v)
        drop_f(Z, F)
        return Z
    factor = b.alloc_qubits(12)
    b.init_and(F, qi, factor)
    drop_f(Z, F)
    m2 = b.alloc_bits(1)
    b.hmr(Z, m2[0])
    b.free(Z)
    ar.mul_xor0(b, v, factor, u)
    ar.mul_xor0(b, p, factor, v)
    Z = guarded.fresh_zero_flag(b, d)
    b.cz(Z, m2[0])
    b.free(m2)
    F = make_f(Z)
    b.del_and(F, qi, factor)
    b.free(factor)
    drop_f(Z, F)
    return Z
'''


def install(km, module):
    mode = os.environ.get('R46_TAIL_MODE', 'fr')
    assert mode in ('fr', 'f', 't3', 't4', 't4m'), mode
    walk = sys.modules['walk_bound']
    t4 = None
    if mode in ('t4', 't4m'):
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            'r46_dqi_syn_t4', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r46_dqi_syn_t4.py'))
        t4 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(t4)
        LOG['t4'] = t4.LOG
    n = 0
    for f in walk.functions(module.decode):
        if f.__name__ == 'term_correction_with_inverse' and f.__code__.co_filename.endswith('/a07/decode_correction_and.py'):
            assert f.__closure__ is None
            assert 'guarded' in f.__globals__ and 'ar' in f.__globals__
            scope = {}
            exec(compile(SRC, __file__ + ':tail', 'exec'), scope)
            f.__globals__['_R46_MODE'] = mode
            f.__globals__['_r46'] = t4
            f.__code__ = scope['term_correction_with_inverse'].__code__
            n += 1
    assert n, 'term_correction_with_inverse (a07/decode_correction_and) not found: is tailctl still in DQI_JUNC_R8P?'
    LOG.update(swapped=n, mode=mode)
    import atexit
    atexit.register(lambda: print('r46_dqi_syn_tail', LOG, file=sys.stderr))
    return dict(LOG)
