"""dqi_junc_r8p_tailctl.py -- the synthesis tail's correction without its 12-qubit factor register.

exp/frontier/a07/decode_correction_and.term_correction_with_inverse (the variant in use) forms factor = F AND qi
(12 qubits, F = [d = 0] and not pole, qi = 1/q) and adds p * factor into v and v * factor into u with two plain
products. That factor register (12) and F are the 13 at synthesis's peak (905 in ip3/ip4). Here each product is
the kit's rank-42 mul_xor0 with every one of its CCX terms controlled by F (one temporary AND per term, erased
free): v ^= F (p qi), u ^= F (v qi). +42 Toffolis per product, -12 in the factor's 12 ANDs: +72 per call, two calls
(forward and reverse). Width: one ancilla at a time instead of 12. Installed by swapping the code of every live
clone (it reads only the globals the original reads: guarded, ar).
"""
import sys

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
    if reverse:
        Z = guarded.fresh_zero_flag(b, d)
        b.cz(Z, token[0])
        b.free(token)
    F = b.alloc_qubits(1)[0]
    b.x(pole)
    b.init_and(Z, pole, F)
    b.x(pole)
    cb = _Ctl(b, F)
    if reverse:
        ar.mul_xor0(cb, v, qi, u)
        ar.mul_xor0(cb, p, qi, v)
    else:
        ar.mul_xor0(cb, p, qi, v)
        ar.mul_xor0(cb, v, qi, u)
    b.x(pole)
    b.del_and(Z, pole, F)
    b.x(pole)
    b.free(F)
    if reverse:
        return Z
    zm = b.alloc_bits(1)
    b.hmr(Z, zm[0])
    b.free(Z)
    return zm
'''


def install(km, module):
    walk = sys.modules['walk_bound']
    scope = {}
    exec(compile(SRC, __file__ + ':tail', 'exec'), scope)
    code = scope['term_correction_with_inverse'].__code__
    n = 0
    for f in walk.functions(module.decode):
        if f.__name__ == 'term_correction_with_inverse' and f.__code__.co_filename.endswith('/a07/decode_correction_and.py'):
            assert f.__closure__ is None
            assert 'guarded' in f.__globals__ and 'ar' in f.__globals__
            f.__code__ = code
            n += 1
    assert n, 'term_correction_with_inverse (decode_correction_and) not found'
    return dict(tail_ctl=n)
