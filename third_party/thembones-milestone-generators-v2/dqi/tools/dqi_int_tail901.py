"""dqi_int_tail901.py -- B5 (17th holder, Nov 16): the synthesis tail (kit `final_decoder/decode.py` tail_inplace and
tail_finish, both directions) at ONE WIRE LESS, exactly.  A DQI_INT_EXT module: it swaps the code of every live clone
of the two kit functions (as r46_dqi_syn_tail.py does for the correction), keeping each clone's globals.

Why the tail holds 902 [M dqi_int_livedump.py on l9r4o, Nov 16 19:53Z]: base = 840 words + G 24 (the junction chain's
scale, `external_extraction > j_in`) + p 12 + pole = 877; every tail kernel then holds an inverse pair (a 12-wire
inverse + its zero flag, 13) and a 12-wire temporary (masked_update's copy, T4's register + F, deferred_multiply's
product, a zero test's 11-AND ladder + its flag): 877 + 13 + 12 = 902.  In every one of the 13 runs the zero flag of
the pair is idle; it is needed only to adjust / unadjust its words (x <- x + [x == 0]) around them.
The construction (exact on every input):
  * a flag whose word is NOT used adjusted inside the run is X-measured before the run; its adjustment of the inverse
    is never undone: the inverse is erased adjusted, and the mask bit of its constant coordinate is paid with the flag's
    own mask when the pair is consumed (forward divide by p's inverse: 0 Toffolis);
  * a product that needs the ADJUSTED non-inverse word (x * q' with q' = q + [q == 0]) runs before the inverse exists, on
    its own zero test (11): forward v <- v q'^pole; reverse d <- d q'^pole and d <- d p'^(not pole); its phase part runs
    after the inverse exists (the word it acts on is unchanged in between);
  * elsewhere the flag is X-measured and tested again (11) when it is next needed; its mask is paid by a CZ on the
    re-tested flag (same function of the same word).
Moves (exact): v <- v q'^pole runs before q's inverse and before T4 (T4 acts only when NOT pole, the masked update only
when pole); in the reverse the undo of S3 runs before T4's undo for the same reason.
Count [P, per decode]: forward +22 (two zero tests of q), reverse +44 (two of q, one of d, one of p); every moment of the
tail at <= base + 24 = 901.
Switch B5_T901 = comma list of pieces (default 'p,d,f'): p = the divide by p's inverse (forward fold + reverse product
first), d = the reverse divide by d's inverse (product with the flag measured, flag tested again), f = tail_finish (both
directions).  Pieces off = the kit's code for that piece.
"""
import os
import sys

PARTS = frozenset(x for x in os.environ.get('B5_T901', 'p,d,f').split(',') if x)
LOG = dict(tail_inplace=0, tail_finish=0, parts=','.join(sorted(PARTS)))

SRC = '''
def tail_finish(b,p,q,d,u,v,pole,reverse=False,token=None,Z=None):
    """B5 tail901 form of the kit's tail_finish (see dqi_int_tail901.py)."""
    if 'f' not in _B5_T901:
        return _B5_T901_KIT['tail_finish'](b,p,q,d,u,v,pole,reverse,token,Z)
    if not reverse:
        Zq0=guarded.fresh_zero_flag(b,q)
        b.cx(Zq0,q[0]);vm=masked_update(b,v,q,pole);b.cx(Zq0,q[0])
        a=b.alloc_bits(1);b.hmr(Zq0,a[0]);b.free(Zq0)
        qi,Zq=guarded.inverse_zero(b,q)
        b.cz(Zq,a[0]);b.free(a)
        z1=b.alloc_bits(1);b.hmr(Zq,z1[0]);b.free(Zq)
        tt=term_correction_with_inverse(b,p,q,d,u,v,pole,qi,Z=Z)
        Zq=guarded.fresh_zero_flag(b,q)
        b.cz(Zq,z1[0]);b.free(z1)
        b.cx(Zq,qi[0])
        z2=b.alloc_bits(1);b.hmr(Zq,z2[0]);b.free(Zq)
        dm=masked_update(b,d,qi,pole)
        m=b.alloc_bits(12);b.hmr(qi,m);b.free(qi)
        b.cswap(pole,u,v);b.cswap(pole,v,d)
        return dict(inverse=(m,z2),term=tt,v=vm,d=dm)
    b.cswap(pole,v,d);b.cswap(pole,u,v)
    m,z2=token['inverse']
    Zq0=guarded.fresh_zero_flag(b,q)
    b.cx(Zq0,q[0]);fd=masked_update(b,d,q,pole);b.cx(Zq0,q[0])
    a=b.alloc_bits(1);b.hmr(Zq0,a[0]);b.free(Zq0)
    qi,Zq=guarded.inverse_zero(b,q)
    b.cz(Zq,a[0]);b.free(a)
    b.cx(Zq,qi[0])
    b.cz(qi,m);b.free(m);b.cz(Zq,z2[0]);b.free(z2)
    dmask=token['d']
    for x,bit in zip(d,dmask):b.ccz(pole,x,bit)
    b.free(dmask);phase_masked_product(b,d,qi,pole,fd)
    z3=b.alloc_bits(1);b.hmr(Zq,z3[0]);b.free(Zq)
    fv=masked_update(b,v,qi,pole)
    Zq=guarded.fresh_zero_flag(b,q)
    b.cz(Zq,z3[0]);b.free(z3)
    vmask=token['v']
    for x,bit in zip(v,vmask):b.ccz(pole,x,bit)
    b.free(vmask)
    b.cx(Zq,q[0]);phase_masked_product(b,v,q,pole,fv);b.cx(Zq,q[0])
    b.cx(Zq,qi[0])
    z4=b.alloc_bits(1);b.hmr(Zq,z4[0]);b.free(Zq)
    Z=term_correction_with_inverse(b,p,q,d,u,v,pole,qi,True,token['term'])
    m=b.alloc_bits(12);b.hmr(qi,m);b.free(qi)
    _B5_T901_PHASE0(b,q,m,z4[0]);b.free(m);b.free(z4)
    return Z


def tail_inplace(b,p,q,d,u,v,pole,reverse=False,token=None):
    """B5 tail901 form of the kit's tail_inplace: (d,u,v) <-> (delta,N0,N1), with p,q the implicit leaders."""
    if not reverse:
        ar.mul_xor0(b,q,v,p)
        square12(b,q);ar.mul_xor0(b,q,u,d);square12(b,q,True)
        if 'p' in _B5_T901:
            pi,Fp=guarded.inverse_zero(b,p)
            b.cx(Fp,pi[0])
            z1=b.alloc_bits(1);b.hmr(Fp,z1[0]);b.free(Fp)
            b.x(pole);m1=masked_update(b,d,pi,pole);b.x(pole)
            m=b.alloc_bits(12);b.hmr(pi,m);b.free(pi)
            tp=((m,z1),m1)
        else:
            masks=[]
            def divide(pi):
                b.x(pole);masks.append(masked_update(b,d,pi,pole));b.x(pole)
            tp=guarded.inverse_apply(b,p,divide,defer=True);tp=(tp,masks[0])
        ar.mul_xor0(b,q,v,p)
        def multiply(di,mm):
            b.x(pole);vm=masked_update(b,v,d,pole);b.x(pole)
            ar.mul_xor0(b,q,u,v);um=deferred_multiply(b,u,d)
            return vm,um
        Z=guarded.fresh_zero_flag(b,d);b.cx(Z,d[0]);td=multiply(None,None);b.cx(Z,d[0])
        tf=tail_finish(b,p,q,d,u,v,pole,Z=Z)
        return dict(p=tp,d=td,finish=tf)
    Z=tail_finish(b,p,q,d,u,v,pole,True,token['finish'])
    vm,um=token['d']
    di=inverse_initialized(b,d)
    if 'd' in _B5_T901:
        b.cx(Z,di[0])
        zz=b.alloc_bits(1);b.hmr(Z,zz[0]);b.free(Z)
        fu=deferred_multiply(b,u,di)
        ar.mul_xor0(b,q,u,v)
        b.x(pole);fv=masked_update(b,v,di,pole);b.x(pole)
        Z=guarded.fresh_zero_flag(b,d)
        b.cz(Z,zz[0]);b.free(zz)
        b.cx(Z,d[0])
        b.cz(u,um);b.free(um);ar.phase_product(b,u,d,fu);b.free(fu)
        b.x(pole)
        for x,bit in zip(v,vm):b.ccz(pole,x,bit)
        b.free(vm);phase_masked_product(b,v,d,pole,fv)
        b.x(pole)
        b.cx(Z,d[0]);b.cx(Z,di[0])
    else:
        b.cx(Z,d[0]);b.cx(Z,di[0])
        reverse_deferred_multiply(b,u,d,di,um);ar.mul_xor0(b,q,u,v)
        b.x(pole);reverse_masked_multiply(b,v,d,di,pole,vm);b.x(pole)
        b.cx(Z,di[0]);b.cx(Z,d[0])
    guarded.erase(b,d,di,Z,False);ar.mul_xor0(b,q,v,p)
    if 'p' in _B5_T901:
        (m,z1),m1=token['p']
        Fp0=guarded.fresh_zero_flag(b,p)
        b.cx(Fp0,p[0]);b.x(pole);f1=masked_update(b,d,p,pole);b.x(pole);b.cx(Fp0,p[0])
        a=b.alloc_bits(1);b.hmr(Fp0,a[0]);b.free(Fp0)
        pi,Fp=guarded.inverse_zero(b,p)
        b.cz(Fp,a[0]);b.free(a)
        b.cx(Fp,pi[0])
        b.cz(pi,m);b.free(m);b.cz(Fp,z1[0]);b.free(z1)
        b.x(pole)
        for x,bit in zip(d,m1):b.ccz(pole,x,bit)
        b.free(m1);phase_masked_product(b,d,pi,pole,f1)
        b.x(pole)
        b.cx(Fp,pi[0])
        guarded.erase(b,p,pi,Fp,False)
    else:
        def multiply(pi):
            b.x(pole);reverse_masked_multiply(b,d,pi,p,pole,token['p'][1]);b.x(pole)
        guarded.inverse_apply(b,p,multiply,incoming=token['p'][0])
    ar.mul_xor0(b,q,v,p)
    square12(b,q);ar.mul_xor0(b,q,u,d);square12(b,q,True)
'''

NEEDED = ('guarded', 'ar', 'masked_update', 'phase_masked_product', 'deferred_multiply', 'reverse_deferred_multiply',
          'reverse_masked_multiply', 'inverse_initialized', 'square12', 'term_correction_with_inverse')


def _phase0_of(guarded):
    """The zero-aware inverse phase (x, mask, zero_mask) that the line's guarded.erase pays with."""
    er = guarded.erase
    fn = er.__globals__.get('phase_inverse36_zero')
    assert callable(fn), 'guarded.erase has no phase_inverse36_zero in its globals'
    assert 'phase_inverse36_zero' in er.__code__.co_names, 'guarded.erase does not call phase_inverse36_zero'
    return fn


def install(km, ns):
    if not PARTS:
        return dict(LOG, on=False)
    walk = sys.modules['walk_bound']
    scope = {}
    exec(compile(SRC, __file__ + ':tail901', 'exec'), scope)
    done = set()
    for fn in walk.functions(ns.decode):
        name = fn.__name__
        if name not in ('tail_finish', 'tail_inplace') or id(fn) in done:
            continue
        if not fn.__code__.co_filename.endswith('/final_decoder/decode.py'):
            continue
        assert fn.__closure__ is None, name
        gl = fn.__globals__
        missing = [k for k in NEEDED if k not in gl]
        assert not missing, (name, missing)
        if '_B5_T901_KIT' not in gl:
            gl['_B5_T901_KIT'] = {}
        kit = gl['_B5_T901_KIT']
        if name not in kit:
            import types
            kit[name] = types.FunctionType(fn.__code__, gl, name, fn.__defaults__, None)
        gl['_B5_T901'] = PARTS
        _phase0_of(gl['guarded'])                     # checked now, resolved again at every call (late patches)
        gl['_B5_T901_PHASE0'] = (lambda b, x, m, z, _g=gl: _phase0_of(_g['guarded'])(b, x, m, z))
        fn.__code__ = scope[name].__code__
        done.add(id(fn))
        LOG[name] += 1
    assert LOG['tail_inplace'] and LOG['tail_finish'], 'tail901: kit tail functions not found (%r)' % (LOG,)
    import atexit
    atexit.register(lambda: print('dqi_int_tail901', LOG, file=sys.stderr, flush=True))
    return dict(LOG)
