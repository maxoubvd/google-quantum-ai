"""r46_dqi_ext_fgret.py -- R46-DQI-EXT: the first homogeneous group's LAST return step measured out (-99 Toffolis).

In the first homogeneous group the chain starts without a scale: after the free first step the N register holds the
public constant (q0, 0) (the kit sets D = 1 and scales it by q0). The forward's step with head 1 computes
N = (q0, 0) + H_1 D (99 Toffolis, needed). The return's mirror step runs N ^= H_1 D again (99 Toffolis) only to
bring N back to that constant: N is then a constant plus a bilinear form of two live words, so it is erased by an
X-measurement and the kit's free bilinear phase (exp.arithmetic.division.phase_product24) instead.

    return step, as built:   scale N by 1/q; swap N, D; scale L by 1/q;  N ^= H_1 D            (99 Toffolis)
    return step, here:       scale N by 1/q; swap N, D; scale L by 1/q;  N ^= (q0, 0);  m = HMR(N);
                             phase (-1)^<m, H_1 D>;  N ^= (q0, 0)                                (0 Toffolis)

Install LAST among the patches that touch chain_step (token `ext`, DQI_INT_EXT=...,r46_dqi_ext_fgret). Every
reachable `chain_step` that has the prefix patch's N-only hook (`_dqb_mul_head` in its globals) gets a guard in front;
the guard fires for (reverse, head 1, N in the first homogeneous group's words) and nowhere else. The count of
firings is printed at exit (it must be 1 per decode built) and returned in the log.
"""
import atexit
import os
import sys
import types

STATS = dict(installed=0, fired=0, calls=0)
GROUP = int(os.environ.get('R46_FGRET_GROUP', '3'))        # the first homogeneous group of the line's recipe

GUARD_SRC = '''
def chain_step(b, head, N, D, L, q, ct, c, reverse=False, first=False):
    if reverse and not first and _r46_fg_hit(head, N):
        return _r46_fg_step(b, head, N, D, L, q, ct, core, _dqb_mul_head)
    return _r46_fg_orig(b, head, N, D, L, q, ct, c, reverse, first)
'''


def _hit(head, N):
    STATS['calls'] += 1
    return head[0][0].id // 24 == 1 and N[0].id // 48 == GROUP


def _bits(v):
    i = 0
    while v:
        if v & 1:
            yield i
        v >>= 1
        i += 1


def _step(b, head, N, D, L, q, ct, core, mul_head):
    G = mul_head.__globals__
    ks, groups, K = G['ks'], G['groups'], G['K']
    kap = groups[GROUP]['kap']
    q0 = ks[0] ^ kap
    assert q == ks[1] ^ kap and q0 == q, (q, q0, ks[0], ks[1])     # heads 0 and 1 share their node: q0 = q1
    phase24 = sys.modules['exp.arithmetic.division'].phase_product24
    core.scalepair(b, N, core.inv(q))
    b.swap(N, D)
    core.scalepair(b, L, core.inv(q))
    a, h = head
    core.cx(b, h, a, ct)
    H = core.combine(head)
    low = [i for i in _bits(q0)]
    for i in low:
        b.x(N[i])                                   # N = H_1 D
    m = b.alloc_bits(24)
    b.hmr(N, m)                                     # N = 0, mask m on its old value
    phase24(b, H, D, m, K)                          # (-1)^<m, H_1 D>: classically conditioned CZ, no Toffoli
    b.free(m)
    for i in low:
        b.x(N[i])                                   # N = (q0, 0), what the built step leaves
    core.cx(b, h, a, ct)
    STATS['fired'] += 1


def install(km, module):
    walk = sys.modules['walk_bound']
    seen = set()
    for f in walk.functions(module.decode):
        if f.__name__ != 'chain_step' or '_dqb_mul_head' not in f.__globals__ or id(f) in seen:
            continue
        seen.add(id(f))
        assert f.__closure__ is None
        orig = types.FunctionType(f.__code__, f.__globals__, 'chain_step_r46_orig', f.__defaults__)
        # one original per function object: the name in globals is unique per function
        name = '_r46_fg_orig_%d' % len(seen)
        f.__globals__[name] = orig
        f.__globals__['_r46_fg_hit'] = _hit
        f.__globals__['_r46_fg_step'] = _step
        src = GUARD_SRC.replace('_r46_fg_orig', name)
        sc = {}
        exec(compile(src, __file__ + ':guard', 'exec'), sc)
        f.__code__ = sc['chain_step'].__code__
        STATS['installed'] += 1
    assert STATS['installed'], 'no chain_step with the N-only hook found'
    return dict(installed=STATS['installed'], group=GROUP)


@atexit.register
def _report():
    print('r46_dqi_ext_fgret: guards installed %d, chain_step calls seen %d, measured-out return steps %d'
          % (STATS['installed'], STATS['calls'], STATS['fired']), flush=True)
