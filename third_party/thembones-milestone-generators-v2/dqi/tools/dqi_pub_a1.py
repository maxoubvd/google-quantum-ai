"""dqi_pub_a1.py -- builder D1, item a1: the un-synthesis body in CRT form (one product fewer per body).

Mechanism credit: the synthesis closing of the public DQI decoder (Syed Asad Ali, Apache-2.0) keeps the full values
of the new polynomial at r points through its division, so the old polynomial costs r products; our own code here.

The live un-synthesis sweep (`exp/round7/root/norm_deferred_flag_return.sweep`, with J's divider patches in its
globals) runs per body: divide_from_norms(x = U, y = W, before = completion, after = body), where the completion adds
F lead co (F = [W = 0]) to [U, T] and the body quotient_add adds Qt floor(l S_c / q) to T in len(S) rank-one products
(S, T q-adic digit words, q = x^2 + x + k, l = the quotient left in U).
Since l S_c = q floor(l S_c / q) + U_pre (U_pre = U before the division), floor(l S_c / q) = (l S_c - U_pre) / q
exactly, and its len(T) = len(S) - 1 values at t_j = 0..len(T)-1 determine it:
  premap   T += Qt sum_j U_pre(t_j) / q(t_j) L_j              free, at the end of `before`
  body     T += Qt sum_j l(t_j) S_c(t_j) / q(t_j) L_j          len(T) rank-one products (kit's rank_affine)
Zero source (F = 1, the division is the identity): the raw residue is h W = 0, so premap(U_pre) is wrong exactly by
F lead premap(co_U); that is cancelled by adding premap(co_U) to the completion's T-constants (no gate).
Forms checked against the floor map by tools/dqi_pub_a1model.py. install(km, ns) with ns.decode; DQI_PUB_A1=0 = off.
"""
import os
import sys
from functools import lru_cache

TOOLS = os.path.dirname(os.path.abspath(__file__))
if TOOLS not in sys.path:
    sys.path.insert(0, TOOLS)
import dqi_pub_a1model as MODEL                                       # noqa: E402  (pure-python GF(2^12) forms)

LOG = dict(bodies=0, products_saved=0)


@lru_cache(None)
def forms(m, n, k):
    return MODEL.crt_forms(m, n, k)


SWEEP_SRC = '''
def sweep(b, C, B, folds, divide):
    previous.qadic(b, C, B, k=previous.K)
    beta = b.alloc_qubits(12)
    b.x(beta[0])
    leaders = [C[-2], C[-1], beta]
    H = b.alloc_qubits(24)
    h0 = [H[:12], H[12:]]
    b.free(H)
    A = h0 + list(C[:-2])
    heads = previous.prior.rr.tail.heads_from_banks(A, B)
    frames, _, _ = previous.frames_for(A, B)
    for j, (a, d, k) in enumerate(frames):
        n = len(a)
        a1, a2, beta = leaders
        U = previous.combine(a[:2])
        W = previous.combine(d[:2])
        ft = list(a[2:]) + [a1, a2]
        fs = list(d) + [beta]
        st = list(d[2:]) + [beta]
        ss = list(a[2:]) + [a1, a2]
        if j:
            previous.qframe(b, a, d, leaders, frames[j-1][2], True)
            previous.qframe(b, a, d, leaders, k)
            fm1 = folds[2*j][1]
            late1 = _D1A1_H2 and not _D1A1_UNFLAGGED(fm1)
            def correct_first(F):
                if not late1:
                    _D1A1_FIX(b, F, a[:2], ft, k)
                target = a + [a1, a2]
                previous.controlled_linear_xor(b, F, a2, target, _D1A1_CO(k, n + 1, len(target), len(ft)))
                if late1:
                    _D1A1_PRE(b, a[:2], ft, k)
            def first():
                _D1A1_BODY(b, a[:2], fs, ft, k)
            # unflagged folds (r6n never calls `before` there): the premap runs on the raw residue first
            if not late1:
                _D1A1_PRE(b, a[:2], ft, k)
            divide_from_norms(b, U, W, fm1, correct_first, first, divide, K=k)
        source = previous.combine(a[2:4]) if n > 3 else previous.combine([a[-1], a1])
        fm2 = folds[2*j+1][1]
        late2 = _D1A1_H2 and not _D1A1_UNFLAGGED(fm2)
        def correct_second(F):
            if not late2:
                _D1A1_FIX(b, F, d[:2], st, k)
            target = d + [beta]
            previous.controlled_linear_xor(b, F, beta, target, _D1A1_CO(k, n, len(target), len(st)))
            if late2:
                _D1A1_PRE(b, d[:2], st, k)
        def second():
            _D1A1_BODY(b, d[:2], ss, st, k)
        if not late2:
            _D1A1_PRE(b, d[:2], st, k)
        divide_from_norms(b, W, source, fm2, correct_second, second, divide, K=k)
    return A, heads, leaders
'''


def helpers():
    sn = next(m for m in list(sys.modules.values())
              if 'rank_affine' in getattr(m, '__dict__', {}) and 'projected_product' in m.__dict__)
    rank_affine = sn.rank_affine
    common = next(m for m in list(sys.modules.values())
                  if 'linear_words' in getattr(m, '__dict__', {}) and 'matrix_ops' in m.__dict__)
    cx = common.cx
    pr = next(m for m in list(sys.modules.values())
              if 'constant_qadic' in getattr(m, '__dict__', {}) and 'complete_zero' in m.__dict__)

    @lru_cache(None)
    def co(k, degree, length, n_target):
        c = list(pr.constant_qadic(k, degree, length))
        if os.environ.get('DQI_PUB_A1_CO', 'orig') == 'orig':
            return tuple(c)             # the completion restores the true residue h W (0 at W = 0): no change
        pre, _ = forms(n_target + 1, n_target, k)
        off = length - n_target
        for i in range(n_target):
            c[off + i] ^= MODEL.mul(pre[i][0], c[0]) ^ MODEL.mul(pre[i][1], c[1])
        return tuple(c)

    def premap(b, Uw, T, k):
        pre, _ = forms(len(T) + 1, len(T), k)
        for i, (c0, c1) in enumerate(pre):
            if c0:
                cx(b, Uw[0], T[i], c0)
            if c1:
                cx(b, Uw[1], T[i], c1)

    ctl = next(m for m in list(sys.modules.values())
               if 'controlled_linear_xor' in getattr(m, '__dict__', {}) and 'divide_safe138' in m.__dict__)

    def fix(b, F, Uw, T, k):
        """T -= F premap(U): at a zero source (F = 1) the raw residue is not h W, and the premap must vanish there.
        Two controlled linear xors (12 Toffolis each, no scratch); called first in the divider's `before`."""
        pre, _ = forms(len(T) + 1, len(T), k)
        for w in (0, 1):
            coeffs = [pre[i][w] for i in range(len(T))]
            if any(coeffs):
                ctl.controlled_linear_xor(b, F, Uw[w], T, coeffs)
        LOG['fixes'] = LOG.get('fixes', 0) + 1

    def body(b, head, source, target, k):
        _, ups = forms(len(source), len(target), k)
        for ac, bc, tc in ups:
            rank_affine(b, head, ac, source, bc, target, tc)
        LOG['bodies'] += 1
        LOG['products_saved'] += len(source) - len(target)

    return co, premap, body, fix


def install(km, ns):
    if os.environ.get('DQI_PUB_A1', '1') != '1':
        return dict(a1='off')
    walk = sys.modules['walk_bound']
    co, premap, body, fix = helpers()
    made = {}
    done = 0
    for fn in walk.functions(ns.decode):
        g = fn.__globals__
        s = g.get('sweep')
        if fn.__name__ != 'consume_synthesis' or s is None or not hasattr(s, '__code__'):
            continue
        if not s.__code__.co_filename.endswith('norm_deferred_flag_return.py'):
            continue
        if id(s) not in made:
            env = dict(s.__globals__)
            w = sys.modules.get('dqi_build_patch_r5w')
            noflag = w.NOFLAG if w is not None else {}
            env.update(_D1A1_CO=co, _D1A1_PRE=premap, _D1A1_BODY=body, _D1A1_FIX=fix,
                       _D1A1_H2=os.environ.get('DQI_PUB_A1_H2', '1') == '1',
                       _D1A1_UNFLAGGED=lambda fm, noflag=noflag: id(fm) in noflag)
            exec(compile(SWEEP_SRC, __file__, 'exec'), env)
            made[id(s)] = env['sweep']
        g['sweep'] = made[id(s)]
        done += 1
    assert done, 'no live consume_synthesis (norm_deferred_flag_return.sweep) found'
    LOG.update(installed=done, sweeps=len(made))
    return dict(LOG)
