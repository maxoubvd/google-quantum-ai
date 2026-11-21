"""dqi_junc_r8p_latem.py -- un-synthesis division with the 6-bit subnorm computed AFTER the half scalings (a pure
reordering of the kit's `divide_from_norms` as patched by the build line's r6n, tools/dqi_build_patch_r6n.py).

In r6n's division the divisor's norm N (12) is computed after the 24-bit product, then M = Norm6(N) (6) at once,
then the two 12-bit halves of x are scaled by conj6(N) through a 12-qubit target; M is first used after that loop
(its inverse, then the 6-bit part scalings). So N, M and the target coexist: 30 (+ the divisor's zero flag) at
un-synthesis's peak (909). M depends only on N, which the loop does not change, so computing it after the loop
gives the same gates in another order: the loop holds N + target (24), the part scalings N + M + a 6-bit target.
Install after r6n (it re-derives r6n's two variants with one more source edit and puts its wrapper in r6n's place).
"""
import sys

MBLOCK = ("    ar.linear_inplace(b, N, to)\n    M = b.alloc_qubits(6)\n    norm(b, N, M, k6)\n"
          "    ar.linear_inplace(b, N, frm)\n")
LOOP_END = "        masks.append(m)\n    ar.linear_inplace(b, N, conj)\n"
LATE_M = ((MBLOCK, ""), (LOOP_END, LOOP_END + MBLOCK))


def install(km, module):
    r6n = sys.modules.get('dqi_build_patch_r6n')
    assert r6n is not None, 'install after r6n'
    walk = sys.modules['walk_bound']
    mod = sys.modules['norm_deferred_flag_return']
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    flagged = r6n._variant(mod, src, r6n.FLAGGED + LATE_M, 'divide_flag_late_m')
    unflagged = r6n._variant(mod, src, r6n.UNFLAGGED + LATE_M, 'divide_noflag_late_m')
    w = sys.modules.get('dqi_build_patch_r5w')
    noflag = w.NOFLAG if w is not None else {}
    stats = dict(flagged=0, unflagged=0)

    def divide(b, x, y, fm, before, after, owner, K=None):
        if id(fm) in noflag:
            stats['unflagged'] += 1
            return unflagged(b, x, y, fm, before, after, owner, K)
        stats['flagged'] += 1
        return flagged(b, x, y, fm, before, after, owner, K)
    fns = walk.functions(module.decode)
    scopes = {id(f.__globals__): f.__globals__ for f in fns}
    scopes[id(vars(mod))] = vars(mod)
    n = 0
    for g in scopes.values():
        cur = g.get('divide_from_norms')
        if cur is not None and getattr(cur, '__code__', None) is not None \
                and cur.__code__.co_filename.endswith('dqi_build_patch_r6n.py') and cur.__name__ == 'divide':
            g['divide_from_norms'] = divide
            n += 1
    assert n, "r6n's division wrapper not found"
    import atexit
    atexit.register(lambda: print('dqb_latem', stats, file=sys.stderr))
    return dict(late_m=True, scopes=n)
