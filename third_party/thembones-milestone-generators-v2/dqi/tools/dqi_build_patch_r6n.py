"""Round 6 patch "N": narrower un-synthesis division (workstream dqi_build).

The kit's `divide_from_norms` computes the divisor's 12-bit norm BEFORE the 24-qubit product and keeps it alive
across the product only to retire its zero flag (measured, re-created afterwards from the 6-bit subnorm). That
norm register is what puts un-synthesis at 914 qubits on the T1 line. This patch computes the norm AFTER the
product (the divisor is not changed by the product):
  flagged folds:   zero flag straight from the 24-bit divisor (23 ANDs instead of 11), kept alive as one qubit
                   across the product, no re-creation from the subnorm (-5): +7 Toffolis, 11 qubits narrower;
  unflagged folds (policy W, patch r5w): norm after the product, no Toffoli change, 12 qubits narrower.
Install after r5w (works without it: then every fold takes the flagged form). DQI_R6N=0 switches it off.
The kit's `early` mask bit is no longer measured, so it is stored as 0 explicitly: allocated bits are not clean,
and leaving it unset gave phase errors on half the shots that have a terminated fold.
"""
import os
import sys

ON = os.environ.get('DQI_R6N', '1') != '0'
STATS = dict(flagged=0, unflagged=0)

HEAD = ("    N = b.alloc_qubits(12)\n    norm(b, y, N, K)\n    F = zero_flag(b, N)\n    before(F)\n"
        "    b.cx(F, y[0])\n    early = b.alloc_bits(1)\n    b.hmr(F, early[0])\n    b.free(F)\n")
PRODUCT = "    oldmask = g['multiply_forward99'](b, x, y, K)\n    b.cx(y[12:], y[:12])\n"
SECOND = "    F = zero_flag(b, M)\n    b.cx(F, N[0])\n    b.cx(F, M[0])\n"
UNPATCH = "    b.cx(F, y[0])\n    b.cx(F, N[0])\n    b.cx(F, M[0])\n    after()\n"
RETIRE = "    b.hmr(F, fmask[0])\n    b.free(F)\n"
LATE_NORM = PRODUCT + "    N = b.alloc_qubits(12)\n    norm(b, y, N, K)\n"

FLAGGED = (
    (HEAD, "    F = zero_flag(b, y)\n    before(F)\n    b.cx(F, y[0])\n    early = b.alloc_bits(1)\n"
           "    b.bit_store(early, False)\n"),
    (PRODUCT, LATE_NORM),
    (SECOND, ""),
)
UNFLAGGED = (
    (HEAD, "    early = b.alloc_bits(1)\n    b.bit_store(early, False)\n"),
    (PRODUCT, LATE_NORM),
    (SECOND, ""),
    (UNPATCH, "    after()\n"),
    (RETIRE, ""),
)


def _variant(mod, src, edits, name):
    for before, after in edits:
        assert src.count(before) == 1, (name, before)
        src = src.replace(before, after)
    src = src.replace('def divide_from_norms(', 'def %s(' % name)
    scope = dict(vars(mod))
    exec(compile(src, __file__ + ':' + name, 'exec'), scope)
    return scope[name]


def install(km, module):
    if not ON:
        return dict(on=False)
    walk = sys.modules['walk_bound']
    fns = walk.functions(module.decode)
    mod = sys.modules['norm_deferred_flag_return']
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    flagged = _variant(mod, src, FLAGGED, 'divide_flag_late')
    unflagged = _variant(mod, src, UNFLAGGED, 'divide_noflag_late')
    w = sys.modules.get('dqi_build_patch_r5w')
    noflag = w.NOFLAG if w is not None else {}
    current = vars(mod)['divide_from_norms']

    def divide(b, x, y, fm, before, after, owner, K=None):
        if id(fm) in noflag:
            STATS['unflagged'] += 1
            return unflagged(b, x, y, fm, before, after, owner, K)
        STATS['flagged'] += 1
        return flagged(b, x, y, fm, before, after, owner, K)
    scopes = {id(f.__globals__): f.__globals__ for f in fns}
    scopes[id(vars(mod))] = vars(mod)
    n = 0
    for g in scopes.values():
        if g.get('divide_from_norms') is current:
            g['divide_from_norms'] = divide
            n += 1
    assert n, 'divide_from_norms not found'
    import atexit
    atexit.register(lambda: print('dqb_r6n', STATS, file=sys.stderr))
    return dict(on=True, scopes=n, over_r5w=w is not None)
