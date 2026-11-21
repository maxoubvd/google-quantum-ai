"""dqi_kern_flagn.py -- the decoder's flagged un-synthesis division with its zero flag from the 12-bit norm
(workstream dqi_kern, builder B6, item D3).

INVALID, do not install (B6, 21:55Z): the decoder's terminated folds have y = 0 with x != 0, and the decoder's
own flagged form then returns x ^ the caller's correction; a flag taken after the product cannot keep x (the
product with y = 0 is 0). Measured with `DQI_KERN_YZERO_X=1 DQI_KERN_BEFORE=5592405 ... dqi_kern_divtest.py ip9g`:
237/256 values, 11 phase flips (ip9f, the decoder's own: 256/256, 0). Kept for the record.

The decoder's flagged form (the kit's divide_from_norms + r6n's FLAGGED edits + latem's late M) tests the 24-bit
divisor for zero BEFORE the product (23 ANDs; r6n moved the norm after the product, so the 12-bit test, 11 ANDs,
was no longer available there). The flag can as well be taken AFTER the product from the norm: when y = 0 the
promise x = 0 makes the product 0, so the caller's before(F) correction (a constant added to words that include
x) lands on the same value as in the kit's order, and every later scaling is by the safe divisor 1. Edits on top
of r6n's unflagged head and latem's late M: N = Norm(y) after the product, F = zero_flag(N), before(F), y and N
made safe there; the old x's mask is paid after y is restored (old x = x_final * y with the unsafe y). Everything
else, the free retirement of F through M t included, is the decoder's own. Width as the decoder's flagged form.
Install after latem (DQI_KERN_PATCHES run last); unflagged folds keep latem's form.
"""
import sys

STATS = dict(flagged=0, unflagged=0)


def edits(r6n, latem):
    head_unflagged = r6n.UNFLAGGED[0]
    flag_after = (r6n.PRODUCT, r6n.LATE_NORM + "    F = zero_flag(b, N)\n    before(F)\n    b.cx(F, y[0])\n"
                  "    b.cx(F, N[0])\n")
    oldpay = ("    g['phase_product24'](b, x, y, oldmask, K)\n    b.free(oldmask)\n\n    b.cx(F, y[0])\n",
              "    b.cx(F, y[0])\n    g['phase_product24'](b, x, y, oldmask, K)\n    b.free(oldmask)\n")
    return (head_unflagged, flag_after, (r6n.SECOND, ""), oldpay) + tuple(latem.LATE_M)


def variant(mod, name='divide_flag_n'):
    r6n = sys.modules['dqi_build_patch_r6n']
    latem = sys.modules['dqi_junc_r8p_latem']
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    return r6n._variant(mod, src, edits(r6n, latem), name)


def install(km, module):
    raise RuntimeError('dqi_kern_flagn is invalid for the decoder (see the module docstring)')
    mod = sys.modules['norm_deferred_flag_return']
    r6n = sys.modules['dqi_build_patch_r6n']
    latem = sys.modules['dqi_junc_r8p_latem']
    text = open(mod.__file__).read()
    start = text.index('def divide_from_norms(')
    src = text[start:text.index('\ndef sweep(', start)].rstrip() + '\n'
    flagged = variant(mod)
    unflagged = r6n._variant(mod, src, r6n.UNFLAGGED + latem.LATE_M, 'divide_noflag_late_m_b6')
    w = sys.modules.get('dqi_build_patch_r5w')
    noflag = w.NOFLAG if w is not None else {}

    def divide(b, x, y, fm, before, after, owner, K=None):
        if id(fm) in noflag:
            STATS['unflagged'] += 1
            return unflagged(b, x, y, fm, before, after, owner, K)
        STATS['flagged'] += 1
        return flagged(b, x, y, fm, before, after, owner, K)
    walk = sys.modules['walk_bound']
    fns = walk.functions(module.decode)
    scopes = {id(f.__globals__): f.__globals__ for f in fns}
    scopes[id(vars(mod))] = vars(mod)
    n = 0
    for g in scopes.values():
        cur = g.get('divide_from_norms')
        if cur is not None and getattr(cur, '__code__', None) is not None \
                and cur.__code__.co_filename.endswith('dqi_junc_r8p_latem.py') and cur.__name__ == 'divide':
            g['divide_from_norms'] = divide
            n += 1
    assert n, "latem's division wrapper not found"
    return dict(kern_flagn=n)
