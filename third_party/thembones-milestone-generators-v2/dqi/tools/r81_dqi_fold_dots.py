#!/usr/bin/env python3
"""r81_dqi_fold_dots.py -- design seat R81-DQI-FOLD, round 2 (Nov 18): THE FUSED DOTS of the converter's kit levels,
as ONE module in B5's DQI_INT_EXT form (install(km, ns); token r81_dqi_fold_dots, LAST).

What it replaces (and nothing else). The kit's post_look (exp/r-main/narrow240/measured_heads.py) makes, at the
levels below the direct threshold (m = 2 .. 5, n = m - 1 = 1 .. 4), in convert and in unconvert,
    def dots():
     for v in payloads:
      for i in range(n):ar.mul_xor0(b,g[i],v[i],v[-1])
= two dot products of ONE n-word vector g with the two payloads: u[-1] ^= sum g[i] u[i], w[-1] ^= sum g[i] w[i]:
2 n K-products of 42 ANDs, consecutive, operands read-only, targets disjoint from every operand. For n = 2, 3, 4 the
bilinear map "two dots of one vector" has schemes of 154 / 224 / 308 rank-one terms over GF(2) against 168 / 252 /
336 (R81-DQI-FOLD round 1: K(2,2) at 22 and K(2,3) at 32 ANDs over GF(8), lifted through the 7 points of
P^1(GF(8)); checked by the designer and, with independent code, by the re-pricer). Each term (A, B, C) is ONE CCX in
the kit's own in-place form (field.mul_xor0): CNOT fan-in of the bits of A onto one wire of g, of B onto one wire of
(u, w), CNOT fan-out of C from one target wire, CCX, undo. 0 work qubits. n = 1 stays the kit's loop (one scalar
times two words is exactly 12 a place).
Saving [P]: 14 + 28 + 28 = 70 Toffolis a stage, 140 a decode.

Exactness (the lane's rule: exact by form or on all inputs).
 1. At install every embedded scheme is checked as a TENSOR IDENTITY against the definition on all basis pairs (a
    bilinear map is fixed by them: this is all inputs), with km.GF2Field(12) agreeing with this file's multiplier
    on the 144 basis pairs; a scheme that fails is not used (assert).
 2. The swap touches a post_look only where the BYTECODE of its nested dots equals the bytecode of the loop above
    compiled by the same interpreter (code, constants, names, free variables); any other dots is refused and
    listed in LOG['refused'].
 3. At a call the fused path is taken only for two payloads, n in the installed sizes, twelve wires a word, and
    2 + 3 n words on pairwise distinct wires; otherwise it returns False and the kit's loop runs unchanged.
 4. `python -I r81_dqi_fold_dots.py test`: the emitter's gate list is simulated SYMBOLICALLY (every wire a GF(2)
    form of the inputs plus products of two forms): operands restored, targets ^= the definition's tensor, for any
    initial target: exact on all inputs; CCX counts 154 / 224 / 308; negative controls.
 5. `selftest` / `convtest` (in the line's build environment, see the usage): the kit's own closure machinery and
    simulator; the converter alone before and after.

Environment: R81_DQI_FOLD_DOTS_OFF=1 installs nothing; R81_DQI_FOLD_DOTS_N=2 (or 2,3 / 2,3,4: default) limits the
sizes (the checkpoint n = 2 alone is 28 Toffolis a decode); R81_DQI_FOLD_DOTS_OUT=<json> receives LOG after every call.
Public mathematics, credited: flip graphs (M. Kauers, J. Moosbauer, arXiv 2212.01175) and plus transitions (Y. Arai,
Y. Ichikawa, K. Hukushima, arXiv 2312.17094) through DQ6's walker (the GF(8) schemes); Toom-Cook evaluation /
interpolation (the lift); the kit's in-place rank-one product form (field.mul_xor0).
Usage:
  python -I tools/r81_dqi_fold_dots.py test
  ~/qcc/.venv/bin/python -B tools/r81_dqi_fold_dots.py kittest <kit copy, e.g. ~/qcc/ws/dqi_zn/src/kv__noa16_cap900_noret>
  bash tools/sh/dqi_zn_lenv.sh <core> ~/qcc/.venv/bin/python -B r81_dqi_fold_dots.py selftest
  bash tools/sh/dqi_zn_lenv.sh <core> ~/qcc/.venv/bin/python -B r81_dqi_fold_dots.py convtest <batches> <seed> [kinds]
"""
import gc
import json
import os
import sys
import types

MOD12 = 0x1053
NAME = 'r81_dqi_fold_dots'
OFF = os.environ.get('R81_DQI_FOLD_DOTS_OFF', '0') == '1'
SIZES = tuple(int(x) for x in os.environ.get('R81_DQI_FOLD_DOTS_N', '2,3,4').split(',') if x)
OUT = os.environ.get('R81_DQI_FOLD_DOTS_OUT')
COUNT = {2: (154, 168), 3: (224, 252), 4: (308, 336)}          # (fused, the kit's 2 n x 42)
LOG = dict(installed=False, off=OFF, sizes=list(SIZES), codes=0, functions=0, already=0, swapped=[], refused=[],
           fused={}, plain={}, sites={}, saved=0)
TERMS = {}

REF_SRC = ("def _tpl(b, g, n, payloads):\n"
           " def dots():\n"
           "  for v in payloads:\n"
           "   for i in range(n):ar.mul_xor0(b,g[i],v[i],v[-1])\n"
           " return dots\n")
NEW_SRC = ("def _tpl(b, g, n, payloads):\n"
           " def dots():\n"
           "  if __import__('" + NAME + "').fused_dots(b,g,payloads,n):return\n"
           "  for v in payloads:\n"
           "   for i in range(n):ar.mul_xor0(b,g[i],v[i],v[-1])\n"
           " return dots\n")


def gmul(a, b, n=12, mod=MOD12):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a >> n:
            a ^= mod
    return r


def _bits(x):
    while x:
        t = x & -x
        yield t.bit_length() - 1
        x ^= t


def expected(n, i, j):
    """The definition on a basis pair: operand a bit i = x^la of g[ja]; operand b bit j = x^lb of word jb of payload
    k (bit 12 (k n + jb) + lb); the output has the product in word k (bits 12 k ..) when ja == jb, else nothing."""
    ja, la = divmod(i, 12)
    k, r = divmod(j, 12 * n)
    jb, lb = divmod(r, 12)
    return (gmul(1 << la, 1 << lb) << (12 * k)) if ja == jb else 0


def check_terms(n, terms):
    """Exact: the scheme's tensor equals the definition's on every basis pair (hence on all inputs)."""
    na, nb = 12 * n, 24 * n
    T = [[0] * nb for _ in range(na)]
    for a, bm, c in terms:
        if not (a and bm and c) or a >> na or bm >> nb or c >> 24:
            return False, ('term', a, bm, c)
        ib = list(_bits(bm))
        for i in _bits(a):
            row = T[i]
            for j in ib:
                row[j] ^= c
    for i in range(na):
        for j in range(nb):
            if T[i][j] != expected(n, i, j):
                return False, (i, j)
    return True, None


def _parse(n):
    lines = [ln.split() for ln in _DATA[n].splitlines() if ln.strip()]
    na, nb, nc, r = (int(v) for v in lines[0])
    assert (na, nb, nc) == (12 * n, 24 * n, 24) and r == COUNT[n][0] == len(lines) - 1, 'header of n = %d' % n
    return [tuple(int(v, 16) for v in ln[:3]) for ln in lines[1:]]


def _load():
    for n in (2, 3, 4):
        if n not in TERMS:
            terms = _parse(n)
            ok, bad = check_terms(n, terms)
            assert ok, 'scheme n = %d NOT exact at %s' % (n, bad)
            TERMS[n] = terms
    return TERMS


def emit(b, A, B, C, terms):
    """One CCX a term, the kit's in-place form (field.mul_xor0) with two different masks and flat wire lists."""
    for ma, mb, mc in terms:
        ia, ib, ic = list(_bits(ma)), list(_bits(mb)), list(_bits(mc))
        pa, pb, pc = A[ia[0]], B[ib[0]], C[ic[0]]
        for i in ia[1:]:
            b.cx(A[i], pa)
        for i in ib[1:]:
            b.cx(B[i], pb)
        for j in ic[1:]:
            b.cx(pc, C[j])
        b.ccx(pa, pb, pc)
        for j in reversed(ic[1:]):
            b.cx(pc, C[j])
        for i in reversed(ib[1:]):
            b.cx(B[i], pb)
        for i in reversed(ia[1:]):
            b.cx(A[i], pa)


def _note(kind, n):
    d = LOG[kind]
    d[n] = d.get(n, 0) + 1
    if OUT:
        try:
            with open(OUT, 'w') as fh:
                json.dump(LOG, fh)
        except OSError:
            pass


def fused_dots(b, g, payloads, n):
    """True when the two dots were emitted here; False = the caller runs the kit's loop (nothing was emitted)."""
    if OFF or n not in SIZES or n not in COUNT:
        _note('plain', n if isinstance(n, int) else -1)
        return False
    try:
        if len(payloads) != 2:
            _note('plain', -2)
            return False
        gw = [g[i] for i in range(n)]
        vw = [[v[i] for i in range(n)] for v in payloads]
        tw = [v[-1] for v in payloads]
        words = gw + vw[0] + vw[1] + tw
        if any(len(w) != 12 for w in words):
            _note('plain', -3)
            return False
        A = [w[l] for w in gw for l in range(12)]
        B = [w[l] for k in range(2) for w in vw[k] for l in range(12)]
        C = [w[l] for w in tw for l in range(12)]
        names = [str(q) for q in A + B + C]
    except (TypeError, IndexError):
        _note('plain', -4)
        return False
    if len(set(names)) != len(names):
        _note('plain', -5)
        return False
    emit(b, A, B, C, _load()[n])
    LOG['saved'] += COUNT[n][1] - COUNT[n][0]
    co = sys._getframe(1).f_code
    site = '%s:%d n=%d' % (os.path.basename(co.co_filename), co.co_firstlineno, n)
    LOG['sites'][site] = LOG['sites'].get(site, 0) + 1
    _note('fused', n)
    return True


def _nested(code, name='dots'):
    for c in code.co_consts:
        if isinstance(c, types.CodeType) and c.co_name == name:
            return c
    return None


def _tpl(src):
    """The nested dots code of a template, compiled but never executed."""
    co = compile(src, '<' + NAME + '>', 'exec')
    return _nested([c for c in co.co_consts if isinstance(c, types.CodeType) and c.co_name == '_tpl'][0])


def _variants():
    """(reference loop, fused loop) as this interpreter compiles them: with `ar` a plain global (the clone scopes of
    exp/r-main/narrow240/measured_heads.py) and with `ar` a module-level import (the kit's safe_scalar/converter.py:
    from Python 3.12 the compiler then loads the attribute without the method form). Same loop text in both."""
    return [(_tpl(pre + REF_SRC), _tpl(pre + NEW_SRC)) for pre in ('', 'import sys as ar\n')]


def _same(a, b):
    return all(getattr(a, k) == getattr(b, k) for k in
               ('co_code', 'co_consts', 'co_names', 'co_varnames', 'co_freevars', 'co_cellvars', 'co_argcount'))


def install(km, ns):
    if OFF:
        return dict(LOG)
    f = km.GF2Field(12)
    assert f.modulus == MOD12, 'field modulus %#x' % f.modulus
    for i in range(12):
        for j in range(12):
            assert f.mul(1 << i, 1 << j) == gmul(1 << i, 1 << j), 'multiplier differs from the kit field'
    _load()
    variants = _variants()
    done = {}
    for fn in [o for o in gc.get_objects() if isinstance(o, types.FunctionType)]:
        co = fn.__code__
        if co.co_name != 'post_look':
            continue
        if id(co) not in done:
            old = _nested(co)
            new = None
            where = '%s:%d' % (os.path.basename(co.co_filename), co.co_firstlineno)
            if old is None:
                pass                                              # a wrapper named post_look: nothing to swap
            elif any(_same(old, nt) for _, nt in variants):
                LOG['already'] += 1
            else:
                hit = [nt for ref, nt in variants if _same(old, ref)]
                if not hit:
                    LOG['refused'].append(where)
                else:
                    assert hit[0].co_freevars == old.co_freevars
                    kw = dict(co_filename=old.co_filename, co_firstlineno=old.co_firstlineno)
                    if hasattr(old, 'co_qualname'):
                        kw['co_qualname'] = old.co_qualname
                    nd = hit[0].replace(**kw)
                    new = co.replace(co_consts=tuple(nd if c is old else c for c in co.co_consts))
                    LOG['codes'] += 1
                    LOG['swapped'].append(where)
            done[id(co)] = (co, new)                              # co kept alive: its id stays unique
        new = done[id(co)][1]
        if new is not None:
            fn.__code__ = new
            LOG['functions'] += 1
    LOG['installed'] = LOG['functions'] > 0
    sys.stderr.write(NAME + ': %s\n' % json.dumps(LOG))
    return dict(LOG)


# ---------------------------------------------------------------------------------------------------------------
# tests


class _Rec:
    def __init__(self):
        self.ops = []

    def cx(self, c, t):
        self.ops.append((c, t))

    def ccx(self, a, b, t):
        self.ops.append((a, b, t))


def _layout(n):
    na, nb = 12 * n, 24 * n
    g = [[12 * j + l for l in range(12)] for j in range(n)]
    pay = [[[na + 12 * (k * n + j) + l for l in range(12)] for j in range(n)]
           + [[na + nb + 12 * k + l for l in range(12)]] for k in range(2)]
    return g, pay


def _symbolic(n, terms=None):
    """Symbolic run of the emitted gate list. Wire w starts as the form x_w. A CX adds forms and product sets; a CCX
    needs two purely linear controls, the first inside operand a, the second inside operand b, and adds their
    product (a set of basis pairs) to the target. Returns (ccx, cx) after asserting: every wire's linear part is its
    own input bit, operand wires carry no product, and each target wire carries exactly the definition's pairs."""
    na, nb = 12 * n, 24 * n
    N = na + nb + 24
    g, pay = _layout(n)
    rec = _Rec()
    if terms is None:
        assert fused_dots(rec, g, pay, n), 'fused path refused n = %d' % n
    else:
        emit(rec, [q for w in g for q in w], [q for k in range(2) for w in pay[k][:n] for q in w],
             [q for k in range(2) for q in pay[k][-1]], terms)
    lin = [1 << w for w in range(N)]
    quad = [0] * N
    ccx = cx = 0
    for op in rec.ops:
        if len(op) == 2:
            c, t = op
            lin[t] ^= lin[c]
            quad[t] ^= quad[c]
            cx += 1
        else:
            a, bb, t = op
            assert quad[a] == 0 and quad[bb] == 0, 'a control carries a product'
            la, lb = lin[a], lin[bb]
            assert la and la >> na == 0, 'first control is not a form of g'
            assert lb and lb & ((1 << na) - 1) == 0 and lb >> (na + nb) == 0, 'second control is not a form of (u, w)'
            lbm = lb >> na
            q = 0
            for i in _bits(la):
                q ^= lbm << (i * nb)
            quad[t] ^= q
            ccx += 1
    bad = []
    for w in range(N):
        if lin[w] != 1 << w:
            bad.append(('linear', w))
    for w in range(na + nb):
        if quad[w]:
            bad.append(('operand product', w))
    for k in range(2):
        for l in range(12):
            want = 0
            for i in range(na):
                for j in range(nb):
                    if expected(n, i, j) >> (12 * k + l) & 1:
                        want |= 1 << (i * nb + j)
            if quad[na + nb + 12 * k + l] != want:
                bad.append(('target', k, l))
    return ccx, cx, bad


def _test():
    import hashlib
    ok = True
    _load()
    here = os.path.dirname(os.path.abspath(__file__))
    for n, fname in ((2, 'r81_dqi_fold_dots2x2_K154.txt'), (3, 'r81_dqi_fold_dots2x3_K224.txt'),
                     (4, 'r81_dqi_fold_dots2x4_K308.txt')):
        p = os.path.join(here, fname)
        same = 'file absent'
        if os.path.exists(p):
            with open(p) as fh:
                lines = [ln.split() for ln in fh if ln.strip()]
            same = 'EQUAL to ' + fname if [tuple(int(v, 16) for v in ln[:3]) for ln in lines[1:]] == TERMS[n] \
                else 'DIFFERS from ' + fname
            ok &= same.startswith('EQUAL')
        ccx, cx, bad = _symbolic(n)
        good = not bad and ccx == COUNT[n][0]
        ok &= good
        print('n = %d: embedded scheme %s; tensor identity on all %d x %d basis pairs PASS; symbolic circuit: %d CCX, '
              '%d CX, 0 allocations, operands restored, targets ^= definition: %s (kit loop %d: saves %d)' % (
                  n, same, 12 * n, 24 * n, ccx, cx, 'PASS' if good else 'FAIL %s' % bad[:4], COUNT[n][1],
                  COUNT[n][1] - ccx))
    # negative controls: the checks must see one flipped bit, in each of the three factors
    seen = 0
    for n in (2, 3, 4):
        for pos, which in ((0, 0), (len(TERMS[n]) // 2, 1), (len(TERMS[n]) - 1, 2)):
            t = list(TERMS[n][pos])
            t[which] ^= 1 << 5
            terms = TERMS[n][:pos] + [tuple(t)] + TERMS[n][pos + 1:]
            c1 = check_terms(n, terms)[0]
            try:
                c2 = not _symbolic(n, terms)[2]
            except AssertionError:
                c2 = False
            seen += (not c1) and (not c2)
    print('negative controls (one bit flipped in a term): caught by the tensor check AND by the symbolic run: %d of 9'
          % seen)
    ok &= seen == 9
    # the guards: n = 1, three payloads, a shared wire, a short word -> the kit's loop (False, nothing emitted)
    rec = _Rec()
    g, pay = _layout(2)
    g1, pay1 = _layout(1)
    guards = [fused_dots(rec, g1, pay1, 1), fused_dots(rec, g, pay + [pay[0]], 2),
              fused_dots(rec, [g[0], g[0]], pay, 2), fused_dots(rec, [g[0][:11], g[1]], pay, 2),
              fused_dots(rec, g, [pay[0], pay[0]], 2)]
    print('guards (n = 1, three payloads, a repeated word, an 11-wire word, the same payload twice): fused %s, gates '
          'emitted %d' % (guards, len(rec.ops)))
    ok &= not any(guards) and not rec.ops
    # the swap on a copy of the kit's post_look text (no kit needed): bytecode equality, patched closure
    src = ("def post_look(b,f,g,payloads,H,parent,reverse=False,incoming=None,defer=False,phase=None):\n"
           " m=len(f)-1;n=m-1;a=g[-1];q0=f[-2];q1=f[-1]\n"
           " def dots():\n"
           "  for v in payloads:\n"
           "   for i in range(n):ar.mul_xor0(b,g[i],v[i],v[-1])\n"
           " dots()\n")
    calls = []
    arns = types.SimpleNamespace(mul_xor0=lambda b, x, y, z: calls.append((x, y, z)))
    fakemod = types.ModuleType('r81_fake_ar')
    fakemod.mul_xor0 = arns.mul_xor0
    sys.modules['r81_fake_ar'] = fakemod
    scope = dict(ar=arns)
    exec(compile(src, 'kit_copy', 'exec'), scope)
    pl = scope['post_look']
    clone = types.FunctionType(pl.__code__, dict(scope), 'post_look', pl.__defaults__, pl.__closure__)
    scope2 = {}
    exec(compile('import r81_fake_ar as ar\n' + src, 'kit_copy_import', 'exec'), scope2)
    pli = scope2['post_look']
    other = {}
    exec(compile(src.replace('v[i],v[-1]', 'v[-1],v[i]'), 'not_the_kit', 'exec'), dict(scope), other)
    keep = other['post_look']
    fake = types.SimpleNamespace(GF2Field=lambda n: types.SimpleNamespace(modulus=MOD12, mul=gmul))
    sys.modules.setdefault(NAME, sys.modules[__name__])
    var = _variants()
    print('bytecode of the two reference variants (plain global / module-level import) differs on this Python %s: %s'
          % (sys.version.split()[0], not _same(var[0][0], var[1][0])))
    before = dict(LOG, fused=dict(LOG['fused']), plain=dict(LOG['plain']))
    log = install(fake, None)
    for fn, nn in ((pl, 3), (clone, 2), (pli, 4), (pli, 2), (keep, 2), (pl, 1), (pli, 1)):
        rec = _Rec()
        g, pay = _layout(nn)
        f = [None] * (nn + 2)
        del calls[:]
        fn(rec, f, g, pay, None, None)
        print('  %s n = %d: %d CCX emitted by the fused path, %d calls of the loop mul_xor0' % (
            fn.__code__.co_filename, nn, sum(1 for op in rec.ops if len(op) == 3), len(calls)))
        want = (COUNT[nn][0], 0) if (fn is not keep and nn > 1) else (0, 2 * nn)
        ok &= (sum(1 for op in rec.ops if len(op) == 3), len(calls)) == want
    print('swap on text copies of the kit post_look: codes %d, functions %d, refused %s' % (
        log['codes'] - before['codes'], log['functions'] - before['functions'], log['refused']))
    ok &= log['codes'] - before['codes'] == 2 and log['functions'] - before['functions'] == 3 \
        and any('not_the_kit' in r for r in log['refused'])
    again = install(fake, None)
    print('second install: codes %d more, already %d' % (again['codes'] - log['codes'], again['already']))
    ok &= again['codes'] == log['codes'] and again['already'] == 2
    print('sha256 of the embedded schemes: %s' % hashlib.sha256(repr(sorted(TERMS.items())).encode()).hexdigest()[:16])
    print('TEST %s' % ('PASS' if ok else 'FAIL'))
    return ok


def _entry():
    """The line's build environment (tools/sh/dqi_zn_lenv.sh): kickmix, S3's entry, B5's converter harness."""
    import importlib.util
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import dqi_int_convtest as CT                        # read-only; sets sys.modules['circ'] = kickmix
    spec = importlib.util.spec_from_file_location('dqi_sw_entry', os.path.join(here, 'dqi_sw_entry.py'))
    mod = importlib.util.module_from_spec(spec)
    sys.modules['dqi_sw_entry'] = mod
    spec.loader.exec_module(mod)
    return CT, mod


def _kit_clones():
    out = []
    for fn in [o for o in gc.get_objects() if isinstance(o, types.FunctionType)]:
        if fn.__code__.co_name == 'post_look' and _nested(fn.__code__) is not None and 'ar' in fn.__globals__:
            out.append(fn)
    return out


def _closure_circuit(km, code, gl, n):
    b = km.CircuitBuilder()
    g = [b.alloc_qubits(12) for _ in range(n)]
    pay = [[b.alloc_qubits(12) for _ in range(n + 1)] for _ in range(2)]
    cells = dict(b=b, g=g, n=n, payloads=pay)
    fn = types.FunctionType(code, gl, 'dots', None, tuple(types.CellType(cells[k]) for k in code.co_freevars))
    fn()
    return b.finish_circuit(), g, pay


def _np_bits(np, vals):
    bits = (vals[:, :, None] >> np.arange(12)) & 1 == 1
    return bits.transpose(1, 2, 0).reshape(-1, vals.shape[0])


def _np_words(np, bits, nwords):
    b = np.asarray(bits).reshape(nwords, 12, -1).astype(np.int64)
    return (b * (1 << np.arange(12))[None, :, None]).sum(axis=1).T


def _selftest(km, ns=None, shots=256, seed=81):
    """T2: for every distinct post_look code with a nested dots: the kit's own closure (the nested dots, the clone's
    globals, its ar.mul_xor0) on kickmix before and after install(), n = 1 .. 4, values on random shots against
    this file's arithmetic, Toffoli-like counts of the built circuits."""
    import numpy as np
    groups = {}
    for fn in _kit_clones():
        co = fn.__code__
        if id(co) not in groups:
            groups[id(co)] = [fn, _nested(co), 0, co]
        groups[id(co)][2] += 1
    assert groups, 'no post_look with a nested dots in this environment'
    variants = _variants()
    for fn, old, cnt, co in groups.values():
        print('post_look %s:%d: %d function objects; nested dots equals the reference loop (plain, import): %s' % (
            os.path.basename(co.co_filename), co.co_firstlineno, cnt, [_same(old, ref) for ref, _ in variants]),
            flush=True)
    log = install(km, ns)
    print('install:', json.dumps(log), flush=True)
    rng = np.random.default_rng(seed)
    ok = True
    for fn, old, cnt, co in groups.values():
        new = _nested(fn.__code__)
        swapped = new is not old
        gl = fn.__globals__
        tot = [0, 0]
        for n in (1, 2, 3, 4):
            res = {}
            vals = rng.integers(0, 4096, size=(shots, 3 * n + 2), dtype=np.int64)
            for tag, code in (('kit', old), ('new', new)):
                c, g, pay = _closure_circuit(km, code, gl, n)
                words = g + pay[0][:n] + pay[1][:n] + [pay[0][-1], pay[1][-1]]
                qa = km.array([q for w in words for q in w])
                sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed)
                sim.use_same_registers_as(c)
                sim.write_across_shots(qa, _np_bits(np, vals))
                sim.do(c)
                res[tag] = (_np_words(np, sim.read_across_shots(qa), 3 * n + 2),
                            np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool),
                            c.operation_statistics().num_toffoli_like, c.num_qubits,
                            sim.operation_statistics().num_toffoli_like / shots)
            want = vals.copy()
            for k in range(2):
                for i in range(n):
                    want[:, 3 * n + k] ^= np.array([gmul(int(a), int(bb)) for a, bb in
                                                    zip(vals[:, i], vals[:, n + k * n + i])], dtype=np.int64)
            e_kit = int((res['kit'][0] != want).any(axis=1).sum())
            e_new = int((res['new'][0] != want).any(axis=1).sum())
            ph = int(res['new'][1].sum()) + int(res['kit'][1].sum())
            tk, tn = res['kit'][2], res['new'][2]
            tot[0] += tk
            tot[1] += tn
            expect = COUNT[n][0] if (swapped and n in COUNT and n in SIZES) else 84 * n
            good = e_kit == 0 and e_new == 0 and ph == 0 and tn == expect and tk == 84 * n \
                and res['new'][3] == res['kit'][3] == 12 * (3 * n + 2)
            ok &= good
            print('  n = %d: Toffoli-like kit %d, new %d (executed a shot %.1f / %.1f); qubits %d / %d; wrong shots of '
                  '%d: kit %d, new %d; phase flips %d: %s' % (n, tk, tn, res['kit'][4], res['new'][4], res['kit'][3],
                                                             res['new'][3], shots, e_kit, e_new, ph,
                                                             'PASS' if good else 'FAIL'), flush=True)
        print('  %s:%d swapped %s: a stage (n = 1 .. 4): kit %d, new %d: %d fewer' % (
            os.path.basename(co.co_filename), co.co_firstlineno, swapped, tot[0], tot[1], tot[0] - tot[1]), flush=True)
        ok &= swapped and tot[0] - tot[1] == sum(COUNT[n][1] - COUNT[n][0] for n in SIZES if n in COUNT)
    print('SELFTEST %s' % ('PASS' if ok else 'FAIL'), flush=True)
    return ok


def _kittest(root):
    """T2 without the entry's eight minutes: kickmix, the kit copy's own field.py (its mul_xor0) and the REAL
    bytecode of the kit copy's two post_look definitions (compiled from the files with their own imports in view,
    never executed), as functions whose `ar` is that field module."""
    import importlib.util
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    spec = importlib.util.spec_from_file_location('r81_kit_field',
                                                  os.path.join(root, 'playbook/kit/decoder/exp/arithmetic/field.py'))
    field = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(field)
    held = []
    for rel in ('playbook/kit/decoder/exp/safe_scalar/converter.py', 'exp/r-main/narrow240/measured_heads.py'):
        p = os.path.join(root, rel)
        with open(p) as fh:
            co = compile(fh.read(), p, 'exec')
        code = [c for c in co.co_consts if isinstance(c, types.CodeType) and c.co_name == 'post_look'][0]
        held.append(types.FunctionType(code, {'ar': field, '__name__': 'r81_kit_copy'}, 'post_look'))
    print('kit copy %s: field.tensor(12) has %d terms' % (root, len(field.tensor(12))), flush=True)
    return _selftest(km, None)


def _convtest(batches, seed, kinds):
    """T3: B5's converter harness (dqi_int_convtest: build, run, charts) with this module installed last: forward
    words against the unpatched converter's, round trip clean, Toffolis, qubits."""
    import numpy as np
    CT, mod = _entry()
    km = CT.km
    conv = mod.decode.__globals__['g']['convert_safe']
    saved = CT.hankel_saved(mod.decode)
    CT.TR['at'] = 10 ** 9
    circ = {}
    for mode in ('fwd', 'round'):
        circ[('ref', mode)] = CT.build(conv, mode, saved)
        print('built ref', mode, flush=True)
    ok = _selftest(km, types.SimpleNamespace(decode=mod.decode))
    base = dict(fused=dict(LOG['fused']), plain=dict(LOG['plain']), saved=LOG['saved'])
    conv = mod.decode.__globals__['g']['convert_safe']
    for mode in ('fwd', 'round'):
        circ[('new', mode)] = CT.build(conv, mode, saved)
        print('built new %s; fused calls so far %s, loop calls %s, sites %s' % (
            mode, {k: v - base['fused'].get(k, 0) for k, v in LOG['fused'].items()},
            {k: v - base['plain'].get(k, 0) for k, v in LOG['plain'].items()}, LOG['sites']), flush=True)
    stat = {}
    for key, (c, _, _, pk, _) in circ.items():
        stat[key] = c.operation_statistics().num_toffoli_like
        print('%-3s %-5s qubits %d, live peak %d, static Toffoli-like %d' % (key[0], key[1], c.num_qubits, pk,
                                                                             stat[key]), flush=True)
    print('static Toffoli-like, new - ref: fwd %d, round %d; module count of Toffolis saved in the four builds: %d' % (
        stat[('new', 'fwd')] - stat[('ref', 'fwd')], stat[('new', 'round')] - stat[('ref', 'round')],
        LOG['saved'] - base['saved']), flush=True)
    rng = np.random.default_rng(seed)
    tot = dict(shots=0, fwd_equal=0, ref_clean=0, new_clean=0)
    tof = {}
    for bi in range(batches):
        An, Bn, label = CT.charts(kinds[bi % len(kinds)], rng)
        out = {}
        for key, (c, A, B, _, _) in circ.items():
            out[key] = CT.run(c, A, B, An, Bn, seed * 1000 + bi)
            tof[key] = tof.get(key, 0) + out[key][4]
        rf, nf = out[('ref', 'fwd')], out[('new', 'fwd')]
        eq = (rf[0] == nf[0]).all(axis=1) & (rf[1] == nf[1]).all(axis=1)
        clean = {}
        for who in ('ref', 'new'):
            r = out[(who, 'round')]
            clean[who] = (r[0] == An).all(axis=1) & (r[1] == Bn).all(axis=1) & ~r[2] & ~r[3]
        tot['shots'] += len(eq)
        tot['fwd_equal'] += int(eq.sum())
        tot['ref_clean'] += int(clean['ref'].sum())
        tot['new_clean'] += int(clean['new'].sum())
        print('batch %d (%s): forward equal %d/%d, round trip clean ref %d new %d; executed Toffolis a shot: '
              'fwd %.1f -> %.1f, round %.1f -> %.1f' % (
                  bi, label, eq.sum(), len(eq), clean['ref'].sum(), clean['new'].sum(), rf[4], nf[4],
                  out[('ref', 'round')][4], out[('new', 'round')][4]), flush=True)
    print('TOTAL', tot)
    ex = {'%s %s' % k: round(v / batches, 1) for k, v in tof.items()}
    print('executed Toffolis per shot:', ex)
    d_f = stat[('new', 'fwd')] - stat[('ref', 'fwd')]
    d_r = stat[('new', 'round')] - stat[('ref', 'round')]
    ok &= tot['fwd_equal'] == tot['shots'] and tot['new_clean'] == tot['ref_clean'] and d_f == -70 and d_r == -140 \
        and circ[('new', 'round')][0].num_qubits == circ[('ref', 'round')][0].num_qubits
    print('CONVTEST %s (forward equal on every shot; new clean = ref clean; static -70 / -140; qubits equal)' % (
        'PASS' if ok else 'FAIL'), flush=True)
    return ok


def _main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'test'
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    import importlib
    me = importlib.import_module(NAME)                   # one module object: the swapped code calls it by name
    if cmd == 'test':
        sys.exit(0 if me._test() else 1)
    if cmd == 'selftest':
        CT, mod = me._entry()
        sys.exit(0 if me._selftest(CT.km, types.SimpleNamespace(decode=mod.decode)) else 1)
    if cmd == 'kittest':
        sys.exit(0 if me._kittest(os.path.expanduser(sys.argv[2])) else 1)
    if cmd == 'convtest':
        kinds = [int(x) for x in (sys.argv[4] if len(sys.argv) > 4 else '0,1,2,3,4').split(',')]
        sys.exit(0 if me._convtest(int(sys.argv[2]), int(sys.argv[3]), kinds) else 1)
    raise SystemExit('usage: test | kittest <kit copy> | selftest | convtest <batches> <seed> [kinds]')


_DATA = {
    2: '''
24 48 24 154
7e07e0 8b5c6b 87e
6a1110 3318b5f5ac6b 8a000
141 f5a7efc6b7ef 88b000
1417b1 7ef7ef8b58b5 8a88b
6f06a1 331331 88b08a
7b16a1 331f5a 7f
110110 4de0008b5000 8a801
141051 8b57ef7efb84 87e8f4
6a16a1 331f5ac6b7ef 88b0f5
6a17e0 7eff5af5a4de 8010f5
7b16f0 331b847efc6b 8a0f5
7b17b1 7efb84000000 7f0f5
51051 c6b4de7efb84 87e07f
7e0051 c6b7ef7ef4de 8a8f4
517b1 331331f5af5a f588b
6f0141 c6b7efb848b5 8f4
516f0 4de4de000000 8f40f5
517b1 4de4def5af5a f587e
7b16a1 4deb84000000 f5000
110000 7ef7ef7efb84 87e88b
51141 8b57ef331331 88b8f4
51051 c6b4def5ac6b 8a07f
260260 86c426 f58
d6f250 eef86cac9426 a6000
f0f ac92a54262a5 f9b000
f0ff3f 2a52a586c86c a6f9b
30d6f eefeef f9b0a6
f3fd6f eefac9 65
250250 c4a00086c000 a6f3d
f0fd5f 86c2a52a5683 f58ffe
d6fd6f eefac94262a5 f9b0c3
d6f260 2a5ac9ac9c4a f3d0c3
f3f030 eef6832a5426 a60c3
f3ff3f 2a5683000000 650c3
d5fd5f 426c4a2a5683 f58065
260d5f 4262a52a5c4a a6ffe
d5ff3f eefeefac9ac9 c3f9b
30f0f 4262a568386c ffe
d5f030 c4ac4a000000 ffe0c3
d5ff3f c4ac4aac9ac9 c3f58
f3fd6f c4a683000000 c3000
250000 2a52a52a5683 f58f9b
d5ff0f 86c2a5eefeef f9bffe
d5fd5f 426c4aac9426 a6065
664664 48a5e8 b87
c5d8ea fcc48aa245e8 561000
a39 a24eae5e8eae c8b000
a394b7 eaeeae48a48a 561c8b
e8ec5d fccfcc c8b561
4b7c5d fcca24 26d
8ea8ea 16200048a000 5619ea
a392d3 48aeaeeaeb46 b87ee6
c5dc5d fcca245e8eae c8b70c
c5d664 eaea24a24162 9ea70c
4b7e8e fccb46eae5e8 56170c
4b74b7 eaeb46000000 26d70c
2d32d3 5e8162eaeb46 b8726d
6642d3 5e8eaeeae162 561ee6
2d34b7 fccfcca24a24 70cc8b
e8ea39 5e8eaeb4648a ee6
2d3e8e 162162000000 ee670c
2d34b7 162162a24a24 70cb87
4b7c5d 162b46000000 70c000
8ea000 eaeeaeeaeb46 b87c8b
2d3a39 48aeaefccfcc c8bee6
2d32d3 5e8162a245e8 56126d
9c09c fb3fc7 de4
fbb872 40afb3bcdfc7 4f6000
f27 bcd47efc747e a8f000
f277c9 47e47efb3fb3 4f6a8f
8eefbb 40a40a a8f4f6
7c9fbb 40abcd 39d
872872 74000fb3000 4f6e79
f27755 fb347e47ebb9 de4912
fbbfbb 40abcdfc747e a8f76b
fbb09c 47ebcdbcd074 e7976b
7c98ee 40abb947efc7 4f676b
7c97c9 47ebb9000000 39d76b
755755 fc707447ebb9 de439d
9c755 fc747e47e074 4f6912
7557c9 40a40abcdbcd 76ba8f
8eef27 fc747ebb9fb3 912
7558ee 74074000000 91276b
7557c9 74074bcdbcd 76bde4
7c9fbb 74bb9000000 76b000
872000 47e47e47ebb9 de4a8f
755f27 fb347e40a40a a8f912
755755 fc7074bcdfc7 4f639d
1fe1fe eed703 ed4
77334 56feed26c703 2c5000
189 26cc81703c81 bf1000
189343 c81c81eedeed 2c5bf1
2ca077 56f56f bf12c5
343077 56f26c 7e0
334334 9ee000eed000 2c5934
1892bd eedc81c81b82 ed4c11
77077 56f26c703c81 bf1525
771fe c8126c26c9ee 934525
3432ca 56fb82c81703 2c5525
343343 c81b82000000 7e0525
2bd2bd 7039eec81b82 ed47e0
1fe2bd 703c81c819ee 2c5c11
2bd343 56f56f26c26c 525bf1
2ca189 703c81b82eed c11
2bd2ca 9ee9ee000000 c11525
2bd343 9ee9ee26c26c 525ed4
343077 9eeb82000000 525000
334000 c81c81c81b82 ed4bf1
2bd189 eedc8156f56f bf1c11
2bd2bd 7039ee26c703 2c57e0
f8ef8e cdfd61 95d
e98fc f9cdfd98d61 3fd000
f67 d98147d61147 da4000
f67815 147147cdfcdf 3fdda4
7720e9 f90f9 da43fd
8150e9 f9d98 704
8fc8fc 1be000cdf000 3fde59
f6779b cdf147147c26 95daa0
e90e9 f9d98d61147 da44f9
e9f8e 147d98d981be e594f9
815772 f9c26147d61 3fd4f9
815815 147c26000000 7044f9
79b79b d611be147c26 95d704
f8e79b d611471471be 3fdaa0
79b815 f90f9d98d98 4f9da4
772f67 d61147c26cdf aa0
79b772 1be1be000000 aa04f9
79b815 1be1bed98d98 4f995d
8150e9 1bec26000000 4f9000
8fc000 147147147c26 95dda4
79bf67 cdf1470f90f9 da4aa0
79b79b d611bed98d61 3fd704
a00a0 198fad 114
888bf0 a6f1985c2fad fe000
828 5c245afad45a 1bb000
828378 45a45a198198 fe1bb
b50888 a6fa6f 1bb0fe
378888 a6f5c2 51
bf0bf0 e35000198000 fe145
8283d8 19845a45abf7 1141ea
888888 a6f5c2fad45a 1bb0af
8880a0 45a5c25c2e35 1450af
378b50 a6fbf745afad fe0af
378378 45abf7000000 510af
3d83d8 fade3545abf7 114051
a03d8 fad45a45ae35 fe1ea
3d8378 a6fa6f5c25c2 af1bb
b50828 fad45abf7198 1ea
3d8b50 e35e35000000 1ea0af
3d8378 e35e355c25c2 af114
378888 e35bf7000000 af000
bf0000 45a45a45abf7 1141bb
3d8828 19845aa6fa6f 1bb1ea
3d83d8 fade355c2fad fe051
''',
    3: '''
36 72 24 224
7e07e0 8b5c6b 87e
517b16f0 4dec6b7eff5ab84331 8a0f5
110000 7ef7ef0007efb84 87e88b
110000141 f5a7ef000c6b7ef 88b000
141051051 8b5f5ac6b4de331331 8f4000
1106a17e0 4deb848b5b848b5000 8010f5
7b10516f0 7ef8b5c6bf5ab84331 7f8f4
51051 8b53318b54de4de8b5 87e07f
1410517b1 7efb84b84331b84b84 f588b
1107b16f0 4def5ac6bf5ac6bf5a 8a0f5
6a1110110 4de0000008b5000 8a801
7e07e0051 4de000000b84000000 8a8f4
7e0051141 8b57ef000331331 88b8f4
6f00516f0 3314deb844dec6b7ef 88b801
517b1 7ef8b58b53314de4de f5
7e0000000 4dec6b7eff5ab84331 87e8f4
7b16a1110 3318b5000f5ac6b 8a000
510516f0 7efc6b8b5f5ab84331 8f40f5
7b1000000 7efc6bc6b000000000 80108a
6f07b16a1 b843313317ef8b58b5 f507f
141051 8b57ef0007efb84 87e8f4
7b16a1 f5af5a4dec6bc6b331 7f
6a16f0110 4dec6b7efb847ef4de 88b87e
110000000 7efb84b844de331331 7f88b
6f0051051 8b53318b54dec6bf5a 8a07f
6a17b16a1 b847ef8b57ef8b58b5 f5000
511417b1 4deb84000b84f5ac6b 8a88b
1417b17b1 7efb84000000000 7f0f5
7b16f06a1 331331 88b08a
1100517b1 7efc6bc6b331b84b84 f587e
1416a16a1 331f5a000c6b7ef 88b0f5
6f0141 c6b7ef000b848b5 8f4
260260 86c426 f58
d5ff3f030 c4a4262a5ac9683eef a60c3
250000 2a52a50002a5683 f58f9b
250000f0f ac92a50004262a5 f9b000
f0fd5fd5f 86cac9426c4aeefeef ffe000
250d6f260 c4a68386c68386c000 f3d0c3
f3fd5f030 2a586c426ac9683eef 65ffe
d5fd5f 86ceef86cc4ac4a86c f58065
f0fd5ff3f 2a5683683eef683683 c3f9b
250f3f030 c4aac9426ac9426ac9 a60c3
d6f250250 c4a00000086c000 a6f3d
260260d5f c4a000000683000000 a6ffe
260d5ff0f 86c2a5000eefeef f9bffe
30d5f030 eefc4a683c4a4262a5 f9bf3d
d5ff3f 2a586c86ceefc4ac4a c3
260000000 c4a4262a5ac9683eef f58ffe
f3fd6f250 eef86c000ac9426 a6000
d5fd5f030 2a542686cac9683eef ffe0c3
f3f000000 2a5426426000000000 f3d0a6
30f3fd6f 683eefeef2a586c86c c3065
f0fd5f 86c2a50002a5683 f58ffe
f3fd6f ac9ac9c4a426426eef 65
d6f030250 c4a4262a56832a5c4a f9bf58
250000000 2a5683683c4aeefeef 65f9b
30d5fd5f 86ceef86cc4a426ac9 a6065
d6ff3fd6f 6832a586c2a586c86c c3000
d5ff0ff3f c4a683000683ac9426 a6f9b
f0ff3ff3f 2a5683000000000 650c3
f3f030d6f eefeef f9b0a6
250d5ff3f 2a5426426eef683683 c3f58
f0fd6fd6f eefac90004262a5 f9b0c3
30f0f 4262a500068386c ffe
664664 48a5e8 b87
2d34b7e8e 1625e8eaea24b46fcc 56170c
8ea000 eaeeae000eaeb46 b87c8b
8ea000a39 a24eae0005e8eae c8b000
a392d32d3 48aa245e8162fccfcc ee6000
8eac5d664 162b4648ab4648a000 9ea70c
4b72d3e8e eae48a5e8a24b46fcc 26dee6
2d32d3 48afcc48a16216248a b8726d
a392d34b7 eaeb46b46fccb46b46 70cc8b
8ea4b7e8e 162a245e8a245e8a24 56170c
c5d8ea8ea 16200000048a000 5619ea
6646642d3 162000000b46000000 561ee6
6642d3a39 48aeae000fccfcc c8bee6
e8e2d3e8e fcc162b461625e8eae c8b9ea
2d34b7 eae48a48afcc162162 70c
664000000 1625e8eaea24b46fcc b87ee6
4b7c5d8ea fcc48a000a245e8 561000
2d32d3e8e eae5e848aa24b46fcc ee670c
4b7000000 eae5e85e8000000000 9ea561
e8e4b7c5d b46fccfcceae48a48a 70c26d
a392d3 48aeae000eaeb46 b87ee6
4b7c5d a24a241625e85e8fcc 26d
c5de8e8ea 1625e8eaeb46eae162 c8bb87
8ea000000 eaeb46b46162fccfcc 26dc8b
e8e2d32d3 48afcc48a1625e8a24 56126d
c5d4b7c5d b46eae48aeae48a48a 70c000
2d3a394b7 162b46000b46a245e8 561c8b
a394b74b7 eaeb46000000000 26d70c
4b7e8ec5d fccfcc c8b561
8ea2d34b7 eae5e85e8fccb46b46 70cb87
a39c5dc5d fcca240005e8eae c8b70c
e8ea39 5e8eae000b4648a ee6
9c09c fb3fc7 de4
7557c98ee 74fc747ebcdbb940a 4f676b
872000 47e47e00047ebb9 de4a8f
872000f27 bcd47e000fc747e a8f000
f27755755 fb3bcdfc707440a40a 912000
872fbb09c 74bb9fb3bb9fb3000 e7976b
7c97558ee 47efb3fc7bcdbb940a 39d912
755755 fb340afb3074074fb3 de439d
f277557c9 47ebb9bb940abb9bb9 76ba8f
8727c98ee 74bcdfc7bcdfc7bcd 4f676b
fbb872872 74000000fb3000 4f6e79
9c09c755 74000000bb9000000 4f6912
9c755f27 fb347e00040a40a a8f912
8ee7558ee 40a074bb9074fc747e a8fe79
7557c9 47efb3fb340a074074 76b
9c000000 74fc747ebcdbb940a de4912
7c9fbb872 40afb3000bcdfc7 4f6000
7557558ee 47efc7fb3bcdbb940a 91276b
7c9000000 47efc7fc7000000000 e794f6
8ee7c9fbb bb940a40a47efb3fb3 76b39d
f27755 fb347e00047ebb9 de4912
7c9fbb bcdbcd074fc7fc740a 39d
fbb8ee872 74fc747ebb947e074 a8fde4
872000000 47ebb9bb907440a40a 39da8f
8ee755755 fb340afb3074fc7bcd 4f639d
fbb7c9fbb bb947efb347efb3fb3 76b000
755f277c9 74bb9000bb9bcdfc7 4f6a8f
f277c97c9 47ebb9000000000 39d76b
7c98eefbb 40a40a a8f4f6
8727557c9 47efc7fc740abb9bb9 76bde4
f27fbbfbb 40abcd000fc747e a8f76b
8eef27 fc747e000bb9fb3 912
1fe1fe eed703 ed4
2bd3432ca 9ee703c8126cb8256f 2c5525
334000 c81c81000c81b82 ed4bf1
334000189 26cc81000703c81 bf1000
1892bd2bd eed26c7039ee56f56f c11000
3340771fe 9eeb82eedb82eed000 934525
3432bd2ca c81eed70326cb8256f 7e0c11
2bd2bd eed56feed9ee9eeeed ed47e0
1892bd343 c81b82b8256fb82b82 525bf1
3343432ca 9ee26c70326c70326c 2c5525
77334334 9ee000000eed000 2c5934
1fe1fe2bd 9ee000000b82000000 2c5c11
1fe2bd189 eedc8100056f56f bf1c11
2ca2bd2ca 56f9eeb829ee703c81 bf1934
2bd343 c81eedeed56f9ee9ee 525
1fe000000 9ee703c8126cb8256f ed4c11
343077334 56feed00026c703 2c5000
2bd2bd2ca c81703eed26cb8256f c11525
343000000 c81703703000000000 9342c5
2ca343077 b8256f56fc81eedeed 5257e0
1892bd eedc81000c81b82 ed4c11
343077 26c26c9ee70370356f 7e0
772ca334 9ee703c81b82c819ee bf1ed4
334000000 c81b82b829ee56f56f 7e0bf1
2ca2bd2bd eed56feed9ee70326c 2c57e0
77343077 b82c81eedc81eedeed 525000
2bd189343 9eeb82000b8226c703 2c5bf1
189343343 c81b82000000000 7e0525
3432ca077 56f56f bf12c5
3342bd343 c8170370356fb82b82 525ed4
189077077 56f26c000703c81 bf1525
2ca189 703c81000b82eed c11
f8ef8e cdfd61 95d
79b815772 1bed61147d98c260f9 3fd4f9
8fc000 147147000147c26 95dda4
8fc000f67 d98147000d61147 da4000
f6779b79b cdfd98d611be0f90f9 aa0000
8fc0e9f8e 1bec26cdfc26cdf000 e594f9
81579b772 147cdfd61d98c260f9 704aa0
79b79b cdf0f9cdf1be1becdf 95d704
f6779b815 147c26c260f9c26c26 4f9da4
8fc815772 1bed98d61d98d61d98 3fd4f9
e98fc8fc 1be000000cdf000 3fde59
f8ef8e79b 1be000000c26000000 3fdaa0
f8e79bf67 cdf1470000f90f9 da4aa0
77279b772 f91bec261bed61147 da4e59
79b815 147cdfcdf0f91be1be 4f9
f8e000000 1bed61147d98c260f9 95daa0
8150e98fc f9cdf000d98d61 3fd000
79b79b772 147d61cdfd98c260f9 aa04f9
815000000 147d61d61000000000 e593fd
7728150e9 c260f90f9147cdfcdf 4f9704
f6779b cdf147000147c26 95daa0
8150e9 d98d981bed61d610f9 704
e97728fc 1bed61147c261471be da495d
8fc000000 147c26c261be0f90f9 704da4
77279b79b cdf0f9cdf1bed61d98 3fd704
e98150e9 c26147cdf147cdfcdf 4f9000
79bf67815 1bec26000c26d98d61 3fdda4
f67815815 147c26000000000 7044f9
8157720e9 f90f9 da43fd
8fc79b815 147d61d610f9c26c26 4f995d
f670e90e9 f9d98000d61147 da44f9
772f67 d61147000c26cdf aa0
a00a0 198fad 114
3d8378b50 e35fad45a5c2bf7a6f fe0af
bf0000 45a45a00045abf7 1141bb
bf0000828 5c245a000fad45a 1bb000
8283d83d8 1985c2fade35a6fa6f 1ea000
bf08880a0 e35bf7198bf7198000 1450af
3783d8b50 45a198fad5c2bf7a6f 511ea
3d83d8 198a6f198e35e35198 114051
8283d8378 45abf7bf7a6fbf7bf7 af1bb
bf0378b50 e355c2fad5c2fad5c2 fe0af
888bf0bf0 e35000000198000 fe145
a00a03d8 e35000000bf7000000 fe1ea
a03d8828 19845a000a6fa6f 1bb1ea
b503d8b50 a6fe35bf7e35fad45a 1bb145
3d8378 45a198198a6fe35e35 af
a0000000 e35fad45a5c2bf7a6f 1141ea
378888bf0 a6f1980005c2fad fe000
3d83d8b50 45afad1985c2bf7a6f 1ea0af
378000000 45afadfad000000000 1450fe
b50378888 bf7a6fa6f45a198198 af051
8283d8 19845a00045abf7 1141ea
378888 5c25c2e35fadfada6f 51
888b50bf0 e35fad45abf745ae35 1bb114
bf0000000 45abf7bf7e35a6fa6f 511bb
b503d83d8 198a6f198e35fad5c2 fe051
888378888 bf745a19845a198198 af000
3d8828378 e35bf7000bf75c2fad fe1bb
828378378 45abf7000000000 510af
378b50888 a6fa6f 1bb0fe
bf03d8378 45afadfada6fbf7bf7 af114
828888888 a6f5c2000fad45a 1bb0af
b50828 fad45a000bf7198 1ea
''',
    4: '''
48 96 24 308
7e07e0 8b5c6b 87e
6a1110 3318b5000000f5ac6b 8a000
141 f5a7ef000000c6b7ef 88b000
1417b1 7ef7ef0000008b58b5 8a88b
6f06a1 331331 88b08a
7b16a1 331f5a 7f
110110 4de0000000008b5000 8a801
141051 8b57ef0000007efb84 87e8f4
6a16a1 331f5a000000c6b7ef 88b0f5
6a17e0 7eff5a000000f5a4de 8010f5
7b16f0 331b840000007efc6b 8a0f5
7b17b1 7efb84000000000000 7f0f5
51051 c6b4de0000007efb84 87e07f
7e0051 c6b7ef0000007ef4de 8a8f4
517b1 331331000000f5af5a f588b
6f0141 c6b7ef000000b848b5 8f4
516f0 4de4de000000000000 8f40f5
517b1 4de4de000000f5af5a f587e
7b16a1 4deb84000000000000 f5000
110000 7ef7ef0000007efb84 87e88b
51141 8b57ef000000331331 88b8f4
51051 c6b4de000000f5ac6b 8a07f
7e07e0000000 8b5c6b000000 87e
6a1110000000 3318b5000000f5ac6b000000 8a000
141000000 f5a7ef000000c6b7ef000000 88b000
1417b1000000 7ef7ef0000008b58b5000000 8a88b
6f06a1000000 331331000000 88b08a
7b16a1000000 331f5a000000 7f
110110000000 4de0000000008b5000000000 8a801
141051000000 8b57ef0000007efb84000000 87e8f4
6a16a1000000 331f5a000000c6b7ef000000 88b0f5
6a17e0000000 7eff5a000000f5a4de000000 8010f5
7b16f0000000 331b840000007efc6b000000 8a0f5
7b17b1000000 7efb84000000000000000000 7f0f5
51051000000 c6b4de0000007efb84000000 87e07f
7e0051000000 c6b7ef0000007ef4de000000 8a8f4
517b1000000 331331000000f5af5a000000 f588b
6f0141000000 c6b7ef000000b848b5000000 8f4
516f0000000 4de4de000000000000000000 8f40f5
517b1000000 4de4de000000f5af5a000000 f587e
7b16a1000000 4deb84000000000000000000 f5000
110000000000 7ef7ef0000007efb84000000 87e88b
51141000000 8b57ef000000331331000000 88b8f4
51051000000 c6b4de000000f5ac6b000000 8a07f
260260 86c426 f58
d6f250 eef86c000000ac9426 a6000
f0f ac92a50000004262a5 f9b000
f0ff3f 2a52a500000086c86c a6f9b
30d6f eefeef f9b0a6
f3fd6f eefac9 65
250250 c4a00000000086c000 a6f3d
f0fd5f 86c2a50000002a5683 f58ffe
d6fd6f eefac90000004262a5 f9b0c3
d6f260 2a5ac9000000ac9c4a f3d0c3
f3f030 eef6830000002a5426 a60c3
f3ff3f 2a5683000000000000 650c3
d5fd5f 426c4a0000002a5683 f58065
260d5f 4262a50000002a5c4a a6ffe
d5ff3f eefeef000000ac9ac9 c3f9b
30f0f 4262a500000068386c ffe
d5f030 c4ac4a000000000000 ffe0c3
d5ff3f c4ac4a000000ac9ac9 c3f58
f3fd6f c4a683000000000000 c3000
250000 2a52a50000002a5683 f58f9b
d5ff0f 86c2a5000000eefeef f9bffe
d5fd5f 426c4a000000ac9426 a6065
260260000000 86c426000000 f58
d6f250000000 eef86c000000ac9426000000 a6000
f0f000000 ac92a50000004262a5000000 f9b000
f0ff3f000000 2a52a500000086c86c000000 a6f9b
30d6f000000 eefeef000000 f9b0a6
f3fd6f000000 eefac9000000 65
250250000000 c4a00000000086c000000000 a6f3d
f0fd5f000000 86c2a50000002a5683000000 f58ffe
d6fd6f000000 eefac90000004262a5000000 f9b0c3
d6f260000000 2a5ac9000000ac9c4a000000 f3d0c3
f3f030000000 eef6830000002a5426000000 a60c3
f3ff3f000000 2a5683000000000000000000 650c3
d5fd5f000000 426c4a0000002a5683000000 f58065
260d5f000000 4262a50000002a5c4a000000 a6ffe
d5ff3f000000 eefeef000000ac9ac9000000 c3f9b
30f0f000000 4262a500000068386c000000 ffe
d5f030000000 c4ac4a000000000000000000 ffe0c3
d5ff3f000000 c4ac4a000000ac9ac9000000 c3f58
f3fd6f000000 c4a683000000000000000000 c3000
250000000000 2a52a50000002a5683000000 f58f9b
d5ff0f000000 86c2a5000000eefeef000000 f9bffe
d5fd5f000000 426c4a000000ac9426000000 a6065
664664 48a5e8 b87
c5d8ea fcc48a000000a245e8 561000
a39 a24eae0000005e8eae c8b000
a394b7 eaeeae00000048a48a 561c8b
e8ec5d fccfcc c8b561
4b7c5d fcca24 26d
8ea8ea 16200000000048a000 5619ea
a392d3 48aeae000000eaeb46 b87ee6
c5dc5d fcca240000005e8eae c8b70c
c5d664 eaea24000000a24162 9ea70c
4b7e8e fccb46000000eae5e8 56170c
4b74b7 eaeb46000000000000 26d70c
2d32d3 5e8162000000eaeb46 b8726d
6642d3 5e8eae000000eae162 561ee6
2d34b7 fccfcc000000a24a24 70cc8b
e8ea39 5e8eae000000b4648a ee6
2d3e8e 162162000000000000 ee670c
2d34b7 162162000000a24a24 70cb87
4b7c5d 162b46000000000000 70c000
8ea000 eaeeae000000eaeb46 b87c8b
2d3a39 48aeae000000fccfcc c8bee6
2d32d3 5e8162000000a245e8 56126d
664664000000 48a5e8000000 b87
c5d8ea000000 fcc48a000000a245e8000000 561000
a39000000 a24eae0000005e8eae000000 c8b000
a394b7000000 eaeeae00000048a48a000000 561c8b
e8ec5d000000 fccfcc000000 c8b561
4b7c5d000000 fcca24000000 26d
8ea8ea000000 16200000000048a000000000 5619ea
a392d3000000 48aeae000000eaeb46000000 b87ee6
c5dc5d000000 fcca240000005e8eae000000 c8b70c
c5d664000000 eaea24000000a24162000000 9ea70c
4b7e8e000000 fccb46000000eae5e8000000 56170c
4b74b7000000 eaeb46000000000000000000 26d70c
2d32d3000000 5e8162000000eaeb46000000 b8726d
6642d3000000 5e8eae000000eae162000000 561ee6
2d34b7000000 fccfcc000000a24a24000000 70cc8b
e8ea39000000 5e8eae000000b4648a000000 ee6
2d3e8e000000 162162000000000000000000 ee670c
2d34b7000000 162162000000a24a24000000 70cb87
4b7c5d000000 162b46000000000000000000 70c000
8ea000000000 eaeeae000000eaeb46000000 b87c8b
2d3a39000000 48aeae000000fccfcc000000 c8bee6
2d32d3000000 5e8162000000a245e8000000 56126d
9c09c fb3fc7 de4
fbb872 40afb3000000bcdfc7 4f6000
f27 bcd47e000000fc747e a8f000
f277c9 47e47e000000fb3fb3 4f6a8f
8eefbb 40a40a a8f4f6
7c9fbb 40abcd 39d
872872 74000000000fb3000 4f6e79
f27755 fb347e00000047ebb9 de4912
fbbfbb 40abcd000000fc747e a8f76b
fbb09c 47ebcd000000bcd074 e7976b
7c98ee 40abb900000047efc7 4f676b
7c97c9 47ebb9000000000000 39d76b
755755 fc707400000047ebb9 de439d
9c755 fc747e00000047e074 4f6912
7557c9 40a40a000000bcdbcd 76ba8f
8eef27 fc747e000000bb9fb3 912
7558ee 74074000000000000 91276b
7557c9 74074000000bcdbcd 76bde4
7c9fbb 74bb9000000000000 76b000
872000 47e47e00000047ebb9 de4a8f
755f27 fb347e00000040a40a a8f912
755755 fc7074000000bcdfc7 4f639d
9c09c000000 fb3fc7000000 de4
fbb872000000 40afb3000000bcdfc7000000 4f6000
f27000000 bcd47e000000fc747e000000 a8f000
f277c9000000 47e47e000000fb3fb3000000 4f6a8f
8eefbb000000 40a40a000000 a8f4f6
7c9fbb000000 40abcd000000 39d
872872000000 74000000000fb3000000000 4f6e79
f27755000000 fb347e00000047ebb9000000 de4912
fbbfbb000000 40abcd000000fc747e000000 a8f76b
fbb09c000000 47ebcd000000bcd074000000 e7976b
7c98ee000000 40abb900000047efc7000000 4f676b
7c97c9000000 47ebb9000000000000000000 39d76b
755755000000 fc707400000047ebb9000000 de439d
9c755000000 fc747e00000047e074000000 4f6912
7557c9000000 40a40a000000bcdbcd000000 76ba8f
8eef27000000 fc747e000000bb9fb3000000 912
7558ee000000 74074000000000000000000 91276b
7557c9000000 74074000000bcdbcd000000 76bde4
7c9fbb000000 74bb9000000000000000000 76b000
872000000000 47e47e00000047ebb9000000 de4a8f
755f27000000 fb347e00000040a40a000000 a8f912
755755000000 fc7074000000bcdfc7000000 4f639d
1fe1fe eed703 ed4
77334 56feed00000026c703 2c5000
189 26cc81000000703c81 bf1000
189343 c81c81000000eedeed 2c5bf1
2ca077 56f56f bf12c5
343077 56f26c 7e0
334334 9ee000000000eed000 2c5934
1892bd eedc81000000c81b82 ed4c11
77077 56f26c000000703c81 bf1525
771fe c8126c00000026c9ee 934525
3432ca 56fb82000000c81703 2c5525
343343 c81b82000000000000 7e0525
2bd2bd 7039ee000000c81b82 ed47e0
1fe2bd 703c81000000c819ee 2c5c11
2bd343 56f56f00000026c26c 525bf1
2ca189 703c81000000b82eed c11
2bd2ca 9ee9ee000000000000 c11525
2bd343 9ee9ee00000026c26c 525ed4
343077 9eeb82000000000000 525000
334000 c81c81000000c81b82 ed4bf1
2bd189 eedc8100000056f56f bf1c11
2bd2bd 7039ee00000026c703 2c57e0
1fe1fe000000 eed703000000 ed4
77334000000 56feed00000026c703000000 2c5000
189000000 26cc81000000703c81000000 bf1000
189343000000 c81c81000000eedeed000000 2c5bf1
2ca077000000 56f56f000000 bf12c5
343077000000 56f26c000000 7e0
334334000000 9ee000000000eed000000000 2c5934
1892bd000000 eedc81000000c81b82000000 ed4c11
77077000000 56f26c000000703c81000000 bf1525
771fe000000 c8126c00000026c9ee000000 934525
3432ca000000 56fb82000000c81703000000 2c5525
343343000000 c81b82000000000000000000 7e0525
2bd2bd000000 7039ee000000c81b82000000 ed47e0
1fe2bd000000 703c81000000c819ee000000 2c5c11
2bd343000000 56f56f00000026c26c000000 525bf1
2ca189000000 703c81000000b82eed000000 c11
2bd2ca000000 9ee9ee000000000000000000 c11525
2bd343000000 9ee9ee00000026c26c000000 525ed4
343077000000 9eeb82000000000000000000 525000
334000000000 c81c81000000c81b82000000 ed4bf1
2bd189000000 eedc8100000056f56f000000 bf1c11
2bd2bd000000 7039ee00000026c703000000 2c57e0
f8ef8e cdfd61 95d
e98fc f9cdf000000d98d61 3fd000
f67 d98147000000d61147 da4000
f67815 147147000000cdfcdf 3fdda4
7720e9 f90f9 da43fd
8150e9 f9d98 704
8fc8fc 1be000000000cdf000 3fde59
f6779b cdf147000000147c26 95daa0
e90e9 f9d98000000d61147 da44f9
e9f8e 147d98000000d981be e594f9
815772 f9c26000000147d61 3fd4f9
815815 147c26000000000000 7044f9
79b79b d611be000000147c26 95d704
f8e79b d611470000001471be 3fdaa0
79b815 f90f9000000d98d98 4f9da4
772f67 d61147000000c26cdf aa0
79b772 1be1be000000000000 aa04f9
79b815 1be1be000000d98d98 4f995d
8150e9 1bec26000000000000 4f9000
8fc000 147147000000147c26 95dda4
79bf67 cdf1470000000f90f9 da4aa0
79b79b d611be000000d98d61 3fd704
f8ef8e000000 cdfd61000000 95d
e98fc000000 f9cdf000000d98d61000000 3fd000
f67000000 d98147000000d61147000000 da4000
f67815000000 147147000000cdfcdf000000 3fdda4
7720e9000000 f90f9000000 da43fd
8150e9000000 f9d98000000 704
8fc8fc000000 1be000000000cdf000000000 3fde59
f6779b000000 cdf147000000147c26000000 95daa0
e90e9000000 f9d98000000d61147000000 da44f9
e9f8e000000 147d98000000d981be000000 e594f9
815772000000 f9c26000000147d61000000 3fd4f9
815815000000 147c26000000000000000000 7044f9
79b79b000000 d611be000000147c26000000 95d704
f8e79b000000 d611470000001471be000000 3fdaa0
79b815000000 f90f9000000d98d98000000 4f9da4
772f67000000 d61147000000c26cdf000000 aa0
79b772000000 1be1be000000000000000000 aa04f9
79b815000000 1be1be000000d98d98000000 4f995d
8150e9000000 1bec26000000000000000000 4f9000
8fc000000000 147147000000147c26000000 95dda4
79bf67000000 cdf1470000000f90f9000000 da4aa0
79b79b000000 d611be000000d98d61000000 3fd704
a00a0 198fad 114
888bf0 a6f1980000005c2fad fe000
828 5c245a000000fad45a 1bb000
828378 45a45a000000198198 fe1bb
b50888 a6fa6f 1bb0fe
378888 a6f5c2 51
bf0bf0 e35000000000198000 fe145
8283d8 19845a00000045abf7 1141ea
888888 a6f5c2000000fad45a 1bb0af
8880a0 45a5c20000005c2e35 1450af
378b50 a6fbf700000045afad fe0af
378378 45abf7000000000000 510af
3d83d8 fade3500000045abf7 114051
a03d8 fad45a00000045ae35 fe1ea
3d8378 a6fa6f0000005c25c2 af1bb
b50828 fad45a000000bf7198 1ea
3d8b50 e35e35000000000000 1ea0af
3d8378 e35e350000005c25c2 af114
378888 e35bf7000000000000 af000
bf0000 45a45a00000045abf7 1141bb
3d8828 19845a000000a6fa6f 1bb1ea
3d83d8 fade350000005c2fad fe051
a00a0000000 198fad000000 114
888bf0000000 a6f1980000005c2fad000000 fe000
828000000 5c245a000000fad45a000000 1bb000
828378000000 45a45a000000198198000000 fe1bb
b50888000000 a6fa6f000000 1bb0fe
378888000000 a6f5c2000000 51
bf0bf0000000 e35000000000198000000000 fe145
8283d8000000 19845a00000045abf7000000 1141ea
888888000000 a6f5c2000000fad45a000000 1bb0af
8880a0000000 45a5c20000005c2e35000000 1450af
378b50000000 a6fbf700000045afad000000 fe0af
378378000000 45abf7000000000000000000 510af
3d83d8000000 fade3500000045abf7000000 114051
a03d8000000 fad45a00000045ae35000000 fe1ea
3d8378000000 a6fa6f0000005c25c2000000 af1bb
b50828000000 fad45a000000bf7198000000 1ea
3d8b50000000 e35e35000000000000000000 1ea0af
3d8378000000 e35e350000005c25c2000000 af114
378888000000 e35bf7000000000000000000 af000
bf0000000000 45a45a00000045abf7000000 1141bb
3d8828000000 19845a000000a6fa6f000000 1bb1ea
3d83d8000000 fade350000005c2fad000000 fe051
''',
}


if __name__ == '__main__':
    _main()
