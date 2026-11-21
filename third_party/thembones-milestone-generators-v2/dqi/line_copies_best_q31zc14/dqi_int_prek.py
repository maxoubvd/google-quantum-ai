#!/usr/bin/env python3
"""dqi_int_prek.py -- B5 (24th holder, Nov 18): the converter's runops products as ONE two-by-n kernel, in the
DQI_INT_EXT form (install(km, ns); token dqi_int_prek, LAST).

What it replaces. The kit's runops (playbook/kit/decoder/exp/safe_scalar/converter.py, cloned into the converter
scopes) executes an op list in order (reversed when undo):
    ('s', x, flag) scale x in place; ('m', x, y, z) ar.mul_xor0(b, x, y, z) = z ^= x y; ('c', ...) a convolution.
The kit's pre (same file, and its clone in exp/r-main/narrow240/measured_heads.py) puts into ONE list
    [('m', g[i], v[-1], v[i]) for v in payloads for i in range(n)]
= K(2, n): the two payload scalars v[-1] times the n words g[i], each product XORed into its own word v[i]: 2 n
K-products of 42 ANDs, consecutive. R81-DQI-FOLD's K(2,2) at 22 ANDs over GF(8) (lifted to K: 154 terms) and K(2,3) at
32 (224 terms) compute exactly this map (files tools/r81_dqi_fold_k22_K154.txt, ..._k23_K224.txt; checked by the
designer and by the re-pricer's independent code); n = 4 is two K(2,2) blocks (308). Saving 14 / 28 / 28 a block.

How. The swapped runops passes its op sequence through plan(): every MAXIMAL run of consecutive 'm' ops is buffered and
flushed before the next non-'m' op runs; a run that is exactly the full product {s_0, s_1} x {w_0 .. w_(n-1)} (two
distinct scalar registers on the y side, n distinct words on the x side, one op for each pair, 2 n distinct targets,
all 2 + 3 n registers of 12 wires on pairwise distinct wires) is emitted as the fused kernel at that moment; any other
run is yielded op by op in its original order (the kit's own code then runs it). Nothing else changes: the 's' and
'c' ops, their order, and every 'm' op that is not in such a run.

Exactness (the lane's rule: exact by form or checked on all inputs).
 1. At install each scheme is checked as a TENSOR IDENTITY against the definition (output word i n + j ^= s_i w_j in
    GF(2^12) mod 0x1053, this file's own multiplier, cross-checked against km.GF2Field(12) on the 144 basis pairs) on all
    basis pairs = all inputs (the map is bilinear). A scheme that fails stops the build (assert).
 2. A runops is swapped only where its source text equals the kit's runops (whitespace-normalised); others refused.
 3. The emitted gate list is the in-place rank-one form of field.mul_xor0 (CNOT fan-in of a term's first mask onto one
    wire of the scalars, of its second mask onto one wire of the words, CNOT fan-out of its output mask from one target
    wire, one CCX, undo): 0 work wires. `python -I tools/dqi_int_prek.py test` simulates the emitted list SYMBOLICALLY
    (every wire a GF(2) form of the inputs plus products of two forms): operands restored, targets ^= the definition,
    for n = 2, 3, 4 and both run orders (forward and the undo order), plus negative controls and pass-through cases.
Environment: DQI_INT_PREK_OFF=1 installs nothing; DQI_INT_PREK_N=2,3,4 (default) limits the sizes;
DQI_INT_PREK_OUT=<json> receives the call log after every fused block.
Credits: R81-DQI-FOLD (the K(2,2) / K(2,3) schemes and their lift; public mathematics: flip graphs, M. Kauers and
J. Moosbauer, arXiv 2212.01175; plus transitions, Y. Arai, Y. Ichikawa, K. Hukushima, arXiv 2312.17094; Toom-Cook); the
emitter follows tools/r81_dqi_fold_dots.py's emit() (R81 round 2) and the kit's field.mul_xor0.
"""
import gc
import inspect
import json
import os
import sys
import types

MOD12 = 0x1053
NAME = 'dqi_int_prek'
OFF = os.environ.get('DQI_INT_PREK_OFF', '0') == '1'
SIZES = tuple(int(x) for x in os.environ.get('DQI_INT_PREK_N', '2,3,4').split(',') if x)
OUT = os.environ.get('DQI_INT_PREK_OUT')
HERE = os.path.dirname(os.path.abspath(__file__))
FILES = {2: os.environ.get('DQI_INT_PREK_FILE2', 'r81_dqi_fold_k22_K154.txt'),       # any K-level K(2,2) / K(2,3)
         3: os.environ.get('DQI_INT_PREK_FILE3', 'r81_dqi_fold_k23_K224.txt')}       # scheme in R81's format
RANK = {2: 154, 3: 224}                    # replaced by the files' headers in _load() (each checked exactly there)
KIT = {2: 168, 3: 252, 4: 336}
LOG = dict(installed=False, off=OFF, sizes=list(SIZES), swapped=[], refused=[], functions=0, fused={}, runs={},
           sites={}, saved=0)
TERMS = {}

REF_SRC = """def runops(b,ops,a,ai,undo=False):
    for op in (reversed(ops)if undo else ops):
        if op[0]=='s':scale(b,op[1],a,ai,op[2]^undo)
        elif op[0]=='m':ar.mul_xor0(b,op[1],op[2],op[3])
        elif op[0]=='c':convolution_xor(b,op[1],op[2],op[3],undo)
        else:raise ValueError(op[0])
"""
NEW_SRC = """def runops(b,ops,a,ai,undo=False):
    for op in __import__('dqi_int_prek').plan(b,(reversed(ops)if undo else ops)):
        if op[0]=='s':scale(b,op[1],a,ai,op[2]^undo)
        elif op[0]=='m':ar.mul_xor0(b,op[1],op[2],op[3])
        elif op[0]=='c':convolution_xor(b,op[1],op[2],op[3],undo)
        else:raise ValueError(op[0])
"""


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
    """Definition on a basis pair: operand a bit i = x^la of scalar sa; operand b bit j = x^lb of word wb; the product
    s_sa w_wb lands in output word sa n + wb (head-major, as tools/r81_dqi_fold_klift.py lays out K(p, n))."""
    sa, la = divmod(i, 12)
    wb, lb = divmod(j, 12)
    return gmul(1 << la, 1 << lb) << (12 * (sa * n + wb))


def check_terms(n, terms):
    na, nb, nc = 24, 12 * n, 24 * n
    T = [[0] * nb for _ in range(na)]
    for a, bm, c in terms:
        if not (a and bm and c) or a >> na or bm >> nb or c >> nc:
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


def _load():
    for n, fn in FILES.items():
        if n in TERMS:
            continue
        with open(fn if os.path.isabs(fn) else os.path.join(HERE, fn)) as fh:
            lines = [ln.split() for ln in fh if ln.strip()]
        na, nb, nc, r = (int(v) for v in lines[0])
        assert (na, nb, nc, r) == (24, 12 * n, 24 * n, len(lines) - 1), 'header %s' % fn
        terms = [tuple(int(v, 16) for v in ln[:3]) for ln in lines[1:]]
        ok, bad = check_terms(n, terms)
        assert ok, '%s is NOT the K(2,%d) map: %s' % (fn, n, bad)
        TERMS[n] = terms
        RANK[n] = r
    KCOST.update({2: RANK[2], 3: RANK[3], 4: 2 * RANK[2]})
    LOG['ranks'] = dict(RANK)
    return TERMS


def emit(b, A, B, C, terms):
    """One CCX a term in the kit's in-place rank-one form (field.mul_xor0), two masks, flat wire lists."""
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


def _note(key, n):
    d = LOG.setdefault(key, {})
    d[str(n)] = d.get(str(n), 0) + 1
    if OUT:
        try:
            with open(OUT, 'w') as fh:
                json.dump(LOG, fh)
        except OSError:
            pass


def _wires(reg):
    return tuple(str(q) for q in reg)


def match(run):
    """run: consecutive ('m', x, y, z) ops. Returns (scalars, words, targets[s][w]) when the run is exactly the full
    product of two scalar registers (y side) with n word registers (x side), else None."""
    if len(run) < 4 or len(run) % 2:
        return None
    try:
        ys, xs, tg = [], [], {}
        for op in run:
            x, y, z = op[1], op[2], op[3]
            if not all(len(r) == 12 for r in (x, y, z)):
                return None
            ky, kx = _wires(y), _wires(x)
            if ky not in [_wires(s) for s in ys]:
                ys.append(y)
            if kx not in [_wires(w) for w in xs]:
                xs.append(x)
            key = (ky, kx)
            if key in tg:
                return None
            tg[key] = z
    except (TypeError, IndexError):
        return None
    n = len(xs)
    if len(ys) != 2 or n < 2 or len(run) != 2 * n or len(tg) != 2 * n:
        return None
    regs = ys + xs + list(tg.values())
    names = [w for r in regs for w in _wires(r)]
    if len(set(names)) != len(names):
        return None
    targets = [[tg[(_wires(s), _wires(w))] for w in xs] for s in ys]
    return ys, xs, targets


def fuse(b, ys, xs, targets):
    """Emit K(2, n) for the matched run; n = 4 as two K(2,2) blocks. Returns the Toffolis saved."""
    n = len(xs)
    blocks = [list(range(n))] if n in TERMS else ([[0, 1], [2, 3]] if n == 4 else None)
    if blocks is None:
        return None
    saved = 0
    for blk in blocks:
        k = len(blk)
        A = [q for s in ys for q in s]
        B = [q for j in blk for q in xs[j]]
        C = [q for si in range(2) for j in blk for q in targets[si][j]]
        emit(b, A, B, C, TERMS[k])
        saved += 84 * k - RANK[k]
    return saved


def plan(b, seq):
    """Yield the ops of seq in order; each maximal run of 'm' ops that is a K(2, n) block is emitted here instead."""
    run = []

    def flush():
        if not run:
            return []
        m = None if OFF else match(run)
        if m is not None and len(m[1]) in SIZES:
            sv = fuse(b, *m)
            if sv is not None:
                n = len(m[1])
                LOG['saved'] += sv
                co = sys._getframe(3).f_code if hasattr(sys, '_getframe') else None
                site = '%s:%s n=%d' % (os.path.basename(co.co_filename) if co else '?', co.co_name if co else '?', n)
                LOG['sites'][site] = LOG['sites'].get(site, 0) + 1
                _note('fused', n)
                del run[:]
                return []
        _note('runs', len(run))
        out = list(run)
        del run[:]
        return out

    for op in seq:
        if op[0] == 'm':
            run.append(op)
            continue
        for o in flush():
            yield o
        yield op
    for o in flush():
        yield o


def _norm(src):
    return ''.join(src.split())


def _field_check(km):
    """The kit's multiplier must equal this file's on the 144 basis pairs: through GF2Field of kickmix (km when a caller
    passes it; the line's entry passes km = None, tools/dqi_sw_entry.py L30: then sys.modules circ / kickmix), else the
    kit's own field module (field(12).mul, the module that carries mul_xor0). Returns which one was used."""
    mul, used = None, None
    for cand in (km, sys.modules.get('circ'), sys.modules.get('kickmix')):
        try:
            if cand is not None and hasattr(cand, 'GF2Field'):
                f = cand.GF2Field(12)
                assert f.modulus == MOD12, 'field modulus %#x' % f.modulus
                mul, used = f.mul, 'GF2Field'
                break
        except AssertionError:
            raise
        except Exception:
            continue
    if mul is None:
        for name, m in list(sys.modules.items()):
            try:
                if m is not None and callable(getattr(m, 'mul_xor0', None)) and callable(getattr(m, 'field', None)):
                    mul, used = m.field(12).mul, 'field module ' + name
                    break
            except Exception:
                continue
    assert mul is not None, 'no field multiplier found to check against'
    for i in range(12):
        for j in range(12):
            assert mul(1 << i, 1 << j) == gmul(1 << i, 1 << j), 'multiplier differs from the kit field (%s)' % used
    return used


KCOST = {2: 154, 3: 224, 4: 308}


def _regs_ok(regs):
    try:
        if not all(len(r) == 12 for r in regs):
            return False
        names = [w for r in regs for w in _wires(r)]
    except (TypeError, IndexError):
        return False
    return len(set(names)) == len(names)


def _make_pair_level(orig):
    """D1's pair_level(b, prog, words, d1, d2, T1, T2): T1[i] ^= words[i] d1, T2[i] ^= words[i] d2 (tools/dqi_pub_pair.py
    L848; callers: its 'P' ops of the converter's runops, its plain loop, B5's plainpair / pairscale). The same map is
    K(2, n) with scalars (d1, d2): for n = 2, 3 it is emitted by R81's scheme when that is cheaper than the loaded
    program's own count (the line's mixed program: 67 n + 34 = 168 / 235 against 154 / 224); else D1's program."""
    def pair_level(b, prog, words, d1, d2, T1, T2):
        n = len(words)
        if not OFF and n in SIZES and n in (2, 3) and n in TERMS:
            try:
                theirs = int(prog.counts(n)['total'])
            except Exception:
                theirs = None
            if theirs is not None and KCOST[n] < theirs and len(T1) == len(T2) == n \
                    and _regs_ok([d1, d2] + list(words) + list(T1) + list(T2)):
                fuse(b, [d1, d2], list(words), [list(T1), list(T2)])
                LOG['saved'] += theirs - KCOST[n]
                LOG['pair_saved'] = LOG.get('pair_saved', 0) + theirs - KCOST[n]
                _note('pair', n)
                return None
        _note('pair_kept', n)
        return orig(b, prog, words, d1, d2, T1, T2)
    pair_level._b5_prek = True
    pair_level._b5_orig = orig
    return pair_level


def _make_pair_ops(orig):
    """D1's _pair_ops (the runops pre-pass that turns K(2, n >= 3) runs into ('P', W, h1, h2, T1, T2)), then the same
    for a remaining run of exactly four 'm' ops that is K(2, 2): it becomes one 'P' op, run by pair_level (above)."""
    def pair_ops(ops):
        out = orig(ops)
        res, run = [], []
        for op in list(out) + [None]:
            if op is not None and op[0] == 'm':
                run.append(op)
                continue
            if run:
                m = match(run) if (len(run) == 4 and 2 in SIZES and not OFF) else None
                if m is not None and len(m[1]) == 2:
                    ys, xs, tg = m
                    res.append(('P', list(xs), ys[0], ys[1], list(tg[0]), list(tg[1])))
                    _note('grouped', 2)
                else:
                    res.extend(run)
                run = []
            if op is not None:
                res.append(op)
        return res
    pair_ops._b5_prek = True
    return pair_ops


def _install_pub_pair():
    """Wraps dqi_pub_pair.pair_level (a module global: its _pair_run / _pair_loop and B5's plainpair look it up at call
    time) and, in every function's globals that hold D1's _D1_PAIR_OPS, the pre-pass. Returns a log dict."""
    pp = sys.modules.get('dqi_pub_pair')
    if pp is None or not hasattr(pp, 'pair_level'):
        return dict(pub_pair='absent')
    out = dict(pub_pair='present')
    if not getattr(pp.pair_level, '_b5_prek', False):
        pp.pair_level = _make_pair_level(pp.pair_level)
        out['pair_level'] = 'wrapped'
    else:
        out['pair_level'] = 'already'
    orig_ops = getattr(pp, '_pair_ops', None)
    wrapped = {}
    for fn in [o for o in gc.get_objects() if isinstance(o, types.FunctionType)]:
        gl = fn.__globals__
        cur = gl.get('_D1_PAIR_OPS')
        if cur is None or id(gl) in wrapped:
            continue
        if getattr(cur, '_b5_prek', False):
            wrapped[id(gl)] = 'already'
        elif orig_ops is not None and cur is orig_ops:
            gl['_D1_PAIR_OPS'] = _make_pair_ops(cur)
            wrapped[id(gl)] = 'wrapped'
        else:
            wrapped[id(gl)] = 'foreign'
    out['pair_ops_globals'] = sorted(wrapped.values())
    return out


def install(km, ns):
    if OFF:
        return dict(LOG)
    LOG['field'] = _field_check(km)
    _load()
    LOG['pub_pair'] = _install_pub_pair()
    ref = _norm(REF_SRC)
    newco = [c for c in compile(NEW_SRC, '<' + NAME + '>', 'exec').co_consts if isinstance(c, types.CodeType)][0]
    done = {}
    for fn in [o for o in gc.get_objects() if isinstance(o, types.FunctionType)]:
        co = fn.__code__
        if co.co_name != 'runops':
            continue
        if id(co) not in done:
            where = '%s:%d' % (os.path.basename(co.co_filename), co.co_firstlineno)
            new = None
            try:
                src = inspect.getsource(co)
            except (OSError, TypeError):
                src = None
            if co is newco or (co.co_filename == '<' + NAME + '>'):
                pass
            elif src is not None and _norm(src) == ref:
                new = newco.replace(co_filename=co.co_filename, co_firstlineno=co.co_firstlineno)
                LOG['swapped'].append(where)
            else:
                LOG['refused'].append(where)
            done[id(co)] = (co, new)
        new = done[id(co)][1]
        if new is not None:
            fn.__code__ = new
            LOG['functions'] += 1
    LOG['installed'] = LOG['functions'] > 0
    sys.stderr.write(NAME + ': %s\n' % json.dumps(LOG))
    return dict(LOG)


# ---------------------------------------------------------------------------------------------------------------
# tests (no kit needed): python -I tools/dqi_int_prek.py test


class _Rec:
    def __init__(self):
        self.ops = []

    def cx(self, c, t):
        self.ops.append((c, t))

    def ccx(self, a, b, t):
        self.ops.append((a, b, t))


def _symbolic(ops, nwires, a_set, b_set, want):
    """Symbolic run: wire w starts as x_w; CX adds forms and product sets; CCX needs two linear controls (one in a_set's
    span, one in b_set's span) and adds their product (a set of (i, j) pairs) to the target. want[w] = the expected set
    of (i, j) pairs on wire w. Returns the list of defects."""
    lin = [1 << w for w in range(nwires)]
    quad = [set() for _ in range(nwires)]
    amask = sum(1 << w for w in a_set)
    bmask = sum(1 << w for w in b_set)
    for op in ops:
        if len(op) == 2:
            c, t = op
            lin[t] ^= lin[c]
            quad[t] ^= quad[c]
        else:
            x, y, t = op
            if quad[x] or quad[y]:
                return ['control carries a product']
            lx, ly = lin[x], lin[y]
            if lx & ~amask or ly & ~bmask or not lx or not ly:
                return ['control outside its operand']
            quad[t] ^= {(i, j) for i in _bits(lx) for j in _bits(ly)}
    bad = []
    for w in range(nwires):
        if lin[w] != 1 << w:
            bad.append(('linear', w))
        if quad[w] != want.get(w, set()):
            bad.append(('quad', w))
    return bad


def _kit_like(n, undo, rec, mutate=None):
    """A runops call shaped like the kit's pre: two scale ops, the 2 n 'm' ops, one more 'm' and a conv; returns the
    ops the fused path left to the kit, and the layout."""
    g = [[12 * j + l for l in range(12)] for j in range(n)]
    base = 12 * n
    pay = [[[base + 12 * (k * (n + 1) + j) + l for l in range(12)] for j in range(n + 1)] for k in range(2)]
    ops = [('s', v[-1], True) for v in pay]
    ops += [('m', g[i], v[-1], v[i]) for v in pay for i in range(n)]
    if mutate:
        ops = mutate(ops)
    ops += [('c', 'conv')]
    seq = reversed(ops) if undo else ops
    left = list(plan(rec, seq))
    return left, g, pay, base + 24 * (n + 1)


def _test():
    ok = True
    _load()
    for n in (2, 3):
        print('K(2,%d): %s: tensor identity on all 24 x %d basis pairs PASS (%d terms against %d)' % (
            n, FILES[n], 12 * n, RANK[n], KIT[n]))
    for n in (2, 3, 4):
        for undo in (False, True):
            rec = _Rec()
            before = dict(LOG['fused'])
            left, g, pay, N = _kit_like(n, undo, rec)
            kinds = [o[0] for o in left]
            # the products: product of scalar k (bits of pay[k][-1]) and word i (bits of g[i]) into pay[k][i]
            want = {}
            a_set = [q for k in range(2) for q in pay[k][-1]]
            b_set = [q for j in range(n) for q in g[j]]
            for k in range(2):
                for i in range(n):
                    for la in range(12):
                        for lb in range(12):
                            prod = gmul(1 << la, 1 << lb)
                            for l in _bits(prod):
                                w = pay[k][i][l]
                                want.setdefault(w, set()).symmetric_difference_update({(pay[k][-1][la], g[i][lb])})
            bad = _symbolic(rec.ops, N, a_set, b_set, want)
            ccx = sum(1 for o in rec.ops if len(o) == 3)
            exp = KCOST[n]
            good = not bad and ccx == exp and 'm' not in kinds and kinds.count('s') == 2 and kinds.count('c') == 1
            ok &= good
            print('n = %d, undo %s: fused %d CCX (kit %d), %d CX, 0 allocations; ops left to the kit %s; symbolic: %s' % (
                n, undo, ccx, KIT[n], len(rec.ops) - ccx, kinds, 'PASS' if good else 'FAIL %s' % bad[:4]))
    # pass-through: n = 1, a run with a foreign op, a repeated target, three scalars
    cases = {
        'n = 1 (a star)': (1, None),
        'one op dropped': (2, lambda ops: ops[:3] + ops[4:]),
        'a target repeated': (2, lambda ops: ops[:2] + [ops[2], (ops[3][0], ops[3][1], ops[3][2], ops[2][3])] + ops[4:]),
        'an op interleaved': (2, lambda ops: ops[:4] + [('s', ops[0][1], True)] + ops[4:]),
    }
    for label, (n, mut) in cases.items():
        rec = _Rec()
        left, g, pay, N = _kit_like(n, False, rec, mut)
        good = not rec.ops and sum(1 for o in left if o[0] == 'm') == sum(
            1 for o in (mut(list(_kit_like_ops(n))) if mut else _kit_like_ops(n)) if o[0] == 'm')
        ok &= good
        print('pass-through, %s: fused gates %d, ops left in order %s: %s' % (
            label, len(rec.ops), ''.join(o[0] for o in left), 'PASS' if good else 'FAIL'))
    # negative control: one flipped bit in a term is caught by the tensor check
    caught = 0
    for n in (2, 3):
        for pos, which in ((0, 0), (len(TERMS[n]) // 2, 1), (len(TERMS[n]) - 1, 2)):
            t = list(TERMS[n][pos])
            t[which] ^= 1 << 3
            caught += not check_terms(n, TERMS[n][:pos] + [tuple(t)] + TERMS[n][pos + 1:])[0]
    print('negative controls (one bit flipped in a term): caught %d of 6' % caught)
    ok &= caught == 6
    # the swap on a text copy of the kit's runops (no kit): source equality, patched function, refusal of others
    fake = types.SimpleNamespace(GF2Field=lambda n: types.SimpleNamespace(modulus=MOD12, mul=gmul))
    calls = []
    scope = dict(ar=types.SimpleNamespace(mul_xor0=lambda b, x, y, z: calls.append('m')),
                 scale=lambda b, x, a, ai, d: calls.append('s'), convolution_xor=lambda b, *a: calls.append('c'))
    path = os.path.join(os.environ.get('TMPDIR', '/tmp'), 'dqi_int_prek_kitcopy.py')
    with open(path, 'w') as fh:
        fh.write(REF_SRC + '\n\ndef other(b):\n    return b\n')
    exec(compile(open(path).read(), path, 'exec'), scope)
    other_scope = dict(scope)
    exec(compile(REF_SRC.replace("op[2]^undo", "op[2]"), path + '.x', 'exec'), other_scope)
    sys.modules.setdefault(NAME, sys.modules[__name__])
    log = install(fake, None)
    rec = _Rec()
    g = [[12 * j + l for l in range(12)] for j in range(2)]
    pay = [[[24 + 36 * k + 12 * j + l for l in range(12)] for j in range(3)] for k in range(2)]
    ops = [('s', v[-1], True) for v in pay] + [('m', g[i], v[-1], v[i]) for v in pay for i in range(2)]
    del calls[:]
    scope['runops'](rec, ops, None, None, False)
    fused_ok = calls == ['s', 's'] and sum(1 for o in rec.ops if len(o) == 3) == 154
    del calls[:]
    rec2 = _Rec()
    other_scope['runops'](rec2, ops, None, None, False)
    other_ok = calls == ['s', 's', 'm', 'm', 'm', 'm'] and not rec2.ops
    print('swap on a text copy of the kit runops: swapped %s, refused %s; patched call: kit ops %s + %d CCX: %s; '
          'a different runops left alone: %s' % (log['swapped'], log['refused'], ['s', 's'],
                                                 sum(1 for o in rec.ops if len(o) == 3), fused_ok, other_ok))
    ok &= fused_ok and other_ok and len(log['swapped']) == 1 and len(log['refused']) == 1
    ok &= _test_pub_pair()
    print('TEST %s' % ('PASS' if ok else 'FAIL'))
    return ok


def _pair_layout(n):
    d = [[12 * k + l for l in range(12)] for k in range(2)]
    w = [[24 + 12 * j + l for l in range(12)] for j in range(n)]
    T = [[[24 + 12 * n + 12 * (k * n + j) + l for l in range(12)] for j in range(n)] for k in range(2)]
    return d, w, T, 24 + 36 * n


def _pair_want(d, w, T, n):
    want = {}
    for k in range(2):
        for i in range(n):
            for la in range(12):
                for lb in range(12):
                    for l in _bits(gmul(1 << la, 1 << lb)):
                        want.setdefault(T[k][i][l], set()).symmetric_difference_update({(d[k][la], w[i][lb])})
    return want


def _test_pub_pair():
    """The D1 path without the kit: a stand-in dqi_pub_pair (pair_level records its calls, _pair_ops returns the list),
    the wrapped pair_level simulated symbolically, the D1 program kept where it is cheaper, the n = 2 grouping."""
    ok = True
    calls = []
    fake = types.ModuleType('dqi_pub_pair')
    fake.pair_level = lambda b, prog, words, d1, d2, T1, T2: calls.append(len(words))
    fake._pair_ops = lambda ops: list(ops)
    saved_mod = sys.modules.get('dqi_pub_pair')
    sys.modules['dqi_pub_pair'] = fake
    try:
        g = {'_D1_PAIR_OPS': fake._pair_ops}
        holder = types.FunctionType(compile('def runops(b, ops):\n    return _D1_PAIR_OPS(ops)\n', '<t>', 'exec')
                                    .co_consts[0], g, 'runops')
        log = _install_pub_pair()
        mixed = types.SimpleNamespace(counts=lambda n: dict(total=67 * n + 34))
        cheap = types.SimpleNamespace(counts=lambda n: dict(total=100))
        for n in (2, 3):
            rec = _Rec()
            del calls[:]
            d, w, T, N = _pair_layout(n)
            fake.pair_level(rec, mixed, w, d[0], d[1], T[0], T[1])
            bad = _symbolic(rec.ops, N, [q for r in d for q in r], [q for r in w for q in r], _pair_want(d, w, T, n))
            ccx = sum(1 for o in rec.ops if len(o) == 3)
            good = not bad and ccx == KCOST[n] and not calls
            ok &= good
            print('pair_level wrapped, n = %d: %d CCX (D1 mixed %d), D1 called %s, symbolic %s' % (
                n, ccx, 67 * n + 34, bool(calls), 'PASS' if good else 'FAIL %s' % bad[:3]))
        for n, prog, label in ((4, mixed, 'n = 4 (D1 302 < 308)'), (2, cheap, 'n = 2, a program cheaper than 154'),
                               (5, mixed, 'n = 5')):
            rec = _Rec()
            del calls[:]
            d, w, T, N = _pair_layout(n)
            fake.pair_level(rec, prog, w, d[0], d[1], T[0], T[1])
            good = not rec.ops and calls == [n]
            ok &= good
            print('pair_level wrapped, %s: left to D1 (%s): %s' % (label, calls, 'PASS' if good else 'FAIL'))
        ops = _kit_like_ops(2) + [('c', 'conv')]
        res = holder(None, ops)
        kinds = ''.join(o[0] for o in res)
        p = [o for o in res if o[0] == 'P']
        good = kinds == 'ssPc' and len(p) == 1 and len(p[0][1]) == 2
        ops3 = _kit_like_ops(3) + [('c', 'conv')]
        kinds3 = ''.join(o[0] for o in holder(None, ops3))
        good &= kinds3 == 'ssmmmmmmc'
        ok &= good
        print('pre-pass wrapped (%s): n = 2 run -> %s; n = 3 run (left to D1\'s own grouping) -> %s: %s' % (
            log, kinds, kinds3, 'PASS' if good else 'FAIL'))
        again = _install_pub_pair()
        good = again['pair_level'] == 'already' and again['pair_ops_globals'] == ['already']
        ok &= good
        print('second install: %s: %s' % (again, 'PASS' if good else 'FAIL'))
    finally:
        if saved_mod is None:
            sys.modules.pop('dqi_pub_pair', None)
        else:
            sys.modules['dqi_pub_pair'] = saved_mod
    return ok


def _kittest(root, shots=256, seed=82):
    """T2 on a kit copy: kickmix, the kit's own field.py (its mul_xor0) and the REAL runops of the kit's
    safe_scalar/converter.py (compiled from the file, never executed as a module): the kit runops and the swapped one on
    the same op lists (the 2 n 'm' ops of pre, forward and undo order), values on random shots against this file's
    arithmetic, phase flips, Toffoli-like counts, qubits."""
    import importlib.util
    import numpy as np
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    spec = importlib.util.spec_from_file_location('prek_kit_field',
                                                  os.path.join(root, 'playbook/kit/decoder/exp/arithmetic/field.py'))
    field = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(field)
    p = os.path.join(root, 'playbook/kit/decoder/exp/safe_scalar/converter.py')
    with open(p) as fh:
        co = compile(fh.read(), p, 'exec')
    code = [c for c in co.co_consts if isinstance(c, types.CodeType) and c.co_name == 'runops'][0]
    gl = {'ar': field, '__name__': 'prek_kit_copy', 'scale': None, 'convolution_xor': None}
    fn = types.FunctionType(code, gl, 'runops')
    log = install(km, None)
    print('install:', json.dumps(log), flush=True)
    swapped = fn.__code__ is not code
    ref = types.FunctionType(code, gl, 'runops')
    rng = np.random.default_rng(seed)
    ok = swapped
    for n in (2, 3, 4):
        for undo in (False, True):
            res = {}
            vals = rng.integers(0, 4096, size=(shots, 3 * n + 2), dtype=np.int64)
            for tag, f in (('kit', ref), ('new', fn)):
                b = km.CircuitBuilder()
                g = [b.alloc_qubits(12) for _ in range(n)]
                pay = [[b.alloc_qubits(12) for _ in range(n + 1)] for _ in range(2)]
                ops = [('m', g[i], v[-1], v[i]) for v in pay for i in range(n)]
                f(b, ops, None, None, undo)
                c = b.finish_circuit()
                words = g + [pay[k][-1] for k in range(2)] + [pay[k][i] for k in range(2) for i in range(n)]
                qa = km.array([q for w in words for q in w])
                sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed)
                sim.use_same_registers_as(c)
                bits = ((vals[:, :, None] >> np.arange(12)) & 1 == 1).transpose(1, 2, 0).reshape(-1, shots)
                sim.write_across_shots(qa, bits)
                sim.do(c)
                out = np.asarray(sim.read_across_shots(qa)).reshape(3 * n + 2, 12, -1).astype(np.int64)
                out = (out * (1 << np.arange(12))[None, :, None]).sum(axis=1).T
                res[tag] = (out, int(np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool).sum()),
                            c.operation_statistics().num_toffoli_like, c.num_qubits)
            want = vals.copy()
            for k in range(2):
                for i in range(n):
                    want[:, n + 2 + k * n + i] ^= np.array([gmul(int(x), int(y)) for x, y in
                                                            zip(vals[:, i], vals[:, n + k])], dtype=np.int64)
            e_kit = int((res['kit'][0] != want).any(axis=1).sum())
            e_new = int((res['new'][0] != want).any(axis=1).sum())
            good = e_kit == 0 and e_new == 0 and res['kit'][1] == 0 and res['new'][1] == 0 \
                and res['kit'][2] == 84 * n and res['new'][2] == KCOST[n] \
                and res['kit'][3] == res['new'][3] == 12 * (3 * n + 2)
            ok &= good
            print('  n = %d undo %s: Toffoli-like kit %d, new %d; qubits %d / %d; wrong shots of %d: kit %d, new %d; '
                  'phase flips %d / %d: %s' % (n, undo, res['kit'][2], res['new'][2], res['kit'][3], res['new'][3],
                                               shots, e_kit, e_new, res['kit'][1], res['new'][1],
                                               'PASS' if good else 'FAIL'), flush=True)
    print('KITTEST %s (runops of %s swapped: %s)' % ('PASS' if ok else 'FAIL', p, swapped), flush=True)
    return ok


def _pairtest(spec, shots=256, seed=83):
    """T2 of the D1 path: the REAL tools/dqi_pub_pair.py with the line's program (spec as DQI_PUB_PAIR_PROG), its
    original pair_level against the wrapped one on kickmix, n = 2, 3, 4 words: values on random shots against this
    file's arithmetic, phase flips, Toffoli-like counts (the wrapped one must read 154 / 224 / D1's own at n = 4)."""
    import numpy as np
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    import dqi_pub_pair as pp
    prog = pp.load_program(spec)
    orig = pp.pair_level
    _load()
    info = _install_pub_pair()
    print('install on the real dqi_pub_pair:', info, flush=True)
    rng = np.random.default_rng(seed)
    ok = getattr(pp.pair_level, '_b5_prek', False)
    for n in (2, 3, 4):
        res = {}
        vals = rng.integers(0, 4096, size=(shots, 2 + 3 * n), dtype=np.int64)
        for tag, f in (('d1', orig), ('new', pp.pair_level)):
            b = km.CircuitBuilder()
            d = [b.alloc_qubits(12) for _ in range(2)]
            w = [b.alloc_qubits(12) for _ in range(n)]
            T = [[b.alloc_qubits(12) for _ in range(n)] for _ in range(2)]
            f(b, prog, w, d[0], d[1], T[0], T[1])
            c = b.finish_circuit()
            words = d + w + T[0] + T[1]
            qa = km.array([q for x in words for q in x])
            sim = km.Simulator(batch_size=shots, debug_print=False, counting_operations=True, seed=seed)
            sim.use_same_registers_as(c)
            bits = ((vals[:, :, None] >> np.arange(12)) & 1 == 1).transpose(1, 2, 0).reshape(-1, shots)
            sim.write_across_shots(qa, bits)
            sim.do(c)
            out = np.asarray(sim.read_across_shots(qa)).reshape(2 + 3 * n, 12, -1).astype(np.int64)
            out = (out * (1 << np.arange(12))[None, :, None]).sum(axis=1).T
            res[tag] = (out, int(np.asarray(sim.read_phase_flipped_across_shots(), dtype=bool).sum()),
                        c.operation_statistics().num_toffoli_like, c.num_qubits)
        want = vals.copy()
        for k in range(2):
            for i in range(n):
                want[:, 2 + n + k * n + i] ^= np.array([gmul(int(x), int(y)) for x, y in
                                                        zip(vals[:, 2 + i], vals[:, k])], dtype=np.int64)
        e_d1 = int((res['d1'][0] != want).any(axis=1).sum())
        e_new = int((res['new'][0] != want).any(axis=1).sum())
        theirs = prog.counts(n)['total']
        expect = KCOST[n] if (n in (2, 3) and KCOST[n] < theirs) else theirs
        good = e_d1 == 0 and e_new == 0 and res['d1'][1] == 0 and res['new'][1] == 0 and res['d1'][2] == theirs \
            and res['new'][2] == expect and res['d1'][3] == res['new'][3]
        ok &= good
        print('  n = %d: Toffoli-like D1 %d (its count %d), new %d; qubits %d / %d; wrong shots of %d: D1 %d, new %d; '
              'phase flips %d / %d: %s' % (n, res['d1'][2], theirs, res['new'][2], res['d1'][3], res['new'][3], shots,
                                           e_d1, e_new, res['d1'][1], res['new'][1], 'PASS' if good else 'FAIL'),
              flush=True)
    print('PAIRTEST %s' % ('PASS' if ok else 'FAIL'), flush=True)
    return ok


def _kit_like_ops(n):
    g = [[12 * j + l for l in range(12)] for j in range(n)]
    base = 12 * n
    pay = [[[base + 12 * (k * (n + 1) + j) + l for l in range(12)] for j in range(n + 1)] for k in range(2)]
    return [('s', v[-1], True) for v in pay] + [('m', g[i], v[-1], v[i]) for v in pay for i in range(n)]


if __name__ == '__main__':
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    import importlib
    me = importlib.import_module(NAME)
    if len(sys.argv) > 1 and sys.argv[1] == 'test':
        sys.exit(0 if me._test() else 1)
    if len(sys.argv) > 2 and sys.argv[1] == 'kittest':
        sys.exit(0 if me._kittest(os.path.expanduser(sys.argv[2])) else 1)
    if len(sys.argv) > 2 and sys.argv[1] == 'pairtest':
        sys.exit(0 if me._pairtest(sys.argv[2]) else 1)
    raise SystemExit('usage: dqi_int_prek.py test | kittest <kit copy root> | pairtest <DQI_PUB_PAIR_PROG spec>')
