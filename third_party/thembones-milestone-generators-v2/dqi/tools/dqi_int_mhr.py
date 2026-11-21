"""dqi_int_mhr.py -- B5 (21st holder, Nov 17): the kit's high_full_restore with its top head word made by the conv.

The kit's restore (narrow240/high_stream.py high_full_restore, every clone of it, every caller):
    embed (q0, q1) as q or x q under the gap flag F (dk.shifted_pair: 24 CSWAP in, 12 CSWAP out);
    conv at m = len(g) points into f[:-2] + [lo]  (m K-products: the product mod prod(x - t), t < m);
    q0 measured into masks[0], lo swapped in;
    hi = sel(g[-1], g[-2]; F) x q1  (selected: 12 ANDs; one K-product: 42);
    q1 measured into masks[1], hi swapped in; the m low words corrected by hi x (x^m mod prod(x - t)) (CNOTs).
Here (variant X): the conv runs at m + 1 points into f[:-2] + [lo, hi]. The product g x q_emb has degree <= m on
every input (its coefficient m + 1 is g[-1] w2, and F = 1 means g[-1] = 0 while F = 0 means w2 = 0), so m + 1 points
give every coefficient exactly, hi included: hi = g[-1] w1 + g[-2] w2 = the kit's sel x q1. Same masks, same order of
the two measurements, same values in f afterwards: the callers' phase algebra is untouched.
Count [P]: -selected (12) - one K-product (42) + one more conv point (42) = -12 a call. Width: hi is live during the
conv (+12 wires at the conv against the kit's restore).
Install: anywhere in DQI_INT_EXT (it swaps the CODE of the kit's function objects, so every wrapper and scope that
holds them, pdl2's hook included, runs the new body). Switches:
    DQI_INT_MHR_OFF=1        install nothing.
    DQI_INT_MHR_LEVELS=a,b   only levels m = len(g) in this set take the new body (default: all).
    DQI_INT_MHR_OUT=<file>   at exit, a JSON record: per call (m, caller, live qubits before / at the conv, mode).
    DQI_INT_MHR_E=1          variant E at the deferred callers (high_stream direct_post, a16 / a22 bottom_post): the
                             swap-out goes too, the three embedded words are measured, and the kit direct_post's
                             reverse pays the F part of their debt by CZ_IF (source swap of its repayment line).
    DQI_INT_MHR_E_LEVELS=..  limits E; DQI_INT_MHR_XOFF=1: callers that do not take E keep the kit's body.
Measured (Nov 17, box 6): 18 restore calls a decode; every X / E conv sits at 902 in the generator (866-867 live
before the restore), level 16 at 914 and level 6 at 910: each new 902 moment needs a park in the 901 chain.
"""
import atexit
import json
import os
import sys
import types

LOG = dict(patched=0, calls=0, new=0, old=0, levels={}, callers={})
CALLS = []
_ORIG_CODE = {}

_SRC = '''
def high_full_restore(b,f,g,F):
 m=len(g);mode=_b5_mhr_on(b,m)
 if not mode:return _b5_mhr_old(b,f,g,F)
 _b5_mhr_cap(12)
 lo=b.alloc_qubits(12);hi=b.alloc_qubits(12)
 if mode==2:
  z=b.alloc_qubits(12)
  b.cswap(F,f[-1],z);b.cswap(F,f[-2],f[-1])
  early=_b5_mhr_early(b,F,m)
  _b5_mhr_mark(b,m)
  dk.conv(b,g,[f[-2],f[-1],z],f[:-2]+[lo,hi],m+1)
  masks=[b.alloc_bits(12),b.alloc_bits(12),b.alloc_bits(12)]
  b.hmr(f[-2],masks[0]);b.hmr(f[-1],masks[1]);b.hmr(z,masks[2]);b.free(z)
  b.swap(f[-2],lo);b.swap(f[-1],hi)
  b.free(hi);b.free(lo)
  _b5_mhr_relist(b,F if early else None)
  _b5_mhr_cap(-12)
  return masks
 with dk.shifted_pair(b,f[-2:],F) as q:
  _b5_mhr_mark(b,m)
  dk.conv(b,g,q,f[:-2]+[lo,hi],m+1)
 masks=[b.alloc_bits(12),b.alloc_bits(12)]
 b.hmr(f[-2],masks[0]);b.swap(f[-2],lo)
 b.hmr(f[-1],masks[1]);b.swap(f[-1],hi)
 b.free(hi);b.free(lo)
 _b5_mhr_relist(b)
 _b5_mhr_cap(-12)
 return masks
'''

# kickmix reuses freed qubits last in, first out, and one alloc / free cycle reverses the order inside a block
# [M kickmix, Nov 17 18:58Z]. The kit's restore cycles each of its two blocks twice (lo then hi on the same wires, z then
# the selected word), which leaves the free list as it found it; this restore cycles three blocks once, which leaves
# each block reversed, so every later allocation (the root's flags among them) lands on other wires than in the base
# and the 901 spec's parks stop fitting (E3, E5: the root's two runs). RELIST (default on) runs one more cycle of the
# three blocks (alloc / free only, no gate) and restores the list exactly. DQI_INT_MHR_RELIST=0 turns it off.
RELIST = os.environ.get('DQI_INT_MHR_RELIST', '1') == '1'


def _relist(b, F=None):
    """Restore the free list (see RELIST). With an early-held F (EARLYF), F's wire sits right below the three blocks:
    it is re-taken here (allocated, at |0>) and freed at the caller's own hold, so the list is the base recipe's
    exactly [M kickmix, scratchpad kmhold.py: same pickup wire, same held wire, same 45-qubit probe]."""
    if not RELIST and F is None:
        return
    x = b.alloc_qubits(12)
    y = b.alloc_qubits(12)
    w = b.alloc_qubits(12)
    v = b.alloc_qubits(1) if F is not None else None
    b.free(w)
    b.free(y)
    b.free(x)
    LOG['relist'] = LOG.get('relist', 0) + 1
    if F is not None:
        EARLY[id(F)] = (v, EARLY_MASK.pop(id(F)))


# EARLYF (variant E only): the level's gap flag F is not read between E's entry and look_solve's hold_flag(b, F) right
# after direct_post returns (direct_post's tail, pdr's pickup and a07's wrapper do not read it), and the kit's hold is
# hmr(F) into a new bit + free(F). So the hold is done at E's entry instead: the conv then needs no 902nd wire, and
# the caller's hold returns that same bit (the wrapper below). Skipped where F is kept (m in KEEP_FULL: held_gap).
EARLYF = os.environ.get('DQI_INT_MHR_EARLYF', '0') == '1'
EARLY = {}
EARLY_MASK = {}


_el = os.environ.get('DQI_INT_MHR_EARLY_LEVELS', '')
EARLY_LEVELS = None if not _el else frozenset(int(x) for x in _el.split(',') if x)


def _early(b, F, m):
    gl = sys._getframe(1).f_globals
    if EARLYF and EARLY_LEVELS is not None and m not in EARLY_LEVELS:
        return False
    if not EARLYF or m in gl.get('KEEP_FULL', frozenset()):
        if EARLYF:
            LOG['early_skipped_full'] = LOG.get('early_skipped_full', 0) + 1
        return False
    mk = b.alloc_bits(1)
    b.hmr(F, mk[0])
    b.free(F)
    EARLY_MASK[id(F)] = mk
    LOG['early'] = LOG.get('early', 0) + 1
    LOG.setdefault('early_levels', []).append(m)
    return True


def _make_hold(orig):
    def hold_flag(b, F, *a, **kw):
        e = EARLY.pop(id(F), None)
        if e is None:
            return orig(b, F, *a, **kw)
        v, mk = e
        b.free(v)
        LOG['early_held'] = LOG.get('early_held', 0) + 1
        _zp_hold_tick(b, F, mk)
        return mk
    hold_flag._b5_mhr = True
    return hold_flag


# 24th holder (Nov 18 21:10Z): the early-held hold never reached zpark's own hold wrapper, so zpark's call clock (its
# oracle's time base, B5_ZP_ORACLE) ran one kit call behind from that hold on: zd7's build read b5_clock 905 against
# 906 and b5_k1_miss 55 against 0, victims 16 against 25 [M build logs] and P1 lost 37 links in unconvert (zero-map diff,
# tools/dqi_int_zdiff.py). ZPTICK (default on) does here exactly what zpark's hold does for a flag it does not park:
# tick its clock, drop the flag from its live map, file the mask under the flag's key (its FOLD / consume bookkeeping).
# No gate is emitted. DQI_INT_MHR_ZPTICK=0 restores the 21st holder's behaviour.
ZPTICK = os.environ.get('DQI_INT_MHR_ZPTICK', '1') == '1'


def _zp_hold_tick(b, F, mk):
    if not ZPTICK:
        return
    for name in ('dqi_int_zpk2', 'dqi_int_zpk1', 'dqi_int_zpk0', 'dqi_int_zpk', 'r46_dqi_conv_zpark'):
        zp = sys.modules.get(name)
        if zp is not None and callable(getattr(zp, 'state', None)):
            break
    else:
        LOG['zp_tick_nomodule'] = LOG.get('zp_tick_nomodule', 0) + 1
        return
    st = zp.state(b)
    if callable(getattr(zp, 'tick', None)):
        zp.tick()
    act = zp.active(st) if callable(getattr(zp, 'active', None)) else True
    key = st['live'].pop(getattr(F, 'id', None), None) if act else None
    if key is not None:
        st['masks'][mk[0].id] = (key, st['epoch'].get(key, 0), mk)
    LOG['zp_ticks'] = LOG.get('zp_ticks', 0) + 1
    LOG['zp_tick_keyed'] = LOG.get('zp_tick_keyed', 0) + (key is not None)

# zpark (dqi_int_zpk2 / r46_dqi_conv_zpark) evicts a parked gap flag whenever an allocation would pass its module CAP;
# the extra word of X / E would trigger evictions at every restore and change zpark's flag hand-overs downstream (E3:
# the root's root_child_post runs lost their idle flag and the 901 spec could not cover them). CAPUP (default on)
# raises that cap by 12 for the restore's own duration: zpark's state stays the base recipe's, and the 901 chain parks
# the restore moments themselves. DQI_INT_MHR_CAPUP=0 turns it off.
CAPUP = os.environ.get('DQI_INT_MHR_CAPUP', '1') == '1'


def _cap(d):
    if not CAPUP:
        return
    for name in ('dqi_int_zpk2', 'dqi_int_zpk1', 'dqi_int_zpk0', 'dqi_int_zpk', 'r46_dqi_conv_zpark'):
        mod = sys.modules.get(name)
        if mod is not None and isinstance(getattr(mod, 'CAP', None), int):
            mod.CAP += d
            LOG['capup'] = LOG.get('capup', 0) + (1 if d > 0 else 0)

# Variant E's repayment: the kit's direct_post (high_stream.py) reverse pays the deferred heads by
#   for word,mask in zip([q0,q1],phase['heads'].pop(m)): b.cz(word,mask);b.free(mask)
# Its source line is replaced by _b5_mhr_pay(...), which pays a 2-mask entry the kit's way and a 3-mask entry
# (the embedded words measured: debt m0 q0 + m1 q1 + F ((m0 + m1) q0 + (m1 + m2) q1)) by CZ_IF gates only.
_POP_OLD = ("  for word,mask in zip([q0,q1],phase['heads'].pop(m)):\n"
            "   b.cz(word,mask);b.free(mask)\n")
_POP_NEW = "  _b5_mhr_pay(b,q0,q1,F,phase['heads'].pop(m))\n"
E_ON = os.environ.get('DQI_INT_MHR_E', '0') == '1'
E_CALLERS = ('direct_post', 'bottom_post')
_E_READY = []


def _pay(b, q0, q1, F, ms):
    if len(ms) == 2:
        for word, mask in zip([q0, q1], ms):
            b.cz(word, mask)
            b.free(mask)
        return
    LOG['e_paid'] = LOG.get('e_paid', 0) + 1
    m0, m1, m2 = ms
    b.cz(q0, m0)
    b.cz(q1, m1)
    for j in range(12):
        b.ccz(F, q0[j], m0[j])
        b.ccz(F, q0[j], m1[j])
        b.ccz(F, q1[j], m1[j])
        b.ccz(F, q1[j], m2[j])
    for mk in ms:
        b.free(mk)


def _levels():
    s = os.environ.get('DQI_INT_MHR_LEVELS', '')
    return None if not s else frozenset(int(x) for x in s.split(',') if x)


LEVELS = _levels()
_es = os.environ.get('DQI_INT_MHR_E_LEVELS', '')
E_LEVELS = None if not _es else frozenset(int(x) for x in _es.split(',') if x)
X_OFF = os.environ.get('DQI_INT_MHR_XOFF', '0') == '1'     # callers that do not take E keep the kit's body


def _live(b):
    try:
        v = b.num_allocated_qubits
        return v() if callable(v) else v
    except Exception:
        return -1


def _caller(depth=3):
    try:
        fr = sys._getframe(depth)
        return '%s:%s' % (os.path.basename(fr.f_code.co_filename), fr.f_code.co_name)
    except Exception:
        return '?'


def _restore_caller():
    """The first frame above the restore that is neither the restore nor pdl2's wrapper of it."""
    d = 2
    while d < 9:
        try:
            co = sys._getframe(d).f_code
        except ValueError:
            return None
        if co.co_name == 'high_full_restore' or co.co_filename.endswith('dqi_zn_pdl2.py'):
            d += 1
            continue
        return co
    return None


def _on(b, m):
    LOG['calls'] += 1
    co = _restore_caller()
    who = '%s:%s' % (os.path.basename(co.co_filename), co.co_name) if co is not None else '?'
    LOG['callers'][who] = LOG['callers'].get(who, 0) + 1
    mode = 1 if (LEVELS is None or m in LEVELS) else 0
    if mode and E_ON and _E_READY and co is not None and co.co_name in E_CALLERS \
            and (E_LEVELS is None or m in E_LEVELS):
        mode = 2
    elif mode == 1 and X_OFF:
        mode = 0
    LOG[('old', 'new', 'newE')[mode]] = LOG.get(('old', 'new', 'newE')[mode], 0) + 1
    key = '%d%s' % (m, ('k', '', 'E')[mode])
    LOG['levels'][key] = LOG['levels'].get(key, 0) + 1
    CALLS.append([m, who, _live(b), None, mode])
    return mode


def _mark(b, m):
    if CALLS:
        CALLS[-1][3] = _live(b)


def _dump():
    out = os.environ.get('DQI_INT_MHR_OUT')
    if not out:
        return
    rec = dict(LOG, calls_list=CALLS, early_left=len(EARLY) + len(EARLY_MASK))
    pdr = sys.modules.get('dqi_z2_pdr')            # diagnostics only: Z2's pdr counters (pickups, rebuilds, ...)
    if pdr is not None:
        try:
            rec['pdr_log'] = json.loads(json.dumps(getattr(pdr, 'LOG', {}), default=str))
            st = getattr(pdr, 'ST', {})
            rec['pdr_pend_left'] = sorted(map(str, st.get('pend', {}).keys())) if isinstance(st, dict) else None
        except Exception as e:
            rec['pdr_err'] = repr(e)[:200]
    try:
        with open(out, 'w') as fh:
            json.dump(rec, fh, default=str)
    except Exception:
        pass


def _is_kit_restore(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'high_full_restore'
            and co.co_filename.endswith('high_stream.py') and not getattr(fn, '_b5_mhr', False))


def _patch(fn):
    """Swap the code of one kit function object; its globals (the kit scope) get the helpers."""
    gl = fn.__globals__
    if fn.__closure__ or fn.__code__.co_freevars:
        LOG.setdefault('skipped_closure', 0)
        LOG['skipped_closure'] += 1
        return
    old = types.FunctionType(fn.__code__, gl, 'high_full_restore_kit', fn.__defaults__, fn.__closure__)
    gl['_b5_mhr_old'] = old
    gl['_b5_mhr_on'] = _on
    gl['_b5_mhr_mark'] = _mark
    gl['_b5_mhr_cap'] = _cap
    gl['_b5_mhr_relist'] = _relist
    gl['_b5_mhr_early'] = _early
    ns = {}
    exec(compile(_SRC, fn.__code__.co_filename, 'exec'), gl, ns)
    new_code = ns['high_full_restore'].__code__.replace(co_firstlineno=fn.__code__.co_firstlineno)
    _ORIG_CODE[id(fn)] = fn.__code__
    fn.__code__ = new_code
    fn._b5_mhr = True
    LOG['patched'] += 1


POSTS = []


def _is_kit_post(fn):
    co = getattr(fn, '__code__', None)
    return (isinstance(fn, types.FunctionType) and co is not None and co.co_name == 'direct_post'
            and co.co_filename.endswith('high_stream.py') and not getattr(fn, '_b5_mhr', False))


def _patch_post(fn):
    """Variant E: swap the code of one kit direct_post for its own source with the heads repayment replaced."""
    import inspect
    LOG['posts_found'] = LOG.get('posts_found', 0) + 1
    if fn.__closure__ or fn.__code__.co_freevars:
        return
    try:
        src = inspect.getsource(fn)
    except Exception as e:
        LOG.setdefault('post_errors', []).append(repr(e)[:120])
        return
    if src.count(_POP_OLD) != 1:
        LOG.setdefault('post_errors', []).append('pop line not found once')
        return
    gl = fn.__globals__
    gl['_b5_mhr_pay'] = _pay
    ns = {}
    pad = '\n' * (fn.__code__.co_firstlineno - 1)
    exec(compile(pad + src.replace(_POP_OLD, _POP_NEW), fn.__code__.co_filename, 'exec'), gl, ns)
    newf = ns.get(fn.__code__.co_name)
    if newf is None or newf.__code__.co_freevars:
        LOG.setdefault('post_errors', []).append('compile')
        return
    fn.__code__ = newf.__code__
    fn._b5_mhr = True
    LOG['posts_patched'] = LOG.get('posts_patched', 0) + 1


def _visit(v, seen):
    if id(v) in seen:
        return
    seen.add(id(v))
    if _is_kit_restore(v):
        _patch(v)
    elif E_ON and _is_kit_post(v):
        POSTS.append(v)


def _cells(fn):
    """Contents of a Python function's closure cells (empty cells skipped)."""
    out = []
    if not isinstance(fn, types.FunctionType):
        return out
    for cl in (fn.__closure__ or ()):
        try:
            out.append(cl.cell_contents)
        except Exception:
            pass
    return out


def install(km, ns):
    if os.environ.get('DQI_INT_MHR_OFF', '0') == '1':
        return dict(LOG, off=True)
    walk = sys.modules['walk_bound']
    seen = set()
    for fn in walk.functions(ns.decode):
        if not isinstance(fn, types.FunctionType):
            continue
        _visit(fn, seen)
        cells = _cells(fn)
        ds = [fn.__globals__] + [c for c in cells if isinstance(c, dict)]
        for c in cells:
            if isinstance(c, types.FunctionType):
                _visit(c, seen)
        for d in ds:
            try:
                vals = list(d.values())
            except Exception:
                continue
            for v in vals:
                if not isinstance(v, types.FunctionType):
                    continue
                _visit(v, seen)
                # a wrapper (pdl2's make_restore) holds the kit's function in its closure as 'orig'
                for c in _cells(v):
                    if isinstance(c, types.FunctionType):
                        _visit(c, seen)
    if EARLYF:                                     # outermost wrapper of every scope's hold_flag (installed last)
        seen2 = set()
        for fn in walk.functions(ns.decode):
            if not isinstance(fn, types.FunctionType):
                continue
            ds = [fn.__globals__] + [c for c in _cells(fn) if isinstance(c, dict)]
            for d in ds:
                if id(d) in seen2:
                    continue
                seen2.add(id(d))
                h = d.get('hold_flag') if isinstance(d, dict) else None
                if isinstance(h, types.FunctionType) and not getattr(h, '_b5_mhr', False):
                    d['hold_flag'] = _make_hold(h)
                    LOG['hold_wraps'] = LOG.get('hold_wraps', 0) + 1
    for fn in POSTS:
        _patch_post(fn)
    if E_ON and POSTS and LOG.get('posts_patched', 0) == len(POSTS):
        _E_READY.append(True)
    atexit.register(_dump)
    LOG['want_levels'] = sorted(LEVELS) if LEVELS is not None else 'all'
    LOG['e_ready'] = bool(_E_READY)
    return dict(patched=LOG['patched'], want_levels=LOG['want_levels'], e_ready=LOG['e_ready'],
                posts=LOG.get('posts_patched', 0), post_errors=LOG.get('post_errors'))
