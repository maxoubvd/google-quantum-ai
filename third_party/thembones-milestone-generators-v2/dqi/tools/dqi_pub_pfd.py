"""dqi_pub_pfd.py -- builder D1: the prefix block (J = 5) with its forward erasure deferred to the return.

The build line's block (tools/dqi_build_patch_prefix5.py, read-only) is called identically on both trips (it is an
involution on the later groups' L words). Each call computes U at the points alpha = 1, 4, infinity into three clean
words and at alpha = 0 in place, applies the 14 group products, and then erases the three clean U words by
recomputing Y and G at each point (2 + 2 + 3 = 7 twelve-bit products, 294 Toffolis).
Deferral: on the FORWARD call the three clean U words are X-measured instead (masks m_p kept as classical bits).
The head words that U is made of (heads 1..4) hold the same values when the RETURN call runs, so the return call,
right after it has recomputed U(alpha_p) into the same kind of word, pays the forward's phase (-1)^(m_p . U_p) with
a classically controlled Z on that word (free), and then erases its own U words as before. Point 0 (in place in the
head words) is still undone on both trips. The skipped recomputation allocated and freed two 12-qubit transients
per point; kickmix reuses freed qubits last-in first-out per qubit, so each such cycle reverses that chunk of the
free stack: the forward replays the same allocate / free sequence without gates, which keeps every later placement
identical (without the replay the final routing failed on q21, build p2pfd).
Measured: p2pfdb = p1a1f + pfd, 360,702.8 executed (quick), -293.9 against p1a1f; static -294; Q 902.
Implementation: the source of tools/dqi_build_patch_prefix5.py is executed in a private module with two exact text
replacements (asserted to match once each); nothing in the original file changes.
Install (from tools/dqi_pub_entry.py): DQI_PUB_P=...,pfd  (after B5's decode has installed prefix5).
"""
import os
import sys
import types

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, 'dqi_build_patch_prefix5.py')

_PAY_OLD = ("            point0(True)\n"
            "            R1, R2, R3 = U[1], U[4], U[None]\n")
_PAY_NEW = ("            if not _p0:\n"
            "                point0(True)\n"
            "                _B1 = b1\n"
            "            _pfd_pay(b, U, b1)\n"
            "            if _PFD_ZJOIN:\n"
            "                _zX, _zL, _za0 = words[68], words[69], words[0]\n"
            "                b.cx(_za0, _zX); square12(b, _zX)\n"
            "                _zc = mul(sq(k0), sq(k1))\n"
            "                scale(b, _zL, _zc); ar.mul_xor0(b, _B1, _zX, _zL); scale(b, _zL, inv(_zc))\n"
            "                square12(b, _zX, True); b.cx(_za0, _zX)\n"
            "                if not _p0:\n"
            "                    cx(b, b4, _zL, inv(sq(k1)))\n"
            "            R1, R2, R3 = U[1], U[4], U[None]\n")
# DQI_PUB_PFD_P0=1: on the FORWARD call only, the point alpha = 0 runs out of place and first: Y(0), G(0), U(0) into
# three clean words (3 products), Y(0) and G(0) erased by measurement at once (free: their repair is bilinear in live
# words), U(0) kept for the group loop and X-measured at the end like the other U words; the return, whose point 0
# stays in place (b_1 holds U(0) there), pays the mask on b_1. Forward point 0: 6 products -> 3 (-126 [P count]).
# The forward block peaks at 60 clean qubits instead of 48 (12 wires created early, at high-water 888 -> 900); at the
# block's exit the free stack is set to the one the unchanged line has there (recorded by a scratch decode at
# install), with the early wires at the bottom, smallest index on top: kickmix pops before it creates and creates in
# index order, so each early wire is taken exactly where the unchanged line creates it (wires 888..901 are created
# later in the same forward extract [M, S3's q902s trace]).
_START_OLD = ("            U = {}\n"
              "            for p in TEMP_POINTS:\n")
_START_NEW = ("            U = {}\n"
              "            _p0 = _pfd_p0(b)\n"
              "            if _p0:\n"
              "                _k0s, _k1s = sq(k0), sq(k1)\n"
              "                _Y0 = b.alloc_qubits(12)\n"
              "                b.cx(a3, _Y0); cx(b, b3, _Y0, k1); sqmul(a3, b4, _Y0)\n"
              "                _G0 = b.alloc_qubits(12)\n"
              "                cx(b, b2, _G0, _k1s); cx(b, b4, _G0, _k1s); sqmul(a2, _Y0, _G0)\n"
              "                _B1 = b.alloc_qubits(12)\n"
              "                cx(b, b1, _B1, mul(_k1s, k0)); cx(b, a1, _B1, _k1s); cx(b, _Y0, _B1, _k0s)\n"
              "                sqmul(a1, _G0, _B1)\n"
              "                if _PFD_ZJOIN:\n"
              "                    cx(b, _G0, words[69], inv(_k1s))\n"
              "                cx(b, b2, _G0, _k1s); cx(b, b4, _G0, _k1s)\n"
              "                square12(b, a2); ar.erase_product(b, _Y0, a2, _G0); square12(b, a2, True)\n"
              "                b.free(_G0)\n"
              "                b.cx(a3, _Y0); cx(b, b3, _Y0, k1)\n"
              "                square12(b, a3); ar.erase_product(b, b4, a3, _Y0); square12(b, a3, True)\n"
              "                b.free(_Y0)\n"
              "            for p in TEMP_POINTS:\n")
_LOOP_OLD = [("                cx(b, b1, R1, lam[0][0]); cx(b, b1, R2, lam[0][1])\n",
              "                cx(b, _B1, R1, lam[0][0]); cx(b, _B1, R2, lam[0][1])\n"),
             ("                cx(b, b1, R2, lam[0][1]); cx(b, b1, R1, lam[0][0])\n",
              "                cx(b, _B1, R2, lam[0][1]); cx(b, _B1, R1, lam[0][0])\n")]
P0 = os.environ.get('DQI_PUB_PFD_P0', '0') == '1'
P0_STATE = dict(active=P0, recording=False, order=None, fixes=0, misses=0, extras=0)
# DQI_PUB_ZJOIN=1 (needs the zt module's matching switch): the node z = 0 joins the block. At point 0 the block
# holds U(0) in b_1 and G(0) in b_4, and z = 0 is a later node with w = 0, kappa_g = 0, so the zero node's L word
# gets U(0) (x_0 + a_0)^2 / (k0^2 k1^2) + G(0) / k1^2 (one 12-bit product) and its chain runs heads 1..4 without
# their L products (4 x 42 per trip). Identity: tools/dqi_pub_zcheck.py, 500 / 500 [M classical].
ZJOIN = os.environ.get('DQI_PUB_ZJOIN', '0') == '1'
_END_OLD = ("            point0(False)\n"
            "            for p in reversed(TEMP_POINTS):\n"
            "                temp_point(p, U[p], False, p == TEMP_POINTS[-1])\n"
            "                b.free(U[p])\n")
_END_NEW = ("            if not _p0:\n"
            "                point0(False)\n"
            "            _pfd_end(b, U, temp_point, p0u=_B1 if _p0 else None)\n")
_TP_OLD = "                temp_point(p, U[p], True, p == TEMP_POINTS[-1])\n"
_TP_NEW = "                _pfd_tp(b, temp_point, p, U[p], True, p == TEMP_POINTS[-1])\n"

# DQI_PUB_PFD_RINF=1: on the RETURN call only, the point at infinity uses a clean G word instead of G in place in
# the b_4 register (saves the in-place undo product in its compute and in its erasure: -2 products). The return
# block runs when the circuit's high-water mark is already 901 (S3's q902s trace), so its 60-qubit peak creates no
# wire; the extra 12-qubit G cycle is undone on the free stack by an exact two-chunk reorder (_fix24).
RINF = os.environ.get('DQI_PUB_PFD_RINF', '0') == '1'

STATS = dict(forward=0, reverse=0, rinf=0)
_CALLS = {}       # id(builder) -> number of block calls so far (even: the next call is a forward call)
_MASKS = {}       # id(builder) -> stack of {point: classical mask bits} from forward calls
_DIR = {}         # id(builder) -> direction of the call in progress


def _pfd_start(b):
    """Block entry: the first call of a builder is the forward trip, the second the return."""
    k = id(b)
    n = _CALLS.get(k, 0)
    _CALLS[k] = n + 1
    _DIR[k] = n % 2 == 0


def _fix24(b):
    """After temp_point with a clean G where the original kept G in place: the original popped Y (x1..x12) and
    pushed it back, leaving the next chunk y untouched; the clean version also popped and pushed G (y1..y12), which
    reverses y. Pop both chunks and push them back in the original order (no gates)."""
    reg = b.alloc_qubits(24)                  # pops x12..x1 then y12..y1
    for i in range(12, 24):                   # push y12, ..., y1  (y back in its original order)
        b.free(reg[i])
    for i in reversed(range(12)):             # push x1, ..., x12  (x12 on top, as the original left it)
        b.free(reg[i])


def _pfd_tp(b, temp_point, p, U, compute, inplace):
    if inplace and RINF and not _DIR[id(b)]:
        temp_point(p, U, compute, False)
        _fix24(b)
        STATS['rinf'] += 1
    else:
        temp_point(p, U, compute, inplace)


def _pfd_p0(b):
    """True on a forward call when the forward point 0 runs out of place (DQI_PUB_PFD_P0=1)."""
    return P0_STATE['active'] and _DIR[id(b)]


def _pfd_pay(b, U, b1):
    """Start of the group loop: on a return call pay the forward's measured U words (key 0: U(0), which the return
    holds in place in b_1)."""
    k = id(b)
    if _DIR[k]:
        STATS['forward'] += 1
        return
    STATS['reverse'] += 1
    masks = _MASKS[k].pop()
    for p, m in masks.items():
        b.cz(b1 if p == 0 else U[p], m)   # (-1)^(m_p . U(alpha_p)): the forward's phase, U recomputed here
        b.free(m)


def _qidx(q):
    return int(str(q)[1:])


def _spare(b):
    n = b.finish_circuit().num_qubits - b.num_allocated_qubits
    return list(b.alloc_qubits(n)) if n > 0 else []


def _p0_snapshot(b):
    got = _spare(b)
    order = [_qidx(q) for q in got]                  # stack top first
    for q in reversed(got):
        b.free(q)
    return order


def _p0_restore(b):
    """Forward block exit under P0: the unchanged line's free stack, the wires P0 created early below it."""
    order = P0_STATE['order']
    got = _spare(b)
    have = {_qidx(q): q for q in got}
    extra = sorted(set(have) - set(order or ()))
    if order is None or not set(order) <= set(have) or (extra and order and extra[0] <= max(order)):
        P0_STATE['misses'] += 1
        for q in reversed(got):
            b.free(q)
        print(f'dqi_pub_pfd: P0 stack MISS (recorded {None if order is None else len(order)}, spare {len(have)})',
              file=sys.stderr, flush=True)
        return
    for i in reversed(extra):                        # bottom: the early wires, smallest index on top of them
        b.free(have[i])
    for i in reversed(order):
        b.free(have[i])
    P0_STATE['fixes'] += 1
    P0_STATE['extras'] = len(extra)
    print(f'dqi_pub_pfd: P0 stack restored ({len(order)} wires, {len(extra)} early wires at the bottom)',
          file=sys.stderr, flush=True)


def _pfd_end(b, U, temp_point, points=None, p0u=None):
    """End of the block: forward call measures the three clean U words; return call erases them as before."""
    k = id(b)
    pts = points if points is not None else _MOD.TEMP_POINTS
    if _DIR[k]:
        masks = {}
        for p in reversed(pts):
            # keep the free stack identical to the original erasure: temp_point(p, U, False, inplace) allocates
            # Y (and G when not in place) and frees them last-in first-out; kickmix reuses freed qubits in that
            # order, so the same allocate / free sequence is replayed here without gates (fix round 1, 12:01Z:
            # without it the final routing failed on q21, as a3 v5 did).
            inplace = p == pts[-1]
            Y = b.alloc_qubits(12)
            if not inplace:
                Gr = b.alloc_qubits(12)
                b.free(Gr)
            b.free(Y)
            m = b.alloc_bits(12)
            b.hmr(U[p], m)
            masks[p] = m
            b.free(U[p])
        if p0u is not None:                          # P0: U(0) out of place, measured like the others
            m = b.alloc_bits(12)
            b.hmr(p0u, m)
            masks[0] = m
            b.free(p0u)
        _MASKS.setdefault(k, []).append(masks)
    else:
        for p in reversed(pts):
            _pfd_tp(b, temp_point, p, U[p], False, p == pts[-1])
            b.free(U[p])


def _load_module():
    src = open(_SRC).read()
    pairs = [(_PAY_OLD, _PAY_NEW), (_END_OLD, _END_NEW), (_TP_OLD, _TP_NEW), (_START_OLD, _START_NEW)] + _LOOP_OLD
    for old, new in pairs:
        assert src.count(old) == 1, 'prefix5 source changed: deferral not applied'
        src = src.replace(old, new)
    mod = types.ModuleType('dqi_pub_pfd_prefix5')
    mod.__file__ = _SRC
    mod._pfd_pay = _pfd_pay
    mod._pfd_end = _pfd_end
    mod._pfd_tp = _pfd_tp
    mod._pfd_p0 = _pfd_p0
    mod._PFD_ZJOIN = ZJOIN
    exec(compile(src, _SRC + ':pfd', 'exec'), mod.__dict__)
    return mod


_MOD = _load_module()


def install(km, ns):
    """Replace the prefix block of every live external_extraction clone with the deferring copy."""
    walk = sys.modules['walk_bound']
    fns = walk.functions(ns.decode)
    need = ('ar', 'cx', 'scale', 'scalepair', 'combine', 'groups', 'ks', 'K', 'mul', 'inv', 'square12',
            'external_square24', 'mul24_xor99')
    mul_head = next(f for f in fns if f.__name__ == 'external_mul_head' and all(k in f.__globals__ for k in need))
    inner = _MOD.make_block(mul_head.__globals__)

    def block(b, words, homogeneous):
        _pfd_start(b)
        out = inner(b, words, homogeneous)
        if _DIR[id(b)]:
            if P0_STATE['recording']:
                P0_STATE['order'] = _p0_snapshot(b)
            elif P0_STATE['active']:
                _p0_restore(b)
        return out

    n = 0
    for f in fns:
        if f.__name__ == 'external_extraction' and '_dqb_block' in f.__code__.co_names:
            f.__globals__['_dqb_block'] = block
            n += 1
    assert n, 'no external_extraction with the prefix block found (is prefix5 installed?)'
    if P0:
        _p0_capture(km or sys.modules['circ'], ns.decode)
    return dict(extraction_clones=n, rinf=RINF, zjoin=ZJOIN, p0=P0,
                p0_recorded=None if P0_STATE['order'] is None else len(P0_STATE['order']))


def _p0_capture(km, decode):
    """Scratch decode of the line with the forward point 0 in place (as it is without P0), recording the free stack
    at the forward block's exit (as tools/dqi_pub_zt.py records its own two points)."""
    a3 = sys.modules.get('dqi_pub_a3')
    st = dict(a3.STACK) if a3 is not None else None   # a3's replay counters must not be spent on this scratch
    CB = km.CircuitBuilder
    tagged = CB.mark
    fv = getattr(getattr(tagged, '__code__', None), 'co_freevars', ())
    if 'done' in fv and 'original' in fv:                # r6tag's once-per-process tag: not on the scratch build
        CB.mark = dict(zip(fv, tagged.__closure__))['original'].cell_contents
    P0_STATE.update(active=False, recording=True)
    b0 = None
    try:
        b0 = CB()
        S0 = [b0.alloc_qubits(12) for _ in range(70)]
        E0 = [b0.alloc_bits(12) for _ in range(4095)]
        decode(b0, S0, E0)
    finally:
        CB.mark = tagged
        if a3 is not None:
            a3.STACK.update(st)
        P0_STATE.update(active=True, recording=False)
        if b0 is not None:
            for d in (_CALLS, _DIR, _MASKS):
                d.pop(id(b0), None)
    del b0, S0, E0
    import gc
    gc.collect()
    print(f'dqi_pub_pfd: P0 capture: forward block exit stack of {len(P0_STATE["order"] or ())} wires',
          file=sys.stderr, flush=True)
