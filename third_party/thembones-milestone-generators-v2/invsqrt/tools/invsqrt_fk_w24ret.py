"""F1 (fifth holder): digit 23 (word bit w24) retired across MAIN's rows 22..0 in a FROZEN FC g8 build (FC's src
folder is read-only: its st_main is patched IN MEMORY, cfg key "_w24ret": 1, next to FC's w25ret).

Forward: right after row 23 (its last reader) digit 23's wire is X-measured (bit m23) and freed for rows 22..0;
after the word's inversion w24 is rebuilt into a fresh wire by F1's exact 3-AND program (k = 2 UNSAT [P])
    g0 = (q21+q20+w23)(q19+q17); g1 = (q21+w23)(q19+w23); g2 = (q21+q20+q17+g0)(q19+q17+g0+g1); w24 = q20+q17+g0+g2
(q_i = wire i = v[i-14]; w23 = the word bit; + = XOR), temps erased by measurement, and the kick of m23 is repaid
by Z(w24) and NEG under m23. Undo (explicit blocks): at the undo's start w24 is X-measured (bit mP, 0 ANDs); before
row 23's undo digit 23 is re-created in its own wire from (v, digit 22) by the same program in the digit convention
(d23 = 1 + f(.., 1 + d22, ..)), on wires that are free there, and the kick of mP is repaid (Z, NEG under mP).
Price: 3 + 3 Toffolis per shot.
    python invsqrt_fk_w24ret.py emit SRC OUTDIR ANGLES.json CFG.json   build, kmx + kmb, sha, counts (no LP)
    python invsqrt_fk_w24ret.py run  SRC WD OUTDIR 'CFG json'           FC's run8 (all inputs + LP) on the patch
CFG: FC's cfg (fc2_src / ep_prog_dir are pointed into SRC here); add "_w24ret": 1 for the patch.
"""
import hashlib
import inspect
import json
import os
import sys

# F_b = F1's minimum programs (box 4 minand): digit b = F_b(q.., X, ..) with X = NOT digit b-1 (= the word bit w_b),
# on every input; measured between the word's inversion of w_b and of w_(b+1) (h1h43a ops 10,351 / 10,350; exact
# also on h1b42a). Digit 23: k = 3 (k = 2 UNSAT); digit 22: k = 4 (k = 3 UNSAT).
PROGS = {
    23: {"inputs": ["q21", "q19", "q20", "X", "q17"],
         "gates": [{"left": [0, 2, 3], "right": [1, 4]},
                   {"left": [0, 3], "right": [1, 3]},
                   {"left": [0, 2, 4, 5], "right": [1, 4, 5, 6]}],
         "outputs": [[2, 4, 5, 7]]},
    22: {"inputs": ["q20", "q21", "q19", "q18", "X", "q17"],
         "gates": [{"left": [0, 2, 3, 4], "right": [1, 2, 3, 4]},
                   {"left": [0], "right": [5, 6]},
                   {"left": [0, 1, 4, 6], "right": [3, 4, 5]},
                   {"left": [0, 3, 4, 5, 6, 7, 8], "right": [1, 4, 5, 7]}],
         "outputs": [[-1, 2, 6, 9]]},
}
PROG_F = PROGS[23]


def _toggle(lst):
    return [i for i in lst if i != -1] if -1 in lst else list(lst) + [-1]


def word_prog(b):
    """after the inversion (word bits): w_(b+1) = NOT d_b = 1 + F_b(q, w_b)."""
    P = PROGS[b]
    ins = [("w%d" % b) if t == "X" else t for t in P["inputs"]]
    return {"inputs": ins, "gates": P["gates"], "outputs": [_toggle(P["outputs"][0])]}


PROG_W24 = word_prog(23)


def digit_prog(dwire, b=23):
    """digit convention (right after row b): d_b = F_b(.., 1 + d_(b-1), ..), d_(b-1) on wire dwire."""
    P = PROGS[b]
    xi = P["inputs"].index("X")
    g = []
    for gt in P["gates"]:
        L, R = list(gt["left"]), list(gt["right"])
        if xi in L:
            L = _toggle(L)
        if xi in R:
            R = _toggle(R)
        g.append({"left": L, "right": R})
    ins = [("q%d" % dwire) if t == "X" else t for t in P["inputs"]]
    return {"inputs": ins, "gates": g, "outputs": [list(P["outputs"][0])]}


class W24:
    def __init__(self, F8, K):
        self.F8, self.K = F8, K
        self.log = []

    def retire(self, c, wire, d22, v, b=23):
        assert [int(q) for q in v] == list(range(14, 22)), v
        if _STATE["uret"]:                             # F1 sixth holder: undo-only form (cfg "_uret")
            prog = digit_prog(d22, b)
            alt = (_STATE.get("cfg") or {}).get("_uprog", {}).get(str(b))   # word-convention file for the re-creation
            if alt:
                cfg_ = _STATE["cfg"]
                p = os.path.expanduser(str(alt))
                if not os.path.isabs(p):
                    p = os.path.join(os.path.expanduser(cfg_.get("fc2_src") or "."), p)
                prog = self.F8.__dict__["_M7Retire"].digit_form(json.load(open(p)), b, d22)
            r = _uret_retire(self.F8, self.K, c, wire, prog, ("m8", b))
            self.log.append(dict(uret_digit=b, wire=int(wire)))
            return r
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        m23 = c.hmr(wire)
        c.free_([wire])
        c.ops.append(("BEND", bid))
        mP = c.bit()
        cs = self.K.Circ(n_fixed=c.next)               # scratch circuit: the re-creation on wires free right here
        cs.free = [w for w in c.free if w != wire]
        rb = self.F8._make_rebuild(digit_prog(d22, b), [])
        rb(cs, wire)
        assert cs.next == c.next, "re-creation needs more free wires than there are after row %d" % b
        c.blocks[bid] = list(cs.ops) + [("Z", wire, mP), ("NEG", mP)]
        self.log.append(dict(digit=b, drop_wire=int(wire), below=int(d22), free_after=len(c.free),
                             recreate_ops=len(cs.ops)))
        return (wire, m23, mP)

    def retire22(self, c, wire, dlow, v, cfg):
        """F1 sixth holder: digit 21 (word bit w22) lent in the UNDO only, any m. cfg "_w22ret": a word-convention
        program file (w22 from q.., w21; relative paths in the src folder), converted to the digit convention
        (d21 from q.., d20 on wire dlow) by FC's _M7Retire.digit_form; re-created before row 21's undo."""
        p = os.path.expanduser(str(cfg["_w22ret"]))
        if not os.path.isabs(p):
            p = os.path.join(os.path.expanduser(cfg.get("fc2_src") or "."), p)
        prog = self.F8.__dict__["_M7Retire"].digit_form(json.load(open(p)), 21, dlow)
        self.log.append(dict(w22_lend=int(wire), dlow=int(dlow), ands=len(prog["gates"])))
        return _uret_retire(self.F8, self.K, c, wire, prog, ("w22", 21), frugal=True)

    def retire21(self, c, wire, dlow, v, cfg):
        """F1 sixth holder: digit 20 (word bit w21) lent in the UNDO only. cfg "_w21ret": a DIGIT-convention program
        file (d20 from q15.. and the lower digit, input named "d"), measured at the lend point; re-created before
        row 20's undo with the frugal emitter (5 wires free there when digits 21..23 are lent too)."""
        p = os.path.expanduser(str(cfg["_w21ret"]))
        if not os.path.isabs(p):
            p = os.path.join(os.path.expanduser(cfg.get("fc2_src") or "."), p)
        P = json.load(open(p))
        prog = dict(P)
        prog["inputs"] = [("q%d" % dlow) if t == "d" else t for t in P["inputs"]]
        self.log.append(dict(w21_lend=int(wire), dlow=int(dlow), ands=len(prog["gates"])))
        return _uret_retire(self.F8, self.K, c, wire, prog, ("w21", 20), frugal=True)

    def rebuild22(self, c, ret):
        return _uret_rebuild(c, ret)

    def rebuild(self, c, word, v, ret, j=24):
        if len(ret) == 4 and ret[3] == "uret":        # F1 sixth holder: undo-only form
            return _uret_rebuild(c, ret)
        wire0, m23, mP = ret
        t = c.alloc(1)[0]
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        self.F8._make_rebuild(word_prog(j - 1), word)(c, t)
        c.z(t, m23)
        c.neg(m23)
        c.ops.append(("BEND", bid))
        c.blocks[bid] = [("HMR", t, mP)]
        self.log.append(dict(rebuild_wire=int(t), live=c.live))
        return t


A_OLD = "            w25 = _retire_digit24(c, digit[24], pos(34), v)\n"
A_NEW = A_OLD + ("        if b == 23 and cfg.get(\"_w24ret\") and m == 8:   # F1: digit 23 retired after row 23\n"
                 "            _F1W24_ret = _F1W24.retire(c, digit[23], pos(33), v)\n"
                 "            digit[23] = None\n"
                 "        if b == 22 and cfg.get(\"_w23ret\") and m == 8:   # F1: digit 22 retired after row 22\n"
                 "            _F1W23_ret = _F1W24.retire(c, digit[22], pos(32), v, 22)\n"
                 "            digit[22] = None\n"
                 "        if b == 21 and cfg.get(\"_w22ret\"):   # F1 (sixth): digit 21 lent in the undo only (any m)\n"
                 "            _F1W22_ret = _F1W24.retire22(c, digit[21], pos(23 + m), v, cfg)\n"
                 "            digit[21] = None\n"
                 "        if b == 20 and cfg.get(\"_w21ret\"):   # F1 (sixth): digit 20 lent in the undo only (any m)\n"
                 "            _F1W21_ret = _F1W24.retire21(c, digit[20], pos(22 + m), v, cfg)\n"
                 "            digit[20] = None\n")
B_OLD = "            if j != 25:\n                c.x(w)\n"
B_NEW = "            if j != 25 and w is not None:\n                c.x(w)\n"
C_OLD = "        word[25] = _rebuild_w25(c, word[24], v, w25)\n"
C_NEW = ("        if cfg.get(\"_w21ret\"):\n"
         "            word[21] = _F1W24.rebuild22(c, _F1W21_ret)\n"
         "        if cfg.get(\"_w22ret\"):\n"
         "            word[22] = _F1W24.rebuild22(c, _F1W22_ret)\n"
         "        if cfg.get(\"_w23ret\") and m == 8:\n"
         "            word[23] = _F1W24.rebuild(c, word, v, _F1W23_ret, 23)\n"
         "        if cfg.get(\"_w24ret\") and m == 8:\n"
         "            word[24] = _F1W24.rebuild(c, word, v, _F1W24_ret)\n") + C_OLD


def load(src):
    """FC's frozen modules first on the path and cached (so no later sys.path change reaches FC's live tools/)."""
    src = os.path.abspath(os.path.expanduser(src))
    sys.path.insert(0, src)
    import invsqrt_fc_kit as K  # noqa: F401
    import invsqrt_fc_g8 as F8
    for f in sorted(os.listdir(src)):
        if f.endswith(".py") and "-" not in f and f.startswith("invsqrt_fc") and "run" not in f:
            __import__(f[:-3])
    return F8, sys.modules["invsqrt_fc_kit"]


def apply(F8, K):
    s = inspect.getsource(F8.st_main)
    for old, new in ((A_OLD, A_NEW), (B_OLD, B_NEW), (C_OLD, C_NEW)):
        assert s.count(old) == 1, old
        s = s.replace(old, new)
    F8.__dict__["_F1W24"] = W24(F8, K)
    exec(compile(s, F8.__file__ + " [F1 w24ret patch]", "exec"), F8.__dict__)
    _install_uret(F8, K)
    return F8.__dict__["_F1W24"]


# ---- F1 sixth holder (2026-11-14 21:45Z): UNDO-ONLY digit retires, cfg "_uret": 1 (off: builds byte-identical) ----
# Measured on FC's builds: NW.block's forward uses at most 2 carry wires at any C, so MAIN's forward rows cost the same
# at C 2 and C 5; only the explicit undo blocks (built at forward time with C = cap - live) use a retired digit's wire.
# So the forward rebuild of a retired digit buys nothing once the rows below run at C >= 2. Undo-only form:
#   forward: the digit stays on its wire (an EMPTY block after its row; the word's inversion X's it at MAIN's end);
#   undo:    right after the inversion's mirror the digit is X-measured (bit mP, 0 ANDs); its wire is LENT to every
#            explicit undo built in the forward between its row and MAIN's end (they run while it is free); before
#            its own row's undo it is re-created by the same program (the existing undo side), then Z under mP
#            (no NEG: it is measured in the digit convention). Price k ANDs per shot instead of 2k.
_STATE = {"uret": False}
ULOG = []


def _rebuild_frugal(prog):
    """re-creation emitter with fewer live temps than FC's _make_rebuild (same AND count): gates read by no later gate
    and whose factors read inputs only ("free" gates) are computed into a temp, XORed into the target and erased at
    once; the other gates are computed in order, the last of them straight into the target, erased in reverse."""
    ins = [int(t[1:]) for t in prog["inputs"]]
    for t in prog["inputs"]:
        assert t.startswith("q"), "frugal rebuild: inputs must be wires (q<i>)"
    nin, gates = len(ins), prog["gates"]
    (out,) = prog["outputs"]
    k = len(gates)
    used_later = set()
    for g in gates:
        for j in g["left"] + g["right"]:
            if j >= nin:
                used_later.add(j)
    free = [s for s in range(k) if (nin + s) not in used_later
            and all(j < nin for j in gates[s]["left"] + gates[s]["right"])]
    nested = [s for s in range(k) if s not in free]

    def rebuild(c, t):
        val = {i: ins[i] for i in range(nin)}

        def side(items, other):
            ws = [val[i] for i in items if i != -1]
            one = -1 in items
            excl = [w for w in ws if w not in {val[j] for j in other if j != -1}]
            if excl:
                acc, tmp = excl[0], None
            else:
                acc = tmp = c.alloc(1)[0]
            done = []
            for w in ws:
                if w != acc:
                    c.cx(w, acc)
                    done.append(("cx", w, acc))
            if one:
                c.x(acc)
                done.append(("x", acc))
            if tmp is not None:
                done.append(("tmp", tmp))
            return acc, done

        def undo(done):
            for d in reversed(done):
                if d[0] == "cx":
                    c.cx(d[1], d[2])
                elif d[0] == "x":
                    c.x(d[1])
                else:
                    c.free_([d[1]])

        def gate(s, into=None):
            gt = gates[s]
            a, da = side(gt["left"], gt["right"])
            b, db = side(gt["right"], gt["left"])
            q = c.and_compute(a, b, into) if into is not None else c.and_compute(a, b)
            undo(db)
            undo(da)
            return q

        def erase(s, q):
            gt = gates[s]
            a, da = side(gt["left"], gt["right"])
            b, db = side(gt["right"], gt["left"])
            c.and_erase(a, b, q)
            undo(db)
            undo(da)
        direct = nested[-1] if nested and (nin + nested[-1]) in out else None
        temps = []
        for s in nested:
            if s == direct:
                val[nin + s] = gate(s, t)
            else:
                val[nin + s] = gate(s)
                temps.append(s)
        for i in out:
            if i == -1:
                c.x(t)
            elif i < nin or (i - nin) in nested and (i - nin) != direct:
                c.cx(val[i], t)
        for s in reversed(temps):                  # nested temps erased first: the free gates then need one temp
            erase(s, val[nin + s])
        for s in free:
            if (nin + s) not in out:               # (a free gate outside the output does nothing)
                continue
            q = gate(s)
            c.cx(q, t)
            erase(s, q)
    return rebuild


def _uret_retire(F8, K, c, wire, prog, tag, frugal=False):
    """empty forward block; explicit undo = re-create the digit on `wire` from `prog` on the wires free there."""
    lend = list(getattr(c, "_fk_lend", None) or [])
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    c.ops.append(("BEND", bid))
    mP = c.bit()
    cs = K.Circ(n_fixed=c.next)                    # scratch: free wires at the re-creation = free now + lent ones
    cs.free = list(c.free) + [w for w in lend if w != wire]
    frugal = frugal or bool((_STATE.get("cfg") or {}).get("_frugal_all"))   # (debug: every re-creation frugal)
    (_rebuild_frugal(prog) if frugal else F8._make_rebuild(prog, []))(cs, wire)
    assert cs.next == c.next, "uret re-creation of %s needs more free wires than there are" % (tag,)
    c.blocks[bid] = list(cs.ops) + [("Z", wire, mP)]
    c._fk_lend = lend + [int(wire)]
    ULOG.append(dict(retire=tag, wire=int(wire), recreate_ops=len(cs.ops), lend=list(c._fk_lend)))
    return (wire, None, mP, "uret")


def _uret_rebuild(c, ret):
    """at MAIN's end (after the inversion loop, which skipped the None word bit): an empty block whose undo
    X-measures the digit, then the X of the inversion; the word bit stays on the digit's own wire."""
    wire, _, mP, _ = ret
    c._fk_lend = [w for w in (getattr(c, "_fk_lend", None) or []) if w != wire]
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    c.ops.append(("BEND", bid))
    c.blocks[bid] = [("HMR", wire, mP)]
    c.x(wire)
    ULOG.append(dict(rebuild=int(wire), live=c.live))
    return wire


def _install_uret(F8, K):
    """lending NW.block, the m 7 _M7Retire hooks and the cfg switch (build_full wrapper); idempotent."""
    if F8.__dict__.get("_fk_uret_installed"):
        return
    F8.__dict__["_fk_uret_installed"] = True
    NW = sys.modules["invsqrt_fc_narrow"]
    orig_block = NW.block

    def block(c, src, tgt, cin):
        lend = list(getattr(c, "_fk_lend", None) or [])
        if not lend:
            return orig_block(c, src, tgt, cin)
        busy = set(w for w in list(src) + list(tgt) + [cin] if isinstance(w, int))
        assert not busy.intersection(lend), ("lent wire read by the block", lend)
        if not hasattr(c, "blocks"):
            c.blocks = {}
        bid = len(c.blocks)
        c.blocks[bid] = None
        c.ops.append(("BBEG", bid))
        debts = NW.fwd(c, src, tgt, cin)
        c.ops.append(("BEND", bid))
        for w in lend:                             # free (X-measured) while this explicit undo runs
            c.free.append(w)
            c.live -= 1
        C = c.cap - c.live
        save = c.ops
        c.ops = []
        for w in tgt:
            c.x(w)
        n = len(tgt)
        if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, cin)), C)[0] is not None:
            cost = NW.game_add(c, src, tgt, cin, debts, C)
        elif C > NW.cuccaro_need(src, tgt, cin):
            cost = float(NW.hybrid_add(c, src, tgt, cin, debts, C))
        else:
            cost = float(NW.cuccaro_add(c, src, tgt, cin, debts))
        for w in tgt:
            c.x(w)
        c.blocks[bid] = c.ops
        c.ops = save
        for w in lend:
            c.free.remove(w)
            c.live += 1
        if not hasattr(c, "block_log"):
            c.block_log = []
        c.block_log.append((len(tgt), C, cost))
    NW.block = block

    def patch_m7():
        """FC's _M7Retire (re-created whenever FC2's w26 patch re-execs st_main's source, which carries the class):
        patched at build time, after the fc2 patches."""
        M7 = F8.__dict__.get("_M7Retire")
        if M7 is None or getattr(M7, "_fk_uret_patched", False):
            return
        o_ret, o_reb = M7.retire, M7.rebuild

        def retire(self, c, wire, dlow, v, b, cfg):
            if not _STATE["uret"]:
                return o_ret(self, c, wire, dlow, v, b, cfg)
            return _uret_retire(F8, K, c, wire, self.digit_form(self.prog(b, cfg), b, dlow), ("m7", b))

        def rebuild(self, c, word, v, ret, b, cfg):
            if len(ret) == 4 and ret[3] == "uret":
                return _uret_rebuild(c, ret)
            return o_reb(self, c, word, v, ret, b, cfg)
        M7.retire, M7.rebuild = retire, rebuild
        M7._fk_uret_patched = True

    obuild = F8.build

    def build(*a, **kw):
        patch_m7()
        return obuild(*a, **kw)
    F8.__dict__["build"] = build

    ob = F8.build_full

    def build_full(cfg, ang=None, feat=None):
        _STATE["uret"] = bool(cfg.get("_uret"))
        _STATE["ladder"] = bool(cfg.get("_ladder"))
        _STATE["cfg"] = cfg
        _install_ladder_hook()                     # outermost A.add wrapper (after every other patch)
        c, regs = ob(cfg, ang, feat)
        assert not getattr(c, "_fk_lend", None), "uret: a lent wire was never returned"
        assert not _STATE["ladder"] or getattr(c, "_fk_ladder_done", False), "ladder: MAIN's first add not found"
        return c, regs
    F8.__dict__["build_full"] = build_full


# ---- F1 sixth holder (2026-11-14 22:05Z): MAIN's first add as an AND ladder, cfg "_ladder": 1 (off: byte-identical) ----
# Before the first add the window is (a, s, s, ..., s) (a = U's bit at position 26, s = its sign at 27, the E wires
# copies of s) and the add is R = win - D with D = 1 + 2v + 2^(m+1), so R = (2^(m+1) - 1) - 2v + a - 2s:
#   R bit 0 = NOT a;  R >> 1 = ~v + (a - s)  (an (m+1)-bit increment / decrement of ~v, or nothing).
# Bit j of R >> 1 = NOT v_j ^ t_j with t_0 = p = a ^ s and t_(j+1) = t_j AND y_j, y_j = v_j ^ NOT s (both directions in
# one ladder). The t's are computed INTO the E wires (cleared first: they held copies of s): m ANDs, no clean wire
# needed, and the mirror erases them by measurement (ANDC -> ANDE: 0 Toffolis). Same window values on every input.
def _ladder_first_add(c, src, tgt):
    m = len(tgt) - 2
    v = list(src[1:m + 1])
    a, s, E = tgt[0], tgt[1], list(tgt[2:])
    for t in tgt:                                  # undo the caller's complement frame: tgt = (a, s, s, ..., s)
        c.x(t)
    for e in E:                                    # the copies of s -> 0
        c.cx(s, e)
    c.cx(a, s)                                     # s-wire := p = a ^ s (t_0)
    prev = s
    lazy = bool((_STATE.get("cfg") or {}).get("_r25lazy"))
    for i in range(m - 1 if lazy else m):          # _r25lazy: t_m (the top, rare) is left to row 25's top step
        c.cx(s, v[i])                              # v_i -> y_i = v_i ^ p ^ a ^ 1 = v_i ^ NOT s
        c.cx(a, v[i])
        c.x(v[i])
        c.and_compute(prev, v[i], E[i])            # t_(i+1) into E[i] (position 28 + i)
        c.x(v[i])
        c.cx(a, v[i])
        c.cx(s, v[i])
        prev = E[i]
    if lazy:
        c._fk_lazy = dict(top=E[m - 1], v=list(v), m=m)
    for i in range(m - 1):                         # bit i+1 of R >> 1 = NOT v_(i+1) ^ t_(i+1); the top bit is t_m
        c.cx(v[i + 1], E[i])
        c.x(E[i])
    c.cx(v[0], s)                                  # bit 0 of R >> 1 = NOT v_0 ^ p
    c.x(s)
    c.x(a)                                         # R bit 0 = NOT a
    for t in tgt:                                  # the caller's complement frame again
        c.x(t)


def _r25_stub_block(c, src, tgt):
    """COUNT-ONLY (cfg "_fk_r25stub"): MAIN row 25 as a narrow block whose forward runs with one more clean wire
    (as if its top target were born at the top step) and whose explicit undo is built at the real cap."""
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    c.cap += 1
    debts = NW.fwd(c, src, tgt, None)
    c.cap -= 1
    c.ops.append(("BEND", bid))
    C = c.cap - c.live
    save = c.ops
    c.ops = []
    for w in tgt:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:
        cost = NW.game_add(c, src, tgt, None, debts, C)
    elif C > NW.cuccaro_need(src, tgt, None):
        cost = float(NW.hybrid_add(c, src, tgt, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, src, tgt, None, debts))
    for w in tgt:
        c.x(w)
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))


# ---- F1 sixth holder (2026-11-14 23:10Z): MAIN row 25 with a LAZY TOP, cfg "_r25lazy": 1 (with "_ladder": 1) ----
# After the ladder, row 25's top target (position 26 + m) holds t_m = t_(m-1) AND (v_(m-1) ^ v_0 ^ 1) [t_(m-1) = 1 forces
# v_0 = s], a 1-AND function of wires intact until row 25's top step (t_(m-1) = R bit m ^ NOT v_(m-1), on target m + 1,
# whose sum bit is written after the top). So the ladder skips t_m and row 25 runs as a narrow block whose top wire is
# born at the top step: its forward uses two carry wires (one-clean-wire form before), its explicit undo is built at
# the real cap (the top is live in the undo), then t_m is erased by measurement. Ladder -1 AND, row 25 +1 AND.
def _load_t_m(c, tgt, top, lz):
    """row 25's top (complemented frame): NOT t_m, t_m = t_(m-1) AND (v_(m-1) ^ v_0 ^ 1); 1 AND."""
    m, v = lz["m"], lz["v"]
    tb = tgt[m + 1]                                # holds v_(m-1) ^ t_(m-1) here (complemented frame)
    c.cx(v[m - 1], tb)
    c.cx(v[0], v[m - 1])
    c.x(v[m - 1])
    c.and_compute(tb, v[m - 1], top)               # t_m
    c.x(v[m - 1])
    c.cx(v[0], v[m - 1])
    c.cx(v[m - 1], tb)
    c.x(top)


def _r24_form(c, tgt, d, v, top, erase):
    """row 24's top in its framed coordinates (inside fadd's frames): top' = q19 ^ q20 ^ T' ^ d ^ 1 ^ A B with
    A = d ^ q19, B = q21 ^ T' ^ d ^ 1 (T' = tgt[-2] framed, d = the row's control, q19..q21 = v[-3], v[-2], v[-1];
    from R25 bit 8 = 1 ^ q19 ^ q20 ^ d ^ T ^ (d ^ q19)(q21 ^ T) before the frame [M, all inputs, m 7 and m 8]).
    erase=False: load top' into the |0> wire `top` (1 AND); erase=True: clear it by measurement (0 Toffolis)."""
    q19, q20, q21 = v[-3], v[-2], v[-1]
    T = tgt[-2]
    if erase:                                      # top = L' ^ A B -> A B
        for w in (q19, q20, T, d):
            c.cx(w, top)
        c.x(top)
    c.cx(d, T)                                     # T -> T' ^ d
    c.cx(q21, T)
    c.x(T)                                         # T -> B = q21 ^ T' ^ d ^ 1
    c.cx(q19, d)                                   # d -> A = d ^ q19
    if erase:
        c.ops.append(("ANDE", d, T, top))
    else:
        c.and_compute(d, T, top)
    c.cx(q19, d)
    c.x(T)
    c.cx(q21, T)
    c.cx(d, T)
    if not erase:                                  # top = A B ^ L'
        for w in (q19, q20, T, d):
            c.cx(w, top)
        c.x(top)


def _fwd_lazy(c, NW, src, tgt, top, lz, loader=None):
    """copy of NW.fwd (tgt += src, at most 2 carry wires) where tgt[-1] = `top` is |0> and NOT live on entry: at the
    top step it is taken (moving the held carry off it if needed) and loaded by `loader` (default: row 25's NOT t_m)."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    debts = {}
    cur = NW._lit(None)
    cur_fresh = False
    deferred = None
    for i in range(n - 1):
        a, bw = src[i], tgt[i]
        if i == n - 2:                             # the top step: take the top wire and load its value
            if cur_fresh and cur[1] == top:
                q = c.alloc(1)[0]
                c.cx(top, q)
                c.cx(q, top)
                cur = ("w", q, cur[2])
            else:
                c.free.remove(top)
                c.live += 1
                c.peak = max(c.peak, c.live)
            (loader or _load_t_m)(c, tgt, top, lz)
        f = NW._form_b(a, bw, cur)
        nxt, nfresh = None, False
        if i < n - 2:
            if f[0] == "const":
                nxt = ("k", f[1])
            elif f[0] == "alias":
                nxt = ("w", f[1], f[2])
            else:
                _, xw, yw, ops = f
                NW._emit(c, ops)
                t = c.alloc(1)[0]
                c.and_compute(xw, yw, t)
                NW._emit(c, ops[::-1])
                NW._xor_lit(c, cur, t)
                nxt, nfresh = ("w", t, 0), True
        else:
            if f[0] == "const":
                if f[1]:
                    c.x(top)
            elif f[0] == "alias":
                c.cx(f[1], top)
                if f[2]:
                    c.x(top)
            else:
                _, xw, yw, ops = f
                NW._emit(c, ops)
                c.ccx(xw, yw, top)
                NW._emit(c, ops[::-1])
                NW._xor_lit(c, cur, top)
            NW._xor_src(c, src[n - 1], top)
        if f[0] == "alias" and i < n - 2:
            new_def = i
        else:
            new_def = None
            NW._xor_src(c, a, bw)
            NW._xor_lit(c, cur, bw)
        if cur_fresh:
            debts[i] = c.hmr(cur[1])
            c.free_([cur[1]])
        elif deferred is not None:
            c.x(tgt[deferred])
        deferred = new_def
        cur, cur_fresh = nxt, nfresh
    assert deferred is None
    return debts


def _r25_lazy_block(c, src, tgt, lz):
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    top, m, v = lz["top"], lz["m"], lz["v"]
    c.x(top)                                       # st_main's complement frame made the placeholder 1: back to |0>
    c.free.insert(0, top)                          # free, and taken last by the carries
    c.live -= 1
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    debts = _fwd_lazy(c, NW, src, tgt, top, lz)
    c.ops.append(("BEND", bid))
    assert top not in c.free, "lazy top: the top wire was not taken"
    C = c.cap - c.live + (1 if (_STATE.get("cfg") or {}).get("_fk_undo41") else 0)   # COUNT-ONLY stub knob
    save = c.ops
    c.ops = []
    for w in tgt:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:
        cost = NW.game_add(c, src, tgt, None, debts, C)
    elif C > NW.cuccaro_need(src, tgt, None):
        cost = float(NW.hybrid_add(c, src, tgt, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, src, tgt, None, debts))
    for w in tgt:
        c.x(w)
    tb = tgt[m + 1]                                # erase t_m (the top's pre-add value is NOT t_m again)
    c.x(top)
    c.cx(v[m - 1], tb)
    c.cx(v[0], v[m - 1])
    c.x(v[m - 1])
    c.ops.append(("ANDE", tb, v[m - 1], top))
    c.x(v[m - 1])
    c.cx(v[0], v[m - 1])
    c.cx(v[m - 1], tb)
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    ULOG.append(dict(r25lazy=int(top), C_undo=C))
    if (_STATE.get("cfg") or {}).get("_r24lazy"):
        c._fk_r24 = dict(ctl=top, v=list(v), m=m)  # row 25's top is row 24's control


def _r24_lazy_block(c, src, tgt, info):
    """MAIN row 24 (inside fadd's frames) as a narrow block with a lazy top: its top target (R25 bit 8, a 1-AND
    function of the control, the next target and q19..q21) is erased by measurement at the block's start and
    re-created at the top step (+1 AND), so the forward uses two carry wires; explicit undo = the plain undo add."""
    NW = sys.modules["invsqrt_fc_narrow"]
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    top, d, v = tgt[-1], info["ctl"], info["v"]
    if not hasattr(c, "blocks"):
        c.blocks = {}
    bid = len(c.blocks)
    c.blocks[bid] = None
    c.ops.append(("BBEG", bid))
    _r24_form(c, tgt, d, v, top, erase=True)
    c.free.insert(0, top)
    c.live -= 1

    def loader(c_, tg, tp, lz):
        _r24_form(c_, tg, d, v, tp, erase=False)
    debts = _fwd_lazy(c, NW, src, tgt, top, info, loader=loader)
    c.ops.append(("BEND", bid))
    assert top not in c.free, "lazy top (row 24): the top wire was not taken"
    C = c.cap - c.live + (1 if (_STATE.get("cfg") or {}).get("_fk_undo41") else 0)   # COUNT-ONLY stub knob
    save = c.ops
    c.ops = []
    for w in tgt:
        c.x(w)
    if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, None)), C)[0] is not None:
        cost = NW.game_add(c, src, tgt, None, debts, C)
    elif C > NW.cuccaro_need(src, tgt, None):
        cost = float(NW.hybrid_add(c, src, tgt, None, debts, C))
    else:
        cost = float(NW.cuccaro_add(c, src, tgt, None, debts))
    for w in tgt:
        c.x(w)
    c.blocks[bid] = c.ops
    c.ops = save
    if not hasattr(c, "block_log"):
        c.block_log = []
    c.block_log.append((n, C, cost))
    ULOG.append(dict(r24lazy=int(top), C_undo=C))


def _install_ladder_hook():
    A = sys.modules["invsqrt_fc_arith"]
    if getattr(A.add, "_fk_ladder", False):
        return
    orig = A.add
    one = getattr(A, "ONE", "1")

    def add(c, src, tgt, cin=None, force=None):
        if (_STATE["ladder"] and not getattr(c, "_fk_ladder_done", False) and getattr(c, "tags", None)
                and c.tags[-1][1] == "MAIN" and cin is None and force is None):
            s_ = list(src)
            if (len(s_) == len(tgt) and len(tgt) >= 4 and s_[0] == one and s_[-1] == one
                    and all(isinstance(w, int) for w in s_[1:-1]) and all(isinstance(w, int) for w in tgt)):
                c._fk_ladder_done = True
                c._fk_r25_next = True
                return _ladder_first_add(c, s_, list(tgt))
        if ((_STATE.get("cfg") or {}).get("_fk_r25stub") and getattr(c, "_fk_r25_next", False)
                and c.tags[-1][1] == "MAIN" and cin is None and force is None):
            c._fk_r25_next = False                 # COUNT-ONLY stub: row 25's forward at cap + 1, its undo at cap
            return _r25_stub_block(c, list(src), list(tgt))
        if (getattr(c, "_fk_lazy", None) and getattr(c, "_fk_r25_next", False)
                and c.tags[-1][1] == "MAIN" and cin is None and force is None):
            c._fk_r25_next = False
            lz = c._fk_lazy
            c._fk_lazy = None
            assert list(tgt)[-1] == lz["top"] and len(tgt) == lz["m"] + 3, "lazy top: row 25 not where expected"
            return _r25_lazy_block(c, list(src), list(tgt), lz)
        if (getattr(c, "_fk_r24", None) and c.tags[-1][1] == "MAIN" and cin is None and force is None):
            info = c._fk_r24
            c._fk_r24 = None
            assert len(tgt) == info["m"] + 3 and info["ctl"] not in list(tgt), "lazy top: row 24 not where expected"
            return _r24_lazy_block(c, list(src), list(tgt), info)
        return orig(c, src, tgt, cin, force)
    add._fk_ladder = True
    A.add = add


def cfg_for(F8, src, cfg_in):
    cfg = dict(F8.CFG)
    cfg.update(cfg_in)
    cfg["fc2_src"] = src
    cfg["ep_prog_dir"] = os.path.join(src, "fcsup")
    if cfg.get("_eplane"):                         # F1: the endpoint lends the constant token's wire as a lane
        import invsqrt_fk_epN as EPN               # (the src folder must hold F1's epN with const_lane)
        if not getattr(EPN.emit_ring, "_f1_eplane", False):
            orig = EPN.emit_ring

            slots = bool(cfg.get("_retslots"))
            w25l = bool(cfg.get("_w25lane"))

            def emit_ring(*a, **k):
                k.setdefault("const_lane", True)
                if slots:
                    k.setdefault("retire_slots", True)
                if w25l:
                    k.setdefault("late_w25", (20, 21))      # q20, q21 wires (v[6], v[7])
                return orig(*a, **k)
            emit_ring._f1_eplane = True
            EPN.emit_ring = emit_ring
    return cfg


def main():
    mode = sys.argv[1]
    src = os.path.abspath(os.path.expanduser(sys.argv[2]))
    F8, K = load(src)
    w = apply(F8, K)
    if mode == "emit":
        out, angp, cfgp = sys.argv[3], sys.argv[4], sys.argv[5]
        a = json.load(open(angp))
        cfg = cfg_for(F8, src, json.load(open(cfgp)))
        K.KMX_CBLOCKS = bool(cfg.get("kmx_cblocks"))
        c, regs = F8.build_full(cfg, [int(t) for t in a["ang"]], [int(t) for t in a["feat"]])
        os.makedirs(out, exist_ok=True)
        p = os.path.join(out, "invsqrt.kmx")
        with open(p, "w") as f:
            f.write(K.kmx_text(c, regs["X"]))
        sys.path.insert(0, os.path.join(os.path.expanduser("~/qcc/data"), "kickmix"))
        import kickmix as km
        kc = km.Circuit.from_file(p, format="kmx")
        kb = os.path.join(out, "invsqrt.kmb")
        kc.to_file(kb, format="kmb")
        sha = hashlib.sha256(open(kb, "rb").read()).hexdigest()
        print("emit %s: sha %s toffolis_exp %.1f static %d wires %d peak %d stage peaks %s; patch log %s" % (
            out, sha[:20], c.toffolis_exp(), c.toffolis(), c.next, c.peak, c.stage_peaks, w.log))
    elif mode == "run":
        wd, out, cfgs = sys.argv[3], sys.argv[4], sys.argv[5]
        cfg = cfg_for(F8, src, json.loads(cfgs))
        tools = os.path.dirname(os.path.abspath(__file__))
        sys.path.append(tools)
        import invsqrt_fc_run8 as R8
        sys.argv = ["invsqrt_fc_run8.py", wd, out, json.dumps(cfg)]
        R8.main()
        print("patch log", w.log, flush=True)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
