#!/usr/bin/env python3
"""r62_inv_word_p2.py -- R62-INV-WORD (designer): THE WORD'S BIT 0 ON THE SLOT WIRE + THE ONE-WIRE LADDER (window W3).
In the parity-slot frame with the table word first (I39: slot "wfirst"), the word's lowest target wire, the SLOT, is
ZERO when the word is added. So (option h1c "r62": {"w0slot": 1}):
  forward  the table's output 0 is emitted INTO the slot wire (5 new table wires, not 6); position 0 of the pass is
           one AND by hand, c1 = w0 OR p (the framed target bit is 1), erased by measurement after the pass (its two
           factors are still there: no debt); the streamed pass (I39's "wstream", the held carry) runs on positions
           1.. with c1 as carry-in;
  undo     the slot already HOLDS w0 xor p: one CNOT gives w0 on its wire; the table is re-emitted for outputs 1..5
           only (same 10 gates, 5 outputs + 2 temporaries: FC's span search has a schedule there [M]); c1 again by
           one AND (free erase); the un-add on positions 1.. ; the slot's w0 is cleared with the table (bit 0's
           measured clear, as today).
So the word's table and un-add hold one wire less in both directions.
  "unadd1": 1  THE ONE-WIRE LADDER: the un-add as Cuccaro MAJ / UMA in the five word wires (carries in place, debts
           repaid there) and, above the word's last literal, ONE temporary A = c_(L+1): c_(L+2) = A t1 goes into its
           target by a CCX, its debt is a conditioned CZ, the write into c0's phantom top is ONE CCZ under m0
           (F3's XCCZ op: half a Toffoli a shot).  Without it the generator's own forms run on the shifted lists.
Needs I39's patches (slot .. slot11) and tools/r62_inv_word_p1.py.  Idempotent (marker).
  r62_inv_word_p2.py [SRC]      Safe to run: edits only MY snapshot."""
import os
import sys

SRC = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/qcc/ws/r62_inv_word/i39/src")
P = os.path.join(SRC, "invsqrt_fc2_patch_h1c.py")
assert "/r62_inv_word/" in P.replace("\\", "/"), "only my own source copy"
s = open(P).read()
MARK = "# R62-P2 (r62_inv_word_p2.py)"
if MARK in s:
    print("already patched:", P)
    sys.exit(0)


def rep(a, b, n=1):
    global s
    assert s.count(a) == n, ("anchor", s.count(a), a[:90])
    s = s.replace(a, b)


# ---- 0. the un-add (module level, before _sword_g0)
FUNC = '''def _r62_unadd(c, NW, src, tgt, cin, debts, m0, Cc, force):   ''' + MARK + '''
    """tgt += src + cin (the caller sets the X frame): word wires on the low positions, zero sources above, the top
    is c0's phantom. force (or no clean wire): Cuccaro MAJ / UMA in the word wires and ONE ladder wire."""
    n = len(tgt)
    src = list(src) + [None] * (n - len(src))
    debts = debts or {}
    L = max(i for i, a in enumerate(src) if NW._w(a)) + 1
    ok = (all(NW._w(a) for a in src[:L]) and all(a is None for a in src[L:]) and n - 1 - L == 3 and NW._w(cin)
          and m0 is not None)
    if not (force and ok):
        if NW.plan(max(0, n - 2 - NW.free_prefix(src, n, cin)), Cc)[0] is not None:
            NW.game_add(c, src, tgt, cin, debts, Cc)
            return "game"
        if Cc >= 2:                                    # F3's mixed form, as the stage uses for the word today
            import invsqrt_fw_c1add as CA
            if CA.NW is None:
                CA.NW = NW
            CA.mixed_add(c, src, tgt, cin, debts, Cc, L - 1, int(os.environ.get("C1CCZ", "1")))
            return "mixed"
        if Cc > NW.cuccaro_need(src, tgt, cin):
            NW.hybrid_add(c, src, tgt, cin, debts, Cc)
            return "hybrid"
        if not ok:
            NW.cuccaro_add(c, src, tgt, cin, debts)
            return "cuccaro"
    import invsqrt_fw_fcgen as _G
    _G._install_ops()
    chain = [cin] + src[:L]
    for i in range(L):                                 # MAJ: src[i] -> c_(i+1), debts repaid there
        x, y, z = chain[i], tgt[i], src[i]
        c.cx(z, y)
        c.cx(z, x)
        c.ccx(x, y, z)
        if (i + 1) in debts:
            c.z(z, debts[i + 1])
    k = src[L - 1]                                     # c_L
    t0, t1, t2 = tgt[L], tgt[L + 1], tgt[L + 2]
    A = c.alloc(1)[0]
    c.and_compute(k, t0, A)                            # c_(L+1) = c_L t0
    if (L + 1) in debts:
        c.z(A, debts[L + 1])
    if (L + 2) in debts:
        c.cz(A, t1, debts[L + 2])                      # c_(L+2) = A t1: its debt is a conditioned CZ
    c.xccz(A, t1, t2, m0)                              # the phantom top ^= c_(L+3) = A t1 t2: a CCZ under m0
    c.ccx(A, t1, t2)                                   # t2 ^= c_(L+2)
    c.cx(A, t1)                                        # t1 ^= c_(L+1)
    c.and_erase(k, t0, A)
    c.cx(k, t0)                                        # t0 ^= c_L
    for i in range(L - 1, -1, -1):                     # UMA
        x, y, z = chain[i], tgt[i], src[i]
        c.ccx(x, y, z)
        c.cx(z, x)
        c.cx(x, y)
    return "r62-ladder1"


def _sword_g0('''
rep("def _sword_g0(", FUNC)

# ---- 1. table(): an undo stage that leaves bit 0 alone
rep("    def table(bw_, st8=None):\n", "    def table(bw_, st8=None, skip0=False):                 " + MARK + "\n")
rep("        k_ = 0\n        done_ = []\n        for b_, given, pg in (st8 or stages):\n",
    "        k_ = 0\n        done_ = [0] if skip0 else []\n        for b_, given, pg in (st8 or stages):\n")

# ---- 2. forward: output 0 on the slot wire
rep('    bw = c.alloc(nbit)\n    _lend = list(F8.__dict__.get("_i39_lend") or [])',
    '    _r62w0 = (bool((o.get("r62") or {}).get("w0slot")) and bool(first) and bool(F8.__dict__.get("_i39_lend"))\n'
    '              and bool(stages) and s_out)             ' + MARK + '\n'
    '    if _r62w0:\n'
    '        assert Rt[0] == F8.__dict__["_i39_lend"][0], "w0slot: the slot is the pass\'s lowest target"\n'
    '        c.x(Rt[0])                                     # the framed slot (1) -> 0: the word\'s bit 0 lives here\n'
    '        bw = [Rt[0]] + c.alloc(nbit - 1)\n'
    '    else:\n'
    '        bw = c.alloc(nbit)\n'
    '    _lend = [w_ for w_ in list(F8.__dict__.get("_i39_lend") or []) if not (_r62w0 and w_ == Rt[0])]')

# ---- 3. forward: position 0 by hand, the streamed pass above it
A3 = ('    _cin(True)\n'
      '    if sw.get("hold"):                                # FC8: the word\'s carry held for the rows\' pass\n'
      '        import invsqrt_fc_s8m3 as _M3h\n'
      '        debts = _M3h.word_hold(c, NW, _fwd_stream, lits_, Rt, _cw, srcw, bid)\n')
R3 = ('    _cin(True)\n'
      '    if _r62w0:                                         ' + MARK + '\n'
      '        assert _cw is not None and sw.get("hold"), "w0slot: needs a carry-in wire and the held carry"\n'
      '        import invsqrt_fc_s8m3 as _M3h\n'
      '        _c1 = c.alloc(1)[0]\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.and_compute(Rt[0], _cw, _c1)                 # NOT (w0 OR p)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.x(_c1)                                       # c1 = w0 OR p: the framed target bit is 1\n'
      '        c.cx(_cw, Rt[0])\n'
      '        c.x(Rt[0])                                     # the framed sum bit 1 ^ w0 ^ p\n'
      '        debts = _M3h.word_hold(c, NW, _fwd_stream, lits_[1:], Rt[1:], _c1, srcw[1:], bid)\n'
      '        _M3h.HOLD["rpos"] += 1\n'
      '        _M3h.HOLD["bits"] = {k_ + 1: v_ for k_, v_ in _M3h.HOLD["bits"].items()}\n'
      '        c.x(Rt[0])\n'
      '        c.cx(_cw, Rt[0])                               # the slot is w0 again: c1 erased by measurement\n'
      '        c.x(_c1)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.and_erase(Rt[0], _cw, _c1)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.cx(_cw, Rt[0])\n'
      '        c.x(Rt[0])\n'
      '    elif sw.get("hold"):                              # FC8: the word\'s carry held for the rows\' pass\n'
      '        import invsqrt_fc_s8m3 as _M3h\n'
      '        debts = _M3h.word_hold(c, NW, _fwd_stream, lits_, Rt, _cw, srcw, bid)\n')
rep(A3, R3)
rep("    mbits = [_mb[w] if w in _mb else c.hmr(w) for w in bw]",
    "    mbits = [None if (_r62w0 and w == bw[0]) else (_mb[w] if w in _mb else c.hmr(w)) for w in bw]")
rep("    c.free_([w for w in bw if w not in _mb])\n",
    "    c.free_([w for w in bw if w not in _mb and not (_r62w0 and w == bw[0])])   " + MARK + "\n")

# ---- 4. undo: table for outputs 1..5, position 0 by hand
A4 = ('    if _mb:                                            # I39-SLOT5 (inv_i39_patch_slot5.py): the same wires back\n'
      '        _take(c, list(bw))\n'
      '    else:\n'
      '        bw2 = c.alloc(nbit)\n'
      '        assert list(bw2) == list(bw), ("table wires", bw2, bw)\n'
      '    k_u = table(bw, _ust8)\n'
      '    for w, mb in zip(bw, mbits):\n'
      '        c.z(w, mb)\n'
      '    Cc = c.cap - c.live\n'
      '    for w in Rt:\n'
      '        c.x(w)\n'
      '    _cin(True)\n'
      '    if NW.plan(')
R4 = ('    if _r62w0:                                         ' + MARK + '\n'
      '        if _mb:\n'
      '            _take(c, list(bw[1:]))\n'
      '        else:                                          # a build without the streamed release\n'
      '            bw2 = c.alloc(nbit - 1)\n'
      '            assert list(bw2) == list(bw[1:]), ("table wires", bw2, bw)\n'
      '        assert len(stages) == 1 and list(stages[0][0]) == list(range(nbit)) and not stages[0][1]\n'
      '        _pg6 = stages[0][2]\n'
      '        k_u = table(bw, [[list(range(1, nbit)), [], {"gates": _pg6["gates"], "outputs": _pg6["outputs"][1:]}]],\n'
      '                    skip0=True)\n'
      '        for w, mb in zip(bw[1:], mbits[1:]):\n'
      '            c.z(w, mb)\n'
      '        _cin(True)\n'
      '        c.x(Rt[0])\n'
      '        c.cx(_cw, Rt[0])                               # the framed slot 1 ^ w0 ^ p -> w0 on its wire\n'
      '        _c1 = c.alloc(1)[0]\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.and_compute(Rt[0], _cw, _c1)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.x(_c1)                                       # c1 = w0 OR p (erased clean in the forward: no debt)\n'
      '        Cc = c.cap - c.live\n'
      '        for w in Rt[1:]:\n'
      '            c.x(w)\n'
      '        how = _r62_unadd(c, NW, srcw[1:], Rt[1:], _c1, debts, (_M3x.PH or {}).get("m0"), Cc,\n'
      '                         bool((o.get("r62") or {}).get("unadd1")))\n'
      '        for w in Rt[1:]:\n'
      '            c.x(w)\n'
      '        c.x(_c1)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)\n'
      '        c.and_erase(Rt[0], _cw, _c1)\n'
      '        c.x(Rt[0])\n'
      '        c.x(_cw)                                       # the slot holds w0: bit 0\'s measured clear takes it\n'
      '        _cin(False)\n'
      '    if _r62w0:\n'
      '        pass\n'
      '    elif _mb:                                          # I39-SLOT5 (inv_i39_patch_slot5.py): the same wires back\n'
      '        _take(c, list(bw))\n'
      '    else:\n'
      '        bw2 = c.alloc(nbit)\n'
      '        assert list(bw2) == list(bw), ("table wires", bw2, bw)\n'
      '    if not _r62w0:\n'
      '        k_u = table(bw, _ust8)\n'
      '        for w, mb in zip(bw, mbits):\n'
      '            c.z(w, mb)\n'
      '        Cc = c.cap - c.live\n'
      '        for w in Rt:\n'
      '            c.x(w)\n'
      '        _cin(True)\n'
      '    if _r62w0:\n'
      '        pass\n'
      '    elif NW.plan(')
rep(A4, R4)
A5 = ('        NW.cuccaro_add(c, srcw, Rt, _cw, debts)\n'
      '        how = "cuccaro"\n'
      '    _cin(False)\n'
      '    for w in Rt:\n'
      '        c.x(w)\n')                                   # (ends here: I39's slot12 puts its correction next)
R5 = ('        NW.cuccaro_add(c, srcw, Rt, _cw, debts)\n'
      '        how = "cuccaro"\n'
      '    if not _r62w0:                                     ' + MARK + '\n'
      '        _cin(False)\n'
      '        for w in Rt:\n'
      '            c.x(w)\n')
rep(A5, R5)
A6 = ('            g_ = c.hmr(bw[b_])\n'
      '            _gone.add(b_)\n'
      '            if sw.get("clrfree"):\n')
R6 = ('            g_ = c.hmr(bw[b_])\n'
      '            _gone.add(b_)\n'
      '            if _r62w0 and b_ == 0:                     ' + MARK + '\n'
      '                c.x(bw[0])                             # the slot back to its framed zero\n'
      '                _fr.add(bw[0])\n'
      '            elif sw.get("clrfree"):\n')
rep(A6, R6)
open(P, "w").write(s)
print("patched:", P)
