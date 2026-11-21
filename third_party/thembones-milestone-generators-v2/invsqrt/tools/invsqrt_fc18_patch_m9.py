#!/usr/bin/env python3
"""invsqrt_fc18_patch_m9.py -- FC eighteenth holder: FC18-M9, digit 23's erase phase in row 24's undo paid at once.
Row 24's undo X-measures its top wire (bit m9) and today repays (-1)^(m9 top) in two parts: Cliffords (A0) at once and
B0 o8' at the narrow game's top step through a hook, two CCZ under m9 on every policy path (8 XCCZ under path & m9 in
R25c's emit = 1.0 a shot).  Simulated on ALL inputs at that HMR (tools/invsqrt_fc_s18m9frame.py), the top wire is a
function of 16 live wires (v0..v6, the control d, the window t0..t7) of degree 3, and
    top = Q + (v4 ^ d)(v4 ^ v5 ^ t7 ^ d ^ 1)(v3 ^ v4)        Q = 6 CZ terms, 2 Z terms and a sign
on every reached pattern (tools/invsqrt_fc_s18cubrank.py, tools/invsqrt_fc_s18m9q.py): ONE CCZ under m9 (0.5 a shot)
and Cliffords pay all of it; the hook pays nothing.  Active with cfg "_fc18m9": the device JSON (m9q's output).
  invsqrt_fc18_patch_m9.py SRCDIR          (a source COPY under ~/qcc/ws/invsqrt_fc/; edits invsqrt_fc17_skipe.py)"""
import os
import sys


def main():
    src = os.path.abspath(os.path.expanduser(sys.argv[1]))
    assert "/invsqrt_fc/" in src, "only a source copy in FC's folder"
    p = os.path.join(src, "invsqrt_fc17_skipe.py")
    s = open(p).read()
    old1 = '''    m9 = c.hmr(top)                                # digit 23 X-measured (P9)
    emit_prod(c, [F(o7, v4, v6, d), F(v4, v5, v6, d, ONE)], m9)        # A0 = P Q ^ v6 o7'
    emit_prod(c, [F(v6), F(o7)], m9)
'''
    new1 = '''    m9 = c.hmr(top)                                # digit 23 X-measured (P9)
    _m9dev = (sys.modules["invsqrt_fk_w24ret"]._STATE.get("cfg") or {}).get("_fc18m9")   # FC18-M9
    _n3 = 0
    if _m9dev:                                     # the WHOLE m9 phase here: one CCZ under m9 and Cliffords
        import json as _json
        _os = __import__("os")                     # (r24_block imports os locally further down)
        _dv = _json.load(open(_os.path.expanduser(_m9dev)))
        _rl = dict(("v%d" % _j, v[_j]) for _j in range(7))
        _rl["d"] = d
        _rl.update(("t%d" % _j, tgt[_j]) for _j in range(8))
        if _os.environ.get("FC18_M9DBG"):
            print("FC18-M9 wires: top %s t8 %s d %s v %s tgt %s" % (top, t8, d, list(v), list(tgt)), file=sys.stderr)
        _n3 = emit_prod(c, [F(*([_rl[_r] for _r in _ws] + ([ONE] if _k else []))) for _ws, _k in _dv["forms"]], m9)
        for _a, _b in _dv["quad"]:
            c.cz(_rl[_a], _rl[_b], m9)
        for _a in _dv["lin"]:
            c.z(_rl[_a], m9)
        if _dv["const"]:
            c.neg(m9)
    else:
        emit_prod(c, [F(o7, v4, v6, d), F(v4, v5, v6, d, ONE)], m9)        # A0 = P Q ^ v6 o7'
        emit_prod(c, [F(v6), F(o7)], m9)
'''
    assert s.count(old1) == 1, "m9 site not found once"
    s = s.replace(old1, new1)
    old2 = '''    def hook(c_, src_, tgt_, carry, n_, f, g=None):
        cond = m9 if g is None else c_.band(((g, 1), (m9, 1)))   # fcgen's policy paths: under the path's bit too
'''
    new2 = '''    def hook(c_, src_, tgt_, carry, n_, f, g=None):
        if _m9dev:                                 # FC18-M9: paid at the measurement
            return
        cond = m9 if g is None else c_.band(((g, 1), (m9, 1)))   # fcgen's policy paths: under the path's bit too
'''
    assert s.count(old2) == 1, "hook not found once"
    s = s.replace(old2, new2)
    old3 = '''    nccz = [0]

    def hook('''
    new3 = '''    nccz = [_n3]

    def hook('''
    assert s.count(old3) == 1, "nccz not found once"
    s = s.replace(old3, new3)
    if "\nimport os\n" not in s.split('"""', 2)[2][:400]:
        s = s.replace("\nimport sys\n", "\nimport os\nimport sys\n", 1)
    open(p, "w", newline="\n").write(s)
    print("FC18-M9 patched", p)


if __name__ == "__main__":
    main()
