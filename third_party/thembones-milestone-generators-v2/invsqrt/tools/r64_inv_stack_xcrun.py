#!/usr/bin/env python3
"""r64_inv_stack_xcrun.py -- R64-INV-STACK (designer): COPY of F4's runner invsqrt_fx_xcrun.py (same modes and
arguments) + ONE knob: env I39_CAPS="H1=39,DIV1=39,DIV2=39,END=39" plans the tags named at their own cap (the form
of I39's slot check: MAIN2 = the MAIN after the stage; other tags: the cfg's cap), so that a FILE can be emitted
with some stages already one wire narrower.
  r64_inv_stack_xcrun.py emit SRC OUT ANGLES.json CFG.json
Safe to run: reads SRC, writes only OUT."""
import json
import os
import sys

src = os.path.abspath(os.path.expanduser(sys.argv[2]))
sys.path.insert(0, src)
import invsqrt_fw_fcpolrun as R  # noqa: E402

_load = R.load


def load(s):
    F8, K, w = _load(s)
    mode = sys.argv[1]
    cfg_in = json.loads(sys.argv[5]) if mode == "run" else json.load(open(sys.argv[5]))
    cfg = R.cfg_for(F8, s, cfg_in)
    if os.environ.get("FX_XC", "1") != "0":
        import invsqrt_fx_xc as XC
        XC.apply(sys.modules["invsqrt_fc_arith"], cfg)
        print("fx xc applied: XC=%r XT=%r XTK=%r" % (os.environ.get("R41_XC", "DIV1"), os.environ.get("R41_XT", ""),
                                                     os.environ.get("R41_XTK", "")), flush=True)
    caps = dict((a_.split("=")[0], int(a_.split("=")[1])) for a_ in os.environ.get("I39_CAPS", "").split(",") if a_)
    if caps:
        otag = K.Circ.tag

        def tag(self, name):
            n2 = sum(1 for t_ in self.tags if t_[1] == "MAIN") + (1 if name == "MAIN" else 0)
            key = "MAIN2" if (name == "MAIN" and n2 >= 2) else name
            if name != "UNDO":
                self.cap = caps.get(key, cfg["cap"])
            return otag(self, name)
        K.Circ.tag = tag
        print("per-tag caps:", caps, "(others %d)" % cfg["cap"], flush=True)
    return F8, K, w


R.load = load
R.main()
