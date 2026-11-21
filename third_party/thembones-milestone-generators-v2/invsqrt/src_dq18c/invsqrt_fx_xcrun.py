#!/usr/bin/env python3
"""invsqrt_fx_xcrun.py -- F4 (ninth holder): F3's runner (invsqrt_fw_fcpolrun.py of a FROZEN src copy) with the
dividers' tail wrapper invsqrt_fx_xc.py (a copy of R41-INV-MAIN's r41_inv_main_xc.py) applied after the seats'
patches. Lives in the src copy next to the runner (same modes, same arguments):
  invsqrt_fx_xcrun.py run  SRC WD OUT CFGJSON
  invsqrt_fx_xcrun.py emit SRC OUT ANGLES.json CFG.json
Env as the designer's wrapper: R41_XC (default DIV1; "" = none), R41_XT, R41_XTK, R41_XC4, ...; FX_XC=0 skips it.
"""
import json
import os
import sys

src = os.path.abspath(os.path.expanduser(sys.argv[2]))
sys.path.insert(0, src)
import invsqrt_fw_fcpolrun as R  # noqa: E402

_load = R.load


def load(s):
    F8, K, w = _load(s)
    if os.environ.get("FX_XC", "1") != "0":
        mode = sys.argv[1]
        cfg_in = json.loads(sys.argv[5]) if mode == "run" else json.load(open(sys.argv[5]))
        cfg = R.cfg_for(F8, s, cfg_in)
        import invsqrt_fx_xc as XC
        XC.apply(sys.modules["invsqrt_fc_arith"], cfg)
        print("fx xc applied: XC=%r XT=%r XTK=%r" % (os.environ.get("R41_XC", "DIV1"), os.environ.get("R41_XT", ""),
                                                     os.environ.get("R41_XTK", "")), flush=True)
    return F8, K, w


R.load = load
R.main()
