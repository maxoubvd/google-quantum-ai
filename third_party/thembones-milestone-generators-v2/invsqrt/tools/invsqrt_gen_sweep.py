"""invsqrt_gen: build several option sets of one variant without the exact check and print, for each, the peak
wire count, the static Toffoli count and the AND count of every lookup block (own tool, workstream invsqrt_gen).

usage: invsqrt_gen_sweep.py <variant> <base key=value ...> :: <key=value ...> :: <key=value ...> ...
  every "::" starts one option set that is laid over the base set. The lookup programs found on the way are
  cached by invsqrt_gen_pool.lookup_block as in a normal build.
  sample=N (in the base or in a set; not passed to the builder): also run the exact simulator on N random
  inputs plus the 2 x 2,048 edge inputs and print pass / fail, the mean executed Toffolis, the rotation count and
  the predicted score Q x (Toffolis + (T + rotations x T per rotation) / 4), with the T per rotation of the
  tables' margin on ALL inputs (model check, a few seconds per table set).
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_gen_model as M
import invsqrt_gen_build as B

_MARGIN = {}


def full_margin(P, tabs):
    """Smallest margin over all inputs (units of 2^-64 turn) of the nominal pipeline with these tables."""
    k = (P.key(), tabs.get("sha"))
    if k not in _MARGIN:
        _MARGIN[k] = M.check(P, tabs)[0]["min_margin"]
    return _MARGIN[k]


def main():
    variant = sys.argv[1]
    parts = " ".join(sys.argv[2:]).split("::")
    base = dict(kv.split("=") for kv in parts[0].split())
    pf = set(M.Params.__dataclass_fields__)
    for var in parts[1:] or [""]:
        kw = dict(base)
        kw.update(dict(kv.split("=") for kv in var.split()))
        tabname = kw.pop("tabs", None)                            # tables from tabs/<name>.json, not from our fit
        kw = {k: int(v) for k, v in kw.items()}
        sample = kw.pop("sample", 0)
        P = M.Params(**{k: v for k, v in kw.items() if k in pf})
        tabs = M.load_tables(tabname) if tabname else M.fit_tables(P)
        t0 = time.time()
        try:
            c = B.VARIANTS[variant](P, tabs, **{k: v for k, v in kw.items() if k not in pf})
        except Exception as ex:                                   # a set that does not build is reported, not fatal
            print(f"{var.strip()}: FAILED {type(ex).__name__}: {ex}", flush=True)
            continue
        st = c.static_counts()
        lut = c.meta.get("lut", {})
        line = (f"{var.strip() or '(base)'}: Q {c.peak}  static CCX {st.get('CCX', 0)}  lookups "
                f"{sum(v['ands'] for v in lut.values())} = " + " ".join(f"{k}:{v['ands']}" for k, v in lut.items())
                + f"  helpers {c.meta.get('n_helpers')}")
        if sample:
            rng = np.random.default_rng(0)
            xs = np.unique(np.concatenate([rng.integers(1, M.MAX_X + 1, sample), np.arange(1, 2049),
                                           np.arange(M.MAX_X - 2048, M.MAX_X + 1)]))
            r = B.check(c, P, tabs, xs)
            mm = full_margin(P, tabs)
            if kw.get("angf"):
                # angf (G2, 2026-11-12): tuned rotation angles; the margin is the circuit's own, measured on ALL
                # inputs (the nominal model margin above does not see the angles)
                mm = B.check(c, P, tabs)["min_margin"]
            line += f"  | sample {'ok' if r['ok'] else 'FAIL'}  Tof {r['tof']:.2f}  rot {r['n_rot']} + {r['n_t']} T"
            if mm > 0 and r["n_rot"]:
                tper = 0.56 * math.log2(1 / (math.pi * mm / 2 ** 64 / r["n_rot"])) + 5.3
                line += f"  score {r['Q'] * (r['tof'] + (r['n_t'] + r['n_rot'] * tper) / 4):,.0f}"
            line += "  stages " + " ".join(f"{n}:{t:.1f}" for n, t, _ in r["stages"] if t)
        print(line + f"  [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
