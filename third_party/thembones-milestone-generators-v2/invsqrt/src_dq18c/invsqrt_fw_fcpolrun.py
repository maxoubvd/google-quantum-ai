#!/usr/bin/env python3
"""invsqrt_fw_fcpolrun.py -- F3 (ninth holder, 2026-11-14): FC's family-C build from a FROZEN src folder with F1's
w24 retire patch (src/invsqrt_fk_w24ret.py, when present) and F3's adaptive undo (tools/invsqrt_fw_fcpol.py) applied
in memory. Same modes as F1's runner:
  emit SRC OUT ANGLES.json CFG.json    build with given angles, write OUT/invsqrt.kmx and OUT/invsqrt.kmb
  run  SRC WD OUT CFGJSON               FC's run8 (all inputs + angle LP + kmx) through the patches
The undo blocks of shape (L 8, C 3) take F3's policy (expected 12.0 against 12.5); nothing else changes, so the
phase function and the angles are those of the base build.
"""
import hashlib
import importlib
import json
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))


def load(src):
    sys.path.insert(0, src)
    import invsqrt_fc_kit as K  # noqa: F401  (from src)
    import invsqrt_fc_g8 as F8
    for f in sorted(os.listdir(src)):
        if f.endswith(".py") and "-" not in f and f.startswith("invsqrt_fc") and "run" not in f:
            __import__(f[:-3])
    w = None
    if os.path.exists(os.path.join(src, "invsqrt_fk_w24ret.py")):
        W = importlib.import_module("invsqrt_fk_w24ret")
        w = W.apply(F8, sys.modules["invsqrt_fc_kit"])
    sys.path.append(TOOLS)
    if os.environ.get("FCPOL", "1") == "1":        # FCPOL=0: the base build (reproduces FC's raw kmb)
        import invsqrt_fw_fcpol as P   # imports invsqrt_fc_kit / invsqrt_fc_narrow: the src copies (already loaded)
        assert P.K is sys.modules["invsqrt_fc_kit"] and P.NW is sys.modules["invsqrt_fc_narrow"]
        P.apply()
    gen = os.environ.get("FCGEN", "")              # e.g. "7,2,1;8,2,1": (L, C, ccz) shapes for the general emitter
    if gen:
        import invsqrt_fw_fcgen as G
        assert G.K is sys.modules["invsqrt_fc_kit"] and G.NW is sys.modules["invsqrt_fc_narrow"]
        G.apply({(int(a), int(b)): int(z) for a, b, z in (t.split(",") for t in gen.split(";"))})
    if os.environ.get("C1ADD", "0") == "1":        # rows at < 2 clean wires (cuccaro_fwd): F3's c1_add (top carry
        import invsqrt_fw_c1add as C1             # straight into the top bit; spare clean wire as a clean position)
        pred = None
        if os.environ.get("C1PRED") == "ctop":     # only rows whose top-position source bit is a constant (MAIN 25, 24)
            def pred(s, t, ci):
                s = list(s) + [None] * (len(t) - len(s))
                return len(t) >= 3 and not isinstance(s[len(t) - 2], int)
        C1.apply(sys.modules["invsqrt_fc_arith"], sys.modules["invsqrt_fc_narrow"], pred)
    if os.environ.get("C1UNDO", "0") == "1":       # hybrid_add / cuccaro_add everywhere -> c1_add, same budgets
        import invsqrt_fw_c1add as C1
        C1.apply_undo(sys.modules["invsqrt_fc_narrow"])
    if os.environ.get("CHIPROG8", "0") == "1" and "invsqrt_fc2_patch_h1c" in sys.modules:   # H1' (m 7): FC2's
        import invsqrt_fw_chiprog as CP            # 17-AND carry rule -> F3's SAT program (8 ANDs; CHIPROG=path)
        CP.apply_h1c(sys.modules["invsqrt_fc2_patch_h1c"])
    if os.environ.get("CHI6K5", "0") == "1":       # FC2's h1c chi mode "f3six6" takes F3's 5-AND minimum (CHI6_K5)
        import invsqrt_fw_chiprog as CP
        CP.CHI6_K6 = CP.CHI6_K5
    if os.environ.get("AFFLEND"):                  # FC's "afflend" (FC4 22:45Z, copied for srcs without the hook):
        _afflend(sys.modules["invsqrt_fc_arith"], json.loads(os.environ["AFFLEND"]))   # outermost add wrapper
    return F8, sys.modules["invsqrt_fc_kit"], w


def _afflend(A, entries):
    """around the matching adder call (stage tag, first target wire, carry-in), a live non-operand wire that is an
    affine function of other live wires before (pre) and after (post) the add is cleared by CNOTs and freed, so the
    add's explicit undo block gets one more clean wire; restored by CNOTs after the add (the mirror repeats them)."""
    orig = A.add

    def add(c, src, tgt, *a, **kw):
        tag = c.tags[-1][1] if getattr(c, "tags", None) else None
        hit = None
        for e in entries:
            if tag == e["stage"] and len(tgt) and tgt[0] == e["tgt0"] and kw.get("cin", a[0] if a else None) == e["cin"]:
                hit = e
                break
        if hit is None:
            return orig(c, src, tgt, *a, **kw)
        w_ = int(hit["wire"])
        for q in hit["pre"]:
            c.cx(q, w_)
        if hit.get("cpre"):
            c.x(w_)
        c.free_([w_])
        r = orig(c, src, tgt, *a, **kw)
        c.free.remove(w_)
        c.live += 1
        c.peak = max(c.peak, c.live)
        for q in hit["post"]:
            c.cx(q, w_)
        if hit.get("cpost"):
            c.x(w_)
        return r
    A.add = add


def cfg_for(F8, src, cfg_in):
    cfg = dict(F8.CFG)
    cfg.update(cfg_in)
    cfg["fc2_src"] = src
    cfg["ep_prog_dir"] = os.path.join(src, "fcsup")
    return cfg


def main():
    mode = sys.argv[1]
    src = os.path.abspath(os.path.expanduser(sys.argv[2]))
    F8, K, w = load(src)
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
        print("emit %s: sha %s toffolis_exp %.2f static %d wires %d peak %d kmb qubits %d; w24 log %s" % (
            out, sha[:20], c.toffolis_exp(), c.toffolis(), c.next, c.peak, kc.num_qubits, getattr(w, "log", None)))
    elif mode == "run":
        wd, out, cfgs = sys.argv[3], sys.argv[4], sys.argv[5]
        cfg = cfg_for(F8, src, json.loads(cfgs))
        import invsqrt_fc_run8 as R8
        sys.argv = ["invsqrt_fc_run8.py", wd, out, json.dumps(cfg)]
        R8.main()
        print("w24 log", getattr(w, "log", None), flush=True)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
