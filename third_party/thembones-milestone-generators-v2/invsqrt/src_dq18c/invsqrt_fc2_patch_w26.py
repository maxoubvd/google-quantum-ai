#!/usr/bin/env python3
"""invsqrt_fc2_patch_w26.py -- FC2 patch "w26" for FC's g8 (own code; FC's and F1's files are read-only, their
source is copied in memory): MAIN's word cut to 26 slots and F1's ring endpoint at N = 26 (the line's ring size).

  MAIN   rows 26 .. blast only (cfg "_fc2_blast": 1): the tail loop stops at row 1; the final sign after row 1 takes
         slot 1, digits 1..25 slots 2..26 (27-slot indexing kept inside st_main: slot 0 = None, so FC's w25 / w26late
         branches and F1's word-bit programs keep their indices); st_main returns the 26 real tokens; two leftover
         window wires instead of one (FC's relations clear what is affine).
  END    F1's ring at N = len(word): base slot angles 2^(p + 54 - N), wrap features 2^(p+54) - 2^(p + 54 - N)
         (N = 27 gives FC's values), the constant token at index ep_const - blast, the root-bit retire programs
         read the 27-slot indexing. Rotations: N + 7 = 33.
  Phase model: the word's bit 0 is dropped (R37 model main_rows(blast 1, slot "zero")): converged LP margin of
  today's arithmetic 4.172e7 on all inputs (tools/invsqrt_fc2_slot26.py).
st_main is rebuilt from FC's file with F1's three replacements (invsqrt_fk_w24ret, when its patch is live) and
FC2's cut, then exec'd into FC's module; apply this patch BEFORE h1b (cfg "fc2": ["norm", "w26", "h1b"]).
  invsqrt_fc2_patch_w26.py check SRC N CFG.json 'OVERRIDES'    word = model (26 tokens), x restored, ancillas 0
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

CUT = [("    for b in range(br - 1, -1, -1):\n",
        "    for b in range(br - 1, int(cfg.get(\"_fc2_blast\", 0)) - 1, -1):\n"),
       ("    for b in range(0, br):\n        c.cx(v[0], digit[b])\n",
        "    for b in range(int(cfg.get(\"_fc2_blast\", 0)), br):\n        c.cx(v[0], digit[b])\n"),
       ("    word = [fin] + [digit[b] for b in range(0, 26)]\n",
        "    _bl = int(cfg.get(\"_fc2_blast\", 0))\n"
        "    word = [None] * _bl + [fin] + [digit[b] for b in range(_bl, 26)]\n"),
       ("    if mbu:\n        _init_mbu(c, us, E, win, dsrc, v)\n",
        "    if cfg.get(\"_k4fuse\"):   # FC2 fourth holder: first add fused with row 25 (invsqrt_fc2_k4fuse)\n"
        "        import invsqrt_fc2_k4fuse as _K4F\n"
        "        _K4F.fused_top(c, cfg, A, E, win, pos(25), dsrc, v, m)\n"
        "    elif mbu:\n        _init_mbu(c, us, E, win, dsrc, v)\n"),
       ("    digit = {}\n    for b in range(25, br, -1):\n",
        "    digit = {}\n    if cfg.get(\"_k4fuse\"):\n        digit[25] = None\n"
        "    for b in range(24 if cfg.get(\"_k4fuse\") else 25, br, -1):\n"),
       ("    return word, [w for w in u if w not in word]\n",
        "    return (([w_ for w_ in word if w_ is not None] if cfg.get(\"_fc2_blast\") else word),\n"
        "            [w for w in u if w not in word])\n")]


def _func_src(text, name):
    i = text.index("\ndef %s(" % name) + 1
    j = text.find("\ndef ", i + 1)
    return text[i:j + 1] if j >= 0 else text[i:]


def _make_st_end(F8):
    def st_end_f1(c, word, F, ang, feat, cfg):
        """FC's st_end_f1 (copied) with the ring size N = len(word) and the 27-slot indexing for the programs."""
        import invsqrt_fk_epN as EPN
        N = len(word)
        bl = int(cfg.get("_fc2_blast", 0))
        cap, prof = int(cfg.get("ep_cap", 6)), cfg.get("ep_profile", "generic")
        if cfg.get("ep_auto"):
            cap = min(int(cfg.get("ep_max", 9)), c.cap - c.live)
        word27 = [None] * bl + list(word)
        retire = []
        pdir = os.path.expanduser(cfg.get("ep_prog_dir", "~/qcc/ws/invsqrt_fc/fcsup"))
        for n in cfg.get("ep_retire", []):
            prog = json.load(open(os.path.join(pdir, "sc8_%s.prog.json" % n)))
            retire.append((int(n[1:]), F8._make_rebuild(prog, word27)))
        const = cfg.get("ep_const")
        ci = None if const is None else int(const) - bl
        if cfg.get("ep_lend") and len(retire) >= 2 and const is not None:
            retire = F8._lend_wires(c, retire, word[ci])
        base_ang = [1 << (p + 54 - N) for p in range(N)]
        base_feat = [((1 << (p + 54)) - (1 << (p + 54 - N))) % (1 << 64) for p in range(10)]
        ang = list(ang)[:N] if ang is not None else base_ang
        feat = feat if feat is not None else base_feat
        if retire:
            cap += len(retire)
            if os.environ.get("R63_CAP"):                # R63-INV-GHOST (my copy only): lanes capped below the retires
                cap = min(cap, int(os.environ["R63_CAP"]))
                print("r63: %d retires, lanes capped at %d" % (len(retire), cap), flush=True)
            net = EPN.best_net(N, cap, profile="line" if const is not None else prof)
            if prof == "clifford" or const is not None:
                feat = list(feat)[:7]
            if const is not None:
                extra = {"const_lane": True} if cfg.get("ep_const_lane") else {}
                if cfg.get("ep_const_early"):              # FC sixth holder: the constant's wire, early lane
                    extra["const_early"] = True
                return EPN.emit_ring(c, word, F, net, ang, feat, cap, const_token=ci, const_one=True,
                                     retire=retire, **extra)
            return EPN.emit_ring(c, word, F, net, ang, feat, cap, retire=retire)
        net = EPN.best_net(N, cap, profile="line" if const is not None else prof)
        if prof == "clifford" or const is not None:
            feat = list(feat)[:7]
        if const is not None:
            return EPN.emit_ring(c, word, F, net, ang, feat, cap, const_token=ci, const_one=True)
        return EPN.emit_ring(c, word, F, net, ang, feat, cap)
    return st_end_f1


def apply(F8, cfg):
    done = F8.__dict__.setdefault("_fc2_done", set())
    if "w26" in done:
        return
    assert "h1b" not in done, "apply w26 before h1b (it rebuilds st_main from source)"
    s = _func_src(open(F8.__file__).read(), "st_main")
    if "_F1W24" in F8.__dict__:                        # F1's w24 / w23 retire patch is live: same replacements
        import invsqrt_fk_w24ret as W
        for old, new in ((W.A_OLD, W.A_NEW), (W.B_OLD, W.B_NEW), (W.C_OLD, W.C_NEW)):
            assert s.count(old) == 1, old
            s = s.replace(old, new)
    for old, new in CUT:
        assert s.count(old) == 1, old
        s = s.replace(old, new)
    exec(compile(s, F8.__file__ + " [FC2 w26 patch]", "exec"), F8.__dict__)
    F8.__dict__["st_end_f1"] = _make_st_end(F8)
    done.add("w26")


def check(src, nsamp, cfgfile, over, seed=11):
    """the cut against today's circuit, both simulated: tokens (blast 1) = word27 >> 1 exactly [P: Y(blast 1) is
    even and Y27 = Y(blast 1) + [R_0 >= 0]] on sampled inputs + edges; x restored; ancillas zero; counts."""
    src = os.path.abspath(os.path.expanduser(src))
    sys.path.insert(0, src)
    import invsqrt_fw_fcpolrun as W
    F8, K, _ = W.load(src)
    cfg_in = json.load(open(os.path.expanduser(cfgfile)))
    cfg_in.update(over)
    rng = np.random.default_rng(seed)
    x = np.unique(np.concatenate([rng.integers(1, F8.MAX_X + 1, nsamp), np.arange(1, 300),
                                  np.arange(F8.MAX_X - 300, F8.MAX_X + 1)]))
    res = {}
    for bl in (0, int(cfg_in.get("_fc2_blast", 1))):
        cfg = W.cfg_for(F8, src, dict(cfg_in, _fc2_blast=bl))
        c, regs = F8.build_full(dict(cfg))
        sim = K.Sim(len(x), c.next + 2)
        sim.set_word(regs["X"], x)
        nf = regs["nfwd"]
        sim.run(c.ops, 0, nf)
        Yc = sim.get_word(regs["word"])
        sim.run(c.ops, nf, None)
        others = [w for w in range(c.next) if w not in regs["X"]]
        res[bl] = dict(Yc=Yc, x_ok=bool(np.array_equal(sim.get_word(regs["X"]), x)),
                       anc_ok=bool(not sim.q[others].any()), exp=c.toffolis_exp(), static=c.toffolis(),
                       peak=c.peak, wires=c.next, tokens=len(regs["word"]), rest=len(regs["rest"]),
                       n_zps=len(regs.get("epinfo", {}).get("zps", [])), stage_peaks=c.stage_peaks)
    b1 = max(res)
    out = dict(n=len(x), tok_bad=int((res[b1]["Yc"] != (res[0]["Yc"] >> b1)).sum()))
    for bl, r in res.items():
        out["b%d" % bl] = {k_: v_ for k_, v_ in r.items() if k_ != "Yc"}
    return out


if __name__ == "__main__":
    if sys.argv[1] == "check":
        print(json.dumps(check(sys.argv[2], int(sys.argv[3]), sys.argv[4], json.loads(sys.argv[5]) if len(sys.argv) > 5 else {})))
