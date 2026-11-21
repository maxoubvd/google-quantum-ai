#!/usr/bin/env python3
"""r36_inv_corr_patch.py -- R36 designer seat INV-CORR (own file): the higher-order stage of
research/panel/R36-INV-CORR.md written with FC's kit operations, and a wrapper that builds FC's generator with it
IN MEMORY (FC's and FC2's files are imported read-only; nothing of theirs is edited or written).

  python -I tools/r36_inv_corr_patch.py FOLDER ['{cfg overrides}'] [NSAMP]
FOLDER: an FC build folder (cfg.json + angles.json), e.g. ~/qcc/ws/invsqrt_fc/best/nw44u8a.  Overrides as in
tools/invsqrt_fc2_count.py; extra keys read here: "r36_form": "hold" (default; c0 held to the mirror, use br 8)
or "dump" (c0 handed to MAIN as dividend position 8 through FC's `_pre_lz` hook, use br 7); "r36_jtop": 7, 6, 5
or 4 (how many row pairs are merged; default 6); "r36_neg2": 1 (variant v2, default) or 0 (v1); "r36_k384": 1
(variant v0: one cubic step at Tf >= 384 for both signs, one AND).
Prints expected Toffolis, peak, stage counts and, with NSAMP, the sample check (word against FC's ref_main on the
stage's integer reference, x restored, ancillas zero).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import invsqrt_fc_g8 as F8          # noqa: E402  (FC, read-only)
import invsqrt_fc_gen as G          # noqa: E402  (FC, read-only)
import invsqrt_fc_arith as A        # noqa: E402  (FC, read-only)
import r36_inv_corr_ref as REF      # noqa: E402

_ORIG = {}
_MARK = {}


def stage(c, U, cfg):
    """U: 19 wires, y in two's complement (LSB first).  Leaves U = y + h on the same wires; returns the wire c0."""
    jtop = int(cfg.get("r36_jtop", 6))
    neg2 = bool(cfg.get("r36_neg2", 1))
    k384 = bool(cfg.get("r36_k384", 0))                # kappa = [Tf >= 384] = T8 AND T7 alone (variant v0)
    k256 = bool(cfg.get("r36_k256", 0))                # kappa = [Tf >= 256] = Tf_8: no AND at all (variant v00):
    if k384 or k256:                                   # +kappa (y >= 0), -kappa (y < 0) = T8u - s, T8u = Tf_8 XOR s
        neg2 = False
    L, T = list(U[:9]), list(U[9:])
    s = T[9]
    Gd, c0, c1 = c.alloc(3)
    R = [Gd] + L + [c0, c1]
    for i in range(9):
        c.cx(s, T[i])                                  # T[i] = Tf_i
    # cubic steps: kappa = [Tf >= 320] = p ^ q ^ k2, k2 = [Tf >= 448]
    if k256:
        t8 = c.alloc(1)[0]
        c.cx(T[8], t8)
        c.cx(s, t8)                                    # t8 = T8u (the unfolded bit 8)
        A.add(c, [T[i] for i in range(9)], R[1:11], cin=t8)    # lin 1 with +T8u at the unit position
        c.cx(s, t8)
        c.cx(T[8], t8)
        c.free_([t8])
        p = q = k2 = kp = None
    else:
        p = c.and_compute(T[8], T[7])
    if k256:
        pass
    elif k384:
        q, k2 = p, None                                # q = kappa
    else:
        q = c.and_compute(T[8], T[6])
        k2 = c.and_compute(p, T[6])
        c.cx(p, q)
        c.cx(k2, q)                                    # q = kappa
    if not k256:
        c.x(s)
        kp = c.and_compute(q, s)
        c.x(s)
        A.add(c, [T[i] for i in range(9)], R[1:11], cin=kp)    # lin 1: Tf << 1, + kappa for y >= 0
        c.x(s)
        c.and_erase(q, s, kp)
        c.x(s)
    A.add(c, [T[i] for i in range(9)], R)                      # lin 2: Tf
    # The subtract words, in an order that lets the cubic wires die early: plain rows 8 .. jtop + 1, then the
    # merged rows from the smallest (j = 5) up, so that the widest OR sets run with no cubic wire alive.
    # Borrow-ins (position 0, half a unit each): words 0, 1 take kn; words 2, 3 take kn2 (variant v2).
    rows = []
    for j in range(8, jtop, -1):
        idx = [i for i in range(j) if i + j >= 8]
        words = [{i + j - 8: T[i] for i in idx}]
        if j >= jtop + 2:
            words.append({i + j - 9: T[i] for i in idx if i + j - 9 >= 0})
        rows.append(("plain", j, idx, words))
    for j in range(5, jtop + 1):
        rows.append(("merged", j, list(range(8 - j, j)), None))
    st = {"k": 0, "kn": None, "kn2": None, "cub": not k256, "q": not k256}

    def sub(word):
        k = st["k"]
        cin = None
        if k < 2 and k256:
            cin = s                                    # -s: two half-unit borrow-ins with the sign wire itself
        elif k < 2:
            if st["kn"] is None:
                st["kn"] = c.and_compute(q, s)
            cin = st["kn"]
        elif k < 4 and neg2:
            if st["kn2"] is None:
                st["kn2"] = c.and_compute(k2, s)
            cin = st["kn2"]
        st["k"] = k + 1
        src = [word.get(pos) for pos in range(len(R))]
        while src and src[-1] is None:
            src.pop()
        for t in R:
            c.x(t)
        A.add(c, src, R, cin=cin)
        for t in R:
            c.x(t)

    def cleanup(final=False):
        """at a row boundary (T unmasked): erase the cubic wires that have no reader left."""
        if st["kn"] is not None and (st["k"] >= 2 or final):
            c.and_erase(q, s, st["kn"])
            st["kn"] = None
        if st["q"] and st["kn"] is None and (st["k"] >= 2 or final):
            if k384:
                c.and_erase(T[8], T[7], p)
                st["cub"] = False
            else:
                c.cx(k2, q)
                c.cx(p, q)
                c.and_erase(T[8], T[6], q)
            st["q"] = False
        done2 = (not neg2) or st["k"] >= 4 or final
        if st["kn2"] is not None and done2:
            c.and_erase(k2, s, st["kn2"])
            st["kn2"] = None
        if st["cub"] and not st["q"] and st["kn2"] is None and done2:
            c.and_erase(p, T[6], k2)
            c.and_erase(T[8], T[7], p)
            st["cub"] = False

    for kind, j, idx, words in rows:
        for i in idx:
            c.cx(T[j], T[i])                           # T[i] = x_ij
        if kind == "plain":
            for w in words:
                sub(w)
        else:                                          # M_j = e ? all ones : 2 X_j
            c.cx(T[j], T[j + 1])                       # T[j + 1] = e_j
            e = T[j + 1]
            c.x(e)
            ors = []
            for i in idx:
                c.x(T[i])
                t = c.and_compute(e, T[i])
                c.x(t)                                 # t = e OR x_ij
                ors.append(t)
            c.x(e)
            word = {0: e}
            for k_, t in enumerate(ors):
                word[1 + k_] = t
            sub(word)
            c.x(e)
            for i, t in reversed(list(zip(idx, ors))):
                c.x(t)
                c.and_erase(e, T[i], t)
                c.x(T[i])
            c.x(e)
            c.cx(T[j], T[j + 1])
        for i in idx:
            c.cx(T[j], T[i])
        cleanup()
    assert st["k"] >= (4 if neg2 else 2), "cubic borrow-ins not all placed"
    cleanup(final=True)
    c.cx(T[5], Gd)                                     # the guard is Tf_5
    c.free_([Gd])
    for i in range(9):
        c.cx(s, T[i])                                  # unfold
    A.add(c, [c0, c1], T)                              # T += C (mod 2^10)
    c.cx(s, T[8])                                      # c1 = NOT c0 AND Tf'8 AND Tf'7 AND NOT L'8
    c.cx(s, T[7])
    c.x(c0)
    c.x(L[8])
    u1 = c.and_compute(T[8], T[7])
    u2 = c.and_compute(c0, L[8])
    u3 = c.and_compute(u1, u2)
    c.cx(u3, c1)
    c.and_erase(u1, u2, u3)
    c.and_erase(c0, L[8], u2)
    c.and_erase(T[8], T[7], u1)
    c.x(L[8])
    c.x(c0)
    c.cx(s, T[7])
    c.cx(s, T[8])
    c.free_([c1])
    return c0


def _st_main(c, U, v, cfg):
    c.tag("HO")
    a0 = len(c.ops)
    if os.environ.get("R36_DEBUG"):                    # log every adder of the stage: wires live, clean, form taken
        real_add = A.add
        log = []

        def dbg_add(c_, src, tgt, cin=None, force=None):
            n0, live0, pk0 = len(c_.ops), c_.live, c_.peak
            real_add(c_, src, tgt, cin, force)
            kinds = {o[0] for o in c_.ops[n0:]}
            form = "narrow" if "BBEG" in kinds else "plain"
            ands = sum(1 for o in c_.ops[n0:] if o[0] == "ANDC")
            log.append((len(tgt), live0, c_.cap - live0, form, ands, c_.peak))
        A.add = dbg_add
        try:
            c0 = stage(c, U, cfg)
        finally:
            A.add = real_add
        _MARK["adds"] = log
    else:
        c0 = stage(c, U, cfg)
    _MARK["ho_ops"] = (a0, len(c.ops))
    _MARK["c0"] = c0
    _MARK["U"] = list(U)
    if cfg.get("r36_form", "hold") == "dump":
        cfg["_pre_lz"] = {8: c0}                       # FC's hook: a wire already holding dividend bit 8
    c.tag("MAIN")
    return _ORIG["st_main"](c, U, v, cfg)


def _plain_y(x, cfg):
    m = cfg["m"]
    n, k, Q, d, delta = F8.ref_chain(x, m)
    a = F8.a_top(cfg)
    sg1, q1 = G.ref_div(delta.astype(np.int64) << cfg["F1"], d, a, F8.sched1_of(cfg))
    sg2, P = G.ref_div(q1 + cfg["c2"] * d, d, a, cfg["sched2"])
    return -((P + 1) >> 1), d, k, q1, P


def _ref_U(x, cfg):
    y, d, k, q1, P = _plain_y(x, cfg)
    if cfg.get("r36_k256"):
        U, c0, c1, C = REF.stage_ref(y, 256, None)
    elif cfg.get("r36_k384"):
        U, c0, c1, C = REF.stage_ref(y, 384, None)
    else:
        U, c0, c1, C = REF.stage_ref(y, 320, 448 if cfg.get("r36_neg2", 1) else None)
    if cfg.get("r36_form", "hold") == "dump":
        return U, d, k, q1, P, c0
    return U, d, k, q1, P


def apply():
    if "st_main" not in _ORIG:
        _ORIG["st_main"] = F8.st_main
        _ORIG["ref_U"] = F8.ref_U
        F8.st_main = _st_main
        F8.ref_U = _ref_U


def main():
    folder = os.path.expanduser(sys.argv[1])
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(os.path.join(folder, "cfg.json"))))
    if len(sys.argv) > 2 and sys.argv[2]:
        for k, v in json.loads(sys.argv[2]).items():
            if v is None:
                cfg.pop(k, None)
            else:
                cfg[k] = v
    a = json.load(open(os.path.join(folder, "angles.json")))
    ang, feat = [int(t) for t in a["ang"]], [int(t) for t in a["feat"]]
    apply()
    out = {}
    if len(sys.argv) > 3:
        x = F8.sample_x(int(sys.argv[3]))
        o = F8.full_check(x, dict(cfg), ang=ang, feat=feat)
        c = o["c"]
        out.update(nsamp=len(x), Y_bad=o["Y_bad"], x_ok=o["x_ok"], anc_ok=o["anc_ok"])
    else:
        c, regs = F8.build_full(dict(cfg), ang, feat)
    out.update(exp=c.toffolis_exp(), static=c.toffolis(), wires=c.next, peak=c.peak, stages=c.stage_counts(),
               stage_peaks=getattr(c, "stage_peaks", {}), form=cfg.get("r36_form", "hold"),
               jtop=cfg.get("r36_jtop", 6), neg2=cfg.get("r36_neg2", 1), cap=cfg.get("cap"), bt=cfg.get("bt"),
               br=cfg.get("br"))
    if _MARK.get("adds"):
        out["stage_adds_n_live_clean_form_ands_peak"] = _MARK["adds"]
    print(json.dumps(out))


if __name__ == "__main__":
    main()
