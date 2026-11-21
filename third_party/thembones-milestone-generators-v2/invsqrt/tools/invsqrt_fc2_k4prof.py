#!/usr/bin/env python3
"""invsqrt_fc2_k4prof.py -- FC2 fourth holder: expected Toffolis per region, before and after P1's affine pass.

  invsqrt_fc2_k4prof.py SRC SET.json ANG.json AFFINE.json ['{overrides}'] [nrows]

Builds SRC exactly as F3's runner `emit` does (env knobs as for the build), labels every op of the final circuit:
  MAIN's forward: "init" (ops of MAIN before row 25's adder call: the first add), "r25", "r24", ... (one label per
  top-level A.add / A.fadd / A.frame call in MAIN, rows in call order), "MAIN" (the rest);
  other stages by their tag; the endpoint "END"; the mirror "U:<label of the forward op it undoes>" (an explicit
  undo block takes the label of its forward block).
Ops are mapped to kmb op indices the way the kit's kmx writer emits them, and the pass's site weights (AFFINE.json
of the SAME build) are summed per label. Prints JSON: {label: [before, removed, after]} for MAIN's top rows and
every stage, and the totals (before must equal the build's toffolis_exp).
"""
import json
import os
import sys
from collections import OrderedDict

TOFW = {"CCX": 1.0, "CCZ": 1.0, "CSWAP": 1.0, "ANDC": 1.0, "CANDC": 0.5}
LADDER = [False]
MAINSTART = {}
STARGS = {}


def build(src, setf, angf, over):
    src = os.path.abspath(os.path.expanduser(src))
    sys.path.insert(0, src)
    import invsqrt_fw_fcpolrun as R
    F8, K, _ = R.load(src)
    A = sys.modules["invsqrt_fc_arith"]
    cfg = R.cfg_for(F8, src, json.load(open(setf)))
    cfg.update(over)
    LADDER[0] = bool(cfg.get("_ladder"))
    K.KMX_CBLOCKS = bool(cfg.get("kmx_cblocks"))
    a = json.load(open(angf))
    calls = []
    depth = [0]

    def wrap(name):
        orig = getattr(A, name)

        def f(c, *aa, **kw):
            tag = c.tags[-1][1] if getattr(c, "tags", None) else "?"
            top = depth[0] == 0
            depth[0] += 1
            i0 = len(c.ops)
            try:
                return orig(c, *aa, **kw)
            finally:
                depth[0] -= 1
                if top and tag == "MAIN":
                    calls.append((id(c), name, i0, len(c.ops)))
        f.__name__ = orig.__name__
        return f
    for nm in ("add", "fadd", "frame"):
        setattr(A, nm, wrap(nm))
    # MAIN's own start (after H1, which runs inside MAIN's tag): wrap st_main once w26 has rebuilt it, before h1b
    # wraps it (build_full's patch loop then finds w26 applied and h1b wraps this marker)
    if "w26" in (cfg.get("fc2") or []):
        import importlib
        sys.path.insert(0, src)
        W26 = importlib.import_module("invsqrt_fc2_patch_w26")
        sys.path.pop(0)
        W26.apply(F8, cfg)
    orig_sm = F8.st_main

    def st_main_mark(c, U, v, cfg_):
        MAINSTART[id(c)] = len(c.ops)
        STARGS[id(c)] = (list(U), list(v))
        return orig_sm(c, U, v, cfg_)
    F8.__dict__["st_main"] = st_main_mark
    c, regs = F8.build_full(cfg, [int(t) for t in a["ang"]], [int(t) for t in a["feat"]])
    calls = [t[1:] for t in calls if t[0] == id(c)]
    return c, regs, calls, K


def labels_for(c, regs, calls, nrows):
    nf = regs["nfwd"]
    ops = c.ops
    lab = [None] * len(ops)
    # stage tags (forward + END)
    marks = sorted(c.tags, key=lambda t: t[0])
    for (i0, nm, _), nxt in zip(marks, marks[1:] + [(len(ops), "END*", 0)]):
        for i in range(i0, min(nxt[0], len(ops))):
            lab[i] = nm
    # MAIN's adder calls: frames of the fadd rows are folded into their row (fadd = frame + add + frame)
    main0 = [t[0] for t in marks if t[1] == "MAIN"][0]
    ms = MAINSTART.get(id(c), main0)
    for i in range(main0, ms):
        lab[i] = "H1f"                              # H1's stage (inside MAIN's tag, before st_main proper)
    rows = []
    for name, i0, i1 in sorted(calls, key=lambda t: t[1]):
        if i0 >= ms:
            rows.append((name, i0, i1))
    if rows:
        for i in range(ms, rows[0][1]):
            lab[i] = "init"                         # ops before the first recorded call (ladder bypasses A.add)
        # first call is the first add unless F1's ladder replaced it (its hook wraps A.add outside this wrapper)
        first_is_init = not LADDER[0]
        k = 0
        for j, (name, i0, i1) in enumerate(rows):
            if j == 0 and first_is_init:
                L = "init"
            else:
                L = "r%d" % (25 - k)
                k += 1
            if k > nrows and L != "init":
                L = "MAINrows"
            for i in range(i0, i1):
                lab[i] = L
    # the mirror: replicate Circ.mirror_blocks with labels
    fwd = list(range(nf))
    out = []
    i = nf - 1
    seg = []

    def flush(seg):
        for jj in seg:                            # Circ.mirror: reversed order, BBEG / BEND dropped, 1:1 otherwise
            if ops[jj][0] in ("BBEG", "BEND"):
                continue
            out.append("U:" + str(lab[jj]))
    while i >= 0:
        o = ops[fwd[i]]
        if o[0] == "BEND":
            flush(seg)
            seg = []
            bid = o[1]
            out.extend(["U:" + str(lab[i])] * len(c.blocks[bid]))
            while ops[i][0] != "BBEG" or ops[i][1] != bid:
                i -= 1
            i -= 1
            continue
        seg.append(i)
        i -= 1
    flush(seg)
    nend = len(ops) - nf - len(out)
    assert nend >= 0, (len(ops), nf, len(out))
    for i in range(nf, nf + nend):
        lab[i] = "END"
    for j, L in enumerate(out):
        lab[nf + nend + j] = L
    return lab


MARK = 9000000


def kmb_labels(c, regs, lab, K):
    """label of every kmb op: the (patched) kit writer renders the circuit with a marker line "X q<MARK+r>" before
    each change of label; every other non-register line is one kmb op, in order."""
    names = []
    rid = {}
    marked = []
    prev = object()
    for o, L in zip(c.ops, lab):
        if L != prev:
            if L not in rid:
                rid[L] = len(names)
                names.append(L)
            marked.append(("X", MARK + rid[L], None))
            prev = L
        marked.append(o)
    saved = c.ops
    c.ops = marked
    try:
        txt = K.kmx_text(c, regs["X"])
    finally:
        c.ops = saved
    out = []
    cur = None
    for ln in txt.splitlines():
        if not ln or ln.startswith("APPEND_TO_REGISTER") or ln.startswith("REGISTER"):
            continue
        p = ln.split()
        if p[0] == "X" and len(p) == 2 and p[1].startswith("q") and int(p[1][1:]) >= MARK:
            cur = names[int(p[1][1:]) - MARK]
            continue
        out.append(cur)
    return out


def main():
    src, setf, angf, affp = sys.argv[1:5]
    over = json.loads(sys.argv[5]) if len(sys.argv) > 5 and sys.argv[5].strip() else {}
    nrows = int(sys.argv[6]) if len(sys.argv) > 6 else 3
    c, regs, calls, K = build(src, setf, angf, over)
    lab = labels_for(c, regs, calls, nrows)
    klab = kmb_labels(c, regs, lab, K)
    sites = json.load(open(affp))
    removed = {}
    bad = 0
    for s in sites:
        i = int(s["op"])
        if i >= len(klab):
            bad += 1
            continue
        removed[klab[i]] = removed.get(klab[i], 0.0) + float(s.get("w", 1.0))
    groups = OrderedDict()
    for o, L in zip(c.ops, lab):
        groups.setdefault(L, []).append(o)
    before = OrderedDict((L, c.toffolis_exp(g)) for L, g in groups.items())
    before = OrderedDict((L, v) for L, v in before.items() if v)
    keys = list(before.keys()) + [k for k in removed if k not in before]
    res = OrderedDict((k, [before.get(k, 0.0), removed.get(k, 0.0), before.get(k, 0.0) - removed.get(k, 0.0)])
                      for k in keys)
    nk = None
    if os.environ.get("K4KMB"):                    # the build's own kmb: its op count must equal the mapped count
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import kmb as KB
        nk = len(KB.read(os.environ["K4KMB"]).ops)
    tot_b = sum(v[0] for v in res.values())
    tot_r = sum(v[1] for v in res.values())
    print(json.dumps({"toffolis_exp": c.toffolis_exp(), "sum_before": tot_b, "removed": tot_r,
                      "after": tot_b - tot_r, "sites": len(sites), "unmapped_sites": bad, "peak": c.peak,
                      "kmb_ops_mapped": len(klab), "kmb_ops_file": nk, "regions": res}))


if __name__ == "__main__":
    main()
