"""F1 (fourth holder): retire root bits across F1's endpoint in a FROZEN FC g8 build (FC's src/ folder, read-only:
imported first on sys.path, so the circuit is FC's byte for byte up to END). Each retired wire is X-measured before
the routing, lent to the walk as a lane, rebuilt after it from a minimum-AND program (pysat minand, files
~/qcc/ws/invsqrt_fk/fcsup/sc8_<name>.prog.json; inputs named q<i> = wire i, w<i> = word bit i, -1 = the constant 1;
each gate is AND(XOR left, XOR right); the last gate is computed straight into the wire, earlier gates into a temp
that is erased by HMR + CZ), then the measurement is repaid with Z.
    NWRET="q20,q21"   wires to retire, rebuilt in this order (default)
    NWCAP             lanes for the endpoint (default: FC's ep_auto cap + number of retired wires)
    python invsqrt_fk_nwret.py check SRC CFG.json [nsamp=0: all inputs]   relations at the start of END
    python invsqrt_fk_nwret.py emit  SRC OUTDIR ANGLES.json CFG.json        build, kit check on samples, kmx + kmb
"""
import json
import os
import sys

SRC = os.path.abspath(sys.argv[2])
sys.path.insert(0, SRC)
import invsqrt_fc_kit as K  # noqa: E402  (FC's frozen copies)
import invsqrt_fc_g8 as F8  # noqa: E402
import invsqrt_fc_gen as FG  # noqa: E402
# The endpoint (frozen copy of tools/invsqrt_fk_epN.py, identical on Nov 14 11:30Z) is imported only inside END, as
# FC's own st_end_f1 does: its helper puts tools/ first on sys.path, and an earlier import would make the forward
# stages' lazy imports pick FC's live tools/ copies instead of the frozen ones (measured: DIV1 336 -> 382).

PROG = os.path.expanduser("~/qcc/ws/invsqrt_fk/fcsup")
NAMES = [t for t in os.environ.get("NWRET", "q20,q21").split(",") if t]
PROBE = {}


def load(name):
    return json.load(open(os.path.join(PROG, f"sc8_{name}.prog.json")))


def wire_of(tok, word):
    return word[int(tok[1:])] if tok.startswith("w") else int(tok[1:])


def make_rebuild(prog, word):
    ins = [wire_of(t, word) for t in prog["inputs"]]
    nin, gates = len(ins), prog["gates"]
    (out,) = prog["outputs"]
    last = nin + len(gates) - 1
    assert last in out, "the last gate must feed the output"

    def rebuild(c, t):
        val = {i: ins[i] for i in range(nin)}

        def side(items, other):
            """XOR of items on one wire; returns (wire, undo list, temp)."""
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

        temps = []
        for g, gt in enumerate(gates):
            a, da = side(gt["left"], gt["right"])
            b, db = side(gt["right"], gt["left"])
            if nin + g == last:
                c.and_compute(a, b, t)
                val[nin + g] = t
            else:
                q = c.and_compute(a, b)
                val[nin + g] = q
                temps.append((g, q))
            undo(db)
            undo(da)
        for i in out:
            if i == -1:
                c.x(t)
            elif i != last:
                c.cx(val[i], t)
        for g, q in reversed(temps):                       # erase the temps (their factors are home again)
            gt = gates[g]
            a, da = side(gt["left"], gt["right"])
            b, db = side(gt["right"], gt["left"])
            c.and_erase(a, b, q)
            undo(db)
            undo(da)
    return rebuild


def joint_j8(word):
    """q19 + q20 + q21 jointly in 6 ANDs (box 4 minand, fcsup/J8_k6.prog.json, exact on FC's sc8 layout):
      g0 = (q16+q18+w23+w24+w25)(q17+w22+w25)   g1 = (w22+w24+w25)(w23+w24)   g2 = (q17+w22+w23+g0+g1)(q18+w24)
      g3 = (q17+w22+w23+w24+g0+g2)(q18+w22+w23+w24)   g4 = (w22+w24+g3)(w23+w24+w25+g1)   g5 = (w23+g4) w24
      q19 = 1+w22+w25+g1+g3   q20 = 1+w23+w24+g4   q21 = 1+w25+g5
    g3, g4, g5 go straight into the three |0> targets, g0..g2 into temps (at most 3 at once), erased by HMR + CZ."""
    W = {i: word[i] for i in range(22, 26)}
    q16, q17, q18, t19, t20, t21 = 16, 17, 18, 19, 20, 21

    def xform(c, acc, terms, one=False):
        for t in terms:
            c.cx(t, acc)
        if one:
            c.x(acc)

    def andc(c, la, lt, ra, rt, tgt=None, lone=False, rone=False, erase=False):
        xform(c, la, lt, lone)
        xform(c, ra, rt, rone)
        if erase:
            c.and_erase(la, ra, tgt)
        elif tgt is None:
            tgt = c.and_compute(la, ra)
        else:
            c.and_compute(la, ra, tgt)
        xform(c, ra, list(reversed(rt)), rone)
        xform(c, la, list(reversed(lt)), lone)
        return tgt

    def rebuild(c, _q):
        w22, w23, w24, w25 = W[22], W[23], W[24], W[25]
        T1 = andc(c, w22, [w24, w25], w23, [w24])                                  # g1
        T0 = andc(c, q16, [q18, w23, w24, w25], q17, [w22, w25])                   # g0
        T2 = andc(c, q17, [w22, w23, T0, T1], q18, [w24])                          # g2
        andc(c, q17, [w22, w23, w24, T0, T2], q18, [w22, w23, w24], tgt=t19)       # g3 -> t19
        xform(c, t19, [w22, w25, T1], one=True)                                    # t19 = q19
        andc(c, q17, [w22, w23, T0, T1], q18, [w24], tgt=T2, erase=True)           # erase g2
        andc(c, q16, [q18, w23, w24, w25], q17, [w22, w25], tgt=T0, erase=True)   # erase g0
        # g4: left w22+w24+g3 = t19+1+w24+w25+g1, right w23+w24+w25+g1
        andc(c, t19, [w24, w25, T1], w23, [w24, w25, T1], tgt=t20, lone=True)
        xform(c, t20, [w23, w24], one=True)                                        # t20 = q20
        c.and_compute(t20, w24, t21)                                               # g5 = (w23+g4) w24 = q20 w24
        xform(c, t21, [w25], one=True)                                             # t21 = q21
        andc(c, w22, [w24, w25], w23, [w24], tgt=T1, erase=True)                  # erase g1
    return [(t19, rebuild), (t20, lambda c, q: None), (t21, lambda c, q: None)]


def st_end_ret(c, word, F, ang, feat, cfg):
    import invsqrt_fk_epN as EPN
    PROBE["op"], PROBE["word"], PROBE["live"] = len(c.ops), list(word), c.live
    if NAMES == ["J8"]:
        retire = joint_j8(word)
    else:
        retire = [(int(n[1:]), make_rebuild(load(n), word)) for n in NAMES]
    cap = min(int(cfg.get("ep_max", 9)), c.cap - c.live) if cfg.get("ep_auto") else int(cfg.get("ep_cap", 6))
    cap = int(os.environ.get("NWCAP", cap + len(retire)))
    PROBE["cap"] = cap
    const = cfg.get("ep_const")
    net = EPN.best_net(27, cap, profile="line" if const is not None else cfg.get("ep_profile", "generic"))
    ang = ang if ang is not None else [1 << (p + 27) for p in range(27)]
    if feat is None:
        feat = [((1 << (p + 54)) - (1 << (p + 27))) % (1 << 64) for p in range(10)]
    if cfg.get("ep_profile") == "clifford" or const is not None:
        feat = list(feat)[:7]
    if const is not None:
        return EPN.emit_ring(c, word, F, net, ang, feat, cap, const_token=int(const), const_one=True, retire=retire)
    return EPN.emit_ring(c, word, F, net, ang, feat, cap, retire=retire)


_orig = F8.st_end_f1
F8.st_end_f1 = st_end_ret


def cfg_of(p):
    cfg = dict(F8.CFG)
    cfg.update(json.load(open(p)))
    return cfg


def check(cfgp, nsamp=0):
    """Each retired wire against its program on every input at the start of END (no retire hook needed)."""
    import numpy as np
    c, regs = F8.build_full(cfg_of(cfgp))
    nsamp = int(nsamp)
    if nsamp:
        rng = np.random.default_rng(7)
        x = np.unique(np.concatenate([np.arange(1, 64), rng.integers(1, FG.MAX_X + 1, nsamp)])).astype(np.int64)
    else:
        x = np.arange(1, FG.MAX_X + 1, dtype=np.int64)
    sim = K.Sim(len(x), c.next + 2)
    sim.set_word(regs["X"], x)
    sim.run(c.ops, 0, PROBE["op"])
    word = PROBE["word"]
    print(f"check: {len(x)} inputs, END at op {PROBE['op']}, live {PROBE['live']}, cap {PROBE['cap']}")
    ok = True
    jobs = ([(json.load(open(os.path.expanduser("~/qcc/ws/invsqrt_fk/fcsup_box/J8_k6.prog.json"))),
              ["q19", "q20", "q21"])] if NAMES == ["J8"] else [(load(n), [n]) for n in NAMES])
    for prog, names in jobs:
        ins = [sim.q[wire_of(t, word)].astype(np.uint8) & 1 for t in prog["inputs"]]
        val = list(ins)
        for gt in prog["gates"]:
            def xs(items):
                v = np.zeros(len(x), dtype=np.uint8)
                for i in items:
                    v ^= 1 if i == -1 else val[i]
                return v
            val.append(xs(gt["left"]) & xs(gt["right"]))
        for n, out in zip(names, prog["outputs"]):
            pred = np.zeros(len(x), dtype=np.uint8)
            for i in out:
                pred ^= 1 if i == -1 else val[i]
            got = sim.q[int(n[1:])].astype(np.uint8) & 1
            bad = int(np.count_nonzero(pred != got))
            ok &= bad == 0
            print(f"  {n}: {prog['inputs']} {len(prog['gates'])} ANDs: mismatches {bad}")
    print("CHECK", "PASS" if ok else "FAIL")


def emit(out, angp, cfgp):
    import numpy as np
    a = json.load(open(angp))
    ang = [int(t) for t in a["ang"]]
    feat = [int(t) for t in a["feat"]] if a.get("feat") is not None else None
    c, regs = F8.build_full(cfg_of(cfgp), ang, feat)
    rng = np.random.default_rng(11)
    x = np.unique(np.concatenate([np.arange(1, 64), rng.integers(1, FG.MAX_X + 1, 40000)])).astype(np.int64)
    sim = K.Sim(len(x), c.next + 2)
    sim.set_word(regs["X"], x)
    sim.run(c.ops)
    xin = np.asarray(sim.get_word(regs["X"]), dtype=np.int64)
    others = [q for q in range(c.next) if q not in set(regs["X"])]
    dirty = [q for q in others if np.any(sim.q[q] & 1)]
    print(f"kit check on {len(x)} inputs: x restored {bool(np.all(xin == x))}; dirty wires {dirty[:8]}")
    os.makedirs(out, exist_ok=True)
    p = os.path.join(out, "invsqrt.kmx")
    with open(p + ".tmp", "w") as f:
        f.write(K.kmx_text(c, regs["X"]))
    os.replace(p + ".tmp", p)
    sys.path.insert(0, os.path.join(os.path.expanduser("~/qcc/data"), "kickmix"))
    import kickmix as km
    kc = km.Circuit.from_file(p, format="kmx")
    kc.to_file(os.path.join(out, "invsqrt.kmb.tmp"), format="kmb")
    os.replace(os.path.join(out, "invsqrt.kmb.tmp"), os.path.join(out, "invsqrt.kmb"))
    print(f"emit {out}: retired {NAMES} cap {PROBE['cap']}; kit toffolis {c.toffolis()} wires {c.next} peak "
          f"{c.peak}; END {dict(c.stage_counts()).get('END')}; kickmix qubits {kc.num_qubits}")


def is_function(target, cols):
    import numpy as np
    if not cols:
        return target.min() == target.max()
    key = np.zeros(len(target), dtype=np.uint64)
    for i, col in enumerate(cols):
        key |= col.astype(np.uint64) << np.uint64(i)
    order = np.argsort(key, kind="stable")
    k, t = key[order], target[order]
    return not np.any((k[1:] == k[:-1]) & (t[1:] != t[:-1]))


def support(cfgp, nsamp=0, maxsup=14):
    """Every non-word, non-flag wire live at the start of END: is it a function of the other live wires, and of
    which small set (greedy removal, other wires first, word bits last)."""
    import numpy as np
    c, regs = F8.build_full(cfg_of(cfgp))
    nsamp = int(nsamp)
    if nsamp:
        rng = np.random.default_rng(7)
        x = np.unique(np.concatenate([np.arange(1, 64), rng.integers(1, FG.MAX_X + 1, nsamp)])).astype(np.int64)
    else:
        x = np.arange(1, FG.MAX_X + 1, dtype=np.int64)
    sim = K.Sim(len(x), c.next + 2)
    sim.set_word(regs["X"], x)
    sim.run(c.ops, 0, PROBE["op"])
    live, word, F = PROBE["livew"], PROBE["word"], PROBE["F"]
    name = {w: "w%d" % i for i, w in enumerate(word)}
    name.update({w: "f%d" % i for i, w in enumerate(F)})
    vals = {w: sim.q[w].astype(np.uint8) & 1 for w in live}
    const = [w for w in live if vals[w].min() == vals[w].max()]
    cand = [w for w in live if w not in word and w not in F and w not in const]
    print(f"END live {len(live)}; constants {[name.get(w, 'q%d' % w) for w in const]}; candidates "
          f"{['q%d' % w for w in cand]}; inputs {len(x)}", flush=True)
    for w in cand:
        sup = [v for v in live if v != w and v not in const]
        if len(sup) > 63 or not is_function(vals[w], [vals[v] for v in sup]):
            print(f"  q{w}: not a function of the other live wires (information)", flush=True)
            continue
        best = None
        for order in ((lambda u: (u in word, u)), (lambda u: (u not in word, u)),
                      (lambda u: (u not in word, -u)), (lambda u: (u in word, -u))):
            s = list(sup)
            changed = True
            while changed:
                changed = False
                for v in sorted(s, key=order):
                    trial = [u for u in s if u != v]
                    if is_function(vals[w], [vals[u] for u in trial]):
                        s, changed = trial, True
                        break
            if best is None or len(s) < len(best):
                best = s
        sup = best
        tag = "SMALL" if len(sup) <= int(maxsup) else "large"
        print(f"  q{w}: function of {len(sup)} wires {[name.get(u, 'q%d' % u) for u in sup]} ({tag})", flush=True)


if __name__ == "__main__":
    if sys.argv[1] in ("check", "support"):
        def _probe(c, word, F, ang, feat, cfg):
            PROBE.update(op=len(c.ops), word=list(word), F=list(F), live=c.live, cap=c.cap - c.live,
                         livew=sorted(set(range(c.next)) - set(c.free)))
            return _orig(c, word, F, ang, feat, cfg)
        F8.st_end_f1 = _probe
        (check if sys.argv[1] == "check" else support)(*sys.argv[3:6])
    else:
        emit(*sys.argv[3:6])
