"""lut_synth: synthesiser for the lookup blocks of our invsqrt generator (spec / program format: header of
tools/invsqrt_gen_pool.py; that file and its folder are read-only for this workstream).

  lut_synth_gen.py <block.spec.json> [--mode cheap|full] [--method auto|A|designA] [--seed S] [--out DIR] [...]

Default (cheap mode, method auto): spec in, one line "RESULT <block> ...: ANDs N, ... check: True" out, in seconds
(method A) to about a minute (design A on the address block); nothing is written unless --out is given.
--mode full uses the search caps of the other options (box, C kernel).

Method A (any block): one-sided products (fixed multiplier u, free co-factor v from one linear solve), any
combination of unfinished outputs may be finished next (target freedom), partial sums are copied into still
empty wires (free pool growth), and the last multiplier of every step is searched (Gray walk, C or numba kernel).
Method designA (blocks whose inputs are the eight address bits): the consumed-algebra design of
lut_synth_inv_a.py. Every program is replayed by an independent interpreter of the JSON format and compared with
the output truth tables on all 256 addresses before it is written to DIR/<block>.prog.json."""
import argparse
import json
import os
import random
import sys
import time

import numpy as np

sys.path.insert(0, "<HOME>/Quantum/tools")
from lut_synth_core import FULL, Proj, Span, bits, solve_cols, t_rem
from lut_synth_fast import KERNEL, fsearch, split2


class GM:
    """Wire machine with an arbitrary base: slot 0 = constant, 1..nb-1 = inputs, nb.. = anonymous work wires."""

    def __init__(self, base_vals, nw):
        self.nb = 1 + len(base_vals)
        self.nw = nw
        self.ns = self.nb + nw
        self.val = [FULL] + list(base_vals) + [0] * nw
        self.prog = []
        self.ands = 0

    def form(self, mask):
        v = 0
        s = 0
        val = self.val
        while mask:
            if mask & 1:
                v ^= val[s]
            mask >>= 1
            s += 1
        return v

    def AND(self, um, vm, tgt, tag=""):
        assert tgt >= self.nb and not ((um | vm) >> tgt) & 1, "target wire inside a control form"
        self.val[tgt] ^= self.form(um) & self.form(vm)
        self.prog.append(["and", um, vm, tgt, tag])
        self.ands += 1

    def CX(self, sm, tgt):
        assert tgt >= self.nb and not (sm >> tgt) & 1
        self.val[tgt] ^= self.form(sm)
        self.prog.append(["cx", sm, tgt])

    def express(self, f, slots):
        sp = Span()
        for s in slots:
            sp.add(self.val[s], 1 << s)
        r, c = sp.reduce(f)
        return c if r == 0 else None


def load_spec(path):
    S = json.load(open(path))
    S["_in"] = [(x["name"], int(x["tt"], 16)) for x in S["inputs"]]
    S["_out"] = [(x["name"], int(x["tt"], 16)) for x in S["outputs"]]
    return S


def independent_outputs(S):
    """Indices of outputs that are independent modulo span(1, inputs, earlier outputs); the rest are dependent."""
    sp = Span()
    sp.add(FULL)
    for n, v in S["_in"]:
        sp.add(v)
    ind, dep = [], []
    for i, (n, v) in enumerate(S["_out"]):
        (ind if sp.add(v)[0] else dep).append(i)
    return ind, dep


class SynthG:
    def __init__(self, S, seed=0, log=print):
        self.S = S
        self.rng = random.Random(seed)
        self.log = log
        self.ind, self.dep = independent_outputs(S)
        nw = len(self.ind) + S["max_scratch"]
        self.M = GM([v for n, v in S["_in"]], nw)
        self.T = [S["_out"][i][1] for i in self.ind]
        self.role = {}
        self.stats = []
        self.tries = 0

    def rem(self):
        return t_rem(self.M.val, self.T)

    def settle(self, tw, tau, Trem, t0=0):
        M = self.M
        t = t0
        for j, v in enumerate(Trem):
            if (tau >> j) & 1:
                t ^= v
        d = M.val[tw] ^ t
        if d:
            mask = M.express(d, [s for s in range(M.ns) if s != tw and M.val[s] != 0])
            assert mask is not None
            M.CX(mask, tw)
        assert M.val[tw] == t

    def step(self, budget, rmax=60):
        M, rng = self.M, self.rng
        Trem = self.rem()
        R = len(Trem)
        emp = [s for s in range(M.nb, M.ns) if M.val[s] == 0 and s not in self.role]
        if emp:
            tw = emp[0]
        else:
            fl = [s for s, r in self.role.items() if r == "fill"]
            assert fl, "no wire left for a target"
            tw = rng.choice(fl)
        pool_slots = [s for s in range(M.ns) if s != tw and M.val[s] != 0]
        pool = [(M.val[s], 1 << s) for s in pool_slots]
        pn = [s for s in pool_slots if s != 0]
        K = len(pn)
        d = K + 1
        J = M.val[tw]
        # an empty target wire gives a homogeneous problem (is some target combination in the span?); make it a
        # membership problem by forcing one target tstar into the combination
        tstar = 0
        Tall = Trem
        if J == 0:
            k = rng.randrange(R)
            tstar = Trem[k]
            Trem = Trem[:k] + Trem[k + 1:]
            J = tstar
        NS = M.ns
        t0 = time.time()
        if d * (d + 1) // 2 + R < 256 - 8:
            return False            # products of this pool cannot span enough: a helper is needed first
        r = max(1, int(np.ceil((256 - d - R - budget) / max(1, K))) - 1)
        rmax = min(rmax, r + 6 + 256 // max(1, K))
        while r <= rmax:
            us = []
            for i in range(r - 1):
                um = 0
                while um == 0 or um in us:
                    um = sum(1 << s for s in pn if rng.random() < 0.5)
                us.append(um)
            cols = []
            colmap = {}
            cid = NS
            for i, um in enumerate(us):
                utt = M.form(um)
                for s in pn:
                    cols.append((utt & M.val[s], 1 << cid))
                    colmap[(i, s)] = cid
                    cid += 1
            P = Proj(Trem + [v for v, _ in pool] + [v for v, _ in cols])
            if P.wd > 2:
                r += 1
                continue
            # measured need: codimension left by a random last multiplier (median of 3)
            needs = []
            for _ in range(3):
                um = sum(1 << s for s in pn if rng.random() < 0.5) or (1 << pn[0])
                utt = M.form(um)
                pk = P.pack([utt & M.val[s] for s in pn])
                sq = Span()
                for row in pk:
                    sq.add(sum(int(row[w]) << (64 * w) for w in range(P.wd)))
                needs.append(P.c - sq.rank())
            need = sorted(needs)[1]
            if need > budget:
                r += 1
                continue
            G = P.pack([M.val[a] & M.val[b] for a in pn for b in pn]).reshape(K, K, P.wd)
            Gl, Gh = split2(G)
            bl, bh = split2(P.pack([J]))
            P0 = Proj([v for v, _ in pool] + [v for v, _ in cols])      # quotient without the targets
            if P0.wd <= 2 and not tstar:
                G0l, G0h = split2(P0.pack([M.val[a] & M.val[b] for a in pn for b in pn]).reshape(K, K, P0.wd))
                c0l, c0h = split2(P0.pack([J]))
            else:
                G0l = G0h = None
                c0l = c0h = [0]
            res = None
            start = 1
            total = 1 << min(K, need + 5, budget + 3)
            while start < total:
                cnt = min(1 << 24, total - start)
                idx = fsearch(Gl, Gh, bl[0], bh[0], start, cnt, G0l, G0h, c0l[0], c0h[0])
                if idx < 0:
                    self.tries += cnt
                    start += cnt
                    continue
                self.tries += idx - start + 1
                x = idx ^ (idx >> 1)
                um = sum(1 << pn[k] for k in range(K) if (x >> k) & 1)
                if um and um not in us:
                    utt = M.form(um)
                    cols2 = list(cols)
                    cm2 = dict(colmap)
                    c2 = cid
                    for s in pn:
                        cols2.append((utt & M.val[s], 1 << c2))
                        cm2[(r - 1, s)] = c2
                        c2 += 1
                    sol = solve_cols(pool, cols2, J, Trem, allow_zero=bool(tstar))
                    if sol is not None:
                        res = (um, sol, cm2)
                        break
                start = idx + 1
            if res is None:
                r += 1
                continue
            um_last, (tau, comb), cm = res
            us.append(um_last)
            done = 0
            for i, um in enumerate(us):
                V = 0
                for s in pn:
                    if (comb >> cm[(i, s)]) & 1:
                        V |= 1 << s
                if V == 0:
                    continue
                M.AND(um, V, tw, "A")
                done += 1
                e = [s for s in range(M.nb, M.ns) if M.val[s] == 0 and s not in self.role and s != tw]
                if e and i < len(us) - 1:
                    M.CX(1 << tw, e[0])
                    self.role[e[0]] = "fill"
            self.settle(tw, tau, Trem, tstar)
            self.role[tw] = "tgt"
            R2 = len(self.rem())
            assert R2 == R - 1, f"step did not finish a target ({R} -> {R2})"
            self.stats.append((done, d, R, need))
            self.log(f"  A r={done} pool={d} R={R}->{R2} need=2^{need} ands={M.ands} t={time.time() - t0:.1f}s")
            return True
        return False

    def helper(self):
        """One random product of two pool forms into an empty wire (pool growth; the wire is consumed later)."""
        M, rng = self.M, self.rng
        emp = [s for s in range(M.nb, M.ns) if M.val[s] == 0 and s not in self.role]
        if not emp:
            return False
        pn = [s for s in range(1, M.ns) if M.val[s] != 0]
        while True:
            um = sum(1 << s for s in pn if rng.random() < 0.5)
            vm = sum(1 << s for s in pn if rng.random() < 0.5)
            p = M.form(um) & M.form(vm)
            if um and vm and p and M.express(p, [0] + pn) is None:
                break
        M.AND(um, vm, emp[0], "H")
        self.role[emp[0]] = "fill"
        return True

    def run(self, budget, helpers=0):
        for _ in range(helpers):
            self.helper()
        while self.rem():
            if not self.step(budget):
                if not self.helper():
                    return False
                self.log(f"  helper product, pool now {sum(1 for v in self.M.val if v)}")
        return True


def export(S, M, role, ind, dep, note=""):
    """Name the anonymous wires and append the free steps that put every output on its own wire.
    Returns the program dict of the generator's format."""
    names_in = [n for n, v in S["_in"]]
    outs = S["_out"]
    tw = sorted(s for s, r in role.items() if r == "tgt")
    assert len(tw) == len(ind), (len(tw), len(ind))
    T = [outs[i][1] for i in ind]
    live = [s for s in range(M.ns) if M.val[s] != 0]
    masks = [M.express(t, live) for t in T]
    assert all(m is not None for m in masks), "an output is not in the span of the final wires"
    idx = {s: j for j, s in enumerate(tw)}
    twmask = sum(1 << s for s in tw)
    rows = [sum(1 << idx[s] for s in bits(m & twmask)) for m in masks]
    rest = [m & ~twmask for m in masks]
    n = len(tw)
    ops = []
    A = list(rows)
    for c in range(n):
        if not (A[c] >> c) & 1:
            p = next(i for i in range(c + 1, n) if (A[i] >> c) & 1)
            A[c] ^= A[p]
            ops.append((c, p))
        for i in range(n):
            if i != c and (A[i] >> c) & 1:
                A[i] ^= A[c]
                ops.append((i, c))
    for (i, j) in reversed(ops):
        M.CX(1 << tw[j], tw[i])
    for o in range(n):
        if rest[o]:
            M.CX(rest[o], tw[o])
    assert all(M.val[tw[o]] == T[o] for o in range(n))
    name = {0: "1"}
    for k, nm in enumerate(names_in):
        name[1 + k] = nm
    for o, i in enumerate(ind):
        name[tw[o]] = outs[i][0]
    nb = getattr(M, "nb", 9)
    left = [s for s in range(nb, M.ns) if s not in tw]
    used = [s for s in left if M.val[s] != 0 or any(((st[1] | (st[2] if st[0] == "and" else 0)) >> s) & 1 or
                                                     st[-1 if st[0] == "cx" else 3] == s for st in M.prog)]
    assert len(used) <= S["max_scratch"], f"{len(used)} scratch wires used, {S['max_scratch']} allowed"
    for k, s in enumerate(used):
        name[s] = f"s{k}"
    steps = []
    for st in M.prog:
        if st[0] == "and":
            steps.append(["and", name[st[3]], [name[s] for s in bits(st[1])], [name[s] for s in bits(st[2])]])
        else:
            steps.append(["xor", name[st[2]], [name[s] for s in bits(st[1])]])
    # dependent outputs: linear in constant, inputs and the other outputs
    sp_slots = [0] + list(range(1, nb)) + tw
    for i in dep:
        m = M.express(outs[i][1], [s for s in sp_slots if M.val[s] != 0])
        assert m is not None
        if m:
            steps.append(["xor", outs[i][0], [name[s] for s in bits(m)]])
    return {"block": S["block"], "spec_sha": S["spec_sha"], "n_scratch": len(used), "steps": steps,
            "ands": sum(1 for s in steps if s[0] == "and"), "by": "lut_synth " + note}


def check(S, prog):
    """Independent interpreter of the program format; returns (ok, ands, message)."""
    val = {"1": FULL}
    for n, v in S["_in"]:
        val[n] = v
    for n, v in S["_out"]:
        val[n] = 0
    for k in range(S["max_scratch"]):
        val[f"s{k}"] = 0
    ands = 0

    def form(f):
        x = 0
        for nm in f:
            x ^= val[nm]
        return x

    known = set(val)
    for st in prog["steps"]:
        names = [st[1]] + list(st[2]) + (list(st[3]) if st[0] == "and" else [])
        if any(nm not in known for nm in names):
            return False, ands, f"unknown wire {[nm for nm in names if nm not in known][0]} (program is for another spec)"
        if st[0] == "and":
            if st[1] in st[2] or st[1] in st[3]:
                return False, ands, "target inside its own form"
            if st[1] not in val or st[1] == "1" or st[1] in dict(S["_in"]):
                return False, ands, f"bad target {st[1]}"
            val[st[1]] ^= form(st[2]) & form(st[3])
            ands += 1
        elif st[0] == "xor":
            if st[1] in st[2]:
                return False, ands, "target inside its own form"
            if st[1] not in val or st[1] == "1" or st[1] in dict(S["_in"]):
                return False, ands, f"bad target {st[1]}"
            val[st[1]] ^= form(st[2])
        else:
            return False, ands, f"unknown step {st[0]}"
    bad = [n for n, v in S["_out"] if val[n] != v]
    if bad:
        return False, ands, f"wrong outputs: {bad[:5]}"
    return True, ands, "all outputs equal their truth tables on 256 addresses"


def run_design_a(S, a):
    """Design A (consumed algebra) on a block whose inputs are the eight address bits."""
    from lut_synth_core import BASE
    from lut_synth_inv_a import SynthA, choose_merged, drive
    from lut_synth_smooth import find_smooth
    by_tt = {v: n for n, v in S["_in"]}
    assert all(BASE[1 + k] in by_tt for k in range(8)), "design A needs the eight address bits as inputs"
    S2 = dict(S)
    S2["_in"] = [(by_tt[BASE[1 + k]], BASE[1 + k]) for k in range(8)]
    ind, dep = independent_outputs(S2)
    T = [S["_out"][i][1] for i in ind]
    nw = len(ind) + S["max_scratch"]
    cache = os.path.expanduser(f"~/qcc/ws/lut_synth/smooth_{S['spec_sha'][:12]}.json")
    merged = None
    if a.split < 16:
        if os.path.exists(cache):
            sm = json.load(open(cache))
        else:
            sm = find_smooth(T, a.smooth_iters)
            json.dump(sm, open(cache, "w"))
        print("smooth elements (16-block weight, transitions):", [(o["bw"], o["tr"]) for o in sm[:10]], flush=True)
        merged = choose_merged(sm, a.split, a.nsmooth)
    A = SynthA(T, seed=a.seed, split=a.split, cubic=a.cubic, nw=nw, merged=merged)
    A.bias = a.need_bias
    A.slack = a.b_slack
    A.workers = a.workers
    A.smooth_fills = bool(a.smooth_fills)
    if not drive(A, a):
        return None, None
    A.M.nb = 9
    prog = export(S2, A.M, A.role, ind, dep, note=f"design A seed {a.seed} split {a.split}")
    return prog, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--budget", type=int, default=24)
    ap.add_argument("--helpers", type=int, default=0)
    ap.add_argument("--out", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--method", default="auto", help="A | designA | auto (designA when the inputs are a0..a7)")
    ap.add_argument("--smooth-iters", type=int, default=600)
    ap.add_argument("--mode", default="cheap",
                    help="cheap (default): small searches, seconds to a minute, for pricing a table; "
                         "full: the search caps given by the other options (minutes to hours, C kernel advised)")
    from lut_synth_inv_a import add_design_args
    add_design_args(ap)
    a = ap.parse_args()
    if a.mode == "cheap":
        a.budget = min(a.budget, 14)
        a.b_log = min(a.b_log, 18)
        a.b_att = min(a.b_att, 16)
        a.late_limit = min(a.late_limit, 14)
        a.extra_b = 0
        a.smooth_iters = min(a.smooth_iters, 150)
        a.workers = 1
    S = load_spec(a.spec)
    t0 = time.time()
    if a.method == "auto":
        from lut_synth_core import BASE as _B
        tts = set(v for n, v in S["_in"])
        a.method = "designA" if len(S["_in"]) == 8 and all(_B[1 + k] in tts for k in range(8)) else "A"
    if a.method == "designA":
        print(f"{S['block']}: design A, outputs {len(S['_out'])} scratch {S['max_scratch']} kernel {KERNEL}", flush=True)
        prog, A = run_design_a(S, a)
        if prog is None:
            print("FAILED", flush=True)
            return
        tries = A.stats["tries"]
        note = f"design A split {a.split}"
    else:
        G = SynthG(S, seed=a.seed)
        print(f"{S['block']}: inputs {len(S['_in'])} outputs {len(S['_out'])} (independent {len(G.ind)}) "
              f"scratch {S['max_scratch']} kernel {KERNEL} budget 2^{a.budget}", flush=True)
        if not G.run(a.budget, a.helpers):
            print("FAILED", flush=True)
            return
        prog = export(S, G.M, G.role, G.ind, G.dep, note=f"method A seed {a.seed} budget {a.budget}")
        tries = G.tries
        note = f"method A budget {a.budget} helpers {a.helpers}"
    good, ands, msg = check(S, prog)
    print(f"RESULT {S['block']} seed {a.seed} {note}: ANDs {ands}, scratch used {prog['n_scratch']}, "
          f"check: {good} ({msg}); tries 2^{np.log2(max(1, tries)):.1f}; time {time.time() - t0:.0f}s", flush=True)
    if good and a.out:
        os.makedirs(a.out, exist_ok=True)
        fn = os.path.join(a.out, f"{S['block']}{a.tag}.prog.json")
        best = None
        if os.path.exists(fn):
            old = json.load(open(fn))
            if old.get("spec_sha") == S["spec_sha"]:
                best = old.get("ands")
        if best is None or ands < best:
            json.dump(prog, open(fn, "w"))
            print("wrote", fn, flush=True)
        else:
            print(f"kept existing {fn} ({best} ANDs)", flush=True)


if __name__ == "__main__":
    main()
