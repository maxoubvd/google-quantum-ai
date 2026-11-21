#!/usr/bin/env python3
"""invsqrt_fc_s8spanemit.py -- FC (eighth holder): emit an F1-format AND program with NO recomputation on few wires,
the wires holding XOR combinations of gate values (CNOT basis changes are free). Search = the span model of
tools/invsqrt_fc_s8span.py (COMPUTE / INTO / ERASE on the span S of the wires' gate parts); emission replays it:
  COMPUTE v : fresh temp t, factors formed on private carriers, ANDC (1 Toffoli)
  INTO v    : CNOTs make one wire t the only one outside the hyperplane S', CCX into t (1 Toffoli)
  ERASE u   : CNOTs make one wire t hold exactly u, ANDE (HMR + CZ if m: 0 Toffolis), t freed
  end       : CNOT basis change so that outs[k] holds output k exactly, temporaries back to 0 and freed.
API: emit_span(c, prog, ins, outs, budget) -> AND count, or None when no schedule fits len(outs) + budget wires.
  python invsqrt_fc_s8spanemit.py PROG.json N_INPUTS BUDGET      self-test on a fresh circuit + all-rows check"""
import os
import sys
import time

sys.setrecursionlimit(10000)


def _part(lst, n):
    v = 0
    for x in lst:
        if x >= n:
            v ^= 1 << (x - n)
    return v


def _red(B, v):
    for b in B:
        if v ^ b < v:
            v ^= b
    return v


def _basis(vecs):
    B = []
    for v in vecs:
        v = _red(B, v)
        if v:
            B.append(v)
            B.sort(reverse=True)
    for i in range(len(B)):
        p = B[i].bit_length() - 1
        for j in range(len(B)):
            if j != i and (B[j] >> p) & 1:
                B[j] ^= B[i]
    return tuple(sorted(B, reverse=True))


def _inspan(B, v):
    return _red(list(B), v) == 0


def _coords(B, v):
    piv = [b.bit_length() - 1 for b in B]
    c = 0
    for i, b in enumerate(B):
        if (v >> piv[i]) & 1:
            v ^= b
            c |= 1 << i
    return c if v == 0 else None


def _hyper(B, U, avoid=None):
    d = len(B)
    Uc = [_coords(B, u) for u in U if u]
    ac = _coords(B, avoid) if avoid is not None else None
    res = []
    for f in range(1, 1 << d):
        if any(bin(f & u).count("1") % 2 for u in Uc):
            continue
        if ac is not None and bin(f & ac).count("1") % 2 == 0:
            continue
        ker, k0 = [], None
        for i in range(d):
            if (f >> i) & 1:
                if k0 is None:
                    k0 = i
                else:
                    ker.append(B[i] ^ B[k0])
            else:
                ker.append(B[i])
        res.append((_basis(ker), B[k0]))
    return res


def search(prog, n, W, tlim=120):
    G = len(prog["gates"])
    facs = [(_part(g["left"], n), _part(g["right"], n)) for g in prog["gates"]]
    target = _basis([_part(o, n) for o in prog["outputs"]])
    FULL = (1 << G) - 1
    seen = set()
    t0 = time.time()

    def within(Bv, mask):
        out = mask ^ FULL
        piv, res = [], []
        for b in Bv:
            o = b & out
            for po, pb in piv:
                if o ^ po < o:
                    o ^= po
                    b ^= pb
            if o:
                piv.append((o, b))
                piv.sort(reverse=True)
            elif b:
                res.append(b)
        return res

    def feasible(done, B):
        for t in within(target, done):
            if not _inspan(B, t):
                return False
        for v in range(G):
            if not (done >> v) & 1:
                for f in facs[v]:
                    if f & ~done == 0 and not _inspan(B, f):
                        return False
        return True

    def dfs(done, B):
        if time.time() - t0 > tlim:
            raise TimeoutError
        key = (done, B)
        if key in seen:
            return None
        seen.add(key)
        if not feasible(done, B):
            return None
        if done == FULL and B == target:
            return []
        moves = []
        for v in range(G):
            if (done >> v) & 1:
                continue
            L, R = facs[v]
            if not (_inspan(B, L) and _inspan(B, R)):
                continue
            if len(B) < W:
                moves.append((0, len(B) + 1, ("C", v, None, None), done | (1 << v), _basis(list(B) + [1 << v])))
            if B:
                for Sp, w in _hyper(B, [L, R]):
                    moves.append((2, len(B), ("I", v, Sp, w), done | (1 << v), _basis(list(Sp) + [w ^ (1 << v)])))
        for u in range(G):
            if not (done >> u) & 1 or not _inspan(B, 1 << u):
                continue
            L, R = facs[u]
            if not (_inspan(B, L) and _inspan(B, R)):
                continue
            for Sp, _w in _hyper(B, [L, R], avoid=1 << u):
                moves.append((1, len(Sp), ("E", u, Sp, None), done, Sp))
        moves.sort(key=lambda m: (m[0], m[1]))
        for _o, _d, mv, nd, nb in moves:
            r = dfs(nd, nb)
            if r is not None:
                return [mv] + r
        return None

    try:
        return dfs(0, ())
    except TimeoutError:
        return "timeout"


def emit_span(c, prog, ins, outs, budget, tlim=120):
    n = len(ins)
    G = len(prog["gates"])
    path = search(prog, n, len(outs) + budget, tlim)
    if path is None or path == "timeout":
        return None
    content = {}
    for i, w in enumerate(ins):
        content[w] = 1 << i
    for w in outs:
        content[w] = 0
    temps = []
    work = list(outs)                          # non-input wires in play (outputs + temporaries)
    k = 0

    def gpart(v):
        return v >> n

    def express(vec, exclude=()):
        """wires (inputs and work wires, not in exclude) whose contents XOR to vec; None if impossible."""
        rows = [(content[w], w) for w in list(ins) + work if w not in exclude and content[w]]
        piv = []
        for v, w in rows:
            s = {w}
            for pv, ps in piv:
                if v ^ pv < v:
                    v ^= pv
                    s = s ^ ps
            if v:
                piv.append((v, s))
                piv.sort(key=lambda t: -t[0])
        s = set()
        for pv, ps in piv:
            if vec ^ pv < vec:
                vec ^= pv
                s ^= ps
        return s if vec == 0 else None

    def form_pair(L, R, exclude):
        sL, sR = express(L, exclude), express(R, exclude)
        assert sL is not None and sR is not None, "factor not formable"
        aL, aR = sL - sR, sR - sL
        mode = "plain"
        if not aL or not aR:                   # L R = L (L ^ R) ^ L: the pair (L, L ^ R) has private carriers
            if not aR:
                L, R = R, L
                sL, sR = sR, sL
            sR = sR ^ sL
            R = L ^ R
            aL, aR = sL - sR, sR - sL
            mode = "lin"
            assert aL and aR, "no private carriers"
        cL, cR = min(aL), min(aR)
        for w in sorted(sL - {cL}):
            c.cx(w, cL)
        for w in sorted(sR - {cR}):
            c.cx(w, cR)
        return (cL, sL, cR, sR, mode, L)

    def unform(fp):
        cL, sL, cR, sR, mode, L = fp
        for w in sorted(sR - {cR}):
            c.cx(w, cR)
        for w in sorted(sL - {cL}):
            c.cx(w, cL)

    def fvec(lst):
        return _fv(lst, n)

    for kind, v, Sp, w in path:
        if kind in ("I", "E"):
            # the functional f with kernel Sp on S: make one work wire the only one with f = 1 (outside Sp)
            ones = [x for x in work if gpart(content[x]) and not _inspan(Sp, gpart(content[x]))]
            assert ones, "no wire outside the hyperplane"
            t = ones[0]
            for u in ones[1:]:
                c.cx(t, u)
                content[u] ^= content[t]
        if kind == "C":
            gt = prog["gates"][v]
            L, R = fvec(gt["left"]), fvec(gt["right"])
            fp = form_pair(L, R, ())
            zero = [x for x in outs if content[x] == 0 and x in work]
            if zero:                           # an output wire still (or again) at 0 holds the new gate
                tw = zero[0]
            else:
                tw = c.alloc(1)[0]
                temps.append(tw)
                work.append(tw)
            content[tw] = 0
            c.and_compute(fp[0], fp[2], tw)
            if fp[4] == "lin":
                c.cx(fp[0], tw)                # t = L (L ^ R) ^ L
            unform(fp)
            content[tw] = 1 << (n + v)
            k += 1
        elif kind == "I":
            gt = prog["gates"][v]
            L, R = fvec(gt["left"]), fvec(gt["right"])
            fp = form_pair(L, R, (t,))
            c.ccx(fp[0], fp[2], t)
            if fp[4] == "lin":
                c.cx(fp[0], t)
            unform(fp)
            content[t] ^= 1 << (n + v)
            k += 1
        else:                                  # ERASE v: wire t holds v + s', s' in Sp + inputs
            rest = content[t] ^ (1 << (n + v))
            s = express(rest, (t,))
            assert s is not None
            for x in sorted(s):
                c.cx(x, t)
            content[t] = 1 << (n + v)
            gt = prog["gates"][v]
            L, R = fvec(gt["left"]), fvec(gt["right"])
            fp = form_pair(L, R, (t,))
            if fp[4] == "lin":                 # t = L (L ^ R) ^ L: take the L part out first
                c.cx(fp[0], t)
            c.and_erase(fp[0], fp[2], t, free=False)
            unform(fp)
            content[t] = 0
            work.remove(t)
            if t in temps:
                temps.remove(t)
                c.free_([t])
            else:                              # an output wire emptied: it stays allocated (it is an output)
                work.append(t)
    # end: outputs[k] = output k exactly, temporaries zero
    tg = [fvec(o) for o in prog["outputs"]]
    for j, wo in enumerate(outs):
        s = express(tg[j], ())
        assert s is not None, "output not in the span"
        s = set(s)
        if wo not in s:
            cand = [x for x in s if x in work and x not in outs[:j]]
            assert cand, "no wire to move into the output"
            u = cand[0]
            c.cx(u, wo); c.cx(wo, u); c.cx(u, wo)          # swap contents
            content[u], content[wo] = content[wo], content[u]
            s = set(express(tg[j], ()))
        for x in sorted(s - {wo}):
            c.cx(x, wo)
            content[wo] ^= content[x]
        assert content[wo] == tg[j]
        if -1 in prog["outputs"][j]:
            c.x(wo)
    for tw in list(temps):
        s = express(content[tw], (tw,))
        assert s is not None, "temporary not clearable"
        for x in sorted(s):
            c.cx(x, tw)
        content[tw] = 0
        c.free_([tw])
        temps.remove(tw)
    return k


def _fv(lst, n):
    v = 0
    for x in lst:
        if x >= 0:
            v ^= 1 << x                        # signal index x: inputs 0..n-1, gates n.. (same bit layout)
    return v


if __name__ == "__main__":
    import json
    import numpy as np
    sys.path.insert(0, os.path.expanduser(os.environ.get("FC8SRC", "~/qcc/ws/invsqrt_fc/s6/src_m3l")))
    import invsqrt_fc_kit as K
    prog = json.load(open(sys.argv[1]))
    prog = {"gates": prog["gates"], "outputs": prog["outputs"]}
    n, budget = int(sys.argv[2]), int(sys.argv[3])
    c = K.Circ(n_fixed=n)
    c.cap = n + len(prog["outputs"]) + budget
    ins = list(range(n))
    outs = c.alloc(len(prog["outputs"]))
    live0 = c.live
    t0 = time.time()
    kk = emit_span(c, prog, ins, outs, budget)
    print("emit_span: %s ANDs, %d ops, peak %d (cap %d), live back %s (%.1f s)" % (
        kk, len(c.ops), c.peak, c.cap, c.live == live0, time.time() - t0))
    if kk is None:
        sys.exit(1)
    # classical check on every row: wire values as integers
    rows = 1 << n
    bad = 0
    for r in range(rows):
        val = {w: 0 for w in range(c.next + 2)}
        for i in range(n):
            val[i] = (r >> i) & 1
        for op in c.ops:
            nm = op[0]
            if nm == "CX":
                if op[3] is None:
                    val[op[2]] ^= val[op[1]]
            elif nm == "X":
                val[op[1]] ^= 1
            elif nm in ("ANDC", "CCX"):
                val[op[3]] ^= val[op[1]] & val[op[2]]
            elif nm == "ANDE":
                assert val[op[3]] == val[op[1]] & val[op[2]], ("ANDE on a wrong value", r, op)
                val[op[3]] = 0
        sig = [(r >> i) & 1 for i in range(n)]
        for g in prog["gates"]:
            a = 0
            for x in g["left"]:
                a ^= sig[x]
            b = 0
            for x in g["right"]:
                b ^= sig[x]
            sig.append(a & b)
        for j, wo in enumerate(outs):
            e = 0
            for x in prog["outputs"][j]:
                e ^= 1 if x == -1 else sig[x]
            bad += val[wo] != e
        bad += any(val[i] != ((r >> i) & 1) for i in range(n))
        bad += any(val[w] for w in range(c.next) if w not in outs and w >= n)
    print("all %d rows: %d mismatches" % (rows, bad))
