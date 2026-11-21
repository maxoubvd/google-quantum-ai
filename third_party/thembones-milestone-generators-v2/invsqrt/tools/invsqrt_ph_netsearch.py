"""invsqrt_ph: search for reflection swap networks of the normaliser (own code, workstream invsqrt_ph).

Idea credited to the best public circuit (manwithacat "HW-25" = gengsr + maoxin123; structure read in
research/ws_invsqrt_decode.md 3.1 / 3.2): a shift is two reflections q -> c - q, exchanging two zero wires is a
don't-care, so a few reflection layers with the right per-exponent subsets realise all 22 shifts. The sums, the
chains and the swap lists come from this search, not from the public file.

Model. Wires 0..21, exponent e, shift s = 21 - e. A network is an ordered list of layers; layer l is the set of
swaps (q, c_l - q). For an exponent every bit takes a path (a subset of the layers); no two bits may meet.
  forward  (NORM):   data bits x_0..x_(e-1) on wires 0..e-1 (the leading one is cleared by CNOTs from the
                     thermometer) must end on wires s..20. Every swap has its own free control over e.
  backward (UNNORM): data bits on wires s..20 must end on 0..e-1 and the leading one goes from wire 21 to e.
                     A swap control must be alpha(e >> 1) xor beta (e & 1) on the exponents where it matters
                     (CNOT sum of the half thermometer and the parity); otherwise it is charged PEN.
Bits of an exponent >= E_INT move as one interval (one common path); below that each bit has its own path. The
leading one always has its own path.
Cost = number of distinct swaps over all exponents (+ penalties).

usage: invsqrt_ph_netsearch.py fwd|bwd seconds "c1,c2,..;c1,c2,.." [extra layers 0..2]
       (candidate sum lists separated by ';'; for each, extra layers are tried at every slot)
"""
import itertools
import json
import sys
import time

NW = 22
PEN = 1
E_INT = int(__import__("os").environ.get("PH_EINT", "8"))


def path_end(sums, mask, q):
    for l, c in enumerate(sums):
        if (mask >> l) & 1:
            q = c - q
            if q < 0 or q > NW - 1:
                return None
    return q


def bit_paths(sums, q0, q1):
    return [m for m in range(1 << len(sums)) if path_end(sums, m, q0) == q1]


MODE = {"lead_in": False, "pen": False}


def options(sums, e, fwd, cap=60):
    """List of (pair masks per layer, occupancy before each layer) for exponent e, cheapest first."""
    s = 21 - e
    L = len(sums)
    if MODE["lead_in"]:                            # forward, leading one moved with the data (0..e -> s..21)
        data, lead = [(i, i + s) for i in range(e + 1)], []
    else:
        data = [(i, i + s) for i in range(e)] if fwd else [(i + s, i) for i in range(e)]
        lead = [] if fwd else [(21, e)]
    groups = []                                    # list of (bits, list of path masks)
    if data:
        if e >= E_INT:
            common = set.intersection(*[set(bit_paths(sums, a, b)) for a, b in data])
            if not common:
                return []
            groups.append((data, sorted(common)))
        else:
            for ab in data:
                p = bit_paths(sums, *ab)
                if not p:
                    return []
                groups.append(([ab], p))
    for ab in lead:
        p = bit_paths(sums, *ab)
        if not p:
            return []
        groups.append(([ab], p))
    if not groups:
        return [((0,) * L, (0,) * L)]
    n = 1
    for _, p in groups:
        n *= len(p)
    if n > 300000:
        return []
    out = {}
    for cand in itertools.product(*[p for _, p in groups]):
        pos, pm = [], []
        for (bits, _), m in zip(groups, cand):
            for a, _b in bits:
                pos.append(a)
                pm.append(m)
        masks, occ = [0] * L, [0] * L
        ok = True
        for l in range(L):
            c = sums[l]
            o = 0
            for q in pos:
                o |= 1 << q
            occ[l] = o
            for k in range(len(pos)):
                if (pm[k] >> l) & 1:
                    q = pos[k]
                    if c - q != q:
                        masks[l] |= 1 << min(q, c - q)
                    pos[k] = c - q
            if len(set(pos)) != len(pos):
                ok = False
                break
        if ok:
            out.setdefault(tuple(masks), tuple(occ))
    res = sorted(out.items(), key=lambda t: sum(bin(x).count("1") for x in t[0]))
    return res[:cap]


def family_ok(on, off):
    b0 = all(not ((2 * m in on and 2 * m + 1 in off) or (2 * m in off and 2 * m + 1 in on)) for m in range(11))
    b1 = all(not ((2 * m in on and 2 * m + 1 in on) or (2 * m in off and 2 * m + 1 in off)) for m in range(11))
    return b0 or b1


def cost_of(sums, choice, opts, fwd, detail=False):
    L = len(sums)
    tot, bad = 0, 0
    for l in range(L):
        c = sums[l]
        u = 0
        for e in range(22):
            u |= opts[e][choice[e]][0][l]
        tot += bin(u).count("1")
        if MODE["pen"]:
            for a in range(NW):
                if (u >> a) & 1:
                    on = {e for e in range(22) if (opts[e][choice[e]][0][l] >> a) & 1}
                    off = {e for e in range(22) if e not in on
                           and ((opts[e][choice[e]][1][l] >> a) & 1 or (opts[e][choice[e]][1][l] >> (c - a)) & 1)}
                    if not family_ok(on, off):
                        bad += 1
    return (tot + PEN * bad, tot, bad) if detail else tot + PEN * bad


def descend(sums, opts, fwd):
    choice = [0] * 22
    cur = cost_of(sums, choice, opts, fwd)
    improved = True
    while improved:
        improved = False
        for e in range(21, -1, -1):
            keep = choice[e]
            for k in range(len(opts[e])):
                if k == keep:
                    continue
                choice[e] = k
                c = cost_of(sums, choice, opts, fwd)
                if c < cur:
                    cur, keep, improved = c, k, True
            choice[e] = keep
    return cur, choice


def evaluate(sums, fwd):
    opts = []
    for e in range(22):
        o = options(sums, e, fwd)
        if not o:
            return None
        opts.append(o)
    cur, choice = descend(sums, opts, fwd)
    return cur, choice, opts


def describe(sums, choice, opts, fwd):
    layers = []
    for l, c in enumerate(sums):
        pairs = []
        u = 0
        for e in range(22):
            u |= opts[e][choice[e]][0][l]
        for a in range(NW):
            if (u >> a) & 1:
                on = [e for e in range(22) if (opts[e][choice[e]][0][l] >> a) & 1]
                off = [e for e in range(22) if e not in on
                       and ((opts[e][choice[e]][1][l] >> a) & 1 or (opts[e][choice[e]][1][l] >> (c - a)) & 1)]
                pairs.append([a, c - a, on, off])
        layers.append({"sum": c, "pairs": pairs})
    return {"dir": "fwd" if fwd else "bwd", "sums": list(sums), "layers": layers}


def main():
    fwd = sys.argv[1] in ("fwd", "fwdl")
    MODE["lead_in"] = sys.argv[1] == "fwdl"       # network for UNNORM, searched in its time-reversed form
    MODE["pen"] = sys.argv[1] in ("bwd", "fwdl")
    secs = float(sys.argv[2])
    n_extra = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    t0 = time.time()
    glob = None
    if sys.argv[3] == "enum":
        # all monotone 5-tuples whose 2- and 4-chains give the shifts 1..13, kept if every exponent >= E_INT
        # has a chain (cheap test), best first by the number of swaps of those exponents
        scored = []
        for tup in itertools.combinations(range(6, 36), 5):
            sums = list(tup) if fwd else list(tup[::-1])
            got = set()
            for i in range(5):
                for j in range(i + 1, 5):
                    d1 = sums[j] - sums[i]
                    got.add(d1)
                    for k in range(j + 1, 5):
                        for l in range(k + 1, 5):
                            got.add(d1 + sums[l] - sums[k])
            if not all((s if fwd else -s) in got for s in range(1, 14)):
                continue
            opts = []
            for e in range(22):
                o = options(sums, e, fwd, cap=8) if e >= E_INT else [((0,) * 5, (0,) * 5)]
                if not o:
                    break
                opts.append(o)
            if len(opts) == 22:
                scored.append((descend(sums, opts, fwd)[0], sums))
        scored.sort()
        print("enum: feasible", len(scored), scored[:12], "t=%.0f" % (time.time() - t0), flush=True)
        cands = [s for _, s in scored[:12]]
    else:
        cands = [[int(v) for v in part.split(",")] for part in sys.argv[3].split(";")]

    def consider(sums):
        nonlocal glob
        r = evaluate(sums, fwd)
        if r and (glob is None or r[0] < glob[0]):
            glob = (r[0], list(sums), r[1], r[2])
            print("best", cost_of(sums, r[1], r[2], fwd, True), sums, "t=%.0f" % (time.time() - t0), flush=True)
        return r

    for base in cands:
        consider(base)
    for _round in range(n_extra):
        start = list(glob[1]) if glob else None
        bases = ([start] if start else []) + (cands if _round == 0 else [])
        for base in bases:
            for slot in range(len(base) + 1):
                for c in range(1, 42):
                    if time.time() - t0 > secs:
                        break
                    consider(base[:slot] + [c] + base[slot:])
    if glob is None:
        print("nothing found")
        return
    d = describe(glob[1], glob[2], glob[3], fwd)
    d["cost"] = cost_of(glob[1], glob[2], glob[3], fwd, True)
    if MODE["lead_in"]:                            # time reversal: the same swaps in reverse order
        d["dir"] = "bwd"
        d["sums"] = d["sums"][::-1]
        d["layers"] = [{"sum": l["sum"], "pairs": l["pairs"][::-1]} for l in d["layers"][::-1]]
    print("NET " + json.dumps(d))


if __name__ == "__main__":
    main()
