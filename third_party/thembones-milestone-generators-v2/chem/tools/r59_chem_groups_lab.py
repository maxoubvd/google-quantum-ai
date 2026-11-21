#!/usr/bin/env python3
"""r59_chem_groups_lab.py: designer R59-CHEM-GROUPS. Exact pricing experiments on dumped groups of a real loader.

Reads a dump of tools/chem_ld2_e19grp.py (builder C2's; imported read-only) and prices relation spaces with the
builder's own exact pricer (Group.price: kernels of every fan-out, supply, Lovasz rank over GF(2^20)).

  rand   random c-dim relation spaces inside W-perp: kernel profile and exact products (the statement's check).
  prof   the profile only (no Lovasz), many samples: how many fan-outs keep a plane at each c.
  hyp    a 6-group (W-perp of dimension a): every hyperplane of W-perp, the best one priced exactly.
  keep   bonus space + the extra relations that keep the planes of as many chosen fan-outs as linear algebra allows.

  python -I r59_chem_groups_lab.py rand g.json --g 0 --cs 13,14,15,16 --samples 3
"""
import argparse
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import chem_ld2_e19grp as eg  # noqa: E402  (read-only import of the builder's pricer)


def load(path, gi, seed=1):
    d = json.load(open(path))
    g = d["groups"][gi]
    G = eg.Group(g, seed=seed + 7)
    rel = eg.nullspace(G.wb, G.N)
    return d, g, G, rel


def cmd_rand(a):
    d, g, G, rel = load(a.dump, a.g, a.rseed)
    rng = random.Random(a.rseed)
    print(f"group {a.g}: b {G.b} a {G.a} N {G.N} dim_w {g['dim_w']} W-perp {len(rel)} engine {g['made']}", flush=True)
    for c in [int(x) for x in a.cs.split(",")]:
        if c > len(rel):
            continue
        for _ in range(a.samples):
            t0 = time.time()
            phis = eg.sub_space(rel, c, rng) if c < len(rel) else list(rel)
            res = G.price(phis, reps=a.reps)
            planes = sum(v for k, v in res["prof"].items() if k >= 2)
            print(f"c {c}: need {res['need']} miss {res['miss']} nu {res['nu']} products {res['prod']} "
                  f"fan-outs with a plane {planes} profile {dict(sorted(res['prof'].items()))} "
                  f"({time.time() - t0:.1f} s)", flush=True)


def cmd_prof(a):
    d, g, G, rel = load(a.dump, a.g, a.rseed)
    rng = random.Random(a.rseed)
    print(f"group {a.g}: b {G.b} a {G.a} N {G.N} dim_w {g['dim_w']} W-perp {len(rel)}", flush=True)
    for c in [int(x) for x in a.cs.split(",")]:
        tot = []
        for _ in range(a.samples):
            phis = eg.sub_space(rel, c, rng)
            fans = G.kernels(phis)
            tot.append(sum(1 for s, Y in fans if len(Y) >= 2))
        need = (G.N - c + 1) // 2
        print(f"c {c}: products at a perfect pairing {need}; fan-outs with a plane: mean {sum(tot) / len(tot):.1f} "
              f"min {min(tot)} max {max(tot)} over {len(tot)} random spaces", flush=True)


def cmd_hyp(a):
    """a group whose W-perp has dimension R = a (a 6-group: 15): every hyperplane Phi of W-perp (c = R - 1), the number
    of fan-outs that keep a plane (R_s meets Phi) and the exact price of the best one: the proof that c = a - 1 is
    out of reach there, measured."""
    d, g, G, rel = load(a.dump, a.g, a.rseed)
    R = len(rel)
    print(f"group {a.g}: b {G.b} a {G.a} N {G.N} dim_w {g['dim_w']} W-perp {R}", flush=True)
    rs = []                                           # per fan-out: basis of R_s in coefficients over rel
    for s in range(1, 1 << G.b):
        rows = []
        for t in range(G.a):                          # condition: sum_j x_j (rel_j restricted by s)[t] = 0
            v = 0
            for j, r in enumerate(rel):
                w = 0
                for o in eg.bits(s):
                    w ^= (r >> (o * G.a)) & G.full_a
                if (w >> t) & 1:
                    v |= 1 << j
            rows.append(v)
        rs.append(eg.nullspace(rows, R))
    dims = {}
    for b in rs:
        dims[len(b)] = dims.get(len(b), 0) + 1
    print(f"dim R_s over the {len(rs)} fan-outs: {dict(sorted(dims.items()))}", flush=True)
    best = (-1, 0)
    hist = {}
    for lam in range(1, 1 << R):
        n = 0
        for b in rs:
            if len(b) >= 2 or (len(b) == 1 and bin(lam & b[0]).count("1") % 2 == 0):
                n += 1
        hist[n] = hist.get(n, 0) + 1
        if n > best[0]:
            best = (n, lam)
    need = (G.N - (R - 1) + 1) // 2
    print(f"hyperplanes: {(1 << R) - 1}; fan-outs with a plane: best {best[0]} (needed {need}); top of the "
          f"histogram {sorted(hist.items())[-4:]}", flush=True)
    lam = best[1]
    piv = lam.bit_length() - 1
    phis = []
    for j in range(R):
        if j == piv:
            continue
        phis.append(rel[j] ^ (rel[piv] if (lam >> j) & 1 else 0))
    res = G.price(phis, reps=2)
    print(f"best hyperplane priced exactly: need {res['need']} miss {res['miss']} nu {res['nu']} products "
          f"{res['prod']} profile {dict(sorted(res['prof'].items()))}", flush=True)


def cmd_keep(a):
    """the best LINEAR construction of a space beyond the bonus: a random bonus space Phi0 (dimension a - 2: every
    fan-out has a plane Y0(s)), then the extra relations G of W-perp that annihilate Y0(s) (x) s for a chosen set of
    fan-outs (2 linear conditions each): as many fan-outs as leave --extra independent solutions outside Phi0.
    Priced exactly. (An 8-group at c 16: 14 fan-outs kept by construction, the rest by chance.)"""
    d, g, G, rel = load(a.dump, a.g, a.rseed)
    rng = random.Random(a.rseed)
    R, extra = len(rel), a.extra
    print(f"group {a.g}: b {G.b} a {G.a} N {G.N} dim_w {g['dim_w']} W-perp {R}; bonus {G.a - 2} + {extra} extra",
          flush=True)
    for _ in range(a.samples):
        phi0 = eg.sub_space(rel, G.a - 2, rng)
        fans = {s: Y for s, Y in G.kernels(phi0)}
        order = list(fans)
        rng.shuffle(order)
        keep, conds = [], []
        for s in order:
            new = conds + [G.feat(y, s) for y in fans[s][:2]]
            # relations of W-perp (as combinations of rel) orthogonal to every kept feature
            rows = []
            for f in new:
                v = 0
                for j, r in enumerate(rel):
                    if bin(r & f).count("1") & 1:
                        v |= 1 << j
                rows.append(v)
            sol = eg.nullspace(rows, R)
            if len(sol) >= G.a - 2 + extra:
                keep.append(s)
                conds = new
                last = sol
            if len(sol) == G.a - 2 + extra:
                break
        phis = []
        for cv in last[:G.a - 2 + extra]:
            v = 0
            for j in eg.bits(cv):
                v ^= rel[j]
            phis.append(v)
        res = G.price(phis, reps=2)
        planes = sum(v for k, v in res["prof"].items() if k >= 2)
        print(f"kept by construction {len(keep)} fan-outs; c {len(phis)}: fan-outs with a plane {planes}, need "
              f"{res['need']}, miss {res['miss']}, nu {res['nu']}, products {res['prod']}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("rand", "prof", "hyp", "keep"):
        p = sub.add_parser(name)
        p.add_argument("dump")
        p.add_argument("--g", type=int, default=0)
        p.add_argument("--cs", default="13,14,15,16")
        p.add_argument("--samples", type=int, default=3)
        p.add_argument("--reps", type=int, default=1)
        p.add_argument("--rseed", type=int, default=1)
        p.add_argument("--extra", type=int, default=3)
    a = ap.parse_args()
    {"rand": cmd_rand, "prof": cmd_prof, "hyp": cmd_hyp, "keep": cmd_keep}[a.cmd](a)


if __name__ == "__main__":
    main()
