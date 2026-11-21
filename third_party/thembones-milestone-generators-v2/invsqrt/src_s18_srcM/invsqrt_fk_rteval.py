"""F1 (third holder): price an arbitrary router network (stages = matchings on the 26 ring positions, applied in
order; class k fires the stages of its codeword) under the endpoint's lane model.

Checks: class k's fired product = shift by k (labels as in tools/invsqrt_fk_rtsat.py: class 0 fires nothing, slot of
base position p = ring position p). Controls: each stage's class vector is affine in (f0..f3, p) or not (then priced
as +2 per product of its ANF over the 11 classes, computed before and after a global stage). Walk: visit order local
search; holds: exact ILP of the interval packing under CAP with the feature helpers and p as fixed demands.
    python invsqrt_fk_rteval.py <network.json> [cap=5] [iters=4000]
"""
import json
import random
import sys

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fk_rtsearch as S  # noqa: E402  (AFF: affine class vectors)
import invsqrt_fk_router as R  # noqa: E402

N = 26


def perm_of(stage):
    p = list(range(N))
    for a, b in stage:
        p[a], p[b] = b, a
    return p


def check(stages, codes):
    bad = 0
    for k, code in enumerate(codes):
        pos = list(range(N))
        for s, bit in enumerate(code):
            if bit == "1":
                pm = perm_of(stages[s])
                pos = [pm[x] for x in pos]
        bad += sum(1 for x in range(N) if pos[x] != (x + k) % N)
    return bad


def windows(stages, order):
    vis = {x: i for i, x in enumerate(order)}
    perms = [perm_of(st) for st in stages]
    out = {}
    for s, st in enumerate(stages):
        for pr in st:
            cone = set(pr)
            for t in range(s + 1, len(stages)):
                cone |= {perms[t][x] for x in cone}
            ts = [vis[x] for x in cone]
            out[(s, tuple(pr))] = (min(ts), max(ts))
    return out


def fixed_demand(order, slots):
    vis = {x: i for i, x in enumerate(order)}
    d = [0] * N
    for j, h in R.FEATURE_HELPERS.items():
        d[vis[slots[j]]] += h
    fv = [vis[slots[j]] for j in range(9)]
    for t in range(min(fv), max(fv) + 1):
        d[t] += 1
    return d


def max_holds(win, base, cap):
    keys = list(win)
    A = np.zeros((N, len(keys)))
    for j, kk in enumerate(keys):
        a, b = win[kk]
        A[a:b + 1, j] = 1
    ub = np.array([cap - x for x in base], dtype=float)
    if (ub < 0).any():
        return -1
    res = milp(-np.ones(len(keys)), constraints=LinearConstraint(A, -np.inf, ub), integrality=np.ones(len(keys)),
               bounds=Bounds(0, 1), options=dict(time_limit=20))
    return int(round(-res.fun)) if res.x is not None else -1


def initial_order(stages):
    """chain the innermost stage's pairs through the next stage (greedy path), then the rest."""
    inner = [tuple(p) for p in stages[-1]]
    nxt = perm_of(stages[-2]) if len(stages) > 1 else list(range(N))
    order, seen = [], set()
    pair_of = {}
    for a, b in inner:
        pair_of[a] = (a, b)
        pair_of[b] = (a, b)
    todo = list(inner)
    while todo:
        pr = todo.pop(0)
        if pr[0] in seen:
            continue
        chain = [pr]
        while True:
            last = chain[-1][1]
            y = nxt[last]
            if y in pair_of and pair_of[y][0] not in seen and pair_of[y] not in chain:
                q = pair_of[y]
                chain.append(q if q[0] == y else (q[1], q[0]))
            else:
                break
        for a, b in chain:
            order += [a, b]
            seen |= {a, b}
    order += [x for x in range(N) if x not in seen]
    return order


def main(path, cap=5, iters=4000):
    cap, iters = int(cap), int(iters)
    net = json.load(open(path))
    stages, codes = [[tuple(p) for p in st] for st in net["stages"]], net["codes"]
    print("shift mismatches:", check(stages, codes))
    ctl_cost = 0
    for s in range(len(stages)):
        vec = tuple(int(codes[k][s]) for k in range(len(codes)))
        aff = vec in S.AFF
        if not aff:
            ctl_cost += 2
        print(f"  stage {s}: {len(stages[s])} pairs, control {''.join(map(str, vec))} affine {aff}")
    npairs = sum(len(st) for st in stages)
    slots = list(range(N))
    order = initial_order(stages)
    best_h = max_holds(windows(stages, order), fixed_demand(order, slots), cap)
    rnd = random.Random(1)
    cur, cur_h = order[:], best_h
    for it in range(iters):
        nxt = cur[:]
        i, j = rnd.randrange(N), rnd.randrange(N)
        if rnd.random() < 0.5:
            nxt[i], nxt[j] = nxt[j], nxt[i]
        else:
            a, b = min(i, j), max(i, j)
            nxt[a:b + 1] = nxt[a:b + 1][::-1]
        h = max_holds(windows(stages, nxt), fixed_demand(nxt, slots), cap)
        if h >= cur_h:
            cur, cur_h = nxt, h
            if h > best_h:
                best_h, order = h, nxt[:]
                print(f"  it {it}: holds {h} -> cost {2 * npairs - h + ctl_cost}", flush=True)
    print(f"pairs {npairs}, best holds {best_h}, control extra {ctl_cost}: cost {2 * npairs - best_h + ctl_cost}")
    print("order", order)


if __name__ == "__main__":
    main(*sys.argv[1:4])
