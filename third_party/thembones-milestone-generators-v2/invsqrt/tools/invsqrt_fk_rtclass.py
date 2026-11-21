"""F1 (fourth holder): class-aware window model of the invsqrt endpoint router.

The third holder's model (tools/invsqrt_fk_router.py) keeps a pair of stage s open while ANY slot of its closure cone
(closure of its two positions under ALL later stages, whatever the class) is visited. That is conservative. Here:
  - for each class k, each token's trajectory is the list of pairs that FIRE in class k and move it; the pair is
    needed at the visit of the token's final slot (need(P) = final slots of its two tokens over its firing classes);
  - nesting is enforced only along trajectories: if a token's class-k trajectory passes P (stage s) then Q (stage
    s' > s), W(P) must contain W(Q) (P opened before Q and closed after Q, so that held ANDs and CSWAPs restore);
  - W(P) = smallest interval holding span(need(P)) and every W(Q) nested in it.
Holds are then packed under the lane cap (earliest-deadline greedy and, if scipy is present, an exact ILP).
Cost = 2 x pairs - held. Our own code; the network and walk conventions are those of invsqrt_fk_router.py.
    python invsqrt_fk_rtclass.py public [cap]          price the public network in both models
    python invsqrt_fk_rtclass.py search [cap] [top]    reflection family of invsqrt_fk_rtsearch.py, class-aware
"""
import itertools
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fk_router as R  # noqa: E402

N = 26


def trajectories(centres, ctl):
    """traj[k][token] = list of (stage, pair) fired on the token's path in class k; final[k][token] = slot."""
    traj, final = [], []
    for k in range(len(ctl[0])):
        pos = list(range(N))
        paths = [[] for _ in range(N)]
        for s, c in enumerate(centres):
            if not ctl[s][k]:
                continue
            for t in range(N):
                x = pos[t]
                y = (c - x) % N
                if y != x:
                    paths[t].append((s, (min(x, y), max(x, y))))
                    pos[t] = y
        traj.append(paths)
        final.append(pos)
    return traj, final


def class_windows(centres, ctl, order):
    vis = {x: i for i, x in enumerate(order)}
    traj, final = trajectories(centres, ctl)
    need = {}
    inner = {}
    for k in range(len(traj)):
        for t in range(N):
            path = traj[k][t]
            v = vis[final[k][t]]
            for i, P in enumerate(path):
                need.setdefault(P, set()).add(v)
                inner.setdefault(P, set()).update(path[i + 1:])
    win = {}
    for P in sorted(need, key=lambda p: -p[0]):          # innermost stages first
        lo, hi = min(need[P]), max(need[P])
        for Q in inner.get(P, ()):
            a, b = win[Q]
            lo, hi = min(lo, a), max(hi, b)
        win[P] = (lo, hi)
    return win


def pack(win, base, cap):
    free = [cap - d for d in base]
    if min(free) < 0:
        return -1, []
    held = []
    for key in sorted(win, key=lambda kk: (win[kk][1], win[kk][0])):
        a, b = win[key]
        if all(free[t] > 0 for t in range(a, b + 1)):
            for t in range(a, b + 1):
                free[t] -= 1
            held.append(key)
    return len(held), held


def pack_ilp(win, base, cap):
    try:
        import numpy as np
        from scipy.optimize import milp, LinearConstraint, Bounds
    except Exception:
        return None
    keys = list(win)
    T = len(base)
    A = np.zeros((T, len(keys)))
    for j, key in enumerate(keys):
        a, b = win[key]
        A[a:b + 1, j] = 1
    ub = np.array([cap - d for d in base], dtype=float)
    res = milp(c=-np.ones(len(keys)), constraints=LinearConstraint(A, -np.inf, ub), integrality=np.ones(len(keys)),
               bounds=Bounds(0, 1))
    if not res.success:
        return None
    return int(round(-res.fun))


def public_base(order):
    return R.fixed_demand(order)


def cmd_public(cap=5):
    cap = int(cap)
    ctl = [[R.affine(R.PUB_CONTROLS[s], k) for k in range(11)] for s in range(5)]
    assert R.check_shifts(R.PUB_CENTRES, ctl) == 0
    order = R.public_walk()
    base = public_base(order)
    wc = R.windows(R.PUB_CENTRES, order)
    wa = class_windows(R.PUB_CENTRES, ctl, order)
    npairs = len(wc)
    for name, w in (("closure", wc), ("class-aware", wa)):
        lens = {s: sorted(w[P][1] - w[P][0] + 1 for P in w if P[0] == s) for s in range(5)}
        h, held = pack(w, base, cap)
        hi = pack_ilp(w, base, cap)
        hs = [sum(1 for P in held if P[0] == s) for s in range(5)]
        print(f"{name}: pairs {npairs} held(greedy) {h} {hs} ilp {hi} cost {2 * npairs - (hi if hi is not None else h)}")
        for s in range(5):
            print(f"   S{s} window lengths {lens[s]}")


AFF = None


def cmd_search(cap=5, top=12):
    import invsqrt_fk_rtsearch as S
    cap, top = int(cap), int(top)
    order = S.walk_order(0)
    vis = {x: i for i, x in enumerate(order)}
    results = []
    nnet = 0
    for cs in itertools.product(range(N), repeat=3):
        centres = (cs[0], cs[1], cs[2], 0, 25)
        if any(centres[i] == centres[i + 1] for i in range(4)):
            continue
        subsets = []
        for bits in itertools.product((0, 1), repeat=5):
            perm = list(range(N))
            for s in range(5):
                if bits[s]:
                    perm = [(centres[s] - x) % N for x in perm]
            subsets.append((bits, perm))
        rot = {}
        for bits, perm in subsets:
            d = perm[0] % N
            if all(perm[x] == (x + d) % N for x in range(N)):
                rot.setdefault(('r', d), []).append(bits)
            c = perm[0]
            if all(perm[x] == (c - x) % N for x in range(N)):
                rot.setdefault(('f', c), []).append(bits)
        for kind in ('r', 'f'):
            for r in (1, N - 1):
                for d0 in range(N):
                    need = [((d0 + r * k) % N) if kind == 'r' else ((d0 - r * k) % N) for k in range(11)]
                    if not all((kind, v) in rot for v in need):
                        continue
                    for pick in itertools.product(*[rot[(kind, v)] for v in need]):
                        ctl = [tuple(pick[k][s] for k in range(11)) for s in range(5)]
                        if any(ctl[s] not in S.AFF for s in range(5)):
                            continue
                        L = [(r * p + d0) % N for p in range(N)] if kind == 'r' else [(-r * p + d0) % N for p in range(N)]
                        base = [0] * N
                        for j, h in R.FEATURE_HELPERS.items():
                            base[vis[L[j]]] += h
                        fv = [vis[L[j]] for j in range(9)]
                        for t in range(min(fv), max(fv) + 1):
                            base[t] += 1
                        win = class_windows(centres, ctl, order)
                        h, held = pack(win, base, cap)
                        if h < 0:
                            continue
                        nnet += 1
                        results.append((2 * len(win) - h, len(win), h, centres, kind, r, d0, ctl, held))
    results.sort(key=lambda t: t[0])
    print(f"class-aware search, cap {cap}: {nnet} networks with affine controls")
    for res in results[:top]:
        cost, npairs, h, centres, kind, r, d0, ctl, held = res
        hs = [sum(1 for P in held if P[0] == s) for s in range(5)]
        print(f"cost {cost} pairs {npairs} held {h} {hs} centres {centres} {kind} r {r} d0 {d0} "
              f"ctl {[''.join(map(str, c)) for c in ctl]}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "public"
    if cmd == "public":
        cmd_public(*sys.argv[2:3])
    elif cmd == "search":
        cmd_search(*sys.argv[2:4])
