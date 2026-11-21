"""F1 (fourth holder): endpoint router search on a ring of any size N (26 = the line, 27 = FC's family-C word).

Network: 5 stages, stage s = reflection x -> (c_s - x) mod N with an affine control over (f0..f3, p). Tokens i = 0..N-1
start at ring position i; class k (0..10) must send token i to ring position L((i + k) mod N), L affine (stride r = +-1,
offset free: the slot binding is a free relabelling). The two innermost stages (S3, S4) define the walk (alternating
S4 / S3 path through all N positions). Windows = closure cones (for reflection networks the class-aware model of
invsqrt_fk_rtclass.py gives the same windows). Holds packed by earliest deadline under the lane cap, with the wrap
features (slots L(0..nfeat-1), helpers per feature) and p as fixed demands. Cost = 2 x pairs - held.
    python invsqrt_fk_ringN.py N cap [top] [nfeat]
"""
import itertools
import sys

HELP = {0: 1, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 1, 7: 1, 8: 1, 9: 0}    # epkit's feature helpers; j = 9: CCZ
PFEATS = list(range(9))                                                  # features whose visits p must cover
# feature profiles (helpers per feature j, features that read p): see invsqrt_fk_epN.py
PROFILES = {
    "line": ({0: 1, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 1, 7: 1, 8: 1}, list(range(9))),
    "generic": ({0: 2, 1: 1, 2: 2, 3: 1, 4: 2, 5: 2, 6: 2, 7: 1, 8: 2, 9: 2}, [0, 1, 2, 4, 6, 9]),
    "clifford": ({0: 2, 1: 1, 2: 2, 3: 1, 4: 2, 5: 2, 6: 2, 7: 1, 8: 1, 9: 0}, [0, 1, 2, 4, 6]),
}


def use_profile(name):
    global HELP, PFEATS
    h, pf = PROFILES[name]
    HELP = dict(h)
    PFEATS = list(pf)


def flagvec(k):
    f = [(k >> i) & 1 for i in range(4)]
    return f + [f[1] & (f[2] ^ f[3])]


AFF = {}
for _mask in range(32):                       # stage controls affine in f0..f3 only (p is live only near features)
    _vals = []
    for _k in range(11):
        _v = flagvec(_k)
        _o = _mask & 1
        for _b in range(5):
            if _mask >> (_b + 1) & 1:
                _o ^= _v[_b]
        _vals.append(_o)
    AFF.setdefault(tuple(_vals), _mask)


def pairs(c, N):
    return [(x, (c - x) % N) for x in range(N) if x < (c - x) % N]


def windows(centres, order, N):
    vis = {x: i for i, x in enumerate(order)}
    out = {}
    for s, c in enumerate(centres):
        for pr in pairs(c, N):
            cone = set(pr)
            for t in range(s + 1, len(centres)):
                cc = centres[t]
                cone |= {(cc - x) % N for x in cone}
            ts = [vis[x] for x in cone]
            out[(s, pr)] = (min(ts), max(ts))
    return out


def walk_order(c3, c4, N):
    s3 = lambda x: (c3 - x) % N
    s4 = lambda x: (c4 - x) % N
    for start in sorted((x for x in range(N) if s3(x) == x), reverse=True):
        order, x, seen = [], start, set()
        while x not in seen:
            seen.add(x)
            order.append(x)
            y = s4(x)
            if y in seen:
                break
            seen.add(y)
            order.append(y)
            x = s3(y)
        if len(order) == N:
            return order
    return None


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


def evaluate(centres, N, cap, order, nfeat):
    out = []
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
    vis = {x: i for i, x in enumerate(order)}
    win = None
    for kind in ('r', 'f'):
        for r in (1, N - 1):
            for d0 in range(N):
                need = [((d0 + r * k) % N) if kind == 'r' else ((d0 - r * k) % N) for k in range(11)]
                if not all((kind, v) in rot for v in need):
                    continue
                for pick in itertools.product(*[rot[(kind, v)] for v in need]):
                    ctl = [tuple(pick[k][s] for k in range(11)) for s in range(5)]
                    if any(ctl[s] not in AFF for s in range(5)):
                        continue
                    L = [(r * p + d0) % N for p in range(N)] if kind == 'r' else [(-r * p + d0) % N for p in range(N)]
                    base = [0] * N
                    for j in range(nfeat):
                        base[vis[L[j]]] += HELP.get(j, 0)
                    fv = [vis[L[j]] for j in PFEATS if j < nfeat]
                    for t in range(min(fv), max(fv) + 1):
                        base[t] += 1
                    if win is None:
                        win = windows(centres, order, N)
                    h, held = pack(win, base, cap)
                    if h < 0:
                        continue
                    out.append((2 * len(win) - h, len(win), h, kind, r, d0, ctl, held, L))
    return out


def search(N, cap, top=10, nfeat=10):
    best = []
    for c4 in (N - 1,):
        for c3 in range(N):
            order = walk_order(c3, c4, N)
            if order is None:
                continue
            for cs in itertools.product(range(N), repeat=3):
                centres = (cs[0], cs[1], cs[2], c3, c4)
                if any(centres[i] == centres[i + 1] for i in range(4)):
                    continue
                for res in evaluate(centres, N, cap, order, nfeat):
                    best.append((res[0], centres, order, res))
    best.sort(key=lambda t: t[0])
    return best[:top], len(best)


if __name__ == "__main__":
    N, cap = int(sys.argv[1]), int(sys.argv[2])
    top = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    nfeat = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    if len(sys.argv) > 5:
        use_profile(sys.argv[5])
    res, n = search(N, cap, top, nfeat)
    print(f"ring {N}, cap {cap}, features {nfeat}: {n} networks with affine controls")
    for cost, centres, order, r in res:
        _, npairs, h, kind, rr, d0, ctl, held, L = r
        hs = [sum(1 for P in held if P[0] == s) for s in range(5)]
        print(f"cost {cost} pairs {npairs} held {h} {hs} centres {centres} {kind} r {rr} d0 {d0} "
              f"ctl {[''.join(map(str, c)) for c in ctl]} walk {order[:4]}..")
