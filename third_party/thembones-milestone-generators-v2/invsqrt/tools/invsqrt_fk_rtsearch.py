"""F1 (third holder): search controlled-reflection networks for the invsqrt endpoint router under the lane model of
tools/invsqrt_fk_router.py.

Network: 5 stages, stage s = reflection x -> (c_s - x) mod 26; the two innermost stages form the walk (S3 centre C3,
S4 centre 25; the walk visits S4 pairs in the S3/S4 path order). Classes k = 0..10 must realise
token i -> ring position L((i + k) mod 26) with L(p) = (r p + l0) mod 26 and the word on positions W(i) = (r i + w0).
For every centre triple (S0, S1, S2) and stride r in {1, 25}: the class -> subset assignment is forced (shift values),
the five controls must be affine in (f0..f3, p) (else the network is skipped, or priced with ctlcost), the feature
slots are L(0..8), and the maximum number of held pairs under the lane cap is found by the earliest-deadline greedy
(optimal for unit intervals under a capacity profile). Cost = 2 x pairs - held.
    python invsqrt_fk_rtsearch.py [cap=5] [top=15] [c3=0]
"""
import itertools
import sys

sys.path.insert(0, "<HOME>/Quantum/tools")
import invsqrt_fk_router as R  # noqa: E402

N = 26


def flagvec(k):
    f = [(k >> i) & 1 for i in range(4)]
    return f + [f[1] & (f[2] ^ f[3])]


AFF = {}
for mask in range(64):           # bit 0 = constant, bits 1..5 = f0, f1, f2, f3, p
    vals = []
    for k in range(11):
        v = flagvec(k)
        o = mask & 1
        for b in range(5):
            if mask >> (b + 1) & 1:
                o ^= v[b]
        vals.append(o)
    AFF.setdefault(tuple(vals), mask)


def walk_order(c3, c4=25):
    """alternating S4 / S3 path through all 26 positions, started at an S3 fixed point."""
    s3 = lambda x: (c3 - x) % N
    s4 = lambda x: (c4 - x) % N
    fixed3 = [x for x in range(N) if s3(x) == x]
    best = None
    for start in sorted(fixed3, key=lambda x: -x if START_HIGH[0] else x):
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
            best = order
            break
    return best


START_HIGH = [True]
ANY_CTL = [False]
ALLOW_REP = [False]
STRIDES = [(1, N - 1)]


def max_holds(win, cand, base, cap):
    free = [cap - d for d in base]
    if min(free) < 0:
        return -1, []
    held = []
    for key in sorted(cand, key=lambda kk: (win[kk][1], win[kk][0])):
        a, b = win[key]
        if all(free[t] > 0 for t in range(a, b + 1)):
            for t in range(a, b + 1):
                free[t] -= 1
            held.append(key)
    return len(held), held


def evaluate(centres, cap, c3order):
    """try strides r = +-1 and all offsets; returns list of (cost, info)."""
    out = []
    subsets = []
    for bits in itertools.product((0, 1), repeat=5):
        perm = list(range(N))
        for s in range(5):
            if bits[s]:
                perm = [(centres[s] - x) % N for x in perm]   # perm[pos0] = final position of the token at pos0
        subsets.append((bits, perm))
    # a rotation by d: perm[x] = x + d; a reflection: perm[x] = c - x
    rot = {}
    for bits, perm in subsets:
        d = (perm[0] - 0) % N
        if all(perm[x] == (x + d) % N for x in range(N)):
            rot.setdefault(('r', d), []).append(bits)
        c = perm[0]
        if all(perm[x] == (c - x) % N for x in range(N)):
            rot.setdefault(('f', c), []).append(bits)
    for kind in ('r', 'f'):
        for r in STRIDES[0]:
            for d0 in range(N):
                # class k needs element ('r', d0 + r k) for rotations; reflections ('f', d0 - r k)... with word stride r
                need = [((d0 + r * k) % N) if kind == 'r' else ((d0 - r * k) % N) for k in range(11)]
                if not all((kind, v) in rot for v in need):
                    continue
                choices = [rot[(kind, v)] for v in need]
                for pick in itertools.product(*choices):
                    ctl = [tuple(pick[k][s] for k in range(11)) for s in range(5)]
                    if ANY_CTL[0] is False and any(ctl[s] not in AFF for s in range(5)):
                        continue
                    # slot layout: token i (at word position W(i)) goes to L((i+k)%26).
                    # rotations, stride r: positions x = r i + w0, final x + d0 + r k = r (i + k) + w0 + d0 -> L(p) = r p + w0 + d0
                    # reflections: final c - x = d0 - r k - r i - w0 = -r (i + k) + d0 - w0 -> L(p) = -r p + d0 - w0
                    for w0 in (0,):
                        if kind == 'r':
                            L = [(r * p + w0 + d0) % N for p in range(N)]
                        else:
                            L = [(-r * p + d0 - w0) % N for p in range(N)]
                        order = c3order
                        vis = {x: i for i, x in enumerate(order)}
                        base = [0] * N
                        for j, h in R.FEATURE_HELPERS.items():
                            base[vis[L[j]]] += h
                        fvis = [vis[L[j]] for j in range(9)]
                        for t in range(min(fvis), max(fvis) + 1):
                            base[t] += 1                    # p from the first to the last feature visit
                        win = R.windows(centres, order)
                        cand = list(win)
                        h, held = max_holds(win, cand, base, cap)
                        if h < 0:
                            continue
                        npairs = len(win)
                        out.append((2 * npairs - h, npairs, h, kind, r, d0, ctl, held, L))
    return out


def main(cap=5, top=15, c3=0):
    cap, top, c3 = int(cap), int(top), int(c3)
    order = walk_order(c3)
    assert order is not None
    results = []
    for cs in itertools.product(range(N), repeat=3):
        centres = (cs[0], cs[1], cs[2], c3, 25)
        if any(centres[i] == centres[i + 1] for i in range(4)) and not ALLOW_REP[0]:
            continue
        if any(centres[i] == centres[i + 1] for i in range(4)):
            continue
        for res in evaluate(centres, cap, order):
            results.append((res[0], centres, res))
    results.sort(key=lambda t: t[0])
    print(f"cap {cap}, c3 {c3}, walk {order[:6]}..., networks evaluated with affine controls: {len(results)}")
    seen = 0
    for cost, centres, res in results[:top]:
        _, npairs, h, kind, r, d0, ctl, held, L = res
        hs = [sum(1 for k in held if k[0] == s) for s in range(5)]
        print(f"cost {cost} pairs {npairs} held {h} {hs} centres {centres} kind {kind} r {r} d0 {d0} "
              f"ctl {[''.join(map(str, c)) for c in ctl]} L0..8 {L[:9]}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
