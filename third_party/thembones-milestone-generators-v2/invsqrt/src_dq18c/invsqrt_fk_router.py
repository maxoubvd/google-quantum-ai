"""F1 (third holder): exact model of the invsqrt endpoint router (26-slot ring, controlled reflection stages,
path-local walk, holds limited by free wires).

Model (read from the public core agents/composer/p/w2core_e2.py, parse only; our own code):
  - ring positions 0..25; word position i holds token i (i = 0..24: signs 36 - i; 25: the constant one);
  - stage s = controlled reflection x -> (c_s - x) mod 26 (pairs (x, c_s - x), x < partner), control = an affine
    function of the class flags f0..f3 (k in binary) and p = f1 (f2 ^ f3), terms as bitmasks (0 = constant 1);
  - class k must send token i to ring position SLOTS[(i + k) % 26] (SLOTS = the binding's basis slots);
  - the walk visits ring positions in a fixed order; a pair of stage s must be applied ("open") while any
    position of its cone (closure of {x, partner} under stages s+1..last) is visited: its window is
    [first visit, last visit] of the cone. A HELD pair costs 1 Toffoli and occupies one free wire over its window;
    an unheld pair costs 2 (applied before the walk, undone after, AND recomputed);
  - free wires during the walk: CAP (46 - 44 occupied = 2 new + 3 retired root bits = 5); feature helpers and
    the p wire also take free wires at their visits.
Commands:
  python invsqrt_fk_router.py check            reproduce the public network (shifts, windows, lanes, cost)
  python invsqrt_fk_router.py lanes [holds]    lane profile of the public schedule (holds: R2 hold list json)
"""
import itertools
import json
import sys

N = 26
SLOTS = [21, 22, 23, 24, 25] + list(range(21))          # basis.json 'slots': slot p (base abs position p)
PUB_CENTRES = (2, 6, 24, 0, 25)
PUB_CONTROLS = ((0, 4), (8,), (0, 2, 4, 8), (0, 1, 2), (0, 1))
FEATURE_HELPERS = {0: 1, 1: 1, 2: 1, 3: 1, 4: 2, 5: 2, 6: 1, 7: 1, 8: 1}   # helpers per feature j (public core)
CAP = 5


def flags(k):
    f = [(k >> i) & 1 for i in range(4)]
    p = f[1] & (f[2] ^ f[3])
    return f + [p]


def affine(terms, k):
    v = flags(k)
    out = 0
    for t in terms:
        if t == 0:
            out ^= 1
        else:
            out ^= v[t.bit_length() - 1]
    return out


def pairs(c):
    return [(x, (c - x) % N) for x in range(N) if x < (c - x) % N]


def compose(centres, ctl):
    """ctl[s][k] in {0,1}; returns perm[k][token] = final ring position."""
    out = []
    for k in range(len(ctl[0])):
        pos = list(range(N))                 # pos[token]
        for s, c in enumerate(centres):
            if ctl[s][k]:
                pos = [(c - x) % N for x in pos]
        out.append(pos)
    return out


def check_shifts(centres, ctl, slots=SLOTS):
    perm = compose(centres, ctl)
    bad = 0
    for k, pos in enumerate(perm):
        for i in range(N):
            if pos[i] != slots[(i + k) % N]:
                bad += 1
    return bad


def public_walk():
    order = []
    for k in range(13):
        order += [13 + k, (12 - k) % N]
    return order


def cones(centres, s, pr):
    cone = set(pr)
    for t in range(s + 1, len(centres)):
        c = centres[t]
        cone |= {(c - x) % N for x in cone}
    return cone


def windows(centres, order):
    vis = {x: i for i, x in enumerate(order)}
    out = {}
    for s, c in enumerate(centres):
        for pr in pairs(c):
            cone = cones(centres, s, pr)
            ts = [vis[x] for x in cone]
            out[(s, pr)] = (min(ts), max(ts))
    return out


def fixed_demand(order, slots=SLOTS, p_from=16):
    """feature helpers at their visits + the p wire from visit p_from (public: step 8) to the end."""
    vis = {x: i for i, x in enumerate(order)}
    d = [0] * len(order)
    for j, h in FEATURE_HELPERS.items():
        d[vis[slots[j]]] += h
    for t in range(p_from, len(order)):
        d[t] += 1
    return d


def lane_profile(win, held, base):
    d = list(base)
    for key in held:
        a, b = win[key]
        for t in range(a, b + 1):
            d[t] += 1
    return d


def cmd_check():
    ctl = [[affine(PUB_CONTROLS[s], k) for k in range(11)] for s in range(5)]
    print("controls per stage (k = 0..10):")
    for s in range(5):
        print(f"  S{s} centre {PUB_CENTRES[s]:2d} terms {PUB_CONTROLS[s]}: {''.join(map(str, ctl[s]))}")
    print("shift mismatches vs SLOTS[(i+k)%26]:", check_shifts(PUB_CENTRES, ctl))
    order = public_walk()
    win = windows(PUB_CENTRES, order)
    for s in range(5):
        ws = sorted((win[(s, pr)][1] - win[(s, pr)][0] + 1, pr) for pr in pairs(PUB_CENTRES[s]))
        print(f"  S{s} window lengths (visits): {[w for w, _ in ws]}")
    return ctl, order, win


def cmd_lanes(r2=None):
    ctl, order, win = cmd_check()
    if r2 is None:
        r2 = [(11, 13), (10, 14), (7, 17), (6, 18)]          # the shipped patch's default (W2_R2 unset)
    held = [(3, pr) for pr in pairs(0)] + [(4, pr) for pr in pairs(25)] + [(2, tuple(pr)) for pr in r2]
    base = fixed_demand(order)
    prof = lane_profile(win, held, base)
    print("lane demand per visit:", prof, "max", max(prof), "cap", CAP)
    npairs = sum(len(pairs(c)) for c in PUB_CENTRES)
    print(f"pairs {npairs}, held {len(held)}, cost {2 * npairs - len(held)}; lane-visits used {sum(prof)} of {CAP * 26}")
    # greedy max extra S2 holds under the cap
    free = [CAP - x for x in prof]
    extra = []
    for pr in sorted(pairs(24), key=lambda pr: win[(2, pr)][1] - win[(2, pr)][0]):
        if pr in [tuple(x) for x in r2]:
            continue
        a, b = win[(2, pr)]
        if all(free[t] > 0 for t in range(a, b + 1)):
            for t in range(a, b + 1):
                free[t] -= 1
            extra.append(pr)
    print("extra S2 holds that would fit:", extra)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "check"
    if cmd == "check":
        cmd_check()
    elif cmd == "lanes":
        cmd_lanes(json.loads(sys.argv[2]) if len(sys.argv) > 2 else None)
