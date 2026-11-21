#!/usr/bin/env python3
"""invsqrt_fc2_lanes.py -- FC2 (own code): swap network that shifts a lane of L positions up by k (k = 0..K) when
its top k positions hold zeros, as a low-bit-first binary butterfly with carry selectors (the mechanism of the
public generator's zero-supported butterfly, Apache-2.0, Syed Asad Ali; re-derived here for pair lanes).
Stage j swaps (p, p + 2^j) for p with bit j clear; its selector is bit_j(k) xor [p mod 2^j < k mod 2^j].
A pair is kept only if, for some k where it fires, one of its two positions holds a live token.
Each selector is a function of k, written as a parity of the thermometer bits P_j = [k >= j] (P_0 = 1).
  invsqrt_fc2_lanes.py L K      -> prints the network and checks it on every k
"""
import json
import sys


def selector(k, p, j):
    m = 1 << j
    return ((k >> j) & 1) ^ int((p % m) < (k % m))


def prefix_terms(truth, K):
    """indices j with truth(k) = parity of [k >= j] over the terms (j = 0 means the constant 1)."""
    terms = [j for j in range(K + 1) if ((truth >> j) & 1) ^ (((truth >> (j - 1)) & 1) if j else 0)]
    for k in range(K + 1):
        assert (sum(k >= j for j in terms) & 1) == ((truth >> k) & 1)
    return terms


def network(L, K):
    nst = max(1, (K).bit_length())
    net = []
    for j in range(nst):
        d = 1 << j
        for p in range(L):
            q = p + d
            if p & d or q >= L:
                continue
            net.append((j, p, q))
    # prune: simulate tokens per k, keep a pair if it ever moves a live token
    keep = []
    for (j, p, q) in net:
        truth = sum(selector(k, p, j) << k for k in range(K + 1))
        keep.append(dict(stage=j, a=p, b=q, truth=truth, live=False))
    for k in range(K + 1):
        st = list(range(L - k)) + [None] * k          # token i at position i (live), zeros above
        for op in keep:
            if (op["truth"] >> k) & 1:
                a, b = op["a"], op["b"]
                if st[a] is not None or st[b] is not None:
                    op["live"] = True
                st[a], st[b] = st[b], st[a]
    net2 = [op for op in keep if op["live"]]
    # verify: every k moves token i to i + k, and the reversed network restores
    for k in range(K + 1):
        st = list(range(L - k)) + [None] * k
        for op in net2:
            if (op["truth"] >> k) & 1:
                st[op["a"]], st[op["b"]] = st[op["b"]], st[op["a"]]
        assert st == [None] * k + list(range(L - k)), (k, st)
    for op in net2:
        op["terms"] = prefix_terms(op["truth"], K)
    return net2


if __name__ == "__main__":
    L, K = int(sys.argv[1]), int(sys.argv[2])
    net = network(L, K)
    print(json.dumps(dict(L=L, K=K, swaps=len(net),
                          per_stage=[sum(1 for o in net if o["stage"] == j) for j in range(4)],
                          net=[(o["stage"], o["a"], o["b"], o["terms"]) for o in net])))
