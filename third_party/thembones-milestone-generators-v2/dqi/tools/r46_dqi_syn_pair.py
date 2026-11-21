#!/usr/bin/env python3
"""r46_dqi_syn_pair.py -- seat R46-DQI-SYN: pairs of hyperplanes of forms that together hold all 42 operand forms of
the kit's GF(2^12) product (exp/arithmetic/field.tensor(12)).  A masked operand register of 11 wires holding
F * lambda_j(qi) for a basis of one hyperplane H1 = {m : par(m & h1) = 0} serves the terms inside H1; ONE Toffoli
CCX(F, g(qi), wire) moves the register to a second hyperplane H2; if H1 u H2 holds every mask, a masked product
costs 42 + 1 instead of 42 + 12 (12-wire register) wires 11 instead of 12.
Usage: python r46_dqi_syn_pair.py [kit decoder dir]
"""
import os
import sys


def masks12(kit):
    sys.path.insert(0, os.path.expanduser('~/qcc/data/kickmix'))
    import kickmix as km
    sys.modules['circ'] = km
    sys.path.insert(0, kit)
    from exp.arithmetic import field as ar
    return [m for m, _ in ar.tensor(12)]


def main():
    kit = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
        '~/qcc/ws/dqi_int/src/source_cap902/playbook/kit/decoder')
    masks = masks12(kit)
    par = [bin(i).count('1') & 1 for i in range(4096)]
    P = [sum(par[h & m] << i for i, m in enumerate(masks)) for h in range(4096)]
    pairs = []
    for h1 in range(1, 4096):
        p1 = P[h1]
        for h2 in range(h1 + 1, 4096):
            if not p1 & P[h2]:
                pairs.append((h1, h2))
    print('terms', len(masks), 'covering pairs (no mask outside both hyperplanes):', len(pairs))

    def score(pr):
        h1, h2 = pr
        both = h1 & h2
        in1 = sum(1 for m in masks if not par[m & h1])
        return (0 if both else 1, -in1)
    pairs.sort(key=score)
    for h1, h2 in pairs[:8]:
        in1 = sum(1 for m in masks if not par[m & h1])
        in2 = sum(1 for m in masks if not par[m & h2])
        print('h1 %4d h2 %4d  inside H1 %d  inside H2 %d  common bit of h1, h2: %s' % (
            h1, h2, in1, in2, 'yes' if h1 & h2 else 'no'))


if __name__ == '__main__':
    main()
