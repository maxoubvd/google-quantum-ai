#!/usr/bin/env python3
"""dqi_int_lawfit.py -- workstream dqi_int (B5, Nov 14): find the missing route term of Z2's pdr on gap charts.

Reads a converter-harness output of dqi_int_pdrdbg (pdr in debug mode + dqi_int_pdrfix with PDRFIX_XDBG=1): per
block U and payload k the residual Delta = D ^ R ^ L on every chart, and the candidate words copied at U's route
(xS<U>_<m>_<k> payload tops, xQ<U>_<m> = q0_m', xH<U>_<m> = q1_m', flags xG<name><U>). Candidate terms: each top
s_m, each product q0_m' s_m'' and q1_m' s_m'' (GF(2^12), x^12 + x^6 + x^4 + x + 1). Solves over GF(2) (all bits of all
charts) for a subset of terms whose XOR equals Delta, and prints it (or 'no law in the span').
  python dqi_int_lawfit.py <out_*.txt> [U]
"""
import ast
import re
import sys

MOD = 0x1053


def gmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
        if a & 0x1000:
            a ^= MOD
    return r


def parse(path):
    delta, x = {}, {}
    for line in open(path, errors='replace'):
        m = re.match(r'pdrdbg D(\d+)_(\d+): .*D\^R\^L (\[.*\])', line)
        if m:
            delta[(int(m.group(1)), int(m.group(2)))] = ast.literal_eval(m.group(3))
            continue
        m = re.match(r'pdrdbg (x\w+): (\[.*\])', line)
        if m:
            x[m.group(1)] = ast.literal_eval(m.group(2))
    return delta, x


def solve(target, cands):
    """GF(2) elimination: rows = candidate bit-vectors (python ints), find subset XOR == target."""
    basis = []                                   # (pivot bit, vector, subset bitmask)
    for i, (_, v) in enumerate(cands):
        s = 1 << i
        for p, bv, bs in basis:
            if v >> p & 1:
                v ^= bv
                s ^= bs
        if v:
            p = v.bit_length() - 1
            basis.append((p, v, s))
            basis.sort(key=lambda t: -t[0])
    v, s = target, 0
    for p, bv, bs in basis:
        if v >> p & 1:
            v ^= bv
            s ^= bs
    if v:
        return None
    return [cands[i][0] for i in range(len(cands)) if s >> i & 1]


def vec(vals):
    out = 0
    for i, w in enumerate(vals):
        out |= (w & 0xFFF) << (12 * i)
    return out


def main(path, Uonly=None):
    delta, x = parse(path)
    for (U, k), d in sorted(delta.items()):
        if Uonly is not None and U != Uonly:
            continue
        if not any(d):
            print('U %d k %d: Delta = 0 on all %d charts' % (U, k, len(d)))
            continue
        n = len(d)
        flags = {nm: x.get('xG%s%d' % (nm, U)) for nm in ('B', 'C', 'F', 'C2')}
        print('U %d k %d: Delta nonzero on %d / %d charts; flags (chart 0) %s' % (
            U, k, sum(1 for t in d if t), n, {a: (b[0] & 1 if b else None) for a, b in flags.items()}))
        S = {int(nm.split('_')[1]): v for nm, v in x.items() if nm.startswith('xS%d_' % U) and nm.endswith('_%d' % k)}
        Q = {int(nm.split('_')[1]): v for nm, v in x.items() if nm.startswith('xQ%d_' % U)}
        H = {int(nm.split('_')[1]): v for nm, v in x.items() if nm.startswith('xH%d_' % U)}
        cands = []
        for m, v in sorted(S.items()):
            cands.append(('s%d' % m, vec(v)))
        for name, src in (('q0', Q), ('q1', H)):
            for m, a in sorted(src.items()):
                for m2, s in sorted(S.items()):
                    cands.append(('%s_%d*s%d' % (name, m, m2), vec([gmul(a[i], s[i]) for i in range(n)])))
        law = solve(vec(d), cands)
        print('   law: %s' % (' + '.join(law) if law else 'no law in the span of %d terms' % len(cands)))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None)
