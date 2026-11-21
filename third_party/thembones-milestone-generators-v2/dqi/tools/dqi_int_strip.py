#!/usr/bin/env python3
"""dqi_int_strip.py -- B5 (16th holder, Nov 16): remove the profiling markers (DEBUG_PRINT 0x07, DEBUG_PRINT_Q 0x17,
DEBUG_PRINT_B 0x47, with or without the condition flag) from a marker build, and optionally compare the result with
a plain build of the same recipe (ops and every operand array).

  python dqi_int_strip.py <marker.kmb> [<out.kmb>|-] [<plain.kmb> to compare]

Prints the number of marker ops and, with a plain file, whether the stripped arrays equal the plain file's.
Also prints the map size: stripped op index = marker op index - (markers before it).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kmb  # noqa: E402
import dqi_int_kmbw as kw  # noqa: E402

MARK = (0x07, 0x17, 0x47)


def marker_mask(ops):
    b = ops & 0x7F
    m = np.zeros(len(ops), bool)
    for x in MARK:
        m |= b == x
    m &= ops != 0x81
    return m


def strip(k):
    ops = np.asarray(k.ops)
    m = marker_mask(ops)
    return kw.edit(k, np.flatnonzero(m).tolist(), {}), int(m.sum())


def same(a, b):
    out = {}
    out['ops'] = np.array_equal(np.asarray(a.ops), np.asarray(b.ops))
    for j in range(3):
        out['q3_%d' % j] = np.array_equal(np.asarray(a.q3[j]), np.asarray(b.q3[j]))
    for j in range(2):
        out['q2_%d' % j] = np.array_equal(np.asarray(a.q2[j]), np.asarray(b.q2[j]))
    out['q1'] = np.array_equal(np.asarray(a.q1), np.asarray(b.q1))
    out['bit_t'] = np.array_equal(np.asarray(a.bit_t), np.asarray(b.bit_t))
    out['cond'] = np.array_equal(np.asarray(a.cond), np.asarray(b.cond))
    out['angles'] = np.array_equal(np.asarray(a.angles), np.asarray(b.angles))
    out['nq'] = (a.num_qubits, b.num_qubits)
    out['nbits'] = (a.num_bits, b.num_bits)
    return out


def main():
    k = kmb.read(sys.argv[1])
    s, nm = strip(k)
    print('markers removed %d of %d ops' % (nm, len(k.ops)), flush=True)
    if len(sys.argv) > 2 and sys.argv[2] != '-':
        s.num_qubits = k.num_qubits
        kw.write(sys.argv[2], s)
        print('written', sys.argv[2])
    if len(sys.argv) > 3:
        p = kmb.read(sys.argv[3])
        r = same(s, p)
        print('compare with %s: %s' % (sys.argv[3], r))
        print('ALL EQUAL' if all(v for kk, v in r.items() if kk not in ('nq', 'nbits')) else 'DIFFERENT')


if __name__ == '__main__':
    main()
