#!/usr/bin/env python3
"""s5_inv_h1end_ndiff.py -- S5-INV-H1END: where do two invsqrt kmx texts differ once measurement-bit NAMES are normalized
(each file's bits renamed in order of first appearance)?  Prints the common prefix length, the common suffix length,
the differing windows of both files and a few lines of each.
  python -I s5_inv_h1end_ndiff.py A.kmx B.kmx"""
import os
import re
import sys

BIT = re.compile(r"\bb(\d+)\b")


def norm(path):
    lines = [l.strip() for l in open(os.path.expanduser(path)) if l.strip()]
    out = [BIT.sub("b", l) for l in lines]          # bit names blanked: the window is located on the op structure
    return lines, out


def main():
    la, a = norm(sys.argv[1])
    lb, b = norm(sys.argv[2])
    p = 0
    while p < min(len(a), len(b)) and a[p] == b[p]:
        p += 1
    s = 0
    while s < min(len(a), len(b)) - p and a[len(a) - 1 - s] == b[len(b) - 1 - s]:
        s += 1
    print("A %d lines, B %d lines; common prefix %d, common suffix %d (bit names normalized)" % (len(a), len(b), p, s))
    print("A differs in [%d, %d) (%d lines); B differs in [%d, %d) (%d lines)" % (
        p, len(a) - s, len(a) - s - p, p, len(b) - s, len(b) - s - p))
    for tag, ll, lo in (("A", la, p), ("B", lb, p)):
        print("--- %s first differing lines:" % tag)
        for l in ll[lo:lo + 6]:
            print("   ", l)


if __name__ == "__main__":
    main()
