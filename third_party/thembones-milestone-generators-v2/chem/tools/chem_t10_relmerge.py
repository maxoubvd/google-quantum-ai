"""T10 seat 3: THE RELATION FILE OF A WHOLE RELATION TABLE: the rounding's fitted two-body targets (rel_used_g.txt: block
constants on blocks 5..; block 4 pinned) merged with the one-body step's fitted constants (rel_ob.txt from
chem_t10_obcar.py mode ge: blocks 0..3; block 4 pinned). Same relations, same order, same base targets (asserted).
    python chem_t10_relmerge.py rel_used_g.txt rel_ob.txt OUT.txt"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r65_chem_all15_relg import fmt_g, parse_g, split_target  # noqa: E402

two = parse_g("@" + sys.argv[1])
one = parse_g("@" + sys.argv[2])
assert len(two) == len(one), (len(two), len(one))
out = []
for (k, o, t2), (k1, o1, t1) in zip(two, one):
    assert k == k1 and tuple(o) == tuple(o1), (k, o, k1, o1)
    b2, g2 = split_target(t2)
    b1, g1 = split_target(t1)
    assert str(b2) == str(b1), (b2, b1)
    assert (g2 & 0x1F) == 0 and (g1 >> 4) == 0, "blocks overlap"
    g = g2 | g1
    out.append((k, o, t2 if g == 0 else f"G/{b2}/{g:x}"))
open(sys.argv[3], "w").write("# t10 seat 3: two-body constants (blocks 5..) + one-body constants (blocks 0..3)\n" + fmt_g(out))
print(f"wrote {sys.argv[3]}: {len(out)} relations")
