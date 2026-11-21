"""T10 seat 2: writes tools/chem_t10_ipolg.py = R8's tools/chem_own2_ipolish.py (read-only, copied) with ONE change:
the wire relations are checked against the INPUT's own bit-0 XORs instead of the plain masks, so a table with R65's
fitted block constants (a G table: target = mask parity XOR g(block)) is accepted. R8's move set keeps every wire's
bit-0 XOR by construction (+-2 on wire columns, (+-1, +-1) pair moves), so the G relations of the input hold on the
output; tools/chem_wr7_g2rel.py --check confirms it outside."""
import hashlib, sys
src = open(sys.argv[1], encoding="utf-8").read()
old = ("    def broken(wv_):\n"
       "        return sum(int(np.sum(((wv_[:, c2] ^ wv_[:, c1]) & 1) != wlin(addr, m_))) for c1, c2, m_ in pairs)\n")
new = ("    REF_ = {(c1, c2): ((wv[:, c2] ^ wv[:, c1]) & 1) for c1, c2, m_ in pairs}    # t10: the input's own bit-0 XORs (G table)\n"
       "    def broken(wv_):\n"
       "        return sum(int(np.sum(((wv_[:, c2] ^ wv_[:, c1]) & 1) != REF_[(c1, c2)])) for c1, c2, m_ in pairs)\n")
assert src.count(old) == 1
src = src.replace(old, new)
src = src.replace('"""chem_own2_ipolish', '"""T10 copy (tools/chem_t10_mkipolg.py): wire relations checked against the input\'s own bit-0 XORs.\nchem_own2_ipolish', 1)
open(sys.argv[2], "w", encoding="utf-8", newline="\n").write(src)
print("written", sys.argv[2], hashlib.md5(src.encode()).hexdigest()[:8])
