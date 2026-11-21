"""chem_t10_mkfrel (T10, Nov 17): writes tools/chem_t10_frel.py = RK7's tools/chem_rk7_frel.py (read-only, copied) plus
ONE option, --relp-at F (THE RELATION WALL): past frozen share F the free words of the bit-1 relation columns of the
relation file (RelPlacer modulus 4) freeze FIRST, while the other free words can still compensate their parity-forced
moves (R8's --slack-at for XOR slack, measured +7 at F 0.61 against +23.5 when the slack froze last [M ws_chem_own2
09:20Z]). Default 1.0 = the original tool's behaviour exactly (the patch is inert).
    python chem_t10_mkfrel.py      (re-run after any change of chem_rk7_frel.py)
"""
import os

T = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(T, "chem_rk7_frel.py"), encoding="utf-8").read()

patches = [
    # 1. the option
    ('    ap.add_argument("--relp", default="",',
     '    ap.add_argument("--relp-boost", type=float, default=0.0, help="chem_t10: with --relp-at: multiply the bit-1 relation "\n'
     '                    "columns\' freeze priority by this factor instead of the wall (0 = the wall)")\n'
     '    ap.add_argument("--relp-at", type=float, default=1.0, help="chem_t10: THE RELATION WALL: past this frozen share the "\n'
     '                    "free words of the bit-1 relation columns (RelPlacer modulus 4) freeze first (1.0 = off)")\n'
     '    ap.add_argument("--relp", default="",'),
    # 2. the wall's columns, once the RelPlacer exists
    ('    FP = fp_place0 = None                           # r50_chem_lbl third hold: the FreezePlacer when --fp is given',
     '    FP = fp_place0 = None                           # r50_chem_lbl third hold: the FreezePlacer when --fp is given\n'
     '    RWALL = None                                    # chem_t10: the relation wall\'s columns (bit-1 relation columns)'),
    ('            log(f"relp: {len(RELS)} relations a row on {len(FP.pl.cc)} columns (moduli {sorted(set(FP.mod[FP.pl.cc].tolist()))}), batches {FP.lo}-{FP.hi}")',
     '            log(f"relp: {len(RELS)} relations a row on {len(FP.pl.cc)} columns (moduli {sorted(set(FP.mod[FP.pl.cc].tolist()))}), batches {FP.lo}-{FP.hi}")\n'
     '            if a.relp_at < 1.0:     # chem_t10: the relation wall\n'
     '                RWALL = np.asarray(FP.mod, np.int64) >= 4\n'
     '                log(f"relation wall (chem_t10): past frozen share {a.relp_at} the free words of {int(RWALL.sum())} bit-1 relation columns freeze first")'),
    # 2b. THE SMART PROTECTION (--smart-protect): frel's gauge protects a frozen word only on columns with two coded
    # bits (pcbits >= 2); coding bit 1 there protects EVERY odd frozen word of a bit-1 relation column and walls
    # 350-710 free entries a stage (+2.3 mHa^2 at 8% frozen against the unprotected twin [M T10 16:25Z]); without it
    # the gauge breaks relations (P50b). A negation flips bit 1 of an odd word and keeps bit 0, so a row stays
    # feasible at its place exactly when the parity change can be made up by the row's FREE bit-1 words: the word's
    # relation-incidence vector lies in the GF(2) span of the free words' vectors. Only words outside it are protected.
    ('    RWALL = None                                    # chem_t10: the relation wall\'s columns (bit-1 relation columns)',
     '    RWALL = None                                    # chem_t10: the relation wall\'s columns (bit-1 relation columns)\n'
     '    SMART = None                                    # chem_t10: the smart protection (a function of frozen, words)'),
    ('            if a.relp_at < 1.0:     # chem_t10: the relation wall\n',
     '            if a.smart_protect:     # chem_t10: the smart protection of the bit-1 relations (names kept local:\n'
     '                SMART, sp_inc = _t10_smart_factory(RELS, N - 1)      # main\'s i_ / j_ feed the JAX objective)\n'
     '                assert not np.any((pcbits >= 2) & (sp_inc != 0)), "--smart-protect: leave bit 1 out of --codes"\n'
     '                log(f"smart protection (chem_t10): {int(np.unique(sp_inc[sp_inc != 0]).size)} incidence classes on "\n'
     '                    f"{int((sp_inc != 0).sum())} bit-1 relation columns")\n'
     '            if a.relp_at < 1.0:     # chem_t10: the relation wall\n'),
    # the factory at module level, before negatable()
    ('def negatable(frozen, wv, protected, NB):\n',
     'def _t10_smart_factory(rels, ncols):\n'
     '    """chem_t10: SMART(frozen, words) = the frozen odd words of the bit-1 relation columns whose negation (bit 1 flips,\n'
     '    bit 0 stays) the row\'s FREE bit-1 words cannot make up: incidence vector outside the GF(2) span of the free ones."""\n'
     '    rel1 = [i for i, (k, o, t) in enumerate(rels) if any(b == 1 for (_j, b) in o)]\n'
     '    inc = np.zeros(ncols, np.int64)\n'
     '    for bi, i in enumerate(rel1):\n'
     '        k, o, _t = rels[i]\n'
     '        for (j, b) in o:\n'
     '            if b == 1:\n'
     '                inc[5 * k + j] ^= 1 << bi\n'
     '                inc[5 * k - 5 + j] ^= 1 << bi\n'
     '    cols = np.flatnonzero(inc)\n'
     '    idx = np.arange(1 << len(rel1))\n'
     '\n'
     '    def smart(fr, wv):\n'
     '        reach = np.zeros((fr.shape[0], len(idx)), bool)\n'
     '        reach[:, 0] = True\n'
     '        for c in cols:                  # the span of the free words\' incidence vectors, row by row\n'
     '            r = ~fr[:, c]\n'
     '            if r.any():\n'
     '                reach[r] |= reach[r][:, idx ^ int(inc[c])]\n'
     '        out = np.zeros(fr.shape, bool)\n'
     '        out[:, cols] = fr[:, cols] & (np.mod(wv[:, cols], 2) == 1) & ~reach[:, inc[cols]]\n'
     '        return out\n'
     '    return smart, inc\n'
     '\n'
     '\n'
     'def negatable(frozen, wv, protected, NB):\n'),
    ('        return negatable(frozen, wv, protected, NB)\n',
     '        if SMART is not None:   # chem_t10: the exact minimal protection of the bit-1 relations\n'
     '            protected = protected | SMART(frozen, wv)\n'
     '        return negatable(frozen, wv, protected, NB)\n'),
    ('    ap.add_argument("--relp", default="",',
     '    ap.add_argument("--smart-protect", action="store_true", help="chem_t10: protect a frozen odd word of a bit-1 "\n'
     '                    "relation column from the sign gauge only when the row\'s free words cannot make up its parity")\n'
     '    ap.add_argument("--relp", default="",'),
    # 3. the wall itself, just before frozen entries are taken out of the ranking
    ('            impact[frozen] = -np.inf\n',
     '            if RWALL is not None and frac > a.relp_at:      # chem_t10: the relation wall\n'
     '                rw_free = int(np.sum(RWALL[None, :] & ~frozen))\n'
     '                if a.relp_boost > 0:    # chem_t10 fix round 1 (18:30Z): a BOOST, not a wall: the wall deferred the\n'
     '                    # stiff non-wall words to its end and gave all back (+7.5 over G9s2 at 88% [M t10_m2bw]); x M moves\n'
     '                    # the soft bit-1 payers to mid-order while stiffer words still go first\n'
     '                    impact = np.where(RWALL[None, :], impact * a.relp_boost, impact)\n'
     '                else:\n'
     '                    rw_off = 2.0 * float(np.max(np.where(np.isfinite(impact), np.abs(impact), 0.0))) + 1.0\n'
     '                    impact = np.where(RWALL[None, :], impact + rw_off, impact)      # above all others, own order kept\n'
     '                log(f"relation wall (chem_t10): {rw_free} free words of the bit-1 relation columns before this stage")\n'
     '            impact[frozen] = -np.inf\n'),
]
for old, new in patches:
    assert src.count(old) == 1, f"anchor not unique or missing: {old[:70]!r}"
    src = src.replace(old, new)
src = src.replace('"""chem_rk7_frel', '"""chem_t10_frel (T10, Nov 17: chem_rk7_frel + --relp-at, written by tools/chem_t10_mkfrel.py)\n'
                  'chem_rk7_frel', 1)
open(os.path.join(T, "chem_t10_frel.py"), "w", encoding="utf-8", newline="\n").write(src)
print("written tools/chem_t10_frel.py")
