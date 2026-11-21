"""T10 seat 3 (Nov 18, ledger row O14): writes tools/chem_t10_frelgc.py =
  tools/r65_chem_all15_frelg.py (R65's relation tool with fitted block constants: RK7's RelPlacer + T10's smart
  protection + --relp-g; read-only)
+ R67-CHEM-PATH's CONDITIONAL FREEZE ORDER (--rank-cond / --rank-geo / --cond-cpus), copied from tools/r67_chem_path_round.py
  exactly as tools/chem_t10_mkgjc.py copies it into the G tool (THE RECIPE), with seat 2's END guard
+ --relp-e MODE (exact | blind | off): R53-CHEM-GAUGE's carrier rule, every relation asked of y = eps x z
  (tools/chem_t10_relge.py RelPlacerGE), and at the END the relations counted on C2's STORED words.
Without --relp-e and --rank-cond the output behaves as r65_chem_all15_frelg.py.
    python tools/chem_t10_mkfrelgc.py tools/r65_chem_all15_frelg.py tools/r67_chem_path_round.py tools/chem_t10_frelgc.py
"""
import hashlib
import sys

src = open(sys.argv[1], encoding="utf-8").read().replace("\r\n", "\n")
r67 = open(sys.argv[2], encoding="utf-8").read().replace("\r\n", "\n")
NL = chr(10)


def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:90])
    src = src.replace(old, new)


def cut(text, start, end):
    i = text.index(start)
    j = text.index(end, i)
    return text[i:j]


# ---- (1) R67's conditional order (as chem_t10_mkgjc.py; the set-up goes after frelg's sdyn set-up)
opts = cut(r67, '    ap.add_argument("--rank-cond"', '    ap.add_argument("--reset-free"')
setup = cut(r67, '    COND = a.rank_cond > 0', '    def rules():')
over = cut(r67, '                if COND:        # r67: the order by the conditional stiffness', '            elif a.rank == "sratio":')
assert opts.count("add_argument") == 3 and "cond_pool" in setup and "impact = val_c" in over
sub('    ap.add_argument("--lm-until",', opts + '    ap.add_argument("--lm-until",')
sub('        sdyn_job, sdyn_from = None, None\n\n    def rules():',
    '        sdyn_job, sdyn_from = None, None\n' + setup.rstrip("\n") + '\n\n    def rules():')
sub('                impact = rank_sk + 1e-30 * c\n            elif a.rank == "sratio":',
    '                impact = rank_sk + 1e-30 * c\n' + over + '            elif a.rank == "sratio":')
sub('        rank_cd = cond_get(cond_start(np.where(frozen, qv, theta), w, frozen.copy()))' + NL,
    NL.join(['        try:',
             '            rank_cd = cond_get(cond_start(np.where(frozen, qv, theta), w, frozen.copy()))',
             '        except Exception as ex_:     # t10: e.g. S_FF not positive definite near the END: S_kk instead',
             '            log(f"t10 cond: WARNING the first conditional table failed ({ex_!r}): S_kk in use until one succeeds")',
             '            rank_cd = rank_sk.copy()']) + NL)
sub('                        rank_cd = cond_get(cond_job)' + NL,
    NL.join(['                        try:',
             '                            rank_cd = cond_get(cond_job)',
             '                        except Exception as ex_:     # t10: the last finished table stays in use',
             '                            log(f"t10 cond: WARNING a conditional table failed ({ex_!r}); the last one stays in use")']) + NL)

# ---- (2) --relp-e: R53's carrier rule in the placer
sub('    ap.add_argument("--relp-g-rounds",',
    '    ap.add_argument("--relp-e", default="", choices=("", "exact", "blind", "off"), help="t10 seat 3: with --relp '
    '--relp-g: R53\'s carrier rule: every relation asked of y = eps x z (tools/chem_t10_relge.py RelPlacerGE; exact = '
    'eps from the rows\' words at every freeze); the END counts the relations on C2\'s stored words")\n'
    '    ap.add_argument("--relp-g-rounds",')
sub('                FP = RelPlacerG(R, B, RELS, n_one=N, ncols=N - 1, rounds=a.relp_g_rounds)\n',
    '                if a.relp_e:     # t10 seat 3: the carrier rule (relations on y = eps x z)\n'
    '                    from chem_t10_relge import RelPlacerGE\n'
    '                    FP = RelPlacerGE(R, B, RELS, edges, n_one=N, ncols=N - 1, rounds=a.relp_g_rounds, mode=a.relp_e, '
    'brot=a.bits)\n'
    '                    log(f"relp-e (t10): the relations asked of y = eps x z, eps mode {a.relp_e} (R53\'s carrier rule)")\n'
    '                else:\n'
    '                    FP = RelPlacerG(R, B, RELS, n_one=N, ncols=N - 1, rounds=a.relp_g_rounds)\n')
# a per-freeze line: eps changes on frozen relation words and rows without a consistent pattern
sub('                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n',
    '                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n'
    '                if a.relp_e:     # t10 seat 3\n'
    '                    log(f"relp-e (t10): rows without a consistent pattern {FP.infeasible}; late eps changes on frozen odd '
    'relation words so far {FP.eps_changes}")\n')
# the END: relations on C2's stored words (the table step applies C2's unchanged carrier gauge)
sub('            log(f"fp: rows placed; relations broken {fpbad} (wires {len(WIRES)}, code chains {len(FCH)})"\n',
    '            if a.relp and a.relp_g and a.relp_e:     # t10 seat 3: the words the loaders will read\n'
    '                from chem_t10_relge import stored_words as _t10_stored\n'
    '                fpst = _t10_stored(np.mod(wv, 128), edges, np.asarray(addr), a.bits)\n'
    '                fpbad_st = int(_rcheck_g(fpst, np.asarray(addr), FP.rels_used()))\n'
    '                log(f"relp-e (t10): relations broken on C2\'s STORED words {fpbad_st} (on the signed-range words '
    '{fpbad}: expected nonzero under the rule); late eps changes {FP.eps_changes}")\n'
    '                fpbad = fpbad_st\n'
    '            log(f"fp: rows placed; relations broken {fpbad} (wires {len(WIRES)}, code chains {len(FCH)})"\n')
# ---- (3) R74-CHEM-ORDER's opportunistic fill (--opp L --opp-wait W --opp-until S) and the re-priced placer (--fp-cost
# cond), copied from tools/r74_chem_order_fround.py (read-only; argv[4]): the same lines, the cond start (--cond-from) = 0
if len(sys.argv) > 4:
    r74 = open(sys.argv[4], encoding="utf-8").read().replace("\r\n", "\n")
    o_opts = cut(r74, '    ap.add_argument("--opp", type=float', '    ap.add_argument("--fp-cost"')
    o_fpc = cut(r74, '    ap.add_argument("--fp-cost"', '\n') + '\n'
    o_blk = cut(r74, '                if a.opp > 1 and COND and rank_cd is not None', '                if a.harv > 0 and COND')
    assert o_opts.count("add_argument") == 3 and "r74_wait" in o_blk
    o_blk = o_blk.replace("a.cond_from - 1e-9 <= float(frozen.mean()) < a.opp_until", "float(frozen.mean()) < a.opp_until")
    sub('    ap.add_argument("--relp-e",', o_opts + o_fpc + '    ap.add_argument("--relp-e",')
    sub('                pick = np.argsort(-impact, axis=None)[:need]\n            if FP is not None:',
        '                pick = np.argsort(-impact, axis=None)[:need]\n' + o_blk + '            if FP is not None:')
    sub('                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n',
        '                if a.fp_cost == "cond" and COND and rank_cd is not None:      # r74: the placer priced by the conditional stiffness\n'
        '                    fpc = np.maximum(rank_cd, 0.0) + 1e-30 * c\n'
        '                fpr, fpm = FP.residues(fpt, fpc, frozen, np.round(qv * scale).astype(np.int64))\n')
src = src.replace('"""r65_chem_all15_frelg = chem_t10_frel.py',
                  '"""chem_t10_frelgc (T10 seat 3, tools/chem_t10_mkfrelgc.py: + R67\'s --rank-cond + --relp-e) of '
                  'r65_chem_all15_frelg = chem_t10_frel.py', 1)
open(sys.argv[3], "w", encoding="utf-8", newline="\n").write(src)
print("written", sys.argv[3], hashlib.md5(src.encode()).hexdigest()[:8])
