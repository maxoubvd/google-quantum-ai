"""r65_chem_all15_mkfrelg (R65-CHEM-ALL15, second hold): writes tools/r65_chem_all15_frelg.py = T10's
tools/chem_t10_frel.py (read-only, copied whole; itself RK7's chem_rk7_frel.py + the relation wall + the smart
protection) plus ONE option, --relp-g: THE BLOCK CONSTANTS OF EVERY RELATION FITTED AT FREEZE TIME
(tools/r65_chem_all15_relg.py: RelPlacerG, a subclass of RK7's RelPlacer). Without --relp-g the written tool is T10's
in behaviour (the patch is inert). With it:
  * the relation file may carry G/<base>/<hex> targets (a rel_used_g.txt of an earlier run);
  * the constants live in state.npz as fp_G (15 x blocks, a bit vector over the loader's relations) next to fp_place
    and are restored on a resume; a state without fp_G (a branch from a kept head of a run without constants) starts
    at 0: every frozen relation of that head holds at 0;
  * at the END the relations are checked WITH the constants, `rel_used_g.txt` (the relation file with the fitted
    targets) is written into the run folder and fp_G into rounded.npz.
    python r65_chem_all15_mkfrelg.py      (re-run after any change of chem_t10_frel.py)
"""
import os

T = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(T, "chem_t10_frel.py"), encoding="utf-8").read().replace("\r\n", "\n")
ONES = "int(sum(bin(int(v)).count(chr(49)) for v in FP.G.reshape(-1)))"

patches = [
    # 1. the options
    ('    ap.add_argument("--relp", default="",',
     '    ap.add_argument("--relp-g", action="store_true", help="r65: with --relp: the block constants of every relation fitted at freeze "\n'
     '                    "time (tools/r65_chem_all15_relg.py RelPlacerG); targets stay in the family f(low 4 bits) XOR g(block)")\n'
     '    ap.add_argument("--relp-g-rounds", type=int, default=1, help="r65: alternations (constants, assignment) a freeze")\n'
     '    ap.add_argument("--relp", default="",'),
    # 2. the state of a resume
    ('    FP = fp_place0 = None                           # r50_chem_lbl third hold: the FreezePlacer when --fp is given',
     '    FP = fp_place0 = None                           # r50_chem_lbl third hold: the FreezePlacer when --fp is given\n'
     '    fp_G0 = None                                    # r65: the fitted block constants of a resumed state'),
    ('            if "fp_place" in s.files:           # r50_chem_lbl third hold\n'
     '                fp_place0 = np.array(s["fp_place"], np.int64)',
     '            if "fp_place" in s.files:           # r50_chem_lbl third hold\n'
     '                fp_place0 = np.array(s["fp_place"], np.int64)\n'
     '            if "fp_G" in s.files:               # r65\n'
     '                fp_G0 = np.array(s["fp_G"], np.int64)'),
    # 3. the placer
    ('            RELS = _rparse(a.relp)\n'
     '            assert not WIRES and not FCH, "--relp replaces --wires / --fp-code"\n'
     '            FP = RelPlacer(R, B, RELS, n_one=N, ncols=N - 1)',
     '            if a.relp_g:        # r65: the same placer with the block constants of every relation fitted\n'
     '                from r65_chem_all15_relg import RelPlacerG, check_g as _rcheck_g, fmt_g as _rfmt_g, parse_g as _rparse_g\n'
     '                RELS = _rparse_g(a.relp)\n'
     '            else:\n'
     '                RELS = _rparse(a.relp)\n'
     '            assert not WIRES and not FCH, "--relp replaces --wires / --fp-code"\n'
     '            if a.relp_g:\n'
     '                FP = RelPlacerG(R, B, RELS, n_one=N, ncols=N - 1, rounds=a.relp_g_rounds)\n'
     '                RELS = list(FP.rels)        # the base targets (the smart protection reads the outputs only)\n'
     '                log(f"relp-g (r65): block constants fitted before every freeze ({a.relp_g_rounds} alternation(s)); block 4 pinned; "\n'
     '                    f"constants set in the relation file: {' + ONES + '}")\n'
     '            else:\n'
     '                FP = RelPlacer(R, B, RELS, n_one=N, ncols=N - 1)'),
    ('        if fp_place0 is not None:\n'
     '            FP.place = fp_place0.copy()\n',
     '        if fp_place0 is not None:\n'
     '            FP.place = fp_place0.copy()\n'
     '        if a.relp and a.relp_g:        # r65\n'
     '            if fp_G0 is not None:\n'
     '                FP.set_state(FP.place, fp_G0)\n'
     '            elif frozen.any():\n'
     '                log("relp-g (r65): a state without fp_G: the constants start from the relation file (0 for a plain file): the frozen "\n'
     '                    "relations of a head that ran without constants hold")\n'),
    # 4. the END: the check with the constants, the fitted relation file, fp_G in rounded.npz
    ('            if a.relp:      # chem_rk7: my relations on the placed rows\n'
     '                fpbad = int(_rcheck(np.mod(wv, 128), np.asarray(addr), RELS))',
     '            if a.relp and a.relp_g:     # r65: the relations with their fitted block constants\n'
     '                fpbad = int(_rcheck_g(np.mod(wv, 128), np.asarray(addr), FP.rels_used()))\n'
     '                open(os.path.join(a.out, "rel_used_g.txt"), "w").write(_rfmt_g(FP.rels_used()))\n'
     '                log(f"relp-g (r65): rel_used_g.txt written; constants set to 1: {' + ONES + '}; "\n'
     '                    f"(loader, block) vectors changed in the run: {FP.gflips}")\n'
     '            elif a.relp:    # chem_rk7: my relations on the placed rows\n'
     '                fpbad = int(_rcheck(np.mod(wv, 128), np.asarray(addr), RELS))'),
    ('                 pairs=np.stack([pm_rows, ps_rows], axis=1))\n'
     '        out = {"codes": a.codes,',
     '                 pairs=np.stack([pm_rows, ps_rows], axis=1),\n'
     '                 **({"fp_G": FP.G} if (FP is not None and a.relp and a.relp_g) else {}))        # r65\n'
     '        out = {"codes": a.codes,'),
    # 5. the state
    ('            extra["fp_place"] = FP.place\n',
     '            extra["fp_place"] = FP.place\n'
     '            if a.relp and a.relp_g:     # r65\n'
     '                extra["fp_G"] = FP.G\n'),
]
for i, (old, new) in enumerate(patches):
    assert src.count(old) == 1, (i, src.count(old), old[:80])
    src = src.replace(old, new)
head = ('"""r65_chem_all15_frelg = chem_t10_frel.py (T10, copied whole) + --relp-g (R65-CHEM-ALL15, second hold; written by\n'
        'tools/r65_chem_all15_mkfrelg.py: read that file). Lines of this patch are marked "r65".\n\n')
i0 = src.index('"""')
src = src[:i0] + head + src[i0 + 3:]
open(os.path.join(T, "r65_chem_all15_frelg.py"), "w", encoding="utf-8", newline="\n").write(src)
print("written tools/r65_chem_all15_frelg.py:", len(src), "bytes;", len(patches), "patches, each matched once")
