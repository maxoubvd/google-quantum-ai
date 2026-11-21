"""T10 seat 2: writes tools/chem_t10_froundgjc.py = tools/chem_t10_froundgj.py (R65's froundg + my seeded order
jitter) + R67-CHEM-PATH's CONDITIONAL FREEZE ORDER copied from tools/r67_chem_path_round.py (read-only): the options
--rank-cond RHO / --rank-geo A / --cond-cpus LIST, the group-worker pool after the sdyn set-up, and the ranking value
in the sdyn branch (the latest finished conditional table; module tools/r67_chem_path_cond.py). Without --rank-cond
the output behaves as chem_t10_froundgj.py. The jitter (if on) multiplies the conditional value as it does S_kk."""
import hashlib, sys
src = open(sys.argv[1], encoding="utf-8").read()
r67 = open(sys.argv[2], encoding="utf-8").read().replace("\r\n", "\n")

def sub(old, new):
    global src
    assert src.count(old) == 1, (src.count(old), old[:70])
    src = src.replace(old, new)

def cut(text, start, end):          # the block of R67's file from the line starting with START to the line before END
    i = text.index(start)
    j = text.index(end, i)
    return text[i:j]

opts = cut(r67, '    ap.add_argument("--rank-cond"', '    ap.add_argument("--reset-free"')
setup = cut(r67, '    COND = a.rank_cond > 0', '    def rules():')
over = cut(r67, '                if COND:        # r67: the order by the conditional stiffness', '            elif a.rank == "sratio":')
assert opts.count("add_argument") == 3 and "cond_pool" in setup and "impact = val_c" in over
sub('    ap.add_argument("--lm-until",', opts + '    ap.add_argument("--lm-until",')
sub('        sdyn_job, sdyn_from = None, None\n    JIT_ = None\n',
    '        sdyn_job, sdyn_from = None, None\n' + setup.rstrip("\n") + '\n    JIT_ = None\n')
sub('                impact = rank_sk + 1e-30 * c\n            elif a.rank == "sratio":',
    '                impact = rank_sk + 1e-30 * c\n' + over + '            elif a.rank == "sratio":')
# seat 2 fix (03:45Z, measured on t10_g10k3 / k4: 'AssertionError: S_FF is not positive definite' in the last stage):
# a failed conditional table no longer stops the run: the last finished table stays in use (S_kk at a resume)
NL = chr(10)
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
src = src.replace('"""T10 copy (tools/chem_t10_mkgj.py)', '"""T10 copy (tools/chem_t10_mkgjc.py: + R67\'s --rank-cond) of T10 copy (tools/chem_t10_mkgj.py)', 1)
open(sys.argv[3], "w", encoding="utf-8", newline="\n").write(src)
print("written", sys.argv[3], hashlib.md5(src.encode()).hexdigest()[:8])
