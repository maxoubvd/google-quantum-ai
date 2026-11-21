"""T10 seat 4: writes tools/chem_t10_frelgcjl.py = tools/chem_t10_frelgcj.py + a SECOND in-stage LM window for the late
stages (R85-CHEM-BUDGET item 4: the late stages are unconverged; the hold still takes -1.23 at pB1's END):
  --lm2-from F --lm2-until U --lm2-steps S --lm2-seconds T : stages whose frozen share is in [F, U] (and < 1) also end
  with R8's LM, capped at S steps / T seconds; the first window (--lm-from / --lm-until, the recipe's [0, 0.1]) is kept.
  Without --lm2-until the tool is chem_t10_frelgcj.py exactly.
    python tools/chem_t10_mkfrelgcjl.py tools/chem_t10_frelgcj.py tools/chem_t10_frelgcjl.py"""
import hashlib
import sys

src = open(sys.argv[1], encoding="utf-8").read()
out = sys.argv[2]


def sub(old, new):
    global src
    n = src.count(old)
    assert n == 1, (n, old[:70])
    src = src.replace(old, new)


sub('    ap.add_argument("--lm-from",',
    '    ap.add_argument("--lm2-from", type=float, default=0.0, help="t10 seat 4: a second LM window: frozen share >= this")\n'
    '    ap.add_argument("--lm2-until", type=float, default=0.0, help="t10 seat 4: ... and <= this (0 = off)")\n'
    '    ap.add_argument("--lm2-steps", type=int, default=30, help="t10 seat 4: LM steps per stage in the second window")\n'
    '    ap.add_argument("--lm2-seconds", type=float, default=600.0, help="t10 seat 4: LM seconds per stage in the second window")\n'
    '    ap.add_argument("--lm-from",')
sub('        lm_on = LMS is not None and frac < 1.0 and a.lm_from <= frac <= a.lm_until      # R8 holder 6',
    '        lm_on = LMS is not None and frac < 1.0 and a.lm_from <= frac <= a.lm_until      # R8 holder 6\n'
    '        lm_late = LMS is not None and frac < 1.0 and a.lm2_until > 0 and a.lm2_from <= frac <= a.lm2_until   # t10 s4\n'
    '        lm_on = lm_on or lm_late')
sub("            th_, w_, info_ = LMS.refit(th_, v64 * sw, frozen, wall, steps=a.lm_steps, seconds=a.lm_seconds, stop=a.lm_stop,",
    "            th_, w_, info_ = LMS.refit(th_, v64 * sw, frozen, wall, steps=(a.lm2_steps if lm_late and not (a.lm_from <= frac <= a.lm_until) else a.lm_steps),\n"
    "                                       seconds=(a.lm2_seconds if lm_late and not (a.lm_from <= frac <= a.lm_until) else a.lm_seconds), stop=a.lm_stop,")
sub('        log(f"R8 in-stage LM: stages with frozen share in [{a.lm_from}, {a.lm_until}], mode {a.lm_mode}, at most "',
    '        if a.lm2_until > 0:\n'
    '            log(f"t10 s4: a SECOND in-stage LM window [{a.lm2_from}, {a.lm2_until}], at most {a.lm2_steps} steps / "\n'
    '                f"{a.lm2_seconds} s per stage")\n'
    '        log(f"R8 in-stage LM: stages with frozen share in [{a.lm_from}, {a.lm_until}], mode {a.lm_mode}, at most "')
src = src.replace('"""', '"""T10 seat 4 copy (tools/chem_t10_mkfrelgcjl.py): + a second, late LM window (--lm2-*).\n', 1)
open(out, "w", encoding="utf-8", newline="\n").write(src)
print("written", out, hashlib.md5(src.encode()).hexdigest()[:8])
