"""T10 seat 4: writes tools/chem_t10_frelgcj.py = tools/chem_t10_frelgc.py (seat 3's relation tool, md5 9bd8a225)
+ seat 2's seeded freeze-order jitter (tools/chem_t10_mkgj.py, moved to the relation tool):
  --order-jitter EPS --order-seed S : every word's ranking value x exp(EPS x N(0,1)), one field per run (seeded by S,
  drawn at the first ranking with that ranking's shape: identical on a resume); EPS 0 = the source tool exactly.
    python tools/chem_t10_mkfrelgcj.py tools/chem_t10_frelgc.py tools/chem_t10_frelgcj.py"""
import hashlib
import sys

src = open(sys.argv[1], encoding="utf-8").read()
out = sys.argv[2]


def sub(old, new):
    global src
    n = src.count(old)
    assert n == 1, (n, old[:70])
    src = src.replace(old, new)


# 1. the options (before the R8 lever options)
sub('    ap.add_argument("--lm-until",',
    '    ap.add_argument("--order-jitter", type=float, default=0.0, help="t10: freeze-order jitter: ranking value x '
    'exp(EPS N(0,1)), seeded, fixed per word (0 = off)")\n'
    '    ap.add_argument("--order-seed", type=int, default=0, help="t10: seed of --order-jitter")\n'
    '    ap.add_argument("--lm-until",')
# 2. the field holder, after the sdyn set-up
sub('        sdyn_job, sdyn_from = None, None\n    COND = a.rank_cond > 0',
    '        sdyn_job, sdyn_from = None, None\n'
    '    JIT_ = [None]      # t10 seat 4: the seeded jitter field, drawn at the first ranking (same seed + shape on a resume)\n'
    '    COND = a.rank_cond > 0')
# 3. the ranking value, after the rank family and before the r47 weight division
sub('            if a.rank_pow:      # r47: the row weight taken out of the ranking value',
    '            if a.order_jitter > 0:      # t10: R67-CHEM-PATH\'s skkj order in the relation tool (ranking value only)\n'
    '                if JIT_[0] is None or JIT_[0].shape != np.shape(impact):\n'
    '                    JIT_[0] = np.exp(a.order_jitter * np.random.default_rng(a.order_seed).standard_normal(np.shape(impact)))\n'
    '                    log(f"t10: freeze-order jitter EPS {a.order_jitter} seed {a.order_seed} on {JIT_[0].shape}: field md5 "\n'
    '                        f"{hashlib.md5(JIT_[0].tobytes()).hexdigest()[:8]}, log-sd {float(np.log(JIT_[0]).std()):.4f}")\n'
    '                impact = impact * JIT_[0]\n'
    '            if a.rank_pow:      # r47: the row weight taken out of the ranking value')
if "\nimport hashlib" not in src and "import hashlib\n" not in src:
    src = src.replace("\nimport numpy as np", "\nimport hashlib\nimport numpy as np", 1)
src = src.replace('"""', '"""T10 seat 4 copy (tools/chem_t10_mkfrelgcj.py): + --order-jitter / --order-seed.\n', 1)
open(out, "w", encoding="utf-8", newline="\n").write(src)
print("written", out, hashlib.md5(src.encode()).hexdigest()[:8])
