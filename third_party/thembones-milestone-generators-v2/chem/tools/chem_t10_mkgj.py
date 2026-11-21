"""T10: writes tools/chem_t10_froundgj.py = R65's tools/r65_chem_all15_froundg.py (G10g's tool, read-only, copied)
+ a seeded freeze-order jitter (R67-CHEM-PATH's simulator policy 'skkj' moved into the real rounding):
  --order-jitter EPS --order-seed S : every word's ranking value x exp(EPS x N(0,1)), one field per run (seeded,
  identical on a resume); the placer's costs (fpc) and everything else untouched. EPS 0 = the source tool exactly."""
import hashlib, sys
src = open(sys.argv[1], encoding="utf-8").read()
out = sys.argv[2]

def sub(old, new):
    global src
    n = src.count(old)
    assert n == 1, (n, old[:60])
    src = src.replace(old, new)

# 1. the options (before the R8 lever options)
sub('    ap.add_argument("--lm-until",',
    '    ap.add_argument("--order-jitter", type=float, default=0.0, help="t10: freeze-order jitter: ranking value x '
    'exp(EPS N(0,1)), seeded, fixed per word (0 = off)")\n'
    '    ap.add_argument("--order-seed", type=int, default=0, help="t10: seed of --order-jitter")\n'
    '    ap.add_argument("--lm-until",')
# 2. the field, drawn once after the sdyn set-up (outside it)
sub('        sdyn_job, sdyn_from = None, None\n\n    def rules():',
    '        sdyn_job, sdyn_from = None, None\n'
    '    JIT_ = None\n'
    '    if a.order_jitter > 0:      # t10: R67-CHEM-PATH\'s skkj order in the real tool (seeded; the same field on a resume)\n'
    '        import hashlib as hl_\n'
    '        JIT_ = np.exp(a.order_jitter * np.random.default_rng(a.order_seed).standard_normal(theta.shape))\n'
    '        log(f"t10: freeze-order jitter EPS {a.order_jitter} seed {a.order_seed} on {JIT_.shape}: field md5 "\n'
    '            f"{hl_.md5(JIT_.tobytes()).hexdigest()[:8]}, log-sd {float(np.log(JIT_).std()):.4f}")\n'
    '\n    def rules():')
# 3. the ranking value, after the rank family and before the r47 weight division
sub('            if a.rank_pow:      # r47',
    '            if JIT_ is not None:        # t10: the seeded jitter of the freeze order (ranking value only)\n'
    '                impact = impact * JIT_\n'
    '            if a.rank_pow:      # r47')
src = src.replace('"""', '"""T10 copy (tools/chem_t10_mkgj.py): + --order-jitter / --order-seed.\n', 1)
open(out, "w", encoding="utf-8", newline="\n").write(src)
print("written", out, hashlib.md5(src.encode()).hexdigest()[:8])
