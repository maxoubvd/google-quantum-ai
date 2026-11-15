# Developing on a Mac? A native kickmix build that runs the official verifiers locally (Apple Silicon, open source)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752659 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-25 19:57 UTC · 5 votes

If you work on a Mac you have probably hit the same wall I did: the competition's
`kickmix/kickmix.so` is an x86-64 Linux, CPython 3.13 binary, so every verifier run
meant a Kaggle CPU notebook round trip. For me that was about 19 minutes for a
full three-task check, plus queue time, and it capped how many candidates I could
try in a day.

kickmix is open source inside Google's Qualtran repo, and it builds natively on
macOS. I packaged that build, with NEON SIMD word types for Apple Silicon, as a
small public repo:

**https://github.com/manwithacat/kickmix-macos** (Apache-2.0)

```bash
git clone https://github.com/manwithacat/kickmix-macos.git
cd kickmix-macos
tools/build.sh          # ~2 min: creates .venv, builds kickmix, runs the 3,403 upstream tests
.venv/bin/pip install scipy   # the verifiers need it
```

Then, with the competition data unpacked somewhere (`verifiers/` and
`sample_submission/` from the official zip), the **unmodified** official verifiers
run locally:

```python
import sys; sys.path.insert(0, "verifiers")
import kickmix as km
from invsqrt.verify_submission_invsqrt import ProblemSpecInvSqrt, verify
circuit = km.Circuit.from_file("sample_submission/invsqrt.kmb", format="kmb")
print(verify(circuit, ProblemSpecInvSqrt.for_n_p(10), max_workers=8)["score"])
```

The DQI and chemistry verifiers work the same way (`dqi.verify_submission_dqi`,
`dfthcsos.verify`). DQI needs about 6 GB per worker, so use `max_workers=2` on a
16 GB machine.

**Does it match Kaggle?** On the reference circuits, yes, to the last digit:
invsqrt 2,076,774.233715282 with Q 260 / T 7,735; DQI 24,439,463,583,672.93 with
the identical sampled mean Toffoli count (so the seeded shots are the same as the
shipped binary's); chemistry 3,355,599,588.3155 with Q 730 / T 93,237. A
third-party portfolio I verified locally at 11.38437 then reproduced the same
costs and the same ZIP hash when run as a Kaggle notebook. The caveat that always
applies: the leaderboard uses a secret salt, so DQI's mean and repetition factor
can move slightly; test on two or more local salts before trusting a DQI gain.

**Speed** on a 10-core Apple Silicon Mac: full invsqrt verification 160 s with 8
workers (versus ~860 s hosted), DQI 158 s, chemistry 41 s. The NEON types give
1.9× on the invsqrt circuit and 1.16× on DQI over the portable build.

**Honesty notes:** the NEON additions were written with an AI coding assistant
under my direction, which is why this is a separate repo rather than an upstream
PR (Qualtran's CLA excludes AI-generated code). Everything else is Google's code,
unmodified, with their license and notices intact. I am not affiliated with
Google, Qualtran or Kaggle. Intel Macs and Linux build the same source with the
upstream code paths. If you find a discrepancy against the hosted verifiers,
please open an issue; the repo's CI runs the upstream Python and C++ suites plus a
NEON-versus-portable parity check on every push.

Credit where it is due: the verifiers and kickmix are Google Quantum AI's work;
the portfolio I used to check parity is other participants' public circuits
(Phuc Thai, anazemcev, Dmitrii Gluzdov and the lineages they credit).


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**hsiaosuan** · 2026-10-26 01:02 UTC · 0 votes

Kudos.
my claude build one for  silicon mac for me last friday

**SpeedSci** · 2026-10-26 04:56 UTC · 0 votes

Very useful!
