# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Rules

- **English for everything written to files.** Code, comments, docstrings, and analysis or report
  documents are always in English, even when the user writes in French.
  Chat replies may be in French.
- **Never commit or push.** Do not run `git commit`, `git push`, or any other command that writes
  to git history or a remote (including `gh` commands that create commits, branches, PRs or
  releases). The user handles all git operations themselves.
- **Read the changelog at the start of every session.** Before any other work, read
  `changelog/README.md` and the most recent entries in `changelog/` (at least the last five, more if
  they belong to the step in progress), then the step table in `Analysis/PLAN.md`, so you know what
  was done, why, and which step comes next.
- **Write a changelog entry for every change to code or files.** Whenever you write or modify code,
  circuits, tools, notebooks, the plan or other project files, add a **new** file to `changelog/`
  (never edit an old entry; correct it with a new one). Name it
  `changelog/YYYY-MM-DD-HHMM-short-slug.md` (UTC) and follow the template in `changelog/README.md`:
  what was done, why, and the results or analysis (measured costs, verifier output, nulls, open
  questions). One entry per plan step or per coherent change, written when it is finished.
- **Optimize run time, never at the cost of correctness or result quality.** Code we write
  (generators, tools, search loops, checkers) must run as fast as reasonably possible, because the
  speed of the try-and-verify loop limits how many ideas we can test. Speed never justifies a less
  exact result: an optimized tool must give the same output as the reference method (same circuit
  bytes, same costs, same verdicts), checked on known cases before it replaces it; a fast proxy
  (static cost, reduced run, sampling) only screens candidates and never replaces the full check at a
  gate; no search is cut short or approximated if that could miss a better circuit, unless the
  trade-off is stated and agreed. Within that limit:
  - choose the algorithm with the best complexity first (state the time and memory complexity of
    hot paths in their docstring), and avoid quadratic or worse scans when a linear or log-linear
    method exists;
  - vectorize with NumPy (or bit-parallel integer tricks) instead of Python loops on hot paths,
    and cache or memoize results that are reused;
  - use the 12 cores: `multiprocessing` with `fork` in WSL; size worker counts to the 8 GB WSL
    memory cap;
  - profile before tuning (`cProfile`, `time.perf_counter`), and record the measured run time and
    peak memory of every tool or build in its changelog entry;
  - iterate on the fastest check that answers the question (static cost, reduced verifier run, the
    generators' own exact proofs); run the full official verifiers only at gates (before a
    submission or a final choice).
- **Work step by step.** `Analysis/PLAN.md` splits the work into steps (`S<phase>.<n>`). Work on one
  step at a time, update its status in the plan when it ends, and end it with a short review in the
  chat before starting the next one.

## What this repository is

A workspace for the Kaggle competition *Google Quantum AI — Verifiable Quantum Circuit Compilation*
(final deadline 2027-01-20 23:59 UTC, 5 submissions/day). The goal is to produce cheaper
**kickmix circuits** for three tasks; there is no application to build. The repo has no build
system, no linter and no test suite: the official **verifiers are the test harness**, and a circuit
is "correct" when its verifier prints `"status": "PASS"`.

- `google-verifiable-quantum-circuit-compilation/`: official competition data, as the leaderboard
  runs it. `verifiers/` is the authoritative specification of each task. Treat this tree, including
  `sample_submission/`, as read-only reference material; write new circuits and scripts elsewhere.
- `Instructions/`: snapshot (2026-11-15) of the Kaggle Overview/Data/Rules tabs and all forum topics.
  `Instructions/discussion/README.md` indexes the topics; official answers are marked **Host**.
  The forum holds most of the known techniques, public frontier scores and "what didn't work" lists.
- `Notebooks/` (versioned): our own notebooks, written here and usually run on Kaggle (CPU, or GPU
  for float64-heavy chemistry refits). Any notebook that reuses a public circuit pins it by SHA-256
  and credits its author and license.
- `Analysis/PLAN.md`: our action plan and its phases. The dates are indicative and move with progress.
- `Analysis/`: third-party analyses of the competition (PDF).
- `Scientific Papers/`: the background papers. The competition paper itself is
  `google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf`.

## Researching competitors and the state of the art

Whenever the user asks what other competitors did, for new ideas, or for an audit comparing our
results with other teams or with state-of-the-art methods, start from the live list of public
Kaggle notebooks for this competition:
https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/code

It is more current than the local snapshots. This link is only the starting point, not the only
source: continue with any other relevant source, such as the live forum, the forum snapshot in
`Instructions/discussion/`, the notebooks in `Notebooks/`, `Analysis/`,
`Scientific Papers/`, and papers or code found on the web. Record the author, license and SHA-256 of any public circuit worth reusing.

## Environment setup

The user works on Windows, with two Python environments built from the same pinned `requirements.txt`
(kickmix comes from PyPI in both; the bundled `kickmix/kickmix.so` is a Linux build for the Kaggle
image and is not used):

| Environment | Python | Use it for |
|---|---|---|
| **WSL Ubuntu**, venv `~/venvs/gqai` (Linux filesystem) | 3.14 | **All verifier runs.** Matches the leaderboard's Linux setup |
| Windows, `.venv/` at the repo root | 3.12 | Editing, notebooks, quick scripts that don't simulate much |

Why verifiers must run in WSL: the Windows kickmix wheel simulates about 4.4× slower per core, and
the invsqrt and dqi verifiers only parallelize when the OS supports `fork`
(`invsqrt/verify_submission_invsqrt.py:221`), so on Windows they silently run on 1 of 12 cores.
Measured on this laptop: the sample chemistry check takes 1 min 6 s in WSL against 5 min 9 s on
Windows. A reduced invsqrt run (`--shots-per-x 1`) takes 54 s in WSL against 23 min 25 s on Windows.

Memory: the laptop has 16 GB, and WSL gets up to 8 GB by default. The official dqi verifier can
need about 9 GB with 4 workers, so use `--max-workers 2` for dqi, or ask the user to raise the cap in
`C:\Users\maxim\.wslconfig`.

To recreate the WSL venv: `python3 -m venv ~/venvs/gqai && ~/venvs/gqai/bin/pip install -r requirements.txt`.
It needs `python3-venv`, which can be installed without a password via `wsl -d Ubuntu -u root -- apt-get install -y python3-venv`.
- Two large files are versioned only as `.xz`. Decompress them once before verifying (`-k` keeps the
  archive; `xz` is available in Git Bash):
  ```bash
  cd google-verifiable-quantum-circuit-compilation
  xz -dk verifiers/dfthcsos/FeMoCo76o113e.fcidump.xz   # required by the dfthcsos verifier (SHA-256 checked)
  xz -dk sample_submission/dqi.kmb.xz                  # 846 MB reference dqi circuit
  ```
  The decompressed originals are git-ignored. Never commit them, because both exceed GitHub's 100 MB limit.

## Running the verifiers

Run all verifier commands from `google-verifiable-quantum-circuit-compilation/verifiers/`, inside
WSL. Inside WSL, the repo is at `/mnt/c/Users/maxim/Documents/Projets/google-quantum-ai`.

From Claude Code's Bash tool (Git Bash), call WSL like this. `MSYS_NO_PATHCONV=1` stops Git Bash
from rewriting the `/mnt/c/...` path, and single quotes keep the inner command intact. From the
PowerShell tool, nested quotes get mangled; put multi-step commands in a script file and run
`wsl.exe -d Ubuntu -- bash /mnt/c/.../script.sh`.

```bash
MSYS_NO_PATHCONV=1 wsl.exe -d Ubuntu \
  --cd /mnt/c/Users/maxim/Documents/Projets/google-quantum-ai/google-verifiable-quantum-circuit-compilation/verifiers \
  -- bash -c '~/venvs/gqai/bin/python dfthcsos/verify.py ../sample_submission/dfthcsos'
```

The verifier commands, with `python` meaning `~/venvs/gqai/bin/python`:

```bash
python invsqrt/verify_submission_invsqrt.py CIRCUIT.kmb   # 64 shots x 3.1M inputs, about 1 h in WSL
python dqi/verify_submission_dqi.py CIRCUIT.kmb           # 102,400 shots
python dfthcsos/verify.py RUN_DIR/                        # angle_k_*.kmb + submission.npz, about 1 min in WSL
```

For a reduced invsqrt run, the total shot count must be a multiple of `--batch-size`.
3,133,452 = 12 × 261,121, so `--shots-per-x 1 --batch-size 12` works; the default batch size of 256
only divides the shot count when `--shots-per-x` is a multiple of 64.

Flags for fast iteration (a reduced run is a smoke test, not a proof):
- invsqrt: `--shots-per-x N`, `--batch-size`, `--max-workers`
- dqi: `--num-shots N` (must be a multiple of `--batch-size`, default 256), `--max-workers` (default 4)
- dfthcsos: `--shots-per-address N` (default 64)
- all three: `--salt STR` and `--seed N`; `--format kmb|kmx` (auto-detected by default)

Each verifier prints a JSON summary with status, cost, Toffolis and qubits. The baseline dfthcsos
score is `3355599588.3155365`, which is a quick check that the setup works.

The scripts that write the reference circuits are also the best examples of the kickmix builder
API (`km.CircuitBuilder`, `create_quantum_register(..., name=)`, `alloc_qubits`, `free`, `mark`,
`z_pow`, `finish_circuit`, `Circuit.to_file`, `operation_statistics().num_toffoli_like`,
`km.Simulator`, `km.GF2Field`, unary iteration):
`invsqrt/nr_reference_circuit.py`, `dqi/dqi_reference_circuit.py` and
`dfthcsos/example_submission/angle_k_reference_circuit.py`.

## How scoring works (big picture)

Score = geometric mean over the three tasks of `reference_cost / submitted_cost`. The reference
circuits score 1.0. Every task uses kickmix circuits made of X/CX/Toffoli, phases and
measurement-based uncomputation (`HMR` plus conditioned fixups). Each basis state maps to one basis
state times a phase, so the verifiers check circuits exactly by classical fuzzing.

| Task | Register contract | Cost (lower is better) |
|---|---|---|
| `invsqrt` | One register `x` (22 qubits); all other qubits are scratch and must return to `0`. Phase `1/sqrt(x)` turns within `2^-24/sqrt(x)`. | `Q × (Toffolis + T/4)`; `Z_POW` rotations are priced in T, and the error budget left over sets their synthesis precision |
| `dqi` | Quantum `S0..S69` and classical `E0..E4094`, 12 bits each; no other registers; no `Z_POW` | `R(ε₊) × Q² × Toffolis`; failures are allowed but raise the Clopper–Pearson bound `ε₊`, which must stay below 1/16 |
| `dfthcsos` | Per circuit: `selection` (exactly `ceil(log2(76+R·B))` bits) and `word_0..word_{J-1}` (`brot` bits each); no `Z_POW`; each circuit ≤ 128 qubits | `λ_eff × T_total × Q_total` over the **whole block encoding**, not just the ANGLE circuits |

Things that are easy to get wrong:
- **Seeding:** test inputs and `HMR` outcomes are seeded from the SHA-256 of the circuit's canonical
  `.kmb` bytes, plus a salt. Any edit to a circuit therefore gets fresh tests. A local run uses an
  empty salt; the leaderboard uses a secret salt for each split. A local PASS on a circuit that fails
  on rare measurement branches can still fail on the leaderboard.
- **Toffoli counts are measured, not static.** The verifier counts the mean of `CCX`/`CCZ`/`CSWAP`
  gates executed per shot. A gate conditioned on a measured bit counts only on the shots where that
  bit is set.
- **dfthcsos is mostly outside the submitted circuits.** `T_total` and `Q_total` add a heuristic
  cost (`scoreinfra/costs_prepare.py`, `SpinFreeData_V1`) that depends on the factorization shape
  `R, B, C` and on `brot`. Temporary workspace is shared (`max`, not `sum`) across the state
  preparations and the ANGLE circuits. `λ_eff` comes from the quantized words the circuits load. A
  submission is rejected unless `Λ > 49.0` and `E_gap > 3.914898` (strict inequalities) and both
  errors are ≤ 20 mHa. `verifiers/dfthcsos/SCORING_README.md` is the full contract.
- **The metric rejects things the standalone verifiers accept** (forum topic 751951):
  - every circuit's `max(num_qubits, num_bits) - 1` must be < 699,050;
  - the zip must contain exactly one `invsqrt`, one `dqi` and one `dfthcsos/`, searched recursively
    (a leftover backup copy anywhere counts);
  - `angle_k_0..n-1` must have no gaps;
  - every register must be named;
  - the `submission.npz` arrays must be little-endian floats (`tree_edges` an integer type).
- The submission is `submission.zip` containing `invsqrt.kmb`, `dqi.kmb` and
  `dfthcsos/{angle_k_*.kmb, submission.npz}`. `.kmx` text circuits are accepted too.

## Licensing of reused circuits

The verifiers and kickmix are Apache 2.0. Under the rules (`Instructions/03_rules.md`), winners
license their submission and the source code that generates it under Apache 2.0. They must also
provide a code repository from which the results can be reproduced. When reusing a public circuit or notebook, record its author, license and SHA-256, as the
existing notebooks do. Several forum topics discuss which public DQI decoders are prize-eligible.
