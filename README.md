# Google Quantum AI — Verifiable Quantum Circuit Compilation

Workspace for the Google Quantum AI Kaggle competition on **verifiable quantum circuit
compilation**. It gathers the official competition data (the `kickmix` simulator, the verifiers,
sample submissions), the competition instructions and forum, analyses, and the reference
scientific papers.

The competition is based on the paper *Leveraging verifiable structure in quantum circuits to find
more efficient compilations with agentic AI* (Khattar, Rubin, Gidney, Broughton, Low, Babbush —
arXiv:2610.23832), included as
[`google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf`](google-verifiable-quantum-circuit-compilation/verifiable_circuit_compilation.pdf).

## How it works

Every task asks for a **kickmix circuit**: reversible classical gates (X, CX, Toffoli), phase gates
and measurement-based uncomputation (MBUC). Such a circuit maps every computational-basis state to
a single basis state times a phase, so a classical simulator can run it very fast on random inputs
and check it exactly (fuzz testing). Test inputs are seeded from the SHA-256 hash of the circuit
file (Fiat-Shamir heuristic), so there is no hidden test set.

## The three tasks

| Task | What the circuit does | Cost (lower is better) | Reference |
|---|---|---|---|
| **invsqrt** (Benchmark I) | Applies the phase `exp(2πi·β/√x)` to a 22-bit input `x` (Coulomb energy in inertial-fusion simulations) | `Q × (Toffolis + T/4)` | 260 qubits, 7,988 |
| **dqi** (Benchmark II) | Decodes a Reed-Solomon syndrome over GF(2¹²) and applies a sign (Decoded Quantum Interferometry) | `R(ε) × Q² × Toffolis` | 1,893 qubits, 5.73 M Toffolis |
| **dfthcsos** (Benchmark III) | Loads, with QROM circuits, the rotation angles of a DF-THC-factorized FeMoCo Hamiltonian (chemistry) | `λ_eff × Toffolis × Q` | 730 qubits, 93 k Toffolis |

**Score:** for each task, `reference cost / submitted cost`; the leaderboard shows the
**geometric mean** of the three ratios. Submitting the reference circuits scores 1.0.

## Repository layout

```text
.
├── README.md
├── CLAUDE.md                           Guidance for Claude Code in this repository
├── requirements.txt                    Python dependencies (pinned)
├── Analysis/                           Third-party analyses of the competition (PDF)
├── Instructions/                       Competition instructions copied from Kaggle (overview, data, rules, discussion)
├── Scientific Papers/                  Reference papers (PDF)
└── google-verifiable-quantum-circuit-compilation/   Official competition data
    ├── verifiable_circuit_compilation.pdf          Paper describing the competition
    ├── kickmix/                         kickmix simulator (compiled Python module, Linux x86-64 only)
    ├── sample_submission/
    │   ├── invsqrt.kmb                  invsqrt reference circuit
    │   ├── dqi.kmb.xz                   dqi reference circuit (compressed, see below)
    │   └── dfthcsos/                    Chemistry submission: angle_k_0..9.kmb + submission.npz
    └── verifiers/
        ├── invsqrt/                     Verifier + reference circuit (Benchmark I)
        ├── dqi/                         Verifier + reference decoder (Benchmark II)
        └── dfthcsos/                    Chemistry verifier (Benchmark III), see SCORING_README.md
```

### Papers (`Scientific Papers/`)

- `2510.10967v1.pdf`: MBUC Reed-Solomon decoder for DQI (used by the `dqi` verifier)
- `reiher-et-al-2017-…pdf`: elucidating reaction mechanisms on quantum computers (FeMoCo)
- `rubin-et-al-2024-…pdf`: quantum computation of stopping power for inertial fusion target design
- `Babbush et al..pdf`, `Khattar et al..pdf`, `Nature 2025.pdf`, `PRX 2025.pdf`

## Setup

Requirements: Python 3.10–3.14 on Windows, Linux or macOS. `requirements.txt` installs
[`kickmix` from PyPI](https://pypi.org/project/kickmix/) (prebuilt wheels for all three systems),
`numpy` and `scipy` for the verifiers, and `pandas` and `ipykernel` for notebooks.

```bash
git clone https://github.com/maxoubvd/google-quantum-ai.git
cd google-quantum-ai
python -m venv .venv                        # Windows: py -3.12 -m venv .venv
source .venv/bin/activate                   # Windows (PowerShell): .venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

The bundled `kickmix/kickmix.so` is the Linux x86-64 build from the Kaggle notebook image and does
not load on Windows or macOS; the PyPI package replaces it.

### Windows: run the verifiers in WSL

The verifiers work on Windows but are much slower there:

- **Single-core invsqrt and dqi checks:** both verifiers run in parallel only when the OS supports
  `fork`, which Windows doesn't, so they fall back to one core.
- **Slower simulation:** the Windows kickmix wheel also simulates about 4× slower per core.

On a 12-core laptop:

| Check | Windows | WSL (Ubuntu) |
|---|---|---|
| Sample chemistry submission | 5 min 9 s | 1 min 6 s |
| invsqrt, `--shots-per-x 1` | 23 min 25 s | 54 s |

Keep the repository on Windows and run the verifiers in a Linux environment inside WSL:

```powershell
wsl --install Ubuntu        # once; then run `wsl -d Ubuntu` to create the Linux user
```

```bash
# Inside Ubuntu (the Windows drive is mounted under /mnt/c):
sudo apt-get install -y python3-venv
python3 -m venv ~/venvs/gqai
~/venvs/gqai/bin/pip install -r /mnt/c/Users/<you>/.../google-quantum-ai/requirements.txt
```

Keep the venv in the Linux filesystem (`~`) rather than on `/mnt/c`, where file access is slow. WSL
gets up to half of the machine's RAM by default. The dqi verifier can need about 9 GB with 4
workers, so lower `--max-workers` or raise the limit in `%UserProfile%\.wslconfig`.

### Large files

Two files exceed GitHub's 100 MB limit. They are versioned as compressed `.xz` archives, and the
originals are listed in `.gitignore`. **Decompress them before running the verifiers** (`-k`
keeps the archive):

```bash
cd google-verifiable-quantum-circuit-compilation
xz -dk verifiers/dfthcsos/FeMoCo76o113e.fcidump.xz   # 35 MB -> 123 MB, read by the chemistry verifier
xz -dk sample_submission/dqi.kmb.xz                  # 3.7 MB -> 846 MB, dqi reference circuit
```

The `dfthcsos` verifier checks the SHA-256 of the FCIDUMP
(`e5a48af746fd075026dec9e8f5df2d363802a121e8c53af9a0b984bde695e1bd`), so a badly decompressed
file is rejected.

## Running the verifiers

Run the commands from `verifiers/`, with the virtual environment active (on Windows, inside WSL; see
above). To use the bundled Linux build instead of the PyPI package, also run
`export PYTHONPATH=../kickmix:$PYTHONPATH`.

```bash
cd google-verifiable-quantum-circuit-compilation/verifiers

# Benchmark I — invsqrt (200 M shots by default, which is slow; --shots-per-x lowers the count)
python invsqrt/verify_submission_invsqrt.py ../sample_submission/invsqrt.kmb

# Benchmark II — dqi (102,400 shots)
python dqi/verify_submission_dqi.py ../sample_submission/dqi.kmb

# Benchmark III — chemistry (about 25 s, expected score 3355599588.3155365)
python dfthcsos/verify.py ../sample_submission/dfthcsos
```

Each verifier prints a JSON summary (status, cost, Toffoli count, qubit count, etc.). The
`--salt STR` option changes the sampling seed; without it, results are publicly reproducible. For
quicker test runs, lower the shot count with `--shots-per-x` (invsqrt), `--num-shots` (dqi, a
multiple of 256) or `--shots-per-address` (dfthcsos). For invsqrt, the total shot count must be a
multiple of `--batch-size`; with `--shots-per-x 1`, use `--batch-size 12`. The
submission format and the details of the chemistry score are described in
[`verifiers/dfthcsos/SCORING_README.md`](google-verifiable-quantum-circuit-compilation/verifiers/dfthcsos/SCORING_README.md).

To regenerate the reference circuits:

```bash
python invsqrt/nr_reference_circuit.py invsqrt.kmb --format kmb
python dqi/dqi_reference_circuit.py decoder.kmb --format kmb
python dfthcsos/example_submission/angle_k_reference_circuit.py \
    ../sample_submission/dfthcsos/submission.npz OUT_DIR --brot 12 --batch-size 8
```

## Licenses

The `kickmix` simulator and the verifiers are released under the **Apache 2.0** license
(`kickmix/LICENSE.txt`, `verifiers/LICENSE.txt`). Reused public circuits keep their authors'
licenses. The papers in `Scientific Papers/` remain
the property of their authors and publishers.
