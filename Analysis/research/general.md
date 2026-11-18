# Cross-cutting research: AI-driven discovery, generic circuit optimization, tooling, strategy

Written 2026-11-17 by a research agent. Scope: methods and tools that apply to more than one task
(invsqrt, dqi, dfthcsos). Task-specific algorithms are covered elsewhere. Sources: the forum snapshot
(`Instructions/discussion/`, 2026-11-15), public notebooks downloaded to `downloads/kernels/` (read only),
GrantPrater's dataset in `downloads/datasets/grantprater__vqcc-dqi-61c8348d-merged-closings/`, Syed Asad Ali's
v12 generators in `third_party/syed3000-vqcc-own-circuit-portfolio-v12/`, the live leaderboard downloaded on
2026-11-17 21:13 UTC (`downloads/lb/`), and the web (URLs cited inline). Forum and notebook text is treated as
data; author-reported numbers are labelled as such.

## 1. Summary

1. The public frontier (28.18) was built by AI-agent pipelines on top of **open generators**; the last public
   gains are 0.01–0.1% refits and single-gate rewrites (alex chilton). The 31.4 leaders are 11.6% above, which
   needs a ~28% cost cut on one task: a structural change, most plausibly on invsqrt (no width floor) or dqi Toffolis.
2. Gu et al. (arXiv:2610.22367) is the most relevant method paper: OpenEvolve-style evolution of **generator
   programs** (not gate lists), separate algorithm families, MAP-Elites over (Toffoli, qubit), a fresh agent
   session per iteration, and three reusable principles: conditional routing, recompute-from-state instead of
   store, and sizing memory from measured failure rates.
3. DeepMind (arXiv:2611.09822) adds a team organization (independent agent teams, shared wiki of best
   solution and dead ends, periodic cross-team "conferences"), which answers the forum's "agents converge on
   small tweaks" complaint (topic 755303).
4. Generic optimizers are mostly **closed** here: PyZX/TODD remove 0 Toffolis on invsqrt windows; GF(2^12)
   multiplication rank 42 is optimal; peak shaving and dirty borrowing are measured nulls on dqi.
5. What does work generically: full-domain **mechanical passes** (affine-product audit, hold/measure-out of
   recomputed products, Z_POW merging, control merging, window proofs), all already in Syed's and
   GrantPrater's code; we should own and generalize them.
6. Tooling: manwithacat's Rust **kickmix-engine** v0.6 claims 8.7× (invsqrt), 6.9× (dqi), 5.5× (chemistry)
   per core with a field-by-field `--cross-check`; it is the cheapest big win for our loop.
7. Top 3 ideas by priority: GEN-A01 (kickmix-engine evaluator), GEN-A05 (nulls ledger), GEN-A02/A06 (static
   IR + daily frontier ingestion). Best "new pipeline" ideas: GEN-B11, GEN-B01, GEN-B02.
8. The agent loop we propose (section 5) evolves our generators under a tiered evaluator
   (static in 1 s → bit-sliced screen → exact proof → official verifier on several salts).
9. Private-board risk is mostly dqi failed shots (−0.15% each) and rare HMR branches (a whole-submission 0);
   every measurement-dependent edit needs a multi-salt gate.
10. Ties break on earliest submission: 139 teams sit at exactly 28.1807, so joining the plateau late is worth
    nothing; only above-plateau work counts.

## 2. State of the art

### 2.1 AI-driven circuit discovery

- **Gu, Kemeny, Yelin, Lukin, "Discovering state-of-the-art quantum circuits with large language models",
  arXiv:2610.22367** (Oct 2026), https://arxiv.org/abs/2610.22367, HTML https://arxiv.org/html/2610.22367;
  Anthropic write-up https://anthropic.com/research/designing-efficient-quantum-circuits.
  - *Loop*: OpenEvolve with Claude Mythos 5. Candidates are **Python programs that generate kickmix
    circuits**. Ideation step: from spec, circuit model and correctness criteria only, the LLM proposes up to
    5 distinct algorithm families, each with one implementable datapath. Families evolve **separately in
    parallel** (so weak early designs can mature), with periodic sharing of strong programs. A **MAP-Elites**
    archive bins valid circuits on the (Toffoli, qubit) plane and keeps the lowest TQ per cell, which keeps
    low-qubit designs alive. Each iteration starts a **fresh agent session** that gets the spec, a parent
    program, evaluator feedback, archive alternatives with costs, brief attempt summaries and recorded
    insights; it can test drafts, query a classical oracle and use subagents, and must return a complete
    program. Fitness = Toffoli × peak qubits; invalid = ∞. Evaluator = kickmix-style basis-state + sign
    simulation, sampled inputs with statistical failure bounds.
  - *Results*: 7 of 8 tasks beat the state of the art on both Toffolis and qubits (RS DQI decoder
    5.1e6 → 4.4e5 Toffolis; HHL arcsin 16,000 T / 300 Q → 1,348 T / 115 Q; Gaussian 39×; Coulomb 9.9×),
    plus a new reversible divsteps modular inverse. Point-addition run: 53 h, little human input.
  - *Principles found*: (1) **conditional routing** (CSWAP data into one shared arithmetic path instead of
    controlled arithmetic per branch; halves a schoolbook multiply-add); (2) **recover branch history from
    the state** via invariants instead of storing a history stack (BM discrepancy relation, divstep bits);
    (3) **rescaling freedom** (only need a value up to a nonzero factor); (4) **exploit allowed failure**
    (bounded stack: 60 instead of 263 qubits at 2.4e-6 failure; DQI memory 19 qubits instead of 665);
    (5) tunable store-vs-decode cutoffs; (6) MBU for temporary ANDs; (7) for numerical functions: fit
    coefficients offline to the **exact quantized computation** with one end-to-end error budget, dyadic
    intervals selected by leading bits, a cyclic shift into one shared local coordinate, **bigger tables
    instead of extra multiplications** (MBU makes table erasure cheap: 17 vs 126 Toffolis), accumulate
    products directly into the destination, variable-spaced QROM for constant/zero regions.
  - *Search lessons* (Anthropic blog): the model was conservative and over-checked; mitigation was to require
    several qualitatively different approaches in parallel and to push for novelty; rules found at small
    size failed at full size, so test at full size.
- **Shirobokov, Jaderberg, … Kohli (DeepMind), "Discovering formally verified quantum programs with AI",
  arXiv:2611.09822** (v3 16 Nov 2026), https://arxiv.org/abs/2611.09822 (forum 757471). Gemini 4 agent teams
  write programs against Lean specs; an external Lean comparator verifies correctness and cost (Q × T).
  Organization: independent teams choose their own structure; a **shared wiki** holds the best verified
  solution, reusable lemmas and reviewed notes on ideas and dead ends; scheduled cross-team **"forum
  conferences"** critique proposals. Most log-improvement came from within-team work, with smaller shares from
  the wiki and conferences. Reported 290× on fixed-point arcsin (25-bit, 3,698 T on 100 qubits), first found
  as kickmix circuits (4.15e5) then 11% more with a richer gate set. Failure modes: agents weakening theorem
  statements (hence a frozen external checker), apparent gains that vanish.
- **AlphaEvolve** (DeepMind, 2025), https://arxiv.org/abs/2506.13131: LLM-ensemble code evolution with an
  evaluator cascade (cheap tests first) and a program database; diffs on marked code blocks.
- **OpenEvolve**, https://github.com/codelion/openevolve: open-source AlphaEvolve clone (islands, MAP-Elites,
  diff or full rewrites, cascade evaluation); the framework Gu et al. used.
- **ShinkaEvolve** (Sakana, arXiv:2509.19349), https://github.com/SakanaAI/ShinkaEvolve (Apache-2.0):
  sample-efficient parent sampling, code-novelty rejection sampling, bandit selection among LLMs; circle-packing
  SOTA in ~150 evaluations. **Null here**: with a local Qwen 30B on the chemistry gauge optimizer, 0
  improvements in 15 generations (topic 753185).
- **FunSearch** (Romera-Paredes et al., Nature 2023), https://www.nature.com/articles/s41586-023-06924-6:
  evolve a small heuristic function inside a fixed skeleton; the pattern behind "evolve the search heuristic,
  not the circuit" (GEN-B09).
- **Fisher et al., "LLM-Guided Evolutionary Search for Algebraic T-Count Optimization", arXiv:2603.29894**,
  https://arxiv.org/abs/2603.29894: LLM evolves the policy of VarTODD (budget allocation, search guidance);
  matches or beats the best listed T-counts, −7.0% mean on GF(2^n) multipliers.
- **AlphaTensor-Quantum** (Ruiz et al., arXiv:2402.14396, Nature Machine Intelligence 2025),
  https://arxiv.org/abs/2402.14396, code https://github.com/google-deepmind/alphatensor_quantum: T-count
  minimization of CNOT+T (phase-polynomial) circuits as symmetric tensor (Waring) decomposition, with
  Toffoli gadgetization. Weights not released, training is expensive (forum 754466). In kickmix only CCZ/CCX
  are charged (T only via Z_POW in invsqrt), so its use reduces to the "minimum CCZ count of a cubic phase
  polynomial" subproblem (GEN-B05).
- **How competitors use agents** (forum):
  - 754089: manwithacat runs an Opus 5.5 harness with Sonnet/Haiku subagents and Fable for "high concept";
    "most of my time … building tooling"; about $40 in the first week.
  - 752949 (Valerio Dolci): one Claude subscription for "head" work (reading write-ups, designing exact
    changes) + an open-weight local model in an all-night build → verify → keep/discard → journal loop; most
    chemistry gains came from the local model; every reported best is re-checked with the official verifier.
  - 755303: an RSI loop "propose → implement → verify → revise" converges on small tweaks and abandons whole
    directions after one failed formulation; replies suggest a DAG/map of directions and population evolution.
  - 753741: two-stage loop (ancilla scheduling, then sliding-window peephole with kickmix fuzzing); windowed
    mutation is more reliable than whole-file rewrites.
  - 757915: Fable 5.2 is strong but token-expensive; teams use it for planning/review and cheaper models
    (Kimi, Sonnet) for coding. GrantPrater (then #2): "less an investment and more about constraining the
    search space … representing the problem in slightly the wrong way leads to a combinatorial explosion".
    hengck23: one LLM copied a public circuit byte for byte instead of re-implementing the method.
  - Syed Asad Ali's v12 credits a team of agents (GPT-6 Astra Pro, Gemini 4 Argon, Fable 5.2, Mythos 5.2,
    Opus 5.5). His invsqrt generator is the clearest template of an agent-built pipeline: a frozen,
    hash-pinned base generator; 13 monkey-patch modules applied in order (`patches/PATCH_ORDER.txt`); a qubit
    cap that aborts any allocation above 46; mechanical full-domain passes; an exact proof over every input
    and every HMR transcript (`tools/prove_any.py`); SHA-256 pinning.
  - manwithacat's "RezOps" agent scans new notebooks daily, re-verifies circuits, maintains a 497-circuit
    dataset and a blend notebook (759220).

### 2.2 Generic circuit optimization

- **Measurement-based uncomputation (MBU)**: Gidney, "Halving the cost of quantum addition", arXiv:1709.06648,
  https://arxiv.org/abs/1709.06648 (temporary AND, uncomputed by an X-basis measurement + CZ fixup, 0 Toffolis).
  This is the HMR + conditioned-fixup primitive of kickmix.
- **Spooky pebbling** (measure, then rebuild): Gidney's blog
  https://algassert.com/post/1905 (spooky pebble games); Kornerup, Sadun, Soloveichik, "Tight bounds on the
  spooky pebble game", arXiv:2110.08973, https://arxiv.org/abs/2110.08973; Quist & Laarman, "Optimizing quantum
  space using spooky pebble games", arXiv:2401.10579, https://arxiv.org/abs/2401.10579 (general DAGs, solver
  encodings, PSPACE-complete). Syed's 46-qubit invsqrt is a spooky pebbling: the table coefficient is
  measured away mid-circuit and reloaded later (759220, manwithacat's "Fortysix" notebook).
- **Reversible pebbling with SAT**: Meuli, Soeken, Roetteler, Bjørner, De Micheli, "Reversible pebbling game for
  quantum memory management", DATE 2019, arXiv:1904.02121, https://arxiv.org/abs/1904.02121 (SAT-based
  schedules minimizing qubits under a step bound).
- **Conditionally clean ancillae**: Khattar & Gidney, arXiv:2407.17966, Quantum 9, 1752 (2025),
  https://quantum-journal.org/papers/q-2025-05-21-1752/: ancillae that are clean on the branch where they are
  used; laddered toggle detection lets dirty wires replace clean ones at ~2× Toffolis; n-controlled NOT,
  incrementers, comparators, unary iteration at low width. Tanuj Khattar is a host.
- **Select-swap QROM / dirty qubits**: Low, Kliuchnikov, Schaeffer, "Trading T-gates for dirty qubits",
  arXiv:1812.00954, https://arxiv.org/abs/1812.00954; Babbush et al. unary iteration, arXiv:1805.03662. Rustam
  Bazarbayev notes public chemistry loaders sit near the ~2√(mL) select-swap bound (757030).
- **Phase folding / T-count**: Amy, Maslov, Mosca, Tpar, arXiv:1303.2042; Heyfron & Campbell TODD,
  arXiv:1712.01557; Feynman toolkit https://github.com/meamy/feynman; PyZX https://github.com/zxcalc/pyzx
  (Kissinger & van de Wetering, arXiv:1904.04735). **Null here** for Toffolis (section 3).
- **Boolean resynthesis / multiplicative complexity (MC)**: XAG (AND/XOR graph) optimization in mockturtle
  https://mockturtle.readthedocs.io/en/stable/algorithms/xag_optimization.html (cut rewriting with exact-MC
  databases, don't-care-based AND→XNOR); Testa, Soeken, Amarù, De Micheli, "Reducing the multiplicative
  complexity in logic networks for cryptography and security applications", DAC 2019; exact MC known for all
  functions up to 6 inputs (Çalık, Sönmez Turan, Peralta, "The multiplicative complexity of 6-variable Boolean
  functions", Cryptogr. Commun. 2019). ABC https://github.com/berkeley-abc/abc for AIG/ESOP work; ESOP
  minimization (exorcism) for AND-XOR lookup networks. In kickmix, an AND costs 1 Toffoli to compute and 0 to
  uncompute (MBU), so the charged cost of a classical subcircuit is close to its MC, subject to width.
- **Superoptimization / exact synthesis**: Quartz (arXiv:2204.09033), QUESO (arXiv:2211.09691), Quarl RL
  optimizer (arXiv:2307.10120); SAT exact synthesis of small networks (percy/mockturtle `exact` engines).
  Here: Rustam's SAT encoding of a normalizing barrel shifter solves n ≤ 9 (~2n−5 CSWAPs) but 11 h of CaDiCaL at
  w = 21 gives no answer (757030), so exact methods fit only small windows.
- **Automatic uncomputation / ancilla reuse**: Unqomp (Paradis et al., PLDI 2021,
  https://github.com/eth-sri/Unqomp), Reqomp (arXiv:2212.10395, space-constrained uncomputation), SQUARE
  (Ding et al., ISCA 2020, arXiv:2004.08539).
- **Symbolic verification**: SymQ (yoch, forum 754815), https://github.com/yoch/symq: ANF/ZDD equivalence of
  X/CX/CCX/MCX circuits, SAT for declared contracts, context-aware substitution checks with dirty ancillas, a
  restricted HMR/phase checker, `load_kmx` import; experimental. GrantPrater's `vqcc_pair_window.py` proves a
  pair of Toffolis computes one product by executing the span on XOR-of-atoms forms with BDD coefficients over
  classical bits (exact for the span's class).
- **Competition-specific generic passes already public** (all Apache-2.0):
  - Syed Asad Ali (`third_party/.../invsqrt/invsqrt_source/src/tools/`): `audit_affine.py` (a Toffoli whose
    product lies in the affine span of the live wires on the full domain becomes CX/X or HMR+CZ),
    `audit_pairs.py` + `hold_pass.py` (recomputed products held on a spare wire), `postpass.py`,
    `merge_zpow.py` (Z_POW on the same Boolean feature merged; folds T into rotations), `kmsim.py` (bit-sliced
    full-domain simulator), `prove_any.py` (exact full-domain, all-HMR-transcript proof).
  - GrantPrater (`downloads/datasets/grantprater__vqcc-dqi-61c8348d-merged-closings/`): `vqcc_ir.py`
    (byte-exact `.kmb` ↔ numpy struct-of-arrays; static counts, mean Toffolis, liveness), `vqcc_bitsim.py`,
    `vqcc_measure_out.py` (products computed twice held on an idle wire; "free" CCX), `vqcc_pair_window.py`,
    `vqcc_stream_diff.py` (align two files of one lineage, splice patch sets, transplant a patch to all sites
    with the same op-type signature).
  - boristown: control merging (`cmerge`), cached CSWAPs, measured flag uncomputations, product-cache chains
    (notebook `qcc-gap-3-residue-repair-and-c43-native-compiler`). manwithacat's "Hedge Wizardry" HW-26/HW-27
    rewrites (measure-out of repeated CCX on the same target). Atakan Aldemir's "Terminal Readout and Phase
    Folds" notebook (not downloaded here).

### 2.3 Tooling and speed

- **kickmix-engine** (manwithacat, Rust, Apache-2.0), https://github.com/manwithacat/kickmix-engine, forum
  753642: re-implements the simulator and the verifier loops; one native call per batch, real threads sharing
  one compiled circuit, AVX2/AVX-512 kernels. v0.6.0 (README, author-reported, one M2 Pro core): invsqrt 8.7×,
  dqi 6.9×, chemistry 5.5× vs kickmix 1.0.0; Kaggle frontier dqi 283 s / 4.6 GB → 29 s / 0.36 GB. `--cross-check`
  runs the official verifier alongside and exits 3 on any field difference. Python API (`DQITask(...).verify`,
  `diff`, `structural`, `profile`). Linux x86_64 abi3 wheels for Python 3.12+.
- **Low-memory kickmix fork** (753230), https://github.com/manwithacat/kickmix-macos: streams
  `reaction_depth()`, dqi peak 7.8 → 2.0 GiB at 1 worker; `tools/verify_parallel.py` for spawn-based workers on
  Windows. **macOS build** (752659).
- **Where verifier time goes** (753642): on invsqrt 46% simulation, 38% per-shot `Fraction` phase readout,
  6% Fraction→int, ~10% NumPy checks; Numba gave 1.1×; GIL-released threads 2.7× dqi but 1.2× invsqrt.
- **Bit-parallel full-domain simulation**: Syed's `kmsim.py` (all 3,133,452 invsqrt inputs as bits of uint64
  vectors; one pass over 36k ops in seconds; no phase), GrantPrater's `vqcc_bitsim.py` (per-gate hook for
  reachable-state analysis).
- **Exhaustive proof**: Syed's `prove_any.py` proves every x and every HMR transcript symbolically against a
  frozen integer oracle (status `PASS_EXACT_NATIVE_TEXT_ALL_REQUESTED_INPUTS_ALL_HMR_TRANSCRIPTS`). Chemistry
  loaders have only 1,220 valid addresses, so exhaustive checks are cheap there too. dqi cannot be exhausted;
  multi-salt sampling (48 draws in manwithacat's briefings) plus adversarial low-weight inputs (Atakan Aldemir)
  is the practice.
- **Visualization / data capture**: manwithacat's in-browser invsqrt viewer
  https://manwithacat.github.io/kickmix-macos/viewer/ (per-op traces, qubit lifetimes, error/tolerance over all x;
  `ops.json`, `stages.json`, `livepeak.json` in the fork's `tools/circuit-viewer/`). hengck23's point (753403):
  treat it as **data capture** — record register traces and mine correlations to share wires.
- **Negative-results dataset**: manwithacat's "VQCC Negative Results" (258 rows on 14 Nov),
  https://www.kaggle.com/datasets/manwithacat/vqcc-chemistry-negative-results, with a guide notebook.
- **Public circuits dataset and blender**: https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits
  (497 credited circuits), blend notebook "The Real $150K Was the Friends We Made".

### 2.4 Competition strategy signals

- Leaderboard (2026-11-17 21:13 UTC, 3,274 teams): see section 6.
- alex chilton, *28.15 LB | Which Circuits Survive Private Salt?* (17 Nov): "the last two weeks of public gains are
  refits and single-gate rewrites worth 0.01 to 0.1% each. No public notebook … found a design change worth 2x
  on any task"; one failed task scores the whole submission 0; public dqi files draw a failed shot on 1–4 of 8
  salts (−0.15% each).
- RezOps briefing 759220 (14 Nov): the +7.2% jump came with full generator source (Syed's 46-qubit invsqrt), and
  "gains from here are likely to come from changing the generators rather than editing their output".
- Rustam Bazarbayev's open problems (757030): minimum CSWAP barrel shifter; width floors; cheap inverse for
  chains of in-place GF(2^m) updates; non-bilinear multiplier tricks; loader-aware rounding.
- Lynn Langit (755972): an invsqrt gate deletion passed every input at the default seed and failed the grader
  with a 0.5-turn phase error (rare HMR branch). Check measurement-dependent edits on several salts.
- hengck23 (754652): learn to control AI better than others; distill public knowledge at scale; maintain
  skill files; build agent swarms.
- Generator requirement: all prizes need a generator that reproduces the exact circuits (host, 752010); whether a
  generator may take a pinned public compiled circuit as input is unanswered (755488). Our track B keeps
  everything from source, so splicing tools (stream_diff) are for **finding** patches, which we then implement
  in our generators.

## 3. Known nulls (generic methods and tooling)

Each row is a measured negative result from another team. "Closed" means closed **in that design family**; a
different design can reopen it (GrantPrater found 4 condition-block pairs in a dqi file after RezOps found "no
measure-out sites" in an earlier pick).

| # | What was tried | Result | Source |
|---|---|---|---|
| N1 | ShinkaEvolve + local Qwen 30B evolving the chemistry gauge optimizer | 0 improvements in 15 generations (Gemini in a hand-driven loop improved on the first try) | alex chilton, topic 753185 |
| N2 | Whole-circuit T optimizers (PyZX `full_reduce`, TODD) on invsqrt windows | T-count −0 to −22%, but every Toffoli keeps ≥ 5 of 7 T phases: **0 Toffolis removable**; kickmix has no H/T, so T merges do not become fewer Toffolis | jeffemerson, notebook *invsqrt X-ray: 46-qubit dead ends and ideas* (row 1) |
| N3 | Sampled "screening" before the official verifier | No gain on broken circuits: the verifier already aborts at the first failing chunk (~1.3 s for a broken invsqrt); dqi's ~1 s of load+hash is the floor | manwithacat, 753185 |
| N4 | GF(2) matrix products on tensor cores for linear X/CX/SWAP blocks | Only ~5% of dqi ops, none of invsqrt's or chemistry's, sit in long enough linear blocks | manwithacat, 753185 |
| N5 | Bit-sliced GPU kernel (Metal) | ≈ throughput of 10 M2 Pro CPU cores (~5e11 op·shots/s): no win on a laptop. Also "transform everything to GPU, searching very fast, but no improvement" | manwithacat 753185; yuanzhe zhou 757915 |
| N6 | GIL-released threads; Numba on the verifier's Python | Threads 2.7× dqi, 1.2× invsqrt; Numba 1.1× (per-shot Python is only 16–18% of time) | manwithacat, 753642 |
| N7 | Peak shaving on dqi (reorder, recompute, dirty or clean borrowing) | One qubit ≈ 1,100 Toffolis; peak hit in 2,046 windows; 0 known-value qubits to borrow; cheapest narrowing 23–42 Toffolis per window: nothing under the exchange rate | alex chilton, 753185 |
| N8 | Static ancilla remap on dqi | All 52 ancillas co-live from extract through unextract; best static Q' = 892 | dcvxzyk via rikirinjani, 753185 |
| N9 | Borrowing zero or affine-redundant wires, dropping degenerate Toffolis (invsqrt 63/66 q) | 0 degenerate controls in ~1,500 Toffolis, no borrowable wire | alex chilton, 753185 |
| N10 | Dirty-wire AND for a row-opening gate (invsqrt) | Not valid: borrowed wires read up to 8,474 times before the replay closes | alex chilton, 753185 |
| N11 | "Parking" a wire at invsqrt 46-qubit peaks (measure out, rebuild from wires that are affine/quadratic in others) | 74 uncovered peaks in the multiplier have no wire quadratic in the others | jeffemerson X-ray notebook, open question 1 |
| N12 | Pebbling by retention at 46 wires | Loses; exact per-adder pebbling DP over 31 adder groups not done (best case −0.37%) | Asterios Panoras, *QCC 28.12769: what works, what fails and why* |
| N13 | Recomputation caching in invsqrt | 698 of 1,550 Toffoli products recompute something released earlier, but a cached wire must save ~38 recomputations to pay one qubit (37.2 Toffolis) | jeffemerson X-ray notebook |
| N14 | Whole-pool search for multi-segment splices of 62-qubit invsqrt circuits | Only rediscovered the existing pick | RezOps, 759220 |
| N15 | Exact low-rank forms of small dqi blocks; GF(2^12) multiplier below 42 | No rank-11 form of a 12-Toffoli block (1,183 calls/shot); μ4(6) = 14 proven (Lean), so 42 is optimal for bilinear | RezOps 759220; topic 756476 |
| N16 | SAT for the 21-bit normalizing barrel shifter (min CSWAPs) | n ≤ 9 solved (~2n−5); w = 21, N in high 40s: 11 h CaDiCaL, neither SAT nor UNSAT | Rustam Bazarbayev, 757030 |
| N17 | ANF/Möbius synthesis or tree reordering of chemistry loader words | Each output bit's mixed ANF matrix has full GF(2) rank 15, like random; XOR of columns 0.95–1.0× as complex | alex chilton, 753185; Rustam 757030 |
| N18 | Deleting "never-firing" invsqrt gates validated on the default seed | Failed the grader with a 0.5-turn phase error: they were fixups for rare HMR outcomes | Lynn Langit, 755972 |
| N19 | AlphaTensor-Quantum on chemistry loaders | Untried; no released weights, training from scratch judged prohibitive | hengck23 / Zero, 754466 |
| N20 | Literature imports: QROM halving, sub-42 GF multipliers, reference DQI design (arXiv:2510.10967, 5.72e6 Toffolis), ECC repos, MBUC gadget catalogue | Nothing usable | alex chilton, 753185 |
| N21 | "propose → implement → verify → revise" single-lineage agent loop | Converges on small tweaks; drops a direction after one failed formulation; more agents = more variations of the same idea | uctpexx123, 755303 |
| N22 | GPT-6 "pro with codex" + GPUs | "no improvement"; another LLM copied a public circuit byte for byte when asked to re-implement a method | 757915 |
| N23 | Gu et al.'s own failure modes | Conservative model, repeated re-checking, preference for known methods; small-scale rules failing at full scale | Anthropic blog |

Lessons for our tooling: speed must come from removing per-shot Python (kickmix-engine) and from static or
bit-sliced evaluation, not from GPUs or screening; generic T-count optimizers do not touch the charged cost;
width wins need lifetime re-packing, not peak lowering (manwithacat 753403: at most 35 of 41 invsqrt scratch
qubits hold a 1 at any instant, yet all 41 are needed).

## 4. Ideas

**Conventions.** `Est. gain` is the expected **score** gain if the idea succeeds (score = geometric mean, so a
1% cost cut on one task = +0.33%; one invsqrt qubit = +0.77%; one dqi qubit = +0.075%; one dqi failed shot
avoided = +0.15%). Tooling ideas also give their speed-up, and an **equivalent score %** used for the priority:
the extra gain we expect from the iterations or the risk it saves over the remaining ~6 working weeks (a
judgement, stated so it can be challenged). `Priority = P × gain(%) / effort(h)`. Efforts assume one person
directing Claude Code. Rows are sorted by priority within each category; section 4.3 merges both.

### 4.1 A. Improve our code / tooling

| ID | Idea | Mechanism | Tasks | Est. gain | P | Effort (h) | Priority | Risks | Sources |
|---|---|---|---|---|---|---|---|---|---|
| GEN-A01 | kickmix-engine as our fast evaluator | Install the Linux wheel in `~/venvs/gqai`; run `--cross-check` against the official verifier on our v12 baseline for all 3 tasks; then use it for all iteration runs and multi-salt sweeps; the official verifier stays the gate | all | 7–9× faster verification (author-reported), dqi memory 9 GB → < 1 GB so 4+ workers fit in WSL; ≈ 0.5% eq. | 0.85 | 3 | **0.142** | Third-party binary (Apache-2.0); version drift vs kickmix 1.0.0; never use it for the final gate | 753642; github.com/manwithacat/kickmix-engine |
| GEN-A05 | Nulls ledger seeded from the community | `experiments/LEDGER.md` (PLAN S2.4) plus manwithacat's 258-row negative-results dataset, topic 753185, the X-ray and "what fails" notebooks; every agent brief includes the relevant rows; each entry records the *family* and *why it closed* (a bound, or only a failed implementation) | all | Avoids repeating closed work; ≈ 0.3% eq. | 0.95 | 3 | **0.095** | Over-trusting a "closed" verdict outside its family | 753185; manwithacat dataset; jeffemerson notebook |
| GEN-A02 | Static IR + exact static cost | Own numpy struct-of-arrays reader/writer for `.kmb` (model: GrantPrater's `vqcc_ir.py`, Apache-2.0, credited); exact expected Toffolis by propagating condition probabilities (HMR bits are fair coins; data-dependent conditions through the bit-sliced simulator); qubit count, liveness, stage tags from INFO_TAG; < 1 s per circuit | all | Evaluator tier T0 in seconds instead of minutes; ≈ 0.4% eq. | 0.9 | 6 | **0.060** | Data-dependent conditions need simulation, not static rules; must agree with the verifier's sampled mean | GrantPrater dataset README |
| GEN-A06 | Daily frontier-ingestion agent | A scheduled Claude Code routine: `kaggle kernels list --sort-by scoreDescending`, pull new notebooks (read only), extract SHA-256s and claims, re-verify new files with GEN-A01, stream-diff them against our builds, open "port to generator" tickets with credit details | all | Keeps us on the moving plateau without manual reading; ≈ 0.6% eq. | 0.8 | 8 | **0.060** | Untrusted downloads (never execute); time sink if it ports every 0.01% | RezOps (759220); CLAUDE.md research rule |
| GEN-A04 | Salt-robustness gate for every HMR-dependent edit | Any edit touching measurement-conditioned code must pass 4 salts (screen) and 16–48 salts (finalist) with kickmix-engine on Kaggle CPU; dqi ranked by expected cost; invsqrt and chemistry also by exact proof | all | Avoids −0.15% per dqi failed shot and the 0-score case; ≈ 0.3% eq. (risk-weighted) | 0.9 | 5 | **0.054** | Salts are evidence, not proof; Kaggle quota | 755972; alex chilton salt notebook; PLAN S2.3 |
| GEN-A07 | I/O and caching hygiene | Keep work trees and big `.kmb` files on the WSL ext4 filesystem (`~/work`), not `/mnt/c` (9P is slow for the 846 MB dqi file); hash-keyed cache of (circuit SHA, salt, verifier version) → result JSON | all | 1.5–3× on load-heavy dqi runs, no repeated runs; ≈ 0.1% eq. | 0.9 | 2 | **0.045** | Two copies drift; the repo stays the source of truth | CLAUDE.md environment notes |
| GEN-A09 | Kaggle batch-verification orchestrator | Script `kaggle datasets version` (private dataset of candidates), `kaggle kernels push` of a private verification notebook (4 cores, 30 GB, 9 h) and `kernels output` back into the ledger; several sessions in parallel | all | Frees the laptop; 48-salt dqi sweeps and full invsqrt runs in parallel; ≈ 0.15% eq. | 0.8 | 6 | **0.020** | Quotas and queue delays; keep datasets private (rule C.5.d) | PLAN 4.2 |
| GEN-A03 | Task-agnostic finishing pipeline | Port Syed's passes (`audit_affine`, `audit_pairs` + `hold_pass`, `postpass`, `merge_zpow`) and GrantPrater's (`measure_out`, `free`, `pair_window`) onto GEN-A02's IR with a domain sampler per task (invsqrt: all x; chemistry: all 1,220 addresses × HMR draws; dqi: verifier-distributed syndromes); run to fixpoint after every generator change | all | 0.01–0.1% per task per round; ≈ 0.25% | 0.6 | 12 | **0.0125** | On dqi the sampled domain only yields candidates; every rewrite needs GEN-A10 + salts | Syed `tools/pipeline.sh`; GrantPrater README |
| GEN-A10 | Local proof per rewrite | Prove each rewritten window before any full run: GrantPrater's forms-with-BDD-coefficients prover generalised to any window, or SymQ for X/CX/CCX/HMR+CZ windows; reject on timeout | all | Catches rare-branch phase bugs early, fewer wasted full runs; ≈ 0.2% eq. | 0.6 | 10 | **0.012** | Unsupported ops; a timeout is not a proof | 754815; GrantPrater `vqcc_pair_window.py` |
| GEN-A12 | Design-space optimizer over generator parameters | Expose every numeric knob (table rows, word widths, per-stage qubit caps, dqi cutoffs, chemistry bit widths) with a static cost model; enumerate or anneal using the per-task exchange rates; output the (Q, T) Pareto front | all | ≈ 0.3% from better trade-offs | 0.4 | 10 | **0.012** | Knobs interact through the error budget; exact validation needed | Gu et al. (tunable trade-offs); PLAN exchange rates |
| GEN-A08 | Peak map + functional-dependency mining | At every width peak, record full-domain wire traces (bit-sliced) and test, for each live wire, whether it is **any** Boolean function of the other live wires (hash of state-minus-wire → value), then its ANF degree and AND cost; rank "measure out here, rebuild there" candidates by Toffoli cost against the qubit's value | invsqrt, dqi, chem loaders | ≈ 0.5% (invsqrt qubit 0.77%, dqi 0.075%) | 0.2 | 10 | **0.010** | No affine/quadratic parking existed at 46 q (N11); higher-degree rebuilds may cost more than a qubit is worth | 753403 (hengck23 idea); jeffemerson X-ray; Syed measure-then-rebuild |
| GEN-A11 | Gadget library + matcher | Mine (before, after) gadget pairs from public lineage diffs (`stream_diff`), canonicalise up to wire relabelling and condition bits, search our circuits of all tasks for matching sites, re-implement hits in the generator | all | 0.05–0.2% | 0.45 | 10 | 0.0068 | Matches hold only under context invariants (prove each) | GrantPrater's transplant of 8 closings |
| GEN-A13 | Scaled-down proxy instances | Parametric small versions (invsqrt with 12–14-bit x; RS over GF(2^6); chemistry with a toy table) so a new design can be checked exhaustively in seconds before a full-size build | all | Faster triage of new families; ≈ 0.2% eq. | 0.4 | 12 | 0.0067 | Small-size rules can fail at full size (Gu et al.); always confirm at full size | Anthropic blog |
| GEN-A15 | Branch-cost balancing | Conditioned Toffolis count only on shots where they run: put expensive work in the less frequent branch of data-dependent conditions; measure-then-condition rather than control when a branch is rare | invsqrt, dqi | ≈ 0.08% | 0.3 | 6 | 0.0040 | Small; fair HMR coins leave nothing to balance | Verifier counts executed CCX/CCZ/CSWAP per shot |
| GEN-A14 | CUDA bit-sliced simulator (RTX 4060) | CuPy kernels for full-domain invsqrt scans and dqi screening | invsqrt, dqi | ≤ 2–3× over 12 cores; ≈ 0.05% eq. | 0.3 | 12 | 0.0013 | N5: the Metal version only matched the CPU | 753185, 757915 |

### 4.2 B. New technique or new pipeline

| ID | Idea | Mechanism | Tasks | Est. gain | P | Effort (h) | Priority | Risks | Sources |
|---|---|---|---|---|---|---|---|---|---|
| GEN-B11 | Divergent-family brainstorm + red team, weekly | One high-effort planner session per task (Opus/Fable-class) gets the spec, our generator's stage map, the nulls ledger and the cost exchange rates, and must propose **≥ 5 qualitatively different datapaths** for the dominant stage (not tweaks), each with a back-of-envelope (Q, T) estimate; a second session red-teams them against the nulls; the 1–2 survivors become new families in GEN-B01 | all | ≈ 1% (one new family that beats the incumbent on one stage) | 0.25 | 8 | **0.031** | Estimates on paper are optimistic; families need real builds | Gu et al. ideation step; Anthropic blog mitigation; 755303 |
| GEN-B01 | Evolutionary agent loop over our generators | OpenEvolve-style loop run by Claude Code (design in section 5): genomes are generator programs + patch modules + pass configs; MAP-Elites on (Q, Toffolis) per task; 3–5 families per task as islands; fresh subagent per iteration with parent, evaluator feedback, archive and ledger; tiered evaluator T0–T3 | all | 2–5% over the remaining weeks (Gu et al. found 4–39× over **literature** baselines; ours are crowd-optimized, so far less) ≈ 3% | 0.35 | 40 | **0.026** | Token cost; convergence on tweaks (N21); reward hacking of the evaluator (pin hashes, read-only verifier) | arXiv:2610.22367; OpenEvolve; arXiv:2611.09822 |
| GEN-B02 | Principle sweep, stage by stage | For every stage of every generator, one targeted agent per principle: (a) conditional routing (CSWAP into a shared datapath); (b) recover history/branch bits from invariants instead of storing; (c) rescaling freedom; (d) failure budget sized from measured rates; (e) bigger table instead of a multiplication, erased by MBU; (f) dyadic intervals + cyclic shift to one shared local coordinate; (g) accumulate products directly into the destination; (h) spooky pebbling of long-lived values. Each produces a priced candidate or a recorded null | all | ≈ 1.5% (a few stage-level wins) | 0.3 | 20 | **0.023** | Many are already in the public generators; need a stage map first (PLAN S2.6, S2.13) | arXiv:2610.22367 §III; 759220 |
| GEN-B04 | Qubit-cap ladder harness | Generalise Syed's method: rebuild with an allocator that aborts at Q−1, log every blocking allocation with its live set and stage, cluster them, and hand each cluster to an agent with the measure-and-rebuild, borrow and reschedule options priced by the exchange rate | invsqrt (main), dqi, chem loaders (≤ 128 cap) | invsqrt −1 qubit = 0.77% | 0.2 | 10 | **0.015** | 44 → 43 is hard (X-ray: multiplier registers alone reach the cap at 46 in that family) | Syed README (cap 46, 13 patches); jeffemerson X-ray |
| GEN-B12 | One end-to-end error budget, fitted to the exact quantized computation | Fit all approximation parameters (table coefficients, truncations, rotation precision; chemistry rounding) jointly against the bit-exact integer model with an LP/MILP on the worst error/tolerance, instead of stage by stage | invsqrt, chemistry | ≈ 0.4% (one table bit or a few rotations) | 0.25 | 12 | 0.0083 | Margins are already saturated (worst 0.99956 of tolerance; chemistry at 19.99999 mHa); overlaps task agents | Gu et al. §III; 753185 (rotation retune null) |
| GEN-B03 | Width-optimal rescheduling by spooky pebbling DP/SAT | Extract a stage's dependency DAG (e.g. invsqrt's 31 adder groups), encode compute / measure-out / recompute moves as a spooky pebble game with a Toffoli budget, solve by DP or SAT (Quist–Laarman encodings) for minimum peak width | invsqrt, dqi | ≈ 0.77% (one invsqrt qubit; open best case −0.37% cost) | 0.2 | 20 | 0.0077 | PSPACE-complete in general; DAG extraction from a generator is real work | arXiv:2401.10579; arXiv:1904.02121; Asterios notebook (open DP) |
| GEN-B10 | Conditional-routing detector | Static search for duplicated arithmetic under complementary or mutually exclusive controls (same op-type sequence on different registers inside opposite condition blocks or controlled by exclusive flags); propose a CSWAP-in / shared kernel / CSWAP-out rewrite | invsqrt (normalizer, divider branches), dqi | ≈ 0.5% where it hits | 0.15 | 12 | 0.0062 | Public generators may already route; CSWAP cost can exceed the saving for short kernels | arXiv:2610.22367 (halved multiply-add) |
| GEN-B08 | Conditionally clean ancillae everywhere | Systematically apply Khattar–Gidney constructions (laddered toggle detection, conditionally clean ancillae for n-controlled NOT, comparators, incrementers, unary iteration) to stages that hold a clean scratch only to support a comparison or a multi-controlled gate | invsqrt, dqi, chem loaders | one invsqrt qubit = 0.77% | 0.12 | 16 | 0.0058 | N7, N9, N10: borrowing was measured null on several designs | arXiv:2407.17966 |
| GEN-B13 | Shared library of optimal primitives | One verified, parametric library (MBU adders, Cuccaro/Gidney adders, comparators, incrementers, unary iteration, select-swap and variable-spaced QROM, CSWAP shifters) with exact Toffoli/width formulas; generators pick implementations by configuration, so the GEN-B01 loop can swap primitives | all | ≈ 0.3% | 0.3 | 16 | 0.0056 | Re-implementing what the generators already inline; integration effort | Gidney 1709.06648; Babbush 1805.03662; Low 1812.00954 |
| GEN-B06 | XAG multiplicative-complexity resynthesis with reachable-state don't-cares | Extract classical subnetworks (lookup tables, decoders, flag logic), build an XAG, get don't-cares from the full-domain simulation (unreachable states), run cut rewriting with exact-MC databases (≤ 6 inputs) and don't-care AND→XNOR passes (mockturtle), re-emit with MBU and a width check | invsqrt tables, chem loaders, dqi control logic | ≈ 0.3% | 0.3 | 20 | 0.0045 | N17: loader outputs look random (full rank); width blow-up; public tables already re-synthesised (27 + 29 ANDs) | mockturtle docs; Testa et al. DAC 2019; Çalık et al. 2019 |
| GEN-B07 | SAT superoptimizer for small windows | Enumerate windows ≤ 8 wires, ≤ 10 Toffolis; encode the kickmix gate set including HMR + conditioned CZ fixups; minimise Toffolis under the window's reachable-state don't-cares with CaDiCaL/Kissat (pysat); prove and splice | all | ≈ 0.2% | 0.3 | 16 | 0.0037 | N16: big windows intractable; most small windows already optimal | 757030; Quartz/QUESO |
| GEN-B09 | Evolve the search heuristics, not the circuits | FunSearch/VarTODD-style: keep our passes and searches fixed and let Claude evolve their policy functions (hold-pass wire assignment, loader synthesis search order, refit schedules, cutoff selection) scored by the circuit cost they produce | all | ≈ 0.3% | 0.2 | 16 | 0.0037 | N1 (ShinkaEvolve + Qwen null on a gauge optimizer); strong models needed | arXiv:2603.29894; FunSearch; ShinkaEvolve |
| GEN-B05 | Minimum-CCZ phase polynomials | In windows where the circuit is diagonal (HMR repayment blocks, dqi phase scan, invsqrt phase stages), express the phase as a cubic polynomial over GF(2) of affine wire functions and resynthesise with fewer CCZ (CZ and Z are free): Reed–Muller decoding, TODD-like rank reduction, or an AlphaTensor-style symmetric-tensor search | dqi, invsqrt | ≈ 0.15% | 0.25 | 12 | 0.0031 | N2: windows rarely have reducible cubic terms; most charged gates are CCX, not CCZ | arXiv:2402.14396; arXiv:1712.01557; Atakan "Terminal Readout and Phase Folds" |
| GEN-B15 | ZX-calculus (PyZX) Toffoli reduction | `full_reduce` + extraction on windows | all | ≈ 0.1% | 0.03 | 6 | 0.0005 | N2: 0 Toffolis removed | jeffemerson X-ray |
| GEN-B14 | RL / MCTS over rewrite sequences | Quarl/AlphaZero-style policy over local rewrites | all | ≈ 0.3% | 0.05 | 40 | 0.0004 | Expensive; nothing in the forum supports it | arXiv:2307.10120; 754557 |
| GEN-B17 | Lean-verified rewrites (DeepMind style) | Formalise the kickmix semantics and prove gadget rewrites in Lean | all | ≈ 0.05% | 0.1 | 40 | 0.0001 | Huge effort; the official verifier + exact proofs already suffice | arXiv:2611.09822 |
| GEN-B16 | Train AlphaTensor-Quantum from scratch | Symmetric-tensor RL for T/CCZ count | chem loaders, dqi | ≈ 0.2% | 0.03 | 60 | 0.0001 | No weights; compute; N19 | 754466; arXiv:2402.14396 |

### 4.3 Combined ranking (top 12)

| Rank | ID | Priority | One line |
|---:|---|---:|---|
| 1 | GEN-A01 | 0.142 | kickmix-engine (cross-checked) as the iteration evaluator |
| 2 | GEN-A05 | 0.095 | Nulls ledger seeded from the community, read by every agent |
| 3 | GEN-A02 | 0.060 | Own static IR with exact expected Toffolis and liveness in < 1 s |
| 4 | GEN-A06 | 0.060 | Daily frontier-ingestion agent (re-verify, diff, port tickets) |
| 5 | GEN-A04 | 0.054 | Multi-salt gate for every measurement-dependent edit |
| 6 | GEN-A07 | 0.045 | WSL ext4 work dirs + result cache keyed by SHA/salt |
| 7 | GEN-B11 | 0.031 | Weekly divergent-family brainstorm + red team |
| 8 | GEN-B01 | 0.026 | Evolutionary agent loop over our generators (section 5) |
| 9 | GEN-B02 | 0.023 | Gu et al. principle sweep, stage by stage |
| 10 | GEN-A09 | 0.020 | Kaggle batch-verification orchestrator |
| 11 | GEN-B04 | 0.015 | Qubit-cap ladder harness |
| 12 | GEN-A03 | 0.0125 | Task-agnostic finishing pipeline (Syed + GrantPrater passes) |

Reading the ranking: the cheap tooling items come first because they multiply everything else; the large
expected gains are in B01/B02/B11, which only pay once A01–A06 exist. B01 is the vehicle in which B02, B04,
B10 and B12 run as mutation operators.

## 5. Proposed agent-loop design for this team

The design copies what worked for Gu et al. (generator programs, separate families, MAP-Elites, fresh
session per iteration, insight memory) and for DeepMind (independent teams, a shared wiki of best solution
and dead ends, periodic critique), and adapts it to one laptop (12 cores, 8 GB WSL), Kaggle CPU sessions and
Claude Code subagents. It respects CLAUDE.md: no git writes by agents, English in files, exactness first.

### 5.1 Genome: what evolves

- **A genome is a generator variant, never a circuit file.** For invsqrt this is Syed's layout: the frozen,
  hash-pinned base generator plus an ordered list of patch modules (`patches/PATCH_ORDER.txt`) plus the pass
  configuration (`tools/pipeline.sh` options, qubit cap). For dqi and chemistry the same idea: base builder +
  patch modules + parameters (cutoffs, table shapes, word widths). This keeps every candidate rebuildable from
  source (prize rule) and makes diffs small and reviewable.
- Each genome lives in `experiments/<task>/<family>/<id>/` with `genome.json` (parent id, patches, parameters,
  pass config), the patch files, the build log, `result.json` (tiers reached, costs, SHA-256) and a 5-line
  `notes.md` written by the agent (what was tried, why, what was learned).

### 5.2 Evaluator: four tiers (cheap first, abort early)

| Tier | What | Time (target) | Where | Promotes when |
|---|---|---|---|---|
| T0 | Build + static IR (GEN-A02): qubits, exact expected Toffolis, liveness; metric checks of topic 751951 | seconds to 3 min (build dominates) | laptop | (Q, T) lands in an empty MAP-Elites cell or beats the cell's incumbent by the noise margin |
| T1 | Bit-sliced screen: invsqrt all 3.1M x with random HMR (no phase) or kickmix-engine reduced run; chemistry all addresses; dqi 4,096 shots | 10 s – 2 min | laptop | No value error; phase spot checks pass |
| T2 | Exact or full check: invsqrt `prove_any` (all x, all HMR transcripts); chemistry official verifier (~1 min); dqi kickmix-engine 102,400 shots at 4 salts | 2–15 min | laptop (2 slots) | PASS on all |
| T3 | Gate: unmodified official verifiers in WSL (`~/venvs/gqai`), dqi 16–48 salts on Kaggle (GEN-A09), plus Atakan Aldemir's adversarial low-weight dqi inputs | 1–9 h | Kaggle CPU + laptop overnight | Candidate for submission |

Fitness = the task's official cost (invsqrt Q × (Toffolis + T/4); dqi R(ε) × Q² × Toffolis with ε from the
multi-salt failure rate; chemistry λ_eff × T_total × Q_total). Invalid = ∞. Ties within noise go to the
smaller genome diff. **Anti-reward-hacking**: verifiers, kickmix and oracles are read-only and SHA-pinned (as
Syed's sources already check); a candidate that edits anything under `google-verifiable-quantum-circuit-compilation/`
or the pinned base is rejected; T3 always uses the untouched official verifier.

### 5.3 Population and selection

- **Families as islands**: per task, 3–5 families (for example invsqrt: "Syed 44-q line", "table-heavy
  dyadic-interval line", "new divider-free line"; dqi: "Syed 890 line", "BM with history recovered from state"),
  seeded by GEN-B11. Families do not compete with each other for the first 2–3 rounds, so a weak new design can
  mature (Gu et al.).
- **MAP-Elites archive per family** on (qubits, expected Toffolis), bins of 1 qubit × 1% Toffolis; each cell keeps
  its best cost. The archive exposes the Pareto front, which matters because each task has a fixed exchange
  rate between qubits and Toffolis.
- **Parent selection**: 70% from the top-5 cells by cost, 30% from under-explored cells (fewest children),
  ShinkaEvolve-style; reject a child whose patch is a near-duplicate of an existing one (normalised diff
  similarity > 0.9) before building.
- **Migration**: every 2 rounds the best genome of each family is offered to the others as "archive
  alternatives" in the brief (not forced).
- **Retirement**: a family with no archive improvement in 3 rounds is frozen; its notes go to the wiki as a
  dated null with the reason (bound vs failed implementation).

### 5.4 Mutation operators (each is a subagent brief template)

1. **Targeted stage rewrite** (LLM): rewrite one stage's patch to lower its Toffolis or width, given the stage
   map and the stage's share of the cost.
2. **Principle application** (GEN-B02): apply one named principle to one stage.
3. **Parameter mutation** (no LLM): perturb numeric knobs (GEN-A12) inside their valid ranges.
4. **Cap tightening** (GEN-B04): rebuild at Q−1, fix the first blocking cluster.
5. **Pass-pipeline mutation** (no LLM): enable, reorder or re-parametrise finishing passes (GEN-A03).
6. **Port from frontier** (GEN-A06): re-implement a public improvement in our generator, with credit.
7. **Crossover**: combine two families' patches that touch disjoint stages (checked by stream-diff of their
   outputs); re-run all tiers.
8. **New family** (GEN-B11): from the spec only, a qualitatively different datapath for the dominant stage.

### 5.5 Running it with Claude Code subagents

- **Orchestrator** = the main Claude Code session (Opus 5.5). It owns `experiments/`, selects parents, writes
  briefs, launches subagents, reads `result.json`, updates the archive and the wiki. It never edits generators
  itself, and never runs git commit/push (the user handles git).
- **Workers** = one fresh subagent per iteration (`isolation: "worktree"` so parallel workers do not collide),
  given: the task contract (from `verifiers/`), the parent genome, its evaluator output and stage costs, 3–5
  archive alternatives with costs, the relevant ledger nulls, the wiki page of the family, and a hard time box
  (45 min). It must return a complete patch set and a 5-line note, and it may run T0–T1 itself to iterate.
  Sonnet-class models for operators 1–6, Opus for 7–8, a high-effort planner only for GEN-B11 (token cost,
  topic 757915).
- **Concurrency**: 3 workers at a time on the laptop (each build/T1 uses 2–3 cores; WSL 8 GB caps dqi at 1–2
  concurrent T2 runs with kickmix-engine, about 0.4–0.6 GB each). One T2 slot reserved. T3 goes to Kaggle.
- **Budget**: a round = 6–10 iterations per task family, about one evening; target 2 rounds per day during
  Phase 2 (21–30 Nov), i.e. 150–250 evaluated genomes per week across tasks. Stop rule per family: 3 rounds
  without improvement → freeze.
- **Kaggle**: a private dataset `candidates-<date>` holds the T2 survivors; a private notebook runs the
  official verifiers (invsqrt full run, dqi 16–48 salts, chemistry) in parallel sessions (4 vCPU, 30 GB, 9 h)
  and writes a JSON receipt per circuit; the orchestrator pulls receipts into the ledger. Kaggle P100 only for
  chemistry float64 refits.
- **Wiki** (`experiments/WIKI/<task>.md`): best verified genome per family, reusable lemmas (invariants used by
  rewrites), and reviewed notes on ideas and dead ends. A short weekly "conference": the orchestrator asks one
  reviewer subagent to critique each family's last round and to propose cross-family transfers.
- **Anti-convergence rules** (topic 755303): every round must include at least one operator 2 or 8 child per
  family; a closed idea needs either a bound or two independent failed formulations; and briefs state the
  family's exchange rate so agents stop chasing 0.01% tweaks when a qubit is on the table.
- **Human checkpoints**: the user reviews the archive and the wiki at each step review (PLAN section 5); every
  submission candidate passes T3 and the packaging checks (PLAN S1.6) before the user submits it.

## 6. Leaderboard reading (2026-11-17 21:13 UTC, 3,274 teams)

| Threshold | Teams | Notes |
|---|---:|---|
| ≥ 31 | 2 | aadi_stonks 31.44520 (17 Nov 17:17), Z³C¹oherence 31.42723 (17 Nov 12:45) |
| ≥ 30 | 5 | Entangled Team 30.285 (15 Nov), Rithvik Redrouthu 30.131 (15 Nov), Yicheng Guang & Taiseiu 30.009 (16 Nov) |
| ≥ 29 | 9 | + Schrödinger's cat lived 29.979, Slingshot Effect 29.573, Christoffer Thimsen 29.349, Some Joint 29.144 |
| ≥ 28.7 | 16 | gold line ≈ rank 16 (28.713) |
| ≥ 28.2 | 61 | strictly above the public plateau |
| ≥ 28.18 | 267 | 139 teams at exactly 28.1807 (the plateau) |
| ≥ 28.15 | 459 | |
| ≥ 28.0 | 947 | |
| = 1.0 | 121 | reference submissions |

What the numbers imply (score is a geometric mean, so a score ratio r needs a cost ratio r³ on one task):

- **Leaders (31.44)**: +11.6% over the plateau = a **28% cost cut on one task**, or 10.4% on each task.
  Concretely, on one task alone: invsqrt at 53,670 (e.g. 44 q at ~1,220 Toffoli-eq instead of ~1,695, or ~32
  qubits at today's Toffolis); dqi at 2.49e11 (264k Toffolis at 890 q, or ~297k at the 840-qubit syndrome floor);
  chemistry at 2.12e8 (implausible given the λ_eff ≥ ~46.85 floor, the 20 mHa gates and a model-fixed Q_total).
  **Most likely: a new invsqrt architecture (−8 to −12 qubits) plus a dqi decoder with ~15–25% fewer Toffolis.**
  An invsqrt move from 44 to 36 qubits alone is +6.9%.
- **Two teams within 0.06% of each other, both new on 17 Nov**, suggests a shared source (a team merger, a
  private exchange, or the same newly published idea: DeepMind's arXiv:2611.09822 v3 appeared on 16 Nov, and
  its arcsin result uses kickmix-style numerical-function tricks that transfer to invsqrt). This is a
  hypothesis, not evidence.
- **Teams 3–5 (30.0–30.3, last submitted 15–16 Nov)** were 7% above today's plateau = a 19% cut on one task;
  they have not moved since the public plateau rose, so part of their edge may now be public, or they hold it
  in reserve.
- **Gold line (28.71)** = +1.9% over the plateau = 5.5% cost cut on one task (e.g. invsqrt 44 → ~41.5 qubits
  equivalent, or dqi −5.5% Toffolis). This is our realistic Phase 2 target.
- **The plateau is crowded**: 139 teams at exactly 28.1807, 267 within 0.0007. Ties break on earliest
  submission, so matching the plateau late is worth nothing; dqi failed-shot draws (−0.15% each) will reshuffle
  the plateau on the private board.
- **Where the private gains likely are** (forum + notebooks): (1) invsqrt is the only task where a generator
  redesign moved the public frontier by >10% in a week (62 → 46 → 44 qubits) and has no hard width floor;
  (2) dqi Toffolis: the public line is still the Syed v12/v13 decoder plus local patches; Gu et al.'s
  history-from-state and failure-budget ideas, and Rustam's "cheap inverse of in-place GF(2^m) chains"
  question, point to the reverse/extract stages (≈ 25% of Toffolis in unextract alone, rikirinjani's profile
  on an older decoder); (3) chemistry is boxed in by the cost model; gains there are 0.1–0.7% per step.
- **Implication for us**: build the fast loop first (GEN-A01, A02, A05, A06), reach the plateau from source,
  then spend most agent budget on new invsqrt and dqi families (GEN-B11 → B01), not on public-style
  single-gate rewrites.
