# LLM-driven search loop, local models and fine-tuning: costing and recommendation (2026-11-21)

Author: ML/infra analyst (research agent), 2026-11-21. Scope: milestone M3.4 ("Evolutionary agent loop") of
`Analysis/PLAN.md`, ideas GEN-B01 and INV-B01 of `Analysis/IDEAS.md`, and the agent-loop design in
`Analysis/research/general.md` §5. No file other than this one was changed; nothing was installed. Facts about
models, prices and tooling were checked on the web on 2026-11-21 (sources in §7). Numbers marked "est." are
this report's estimates, not measurements. `Analysis/research/public_frontier_2026-11-21.md` (the ThemBones
generators) did not exist yet when this was written; §1.5 says how it changes the picture.

## 0. Summary

1. **The loop is cheap to run; the bottleneck is the hit rate, not the tokens.** One "diff-style" proposal
   (task contract + one stage's source + parent and 3 archive diffs + feedback ≈ 15k tokens in, ≈ 5k out)
   costs about **$0.004 on Claude Haiku 5.5**, $0.08 on Sonnet 5.5 and $0.16 on Opus 5.5 (API list prices,
   prompt caching on). 1,000 proposals cost $4 / $80 / $160. **API usage is billed separately from the Claude
   subscription quota**, so the quota problem disappears with an API key and a hard spend cap. Our evaluators
   allow 60–200 screened candidates per hour on the laptop (T0 `tools/kmb.py` < 1 s, T1 kickmix-engine
   15–18 s, builds seconds to 3 min); the proposer is the bottleneck only with a local model.
2. **Published evolve-style wins need 150 (ShinkaEvolve, frontier ensemble) to thousands (AlphaEvolve) to
   millions (FunSearch) of proposals,** on objectives where the baseline was not crowd-optimised. Our circuits
   are locally tight (0 removable Toffolis in several scans). Est. 0.1–0.5% of proposals improve an incumbent
   at the stage level; a January campaign of 2,000–4,000 proposals would give **≈ +0.3 to +1% score, with
   P ≈ 0.25** (lower than GEN-B01's register estimate of 3% at P 0.35, which assumed a frontier proposer
   discovering new families). Good targets: small sub-programs with an exhaustive millisecond checker and a
   scalar cost (loaders/QROM, GF(2^12) gadgets, chemistry loader-compiler heuristics, pass orderings). Bad
   targets: whole-generator redesigns (2.2 MB of Python per generator, 13 interlocked invsqrt patches),
   dqi failure-rate trade-offs (noisy multi-salt cost), anything pure-numeric (use SAT/LP/gradient instead).
3. **Best local model for 8 GB VRAM + 16 GB RAM: Qwen3.6-35B-A3B** (Apache-2.0, SWE-bench Verified 73.4,
   LiveCodeBench v6 80.4) at UD-IQ3_XXS (13.2 GB) with llama.cpp `--n-cpu-moe`, est. 15–30 tok/s.
   Runners-up: gpt-oss-20b (MXFP4 ≈ 13 GB, SWE-bench Verified 60.7) and Qwen3.5-9B (Q4_K_M ≈ 6 GB, fully on
   GPU, 33 tok/s measured on an RTX 4060, LiveCodeBench v6 65.6). The best open models (Qwen3.8-27B,
   Qwen3.8-Flash-Next, GLM-5.3-Flash, DeepSeek-V4-Flash) do not fit or run at ≈ 3.5 tok/s. **A local model is
   not cheaper than Haiku 5.5** (laptop electricity ≈ €0.03/h ≈ €0.002–0.003 per local proposal), is 5–10×
   slower, weaker on coding benchmarks (Haiku 5.5: SWE-bench Pro 64.8), and competes with the evaluator for
   CPU and RAM. Use it only as a $0 overnight fallback.
4. **Fine-tuning does not pay off before 2027-01-20.** bf16 LoRA needs ≈ 22 GB for a 9B model and ≈ 74 GB for
   35B-A3B (Unsloth), QLoRA is discouraged for Qwen3.5/3.6, the laptop has 8 GB and Kaggle T4/P100 have no
   bf16; GRPO needs ≥ 10k verifier-scored rollouts with an almost always zero reward; we would have tens of
   accepted edits, not the ≥ 1,000 that rejection-sampling fine-tuning needs. A retrieval bank of verified
   edits, a kickmix API cheat-sheet and the nulls in the prompt give most of the benefit at no training cost.
5. **Recommendation:** a minimal loop (OpenEvolve or a ≈ 400-line custom loop; the two are about equivalent
   per Gideoni et al. 2026) with **Haiku 5.5 for 80% of proposals and Sonnet 5.5 for 20%** (ideation, hard
   operators), our tiered evaluator, a MAP-Elites archive on (width, Toffolis), and a **$25 cap** for a
   **1-day go/no-go: "loader golf"** on the invsqrt K = 6 loaders of S2.8 (exact 64-address check in
   milliseconds, our constructive baseline 44 ANDs per loader, target ≈ 27–29), raced against a non-LLM
   baseline of the same wall time. Place it as M3.4a at the start of Phase 3 (or on a free day before
   1 Dec); production (M3.4b, ≤ $100) only on a go.

## 1. Pricing an LLM-driven search loop on paper

### 1.1 What one proposal costs in tokens

Our generators are far too big to fit in a prompt: measured today, `generators/invsqrt` has 258 Python files,
2.18 MB, 43,648 lines; `generators/dqi` 259 files, 2.24 MB, 46,249 lines; `generators/dfthcsos` 59 files,
0.29 MB, 6,604 lines (≈ 550k tokens per generator for the first two). Every proposal must therefore target
**one marked block** (AlphaEvolve's `EVOLVE-BLOCK` pattern): one patch module, one emitter function, one
heuristic.

| Prompt part | Tokens (est.) | Cacheable across proposals |
|---|---:|---|
| System prompt + task contract (register rules, cost formula, kickmix builder cheat-sheet, MBU idioms) | 2,500–4,000 | yes |
| Stage map with each stage's share of the cost | 500–1,000 | yes |
| Target block source (5–30 KB) | 1,500–8,000 | yes per parent |
| Parent result + evaluator feedback (failing input, Toffolis per stage) | 500–1,500 | no |
| 3–5 archive alternatives as diffs with their costs | 2,000–4,000 | partly |
| Relevant nulls from `experiments/NULLS.md` + family notes | 1,000–2,000 | yes |
| **Input total** | **≈ 10k–20k (use 15k)** | ≈ 70–80% cache hits |
| Output: diff or rewritten block (0.5–3k) + reasoning (2–8k) | **≈ 3k–10k (use 5k)** | — |

Agentic workers (a Claude Code subagent that edits, builds and iterates for 45 min, as `general.md` §5.5
proposed) are 20–30× more expensive: est. 200k–500k input tokens (mostly cached) and 15k–30k output per
iteration. This mode is what burns the subscription quota.

Cost per proposal at list prices (Anthropic pricing page, 2026-11-21; Claude 4.7+ models use a tokenizer that
gives ≈ 30% more tokens for the same text, included below):

| Proposer | $/MTok in / cache hit / out | One diff proposal (15k in, 75% cached, 5k out, +30% tokens) | 1,000 proposals | One agentic iteration (350k in mostly cached, 25k out) |
|---|---|---:|---:|---:|
| Claude Haiku 5.5 (prompt ≤ 100k) | 0.10 / 0.01 / 0.50 | $0.004 | **$4** ($2 batch) | $0.02–0.05 |
| Claude Sonnet 5.5 | 2 / 0.10 / 10 | $0.08 | **$80** ($40 batch) | $0.35–0.60 |
| Claude Opus 5.5 | 4 / 0.20 / 20 | $0.16 | $160 | $0.70–1.20 |
| DeepSeek V4.1 Flash (first-party, peak) | 0.30 / ≈ 0.03 / 1.20 | $0.009 | $9 | $0.05–0.10 |
| Local Qwen3.6-35B-A3B on the laptop | electricity ≈ €0.03/h | ≈ €0.002–0.003 | ≈ €2.5 (and 60–100 h) | not practical (est. 1–2 h each) |

The money is small; the binding constraints are **wall time** and **hit rate**. A $25 cap buys ≈ 5,000 Haiku
proposals or ≈ 300 Sonnet proposals.

### 1.2 Evaluations per hour (our evaluators, measured where stated)

| Task | T0 static (`tools/kmb.py`) | T1 screen | T2 exact / multi-salt | Screens per hour, 3 slots on 12 cores (est.) |
|---|---|---|---|---:|
| invsqrt | build (seconds to 3 min, generator-dependent) + 0.24 s | kickmix-engine full verify **15 s** (S2.1) | `prove_any` all x × all HMR transcripts, 2–15 min | 60–300 (build-bound) |
| dqi | build + 0.85 s (static only; the mean needs the engine) | kickmix-engine **16–18 s**, 0.55 GB, 1 salt | 16 salts ≈ 5 min engine; 48 salts + official at gates | 100–200 |
| chemistry | build + 0.1 s (exact per pass) | official verifier ≈ 45–66 s (engine gives no gain) | all-branch proof + verifier | 40–80 |
| Sub-program benchmarks (loader truth table, GF(2^12) gadget, pass order on a fixed circuit) | — | exhaustive bit-parallel check, **ms to 1 s** | splice into the generator + T1/T2 | 1,000s (proposer-bound) |

With an API proposer (4–8 concurrent calls, 10–40 s each) the evaluator is the bottleneck on whole-task
candidates; with a local proposer (≈ 3–6 min per proposal, §2.3) the proposer is the bottleneck at 10–20
proposals per hour. A January window of ≈ 10 working days × 12 h gives at most ≈ 1,200 h·slot, i.e. 2,000–4,000
realistic whole-task evaluations, or far more on sub-program benchmarks.

### 1.3 How many proposals published evolve-style systems needed

| System / result | Proposer | Proposals or samples | Baseline quality | Gain |
|---|---|---:|---|---|
| FunSearch, cap sets / bin packing (Nature 2023) | PaLM 2 code model | ≈ 10^6 samples per run | academic heuristics | new bounds |
| AlphaEvolve (2025) | Gemini ensemble | thousands to tens of thousands per problem | literature | 4×4 complex matmul, kernels +23%, ~20% of open math problems improved |
| ShinkaEvolve (Sakana 2025) | frontier ensemble, bandit | **≈ 150** (circle packing) | AlphaEvolve's own result | new SOTA |
| CodeEvolve (arXiv:2510.14150) | **open Qwen3-Coder-30B** | matched budgets (≈ 10× cheaper than frontier) | AlphaEvolve circle packing | surpasses AlphaEvolve on 2 instances; beats OpenEvolve/Shinka on 6 of 9 tasks |
| Gideoni, Risi, Gal (arXiv:2602.16805) | various | various | — | simple baselines match or beat evolve pipelines; **the search space and the evaluator matter more than the pipeline** |
| Fisher et al. (arXiv:2603.29894), LLM-guided VarTODD | frontier | not stated | best known T-counts | −7.0% mean T-count on GF(2^n) multipliers |
| Gu et al. (arXiv:2610.22367), OpenEvolve + Claude Mythos 5 | frontier | not reported; point-addition run 53 h | literature (Khattar et al. etc.) | 4–39×; RS DQI decoder 5.1e6 → 4.4e5 Toffolis |
| Forum 753185 (alex chilton): ShinkaEvolve + local Qwen 30B on the chemistry gauge optimiser | local 30B | 15 generations | crowd-optimised | **0 improvements** (Gemini by hand improved on its first try) |
| Forum 752949 (Valerio Dolci): own harness + local Qwen3.8-Flash-Next on 2× RTX 5070 Ti | local 125B-A6B MoE | "session after session through the night" | public chemistry (Oct 2026) | ≈ −5% chemistry over the best public circuit |

Reading: the large published gains are against literature baselines, from frontier proposers, and need
hundreds to thousands of proposals. On a crowd-optimised incumbent, a small proposer either finds nothing
(753185) or finds gains only when the harness is narrow and exact (752949). Gu et al.'s principles (conditional
routing, history recovered from the state, allowed failure, bigger tables instead of multiplications) were
already harvested by the crowd within weeks of the paper.

### 1.4 Expected gain and P(success) per task (est., January campaign of 2,000–4,000 proposals)

Hit-rate model (est.): 30–60% of diff proposals build and pass T1 with a good cheat-sheet (10–30% for a 9B
local model); 0.1–0.5% of all proposals beat the incumbent on their stage, each by 0.05–0.5% of task cost.
Score gain ≈ task gain / 3.

| Task | Good targets (local objective, fast exact check) | Poor targets | Task-cost gain if it works | P | Score gain, expected |
|---|---|---|---|---:|---:|
| invsqrt | K = 6 loaders and other QROM/loader programs (exact over 64–128 addresses, ms); one patch module at a time under the qubit cap with `prove_any`; finishing-pass order | the whole datapath (13 patches rest on full-domain facts of the 128-row table: S2.9 pause); angle refit (exhausted per Asterios); table search (numeric: use the exact integer model + search, not an LLM) | 0.5–3% | 0.25 | +0.1 to +0.25% |
| dqi | GF(2^12)/GF(2^24) multiplier and closing gadgets (exhaustive over 2^24 pairs bit-parallel in < 1 s; DQI-A02 shows 360 → 340 Toffolis per site exists); unary-iteration / extractor emitters | decoder control flow and cutoffs (cost depends on multi-salt failure rates: noisy, 16+ salts per decision); qubit layout of the 890-qubit circuit | 0.5–2% | 0.20 | +0.05 to +0.15% |
| chemistry | loader-compiler heuristics (4,022 Toffolis per pass, exact all-branch proof ≈ 1 min); codebook/address-split heuristics (FunSearch-style: evolve the heuristic, GEN-B09) | λ_eff (0.25% above a hard floor), the cost model itself (fixed), float64 refits (GPU numeric) | 0.3–1.5% | 0.30 | +0.05 to +0.15% |
| **All three** | | | | ≈ 0.25 that any lands | **≈ +0.3 to +1% score if it lands; +0.1 to +0.3% expected** |

A frontier proposer (Opus/Fable class) for the ideation operator (new family per stage) is the only route to
the large steps (the 28% cut the leaders need on one task); its cost is ≈ $1 per agentic iteration and its P is
low (INV-B01: 0.08). That is GEN-B11 + human review, not a bulk loop.

### 1.5 Effect of the ThemBones generators (Apache-2.0)

Public generators of the #2 team give (a) new **families** to seed islands with (the loop's weakest step is
ideation, and a strong seed removes it), (b) new stage-level **baselines** that are less crowd-optimised than
v12's oldest stages, and (c) more "archive alternatives" for prompts. They also raise the bar: a stage that a
top team already tuned is a worse target than one nobody touched. Recommendation: once the frontier report
lands, pick loop targets in the stages where ThemBones and v12 **differ most** in cost (diff the stage maps),
and credit ThemBones with SHA-256 in every derived genome.

### 1.6 Comparison with the current human + Claude workflow

| | Human + Claude Code (today) | Bulk loop (Haiku/Sonnet API) | Local loop (Qwen3.6-35B-A3B) |
|---|---|---|---|
| Attempts per day | 3–10 deep attempts (a plan step takes hours) | 500–2,000 shallow proposals | 150–300 shallow proposals (overnight) |
| Depth per attempt | high: reads papers, derives invariants, writes proofs | low: one block, one idea | lower |
| Cost | subscription quota (binding) | $2–40/day, separate from the quota | ≈ €1/day electricity |
| What it finds | new families, re-orientations (e.g. S2.9 pause), exact proofs | local improvements on narrow blocks; rediscovery of known tricks | same, at lower valid rate |
| Failure mode | quota exhaustion; single-lineage convergence (forum 755303) | reward hacking, near-duplicate edits, noise chasing on dqi | 0 hits (forum 753185) |

The two are complements: the loop replaces the mechanical "try the 20 obvious variants" part of a step, which
is exactly where the quota goes today; the human + Claude session keeps the ideation, the gates and the
decisions.

## 2. Best local open-weight model for this laptop

### 2.1 Constraints

- GPU: RTX 4060 Laptop, 8,188 MiB (checked with `nvidia-smi`, driver 595.79), ≈ 256 GB/s.
- RAM: 15.7 GB total. Windows needs ≈ 3–4 GB; WSL (evaluator) needs 2–4 GB while a loop runs. **Run the LLM
  server natively on Windows** (llama.cpp CUDA build or LM Studio), not inside WSL: the WSL cap is 8 GB and
  raising it starves Windows. The WSL harness calls `http://<windows-host>:8080/v1` (OpenAI-compatible).
  Model budget: ≈ 7 GB VRAM + ≈ 8–10 GB RAM → **≤ 15–16 GB of weights including context**.
- CPU: i5-12450HX, 8 cores / 12 threads, DDR5. MoE expert offload uses CPU threads and memory bandwidth;
  every thread given to the LLM is taken from kickmix-engine. Pin the LLM to 6 threads, the evaluator to 6.

### 2.2 Candidates

| Model (release, licence) | Total / active | Coding benchmarks (vendor-reported) | Quant that fits | File size | VRAM fit | Expected tok/s gen (laptop) | Verdict |
|---|---|---|---|---:|---|---|---|
| **Qwen3.6-35B-A3B** (Apr 2026, Apache-2.0) | 35B / 3B, 262k ctx | SWE-bench Verified **73.4**, LiveCodeBench v6 **80.4**, Terminal-Bench 2.0 51.5 | UD-IQ3_XXS (UD-IQ3_S 13.7, UD-Q3_K_S 15.4 if RAM allows) | **13.2 GB** | attention + shared + some experts on GPU, rest in RAM (`--n-cpu-moe`) | **15–30** (8 GB laptop reports: 33–40 tok/s with 32 GB RAM; 8.6 tok/s without the flag) | **#1** |
| **gpt-oss-20b** (Aug 2025, Apache-2.0) | 21B / 3.6B, 128k ctx | SWE-bench Verified 60.7 (high reasoning); LiveCodeBench v6 74.6 (Qwen's table) | native MXFP4 | ≈ 13 GB | partial expert offload | 25–40 | **#2** (mature tool calling, older) |
| **Qwen3.5-9B** (2026, Apache-2.0) | 9B dense (hybrid attention) | LiveCodeBench v6 65.6 | Q4_K_M | ≈ 6 GB | **fully on GPU** (7.1 GB used) | **33** measured on RTX 4060 (desktop) | **#3** (fast, leaves the CPU to the evaluator) |
| Gemma 4 26B A4B (2026) | 26B / 4B | LiveCodeBench v6 77.1 | Q3/IQ3 | ≈ 12–14 GB (est.) | partial offload | 15–25 (est.) | alternative to #1 |
| Qwen3.8-27B (Aug 2026) | 27B dense | LiveCodeBench v6 **90.3**, SWE-bench Pro **61.7** (best small open model) | Q4_K_M 17–18 GB; Q3 ≈ 13 GB | 13–18 GB | no: 24 of 58 layers on GPU | **≈ 3.5** (Qwen3.5-27B measured 3.57 on RTX 4060) | too slow locally; **Kaggle P100/T4×2 candidate** |
| Devstral Small 2 24B | 24B dense | SWE-bench Verified 68.0 | IQ4_XS | 12.8 GB | no | 3–5 | too slow |
| Qwen3.8-Flash-Next (Aug 2026) | 125B / 6B (+51B N-gram embeddings in RAM) | SWE-bench Pro 62.5, LiveCodeBench v6 91.9 | Q4 ≈ 70 GB (est.) | — | no | — | does not fit (needs ≥ 64 GB RAM) |
| GLM-5.3-Flash / DeepSeek-V4(.1)-Flash | 320B / 18B; 284B / 13B | Terminal-Bench 2.1 84.3 / — | — | 150+ GB | no | — | API only |

Reference point: **Claude Haiku 5.5** (API, $0.10/$0.50) reports SWE-bench Pro 64.8 and SWE-bench Multilingual
83.7 (system card via apidog), i.e. at least the level of the best open models, at 5–10× the local speed.

### 2.3 Throughput per proposal (est.)

| Setup | Prompt processing | 15k-token prompt | 5k-token output | Proposals/hour (one stream) |
|---|---:|---:|---:|---:|
| Qwen3.6-35B-A3B IQ3_XXS, `--n-cpu-moe`, 32k ctx | 300–800 tok/s | 20–50 s | 170–330 s | **10–18** |
| gpt-oss-20b MXFP4, partial offload | 400–1,000 tok/s | 15–40 s | 125–200 s | 15–25 |
| Qwen3.5-9B Q4_K_M, all on GPU | 1,500–2,500 tok/s | 6–10 s | ≈ 150 s | 20–25 |
| Haiku 5.5 API, 6 parallel streams | — | — | — | 300–600 |
| Kaggle T4×2, Qwen3.6-35B-A3B UD-Q4_K_XL (22.4 GB split) | 1,000+ tok/s (est.) | ≈ 15 s | 40–60 tok/s → ≈ 100 s | 25–35 |
| Kaggle P100 16 GB, Qwen3.8-27B Q3 (≈ 13 GB) | 500–800 tok/s (est.) | ≈ 25 s | 15–25 tok/s → 200–330 s | 10–15 |

Reasoning ("thinking") output dominates; cap it (`max_tokens`, or Qwen's non-thinking mode for mechanical
operators) and use grammar-constrained output (llama.cpp GBNF / JSON schema) so the reply is a parseable diff.

### 2.4 Serving stacks

| Stack | MoE CPU offload | Fit here | Notes |
|---|---|---|---|
| **llama.cpp `llama-server`** (Windows CUDA build) | **yes** (`--n-cpu-moe N`, `-ot` regex) | best | OpenAI-compatible API, GBNF/JSON grammars, KV quantization (`-ctk q8_0`), speculative/MTP support |
| LM Studio | yes (llama.cpp backend, GUI toggle) | good | easiest setup; same engine; headless server mode |
| Ollama | partial (automatic split, less control) | acceptable | simplest CLI; worse tok/s on 8 GB MoE in reports |
| vLLM | no useful CPU offload | no (laptop); T4 only | needs weights in VRAM; AWQ/GPTQ on T4 (compute 7.5); **P100 unsupported** (compute 6.0) |
| ExLlamaV2/V3 (EXL2/EXL3) | no | only for ≤ 8 GB models (9B) | fastest GPU-only decoding; no help for 35B-A3B |

Quantization formats: GGUF K-quants/IQ-quants (llama.cpp; Unsloth "UD" dynamic quants keep sensitive layers at
higher precision, the best choice at 3 bits); MXFP4 (gpt-oss native); AWQ/GPTQ (vLLM, GPU-only); EXL2/EXL3
(GPU-only). At 8 GB VRAM only GGUF with expert offload makes a 30B-class MoE usable.

### 2.5 Top-3 and the honest trade-off

1. **Proposer of choice: Claude Haiku 5.5 via API** (not local): $0.004 per proposal, best coding scores in
   its price class, fast, prompt caching; Sonnet 5.5 for the 20% of proposals that need more reasoning.
   Billing is per token on an API key with a spend limit, independent of the Claude Code subscription.
2. **Local fallback #1: Qwen3.6-35B-A3B UD-IQ3_XXS (13.2 GB download)** on Windows llama.cpp, 15–30 tok/s,
   10–18 proposals/hour, for unattended nights at $0.
3. **Local fallback #2: Qwen3.5-9B Q4_K_M (≈ 6 GB)**, all on GPU, ≈ 33 tok/s, lowest interference with the
   evaluator; use for mechanical operators (parameter tweaks, pass orderings) where code quality matters less.
   gpt-oss-20b (≈ 13 GB) is an equal alternative to #2 if Qwen's 3-bit quality disappoints.

Kaggle GPUs (30 h/week) can host a better open model (Qwen3.8-27B on P100, or Qwen3.6-35B-A3B at Q4 on T4×2) and
the whole loop (kickmix-engine on 4 vCPU), but the P100 is already booked for chemistry refits (M3.5), sessions
end after 12 h, and the gain over Haiku 5.5 is nil on quality. Not recommended except as a $0 alternative if
the user does not want any API spend.

## 3. Fine-tuning analysis

### 3.1 Feasibility by method and hardware

| Method | What it needs | Laptop (8 GB, no room in 16 GB RAM) | Kaggle (P100 16 GB / T4×2 16 GB each, fp16 only, 30 h/week) | Data we would have by mid-January |
|---|---|---|---|---|
| Full fine-tune 9B | ≈ 16 bytes/param ≈ 150 GB | no | no | — |
| bf16 LoRA | Qwen3.5-9B ≈ 22 GB; Qwen3.5/3.6-35B-A3B ≈ 74 GB (Unsloth) | no | no (no bf16 on T4/P100; 16 GB per card) | — |
| QLoRA 4-bit | 7–9B dense ≈ 6–10 GB | 7–8B older dense only (Qwen2.5-Coder-7B, Qwen3-8B) | 7–14B dense on one T4 (fp16 compute) | — |
| QLoRA on Qwen3.5/3.6 | — | **discouraged by Unsloth** (quantization error larger than usual; bitsandbytes limits MoE 4-bit) | same | — |
| SFT / rejection-sampling fine-tuning (expert iteration) on verified improvements | ≥ 1,000 good (prompt, accepted diff) pairs | as QLoRA | as QLoRA | **tens** of accepted improvements; a few hundred "valid but not better" edits |
| DPO on accepted vs rejected edits | ≥ 1,000–5,000 pairs | as QLoRA | as QLoRA | pairs are plentiful for *validity* (build fails vs passes), scarce for *improvement* |
| GRPO / RL from verifier reward | group of 8–16 rollouts × thousands of prompts; each rollout = one build + T1 (15 s – 3 min) | no (generation + training + evaluator do not fit) | ≥ 10k rollouts ≈ 50–150 h of evaluator alone, plus slow fp16 generation | reward is 0 for > 99% of rollouts (no improvement), so the gradient signal is near zero |

### 3.2 Expected benefit versus cost

- **What fine-tuning could teach:** the kickmix builder API, the register contract, MBU/HMR idioms, our patch
  conventions. That raises the *valid* rate (est. from 30% to 50% for a small model). It cannot teach the
  *improving* moves, because we have almost no examples of them; those are exactly what we are searching for.
- **What in-context examples give for free:** a cached prefix with a 2–3k-token kickmix cheat-sheet, 3–5
  verified edits retrieved from the archive (same stage, nearest cost), and the relevant nulls. With prompt
  caching on Haiku 5.5 this prefix costs about $0.00005 per proposal. Forum experience (753741: windowed edits
  are more reliable than whole-file rewrites) points to the same lever: narrow the edit window, enrich the
  context.
- **Calendar:** 9 days left before the break (1 Dec), then ≈ 12 working days (4–16 Jan) before the freeze.
  Collecting a dataset, training, evaluating the tuned model and running it would consume the whole January
  window for, at best, a higher valid rate on a model that is still weaker than Haiku 5.5.

**Verdict: no fine-tuning before 2027-01-20.** Reopen only if (a) the loop runs on a local model by necessity,
(b) it has produced ≥ 2,000 valid edits and ≥ 100 accepted ones, and (c) the local valid rate is the measured
bottleneck. None of these is likely before the deadline. Record it as a dated exclusion (GEN-B22 below).

## 4. Recommendation

### 4.1 Minimal architecture

| Component | Choice | Why |
|---|---|---|
| Loop framework | OpenEvolve (Apache-2.0; islands, MAP-Elites, diff mode, cascade evaluation; used by Gu et al.) **or** a ≈ 400-line custom loop on top of `tools/gate.py` and `tools/kmb.py` | Gideoni et al.: pipelines are about equivalent; what matters is the search space and the evaluator. Use OpenEvolve for the 1-day test (less code), keep the option to replace it |
| Proposers | **Haiku 5.5** for 80% of proposals (mutations, rewrites of one block), **Sonnet 5.5** for 20% (ideation of a new block design, crossover, any proposal after 3 failed children of the same parent); ShinkaEvolve-style bandit between them on "improvements per $" | cost: $0.004 vs $0.08 per proposal; escalation keeps quality where it pays |
| Local fallback | Qwen3.6-35B-A3B UD-IQ3_XXS on Windows llama.cpp, used only for unattended nights and only after the go/no-go shows its valid rate is ≥ half of Haiku's | $0, but slow and CPU-hungry |
| Genome | one marked block (emitter function, patch module or heuristic) inside a frozen, SHA-pinned base; never a circuit file | prize rule (rebuild from source), small diffs, cacheable prompts |
| Evaluator | T0 `tools/kmb.py` static cost + metric checks (< 1 s) → T1 exhaustive sub-program check or kickmix-engine 1 salt (ms – 18 s) → T2 `prove_any` / 16 salts / all-branch proof for archive entries only → T3 official verifiers at gates only | the CLAUDE.md speed rule: a proxy screens, the full check gates |
| Archive and selection | MAP-Elites per family on (width, Toffolis), 1 qubit × 1% bins; 70% exploit (top cells) / 30% explore; reject near-duplicate diffs (normalised similarity > 0.9) before building | as `general.md` §5.3 |
| Anti-reward-hacking | verifiers, kickmix, oracles and the base generator read-only and SHA-checked; any candidate touching them = invalid; T3 always on the untouched official verifier | Shirobokov et al. (agents weakening specs), hengck23 (byte-for-byte copy) |
| Budget | experiment: **$25 hard cap** on the API key; production: ≤ $100 total for January; stop rule: 3 rounds without archive gain → freeze the family | spend cap enforced by the provider console, not by the loop |
| Logging | one row per proposal batch in `experiments/LEDGER.md`, a changelog entry per campaign, run time and peak memory recorded | CLAUDE.md |

### 4.2 First experiment (≈ 1 day): "loader golf", go/no-go

- **Target:** the two 6-address-bit loaders of the invsqrt K = 6 table (S2.8). Truth table: 64 addresses ×
  ≈ 20 output bits, fixed. Genome: a Python function that emits an AND/XOR program (kickmix builder calls or a
  tiny IR). Fitness: number of ANDs (Toffoli-like gates) per loader, then peak extra width. Check: all 64
  addresses exactly, bit-parallel, in milliseconds. Known landmarks: our constructive baseline **44 ANDs per
  loader**; lower bound 19–20 products beyond the monomials; target from S2.8 ≈ 27 + 29 ANDs.
- **Why this target:** exact and instant evaluator, a scalar local objective, a known baseline, a known bound
  and a known "good" region, and the result is reusable (M3.2 own invsqrt generator, chemistry/QROM loaders).
  It tests the loop itself, not our patience with a 3-minute build.
- **Protocol:**
  1. Morning (≈ 4 h): evaluator wrapper (truth table + checker + AND count), OpenEvolve config with the
     cheat-sheet and the 44-AND seed, API key with a $25 limit.
  2. Afternoon (≈ 3 h, unattended): 500 Haiku 5.5 proposals + 100 Sonnet 5.5 proposals (≈ $10), 4 islands.
  3. Same wall time, **non-LLM baseline**: our greedy/SAT row-wise synthesis with random restarts on the 12
     cores (Gideoni et al.: always race a simple baseline).
  4. Night (optional): 150 proposals from local Qwen3.6-35B-A3B (or Qwen3.5-9B) to measure the local valid
     rate and hit rate for free.
- **Go** if the loop reaches **≤ 36 ANDs** on at least one loader (−18%) within 600 API proposals **and** beats
  the non-LLM baseline's best by ≥ 2 ANDs. **No-go** if it stays above 40 or does not beat the baseline:
  record the null, keep GEN-B01 frozen, and spend January on the human + Claude workflow. Either way record the
  valid rate, hit rate, $ per accepted improvement and tok/s of each proposer (the numbers this report could
  only estimate).
- **Second calibration (only on a go, ≈ 1 day):** rediscover a known public improvement in our source, e.g.
  DQI-A02 (norm-preserving GF(2^24) closing, 360 → 340 Toffolis per site) on an isolated gadget with an
  exhaustive check. A loop that cannot rediscover a known 5.6% gadget gain will not find unknown ones.

### 4.3 Where it sits in the plan

- Phase 2 tail (before 1 Dec): only if a day is free after S2.14–S2.16; otherwise do not start it, because the
  break would leave it unattended.
- Break: nothing required. An unattended local-model run during the break is possible but not recommended
  (no review, laptop heat, nothing to gate it against).
- **Phase 3: M3.4a on the first day after M3.1 (≈ 5–6 Jan)**, then **M3.4b** alongside M3.2–M3.3 on a go,
  targeting the blocks of the retained families with the largest ThemBones-vs-v12 stage-cost gap.

## 5. Proposed IDEAS rows (not applied)

Scale as in `IDEAS.md` §1 (score gain %, P, effort h, priority = P × gain ÷ effort). The highest existing
cross-cutting B ID is GEN-B17, so the new rows start at GEN-B18.

| ID | Task | Idea | Gain as stated | Score gain % | P | Effort h | Priority | Status |
|---|---|---|---|---:|---:|---:|---:|---|
| GEN-B18 | cross-cutting | Cheap API proposer tier for the loop: Haiku 5.5 bulk (≈ $0.004/proposal) + Sonnet 5.5 escalation, prompt caching, hard spend cap; billed outside the subscription quota | makes GEN-B01 affordable: ≈ 0.3% eq. | 0.30 | 0.60 | 4 | 0.0450 | TODO |
| GEN-B23 | cross-cutting | Exemplar bank for prompts: kickmix cheat-sheet + retrieved verified edits + nulls, cached prefix (replaces fine-tuning) | higher valid rate (est. 30% → 50%): ≈ 0.2% eq. | 0.20 | 0.60 | 3 | 0.0400 | TODO |
| GEN-B21 | cross-cutting | Rediscovery calibration: the loop must rediscover a known public gain (DQI-A02 closing 360 → 340) on an isolated gadget before production | avoids a wasted campaign: ≈ 0.2% eq. | 0.20 | 0.70 | 6 | 0.0233 | TODO |
| GEN-B20 | invsqrt (+ chemistry loaders) | Loader golf: evolve AND/XOR emitters for fixed loader truth tables with an exhaustive ms checker; go/no-go of the whole loop (§4.2) | loaders 44 → ≤ 36 ANDs; ≈ 0.3% score once the K = 6 path is integrated | 0.30 | 0.30 | 8 | 0.0113 | TODO |
| GEN-B19 | cross-cutting | Local open-weight proposer: Qwen3.6-35B-A3B UD-IQ3_XXS on Windows llama.cpp (`--n-cpu-moe`), Qwen3.5-9B as fast variant; $0 overnight fallback | ≈ 0.1% eq. (same hits as Haiku at a 5–10× lower rate) | 0.10 | 0.40 | 4 | 0.0100 | TODO |
| GEN-B22 | cross-cutting | Fine-tuning or RL on verifier reward (LoRA/QLoRA, RFT, DPO, GRPO) | none before the deadline (§3) | 0.10 | 0.05 | 40 | 0.0001 | EXCLUDED (analysis 2026-11-21: hardware, data and calendar) |

Revisions proposed for existing rows:

| ID | Change | Reason |
|---|---|---|
| GEN-B01 | Gain 3% → **1%** (stage-level hits on crowd-optimised incumbents), P 0.35 → **0.25**, effort 40 → **30 h** (OpenEvolve + existing gate tools); priority 0.0262 → 0.0083; status "TODO, gated by GEN-B20" | §1.3–1.4: published gains are against literature baselines with frontier proposers; local 30B null on chemistry (753185) |
| INV-A08 | merge into GEN-B01 (one patch module = one evolve block) | same loop, same evaluator |
| INV-B01 | keep P 0.08; add: ideation step only with a frontier model (Sonnet/Opus) plus human review, not the bulk proposer | family-level redesign is out of reach of Haiku/local models |

## 6. Proposed plan changes (not applied)

| Where | Change |
|---|---|
| `PLAN.md` §9, M3.4 | Split into **M3.4a "Loop go/no-go"** (GEN-B18 + B23 + B20, 1 day, ≤ $25, go criterion of §4.2) and **M3.4b "Loop production"** (GEN-B01 with GEN-B02 operators on blocks of the retained families, ≤ $100, only on a go, frozen after 3 rounds without gain). Replace "fresh subagent per iteration" by "API diff proposals (Haiku/Sonnet); Claude Code subagents only for ideation and gates", to protect the subscription quota |
| `PLAN.md` §7 (Phase 2) | Optional extra step: run M3.4a before 1 Dec if a full day is free after S2.14–S2.16 |
| `PLAN.md` §12 (Risks) | Add: "API spend runaway" → provider-side spend cap; "evaluator reward hacking" → read-only SHA-pinned verifiers, T3 on the untouched official verifier; "local LLM starves the evaluator" → thread pinning, LLM on Windows, evaluator in WSL |
| `general.md` §5.5 | Superseded in part: workers become API diff proposals; the agentic subagent mode is kept for operators 7–8 (crossover, new family) only |
| `IDEAS.md` | Add rows GEN-B18–B23, revise GEN-B01, merge INV-A08 (§5) |

## 7. Sources

Project files: `CLAUDE.md`; `Analysis/PLAN.md` (§7, §9 M3.4); `Analysis/IDEAS.md` (GEN-B01, INV-A08, INV-B01);
`Analysis/research/general.md` §2.1 and §5; forum snapshots `Instructions/discussion/752949-*.md` (Valerio
Dolci, local model setup) and `753185-*.md` (ShinkaEvolve + local Qwen 30B null); generator sizes measured with
`find`/`wc` on 2026-11-21; GPU checked with `nvidia-smi`.

Papers and frameworks:
- Gu, Kemeny, Yelin, Lukin, "Discovering state-of-the-art quantum circuits with large language models",
  arXiv:2610.22367: https://arxiv.org/html/2610.22367
- AlphaEvolve, arXiv:2506.13131: https://arxiv.org/abs/2506.13131
- FunSearch, Nature 2023: https://www.nature.com/articles/s41586-023-06924-6
- ShinkaEvolve, arXiv:2509.19349: https://arxiv.org/html/2509.19349v1 and https://github.com/SakanaAI/ShinkaEvolve
- OpenEvolve: https://github.com/codelion/openevolve
- CodeEvolve, arXiv:2510.14150: https://arxiv.org/html/2510.14150
- Gideoni, Risi, Gal, "Simple Baselines are Competitive with Code Evolution", arXiv:2602.16805: https://arxiv.org/abs/2602.16805
- "Compute Allocation for Self-Evolving LLMs" (BaSE), arXiv:2605.29268: https://arxiv.org/abs/2605.29268
- Fisher et al., LLM-guided VarTODD, arXiv:2603.29894: https://arxiv.org/abs/2603.29894

Models, quantization, throughput:
- Qwen3.6-35B-A3B GGUF (Unsloth; sizes, licence, benchmarks): https://huggingface.co/unsloth/Qwen3.6-35B-A3B-GGUF
- Qwen3.6 SWE-bench write-up: https://dev.to/purpledoubled/qwen36-scores-734-on-swe-bench-with-only-3b-active-parameters-heres-why-that-matters-2fmp
- RTX 4060 8 GB measurements (Qwen3.5-9B 33 tok/s, 27B 3.57, 35B-A3B 8.61 without expert offload): https://dev.to/plasmon_imp/moe-beat-dense-27b-by-24x-on-8gb-vram-the-35b-a3b-benchmark-nobody-expected-24n
- `--n-cpu-moe` on 8 GB laptop GPUs: https://medium.com/coding-nexus/one-flag-tripled-my-tokens-per-second-on-an-8gb-gpu-238f55d1e3bd , https://windowsforum.com/news/llama-cpp-n-cpu-moe-run-qwen3-35b-faster-on-8gb-rtx-laptop-gpus.449607/?amp=1 , https://huggingface.co/blog/Doctor-Shotgun/llamacpp-moe-offload-guide
- VRAM tiers, Qwen3.8-27B / Devstral Small 2 / Gemma 4 figures (Sept 2026): https://shop.zimaspace.com/en-ca/blogs/tech-ai-hub/best-local-coding-models-8gb-12gb-16gb-24gb-vram
- Qwen3.5-9B model card (LiveCodeBench v6 65.6; gpt-oss-20b 74.6 in the same table): https://huggingface.co/QWEN/Qwen3.5-9B
- Qwen3.8 releases: https://raw.githubusercontent.com/QwenLM/Qwen3.8/main/README.md
- Qwen3.8-Flash-Next (125B-A6B, benchmarks): https://www.datacamp.com/blog/qwen3-8-flash-next
- GLM-5.3-Flash (320B-A18B): https://rits.shanghai.nyu.edu/ai/glm-5-3-flash-320b-multimodal-moe/
- DeepSeek V4/V4.1 Flash pricing: https://api-docs.deepseek.com/quick_start/pricing , https://venturebeat.com/ai/deepseek-v4-1-flash-debuts-with-0-003-1m-off-peak-cached-input-rate-and-benchmarks-eclipsing-gpt-5-6-sol-claude-opus-5
- gpt-oss-20b (SWE-bench Verified 60.7, MXFP4; from the model card, not re-checked today): https://openai.com/index/introducing-gpt-oss/
- Unsloth fine-tuning memory (bf16 LoRA 9B ≈ 22 GB, 35B-A3B ≈ 74 GB; QLoRA discouraged for Qwen3.5): https://unsloth.ai/docs/de/modelle/qwen3.5/fine-tune , https://unsloth.ai/docs/models/qwen3.6

Prices:
- Anthropic API pricing (Haiku 5.5 $0.10/$0.50, Sonnet 5.5 $2/$10, Opus 5.5 $4/$20; cache and batch
  multipliers; tokenizer note): https://platform.claude.com/docs/en/about-claude/pricing
- Claude Haiku 5.5 benchmarks (SWE-bench Pro 64.8, Multilingual 83.7): https://apidog.com/blog/claude-haiku-5-5-benchmarks/ , https://anthropic.com/claude-haiku-5-5

Caveats: most third-party model figures are vendor-reported or from secondary blogs and were not reproduced
here. Laptop tok/s for this exact i5-12450HX + 16 GB configuration is interpolated from 8 GB laptop reports
with 32 GB RAM (expect the lower end). Per-proposal hit rates are estimates that the 1-day experiment is
designed to measure.
