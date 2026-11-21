# Public frontier, 20–21 November 2026: the 31.68 → 33.10 jump, its methods, and why our 19 Nov studies missed it

Written 21 Nov 2026 (research step, no code or register changed). Sources, all read as untrusted text and never
executed: the notebooks pulled into `downloads/kernels-2026-11-21/` (one folder per notebook:
ThemBones' *Gift*, Almond Labs 31.70105, Rithvik Redrouthu's *Shut up and Calculate*, tonylica 32.85306,
goodpjw2008 and mizeroluckygall *Dead-Word Loader Chain*, boristown *DAG splicing*, Asterios Panoras 31.69052
*What we tried* and 33.06908 *credits*, mizeroluckygall 37Q split-block refit, ghazarosbarseghyan91 *error
clusters*, sjhageman *decoder phase*, cosmorin milestone/tie reproductions, arsen123 *independent audit*, Syed Asad
Ali's portfolio notebook), the ThemBones generator dataset (`dataset__thembones-milestone-generators/`, 29 MB,
README, NOTICE, build records), the metadata of Rithvik's, tonylica's and Syed's datasets
(`dataset__*-meta/`), the receipts in `downloads/outputs-2026-11-21/` (another process), and one static read of
the 37-qubit invsqrt with `tools/kmb.py` (< 2 s). Numbers are the authors' official-verifier receipts unless
marked *(estimate)* or *(ours)*.

Reference costs: invsqrt 2,076,774.23; dqi 24,439,463,583,672.9; dfthcsos 3,355,599,588.32.
Our current submission (28.234, 19 Nov): invsqrt `ecf908fe` 74,383.7 (44 Q × 1,690.5), dqi `e39b1396`
3.4526e11 (890 Q, 365,893 T), chemistry 294,642,170 (λ 46.9216 × T 14,706 × Q 427).

## 0. Summary

1. **The best public per-task files give 33.10** (geometric mean 33.1019, goodpjw2008 / boristown /
   arsen123 receipts): invsqrt `d2148a7c` 60,955.33 (37 Q), dqi `4b572141` 3.0504e11 (900 Q), chemistry
   252,534,925 (405 Q, T 13,290, λ 46.91821). Against our submission that is **−18.1% invsqrt, −11.6% dqi,
   −14.3% chemistry**, i.e. +17.2% score. The board leader is at **33.348**; one reading that fits it to the
   third decimal is the same set with the dqi Toffolis of `4b572141` at **890 instead of 900 qubits**
   ((900/890)^(2/3) × 33.10 = 33.35), section 4.3.
2. **All of it is in-family.** No new datapath appeared on any task. invsqrt is still Syed Asad Ali's
   root + divider "family C"; dqi is still the in-place Berlekamp–Massey / CRT decoder; chemistry is still
   DF-THC with QROAM angle loaders. The gains came from (a) width engineering on invsqrt (44 → 39 → 37 Q) with
   measurement-based undo and recomputed digits, (b) a chemistry shape that crosses a power-of-two boundary of
   the scorer (B 143 → 127) plus brot 8 → 7 made accurate by a rounding-aware fit, and loaders co-designed with
   the table ("wired" bits, dead words), (c) cheaper multi-output field-multiplication kernels in dqi.
3. **Why our studies missed it** (section 3): our invsqrt floor counted owners at the phase endpoint as if the
   quotient had to be held whole; our chemistry Q proof was stated only for B ∈ [129, 256] and never ran the
   1-minute scorer grid it proposed; brot 7 was declared accuracy-dead from roundings of 8-bit-optimal fits;
   dqi operand-sharing kernels were priced at ≤ 0.5% with P 0.10. Each null had a domain restriction that the
   public work stepped outside of.
4. **Source availability is uneven** (section 1.3): the ThemBones generators (all three tasks, Apache-2.0, with
   byte-for-byte rebuild records) are the only complete public source line. The 37-qubit invsqrt and Rithvik's
   890-qubit dqi have **no public generator**. tonylica's dqi and chemistry rebuild from *recovered KMX assembly*
   plus small generators. Almond Labs' edits ship as binaries and descriptions.
5. **Recommendation** (section 5): resubmit the 33.10 public set at once (credited), and refactor the plan onto
   the ThemBones source line per task: rebuild it in WSL, re-implement the documented public edits as generator
   passes, then spend the remaining Phase 2 days on the two levers with the best price per point: invsqrt width
   (one qubit = 2.7% task cost ≈ 44 Toffolis) and dqi width at the tonylica Toffoli level (one qubit ≈ 0.22%).

## 1. Timeline and the best public files

### 1.1 Timeline (UTC)

| When | Who | What | Board / local score |
|---|---|---|---|
| 19 Nov (our snapshot) | public frontier | 44 Q invsqrt, 890 Q dqi `e39b1396`, Lera/maoxin123 chemistry at Q 427 | 28.234 (us) |
| 20 Nov 04:24 | ThemBones (Danish Ali) | set s108 uploaded: 39 Q / 901 Q / 405 Q | **31.67924** |
| 20 Nov 14:05–16:03 | ThemBones | circuits dataset, generator dataset v2 (`danishlynxthembones/vqcc-thembones-milestone-generators`), *Gift* notebook | — |
| 20 Nov 17:33 | Rithvik Redrouthu | 37-qubit invsqrt `d2148a7c` first verified (cosmorin's log) | — |
| 20 Nov 19:50–23:57 | Almond Labs | LEND29 + T-fusion (invsqrt), livescan/LDH + INV17 + W=3 scan (dqi), recompiled loaders + GRG (chemistry) | 31.69054 → 31.69692 → 31.70105 |
| 20 Nov 22:31 | Syed Asad Ali | portfolio v17 with generators: 43 Q invsqrt 72,440; dqi 890 Q 3.368e11; chemistry 295.0M | 28.71 (local) |
| 20 Nov 23:39–23:59 | Rithvik Redrouthu | `qcc-milestone-set` (two packs, no source): 37 Q invsqrt, 890 Q dqi, chemistry 253.29M | 32.86114 (board) |
| 21 Nov 00:24–00:32 | tonylica (Tony Li) | 32.85306 *task generators* (KMX + builders): 39 Q invsqrt `68eb8c13`, 900 Q dqi `4b572141` | 32.85306 |
| 21 Nov 01:32 | Asterios Panoras | best file per task (Rithvik invsqrt + chemistry, tonylica dqi) | 33.06908 (tie with 3 teams) |
| 21 Nov 06:00–09:00 | boristown, rupam2020, ozertuu, goodpjw2008, areebbb, adilshamim8, alexchilton | dead-word loader chain (13,328 → 13,292 → 13,290), loader 6 recompiled, ozertuu's table | 33.1002 → 33.10186 / 33.10189 |
| 21 Nov ~10:00 | board | leader 33.348; tie of ≈ 100 teams at 33.10 | — |

### 1.2 Per-task cost of the best public files

**invsqrt** (cost = Q × (Toffolis + T/4 + rotation price); 33 rotations ≈ 237 Toffoli-eq in every file):

| SHA-256 | Author | Q | Mean Toffolis | Static | T | Toffoli-eq C | Cost | Generator source |
|---|---|---:|---:|---:|---:|---:|---:|---|
| `ecf908fe` (ours, 19 Nov) | public 44 Q line | 44 | 1,450.9 | — | — | 1,690.5 | 74,383.7 | compiled only |
| `70fe5f79` | Syed Asad Ali v17 | 43 | 1,436.09 | — | — | 1,684.7 | 72,440.35 | **yes** (`invsqrt_generator.zip`, Apache-2.0) |
| `c639a045` | ThemBones (b196lend) | 39 | 1,379.375 | 1,493 | 3 | 1,616.80 | 63,055.04 | **yes**, rebuilt byte for byte twice (`invsqrt/tools/f9_i1_stack.sh`, ~20 min, 1 core) |
| `0612d079` | Almond Labs (LEND29 + T-fusion) | 39 | 1,378.875 | 1,492 | 1 | 1,615.80 | 63,016.05 | no (binary in the notebook; both edits described exactly) |
| — | jeffemerson (cited by boristown) | 39 | — | — | — | — | 62,977 | not checked |
| `68eb8c13` | tonylica (v90 → Geremie Yeo v50 divider → jackvd v220 carry, LEND29) | 39 | 1,339.9 | — | — | 1,581.9 | 61,696.23 | recovered KMX `base_v90.kmx` + 4 small Python patches |
| `d2148a7c` | Rithvik Redrouthu | **37** | 1,406.31 | 1,553 | 8.5 | 1,647.44 | **60,955.33** | **no** (dataset carries zips only; NOTICE credits Syed + ThemBones) |

Static read of `d2148a7c` *(ours, `tools/kmb.py`)*: 37 qubits, 78,869 bits, 60,767 ops (the 39 Q files have
≈ 24,600), 1,468 CCX + 39 CCZ + 42 CSWAP, **1,278 HMR**, 4,054 Z_IF + 3,280 CZ_IF, 2,271 condition frames,
55 Z_POW. The 37 Q file pays +27.4 Toffolis a shot over the 39 Q one (+1.9%) to save 2 qubits (−5.1%).
Worst error / tolerance 0.99947 (arsen123): no precision slack.

**dqi** (cost = R × Q² × T; R = 1.19126 at 0 failed shots, +0.46% per failed shot):

| SHA-256 | Author | Q | Mean T | Static T | Failures (local) | Cost | Generator source |
|---|---|---:|---:|---:|---|---:|---|
| `e39b1396` (ours) | public 890 Q line | 890 | 365,893 | — | 0 | 3.4526e11 | compiled only |
| `a31716fa` | Syed Asad Ali v17 | 890 | 356,981 | — | 0 | 3.3685e11 | **yes** (`dqi_generator.zip`) |
| `2c043a9a` (file A), tickets `859ff092`… | ThemBones (chain zd20) | 901 | 343,998.0 | 371,577 | 0/0 on 2 salts | 3.32669e11 | **yes** (kit + 11 chain modules; file A built once, rebuild not re-run) |
| `c4cbeafd` | Almond Labs | 901 | 343,806.53 | 373,287 | 0/0 | 3.32484e11 | partial (C++17 rewriter embedded in the notebook, applied to file A) |
| (not pinned) | Rithvik Redrouthu | **890** | ≈ 328,000 *(estimate from cost)* | — | — | 3.09449e11 | **no** |
| `4b572141` | tonylica (compact PA001 parent + 4 D16 sites) | 900 | 316,131 | 318,972 | 0 on salt ""; 2/102,400 on arsen123's salt (3.08e11) | **3.05042e11** | build from recovered KMX + `selector.py` + `arithmetic.py`; "the upstream optimizer that chose the formulas was not recovered" |

**chemistry** (cost = λ_eff × T_total × Q_total; T_total = 6,500 fixed + 2 × loader Toffolis per pass):

| Files | Author | (R, B, C), brot, J | Q | Loaders / pass | T | λ_eff | Cost | Source |
|---|---|---|---:|---:|---:|---:|---:|---|
| ours (19 Nov) | Lera table + maoxin123 loaders | (8, 143, 44), 8, 5 | 427 | 4,022 | 14,706 | 46.9216 | 294,642,170 | compiled loaders |
| npz `e05fb18f` | ThemBones | **(9, 127, 46), 7, 5** | **405** | 3,468 | 13,436 | 46.932071 | 255,384,621 | loaders: `.kmx` + compiler tools (randomized, carried); table: GPU fit, carried |
| — | Almond Labs | same | 405 | 3,464 | 13,428 | 46.918298 | 255,157,659 | 4 loaders recompiled (binaries), weights refit |
| — | tonylica | same | 405 | 3,465 | 13,430 | — | 255,214,794 | `build.py` from KMX (9 ThemBones, 6 disassembled Almond) + TSV tables |
| npz `e259454d` + angle_0 `8df36485` | Rithvik Redrouthu | same | 405 | 3,414 | 13,328 | 46.923700 | 253,286,627 | no |
| npz `fc4b4c01` + chain | rupam2020 (9 loaders) + Rithvik (5) + boristown (k6 `64c95e3c`) + ozertuu (table) | same | 405 | **3,395** | **13,290** | 46.918211 | **252,534,925** | loaders as `.kmb`; chain DP in goodpjw2008's notebook |

Fixed part of T at the new shape *(scorer receipt)*: outer prep 337, inner prep 4,213, Givens 1,950
(old shape: 317, 4,095, 2,250). Qubits: persistent 118 (outer 40, inner 43, words J·brot = 35), temporary 287
(old: 84 + 40 + 303).

### 1.3 Licences and generator sources (prize eligibility)

Every file above is published under Apache-2.0 with NOTICE files (ThemBones `NOTICE` + `third_party__*`;
Rithvik's `NOTICE`; tonylica's `notices/RELEASE_NOTICES.txt` 81 kB; Syed's NOTICE, which additionally asks
that he be named first as original builder in any derivative's credit line: cheap to honour, we should). The
rules require winners to license the submission **and the source code that generates it** under Apache-2.0 and
to provide a repository that reproduces the results. Status per base:

| Base | What is public | Reproducible from source? | Risk for a prize submission |
|---|---|---|---|
| ThemBones invsqrt `c639a045` | full generator tree `src_s18_srcM` + patch + staged program + config + angles; build record of two identical rebuilds | **yes** (srcM is carried in its final state; the patch history from Syed's generator is not replayed) | low |
| ThemBones dqi 901 Q | kit `kv_unfused5` (rebuilt byte for byte from `kv_unfused4_z2`), 29 extension modules, oracle, re-roll tool | **mostly**: file A was built once (40 min, 2 cores) and never rebuilt; the parent chain's 693 MB file is not included (used only for the pair check) | medium (we must rebuild it ourselves once) |
| ThemBones chemistry | 15 loader `.kmx`, table `all_c.npz`, `asm_words.npz`, fit and compiler tools (JAX/GPU) | loaders from `.kmx`: yes; `.kmx` and table: carried as data (randomized searches, GPU race) | medium (a carried table + KMX is "data + builder", like Syed's chemistry) |
| Syed v17 (all three) | generator zips with build receipts | yes | low, but 3–10% behind |
| tonylica dqi `4b572141` / invsqrt `68eb8c13` / chemistry | builder + recovered KMX assembly + formula tables | builds byte for byte, but the KMX is disassembly, not generated; the formula optimizer is missing | medium–high |
| Rithvik 37 Q invsqrt, 890 Q dqi, chemistry | binaries only | **no** | high: usable on the board with credit, not as a final prize base unless we regenerate it |
| Almond Labs edits | binaries + exact descriptions (affine relations, op indices) | re-implementable from the description | low once re-implemented |
| dead-word loader chain (rupam2020, boristown) | `.kmb` loaders, chain DP | loaders as binaries | medium (disassemble to KMX like tonylica did) |

Several forum topics already discuss which public DQI decoders are prize-eligible (`Instructions/discussion/`);
a host answer on "recovered assembly as source" would settle the medium cases and is worth asking for.

## 2. The methods, per task

### 2.1 invsqrt: 44 → 39 → 37 qubits

**Datapath (unchanged family).** ThemBones' "family C" is a fork of Syed Asad Ali's public generator: shift
normalisation by an even amount, a short **in-place (consuming) square root**, **two correction divisions**, a
higher-order stage built on small table words, the main divider, and a **cyclic endpoint** whose router feeds
the 33 Z-rotations (angles from an LP that maximises the worst-case margin; 34 binding points between x = 89
and 2.9e6, so 33 is a proven floor for this datapath: Almond Labs' LP is infeasible for every single drop).
Then everything is undone in reverse. Syed's own public v17 already had the three structural ideas (in-place
root, quotient signs kept in freed wires, cyclic routing into the rotations) and reached 43 Q at 1,436 T.

**What took 43/44 Q to 39 Q** (ThemBones README + build log `invsqrt/record/stack_dq18c.log`):

1. **Undo paid by measurement wherever it is exact.** A temporary AND or a quotient digit is erased by an
   X-basis measurement (HMR); its phase is repaid by CZ/CCZ conditioned on the outcome, placed *inside a
   later block that computes the needed product anyway*, so the repair costs no extra Toffoli.
2. **Digits of the main divider are not stored as history.** In the undo they are "re-made by short exact
   programs that borrow no wire" (per-block programs found by search: the log's `game` DP at `L 8, Cc 3`,
   i.e. 8 live, 3 clean carries). This is the general form of our INV-B10 parity-row idea, applied to every row
   whose digit can be recomputed from live state, not only to the 12 feed-free rows.
3. **Wire lending.** A wire whose value is a fixed Boolean function of other live wires over a window
   (`q11 = T5 + W5 + MAJ(T5 + T7, W3, W4)`, checked on all 3,133,452 inputs) is X-measured, its phase repaid by
   Cliffords, lent to a block as a clean carry, and re-made by **one AND** afterwards (−1.125 T a shot in that
   block). Almond's LEND29 is the same move with an affine relation (`q34 = q5 ⊕ q7 ⊕ q9 ⊕ q10 ⊕ 1`, 0 AND,
   −0.5 T a shot); a second lend on the same relation fails (1,417,293 failing shots). Almond's T-fusion uses
   `a + b − (a ⊕ b) = 2(a ∧ b)` with `a ∧ b` available as an early parity to turn two T gates into one S (T 3 → 1).
4. **Every block is fitted onto a fixed wire budget, then a global relabel.** The pipeline emits at 40 wires
   with a busy peak of 39, then a width step relabels to 39 (`wires before 40, peak busy 39, labels after 39`).
5. **Exact all-input checking in seconds.** Each emitted variant is compared with its parent on all 3,133,452
   inputs under four HMR-outcome modes (all zeros 1,327.0 T, all ones 1,425.0 T, two random seeds 1,379.4 T)
   in 5–8 s, then an exact affine pass and a rewrite chain. The full pipeline runs in 401–1,845 s on one core.

The result costs *fewer* Toffolis than the 44 Q line (1,379 vs 1,451) while removing 5 qubits: width and
Toffolis were not in tension, because the stored history was both a width cost and an undo cost.

**What took 39 Q to 37 Q** (Rithvik Redrouthu, no source, no write-up). What the file shows *(ours)*: the op
count grows from ≈ 24,600 to 60,767, with 1,278 HMR, 4,054 classically conditioned Z and 2,271 condition
frames; mean Toffolis +27.4 (+1.9%), static 1,553 (+60), T 8.5 instead of 1–3. That signature is
"recompute instead of keep" pushed further: more digits and temporaries measured out early and re-made later,
more phase repair, two wires fewer. At 37 Q one qubit is worth 1,647/37 ≈ 44.5 Toffolis, so +27 Toffolis for
−2 qubits is a net −3.3%. Asterios Panoras measured that the *minimal* way to delete the 39th wire of the
ThemBones file costs ≈ 449 Toffolis (a loss of 25%), so the 37 Q file is a re-planned schedule, not a local
deletion. tonylica's parallel line (`68eb8c13`, 39 Q, 1,339.9 T) shows a further −39 Toffolis is available at
39 Q from a different divider (Geremie Yeo's v50 rewrite: one rotation moved into the even remainder, the low
carry reinstated in both directions) and a reconstructed carry schedule (q29 borrowed for the carry chain).

### 2.2 dqi: 890 → 901 → 900 qubits, 365.9k → 344.0k → 316.1k Toffolis

**ThemBones 901 Q** (from Atakan Aldemir's public bundle, whose 892 Q carry decoder they rebuilt byte for
byte; lineage SomeJoint / yoyobar / manwithacat / Gidney / Busya PRIME / James Day):

- projective finish with CRT closing (after Syed Asad Ali's v13+ decoder);
- cached 12-bit field inverses at the converter's pivot sites;
- converter steps in which **one field product serves two terms**;
- **direct-product kernels for the two dot products of one word vector with both payloads** (bilinear
  schemes over GF(8) lifted through the points of its projective line): an operand-shared multi-output kernel;
- measured quotient heads at the restore;
- exact re-assignment of wires through proved-zero intervals, then an exact post-pass.

They accepted **+11 qubits** (+2.5% through Q²) for −6.0% Toffolis (365.9k → 344.0k): net −3.6% vs our
frontier file. Every kernel was checked as an identity on all inputs, every step paired against its parent on
12 salts. Failures: 0/102,400 on both local salts for three tickets.

**Almond Labs, −191.5 T a shot at 901 Q:** livescan + LDH (58,431 proved-zero idle intervals, repeated
quadratic products held on them, −32.8); ozertuu's INV17 subfield disjunction (two GF(2²) norm indicators
cannot both vanish on a nonzero GF(2⁶) input, so `a·b = a ⊕ b ⊕ 1`, 119 sites, −117.65); a W = 3 cubic-trace
Chien scan with 61 packed parity checks (−32.1); 19 clean-ladder guards + 4-pass condition-frame deferral
(−8.35); 2 degenerate CCZ → CZ (−2). Static anatomy (Almond's figure): extract 17.7k, synthesis 31.6k, convert
50.8k, scan 44.2k, unconvert 66.9k, pad_unsynthesis 56.0k, unextract 60.7k.

**tonylica 900 Q, 316,131 T** (−8.1% T and −1 Q vs Almond; −11.6% task cost vs ours). From the README of
`tonylica/quantum-ai-32-85306-task-generators`:

- a **compact "lambda-only" PA001 parent** (from submission 58622653): `selector.py` maps the input bits to
  4,024 twelve-bit residuals, sorts them with an odd-even sorting network, selects λ and clears μ, then applies
  the gauge, the root-1330 phase contribution and the final gauge restoration;
- **four D16 substitutions** in the field stages: each row of an input-form table defines two XOR parities
  whose product toggles the outputs chosen by an output mask (a bilinear multi-output formula); all four sites
  share one **70-term** formula plus an **84-term** companion block with site-specific output masks;
- everything else is the inherited decoder as recovered KMX assembly (stages listed in `stages.json`).

Release: 69,782,188 instructions, 318,972 static Toffoli-likes, so static ≈ mean (316,131): far fewer
measurement-conditioned gates than the ThemBones line (371,577 static vs 344.0k mean).

**Rithvik Redrouthu, 890 Q, ≈ 328k T *(estimate)*:** the Syed / Atakan / Areeb Ali 890 Q line with
ThemBones-level kernels; no source and no write-up.

Closed publicly: a trace-based (Berlekamp-trace) root splitting instead of the 4,026-element Chien scan costs
1.10M T vs 160k for convert + scan; a ≤ 19-AND GF(64) inverse and a rank-< 45 GF(16)[t]/t³ multiplication were
not found (160k constructive trials, 6.4e9 flip-graph steps); the two "dead tail" removals fail the verifier
(12/18 and 20/14 value errors); all 900 wires are linearly independent at 3,208 of 3,925 peak moments and live
≈ 91.7% of the run, so no wire is freed by CNOT rewrites alone.

**Tickets.** A DQI "ticket" is the same file plus one trailing `INFO_TAG reroll:<label>`: a new SHA-256, so a
new verifier draw, no gate change. Teams verify several and ship one with 0 failed shots. Locally it only picks
a lucky draw for salts "" and alt-salt-1; on the board the draw comes from the secret salt and the file's hash,
so a ticket is a ±0.03-point lottery, not a lever (Asterios Panoras).

### 2.3 Chemistry: Q 427 → 405, T 14,706 → 13,290, λ 46.92 → 46.918

**The shape moved the scorer's terms** (`scoreinfra/costs_prepare.py`, `SpinFreeData_V1`):

| Term | (8, 143, 44), brot 8 | (9, 127, 46), brot 7 | Why |
|---|---:|---:|---|
| b2 = ⌈log2 B⌉ + χ + 1 | 8 + χ + 1 | **7 + χ + 1** | B ≤ 128 crosses a power of two |
| inner temporary = ⌈log2(RC·⌈(B+1)/16⌉)⌉ + 15·b2 + uniform | 303 | **287** | −15 from b2, −1 from ⌈log2(B+1)⌉ − 2 |
| persistent (outer + inner) | 84 | 83 | b2 −1, ⌈log2 R⌉ +1, … |
| word registers J·brot | 40 | **35** | brot 8 → 7 |
| **Q** | **427** | **405** (−5.2%) | |
| Givens 2(N−1)(2·brot−1) | 2,250 | 1,950 | brot |
| inner prep | 4,095 | 4,213 | RC 352 → 414 |
| outer prep | 317 | 337 | RC |
| selection bits ⌈log2(76 + RB)⌉ | 11 (1,220 rows) | 11 (1,219 rows) | R·B kept ≈ constant |

So R·B was held at ≈ 1,143 rows (the selection register, the loader address space and the loader row count do
not change), B dropped below 128 and R rose to 9 to keep the accuracy, and C rose 44 → 46 to absorb brot 7.
Net of the shape alone: Q −5.2%, fixed T +138 −300 = −162.

**The fit made brot 7 accurate.** ThemBones' table is "a fit with hard sign constraints, our own Givens-tree
labels and staged rounding", from a JAX/optax GPU race (many GPU-hours) then a rounding-with-refit step
(`chem_rk7_table7c.py`, ~4 h): error_1 19.9959 mHa, error_2 19.9992 mHa, E_gap 3.9149010, λ_eff 46.932.
Balram Srinivasan's earlier null (brot 7 rounds to 27.1 mHa) and our NR-0124 (δ7 ≥ 19.3 mHa) were roundings
of 8-bit-optimal continuous fits after the fact.

**Loaders co-designed with the table.** 3,468 Toffolis a pass (−13.8% vs 4,022):

- fewer bit planes (brot 7: 35 output bits per loader instead of 40);
- **"wired" loaders**: during the rounding, 10 of the 15 loaders have chosen output bits forced to be parities
  of address bits (so those bits need no product at all), and 13 more bit relations among the words are kept;
  a finishing chain keeps every word in the signed range while both gates hold;
- exact plans for the first output groups (every hyperplane of a group's relation space, priced by Lovász's
  randomized matroid-parity rank), Syed Asad Ali's 68-AND high-cache form, a rank-additive cubic block in the
  last group of some loaders.

**After ThemBones (3,468 → 3,395 a pass):** Almond Labs recompiled 4 loaders (−1 each) using dead cells;
Rithvik Redrouthu and rupam2020 recompiled loaders against **dead words**: 2,400 of the 91,425 angle words sit
below a Givens rotation whose amplitude is 0 (exactly, or 1e-16 from cos(π/2)), so no metric can see them; the
amplitudes split eleven decades apart, so the criterion is unambiguous. Since loader k XORs batch k−1 ⊕ batch k
into the word registers, loaders from different families can be **chained** if every register error lands on a
dead word of every batch the verifier reads; goodpjw2008's shortest-path DP over 118 loader options from 41
public families finds −38 (13,328 → 13,290) and no more. Tables are interchangeable because every competitive
family loads the same live words: ozertuu's `fc4b4c01` gives λ 46.918211 on them.

Closed publicly: no loader Toffoli writes only dead words; none is replaceable by CNOTs (the delta is never in
the GF(2) span of the address bits and the held products); C 46 → 45 stalls at error_2 21.6 mHa with λ > 48.2
(goodpjw2008) and at a 20.351 mHa discrete plateau (Almond Labs); an augmented-Lagrangian refit of (W, H2_I)
gains nothing measurable (three active limits); a 257-Toffoli k0 fails because the verifier runs the loaders in
reverse as an involution (dirty carrier wire on the reverse pass).

Decomposition of our gap (294.64M → 252.53M, −14.3%): Q −5.15%, T −9.63% (Givens −300, fixed preps +138,
loaders −1,254), λ −0.008%.

## 3. Why our 19 Nov studies missed them

### 3.1 Claim by claim

| Our claim (file, section) | What happened | The assumption that was wrong |
|---|---|---|
| invsqrt: "best plausible target from our own source is F1 on a 44 Q base: −3% to −9%"; "no family reaches the −30% class" (`invsqrt_families` §1.5) | −15.2% at 39 Q (same family), −18.1% at 37 Q | We searched for **new datapaths** and treated in-family width moves (INV-A21–A24, GEN-B03/B04/B08) as a small residue at P 0.10–0.15. The family was not exhausted; its schedule was. |
| invsqrt floor 3: "the rotations need the whole 27-digit quotient at once; a divider design holds 57 owners raw, 43 after hiding v10/v11/e2/e3 and retiring history" (§2.3) | 39 and 37 Q with the same divider | Owners were counted as **stored** values. The public line stores almost no digit history: digits are measured out (HMR) and re-made in the undo by short exact programs; wires are lent when their value is a function of live wires. The endpoint holds what the router needs at each moment, not every digit with its history. |
| invsqrt: "only feed-free rows can be history-free" (impossibility lemma §5.2) | digits of feed rows are recomputed, not stored | The lemma is correct for *reversible rows that keep the remainder in an interval*; it does not forbid **recomputing** a digit from live state at undo time, which costs a few Toffolis but no wire. We bounded F1 by the lemma instead of pricing recomputation. |
| invsqrt rotation floor ≈ 33 (§2.3) | 33 in every public file; LP infeasible at 32 | **Held.** |
| invsqrt exchange rate in our pricing | one qubit = 44 Toffolis at 37 Q | We priced families in Toffolis and quoted width as a side constraint; at Q ≈ 40 width is the dominant lever (Asterios: +0.27 points per qubit vs +0.0065 per Toffoli). |
| chemistry: "Q is pinned at 427 for every viable shape" (`chemistry_families` §1.1, F2) | Q 405 | The proof was stated for **B ∈ [129, 256]** only (F2 §4). b2 = ⌈log2 B⌉ + χ + 1 drops by one at B ≤ 128, taking 15 qubits out of the inner temporary register. The 1-minute `SpinFreeData_V1` grid that F2 itself proposed as proof was never run; it would have shown the step. |
| chemistry: "B 143 → 129 ≈ −1.6%"; "R 8 → 7 ≈ −2%" (§0.4) | (9, 127, 46): Q −5.2% | We moved R and B **one at a time** and priced them in Toffolis (row caches). The move that pays keeps **R·B constant** (same 1,219 rows, same selection width, same loaders) and moves B across a power of two, priced through **Q**. |
| chemistry: "loader products do not scale with rows; a row is worth 1.9 T" (§1.2) | loaders −15.6% a pass | Probably still true, but irrelevant: the loader savings came from **brot** (fewer planes), **table/loader co-design** (wired parity bits enforced during rounding) and **dead-word don't-cares**, none of which change the row count. |
| chemistry: brot 7 "accuracy-killed" (NR-0124 δ7 ≥ 19.3 mHa; Balram 27.1 mHa; CHEM-B11 NULL) | brot 7 at e1 19.9959, e2 19.9992 mHa | The nulls rounded fits that were optimal at 8 bits. A fit with hard sign constraints, Givens-tree labels and **staged rounding with refit** (GPU race) at a larger shape (R 9, C 46) lands under both gates. Our own CHEM-B11 entry said "with codes and QAT": the right idea, killed by evidence from the wrong procedure. |
| chemistry: "realistic ceiling for January −2% to −5%" (§1.6) | −14.3% in one day | Consequence of the three rows above. |
| chemistry: "λ_eff is 0.28% above a hard floor" | 46.918 | **Held** (refits are saturated; public agrees). |
| dqi: "no family reaches the 1.34× target (≤ 2.577e11)" (`dqi_families` §1.3) | best public 3.05e11 | **Held** for families; the gain came in-family. |
| dqi: "GF(2^24) rank-below-99 lottery ≤ −1%" (DQI-B08); "operand-sharing block decompositions 0–0.5%, P 0.10" (DQI-A19) | ThemBones' shared-operand direct-product kernels and converter steps where one product serves two terms (−6% T); tonylica's D16 70/84-term shared formulas (−8% T more) | We priced **single** products (rank of one GF(2^24) multiplication). The public kernels are **multi-output bilinear maps** (one vector times two payloads, one formula reused at four sites with site-specific output masks), whose joint rank is well below the sum of the single ranks. |
| dqi law L3 "every extra word of width doubles the cost" | ThemBones +11 Q for −6% T | L3 is right for words (12+ qubits each); single-qubit trades were left unexplored (one qubit ≈ 0.22% ≈ 780 T at 900 Q). |
| each study: "the leaders' gap is more likely on the other tasks" | the gap was spread: invsqrt −18%, dqi −12%, chemistry −14% | We used the leaders' total gap to argue each task away. Each task had its own in-family lever. |

### 3.2 The common causes

1. **Domain restrictions on nulls.** Each wrong null was correct on its stated domain (B ≥ 129; stored-history
   dividers; single-output products; roundings of an 8-bit fit) and was then quoted without the restriction.
   Rule to add: every structural null in `NULLS.md` carries its domain, and the register entry it kills names
   the domain too.
2. **Families over schedules.** M2.B/M2.C asked for "≥ 5 qualitatively different datapaths". The public gain
   came from the schedule and kernel level of the existing datapaths, which our A-track left to "ports of public
   micro-edits".
3. **No fast exact harness on a generator we control.** ThemBones re-checks every invsqrt variant on all
   3.1M inputs under four HMR modes in 5–8 s, and every dqi kernel as an identity on all inputs. We had the
   evaluation toolkit (S2.1–S2.3) but no generator at the frontier (S2.12 paused), so width moves could only be
   priced on paper.
4. **Pricing in the wrong unit.** At the 39–44 Q level, invsqrt and chemistry moves should have been priced
   first in qubits; we priced them first in Toffolis.

### 3.3 Our ideas against the new base (`Analysis/IDEAS.md`)

**Now IN PUBLIC FRONTIER (from source in the ThemBones line, or exactly described):** INV-A07 (SAT/exact
rebuild programs: ThemBones' digit re-make programs; Almond's 5-AND UNSAT), INV-A12 (angle LP: saturated, 33
rotations proven), INV-A15 (measurement-conditioned repair instead of reverse computation), INV-A24 and
GEN-B04 (qubit-cap ladder: 44 → 39 → 37), GEN-B08 (conditionally clean / lent ancillae: LEND29, q11 lend),
GEN-A08 (functional-dependency mining: the lend relations are exactly that), GEN-B03 (width-optimal
rescheduling: the `game` DP + width step), DQI-A06 (product retention: livescan/LDH), DQI-A09 (proved-zero
interval reassignment), DQI-A13 (cached inverses; INV17 subfield rewrite), DQI-A19 (operand-sharing kernels),
DQI-B03 (projective finish with CRT closing), DQI-B08 (in its multi-output form: D16), CHEM-A02 (parts market:
the chain DP shows −38 is the exchange limit), CHEM-A12 (exact re-planning: hyperplane plans priced by
matroid-parity rank), CHEM-B01 and CHEM-B09 (loader-aware coding: wired bits), CHEM-B11 (brot 7).

**Obsolete (tied to the v12 or 44 Q / Lera bases):** INV-A02, A03, A04, A05, A09, A11, A13, A14, A16, A17,
A18, A19, A21, A22, A23, A25 (stage moves on the v12 datapath; re-derive on the new generator if at all);
DQI-A01, A02, A03, A04, A10, A12, A16 (v12/v13 deltas, now inside the ThemBones kit); CHEM-A01, A03, A07, A09
(Lera/maoxin123 base), CHEM-A06 (the ANGLE circuits now use the full 128-qubit cap), CHEM-B07, B10 (row-axis
shapes at R·B ≠ const: superseded by the (9, 127, 46) shape).

**Still open on the new base:**

- **INV-B10** (parity divider): partly absorbed (digits are re-made, not stored). Worth one cheap screen on the
  ThemBones generator: if some feed-free rows still store or re-make their digits with Toffolis, the parity row
  re-makes them for 0 Toffolis (our S2.18 estimate −35 to −60 T a shot was against a stored-history baseline,
  so the real residue is smaller; re-price it from the generator before building).
- **INV-B01 / GEN-B01** (evolutionary search) now has a concrete search space: ThemBones' per-block programs
  under a wire cap, with the 5–8 s exact check as the judge.
- **DQI-B13** (convert-free Forney): unchanged in principle; convert + unconvert are still ≈ 117.7k of the
  ≈ 328k static Toffolis (Almond's anatomy), so its ceiling is intact; P unchanged (≈ 0.03).
- **DQI-A07/A11/A17** (failure budget, salt-robust selection): more valuable, since `4b572141` fails 2/102,400
  on one foreign salt.
- **DQI-A08** (width/Toffoli Pareto): the 901 → 900 → 890 evidence makes this the main dqi lever (section 4).
- **CHEM-B19** (non-bilinear loaders): still open; the base now has wired bits and 2,400 dead words as
  don't-cares, which helps a synthesis-based loader. Loaders are 51% of T (6,790 of 13,290).
- **CHEM-B21** (C-axis retrain): weakened. At fixed words C 46 → 45 fails twice publicly (21.6 mHa; 20.351 mHa
  plateau); only a joint re-rounding (new words, new loaders) at C 45 remains, worth −0.59% T (−78).
- **CHEM-A10** (helper parking in free word bits) and **CHEM-B14** (word registers as dirty workspace): open on
  the new loaders, small.
- **GEN-A13** (scaled proxies) and **GEN-A10** (local proof per rewrite): the ThemBones harnesses show these are
  the enablers; adopt theirs rather than build ours.

## 4. Where the remaining headroom is

### 4.1 What the public authors say

Asterios Panoras' floors (section 5 of his 31.69052 notebook, on the ThemBones files): the only hard floor is
the 22 input qubits of invsqrt; single-gate deletion is exhausted on all three tasks (none of the 1,493 invsqrt
Toffolis removable by five exact rewrite families; no chemistry loader Toffoli removable over all 1,219
addresses, a deletion costs 4.2–5.6 Ha against a 0.020 limit; dqi wires live 91.7% of the run); λ_eff is
polished out (+0.002 at most); board noise ±0.03 points. His price list at the 31.69 level (points per unit):
**invsqrt qubit +0.27**, invsqrt Toffoli-eq +0.0065, **dqi qubit +0.023**, dqi Toffoli +0.000031, chemistry
Toffoli +0.0008, λ +0.001 per 0.005. Equivalently, +0.1 point needs 0.4 invsqrt qubit, 15 invsqrt Toffolis,
127 chemistry Toffolis, 4 dqi qubits or ≈ 3,300 dqi Toffolis. His verdict: "floors move when someone changes
the design, not when someone deletes a gate"; invest in redesigns that remove structure (a wire, a repeated
block, a table stored twice).

Almond Labs' failed list (section 7 of 31.70105) closes: a second LEND on the same affine relation, dropping a
rotation, a 5-AND zero-scratch block (UNSAT), trace-based root splitting (1.10M T), a ≤ 19-AND GF(64) inverse
and a rank-< 45 GF(16)[t]/t³ product, C 46 → 45 at 7 bits, a 257-T k0. goodpjw2008 adds: −38 is the limit of
loader exchange among 41 public families; C = 46 is the floor for these words; invsqrt (37) and dqi (900) are
full GF(2) rank at their peaks, so no CNOT-only width move exists.

ThemBones' own reading (*Gift* notebook): top-10 teams have much more compute than a single author; the
chemistry table alone took "many GPU-hours".

### 4.2 Our estimates on top of the public set *(estimates)*

| Lever | Size | Basis | Confidence |
|---|---|---|---|
| dqi 900 → 890 Q at the `4b572141` Toffoli level | −2.2% task (+0.25 points) | (890/900)²; Rithvik is at 890 Q with ≈ 328k T, tonylica at 900 Q with 316k | medium: both halves exist publicly, not together |
| invsqrt: tonylica's divider/carry changes (−39 T at 39 Q) on the 37 Q schedule | ≈ −2.4% (37 × 1,608 ≈ 59,500) | additive only if the two lines compose | low–medium |
| invsqrt 37 → 36 Q at +≤ 30 T | −2.0% to −2.7% per qubit | the 39 → 37 step paid +27 T for 2 Q | medium for one qubit, unknown beyond |
| dqi multi-output kernels at more sites (D16-style formulas beyond the four field sites) | −1% to −4% | tonylica replaced four sites; ThemBones' kernels cover the dot products only | low–medium |
| chemistry: brot 6 with the rounding-aware fit | ≈ −7% if accurate (Q −5, Givens −300, loaders −1/7 of planes) | same mechanism as brot 8 → 7 | low (accuracy) |
| chemistry: other power-of-two boundaries in the scorer | unknown; B ≤ 64 at R·B const looks net worse (inner prep +≈ 830 T for Q −14) | hand computation; needs the scorer grid | to check (1 min) |
| chemistry: non-bilinear loaders on dead-word don't-cares (CHEM-B19) | −1% to −3% | S2.23 F3 re-based | low |
| chemistry: joint re-rounding at C 45 | ≤ −0.59% | goodpjw2008's break-even | low |

### 4.3 What the 33.35 leader may have

33.348 is +0.74% over 33.10, i.e. a −2.2% cost product. **The same public set with the dqi at 890 instead of
900 qubits gives 33.1019 × (900/890)^(2/3) = 33.349**, within 0.001 of the leader. Other single moves that fit:
invsqrt at ≈ 36.1 Q-equivalent (−2.2%), or chemistry at T ≈ 13,000. A 36-qubit invsqrt alone would give
33.40 (too high), so if it is one move it is most likely the dqi width; it could also be several smaller moves.
This is inference from one number, not evidence.

## 5. Recommendations for the plan

### 5.1 Base per task

| Task | Board file now (credited) | Source base to adopt | Rebuild cost | Then rebuild from source |
|---|---|---|---|---|
| invsqrt | `d2148a7c` (Rithvik, 60,955) | **ThemBones generator** (`src_s18_srcM` + patch + `f9_i1_stack.sh`; byte-for-byte record) | ≈ 20 min, 1 core, Python 3.13 + numba | Almond's LEND29 + T-fusion as passes (63,016); tonylica's v50 divider and v220 carry if they port; then the width work toward 37 Q |
| dqi | `4b572141` (tonylica, 3.050e11) | **ThemBones kit + chain** (`kv_unfused5`, 29 modules, oracle) | ≈ 40 min, 2 cores; file A never rebuilt by them, so our rebuild is itself a result | Almond's five edits; then tonylica's D16 formula tables (data + `arithmetic.py`) into our chain; then width to 890 |
| chemistry | 13,290 chain + `fc4b4c01` (252.53M) | **ThemBones chem tools + table + `.kmx`**, with tonylica's `build.py` route (KMX + TSV) | seconds for loaders; no GPU needed while the table is carried | our own dead-word criterion and chain DP; disassemble the rupam2020 / Rithvik / boristown loaders to KMX with credit; then CHEM-B19 |

Our v12 generators (invsqrt 46 Q, dqi v12) are retired as bases; their tools (exact integer model, proof
harness, `tools/kmb.py`, kickmix-engine gates) carry over. Environment: the ThemBones and tonylica code
targets **CPython 3.13** (tonylica bundles a cp313 `kickmix.so`; ThemBones uses numba); our WSL venv is 3.14,
so add a second venv `~/venvs/gqai313` with numpy, scipy, numba (and JAX only if we refit the table).

### 5.2 Rest of Phase 2 (21–30 Nov)

| Step | Content | Gate / go-no-go | Effort |
|---|---|---|---|
| S2.25 | Adopt the 33.10 public set: pin the 18 members by SHA-256 (`receipt.json` of goodpjw2008), credits and NOTICEs (Syed Asad Ali first), our gates in WSL (invsqrt full 3 min; dqi 16 salts via kickmix-engine; chemistry 45 s), submit | PASS on all gates; board ≥ 33.10 | 0.5 day |
| S2.26 | Reproduce the ThemBones line in a 3.13 venv: invsqrt `c639a045`, dqi file A `2c043a9a`, chemistry loaders and npz `e05fb18f`; record run time and memory | byte-for-byte SHA match on each (dqi: first independent rebuild) | 1 day |
| S2.27 | invsqrt from source to the public 39 Q level: LEND29, T-fusion, then tonylica's divider/carry ideas as generator passes | ≤ 63,016 (then ≤ 61,696) from source, exact proof on all inputs | 1 day |
| S2.28 | invsqrt width 39 → 38 → 37 from source: liveness diff of `d2148a7c` vs `c639a045` (static, ours), ThemBones' `game` DP and width step at cap 38/37, INV-B10 re-screen as one program family | 38 Q at C ≤ 1,604 or 37 Q at C ≤ 1,647 (both beat 60,955); else keep the public binary on the board | 2–3 days |
| S2.29 | dqi from source: Almond's five edits, then the D16 formulas on our rebuilt chain | ≤ 3.3248e11 (then ≤ 3.05e11) on 16 salts, no more failures than the parent | 1.5 days |
| S2.30 | dqi width 900 → 890 at the 316k level (DQI-A08/A09): proved-zero reassignment, inverse-helper re-emission at the non-multiply peaks | each qubit −0.22%; go if ≥ 4 qubits by 30 Nov | 1 day |
| S2.31 | chemistry from source: dead-word criterion, our chain DP, table exchange; then CHEM-A10 on the new loaders | 13,290 at λ ≤ 46.9182 from our repo | 1 day |
| S2.32 | Freeze (M2.D): `STATUS.md`, build commands per task, Phase 2 review | a fresh session rebuilds and resubmits | 0.5 day |

### 5.3 January (Phase 3)

1. **M3.1 catch-up** (1–2 days): pull the board and new notebooks (GEN-A06), re-run S2.25 on whatever moved.
2. **M3.2 invsqrt ≤ 36 Q** on the ThemBones generator: evolutionary search (GEN-B01) over per-block programs
   under a wire cap, judged by the 5–8 s exact all-input check; target −2% to −5%.
3. **M3.3 dqi multi-output kernels and width**: search shared-operand bilinear formulas (flip graphs with
   plus transitions, Kauers–Moosbauer / Arai et al., as ThemBones credit) for every field site, not only the
   four D16 sites; keep 890 Q; DQI-B13 screen once (1 day) since its ceiling is unchanged.
4. **M3.4 chemistry fit on the GPU** (Kaggle P100/T4 or the RTX 4060): first the 1-minute scorer grid over
   (R, B, C, brot, J) to list every power-of-two step; then a rounding-aware, loader-aware fit (wired bits as
   constraints, dead words as don't-cares) at brot 7 / C 45 and a brot 6 attempt; go if e1, e2 < 20 mHa with
   λ ≤ 46.95.
5. **M3.5 CHEM-B19** non-bilinear loaders with don't-cares, CPU.
6. **M3.6 team decision** before 13 Jan, as planned; the prize-eligibility question (recovered assembly as
   "source") asked on the forum early in January.

## 6. Proposed register updates (not applied)

### 6.1 `Analysis/IDEAS.md` statuses

| Entry | Proposed status |
|---|---|
| INV-A07, INV-A12, INV-A15, INV-A24, GEN-B03, GEN-B04, GEN-B08, GEN-A08 | IN PUBLIC FRONTIER (ThemBones 39 Q generator, Almond LEND29; Rithvik 37 Q compiled only) |
| INV-A02, A03, A04, A05, A09, A11, A13, A14, A16, A17, A18, A19, A21, A22, A23, A25 | OBSOLETE (v12 / 44 Q base retired, S2.25) |
| INV-B10 | OPEN, re-price on the ThemBones generator (partly absorbed: digits are re-made, not stored) |
| INV-B01, GEN-B01 | OPEN, raised: search space = ThemBones per-block programs under a wire cap |
| DQI-A06, A09, A13, A19, B03, B08 | IN PUBLIC FRONTIER (ThemBones kit; Almond livescan/LDH; tonylica D16) |
| DQI-A01, A02, A03, A04, A10, A12, A16 | OBSOLETE (inside the ThemBones kit) |
| DQI-A08 | OPEN, raised: main dqi lever (900 → 890 ≈ −2.2%) |
| DQI-A07, A11, A17 | OPEN, raised (`4b572141` fails 2/102,400 on a foreign salt) |
| DQI-B13 | PARKED, unchanged ceiling |
| CHEM-A02, A12, B01, B09, B11 | IN PUBLIC FRONTIER (dead-word chain DP; hyperplane plans; wired loaders; brot 7) |
| CHEM-A01, A03, A06, A07, A09, B07, B10 | OBSOLETE (Lera / maoxin123 base and row-axis shapes superseded) |
| CHEM-B19 | OPEN, re-based on wired loaders + 2,400 dead words |
| CHEM-B21 | OPEN, weakened (fixed-word C 45 fails twice publicly; joint re-rounding only, ≤ −0.59%) |
| new INV-A26 | Port tonylica's v50 divider rewrite and v220 carry schedule onto the ThemBones generator (−39 T at 39 Q) |
| new INV-A27 | 37 Q from source: reproduce the 39 → 37 schedule (static diff of `d2148a7c`, cap-37 `game` DP) |
| new DQI-A20 | D16 multi-output formulas on our rebuilt chain, then at every field site |
| new CHEM-B23 | Scorer grid over (R, B, C, brot, J) at constant R·B: every power-of-two step of b2, ⌈log2 R⌉, uniform terms |
| new CHEM-B24 | brot 6 with a rounding-aware, loader-aware fit (wired bits, dead words), GPU |
| new GEN-A16 | Second WSL venv on CPython 3.13 (numba) for the public generators |

### 6.2 `experiments/NULLS.md`

**Corrections (new rows that supersede ours; old rows stay):**

- Correct the invsqrt width floor (`invsqrt_families` §2.3 floor 3, "43–44 owners at the endpoint"): broken by
  measured-out digits re-made in the undo and lent wires (ThemBones 39 Q, Rithvik 37 Q). Domain of the old
  claim: stored-history dividers.
- Correct "Q pinned at 427" (`chemistry_families` F2): true only for B ∈ [129, 256]; B ≤ 128 gives inner
  temporary 287 and Q 405 at (9, 127, 46), brot 7.
- Correct NR-0124 / CHEM-B11 (brot 7 accuracy null): true for roundings of 8-bit-optimal fits; a staged
  rounding-with-refit at (9, 127, 46) passes at e1 19.9959 / e2 19.9992 mHa.
- Correct the DQI-B08 / DQI-A19 pricing (≤ −1%, ≤ −0.5%): multi-output shared-operand kernels gave −6% (ThemBones)
  and −8% (tonylica D16).

**New public nulls (with their domains):**

| Task | Null | Number | Source |
|---|---|---|---|
| invsqrt | second −0.5 T lend on the same affine relation as LEND29 (H1 row 7) | 1,417,293 failing shots | Almond Labs |
| invsqrt | drop one of 33 rotations, or snap rotation 13 to the T grid (family-C datapath) | LP infeasible on all 33 | Almond Labs |
| invsqrt | 5-AND zero-scratch rebuild of q18 (6 ANDs today) | UNSAT on the staged templates | Almond Labs |
| invsqrt | single-gate deletion on `c639a045` (five exact rewrite families) | 0 of 1,493 removable | Asterios Panoras |
| invsqrt | minimal deletion of the 39th wire of `c639a045` | +449 T vs a 42.5 T break-even | Asterios Panoras |
| invsqrt | CNOT-only width reduction on `d2148a7c` | full GF(2) rank 37 at 9,763 peak steps | goodpjw2008 |
| dqi | trace-based root splitting instead of the Chien scan | 1.10M T vs 159.8k | Almond Labs |
| dqi | ≤ 19-AND GF(64) inverse; rank-< 45 GF(16)[t]/t³ product (flip-graph neighbourhoods) | 20 ANDs, rank 45 hold | Almond Labs |
| dqi | "dead tail" removals on the 901 Q file | 12/18 and 20/14 value errors | Asterios Panoras |
| dqi | CNOT-only width reduction on `4b572141` | full rank 900 at 3,208 of 3,925 peaks | goodpjw2008 |
| dqi | tickets (INFO_TAG re-rolls) as a cost lever | ±0.03 points, board draw uses the secret salt | Asterios Panoras |
| chemistry | C 46 → 45 at 7 bits, fixed words | e2 21.6 mHa, λ > 48.2; discrete plateau 20.351 mHa | goodpjw2008, Almond Labs |
| chemistry | refit of (W, H2_I) on the 13,290 words | no measurable gain (three active limits) | goodpjw2008 |
| chemistry | loader exchange across 41 public families | −38 is the limit (13,328 → 13,290) | goodpjw2008 |
| chemistry | loader Toffolis writing only dead words; Toffolis replaceable by CNOTs | none of 1,856 | goodpjw2008 |
| chemistry | cheaper public loaders at positions 6 and 14 that change live words | e2 58–1,448 mHa | goodpjw2008 |
| chemistry | 257-T k0 by a fake-carrier dead-cell rewrite | invalid on the reverse (involution) pass | Almond Labs |

### 6.3 `experiments/LEDGER.md`

One row for this research step: "public frontier audit 21 Nov; verdict: adopt the 33.10 set, move every base to
the ThemBones source line; four of our 19 Nov nulls corrected with their domains".

