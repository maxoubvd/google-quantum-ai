# Three-task circuit compilation: methods, results and AI assistance

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759333 — copied 2026-11-15 (raw Markdown, unedited).

**levent erkul** · 2026-11-15 03:45 UTC · 0 votes

Methods and results report, 15 November 2026. This post describes the submitted
circuits and our AI-assisted workflow. It is not a released executable
reproduction bundle.

## Summary

We investigated the three circuit-compilation tasks using published research,
public participant work, bounded experiments and the supplied verifiers.
Our selected three-task submission completed Kaggle evaluation with displayed
public score **6.24934** on 15 November 2026.
For inverse-square-root phasing, we combined direct quadratic approximation with
exact lookup and arithmetic/workspace improvements.
For DQI, we adopted and independently source-built the credited SomeJoint P04
V21 decoder without algorithmic changes.
For chemistry, weight compensation and discrete angle polishing made an 11-bit
representation acceptable at the original factor dimensions and Givens tree.
The work was substantially AI-assisted under the participant's direction, with
claims tied to saved circuit identities and verification results.

## Selected results

| Task | Qubits | Recorded complete task cost | Reference / candidate |
| --- | ---: | ---: | ---: |
| invsqrt | 154 | 544,507.1702641873 | 3.81404 |
| DQI, adopted P04 V21 | 889 | 412,317,672,621.7634 | 59.2734 |
| Chemistry | 733 | 3,108,228,989.5773115 | 1.07959 |

These task costs come from the identified private verifier runs. The reference
comparisons use preserved official-reference records: chemistry's published
worked example and the earlier invsqrt/DQI baseline records. Different circuit
hashes and salts produce different sampled trajectories, especially relevant to
DQI's error-dependent repetition factor. The ratios are not matched-cohort
ablations or claims to beat the public frontier. **6.24934 is the combined
official public score**, not a separate score for each task or a final private rank.

## Inverse-square-root phasing

The circuit phases a 22-bit input by the prescribed inverse-square-root value,
then restores the input and scratch. Instead of a cubic seed followed by Newton
refinement, our selected construction directly evaluates a piecewise quadratic
with two reversible fixed-point Horner products. Shared algebraic-normal-form
coefficient loading supplies the exact finite table.

Four binary normalization flags replace ten unary flags. An adapted
borrowed-padding carry routine reduces peak workspace; a proven constant sign
simplifies repeated arithmetic rows. Finally, assigning suitable values to
unreachable lookup addresses removes nonlinear products without changing any
reachable coefficient. The complete cost is **73.7811% below** the recorded
reference, including rotation-synthesis charges, with 3,284 Toffoli-like
operations and 154 qubits.

The isolated adder adaptation comes from Alex Chilton's pinned source, with
credit to Rudra Bhaskar, anazemcev/Artem1410 and Sumit Chongder. We preserve that
lineage and its Apache notice. Normalization, interpolation, ANF and borrowed
carry are established methods; our contribution is the particular implementation,
integration and validated simplifications. Two private runs each covered all
3,133,452 prescribed inputs at 64 shots per input: 401,081,856 shots combined.

## Decoded quantum interferometry

We selected **SomeJoint/Agniv Banerjee's P04 V21** reversible decoder for the
RS(4095,4025), GF(2^12) task. Its low-width Schur/continued-fraction pipeline
extracts and synthesizes the decoding polynomials, performs the required phase
scan, then reverses computation with measurement-phase repayment. We rebuilt
the authenticated source unchanged; the resulting circuit is not our original
decoder algorithm.

Our 102,400-shot private qualification observed zero value or phase errors,
mean 437,947.14455078123 Toffoli-like operations, 889 qubits and repetition factor
1.1912604648115832. This is one salt, not an all-input or repeated-salt guarantee;
the upstream author reports rare phase failures on other samples. Our earlier
decoder modifications and arithmetic studies are retained separately and are
not added as further savings to P04.

Credits include SomeJoint, GPT-6 Astra Pro as identified upstream, Ahmed Berat
Özer, Alex Chilton, James Day and the wider inherited chain in the detailed
report. P04 contains an Apache license, but the lineage also refers to an earlier
repository without an established explicit license. James Day publicly expressed
permission for use and conditionally discussed a future Apache release. We
disclose that source-delivery scope rather than presenting it as settled; the
competition hosts determine eligibility.

## Chemistry: DFTHC-SOS angle loading

We retained the supplied FeMoCo76o113e factor dimensions `(65,60,60)`, the
75-edge Givens tree and fixed offset/symmetry arrays. Direct 11-bit rounding
initially failed the 20 mHa two-body error gate. Progressive rounding with
continuous compensation, followed by weight-only fitting and bounded discrete
angle changes, reduced the residual until a fully decoded checkpoint passed.

The selected checkpoint has one-body error **6.777859 mHa** and two-body error
**19.924853 mHa**. It uses ordinary loaders with nine 11-bit angles per batch:
nine circuits instead of the reference's ten. Complete resources are 84,989
Toffoli-like operations and 733 qubits, with maximum loader width 122 below the
128-qubit cap. Complete cost is **7.3719% below** the published reference.
The selected two-body margin is only 0.075147 mHa.

The approach was informed by public compensation/polishing descriptions,
including Marculera and Michael Rustad. Our candidate retains the original
large factor shape and ordinary generator, rather than adopting their compact
factorizations or specialized loaders. Local and private Kaggle qualification
each passed 254,464 shots on distinct salts, recovered the intended words exactly
and agreed on complete cost. Dense-QROM, shared-product and exact local-synthesis
studies did not supply a competitive full candidate; their scoped outcomes remain
in the detailed record.

## AI assistance and human direction

OpenAI Codex and GPT-based coding agents substantially assisted literature and
source inspection, experiment design, implementation, optimization, analysis,
testing and drafting. Agent roles included coordination, bounded implementation
and review where completed; automated checks and saved verifier receipts support
the reported results. We do not claim that every proposed review completed or
that the human participant manually audited every gate.

The human participant set the research objectives, directed follow-up questions,
approved experimental progression and made the publication/submission decisions.
This disclosure concerns our workflow; it does not transfer ownership of adopted
participant algorithms or make AI-generated proposals verified results.

## Reproduction, evidence and sources

The submitted ZIP contains `invsqrt.kmb`, `dqi.kmb`, nine chemistry KMB files and
the chemistry `submission.npz`. Its SHA-256 is
`8379394b65b4e129b7bb6cbb28a631b2250ce3a766995aaf82d3cb390b31d49f`.
Our retained manifest pins all 12 selected files.
Reproduction requires the identified source/checkpoints, competition inputs,
dependency versions and CPU environment. The hosted checks used the supplied
Kickmix binary; a newer local SDK is not assumed equivalent. Private notebooks
and local evidence links are not yet publicly accessible reproduction resources.

Detailed methods, commands, runtime/seed records and negative experiments are
retained in our local reports. Those reports, source files and private validation
notebooks are not attached to this post. Accordingly, the post alone is not a
self-contained reproduction release. The adopted participant sources are linked
below; a future source/evidence release must preserve their original notices.

Primary foundations include the [competition companion paper](https://arxiv.org/abs/2610.23832),
[Qualtran/Kickmix](https://github.com/quantumlib/Qualtran),
[DQI](https://www.nature.com/articles/s41586-025-09527-5), and
[spectral-amplification electronic-structure simulation](https://journals.aps.org/prx/abstract/10.1103/pb2g-j9cw).
Implementation provenance includes [Alex Chilton's pinned source](https://www.kaggle.com/code/alexchilton/24-94-lb-verified-three-task-portfolio?scriptVersionId=361918908),
[SomeJoint P04 V21](https://www.kaggle.com/code/somejoint/quantum-circuits-24-65-local-gm-chem-refit?scriptVersionId=363462588),
[Marculera's pipeline](https://www.kaggle.com/datasets/marculera/dfthc-merge-pipeline),
[Rustad's refitting work](https://www.kaggle.com/datasets/michaelrustad/vqcc-dfthc-r10-merged-groups), and
[James Day's permission statement](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753660).
The linked participant publications provide the inherited credit chains. Their
reported results are distinct from the private verification results we measured
on the selected files.
