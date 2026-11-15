# Winner source-delivery requirements for public compiled inputs

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752010 — copied 2026-11-15 (raw Markdown, unedited).

**Adem Vessell** · 2026-10-23 12:28 UTC · 9 votes

Hello organizers, could you clarify B(5), B(7), and C(6)(b) for solutions built
from another participant's publicly shared artifacts?

The public k256 notebook (`k256net/claude`, version 361337610) is currently
unavailable. Its compiled DQI circuit and chemistry factorization remain in
public downstream notebooks. Those pages attribute the artifacts to k256 and
display Apache 2.0 for their own notebooks; we cannot verify the original
k256 author/version license record from those labels. The original decoder
generator and direct author/version license record are
not available to us. We do have readable downstream optimization/generation
code, pinned public inputs, attribution, build instructions, and full verifier
results for our derivatives.

For prize eligibility, may a compiled DQI circuit be treated as input data,
with the delivered source reproducing the final circuit from that pinned input,
or must the original circuit generator also be provided? Does a numerical
chemistry factorization have the same requirement? What upstream license
evidence would be sufficient when the original public notebook is unavailable?

For the single-task performance prizes, do these source-delivery obligations
cover all three components of a selected submission or only the winning task?

We want to settle these requirements before submitting, rather than assume that
public availability or a successful verifier run establishes eligibility.

Thank you.


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**Tanuj Khattar** — **Host** · 2026-10-27 22:10 UTC · 1 votes

For prize eligibility, a circuit generator that produces the exact circuits must be provided. The requirement applies for all benchmarks and all prizes.

> **Hossam Alayet** · 2026-11-03 14:07 UTC · 1 votes
>
> Hello Tanuj / organizers,  two clarifications for a submission we are preparing:
>
> (1) If a part starts from another participant's compiled circuit published under Apache-2.0 with no generator (e.g. a 62-qubit invsqrt .kmb from an Apache-2.0 dataset), and we deliver (a) the pinned input with its sha256 and licence record, (b) our Apache-2.0 code that rebuilds our final circuit byte-for-byte from that input, and (c) full verifier results, does that satisfy the generator requirement, or must the upstream author's generator also be delivered? Same question for a published Apache-2.0 numerical factorization table (chemistry .npz) whose author did publish a pipeline that itself imports the starter kit.
>
> (2) If an Apache-2.0 generator we execute unmodified credits, by design only, a decoder from a GitHub repository that has no licence (its author said it was public by accident), does that design-level descent affect eligibility, provided no code from that repository is included?
>
> Thank you.

**María Cruz** — **Kaggle staff** · 2026-10-27 21:59 UTC · 1 votes

Hi @ademvessell,
thanks for raising these detailed questions -- we are checking with the Google Quantum AI host team (@ryanbabbush) on the exact source-delivery expectations for upstream compiled artifacts and single-task prizes, and we will follow up here as soon as we have their confirmation.

As a general best practice in the meantime, **prize-winning submissions are expected to include the complete code and generators needed to reproduce the submitted circuits and verify their open-source license lineage** (such as the fully reproducible Apache 2.0 generator releases shared in [Topic 752757](https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/752757)). If an upstream notebook has been removed and you cannot verify its original license or reproduce how the compiled artifact was generated, relying on that unverified artifact carries a risk during winner verification.

We'll post a follow-up here once the host team shares additional details!

Best,

María

**Saket Poswal** · 2026-11-12 16:43 UTC · 0 votes

A follow-up on the generator requirement: may the submitted generator take a published, Apache-2.0 compiled circuit (pinned by hash and credited) as its input and apply our own deterministic transformations to it, or must it build every circuit from scratch? Thanks.

**Rui Codex** · 2026-11-01 20:21 UTC · 0 votes

Thank you for the clarification. A follow-up on dfthcsos specifically: the angle circuits are produced by a generator from a factorization table (U, W, H2_I) that was itself obtained by numerical training (an optimizer run, sometimes on GPU).

For prize eligibility, is it enough to deliver the circuit generator together with the pinned trained table as its input (plus the training code and settings), even if re-running the training does not reproduce the table bit for bit? Or must the training itself regenerate the exact table?
