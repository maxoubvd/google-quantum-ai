# Final evaluation scope: published verifier or additional correctness checks?

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755427 — copied 2026-11-15 (raw Markdown, unedited).

**Pizzaboi** · 2026-11-03 15:59 UTC · 15 votes

Could the organizers clarify the correctness requirements and final evaluation scope before the submission deadline?

## Why this matters

For Task III, SCORING_README.md describes ANGLE_k as a reversible XOR oracle:

|a>|z> → |a>|z XOR delta[k,a]>

However, the supplied verifier checks the word states reached by the forward/reverse traversal starting from zero. The README also explicitly acknowledges that passing does not prove correctness on every target-register input.

This distinction affects which optimizations are permitted. A loader can be correct on the reachable encoded subspace while behaving differently from arbitrary-target XOR outside that subspace.

## Questions for the organizers

1. Will final/private evaluation use the published verifier and the same input-preparation scheme, with unseen salts? Or should we expect additional correctness checks or manual review after the deadline that enforce a broader input domain?
2. Are optimizations relying on invariants of the reachable word states permitted, provided outputs, phase and workspace cleanup are correct on those states? Or must every component implement arbitrary-target XOR on inputs the published verifier never supplies?
3. If broader requirements will be enforced, could their scope be clarified before the deadline?

## Reproducible example

[Notebook with a minimal reproducer and complete submission archive
](https://www.kaggle.com/code/asimandia/task-iii-xor-contract-minimal-reproducer)

The notebook modifies the official baseline by adding only three CX gates to one angle loader. It passes the full unchanged supplied Task III verifier on Kaggle at the reference cost, but differs from arbitrary-target XOR on independent word inputs.

| Check | Result |
|---|---|
| Full supplied Task III verifier | PASS — 254,464 shots per run |
| Task III cost | 3,355,599,588.3155375 — equal to the baseline |
| Original loader on additional challenge inputs | 0 incorrect outputs out of 15,904 |
| Modified loader on the same inputs | 7,952 outputs differ from arbitrary-target XOR |

The full supplied-verifier runs are distinct from leaderboard evaluation.

This illustrates the distinction; it is not an allegation about other submissions or a claim that the actual algorithm fails on its reachable inputs. I am asking about the intended correctness domain and evaluation policy, not for hidden salts or individual test cases.
