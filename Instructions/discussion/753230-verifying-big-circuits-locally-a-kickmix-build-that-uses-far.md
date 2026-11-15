# Verifying big circuits locally: a kickmix build that uses far less memory (same results)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753230 — copied 2026-11-15 (raw Markdown, unedited).

**manwithacat** · 2026-10-27 16:56 UTC · 2 votes

If you verify circuits locally, especially large DQI circuits, this may save you memory and time.

All three official verifiers call `circuit.reaction_depth()` after the shots have run. In upstream kickmix that
call first copies every operation of the circuit into a list, about 135–170 bytes per op, just to walk it once.
On a 55-million-op DQI circuit that is a 4–7 GiB spike, and the whole verification peaks at about 8–9 GiB whatever
`max_workers` you choose. The cost grows with the circuit: by my estimate `reaction_depth()` alone passes 16 GiB at
roughly 100 M ops.

I maintain an openly AI-assisted fork of kickmix that computes the same depth while streaming over the circuit, so
the spike is gone: https://github.com/manwithacat/kickmix-macos

Measured on my Mac (all processes of the verification together; same 55 M-op DQI circuit):

| DQI verification | upstream | fork |
| --- | --- | --- |
| 1 worker | 7.8 GiB | 2.0 GiB |
| 4 workers | 7.8 GiB | 4.4 GiB |
| 6 workers | 9.0 GiB | 5.6 GiB |

With the fork, peak memory is roughly 1.9 GiB plus 0.6 GiB per worker, so you can pick the worker count for your
machine. invsqrt and chemistry use little memory either way (under 1.5 GiB).

- **Same results.** The official verifiers give identical summaries, field for field, with and without the
  changes, on a full three-task portfolio. CI checks simulation parity against unmodified upstream on every run, on
  macOS, Linux and Windows.
- **Also faster.** A handful of profiling-driven changes (cheaper per-shot operation counting, word-at-a-time
  `write/read_across_shots`, phase readout without a per-call `Fraction` GCD) take full local verification of all three
  tasks on my Mac from about 115 s to about 77 s.
- **Install:** `pip install "git+https://github.com/manwithacat/kickmix-macos.git#subdirectory=kickmix"`
  (macOS, Linux x86-64/arm64, Windows).

Two tips that apply whichever build you use:
1. The remaining ≈ 1.9 GiB in the parent is mostly the verifier itself: it serialises the whole circuit into memory to
   hash it. That is verifier code, so the fork doesn't change it.
2. The official verifiers parallelise with `fork`, which Windows doesn't have, so on Windows they silently run in one
   process. Either use WSL2, or the repo's `tools/verify_parallel.py`, which runs the unchanged verifiers' own seeded
   batches on `spawn` workers and gives the same report (works with upstream kickmix too).

Upstream kickmix is what Kaggle runs and the official, reviewed source; please treat the fork as a local
convenience and check its claims with the repo's `tools/verify_parity.sh`. The code was written by AI agents at my
direction and is disclosed as such; it is not offered upstream under Qualtran's current policy. Thanks to the
Qualtran team for kickmix and the verifiers.


---

## Comments (2 threads)

Replies are indented with `>`, one level per reply depth.

**manwithacat** — **Topic author** · 2026-10-27 20:15 UTC · 0 votes

A note on how the memory figures above were measured, since it matters for reproducing them. They are the peak of the **summed physical footprint** of every process in the verification (macOS `phys_footprint`; on Linux the equivalent is the summed **PSS**). The DQI verifier forks its workers after loading a ~1.5 GB circuit, so the workers share those pages with the parent; summing **RSS** counts the shared circuit once per worker and can make both builds look about the same, while `/usr/bin/time`'s max RSS sees only the parent. Measured the same way on a Kaggle CPU notebook (Linux, summed PSS, organisers' sample submission, default workers), DQI peaks at 8.95 GB with the shipped kickmix.so and 3.35 GB with the fork, with every result field identical across all three tasks. On Kaggle the fork was also 0.82× the wall time on DQI and 0.87× on chemistry, but about the same on invsqrt, where per-shot Python work dominates. A Kaggle notebook that runs this comparison on Kaggle's own hardware is on its way; I'll link it here.

> **manwithacat** — **Topic author** · 2026-10-27 20:44 UTC · 0 votes
>
> The notebook is up: [kickmix fork on Kaggle: same results, less memory](https://www.kaggle.com/code/manwithacat/kickmix-fork-on-kaggle-same-results-less-memory). It installs the fork from a prebuilt wheel (hash-checked, a couple of seconds) into its own directory, then runs the unchanged official verifiers with the shipped kickmix.so and with the fork, one after the other, and compares every result field, wall time and peak summed PSS. Set `MODE` and `SUBMISSION_DIR` to point it at your own circuits. The same helper is also a utility script you can attach and import: [kickmix_tools](https://www.kaggle.com/code/manwithacat/kickmix-tools). Latest run on the organisers' sample submission: all fields identical; DQI 9.07 GB → 3.20 GB peak.

**Sohier Dane** — **Kaggle staff** · 2026-10-27 18:06 UTC · 0 votes

Nice! Can I suggest filing a PR in the Qualtran repo?

> **manwithacat** — **Topic author** · 2026-10-27 18:14 UTC · 0 votes
>
> i'm afraid they're not accepting AI assisted code. I can understand their reasons. If you know someone at a high-level within the Qualtran team who would be willing to do the code review, I've structured my repo as a collection of ongoing patches that they could cherry pick. All my code is open source, but I'll happily take some google merch if it's offered.

> > **Tanuj Khattar** — **Host** · 2026-10-27 20:27 UTC · 1 votes
> >
> > We are happy to accept AI-generated PRs in Qualtran. Please follow the general best practices of sending small PRs with human readable short descriptions so it's easy to review.

> > **manwithacat** — **Topic author** · 2026-10-27 20:35 UTC · 0 votes
> >
> > Thanks, that's great to hear. One thing I'd like to get right first: Qualtran's CONTRIBUTING.md points to CLA clauses 5 and 7 and says AI-generated code doesn't qualify as the contributor's original creation, so I'd be uncomfortable signing the CLA for these patches as things stand. Is there a way you'd like that handled (e.g. a disclosure in each PR description, or an updated note in CONTRIBUTING.md)? If so, I'll send the reaction_depth streaming change first as a small PR against issue #1993, with a short human-readable description and the AI assistance stated up front.
