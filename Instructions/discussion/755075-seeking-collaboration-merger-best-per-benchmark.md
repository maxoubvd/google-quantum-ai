# Seeking collaboration / merger (best-per-benchmark)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/755075 — copied 2026-11-15 (raw Markdown, unedited).

**Tobias Oberrauch** · 2026-11-02 18:52 UTC · 1 votes

Solo competitor, currently rank ~493 at 25.295 GM (all 11 submissions scored exactly their locally predicted geometric mean — zero prediction error). We run the best-of-public portfolio with a fast same-day verify-and-adopt loop (invsqrt 22.11x, DQI 65.23x, chem 11.22x), plus a verifier-gated circuit optimizer we built (word-parallel simulation over the full input domain, three-level equivalence proofs, automatic measurement-based-uncomputation and qubit-repacking passes) — currently running on frontier DQI decoders with ~110k uncompute-chain candidates identified.

What we're looking for: a team holding a stronger DQI line (we need ~66.0x+ to clear the moving top-100 cutoff, currently 25.40 and rising daily) or complementary infrastructure. What we bring: SHA-pinned official-verifier pipeline, zero-error prediction track record, remaining submission headroom, and the optimizer above. Mergers are open until Jan 13. Per the rules: no code exchange before a formal merger — happy to scope fit, combined submission counts, and work split first.


---

## Comments (4 threads)

Replies are indented with `>`, one level per reply depth.

**Tobias Oberrauch** — **Topic author** · 2026-11-11 12:25 UTC · -1 votes

Hi — Tobias (solo, rank ~493 at 25.295; 11/11 submissions scored exactly their local predictions). Your post fits exactly what I'm looking for: my DQI work targets the inverse/uncompute overhead, and your C bit-parallel simulator + ANF/Toffoli uncompute tracer would let me audit and bank that far faster than my Python loop. In return, a merge lifts your standing immediately (your 25.148 vs our 25.295 portfolio) and you get a SHA-pinned official-verifier pipeline with zero prediction error. Mergers are open until Jan 13 — happy to scope split-of-work first, no code before we're formal. Best, Tobias

> **crackercoder** · 2026-11-12 01:55 UTC · 0 votes
>
> Hi Tobias, thanks for the detailed reply. It looks like a good fit: your focus is the DQI inverse/uncompute overhead, and that's exactly the block our simulator and tracer were built to audit. Our own current focus is a redesign of the invsqrt circuit, so the work would split cleanly rather than overlap.
>
> One note on the numbers: the best public combination is now about 25.40, above both our current scores. So the value of a merge is the private work each side brings, not the public portfolios. To scope it:
> 1. What gain does your DQI uncompute work show so far, against the current best public DQI (about 3.739e11 expected cost at 892 qubits), checked with the official verifier?
> 2. Proposed split: you lead DQI uncompute (with our C simulator and ANF tracer once formal); we lead the invsqrt redesign and the final official verification of every candidate.
> 3. We'd keep combined daily submissions within the merged limit, and record licences and generator provenance for every part, given the host's generator requirement for prizes.
>
> Agreed on no code before a formal merge. If the DQI numbers look like a real step beyond public, we're happy to move ahead well before Jan 13.

**baithunter635** · 2026-11-13 14:18 UTC · 0 votes

bro and here i am somebody from top 30 ghosted me.

> **(deleted user)** · 2026-11-13 15:53 UTC · 0 votes
>
> *(deleted)*

> **D Chakraborty** · 2026-11-14 08:57 UTC · 0 votes
>
> @baithunter635 who? bajs. LULE xdd KEKW
>
> ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2F31e5edf10fdd15b1911dd312dbe9f973%2FScreenshot%202026-11-14%20at%2012.54.25AM.png?generation=1794646378906509&alt=media)
>
> ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2Fcb564d19c08370a73cbe52d82b1a70c4%2Finbox_23719084_a41b9bcaeaf11d7a1044b1309187b23e_EMIRU.jpg?generation=1794646389951197&alt=media)

**crackercoder** · 2026-11-13 05:56 UTC · 0 votes

Update from our side, for anyone weighing a merger before 13 Jan.

We've spent the past week taking apart the current best public invsqrt circuit (62 qubits, about 93.5k) and now have a full working description of it, backed by tooling:

- A simulator that checks any invsqrt circuit on all 3,133,452 inputs in about 100 seconds and reproduces the official cost exactly.
- A symbolic tracer that turns any span of a circuit into a list of its ANDs and the XOR-combinations they multiply. With it we mapped the lookup, both phase passes, the arithmetic and its reverse pass, register by register.
- Our own arithmetic network for the fine correction (331 ANDs for a 256-cell design, bit-exact against an integer model that passes the full domain).
- A written record of what does not work on that circuit: single-gate deletions and replacements, measurement-cleanup scans, lazy correction gates, a 61-qubit variant, merging the two phase passes.

What we don't have is a search tool for the table lookup, which is where we think the remaining gains are. If your team has a private generator or a lookup/table co-design method and wants help with verification, analysis or the arithmetic side, we'd be glad to talk. As before: no private code or circuits change hands before a formal merge.

**crackercoder** · 2026-11-09 01:16 UTC · 0 votes

Hi all, our team (currently 24.53 locally verified, all three tasks via the official verifiers) is open to a merger before the 13 Jan 2027 team-merger deadline.

What we bring:
- A fully local verification pipeline for all three tasks (official verifiers, reproducible builds, SHA-tracked submissions).
- A fast C bit-parallel simulator for the DQI circuit (full 102,400-shot checks in seconds), plus an ANF/Toffoli tracer used to audit uncomputation.
- Fast evaluators for the chemistry metrics (Lambda, E_gap, lambda_eff, error_1/2) and a reproducible, code-only chemistry generator (byte-identical rebuilds, Apache-2.0 provenance) suitable for prize eligibility.
- A written record of leads we have measured and closed (DQI failure-budget and measure-and-repay variants, chemistry Lambda/Q-lock levers, loader decompositions), so a merged team does not repeat them.

What we are looking for: a team with complementary strengths, especially a different DQI construction or a new chemistry loader/factorization idea, where combining work would move the overall score rather than duplicate it.

Per the rules we will not share private code or submissions before a formal merger. Happy to talk about fit, combined submission counts (they must fit within the merged daily limit), and how work would be split. Reply here if interested.
