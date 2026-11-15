# QVC v0.4.0: complete source generators, 46Q inverse / 889Q DQI / FeMoCo

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/759244 — copied 2026-11-15 (raw Markdown, unedited).

**Francisco Angulo de Lafuente** · 2026-11-14 13:22 UTC · 0 votes

I am sharing the complete source and reproducibility record for QVC Chemistry Compiler v0.4.0 with all competition participants.

Repository and selected release: https://github.com/Agnuxo1/qvc-chemistry-compiler/releases/tag/v0.4.0

The release provides generators and mathematical inputs for all 18 submitted files: a 46-qubit inverse-square-root circuit, an 889-qubit DQI decoder, and the deep512 FeMoCo chemistry recipe. The exact submission ZIP, source bundle, SHA-256 checksums, build instructions, licenses, full local verifier receipts and scientific figures are included. The original competition data and native runtime are obtained from this competition, rather than redistributed in the repository.

Local checks give a composite of 26.97516–27.01640 across all 16 measured component-cost combinations. That range is not a confidence interval or a leaderboard score. The inverse has two full 200,540,928-shot checks and an exact full-domain/HMR proof; DQI has four 102,400-shot checks, all PASS, with one value error retained in its cost; chemistry has two full 78,080-shot checks. Inverse and chemistry receipts are reused only because the newly generated files match their exact hashes.

The inverse and DQI generators are Syed Asad Ali's Apache-2.0 work, credited with upstream notices preserved. This project's contributions are chemistry optimization, integration, measurements and documentation by Francisco Angulo de Lafuente. Apache-2.0 is the selected competition handover option; separable original work is also offered under MIT. This does not claim invention of the upstream generators.

The older, different 26ebf39d inverse and e9f54be3 DQI generators remain a separate historical recovery item. They are not the v0.4.0 submission inputs. Local verification and public source availability do not certify prize eligibility or an official result.

Francisco Angulo de Lafuente · https://orcid.org/0009-0001-1634-7063


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**Syed Asad Ali** · 2026-11-14 17:32 UTC · 0 votes

Thanks, Francisco, for the clear credit in this post. Two small fixes so the lineage survives forks of your repo: (1) your NOTICE lists me only among several circuit sources, and (2) docs/ATTRIBUTION.md still names “gengsr on maoxin123’s circuit” for the inverse square root and “Atakan Aldemir” for DQI, with me only as upstream.

As this post says, the 46-qubit inverse-square-root and 889-qubit DQI generators are my Apache-2.0 work. Could you add near the top of NOTICE, and in those two ATTRIBUTION rows: “The 46-qubit inverse-square-root generator and the 889-qubit DQI generator are Syed Asad Ali’s Apache-2.0 work (dataset syed3000/vqcc-own-circuit-portfolio); derived copies must keep this notice.” The chemistry is rightly credited to Lera’s pipeline and your own work.

Thanks again.
