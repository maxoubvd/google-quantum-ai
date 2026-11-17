# VQCC Own-Circuit Portfolio v12 · a complete from-source package

**Public score 27.00209.** Geometric mean **27.0784** at zero failed DQI shots.

This is our 100% original setup, released for the first time. All three circuits are produced by builders we wrote
ourselves, and nothing in them is taken from any public code.

## Fork it and make it better

You are warmly invited to use this dataset and its companion notebook, *VQCC Own-Circuit Portfolio Notebook*, to
improve any of the three circuits in any way you can. Everything here is Apache-2.0 and free to use. If you build on
it, please credit the lineage starting from Syed Asad Ali.

## What is inside

| File | What it is |
|---|---|
| `invsqrt.kmb` | the fusion (inverse square root) circuit |
| `dqi.kmb` | the DQI circuit |
| `dfthcsos__submission.npz`, `dfthcsos__angle_k_0.kmb` … `dfthcsos__angle_k_14.kmb` | the chemistry submission (`__` stands for `/` in the submission zip) |
| `invsqrt_source_zip`, `dqi_source_zip`, `dfthcsos_source_zip` | the complete source code of the three builders, one zip archive each (saved without the `.zip` extension so that Kaggle keeps each archive as one file) |
| `SHA256SUMS` | the SHA-256 of every file above |

**Rebuilding a circuit yourself:** unzip one of the source archives and run
`bash build.sh <competition kickmix folder> <competition verifiers folder> <output folder>`. Its `README.md` lists the
expected SHA-256 of every output file. The companion notebook does this for all three tasks and checks every byte.

## Results (measured with the official verifiers)

| Task | Qubits | Mean Toffolis per shot | Cost | Ratio to reference |
|---|---:|---:|---:|---:|
| Fusion (inverse square root) | 46 | 1,511.57 | 80,162.76 | **25.91×** |
| DQI (Reed–Solomon decoding) | 889 | 383,156 | 3.6073 × 10¹¹ | **67.75×** |
| Chemistry (DF-THC, FeMoCo) | 427 | 14,810 | 296,635,967.59 | **11.31×** |

The DQI cost shown is at zero failed shots. Kaggle's own draw for this file contains two failed shots, which is why the
public score is 27.00209 rather than 27.0784.

## Credits and licence

Copyright 2026 Syed Asad Ali. Licensed under Apache-2.0.

**Credits:** Syed Asad Ali, working with a team of AI agents (GPT-6 Astra Pro, Gemini 4 Argon, Fable 5.2, Mythos 5.2
and Opus 5.5), with Claude Science supporting the research.
