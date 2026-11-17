# Licensing files and source inventory

- **Date:** 2026-11-17 20:10 UTC (approximate; entry written at 20:27 UTC)
- **Plan step:** S1.1 and S1.2
- **Files:** `LICENSE`, `NOTICE`, `third_party/PROVENANCE.md` (created); `downloads/` (git-ignored)

## What was done

- `LICENSE`: the Apache-2.0 text, byte-identical to `verifiers/LICENSE.txt`.
- `NOTICE`: "Copyright 2026 Maxime Bendavid (Kaggle: maxouUu)", the scope of the license (our own code
  only), and the third-party material it does not cover (competition data, `third_party/`,
  `Scientific Papers/`, `Instructions/`, the Gemini PDF).
- Downloaded into `downloads/` the notebooks and datasets behind the best public files: Syed Asad Ali
  (dataset versions 13 and 12, notebook), maoxin123 (two notebooks), GrantPrater (dataset with patch
  tools, notebook), Balram Srinivasan, boristown, Atakan Aldemir, Lera (marculera), Ali Haydar Özdağ,
  Asterios Panoras, alex chilton.
- `third_party/PROVENANCE.md`: one row per source with author chain, Kaggle ref and version, license,
  SHA-256 and whether a generator is published.

## Why

The rules license the winning submission and its generating code under Apache-2.0; reused Apache-2.0
work must keep its LICENSE and NOTICE (section 4). The hosts require, for every prize, a generator that
reproduces the exact circuits (forum topic 752010).

## Results and analysis

- **Only Syed Asad Ali's dataset version 12 ships full builder sources** (`invsqrt_source_zip`,
  `dqi_source_zip`, `dfthcsos_source_zip`, each with a `build.sh`; SHA-256 match the dataset's
  `SHA256SUMS`). They build the 12 Nov state: invsqrt 46 Q / 80,162.76, dqi 889 Q / 3.6073e11,
  chemistry 2.9664e8, geometric mean ≈ 27.08.
- Version 13 (44 Q invsqrt, 890 Q dqi) ships compiled files only.
- The 28.18 frontier parts are compiled files (maoxin123's invsqrt `563d0583` and chemistry loaders) or
  patches of compiled files (GrantPrater's dqi `575b043f`, rebuilt by his tools from two compiled
  inputs that have no generator).
- Syed Asad Ali's NOTICE requires naming him first in the lineage of any derived work.
- Medals do not require generators; prizes do. Prizes need well above the plateau anyway.

## Next

S1.3: user decision on the two tracks (A: submit the credited public frontier now; B: build on the v12
generators).
