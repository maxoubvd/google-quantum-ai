# Provenance of third-party code and circuits

Every file reused from another participant is listed here before it is used in a submission.
A file is **prize-eligible** only if we can rebuild it byte-for-byte from source we hold
(host answer, forum topic 752010).

Each upstream work lives in `third_party/<author>-<work>/` with its original `LICENSE` and `NOTICE`
kept verbatim. Files we change carry a header line saying so (Apache License 2.0, section 4(b)).

| Work | Author (credit chain) | Source (Kaggle ref, version, date pulled) | License | Files and SHA-256 | Generator held? | Rebuild verified? |
|---|---|---|---|---|---|---|
| Builders v12: invsqrt 46 Q, dqi 889 Q, chemistry | Syed Asad Ali (must be named first in any lineage, NOTICE §4(d)) | dataset `syed3000/vqcc-own-circuit-portfolio` **version 12**, pulled 2026-11-17 | Apache-2.0 + NOTICE | `invsqrt_source_zip` 7d9801ae20e942dba2211a656bba6014a0a16d57da534c0b60328d1c5d9d8254; `dqi_source_zip` 49e2b331d4008afb724985caf1038765b7a910f08838b590bf6e8fb692772763; `dfthcsos_source_zip` a3e5747c6211c465cc211df4608beab59898e806d1dbddc5624b674dae8d3fea (match the dataset's SHA256SUMS). Outputs: invsqrt 321aff83…, dqi 7cbe8ecf… | **Yes** (`build.sh` per task) | **Yes**, 2026-11-17: all 17 outputs byte-identical (invsqrt 321aff83…, dqi 7cbe8ecf…, chemistry 15 loaders + npz). The chemistry `submission.npz` (trained weights) is a frozen input copied by the build, not regenerated |
| Portfolio v13: invsqrt 44 Q (75,170.66), dqi 890 Q, chemistry 4,040/pass | Syed Asad Ali | same dataset, version 13 (current) | Apache-2.0 + NOTICE | compiled `.kmb` only | **No** (v13 ships no source; README says sources stay in v12) | — |
| invsqrt `563d0583` (44 Q, 74,560.91), current public best | maoxin123, on Syed Asad Ali's generator line, Ali Haydar Özdağ's 64-row table, Ertuğrul Özer's angle refit | notebook `maoxin123/qcc-28-17-4-022-toffoli-chemistry-loaders`, pulled 2026-11-17 | Apache-2.0 (notebook) | embedded as a compressed blob | **No** (compiled file only) | — |
| dqi `575b043f` / `61c8348d` (890 Q), current public best | GrantPrater; inputs by boristown (`d54d3387`) and Balram Srinivasan (`66b4f18b`), on Syed Asad Ali's v13 decoder | dataset `grantprater/vqcc-dqi-61c8348d-merged-closings`, pulled 2026-11-17 | Apache-2.0 + Syed Asad Ali NOTICE | `dqi_575b043f.kmb.xz` + patch tools (`vqcc_*.py`, chains/pairs JSON) | **Partly**: patch tools rebuild it from two compiled inputs, which have no generator | Not yet |
| Chemistry: 15 loaders at 4,022 Toffolis/pass + weights `703843d3` | maoxin123 (loaders, CC0), Atakan Aldemir (weights), on Lera's (8,143,44) table | `atakanaldemir/qcc-shared-4022-loaders-and-own703-weights`; `marculera/dfthc-c44-build`, pulled 2026-11-17 | Apache-2.0 (Lera NOTICE §4) / CC0 (maoxin123 loaders) | `science17.xz` | **No** for the loaders (compiled); weights come from a refit run | — |
