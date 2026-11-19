# Source Credits, Provenance, and Acknowledgments

This submission and compilation pipeline by **boristown** (team **暗黑AGI**, https://www.kaggle.com/code/boristown/qcc-gap-3-residue-repair-and-c43-native-compiler) builds upon public Apache-2.0 and CC0 work shared in the `google-verifiable-quantum-circuit-compilation` competition. We sincerely thank every author below whose circuits, loaders, mathematical reductions, verification frameworks, and analyses enabled this portfolio. All upstream licenses and `NOTICE` texts are preserved; unmodified official verifiers and `kickmix` binaries remain Google Quantum AI's competition runtime (https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation).

---

## Apache-2.0 Section 4(d) Notice (`syed3000/vqcc-own-circuit-portfolio`)

```text
VQCC own-circuit portfolio
Copyright 2026 Syed Asad Ali. Licensed under the Apache License, Version 2.0.

Credits: Syed Asad Ali, working with a team of AI agents (GPT-6 Astra Pro, Gemini 4 Argon,
Fable 5.2, Mythos 5.2 and Opus 5.5), with Claude Science supporting the research.

Lineage notice (Apache-2.0, Section 4(d): this NOTICE must be carried by every redistribution):
- The 46-qubit inverse-square-root generator, the 889-qubit DQI decoder generator (v12 and
  later), and all circuits built from them are the original work of Syed Asad Ali
  (Kaggle: syed3000; dataset syed3000/vqcc-own-circuit-portfolio).
- This version contains: invsqrt.kmb (44-qubit inverse square root, from that generator line),
  dqi.kmb (890-qubit DQI decoder, from that decoder line), and the chemistry files
  dfthcsos__angle_k_0..14.kmb and dfthcsos__submission.npz (our own loader compiler and word tables).
- Any redistribution of these files, or of works derived from them at any number of steps
  removed, must include this NOTICE (or its full attribution text) in a readable place, and
  must name "Syed Asad Ali" as the original builder in its credit or lineage line, before
  any later contributors.
- Improvements by others are welcome; please credit them after the original builder.
```

---

## Component Provenance & Absorbed Solutions

### 1. Inverse Square Root (`invsqrt.kmb`)
- **Original builder & direct source (Apache-2.0)**: **Syed Asad Ali** (`syed3000`), original builder of the 46-qubit and 44-qubit inverse-square-root generator line (`47b31838...`, `44` qubits), **Ali Haydar Özdağ** (`ahoirg`, 64-row table), **maoxin123** (`864f81cd20489d80e302f48014fcedec3b3fb703122fbcbac41b0bb2a395a859`, `44` qubits, cost `74,490.95`, 1,453.375 mean Toffolis, 33 fine rotations, `worst error/tolerance 0.99946`; prior builds `563d0583...` at `74,560.91` and `f2cdab0b...` at `75,006.38`, Apache-2.0):
  - Original Notebook / Dataset: https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook, https://www.kaggle.com/datasets/syed3000/vqcc-own-circuit-portfolio
  - Absorbed 44Q `864f81cd` & `563d0583` / `f2cdab0b` Builds: https://www.kaggle.com/code/maoxin123/qcc-28-18-4-022-toffoli-chemistry-loaders, https://www.kaggle.com/code/maoxin123/qcc-v27-a-44-qubit-invsqrt-on-the-28-09-frontier
- **Prior & concurrent inverse-square-root contributors thanked & preserved in history**:
  - **jeffemerson**: `9339_plus_angles` (`cmerge` + stack-AND merge cut, `2c125b39...`, https://www.kaggle.com/code/jeffemerson/28-071-44-qubit-stack-our-and-merge-cut)
  - **Leo Tang** (`leotang09`): https://www.kaggle.com/code/leotang09/vqcc-t7421-submission (`a4d46b17...`)
  - **James Barlow** (`manwithacat`): Hedge Wizardry `HW-26` & `HW-27` / `qhdl` (https://www.kaggle.com/code/manwithacat/vqcc-hedge-wizardry-notebook, https://www.kaggle.com/datasets/manwithacat/vqcc-hedge-wizardry, https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
  - **Asterios Panoras** (`asteriospanoras`): measured-CZ dead-write, `q21` HMR analysis, & LP angle refit (https://www.kaggle.com/code/asteriospanoras/qcc-27-3122-frontier-finer-chemistry-refit, https://www.kaggle.com/code/asteriospanoras/qcc-top-1-public-solution-until-not-lb-28-075)
  - **Ertuğrul Özer** (`ozertuu`): fast runner, rotation-angle refits, and verified SOTA portfolios (`8b410e77...`, https://www.kaggle.com/code/ozertuu/25-295-i-m-better-verified-sota-portfolio, https://www.kaggle.com/code/ozertuu/vqcc-verified-portfolio-28-056, https://www.kaggle.com/code/ozertuu/vqcc-verified-portfolio-28-082-sota, https://www.kaggle.com/code/ozertuu/vqcc-super-sota-verified-portfolio-28-091)
  - **Atakan Aldemir** (`atakanaldemir`): base notebook & C43/residue-repair lineage (https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler)
  - **Ali Haydar Özdağ** (`ahoirg`): 64-row-table build

### 2. DQI Decoder (`dqi.kmb`)
- **Original builder & direct source (Apache-2.0)**: **Syed Asad Ali** (`syed3000`), original builder of the 889-qubit (v12+) and 890-qubit (v13/v14) DQI decoder generator line (`46446af1...`, `890` qubits), **Hayato Hikawa** (`hayatohikawa`, `qcc-c8-milestone-portfolio` `4b350932b03931af6040c871209da62aa0ebc74cfac3e77efcb6384e027285cc`), **Balram Srinivasan** (`balramsrinivasan/vqcc-28-10-norm-preserving-dqi-closings`, `66b4f18b`), **Grant Prater** (`grantprater/vqcc-dqi-61c8348d-merged-closings`, 49-site GF24 closing splice + 4-pair window-proven product retention `575b043f`), **Akmal Xodarev** (`busyaprime`, product-chain caching idea: https://www.kaggle.com/code/busyaprime/lb-25-33133-grimoire-33-every-gate-a-pixel), and **boristown** `cmerge` + `8 x HW-27` native patch (`ccc2ab84...`) + 19-chain / 144-chain pure-NumPy product-cache patches + 3-site `02_chain` & condition-block `pair_window` measurement-out patches (`bb9bf415...` `-3` Toffoli-likes & `5f901c556cafd02dc11641eee7766171a454e3a0b6687d890ff9f99944e1934c` `-6` Toffoli-likes from `575b043f`, `890` qubits, verified `ve=0` on official shots & Kaggle public LB):
  - Original Notebook / Dataset: https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook, https://www.kaggle.com/datasets/syed3000/vqcc-own-circuit-portfolio
  - Absorbed C8 Milestone Portfolio: https://www.kaggle.com/datasets/hayatohikawa/qcc-c8-milestone-portfolio
- **Prior DQI lineage contributors thanked**:
  - **James Barlow** (`manwithacat`): `HW-27` / `qhdl` `(CCX -> HMR)` elimination & Friends overlay (https://www.kaggle.com/code/manwithacat/the-real-150k-was-the-friends-we-made)
  - **Heliosli**: https://www.kaggle.com/code/heliosli/quantum-ai-v2 & `qcc-circuits-live`
  - **Kirill Tymchishin** (`kirilltym`): https://www.kaggle.com/code/kirilltym/quantum-circuits-retained-product-audit
  - **Areeb Ali** (`areebbb`): https://www.kaggle.com/code/areebbb/qai-fewer-toffolis-same-qubits-lb-27-31224
  - **Atakan Aldemir**, **Md. Hamid Hosen**, **Agniv Banerjee (`somejoint`)**, **Ali Haydar Özdağ**

### 3. Chemistry (`dfthcsos/angle_k_0..14.kmb` and `dfthcsos/submission.npz`)
- **Direct 4,022-Toffoli/pass & 4,025-Toffoli/pass loaders & base factorization (`Apache-2.0` / `CC0`)**: **maoxin123** (`4,022` AND/pass sequential XOR-difference loaders `k_0` 292 T & `k_1` 287 T in [`maoxin123/qcc-28-17-4-022-toffoli-chemistry-loaders`](https://www.kaggle.com/code/maoxin123/qcc-28-17-4-022-toffoli-chemistry-loaders) / [`maoxin123/qcc-28-18-4-022-toffoli-chemistry-loaders`](https://www.kaggle.com/code/maoxin123/qcc-28-18-4-022-toffoli-chemistry-loaders), and **ghazarosbarseghyan91** packaging), **Hayato Hikawa** (`hayatohikawa/qcc-c8-milestone-portfolio` for `k_2..k_14`), and **Lera (`marculera`)** (`DF-THC 11.38x: (8, 143, 44)` coded-word table `906fa220` and XOR-of-five-words `4,025` loaders, upstream Apache-2.0 `NOTICE` preserved as `CHEM_SOURCE_NOTICE`):
  - Notebook: https://www.kaggle.com/code/marculera/df-thc-11-38x-codes-on-the-xor-of-five-words
  - Dataset: https://www.kaggle.com/datasets/marculera/dfthc-c44-build
  - Prior `4,028`/`4,030`/`4,033` releases built upon **maoxin123**: https://www.kaggle.com/code/maoxin123/qcc-v27-4028-toffoli-coded-word-chemistry, https://www.kaggle.com/datasets/maoxin123/qcc-chemistry-4028-toffoli-loaders, https://www.kaggle.com/code/maoxin123/qcc-v27-4030-toffoli-coded-word-chemistry, https://www.kaggle.com/datasets/maoxin123/qcc-chemistry-4030-toffoli-loaders, https://www.kaggle.com/code/maoxin123/qcc-v27-4033-toffoli-coded-word-chemistry, https://www.kaggle.com/datasets/maoxin123/qcc-chemistry-4033-toffoli-loaders
- **Upstream chemistry lineage & foundational methods**:
  - **Syed Asad Ali** (`syed3000`): https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook
  - **Lera** (`marculera`): compact DF-THC & 5-word XOR loader codes (https://www.kaggle.com/code/marculera/df-thc-11-2x-codes-on-the-xor-of-five-words, https://www.kaggle.com/code/marculera/dfthc-compact-pipeline)
  - **Ertuğrul Özer** (`ozertuu`): *quantum-sim-linalg-cache* (https://www.kaggle.com/datasets/ozertuu/quantum-sim-linalg-cache) and `qccscope` diagnostics
  - **James Barlow** (`manwithacat`): `dfthc_refit` & public circuits (https://www.kaggle.com/datasets/manwithacat/vqcc-apache-public-circuits)
  - **Atakan Aldemir** (`atakanaldemir`): base compilation notebook & pi-gauge / phase folds (https://www.kaggle.com/code/atakanaldemir/qcc-gap-3-residue-repair-and-c43-native-compiler, https://www.kaggle.com/code/atakanaldemir/qcc-terminal-readout-and-phase-folds)
  - **Kirill Tymchishin** (`kirilltym`): exact gauge + pi-gauge weight refinement (https://www.kaggle.com/code/kirilltym/quantum-circuits-exact-gauge-5-word-qrom, https://www.kaggle.com/code/kirilltym/quantum-circuits-pi-gauge-weight-refinement)
  - **prvsiyan**: verifier-guided optimization & salt-surviving circuits (https://www.kaggle.com/code/prvsiyan/verifier-guided-quantum-circuit-optimization, https://www.kaggle.com/datasets/prvsiyan/vqcc-verified-best-of-circuits-1106a)
  - **Ahmed Berat Özer** (`ahmedberatozer`): DF-THC group compression (https://www.kaggle.com/code/ahmedberatozer/quantum-ai-df-thc-group-compression)
  - **Asterios Panoras** (`asteriospanoras`): https://www.kaggle.com/code/asteriospanoras/qcc-27-3122-frontier-finer-chemistry-refit
  - **Areeb Ali** (`areebbb`): https://www.kaggle.com/code/areebbb/qai-fewer-toffolis-same-qubits-lb-27-31224
  - **Akmal Xodarev** (`busyaprime`): https://www.kaggle.com/code/busyaprime/lb-25-33133-grimoire-33-every-gate-a-pixel
  - **Agniv Banerjee** (`somejoint`): https://www.kaggle.com/code/somejoint/quantum-circuits-24-65-local-gm-chem-refit
  - **alex chilton** (`alexchilton`): https://www.kaggle.com/code/alexchilton/24-94-lb-verified-three-task-portfolio, https://www.kaggle.com/code/alexchilton/27-31-lb-which-circuits-survive-private-salt
  - **Heliosli** (`heliosli`): https://www.kaggle.com/code/heliosli/quantum-ai-v2
- **Continuation optimization (`boristown` / 暗黑AGI)**:
  - Exhaustive forward/reverse symbolic HMR word recovery (`src/prepare_4030_base.py`) and Stage A/B projected bounded reduced-gradient L-BFGS with explicit Schur zero-eigenvalue nullspace chart (`src/refit_dfthc_v2.py`), verified against unmodified `dfthc_math` and stock Kaggle `dfthcsos` N4 verifiers.


## This version
Combines maoxin123 44-qubit invsqrt ecf908fe (74,383.70, Apache-2.0), Manvender Raparia chemistry weight refit 6ef06081 + maoxin123 4,022 loaders (via leotang09/vqcc-t7421-submission), and boristown board-proven DQI e39b1396 (-1 Toffoli pair-window AND retention).
