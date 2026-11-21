# Our base on the ThemBones source line (plan step S2.26)

Sources (pristine, never edited):

- `third_party/thembones-milestone-generators-v2/`: team ThemBones (Danish Ali), Kaggle dataset
  `danishlynxthembones/vqcc-thembones-milestone-generators` version 2 (2026-11-20), Apache-2.0, with its NOTICE and the
  `third_party__*` upstream licences (Syed Asad Ali first). Every file matches its `GENERATORS_MANIFEST.sha256` (only
  Kaggle's `dataset-metadata.json` is absent). `SET_*` files and `set_invsqrt/invsqrt.kmb` (c639a045) come from the
  companion dataset `danishlynxthembones/vqcc-thembones-milestone-set`.
- dqi: Tony Li (tonylica), Kaggle dataset `tonylica/quantum-ai-32-85306-task-generators` (2026-11-21, Apache-2.0,
  `notices/RELEASE_NOTICES.txt`), on Team ThemBones' designs, upstream Syed Asad Ali. Not yet copied into the repo
  (101 MB, of which 90 MB is the dqi stages): it sits in `downloads/datasets-2026-11-21/` (git-ignored); the files
  our build uses are pinned in `third_party/PROVENANCE.md`.

Our scripts (run inside WSL from the repo root):

| Script | What it does | Result (2026-11-21) |
|---|---|---|
| `setup_sandbox.sh` | Lays the ThemBones folders out in a sandbox home `~/tb` as their README asks (the `<HOME>` path redaction undone, nothing else changed); creates `~/venvs/tb313` (CPython 3.13.18, numpy 2.5.4, scipy 1.18.1, numba 0.68.1) | about 1 min |
| `check_chem_loaders.py` | The 15 chemistry loaders rebuilt from their `.kmx` and compared with the set's manifest | **15/15 byte-identical**, 0.6 s |
| `build_invsqrt.sh` | ThemBones' one-command invsqrt pipeline (`f9_i1_stack.sh`), with two inputs supplied (see below) | `e2787831…`, 39 Q, 1,379.877 Toffolis a shot; **exact proof PASS against c639a045 on all 3,133,452 inputs**; 322 s, 1.8 GB; deterministic (same bytes on two runs) |
| `build_dqi_tonylica.sh` | tonylica's dqi builder on an ext4 copy of his dataset | **`4b572141…` byte-identical** (the board dqi file), 19 s, 1.3 GB |

What is not reproducible from the published material, and why:

- invsqrt: the FC18-M9 device file `s18/cub/m9dev.json` (config key `_fc18m9`) and the tree `s18/srcR` are not
  published. `build_invsqrt.sh` removes the key (the generator emits its built-in m9 branch: +0.5 Toffoli a shot,
  39 × 0.5 = +19.5 cost units, +0.03% on invsqrt) and points `fc2_src` at `srcM`, as the pipeline's q18 recipe does.
  It also makes fc18_M's emit itself (used only by the first all-input check). Cost of our file:
  39 × (1,379.877 + 237.42) ≈ 63,075 against 63,055 for c639a045. The device can be re-derived with their tools
  (`invsqrt_fc_s18m9frame.py` + `invsqrt_fc_s18m9q.py`; the three affine forms are not on record): S2.27.
- dqi (ThemBones chain): three inputs are missing (`idle_l9r4dp_902.txt`, `spec_k1.txt`, `jt_tensors_r54c.json`), and
  the team never rebuilt its file A. Our dqi base is therefore tonylica's builder, whose inherited decoder is recovered
  KMX assembly (prize status to be asked, M3.6).
- chemistry: the loaders' `.kmx` and the table are carried as data (randomized searches, GPU fit); the assembly and
  gauge tools named in their README (`chem_c3_entry_raw.py`, `lead_gauge.py`) are not in the dataset. The board
  chemistry (goodpjw2008's chain, Ozer's table) is beyond this set anyway: S2.31.
