# Team changes to the dqi generator

This folder is our working copy of Syed Asad Ali's v12 dqi decoder generator (Apache-2.0; `LICENSE` and `NOTICE`
kept verbatim; Syed Asad Ali is the original builder of this decoder line and is credited first). The pristine
source is `third_party/syed3000-vqcc-own-circuit-portfolio-v12/dqi/dqi_source/`. No file of the original tree is
edited: every change is a new file that patches the pinned pipeline at run time, after its SHA-256 checks.

| File | What it does | Added |
|---|---|---|
| `team_build.sh` | Build script derived from `build.sh`: same host-file pins, runs `src/team/build_team.py` with a chosen extractor cutoff `J_MIN` | 2026-11-18 (S2.5) |
| `src/team/build_team.py` | Build with a configurable cutoff: replaces `integrated_pipeline.emit_key_prefix` at run time so the extractor receives `j_min`; records instead of asserting the three j_min-15 pins (static T 384,003, source ledger, event inventory) | 2026-11-18 (S2.5) |

## Builds

| J_MIN | SHA-256 | Static T | Notes |
|---|---|---:|---|
| 15 | `7cbe8ecf…` | 384,003 | byte-identical to v12 (`build.sh`): neutrality check of the team layer |
| 16 | `bec4d1a7…` | 383,887 | |
| 17 | `cc0d7d25…` | 383,771 | **selected**; byte-identical to Ertuğrul Özer's public cutoff-17 file |
| 18 | `32b0a8b4…` | 383,655 | same expected cost as 17, more variance |

Rebuild the selected one (WSL, CPython 3.13 venv, ~95 s, 1 GB):

```bash
PY=~/venvs/gqai313/bin/python WORK=~/work/dqi-team bash generators/dqi/team_build.sh 17 \
  <repo>/google-verifiable-quantum-circuit-compilation/kickmix <repo>/google-verifiable-quantum-circuit-compilation/verifiers \
  OUT_DIR cc0d7d25db72e90130c7c739f874a7ac611751791ee6e6c765d4d881c01b7a52
```
