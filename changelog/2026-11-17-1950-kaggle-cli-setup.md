# Kaggle API client in both Python environments

- **Date:** 2026-11-17 19:50 UTC (approximate; entry written at 20:27 UTC)
- **Plan step:** Phase 1 prerequisites
- **Files:** `requirements.txt`; packages installed in `.venv/` (Windows) and `~/venvs/gqai` (WSL)

## What was done

- Installed `kaggle==2.2.5` in both environments. It failed at import with
  `ModuleNotFoundError: No module named 'six'`, a dependency it does not declare, so `six==1.17.0` was
  installed too. Both are pinned in `requirements.txt`.
- The user saved their API token to `~/.kaggle/access_token` on Windows and in WSL (the token was never
  shown to Claude).
- `.gitignore`: `downloads/` added for raw Kaggle downloads (untrusted, re-downloadable).

## Why

The Kaggle web pages render client-side only and cannot be read from this machine. The API gives the
leaderboard, the notebook list sorted by score, and the notebooks and datasets themselves.

## Results and analysis

`kaggle competitions list`, `kaggle competitions leaderboard -d` and
`kaggle kernels list --competition google-verifiable-quantum-circuit-compilation --sort-by scoreDescending`
work on both sides. The account is entered in the competition, with no submission yet.
Dataset versions are addressed as `owner/slug/<version>` (not `.../versions/<n>`).

## Next

Licensing and source inventory (S1.1, S1.2).
