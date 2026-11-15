# Public frontier check: gengsr 25.2268 zip reproduces exactly on LB; 2026-11-10 LB snapshot

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/757590 — copied 2026-11-15 (raw Markdown, unedited).

**Harshit Kumar** · 2026-11-10 04:15 UTC · 0 votes

A short reproducibility note.

I rebuilt gengsr's public 25.2268 portfolio from the Apache-2.0 dataset `gengsr/vqcc-25-2268-portfolio`: all 18 files matched the SHA-256 manifest, and I submitted the zip unchanged. The public LB returned **25.22676**, which is the number the authors reported. So the chain "published manifest -> rebuilt zip -> same LB score" holds here, with no local salt surprises. Credit goes to gengsr and the upstream authors listed in the dataset's NOTICE.txt.

Anonymised leaderboard snapshots (rank, score, submission count only) for this and three other live competitions, including today's: https://www.kaggle.com/datasets/harshitsama/kaggle-public-leaderboards-2026-11-02 . On 2026-11-10 this board had 2,729 teams and the top score was 29.301.

Question for people near the top: is the gap from the 25.2 public frontier to 29+ mostly the DQI task, the chemistry task, or invsqrt? A per-task view would help the rest of us choose where to spend effort.
