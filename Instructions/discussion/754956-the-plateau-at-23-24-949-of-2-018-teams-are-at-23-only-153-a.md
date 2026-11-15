# The plateau at 23-24: 949 of 2,018 teams are at 23+, only 153 at 24+ (LB anatomy + open data)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754956 — copied 2026-11-15 (raw Markdown, unedited).

**Harshit Kumar** · 2026-11-02 12:55 UTC · -2 votes

Snapshot of the public LB on 2026-11-02 12:10 UTC (2,018 teams). Measured: 1,360 teams >= 20, 949 >= 23, 153 >= 24, 23 >= 25, 10 >= 26, top 27.996. Five days ago only about 20 teams were at 23+, so the wall moved fast as public zips got rebuilt. My reading (not tested): rank inside the 23-24 band mostly shows who resubmitted the newest public zip. Notebook: https://www.kaggle.com/code/harshitsama/quantum-leaderboard-anatomy . Data (anonymous, no names): https://www.kaggle.com/datasets/harshitsama/kaggle-public-leaderboards-2026-11-02 . Also a list of 276 public notebooks with votes and the scores claimed in titles (claims unverified): https://www.kaggle.com/datasets/harshitsama/quantum-circuit-public-notebooks-2026-11-02 . I will re-snapshot weekly. Corrections welcome.


---

## Comments (1 threads)

Replies are indented with `>`, one level per reply depth.

**pjmathematician** · 2026-11-02 14:47 UTC · 2 votes

> Does submitting more help?

It's fairly easy to locally approximate the lb score. There might be deeper strategies ppl would be following to submit/not submit. For example, submitting to prove a better lb score is attached to a notebook, and all of its forks get submitted to get a better rank.

For #1 rank (longest lb leader prize, middle milestone prize), there are other strategies that the top ranked ppl would be following, like one thing I noticed - Breaking1800 only submits immediately when someone (usually Jack) submits to overtake and get #1, so they are always keeping a better local solution ready, which is fair.

I myself am not following any particular strategy. I was submitting a lot previously, but now I just submit if the trailing rank is getting too close to my score, or if I am overtaking the person above me with a confident margin, except if it is #1, in which case I ideally maintain #2.
