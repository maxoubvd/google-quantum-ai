# First public 27+ notebook: a complete from-source package (Apache-2.0). Fork it and make it better

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/758815 — copied 2026-11-15 (raw Markdown, unedited).

**Syed Asad Ali** · 2026-11-12 21:39 UTC · 17 votes

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F30764058%2F642c10b62adae8a8edfa95e62f254375%2Fexec-6d46c2bd-11db-418e-b0c0-3d1dabdf859f.png?generation=1794519544744189&alt=media)

Hi everyone,

Our notebook **VQCC Own-Circuit Portfolio v12** is now public, together with its dataset. As far as we know, it is the first public notebook above 27, and everything in it is built from our own source.

| Task | Qubits | Cost | Ratio |
|---|---:|---:|---:|
| Fusion (inverse square root) | 46 | 80,162.76 | 25.91× |
| DQI (Reed–Solomon decoding) | 889 | 3.6073 × 10¹¹ | 67.75× |
| Chemistry (DF-THC, FeMoCo) | 427 | 296,635,967.59 | 11.31× |

**Public score 27.00209** (27.0784 at zero failed DQI shots).

### What "from source" means here

- **The complete source code** of the three builders is in the dataset, not only the finished circuit files.
- **The notebook rebuilds** every circuit from that code on Kaggle and checks each file byte for byte before it writes `submission.zip`.
- **Nothing is hidden.** Fork it and you get the whole working package: read it, change it, run it and submit your own improved version.

This is our 100% original setup, released for the first time, and nothing in it is taken from any public code. It is Apache-2.0: you are warmly invited to fork it and improve any of the three tasks in any way you can. If you build on it, please credit the lineage starting from Syed Asad Ali.

- **Notebook:** https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook
- **Dataset:** https://www.kaggle.com/datasets/syed3000/vqcc-own-circuit-portfolio

Questions and improvements are very welcome in the comments.


---

## Comments (3 threads)

Replies are indented with `>`, one level per reply depth.

**Syed Asad Ali** — **Topic author** · 2026-11-13 16:34 UTC · 2 votes

@manwithacat  big scoreboard shakeups since I posted last night. Lets see what your daily summary Rezops has to say today 🤪

> **manwithacat** · 2026-11-14 09:46 UTC · 1 votes
>
> The cat was sleeping on top of the Mac mini, so yesterday's briefing was delayed due to thermal shutdown. Today's should be out in a hour or two.

**HARSHA Lekkala1** · 2026-11-13 18:09 UTC · 0 votes

Incredible work on breaking the 27 barrier and open-sourcing the complete builders, Syed! I'm really impressed by the from-source approach. I’m currently looking to merge teams to push the optimizations even further. Let me know if you have any open spots or interest in teaming up!

> **Syed Asad Ali** — **Topic author** · 2026-11-14 10:30 UTC · 0 votes
>
> sure drop me a msg and we can work out something!

**Giovanni Baldon** · 2026-11-14 12:25 UTC · 1 votes

Thanks, this is really helpful
