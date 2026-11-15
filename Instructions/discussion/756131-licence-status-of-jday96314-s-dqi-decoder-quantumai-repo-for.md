# Licence status of jday96314's DQI decoder (QuantumAI repo) for prize-eligible submissions

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/756131 — copied 2026-11-15 (raw Markdown, unedited).

**i4nLe3** · 2026-11-05 14:18 UTC · 1 votes

Hi all, and a question for the hosts in particular,
  
Many public DQI circuits trace back to @jday96314's continued-fraction decoder(build_raw_cf_destructive.py, from github.com/jday96314/QuantumAI). This includes the 891-qubit line distributed in several Apache-2.0 Kaggle datasets. As far as I can tell, the original repository states no licence and is no longer reachable.
The rules require winners to license their submission and its source under Apache 2.0 (B5/B7) and to warrant it is their own work and free of third-party rights (C14a). 

Is a submission derived from this lineage eligible for prizes, including the Nov 20 milestone, provided it is fully credited? Or should prize candidates use an independent implementation?
  
A clear answer would help everyone near the top, since most leading DQI entries seem to share this ancestry.
  
Thanks!


---

## Comments (3 threads)

Replies are indented with `>`, one level per reply depth.

**James Day** · 2026-11-05 23:19 UTC · 7 votes

I believe the public reimplmentations are fine to use/relicense as Apache 2.0. I will make the original repository public under Apache 2.0 if the hosts are concerned.

I only made the repo private to prevent people from copying me as I commit additional optimizations. I don't intend to get anyone disqualified for my mistake (unintentionally configuring it to be public at first).

> **D Chakraborty** · 2026-11-06 04:56 UTC · 2 votes
>
> ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2F82900c373311ec3c9219b8e03a358155%2FApache-2.0%20Under%20Pressure.png?generation=1793940992448900&alt=media)

> > **Syed Asad Ali** · 2026-11-12 21:32 UTC · 1 votes
> >
> > ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F30764058%2F9a202ec2db8a9b1dc9e3e1473f745cb5%2Fexec-ebb3f712-279a-42bb-b2b4-aaf37e42e83e.png?generation=1794519142012711&alt=media)
> > No guns needed on our side 😄 Our three generators (fusion, DQI and chemistry) are Apache-2.0 from day one, and the complete source is public. Fork it, rebuild it, improve it.
> > https://www.kaggle.com/code/syed3000/vqcc-own-circuit-portfolio-notebook

> **Varnavas Tryfonos** · 2026-11-09 17:31 UTC · 0 votes
>
> Hi James, thanks for clarifying that the public reimplementations of your decoder are fine to use and relicense as Apache 2.0. For prize-eligibility paperwork, would you be willing to publish a snapshot of the original repository as it was when it was public (just that commit, none of your later work) under Apache 2.0, as you offered? Or else state the licence of that version in a LICENSE file somewhere public? That would let everyone building on the public line cite a clear licence. Thanks for sharing the decoder in the first place.

> > **James Day** · 2026-11-10 01:14 UTC · 0 votes
> >
> > Done: https://github.com/jday96314/QuantumAI-archive

**John K** · 2026-11-10 02:58 UTC · 0 votes

Thanks for raising this before the Nov 20 milestone. James's note that the public reimplementations can be treated as Apache 2.0, and that he will open the original under Apache 2.0 if needed, is useful. I would still like a one-line host confirmation that a prize submission which credits that lineage stays eligible under B5/B7 and C14a, so prize candidates are not left guessing.

**(deleted user)** · 2026-11-08 05:11 UTC · 0 votes

*(deleted)*
