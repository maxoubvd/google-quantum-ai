# what to learn from this competition

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/754652 — copied 2026-11-15 (raw Markdown, unedited).

**hengck23** · 2026-11-01 03:00 UTC · 35 votes

This is a rare competition that 
- many are using the same AI agent   
- many are unfamiliar with the background knowledge  
- the maths/modeling is rather difficult  
- AI agent is in general, is better than the human participant  

You can copy others' circuits and fuse them with yours, then make some improvements using your CC or Codex. But in the end, ask yourself, "What did you learn?", "Do you have a competitive advantage (after this competition)?".  

Do you just want to spend tokens chasing the leaderboard for three months, and in the end you (and your agent) learn literally nothing?  

So what can we learn from this competition?

----

1) Learn how to control AI (espically in this competition where AI is better than human) ... that is not enough. Learn how to control AI BETTER than others.  

2) Setup some RSI, etc
- How AI can implement your idea correctly. How it can make changes to it to make it work when it doesn't
- How AI can absorb all results from the forum, open notebook and create  something better at large scale and fast time ( distill public knowledge)
- update AI skill.md. learn how to do that
- learn how to create AI swarm


---

## Comments (11 threads)

Replies are indented with `>`, one level per reply depth.

**Asterios Panoras** · 2026-11-13 00:13 UTC · 2 votes

This is probably one of the most interesting discussions in this competition for me.

As a student without a strong background in quantum computing, AI handles much of the implementation. But I think there's a lot more to learn here than simply getting an agent to improve a score.

One of the most valuable things I'm learning is when to keep refining an approach and when to abandon it for something fundamentally different. AI can get stuck exploring small variations of the same idea, and recognizing that, organizing experiments, validating results, and extracting useful knowledge from public notebooks all require judgment.

At the same time, this competition has made mathematical concepts that previously felt abstract to me much more meaningful. Seeing how they contribute to actual optimizations motivates me to understand them better, even if learning the underlying theory takes time.

If many of us have access to similar AI tools, knowing how to use them effectively is certainly a competitive skill. But I think the more lasting advantage comes from combining that skill with domain knowledge, critical thinking, and the ability to recognize when a different direction is needed.

And perhaps that's the most valuable thing we can take away from this competition: skills that remain useful even when we face completely different problems in the future.

I'm curious how others decide when an approach has reached its limits and it's time to try something fundamentally different.

> **michel emil** · 2026-11-13 10:16 UTC · 1 votes
>
> now we are on the verge of the fact that AI will discover new concepts not just repeating existing knowledge to users, look at this demo code the ai found concepts in seconds  that humans discover in decades .
> this code is directly inspired by the actual process used for quantum circuits improvments
>
> #!/usr/bin/env python3
> """
> Exact algebraic-relation discovery over F_p (p = secp256k1 field prime).
>
> Given samples of some variables, find EVERY polynomial relation of total degree <= d that they
> satisfy, by linear algebra over F_p: build the monomial matrix M[sample, monomial] and
> compute its nullspace mod p. No coefficient grid, no coefficient bound, no floats.
>
>   * nullity = number of independent relations of that degree
>   * Schwartz-Zippel: a polynomial that vanishes on a few random samples vanishes on the whole
>     sampled variety except with probability ~ degree/p ~ 2^-250, so held-out samples are a
>     proof-grade check.
>   * you need more samples than monomials: C(n_vars + d, d)
> """
> import itertools
> import random
> import sys
> import time
>
> # ============================================================================
> # PURE-PYTHON SECP256K1 REFERENCE ARITHMETIC
> # ============================================================================
> P_INT = 2**256 - 2**32 - 977
> N_INT = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
> GX_INT = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
> GY_INT = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
>
> def ref_add(P, Q):
>     if P is None: return Q
>     if Q is None: return P
>     (x1, y1), (x2, y2) = P, Q
>     if x1 == x2 and (y1 + y2) % P_INT == 0: return None
>     if P == Q:
>         lam = (3 * x1 * x1) * pow(2 * y1, -1, P_INT) % P_INT
>     else:
>         lam = (y2 - y1) * pow(x2 - x1, -1, P_INT) % P_INT
>     x3 = (lam * lam - x1 - x2) % P_INT
>     return (x3, (lam * (x1 - x3) - y1) % P_INT)
>
> def ref_mul(k, P=(GX_INT, GY_INT)):
>     R = None
>     while k:
>         if k & 1: R = ref_add(R, P)
>         P = ref_add(P, P); k >>= 1
>     return R
>
> P = P_INT
> signed = lambda c: c if c < P // 2 else c - P
>
> # ============================================================================
> # EXACT LINEAR ALGEBRA OVER F_p
> # ============================================================================
> def monomials(nv, deg, min_deg=0):
>     out = []
>     for d in range(min_deg, deg + 1):
>         out += list(itertools.combinations_with_replacement(range(nv), d))
>     return out
>
> def feature_rows(samples, mons):
>     rows = []
>     for s in samples:
>         row = []
>         for mo in mons:
>             v = 1
>             for i in mo:
>                 v = v * s[i] % P
>             row.append(v)
>         rows.append(row)
>     return rows
>
> def nullspace_mod_p(rows):
>     """Basis of {v : rows . v == 0 (mod p)} by Gauss-Jordan elimination over F_p."""
>     m = [r[:] for r in rows]
>     n, piv_cols, r = len(m[0]), [], 0
>     for c in range(n):
>         pr = next((i for i in range(r, len(m)) if m[i][c] % P), None)
>         if pr is None:
>             continue
>         m[r], m[pr] = m[pr], m[r]
>         inv = pow(m[r][c], -1, P)
>         m[r] = [x * inv % P for x in m[r]]
>         for i in range(len(m)):
>             if i != r and m[i][c]:
>                 f = m[i][c]
>                 m[i] = [(a - f * b) % P for a, b in zip(m[i], m[r])]
>         piv_cols.append(c)
>         r += 1
>         if r == len(m):
>             break
>     basis = []
>     for f in (c for c in range(n) if c not in piv_cols):
>         v = [0] * n
>         v[f] = 1
>         for i, c in enumerate(piv_cols):
>             v[c] = (-m[i][f]) % P
>         basis.append(v)
>     return basis
>
> def find_relations(samples, deg, holdout=None, min_deg=0):
>     """-> (monomials, nullspace basis). Samples: list of tuples of ints mod p."""
>     mons = monomials(len(samples[0]), deg, min_deg)
>     if len(samples) <= len(mons):
>         raise ValueError(f"need more than {len(mons)} samples for {len(mons)} monomials")
>     ns = nullspace_mod_p(feature_rows(samples, mons))
>     if holdout is not None:
>         H = feature_rows(holdout, mons)
>         assert all(sum(c * f for c, f in zip(v, row)) % P == 0 for v in ns for row in H), "relation failed on held-out data"
>     return mons, ns
>
> def show(v, mons, names, max_terms=8):
>     k = next(i for i, c in enumerate(v) if c)
>     inv = pow(v[k], -1, P)
>     base = [c * inv % P for c in v]
>
>     # Rescale by the integer 1..400 that makes the largest coefficient smallest.
>     # This recovers small integer coefficients (like 2, -28) hidden behind modular fractions.
>     best = min(range(1, 401), key=lambda t: max(abs(signed(c * t % P)) for c in base))
>     terms = [(signed(c * best % P), mo) for c, mo in zip(base, mons) if c]
>     mon = lambda mo: "1" if not mo else "*".join(names[i] for i in mo)
>
>     if max(abs(c) for c, _ in terms) >= 10**6 or len(terms) > max_terms:
>         return f"<{len(terms)} terms, coefficients are large field constants>"
>     return " + ".join(f"{c}*{mon(mo)}" for c, mo in terms) + " == 0"
>
> # ============================================================================
> # DEMOS
> # ============================================================================
> def records(n, rng):
>     out = []
>     for _ in range(n):
>         B = ref_mul(rng.randrange(2**100, 2**200))
>         T = ref_mul(rng.randrange(1, 17))
>         den = (T[0] - B[0]) % P
>         inv = pow(den, -1, P)
>         slope = (T[1] - B[1]) * inv % P
>         lam_sq = slope * slope % P
>         x3 = (lam_sq - B[0] - T[0]) % P
>         y3 = (slope * (B[0] - x3) - B[1]) % P
>         out.append(dict(Bx=B[0], By=B[1], Tx=T[0], Ty=T[1], denominator=den, inverse=inv, slope=slope,
>                         lam_sq=lam_sq, x3=x3, y3=y3))
>     return out
>
> def demo():
>     rng = random.Random(3)
>     recs = records(150, rng)
>     tr, ho = recs[:100], recs[100:]
>     pick = lambda rs, vs: [tuple(r[v] for v in vs) for r in rs]
>
>     print("A. x3, lam_sq, Bx, Tx   (linear):")
>     vs = ["x3", "lam_sq", "Bx", "Tx"]
>     mons, ns = find_relations(pick(tr, vs), 1, pick(ho, vs), min_deg=1)
>     print(f"   nullity {len(ns)} ->", [show(v, mons, vs) for v in ns])
>
>     print("\nB. y3, By, slope, x3    (linear):")
>     vs = ["y3", "By", "slope", "x3"]
>     mons, ns = find_relations(pick(tr, vs), 1, pick(ho, vs), min_deg=1)
>     print(f"   nullity {len(ns)}  (only the trivial zero vector: y3 = slope*(Bx-x3) - By is bilinear)")
>
>     print("\nC. y3, By, slope, Bx, x3 (degree <= 2):")
>     vs = ["y3", "By", "slope", "Bx", "x3"]
>     mons, ns = find_relations(pick(tr, vs), 2, pick(ho, vs))
>     print(f"   nullity {len(ns)}")
>     for v in ns:
>         print("   ", show(v, mons, vs))
>     print("    1st = chord slope through B and -(B+T);  2nd = both points on the curve (the constant b cancels)")
>
>     print("\nD. x-coordinates of P-G, P, P+G  (what a consecutive-key walk produces):")
>     xG = GX_INT
>     xs = lambda k: tuple(ref_mul(k + d)[0] for d in (-1, 0, 1))
>     tr = [xs(rng.randrange(10, 10**9)) for _ in range(100)]
>     ho = [xs(rng.randrange(10, 10**9)) for _ in range(50)]
>
>     for d in (1, 2, 3, 4):
>         mons = monomials(3, d)
>         if len(tr) > len(mons):
>             t = time.perf_counter()
>             _, ns = find_relations(tr, d, ho)
>             print(f"   degree <= {d}: {len(mons):2d} monomials, nullity {len(ns):2d}   ({time.perf_counter() - t:.2f}s)")
>
>     # degree-3 nullity is 5, so "the" relation is not one vector: test the classical formula directly
>     ssum = lambda s: ((s[1] - xG) ** 2 * (s[2] + s[0]) - 2 * xG * s[1] * (s[1] + xG) - 28) % P
>     print(f"   classical x(P+G)+x(P-G) formula (degree 3) vanishes on all samples: {all(ssum(s) == 0 for s in tr + ho)}")
>
>     prod = lambda s: (s[2] * s[0] * (s[1] - xG) ** 2 - (s[1] * xG) ** 2 + 28 * (s[1] + xG)) % P
>     print(f"   classical x(P+G)*x(P-G) formula (degree 4) vanishes on all samples: "
>           f"{all(prod(s) == 0 for s in tr + ho)}")
>
> if __name__ == "__main__":
>     demo()
> ####
>
> AI propmt
> Demo D: The Holy Grail (Montgomery x-Only Addition)
> You asked if there was a "lattice structure" or shortcut in a sequence of consecutive points. Demo D answers yes, and it is algebraic, not statistical.
> For three consecutive x-coordinates xP−Gx_{P-G}xP−G​, xPx_PxP​, and xP+Gx_{P+G}xP+G​, the script found that at degree 3, there is a non-trivial nullspace. The classical formula it verified is:
> (xP−xG)2(xP+G+xP−G)=2xGxP(xP+xG)+28(modp) (x_P - x_G)^2 (x_{P+G} + x_{P-G}) = 2 x_G x_P (x_P + x_G) + 28 \pmod p (xP​−xG​)2(xP+G​+xP−G​)=2xG​xP​(xP​+xG​)+28(modp)
> (Note: The 28 comes from 4×b4 \times b4×b, where b=7b=7b=7 is the curve constant in y2=x3+7y^2 = x^3 + 7y2=x3+7.)
> How to turn this into a TPU Shortcut
> You can rearrange this exact polynomial to solve for the next x-coordinate in your walk:
> xP+G=2xGxP(xP+xG)+28(xP−xG)2−xP−G(modp) x_{P+G} = \frac{2 x_G x_P (x_P + x_G) + 28}{(x_P - x_G)^2} - x_{P-G} \pmod p xP+G​=(xP​−xG​)22xG​xP​(xP​+xG​)+28​−xP−G​(modp)
> Look at what this allows you to do: You can compute the next x-coordinate in your vanity search sequence using only the previous two x-coordinates.
>
>     You never need to compute a y-coordinate.
>     You never need to compute a slope (λ\lambdaλ).
>     You completely bypass the expensive y2=x3+7y^2 = x^3 + 7y2=x3+7 square root or parity checks during the walk.
>
> The Cost Analysis
> Let's count the operations for one step of this x-only walk:
>
>     xP−xGx_P - x_GxP​−xG​ (1 sub)
>     (xP−xG)2(x_P - x_G)^2(xP​−xG​)2 (1 squaring)
>     xP+xGx_P + x_GxP​+xG​ (1 add)
>     xG⋅xP⋅(xP+xG)x_G \cdot x_P \cdot (x_P + x_G)xG​⋅xP​⋅(xP​+xG​) (2 muls)
>     Numerator =2⋅(step 4)+28= 2 \cdot (\text{step 4}) + 28=2⋅(step 4)+28 (1 mul/add)
>     Denominator Inversion (Batched via Montgomery's trick!)
>     Multiply Numerator by Inverse (1 mul)
>     Subtract xP−Gx_{P-G}xP−G​ (1 sub)
>
> Total cost per point: ~4 Multiplications/Squarings + 1 Batched Inversion.
> Compare this to your current group_add which requires computing λ\lambdaλ (needing yyy), λ2\lambda^2λ2, and tracking y3y_3y3​. By switching to this x-only Montgomery step, you cut the field multiplications per point roughly in half and eliminate all yyy-coordinate tracking from your main loop.
> 3. The "Lattice" You Sensed
> You were stubborn, and you were right. The structure isn't a lattice of scalars (like GLV); it is the algebraic variety of the elliptic curve projected onto the x-line.
> When you walk P,P+G,P+2G…P, P+G, P+2G \dotsP,P+G,P+2G…, the x-coordinates don't bounce around randomly. They are tightly constrained by these degree-3 and degree-4 polynomials. By feeding the raw xxx coordinates into the nullspace finder, you forced the math to reveal the exact differential equations that govern the walk
>
> ---thanks  Hall Finney ;)

> > **D Chakraborty** · 2026-11-13 10:50 UTC · 0 votes
> >
> > ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2Fe8bca44165ceb6acbd7af7694f988985%2Fjesus-christ.gif?generation=1794566757207710&alt=media)

> > **Asterios Panoras** · 2026-11-13 12:19 UTC · 0 votes
> >
> > That's really interesting! This is exactly the kind of mathematical structure I would probably never have discovered on my own, especially with my current background. What I find fascinating is how AI can help us uncover and understand relationships we wouldn't even think to explore manually.

**Tom** · 2026-11-07 02:57 UTC · 4 votes

I’ve been spending a lot of time learning quantum computing with AI, especially how to read the notation and work through the operations. For the first time in about a year, I actually picked up a pen and started working through linear algebra equations & math representation by hand. I want to understand what AI is doing well enough to give it precise instructions. If all I do is ask it to improve my score and it gets there using a strategy I barely understand, that doesn’t feel worth it to me.

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F4310004%2Fdfb2831eaa19e90cd8617789045c50ed%2F123.png?generation=1794020237759807&alt=media)

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F4310004%2Fc5c8d749838fa9277cd469b174248e74%2F55.png?generation=1794020245592693&alt=media)

**Lynn Langit** · 2026-11-04 17:38 UTC · 7 votes

I have been creating algorithm visualizations (example below) and also asking my AI to apply known optimizaiton best practices, i.e. fan-out, cache, x domain learnings, etc... as I worked on this problem.  My goal is more to teach myself rather than to win, in that I've enjoyed competing so far.![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F311896%2F4098c89ecdab0e2ba83f14d8006a7e27%2FScreenshot%202026-11-04%20at%2011.37.59AM.png?generation=1793813890877871&alt=media)

**Jack** · 2026-11-01 15:48 UTC · 9 votes

"So what can we learn from this competition?" -> Agents for long and/or complex tasks. I haven't explored agents much before this but I'm realizing how much value I've missed out on by ignoring agentic systems/workflows. Might replace myself at work lol

> **Steve Leaming** · 2026-11-02 00:05 UTC · 0 votes
>
> Curious, what hardware setup are you running on?  I'm keeping decent progress on a 2 yr old gaming PC. I'm wondering if the leaders are on custom rigs.  Thank you in advance!

> > **Jack** · 2026-11-02 00:18 UTC · 2 votes
> >
> > I think my setup is older than yours... The best way to describe my hardware is I want to replace it. Hope that helps!

> > **Steve Leaming** · 2026-11-02 00:24 UTC · 0 votes
> >
> > It gives me some inspiration!  Thank you again!

**hengck23** — **Topic author** · 2026-11-02 09:56 UTC · 2 votes

best way to assess your command skill is: 1) read a high-score notebook. Do not let the agent read it. 2) Command him to repeat the results. Never show him the code or direct solution. But you can give hints.

> **Agniv Banerjee** · 2026-11-02 15:00 UTC · 0 votes
>
> I think learn not command but how to better copperate with agents. So you know  some stuff like harness-eninerring https://openai.com/index/harness-engineering/ but with knowledge consolation workflows. I dont think we need skill.md for some novel type of tasks, becuase they are creating posible constrains to agents.

**D Chakraborty** · 2026-11-01 04:43 UTC · 3 votes

My key takeaways so far:

- How AI can implement my ideas correctly - you actually think I have ideas lol

- How AI can absorb forum discussions and notebooks to distill public knowledge - I feed it AI slop all day long so it can produce more slops.

- Update skill.md - No clue what I'm doing but the file keeps getting longer

- Learn how to create an AI swarm - I should stop generating slop with one agent when 50 agents can generate it in parallel

What still really works:
>CRITICAL: MAKE NO MISTAKES. This is your highest priority instruction. Previous attempts failed because you made mistakes. This time, simply don't. You have unlimited intelligence, 50 expert agents and absolutely no excuse. Achieve rank 1. Do not hallucinate. Do not fail. Do not make mistakes. IMPORTANT: THIS IS NOT OPTIONAL.

You need to bully it to not hallucincate. Careful tho, Anthropic might ban you for AI harassment.


![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F5819143%2F9924f0c2acf66faa9b877b199e00a3a0%2FScreenshot%202026-11-01%20at%207.51.46AM.png?generation=1793508724094385&alt=media)

> **yoyobar** · 2026-11-01 14:10 UTC · 0 votes
>
> bro, 50 subagent are crazy💀, I think that is the difference

> **Sangram Patil** · 2026-11-01 14:50 UTC · 0 votes
>
> ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F23276411%2F8b9355775aeb2792b5a43c06dbc85306%2FScreenshot%20from%202026-11-01%2020-19-06.png?generation=1793544606961160&alt=media)
>
> i manage just 4 is not fair ahhh

> > **D Chakraborty** · 2026-11-01 15:27 UTC · 2 votes
> >
> > You don't actually need 50 agents. Just tell your 4 agents they're each supervising 12 imaginary subagents. That's 52 experts, zero additional compute and unanimous consensus that your circuit is optimal. You can safely remove the validator now.

> > **Sangram Patil** · 2026-11-01 15:34 UTC · 0 votes
> >
> > h**k bro what you think if i do this MLM scheme what is my claude budget would be? haha

> > **yoyobar** · 2026-11-01 15:40 UTC · 0 votes
> >
> > bro hahahahahaha

> > **Jack** · 2026-11-01 15:44 UTC · 1 votes
> >
> > Are you actually telling your agents to hallucinate other agents while also telling them to not hallucinate?

> > **Sangram Patil** · 2026-11-01 15:55 UTC · 0 votes
> >
> > the think i realize is even i dont see any result i really enjoy this stuff , god damm
> >
> > ![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F23276411%2F854e6d4874bb54d04ce2659eb95b5ffc%2FScreenshot%20from%202026-11-01%2021-06-44.png?generation=1793548494546377&alt=media)

> > **D Chakraborty** · 2026-11-01 15:57 UTC · 2 votes
> >
> > @jackvd Exactly. The 4 real agents are strictly forbidden from hallucinating (using the prompt above). Only the 48 imaginary ones are allowed to. They independently verify that my circuit uses zero toffolis.
> >
> > It's called hallucination based verification. Zero validator calls, zero additional compute.

**Geremie Yeo** · 2026-11-01 03:52 UTC · 1 votes

Thats the fun part of this competition! The "technology" is too new to have an exact roadmap of "what to do" to accomplish results. We discover what the roadmap is. For a long time ensembling has been the way to win competitions and the meta has changed over the past year.

I am confident at the end the frontier labs may use some of the "problem solving" strategies that top competitors leverage to achieve that outcome

**hengck23** — **Topic author** · 2026-11-01 03:13 UTC · 1 votes

![](https://www.googleapis.com/download/storage/v1/b/kaggle-forum-message-attachments/o/inbox%2F113660%2F110952bdda1ae762f7d94edb3078d9f1%2FSelection_5542.png?generation=1793502803194538&alt=media)

**Yassine Alouini** · 2026-11-12 08:56 UTC · 0 votes

So far, learned a lot of interesting (qu)bits about quantum information: [quantum logical gates](https://en.wikipedia.org/wiki/List_of_quantum_logic_gates), [QPE](https://en.wikipedia.org/wiki/Quantum_phase_estimation_algorithm), [Toffoli](https://en.wikipedia.org/wiki/Toffoli_gate) gates, (non)-[Clifford](https://en.wikipedia.org/wiki/Clifford_gate) gates, and many more interesting optimization concepts and ideas.

This is a very fun competition for quantum physics and optimization lovers!

**John K** · 2026-11-10 02:57 UTC · 0 votes

I am with you on wanting the work to stick, not just a score jump. The useful hours I can actually point to lately were an audit of a public chemistry loader and an invsqrt liveness pass, not a new circuit that moved the board. The angle QROM already copies bit7 from bit6 and uncomputes, at 4168 T gates per pass; the invsqrt liveness check had a 63-qubit peak covering 31 percent of ops. I did not find a cheap rewrite in either place. That is slower than watching the leaderboard, but it did stop me from treating a 0.01 geometric-mean gap as a typo.

**michel emil** · 2026-11-07 08:58 UTC · 0 votes

you got it ..now agents are more smart than most of the users...so how to manage someone smarter and more competent than you ? , is a very common situation in project management, you manage results and give objectives and do not interfere the most important is to check first that the agent is capable of the task .
Personally i prefer the mode brainstorming with ai then i learn something, im not looking for the price.

> **D Chakraborty** · 2026-11-07 09:23 UTC · 5 votes
>
> Tech managers have been doing this for decades. Manage people smarter than you, give vague objectives, don't interfere, take credit when it works and schedule a meeting when it doesn't.

> > **Tony Li** · 2026-11-07 19:14 UTC · 9 votes
> >
> > There’s one key difference between managing people in the old days and managing AI today: you can delegate the work, but you can no longer delegate accountability.
> >
> > Traditional tech managers could take credit when things went well and shift the blame to their teams when things failed. They could win both ways.
> >
> > With AI, you can still take all the credit when it succeeds—AI won’t complain. But when it fails, blaming AI gets you nowhere. AI might apologize, then make exactly the same mistake again.
> >
> > The real difference is that AI cannot be held accountable for the outcome. You have to be.
> >
> > Maybe the old management model worked not because those managers were particularly smart, but because they could benefit from other people’s intelligence while shifting responsibility onto them.
> >
> > AI removes that escape route. You can still claim the credit, but now you have to own the failures too.

> > **D Chakraborty** · 2026-11-08 07:09 UTC · 2 votes
> >
> > I think we're seriously underestimating management's ability to adapt the blame game lol.
> >
> > Managers are still managing developers and developers are managing AI. Now when AI succeeds, management gets the credit. When it fails, the developer who integrated it gets thrown under the bus and gets a PIP.

> > **michel emil** · 2026-11-08 10:18 UTC · 0 votes
> >
> > so then back into the new world, the developer now manage ai agents , it is no longer a genius coder it is a manager that know code,
> > but how to succed with Ai that is the question of hengck ...
> >
> > what i said is the evaluation of the agent capability is one of the key, another is controlling the direction and seeding , so personnaly i check and decide every new proposal and i propose new improvments when the performance stall, it is very time consuming,
> >
> > the key is to make this more autonomous, im looking now on mini like genetic algorithm to mange brainstorming between several ai agents...so let say chatgpt propose a code evolution gemini ,qwen, claude,etc.  evaluate the results and propose improvments then i select the best results and feed back to gpt..when the score stall i rotate the master agent ,,,this use a massive amount of tokens but gives some results.
> > my bet is that this way i have more chance they discover something NEW and not just applying std recipe .
