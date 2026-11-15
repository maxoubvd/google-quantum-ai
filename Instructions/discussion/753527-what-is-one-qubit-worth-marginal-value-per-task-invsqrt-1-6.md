# What is one qubit worth? Marginal value per task (invsqrt 1.6%, DQI 0.23%, chem 0.23%)

> Kaggle discussion topic: https://www.kaggle.com/competitions/google-verifiable-quantum-circuit-compilation/discussion/753527 — copied 2026-11-15 (raw Markdown, unedited).

**Harshit Kumar** · 2026-10-28 12:50 UTC · 3 votes

The score is a geometric mean, so 1% off any task's cost gives +0.336% on the score. But qubits are priced very differently:
- invsqrt (cost = Q x C, Q=63): one qubit = 1.61% of cost -> +0.54% score
- DQI (cost = R x Q^2 x T, Q=889): one qubit = 0.23% of cost -> +0.08% score
- chemistry (cost = lambda x T x Q, Q=430): one qubit = 0.23% -> +0.08% score
Cost-neutral trade rates from the public 22.63 table: about 26 Toffoli-equivalents per qubit in invsqrt, about 990 Toffolis per qubit in DQI.

Notebook (arithmetic on published numbers, no circuits redistributed; credits to the public authors): https://www.kaggle.com/code/harshitsama/quantum-circuits-marginal-value-calculator

Leaderboard shape (1,403 teams, 2026-10-28): 689 teams >= 20, 451 >= 22, 272 >= 22.5, only 20 >= 23. The plateau is where public circuits meet. With ratios of about 20x, 59x and 10x, the chemistry ratio is the lowest, so a new chemistry shape has the most room; DQI qubits are cheap to ignore, invsqrt qubits are worth 7x more.

Tell me if any of my formulas are off.
