# Manuscript correction — Section V.B half-life (spec line 202)

**Spec text:** "t1/2 (half-life): ln(2) / k_elim = 8.7 hours (mRNA + LNP complex)" — stated for
a model the same section defines as two-compartment.

**Problem:** category error. `ln(2)/k_elim` is a micro-constant. It equals the *terminal*
half-life only in a one-compartment model. For a two-compartment system, the observable
half-lives come from the eigenvalues of the system matrix, not from k_elim alone.

**Correct values** (k12=0.30/h, k21=0.05/h, k_elim=0.08/h, all as given by the spec):
- alpha (distribution, fast) = 0.420487 /h → t_half = 1.648 h
- beta (terminal, slow) = 0.009513 /h → t_half = **72.86 h** — 8.4x longer than the spec's 8.7 h

Verified two ways: analytically (eigenvalues of the 2x2 system matrix) and empirically (linear
fit of ln(A1) over the solved curve's terminal segment) — both agree to 0.0000%.

**Affects Section V.C — magnitude corrected below.** Day 21 (504 h) is 6.92 terminal
half-lives. Solving the full two-exponential system exactly (not the single-exponential
approximation) gives a day-0 residual at day 21 of **0.686% of the prime dose** (A1 0.0082 μg +
A2 0.0604 μg = 0.0686 μg of 10 μg). Day 56 (1344 h, 18.4 terminal half-lives) residual is
**0.0002%** — negligible.

**Correction to an earlier draft of this note:** an initial estimate of ~0.83% used the
single-exponential approximation `exp(-beta*t)` alone, which ignores the two-exponential
system's actual amplitude split between A1 and A2. The exact figure is 0.686%, not ~0.83%.

**Conclusion, precisely stated:** the residual is sub-1% at day 21 and negligible at day 56 —
**not a substantial accumulation.** Step 5C must still show a non-zero trough and a
non-decreasing C_max across doses (a linear time-invariant system cannot produce a smaller peak
from an identical dose plus a positive residual), but the correct reason is that C_max *cannot
decrease*, not that it rises appreciably. The spec's table (lines 211-215, C_max 4.2→3.8→3.5,
troughs 0→0.15→0.08) is still internally impossible for the reason above, independent of how
large the residual turns out to be.

**Source:** `Step_Outputs/Phase5/StepB/Phase5B_Ledger_*.csv` keys `5B_alpha_per_h`,
`5B_beta_per_h`, `5B_t_half_distribution_h`, `5B_t_half_terminal_h`,
`5B_empirical_terminal_slope_per_h`.
