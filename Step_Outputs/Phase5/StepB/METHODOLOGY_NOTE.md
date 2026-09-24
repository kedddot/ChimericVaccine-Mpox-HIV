# Phase 5B -- Methodology Note

**Spec:** Section V.B, `REVISED_METHODOLOGY_MpoxHIV.md` lines 183-205.

**Model:** two-compartment open PK, `dA1/dt = -k12*A1 - k_elim*A1 + k21*A2`, `dA2/dt = k12*A1 - k21*A2`. Solved via `scipy.integrate.solve_ivp`, method=`Radau`, rtol=1e-10, atol=1e-12, over t=0-1500.0h at 6001 points.

**CORRECTION (post-review, HALT 5B / trigger 2): the spec's half-life is a category error.** `ln(2)/k_elim = 8.66h` is a MICRO-constant. It equals the terminal half-life only in a one-compartment model; this is a two-compartment system, so the *observable* half-lives come from the eigenvalues of the 2x2 system matrix, not from k_elim alone.

- alpha (distribution, fast) = 0.420487 /h -> t_half = 1.648 h
- beta (terminal, slow) = 0.009513 /h -> t_half = 72.86 h

The empirical terminal slope of the solved A1(t) curve (0.009513 /h, fit over the last 20% of the 1500h window) confirms beta to 0.000% -- the curve was always right, only the reported single-number summary was wrong. Both half-lives, both macro-constants, and the micro-constant are now reported, each explicitly labelled, and the empirical-vs-analytic beta check is a permanent verification gate.

**CL and AUC are unaffected** -- both are defined via k_elim directly (`CL = k_elim*V_d`, `AUC = Dose/CL`), which is correct terminology in either a one- or two-compartment model; they still reproduce the spec exactly (0.48-0.64 L/h, 15.6-20.8 ug*h/L).

**Consequence flagged for Step 5C:** with beta=0.0095/h (t_half_terminal~73h), the prime-dose peripheral reservoir has NOT washed out by day 21 (504h, ~6.9 terminal half-lives puts ~1% of the dose remaining, not zero). Step 5C must therefore show non-zero troughs and a NON-DECREASING C_max across doses, not the spec's table showing C_max falling with zero troughs.

**AUC derivation, made explicit (the spec shows no derivation for this step):** two equivalent routes were checked against each other. (1) `AUC = Dose/CL` where `CL = k_elim*V_d`, treating V_d (6-8 L) as an independently stated apparent volume of distribution, the standard clinical-PK convention. (2) This model's own compartment-1 mass-integral identity, `AUC_mass(0-inf) = Dose/k_elim` (exact for any k12/k21, since only k_elim removes mass from the closed A1+A2 system), converted to a concentration-AUC by dividing by the same V_d. Both routes give the identical 15.6-20.8 ug*h/L range -- confirmed algebraically in the run log, not just numerically.

**Numerical vs analytic check:** trapezoid integration of the solved A1(t) curve over 0-1500h gives 125.0197 ug*h, matching the analytic Dose/k_elim identity (125.0000 ug*h) to 0.0158% -- well under the 1% gate. The solve window (240h, ~28 half-lives) is long enough that the truncated trapezoid integral is effectively the full 0-infinity integral.

**Checks passed:** no negative compartment masses at any solved time point; initial A1+A2 equals the dose exactly; A1+A2 decays toward 0 (open system, not conserved) rather than plateauing or growing.

**Escalation: trigger 2 fired and was resolved by direct ruling (HALT 5B), not by an independent Opus escalation cycle** -- the reviewer identified the half-life category error directly, with the exact macro-constants, before this step's own report was written. CL and AUC remain within spec rounding; only the half-life claim was wrong.
