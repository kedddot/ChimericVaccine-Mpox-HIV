# Phase 5D -- Methodology Note

**Spec:** Section V.D, `REVISED_METHODOLOGY_MpoxHIV.md` lines 222-258.

**Ag(t) source:** `Phase5C_FullCurve_IMdepot_*.csv` (never the IV-bolus curve -- the depot fix's whole purpose was a realistic ~4h absorption delay, and an instantaneous-spike-driven T-cell model would give a materially different early response). Concentration_ug_per_L column used as Ag(t), linearly interpolated (`numpy.interp`) at arbitrary t. Primary ka=0.2/h; sensitivity ka=(0.15, 0.25).

**Solver:** `scipy.integrate.solve_ivp`, method=`LSODA`, rtol=1e-10, atol=1e-12. Radau (used in Steps 5B/5C) struggled with this interpolated-forcing right-hand side, whose derivative has kinks at each dose event and at the underlying curve's own hourly sample points; LSODA (variable-order, automatically switches stiff/non-stiff) handled it robustly. This is a solver-choice deviation from Steps 5B/C, made because the RHS itself changed character (interpolated external forcing vs a pure closed-form linear system), not a loosening of accuracy standards -- rtol/atol are unchanged from 5B/5C.

**State vector:** [T_eff, G1, GC_output, B_mem, Plasma_cell, Ab_titer], all arbitrary units, all start at 0.

**T-cell equation, spec ambiguity resolved:** the spec writes `dT_eff/dt = a*Ag(t) - b*T_eff - g*T_eff*T_reg` where g's own stated units already carry a `(T_reg/T_total)` factor, and no T_reg(t) trajectory is defined anywhere in the spec. Rather than invent an unstated T_reg(t) state, b and g were folded into a single effective per-day suppression rate, `-(b+g)*T_eff`; both are still swept across their full stated ranges (b: 0.01-0.02/day, g: 0.001-0.005/day), as a 3-point (low/mid/high) grid.

**GC_output coupling and B-cell/antibody parameters:** see manuscript-corrections items 2-3 above (same file, `MANUSCRIPT_CORRECTIONS_5D.md`) for the full reasoning; not repeated here.

**Alpha scenarios:** three, per the brief, none picked as a winner. Scenario 2 ("BigMHC-weighted") uses the construct-wide BigMHC-IM mean (0.3502) uniformly (tests whether weighting by overall immunogenicity strength vs. an arbitrary uniform 1.0 changes T-cell kinetics shape). Scenario 3 ("per-pathogen") uses BigMHC-IM's own HIV/Mpox means (0.382/0.318) -- BigMHC-IM was chosen over PRIME for this split because it is the TCR ladder's primary rung and a direct immunogenicity probability, not because it gave an expected-looking answer (the PRIME-based cross-check, which disagrees in direction, is reported plainly, not hidden -- see manuscript corrections item 4).

**Direction checks (re-cited from Phase IV Step 3, not re-derived):** BigMHC-IM best score vs Phase4B EL rank, Pearson r=-0.2984, p=0.4024, n=10 (expected negative sign, correct, not significant at n=10). PRIME %Rank vs Phase4B EL rank, Pearson r=+0.9159, p=0.0002, n=10 (expected positive sign, meaningful). Both tools' scores were used only after confirming these signs match their stated conventions.

**Gates, all passed:** no negative state values anywhere; antigen peak <= T_eff peak <= antibody peak in every one of the 54 (ka, scenario, pathogen, b/g) combinations simulated; GC_output's prime-induced peak (before boost 1, day 21) falls in [11.38, 11.50] days, inside the spec's stated 10-14 day window for every combination; T_eff and Ab_titer both decay well below their peak values by end of study (no divergence).

**Determinism:** no stochastic elements anywhere in this step (deterministic ODEs, deterministic solver, no RNG seed needed). Re-running reproduces output exactly.

**Figures:** written (Phase5D_Scenario3_Curves_20260922_2102.png).
