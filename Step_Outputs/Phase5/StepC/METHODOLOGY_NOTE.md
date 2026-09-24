# Phase 5C -- Methodology Note

**Spec:** Section V.C, `REVISED_METHODOLOGY_MpoxHIV.md` lines 207-220.

**Model:** same 2-compartment ODE as Step 5B, integrated as ONE continuous run 0-84 days with instantaneous +10ug dose events at days 0, 21, 56 -- not three independent single-dose runs. State (A1, A2) carries over exactly across each dose boundary. Solver: `scipy.integrate.solve_ivp`, method=`Radau`, rtol=1e-10, atol=1e-12, hourly output resolution (24/day). Concentration = A1/V_d with V_d=7L (midpoint of Step 5B's 6-8L range).

**Consistent with the Step 5B correction, with the invariant stated precisely:** pre-dose residual is always >=0 in this model, so EVERY dose's C_max is >= dose 1's C_max (1.42857 ug/L) -- an unconditional invariant, verified as a gate here. This does NOT require C_max to rise pairwise from one dose to the next: this schedule's gaps grow (21 days, then 35 days), so more washout occurs before dose 3 than before dose 2, and our own model shows dose 3's C_max (1.42862) dip slightly below dose 2's (1.42975) while staying above dose 1's -- legitimate model behaviour, not a bug. The residual is SUB-1% of a single dose throughout (0.081% at the day-21 trough) -- real but small, not a substantial accumulation.

**EXPECTED DISAGREEMENT with the spec's table (lines 211-215), logged not reproduced:** the spec's dose-3 C_max (3.5) is LOWER than its own dose-1 C_max (4.2) -- impossible under the unconditional invariant above, regardless of interval spacing or residual size, since pre-dose-3 residual cannot be negative. A same-magnitude or slightly lower dose-2-to-dose-3 step is not on its own impossible (we reproduce a small instance of it ourselves, for the reason above); the decisive problem is falling below the dose-1 baseline. We did not reverse-engineer parameters to reproduce the spec's table; our own values are reported and the contradiction is logged in `Phase5C_SpecContradiction_IVbolus_20260922_2101.csv`.

**No escalation.** This is one of the brief's explicitly pre-identified non-triggers ("The C_max-decreasing table (5C). Already identified as internally inconsistent.").

---

# Phase 5C -- Depot (IM Absorption) Revision

**Source:** `STEP_5C_DEPOT_FIX_PROMPT.md`. The IV-bolus model above is the faithful reading of Section V.B as literally written; it gives T_max=0h for every dose, which is correct FOR THAT MODEL but contradicts Section III.C's stated T_max of 4-8h (line 69) and Section V.C's own table (6/5/4h) -- because Section V.B never defines an absorption route for the IM administration Section V.A specifies (line 179).

**Fix:** added a depot compartment using `k_IM` from Section III.B line 44 (`CL_depot = k_IM * LNP_dose`, k_IM=0.15-0.25/h): `dD/dt=-ka*D`, `dA1/dt=ka*D-k12*A1-k_elim*A1+k21*A2`, `dA2/dt=k12*A1-k21*A2`, D(0)=Dose, A1(0)=A2(0)=0. A 4th state `E` (cumulative eliminated mass, `dE/dt=k_elim*A1`) was added purely so mass balance is an independent numerical check, not an assumed algebraic identity.

**Validated against the fix prompt's expected single-dose table before trusting the full schedule** (all three ka, T_max within 0.1h, C_max/dose within 0.002):

- ka=0.15/h: T_max=4.5700h (expected 4.57h), C_max/dose=0.2264 (expected 0.2264)
- ka=0.2/h: T_max=3.9200h (expected 3.92h), C_max/dose=0.2685 (expected 0.2685)
- ka=0.25/h: T_max=3.4800h (expected 3.48h), C_max/dose=0.3040 (expected 0.304)

**All gates passed for all three ka values:** T_max>0 for every dose; depot empties by end of study; mass balance (D+A1+A2+E = administered dose) holds to <1e-4 ug; no negative compartments; every dose's C_max >= dose 1's C_max (the corrected invariant, not pairwise monotonicity, per the fix prompt's explicit correction to the original Opus guidance); cumulative AUC within 1% of the IV-bolus run (bioavailability out of the depot = 1.0, so absorption reshapes the curve without changing its area).

**Manuscript corrections, 4 items (1-2 new, 3-4 carried from the IV-bolus run):**

1. Section V.B models an IV bolus, not the IM route the study specifies (T_max=0h vs the spec's own stated 4-8h and table values); adding the depot with k_IM from Section III.B line 44 yields T_max 3.48-4.57h, consistent with the stated range. Recommend the model be respecified with three compartments.
2. The spec's C_max (4.2 ug/L at 10ug dose) implies V_d=2.38L, inconsistent with Section V.B line 203's stated V_d=6-8L (which gives 1.25-1.67 ug/L, IV-bolus; our depot value at V_d=7L spans 0.3235-0.4343 ug/L across ka). A ~3x inconsistency internal to the spec.
3. The spec's dose-3 C_max (3.5) is unconditionally impossible: it is below the spec's own dose-1 C_max (4.2), and pre-dose-3 residual cannot be negative in either model, regardless of dose spacing.
4. Accumulation is negligible in both models, not substantial -- do not describe this schedule as producing meaningful accumulation.

**Both model variants are kept and clearly labelled** (`*_IVbolus_*`, `*_IMdepot_*`) -- the manuscript needs both to show why the change was made.
