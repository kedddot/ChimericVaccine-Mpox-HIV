# Step 5C Report — Multi-Dose Schedule (Section V.C, spec lines 207–220)

> **Superseded update (2026-09-22, later same day):** Finding 4 below (T_max=0h) triggered a
> follow-up revision, `STEP_5C_DEPOT_FIX_PROMPT.md`, adding the missing IM absorption (depot)
> compartment. Both model variants are now on disk, clearly labelled:
> `Phase5C_*_IVbolus_20260922_2013.csv` (the model as originally written, described below,
> unchanged) and `Phase5C_*_IMdepot_20260922_2013.csv` (the fix, ka swept 0.15/0.20/0.25/h per
> `STEP_5C_DEPOT_FIX_PROMPT.md`, all gates passed: T_max lands at 4.57/3.92/3.48h matching the
> prompt's expected table exactly, T_max>0 for every dose, mass balance holds, C_max invariant
> holds, cumulative AUC within 0.36–0.42% of the IV-bolus run). Full detail is in the appended
> "Depot (IM Absorption) Revision" section of `METHODOLOGY_NOTE.md`. Findings 1–3 below are
> unaffected (both models show the same invariant and the same spec contradiction). This
> original report is kept intact below as the record of what motivated the fix.

Run: `Phase5C_DoseMetrics_20260922_1955.csv`, `Phase5C_FullCurve_20260922_1955.csv`,
`Phase5C_SpecContradiction_20260922_1955.csv`, `Phase5C_Ledger_20260922_1955.csv`
(16 numbers). Rule 1: 1810/1810 OK before and after.

## Model

Same 2-compartment ODE as Step 5B (`k12=0.30/h, k21=0.05/h, k_elim=0.08/h`), integrated as
**one continuous run** over 0–84 days with an instantaneous +10 μg jump to A1 at days 0, 21, 56.
State carries over exactly across dose boundaries — this is what lets a prior dose's residual
reach the next dose's peak, which three independent single-dose runs could not show.

## Results

| Dose | Day | C_max (μg/L) | T_max (h post-dose) | AUC, this interval | C_min (trough before next dose) | Cumulative AUC |
|---|---|---|---|---|---|---|
| 1 | 0 | 1.42857 | 0.0 | 17.780 | 0.00116 | 17.780 |
| 2 | 21 | 1.42975 | 0.0 | 18.021 | 0.00005 | 35.801 |
| 3 | 56 | 1.42862 | 0.0 | 17.883 | 0.00024 | 53.684 |

(Concentration = A1/V_d, V_d=7 L, the midpoint of Step 5B's 6–8 L range.)

## Finding 1 — the correct invariant, and a self-correction

**Every dose's C_max is ≥ the first dose's C_max.** Pre-dose residual is always ≥0 in this model
(no negative compartments, no absorption loss route other than k_elim), so C_max at dose *i*
= (Dose + residual_i)/V_d ≥ Dose/V_d = C_max at dose 1. Verified as a hard gate: 1.42857 ≤ all
three values.

**This does *not* mean C_max rises pairwise from dose to dose.** My first implementation
asserted exactly that (stronger claim) and the run failed its own gate — dose 3's C_max
(1.42862) is *below* dose 2's (1.42975). That is legitimate model behaviour, not a bug: this
schedule's inter-dose gaps grow (21 days, then 35 days), so more washout time elapses before
dose 3 than before dose 2, and the pre-dose-3 residual (0.0028 μg total, both compartments) is
smaller than the pre-dose-2 residual (0.0686 μg) despite dose 3 having two prior doses behind it.
The gate was corrected to the true unconditional invariant (≥ dose-1 baseline) rather than the
stronger, interval-dependent pairwise claim. Both are ≥ the dose-1 baseline, satisfying the real
requirement.

## Finding 2 — the spec's table is still wrong, restated precisely

Spec (lines 211–215): C_max 4.2 → 3.8 → 3.5 μg/L, C_min 0 → 0.15 → 0.08 μg/L.

The decisive, unconditional problem: **the spec's own dose-3 C_max (3.5) is lower than its own
dose-1 C_max (4.2).** That specific relationship is impossible regardless of interval spacing,
because pre-dose-3 residual cannot be negative. (A dose-2-to-dose-3 dip on its own is *not*
automatically impossible, per Finding 1 — our model reproduces a small one honestly.) We did not
reverse-engineer parameters to reproduce the spec's table; our own values are reported and the
contradiction is logged (`Phase5C_SpecContradiction_20260922_1955.csv`).

## Finding 3 — magnitude, corrected from Step 5B's HALT

Per the Step 5B correction (exact two-exponential solution, not the single-exponential
approximation): the day-21 trough is 0.00116 μg/L against a 1.42857 μg/L peak —
**0.081% of the dose-1 peak.** Total C_max growth from dose 1 to dose 3 is **0.0035%.** This is
real, non-zero, and directionally consistent with a positive residual, but it is **not a
substantial accumulation** — do not report it as one.

## Finding 4 — new, not yet flagged anywhere: T_max is always 0 in this model

**T_max = 0 h post-dose for all three doses.** This is a direct consequence of how the spec's
own two-compartment equations are written: the 10 μg dose is added directly to the central
compartment (A1) with no absorption/depot compartment or first-order uptake term, so
concentration is highest at the instant of dosing and decreases monotonically afterward
(`dA1/dt|_{t=dose+} = -(k12+k_elim)·A1 + k21·A2`, and the first term dominates whenever A1 is
freshly large, which it always is immediately post-dose).

This contradicts two things the spec itself states: Section III.C's "Time-to-Peak Concentration
(T_max): 4-8 h", and Section V.C's own table, which lists T_max as 6h, 5h, 4h for doses 1-3.

**Reading:** the spec's narrative describes absorption/distribution behaviour (a delay to peak)
that its own stated equations do not produce. The equations as given model an instantaneous IV
bolus directly into the central compartment; a nonzero T_max requires an explicit absorption
term (e.g. a depot/injection-site compartment with its own first-order rate into A1), which
Section V.B/V.C never defines. This is logged as a finding, not silently patched by adding an
unstated absorption compartment.

## No escalation

This step's own pre-identified non-trigger (spec's C_max-decreasing table) is confirmed, restated
more precisely per Finding 1/2 above. Finding 4 (T_max=0) is new and not one of the brief's
pre-identified items; it does not meet any of the 5 escalation triggers (no verification gate
failed once the gate itself was corrected; no >2x contradiction on a headline PK figure — T_max
was never gated; no missing/underivable parameter; no Phase I/II/IV edit needed; not part of the
5D alpha-scenario trigger). Logged for the manuscript-corrections file, not escalated.

## Files

- `Phase5C_DoseMetrics_20260922_1955.csv` — per-dose C_max/T_max/AUC/C_min/cumulative AUC
- `Phase5C_FullCurve_20260922_1955.csv` — full 0-84 day concentration curve, hourly, with dose-event jumps
- `Phase5C_SpecContradiction_20260922_1955.csv` — side-by-side spec vs. model values
- `Phase5C_Ledger_20260922_1955.csv` — 16 numbers, each traced to `Phase5C_MultiDose`
- `METHODOLOGY_NOTE.md` — solver settings, model definition, deviations

Holding here per your HALT — not starting Step 5D.