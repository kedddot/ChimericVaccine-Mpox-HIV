# Phase 5A -- Methodology Note

**Spec:** Section V.A, `REVISED_METHODOLOGY_MpoxHIV.md` lines 166-181.

**What this step does:** reverse-PK dose calculation. `Required_mRNA_Dose = Protein_Target / (T_eff * Bioavailability)`; `Injected_Dose = Required_Serum_Level * V_d`.

**Deviation from the spec:** the spec reports only the midpoint (T_eff=0.30, target=75 ng/mL, BA=0.80, V_d=5-10 L -> 1.56-3.1 ug). This script reproduces that midpoint exactly as a direction check, then sweeps the full stated ranges as a 9-step-per-axis, 4-axis grid (6561 points: T_eff 20-40%, target 50-100 ng/mL, bioavailability 80-90%, V_d 5-10 L) and reports min/median/max rather than only the midpoint (Rule 2: the methodology's numbers are targets, reproduced here, not final results).

**Correction (post-review): bioavailability is now swept, not fixed at 0.80.**

Section V.A's own worked example uses a single BA=0.80, reproduced exactly by the midpoint check above. Section III.C (line 68) separately gives bioavailability as **80-90%** of injected mRNA remaining functional at 2 h post-injection. That range, not the single V.A figure, is used for the sweep below.

**Provenance note:** Section III (new-scheme LNP delivery) has not been executed in this project -- it is out of scope for Section V and separately blocked (`gmx`/`orca`/`RNAfold` not installed). The 80-90% figure is therefore carried in as a **stated value from the spec text**, not something this pipeline derived or measured.

**Effect on the dose gap:** with BA fixed at 0.80 (prior version), the grid gave 0.781-6.250 ug (median 2.314), gap 1.60x-12.80x vs the proposed 10 ug (median 4.32x). With BA swept 0.80-0.90, the grid gives 0.694-6.250 ug (median 2.161), gap 1.60x-14.40x (median 4.63x). Higher BA lowers the required dose, so the fix widens the low end of the grid and raises the gap's upper bound; the high end of the grid (BA=0.80) is essentially unchanged, so the gap's lower bound barely moves.

**No escalation.** This is one of the brief's explicitly pre-identified non-triggers; the fix changes the size of the gap, not whether it exists or whether it needs escalating.
