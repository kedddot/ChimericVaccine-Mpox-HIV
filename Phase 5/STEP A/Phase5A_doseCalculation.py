import os
import sys
import csv
import statistics

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_COMMON = os.path.normpath(os.path.join(_THIS_DIR, "..", "_common"))
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

import phase5_common as common

# =============================================================================
# STEP 5A -- mRNA dose calculation & translation efficiency (spec lines 166-181)
#
# Reverse-PK chain: given a target serum protein level and assumptions about
# translation efficiency and bioavailability, back-calculate the required
# mRNA dose. The spec gives a single midpoint example; this script sweeps
# the full stated ranges and reports the grid, not just the midpoint.
# =============================================================================

PARAMS = {
    "T_EFF_RANGE": (0.20, 0.40),          # translation efficiency, fraction (Section V.A)
    "PROTEIN_TARGET_RANGE_NG_ML": (50.0, 100.0),   # serum therapeutic window (Section V.A)
    "BIOAVAILABILITY_MIDPOINT": 0.80,      # value the spec's own V.A worked example uses
    "BIOAVAILABILITY_RANGE": (0.80, 0.90), # Section III.C line 68 -- III.C was NOT executed;
                                            # this is a stated literature value, not a derived one
    "VD_RANGE_L": (5.0, 10.0),             # volume of distribution, adult IM (Section V.A)
    "GRID_STEPS": 9,                       # steps per swept axis (0,12.5,...,100%)
    "PROPOSED_DOSE_UG": 10.0,              # the spec's Group-2 regimen dose
}


def _grid(lo, hi, n):
    if n == 1:
        return [lo]
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def midpoint_reproduction():
    """Reproduce the spec's own worked example exactly, as a direction check
    before sweeping. Spec (Section V.A): 75 ng/mL / (0.30*0.80) = 312.5 ng/mL
    serum; x V_d(5-10L) = 1.56-3.1 ug. The spec's own V.A example holds BA at
    a single fixed 0.80, even though Section III.C separately states BA as an
    80-90% range -- this check uses the spec's own fixed value so the
    midpoint reproduces exactly, before the sweep below applies the range."""
    t_eff = 0.30
    protein_target = 75.0
    bioavail = PARAMS["BIOAVAILABILITY_MIDPOINT"]
    required_serum_ng_ml = protein_target / (t_eff * bioavail)
    injected_lo = required_serum_ng_ml * 5.0 / 1000.0   # ng/mL * L -> ng -> ug (/1000)
    injected_hi = required_serum_ng_ml * 10.0 / 1000.0
    return required_serum_ng_ml, injected_lo, injected_hi


def sweep_fixed_ba_for_comparison():
    """Re-run the ORIGINAL 3-axis sweep (BA fixed at the V.A midpoint's 0.80)
    purely so this run can report how the dose gap changes once BA is
    correctly swept as a 4th axis. Not written as the step's own grid output
    -- comparison only."""
    t_effs = _grid(*PARAMS["T_EFF_RANGE"], PARAMS["GRID_STEPS"])
    targets = _grid(*PARAMS["PROTEIN_TARGET_RANGE_NG_ML"], PARAMS["GRID_STEPS"])
    vds = _grid(*PARAMS["VD_RANGE_L"], PARAMS["GRID_STEPS"])
    bioavail = PARAMS["BIOAVAILABILITY_MIDPOINT"]
    doses = []
    for t_eff in t_effs:
        for target in targets:
            required_serum = target / (t_eff * bioavail)
            for vd in vds:
                doses.append(required_serum * vd / 1000.0)
    return min(doses), statistics.median(doses), max(doses)


def sweep():
    """Four-axis sweep: T_eff x protein target x V_d x bioavailability. BA is
    included per correction -- Section III.C (line 68) states BA as 80-90%,
    not the single 0.80 the V.A worked example implies. III.C was never
    executed in this project, so this range is a stated literature value
    carried in from the spec text, not something derived by this pipeline."""
    t_effs = _grid(*PARAMS["T_EFF_RANGE"], PARAMS["GRID_STEPS"])
    targets = _grid(*PARAMS["PROTEIN_TARGET_RANGE_NG_ML"], PARAMS["GRID_STEPS"])
    vds = _grid(*PARAMS["VD_RANGE_L"], PARAMS["GRID_STEPS"])
    bas = _grid(*PARAMS["BIOAVAILABILITY_RANGE"], PARAMS["GRID_STEPS"])

    rows = []
    for t_eff in t_effs:
        for target in targets:
            for bioavail in bas:
                required_serum = target / (t_eff * bioavail)   # ng/mL
                for vd in vds:
                    injected_ug = required_serum * vd / 1000.0   # ng/mL * L = ng; /1000 -> ug
                    rows.append({
                        "T_eff": round(t_eff, 4),
                        "Protein_Target_ng_mL": round(target, 3),
                        "Bioavailability": round(bioavail, 4),
                        "V_d_L": round(vd, 3),
                        "Required_Serum_ng_mL": round(required_serum, 4),
                        "Injected_Dose_ug": round(injected_ug, 4),
                    })
    return rows


def run():
    common.print_banner("PHASE 5A -- mRNA DOSE CALCULATION (Section V.A, spec lines 166-181)")

    ok_pre = common.verify_rule1("Step 5A pre-check")
    if not ok_pre:
        print("[ERROR] Rule 1 dirty before Step 5A even ran -- refusing to proceed.")
        sys.exit(2)

    ledger = common.Ledger()

    print("\n" + "-" * 90)
    print("DIRECTION CHECK -- reproduce the spec's own midpoint worked example first")
    req_serum, inj_lo, inj_hi = midpoint_reproduction()
    print(f"  Required serum level: {req_serum:.1f} ng/mL (spec: 312.5 ng/mL)")
    print(f"  Injected dose range (V_d 5-10 L): {inj_lo:.2f}-{inj_hi:.2f} ug (spec: 1.56-3.1 ug)")
    assert abs(req_serum - 312.5) < 0.1, "Midpoint reproduction failed -- required serum level does not match spec"
    assert abs(inj_lo - 1.5625) < 0.01 and abs(inj_hi - 3.125) < 0.01, \
        "Midpoint reproduction failed -- injected dose range does not match spec"
    print("  [SUCCESS] Midpoint reproduces the spec exactly. Formula is correct; proceeding to the full sweep.")
    ledger.N("5A_midpoint_required_serum_ng_mL", round(req_serum, 2), "Phase5A (formula check)")
    ledger.N("5A_midpoint_injected_dose_lo_ug", round(inj_lo, 3), "Phase5A (formula check)")
    ledger.N("5A_midpoint_injected_dose_hi_ug", round(inj_hi, 3), "Phase5A (formula check)")

    print("\n" + "-" * 90)
    print(f"SWEEP -- T_eff {PARAMS['T_EFF_RANGE']}, protein target {PARAMS['PROTEIN_TARGET_RANGE_NG_ML']} ng/mL, "
          f"bioavailability {PARAMS['BIOAVAILABILITY_RANGE']} (Section III.C, not executed -- literature value), "
          f"V_d {PARAMS['VD_RANGE_L']} L, {PARAMS['GRID_STEPS']} steps/axis "
          f"({PARAMS['GRID_STEPS']**4} grid points, 4 axes)")
    rows = sweep()
    doses = [r["Injected_Dose_ug"] for r in rows]
    dose_min, dose_max = min(doses), max(doses)
    dose_median = statistics.median(doses)
    print(f"  Injected dose across full grid: min {dose_min:.3f} ug, median {dose_median:.3f} ug, max {dose_max:.3f} ug")
    ledger.N("5A_grid_dose_min_ug", round(dose_min, 3), "Phase5A_DoseGrid")
    ledger.N("5A_grid_dose_median_ug", round(dose_median, 3), "Phase5A_DoseGrid")
    ledger.N("5A_grid_dose_max_ug", round(dose_max, 3), "Phase5A_DoseGrid")
    ledger.N("5A_bioavailability_range_lo", PARAMS["BIOAVAILABILITY_RANGE"][0], "Section III.C line 68 (not executed)")
    ledger.N("5A_bioavailability_range_hi", PARAMS["BIOAVAILABILITY_RANGE"][1], "Section III.C line 68 (not executed)")

    out_dir = common.step_output_dir("A")
    ts = common.timestamp()
    grid_path = os.path.join(out_dir, f"Phase5A_DoseGrid_{ts}.csv")
    with open(grid_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Grid written: {grid_path} ({len(rows)} rows)")

    print("\n" + "-" * 90)
    print("FLAGGED FOR THE REPORT -- gap between the reverse-PK dose and the proposed regimen dose")
    proposed = PARAMS["PROPOSED_DOSE_UG"]
    ratio_lo = proposed / dose_max     # smallest multiple (using the largest reverse-PK estimate)
    ratio_hi = proposed / dose_min     # largest multiple (using the smallest reverse-PK estimate)
    ratio_mid = proposed / dose_median
    print(f"  Reverse-PK dose range (4-axis grid, BA 0.80-0.90): {dose_min:.2f}-{dose_max:.2f} ug; spec midpoint 1.56-3.1 ug")
    print(f"  Proposed regimen dose: {proposed:.1f} ug")
    print(f"  Proposed / reverse-PK ratio: {ratio_lo:.2f}x (vs grid max) to {ratio_hi:.2f}x (vs grid min); "
          f"{ratio_mid:.2f}x vs grid median")

    old_min, old_median, old_max = sweep_fixed_ba_for_comparison()
    old_ratio_lo = proposed / old_max
    old_ratio_hi = proposed / old_min
    old_ratio_mid = proposed / old_median
    print(f"\n  How the gap changes once BA is correctly swept (was fixed at 0.80, now 0.80-0.90, Section III.C):")
    print(f"    BA fixed at 0.80 (previous run): dose {old_min:.3f}-{old_max:.3f} ug (median {old_median:.3f}); "
          f"gap {old_ratio_lo:.2f}x-{old_ratio_hi:.2f}x (median {old_ratio_mid:.2f}x)")
    print(f"    BA swept 0.80-0.90 (this run):   dose {dose_min:.3f}-{dose_max:.3f} ug (median {dose_median:.3f}); "
          f"gap {ratio_lo:.2f}x-{ratio_hi:.2f}x (median {ratio_mid:.2f}x)")
    print(f"    Effect of the BA fix: widens the grid's low end (higher BA -> lower required dose), so the gap's "
          f"UPPER bound grows from {old_ratio_hi:.2f}x to {ratio_hi:.2f}x; the median gap moves from "
          f"{old_ratio_mid:.2f}x to {ratio_mid:.2f}x. The grid's high end (low BA=0.80, same as before) is "
          f"essentially unchanged, so the LOWER bound of the gap barely moves ({old_ratio_lo:.2f}x -> {ratio_lo:.2f}x).")
    print("  Still NOT flagged as an error (immunogenicity, not serum protein concentration, typically drives "
          "mRNA vaccine dosing) -- it is an unexplained gap between two adjacent paragraphs of the spec, reported "
          "here rather than silently adopting either number.")
    ledger.N("5A_proposed_dose_ug", proposed, "Section V spec, Group 2 regimen")
    ledger.N("5A_dose_gap_ratio_vs_grid_min", round(ratio_hi, 3), "Phase5A_DoseGrid")
    ledger.N("5A_dose_gap_ratio_vs_grid_max", round(ratio_lo, 3), "Phase5A_DoseGrid")
    ledger.N("5A_dose_gap_ratio_vs_grid_median", round(ratio_mid, 3), "Phase5A_DoseGrid")
    ledger.N("5A_dose_gap_ratio_median_BA_fixed_080", round(old_ratio_mid, 3), "Phase5A (comparison, not written to grid)")
    ledger.N("5A_dose_gap_ratio_median_BA_swept", round(ratio_mid, 3), "Phase5A_DoseGrid")

    summary_path = os.path.join(out_dir, f"Phase5A_Summary_{ts}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Metric", "Value", "Unit"])
        w.writerow(["Midpoint_required_serum", round(req_serum, 2), "ng/mL"])
        w.writerow(["Midpoint_injected_dose_lo", round(inj_lo, 3), "ug"])
        w.writerow(["Midpoint_injected_dose_hi", round(inj_hi, 3), "ug"])
        w.writerow(["Grid_dose_min", round(dose_min, 3), "ug"])
        w.writerow(["Grid_dose_median", round(dose_median, 3), "ug"])
        w.writerow(["Grid_dose_max", round(dose_max, 3), "ug"])
        w.writerow(["Proposed_regimen_dose", proposed, "ug"])
        w.writerow(["Gap_ratio_vs_grid_min", round(ratio_hi, 3), "x"])
        w.writerow(["Gap_ratio_vs_grid_max", round(ratio_lo, 3), "x"])
        w.writerow(["Gap_ratio_vs_grid_median", round(ratio_mid, 3), "x"])
        w.writerow(["Bioavailability_range_lo", PARAMS["BIOAVAILABILITY_RANGE"][0], "fraction"])
        w.writerow(["Bioavailability_range_hi", PARAMS["BIOAVAILABILITY_RANGE"][1], "fraction"])
        w.writerow(["Gap_ratio_median_BA_fixed_0.80_superseded", round(old_ratio_mid, 3), "x"])
        w.writerow(["Gap_ratio_median_BA_swept_0.80-0.90", round(ratio_mid, 3), "x"])
    print(f"  Summary written: {summary_path}")

    ledger_path = os.path.join(out_dir, f"Phase5A_Ledger_{ts}.csv")
    ledger.write(ledger_path)
    print(f"  Ledger written: {ledger_path} ({len(ledger.rows)} numbers)")

    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write("# Phase 5A -- Methodology Note\n\n")
        f.write("**Spec:** Section V.A, `REVISED_METHODOLOGY_MpoxHIV.md` lines 166-181.\n\n")
        f.write("**What this step does:** reverse-PK dose calculation. `Required_mRNA_Dose = "
                "Protein_Target / (T_eff * Bioavailability)`; `Injected_Dose = Required_Serum_Level * V_d`.\n\n")
        f.write("**Deviation from the spec:** the spec reports only the midpoint (T_eff=0.30, target=75 ng/mL, "
                "BA=0.80, V_d=5-10 L -> 1.56-3.1 ug). This script reproduces that midpoint exactly as a direction "
                f"check, then sweeps the full stated ranges as a {PARAMS['GRID_STEPS']}-step-per-axis, "
                f"{4}-axis grid ({PARAMS['GRID_STEPS']**4} points: T_eff 20-40%, target 50-100 ng/mL, "
                f"bioavailability 80-90%, V_d 5-10 L) and reports min/median/max rather than only the midpoint "
                "(Rule 2: the methodology's numbers are targets, reproduced here, not final results).\n\n")
        f.write("**Correction (post-review): bioavailability is now swept, not fixed at 0.80.**\n\n")
        f.write("Section V.A's own worked example uses a single BA=0.80, reproduced exactly by the midpoint "
                "check above. Section III.C (line 68) separately gives bioavailability as **80-90%** of "
                "injected mRNA remaining functional at 2 h post-injection. That range, not the single V.A "
                "figure, is used for the sweep below.\n\n")
        f.write("**Provenance note:** Section III (new-scheme LNP delivery) has not been executed in this "
                "project -- it is out of scope for Section V and separately blocked (`gmx`/`orca`/`RNAfold` "
                "not installed). The 80-90% figure is therefore carried in as a **stated value from the spec "
                "text**, not something this pipeline derived or measured.\n\n")
        f.write(f"**Effect on the dose gap:** with BA fixed at 0.80 (prior version), the grid gave "
                f"{old_min:.3f}-{old_max:.3f} ug (median {old_median:.3f}), gap {old_ratio_lo:.2f}x-{old_ratio_hi:.2f}x "
                f"vs the proposed 10 ug (median {old_ratio_mid:.2f}x). With BA swept 0.80-0.90, the grid gives "
                f"{dose_min:.3f}-{dose_max:.3f} ug (median {dose_median:.3f}), gap {ratio_lo:.2f}x-{ratio_hi:.2f}x "
                f"(median {ratio_mid:.2f}x). Higher BA lowers the required dose, so the fix widens the low end of "
                f"the grid and raises the gap's upper bound; the high end of the grid (BA=0.80) is essentially "
                f"unchanged, so the gap's lower bound barely moves.\n\n")
        f.write("**No escalation.** This is one of the brief's explicitly pre-identified non-triggers; the fix "
                "changes the size of the gap, not whether it exists or whether it needs escalating.\n")
    print(f"  Methodology note written: {note_path}")

    ok_post = common.verify_rule1("Step 5A post-check")
    if not ok_post:
        print("[ERROR] Rule 1 violated during Step 5A.")
        sys.exit(2)

    print("\n[SUCCESS] Step 5A complete.")
    return ledger.rows


if __name__ == "__main__":
    run()
