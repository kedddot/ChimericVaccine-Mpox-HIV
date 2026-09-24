import os
import sys
import csv

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_COMMON = os.path.normpath(os.path.join(_THIS_DIR, "..", "_common"))
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

import phase5_common as common
import numpy as np
from scipy.integrate import solve_ivp

# =============================================================================
# STEP 5C -- multi-dose schedule (spec lines 207-220)
#
# Days 0, 21, 56 at 10 ug each. Simulated as ONE continuous integration with
# dose events (state jumps at each dose time), not three independent runs --
# this is the only way a prior dose's residual can carry into the next peak.
#
# Per Step 5B (HALT, trigger 2): k_elim alone is NOT the terminal decay rate
# for this 2-compartment system. beta=0.009513/h (t_half_terminal=72.86h)
# governs long-run washout. Day 21 is 6.92 terminal half-lives after the
# prime -- residual 0.686% of dose (exact 2-exponential solution, corrected
# from an earlier ~0.83% single-exponential approximation).
#
# Conclusion (corrected from an initially over-strong version): pre-dose
# residual is always >=0 in this model, so EVERY dose's C_max is >= the
# first dose's C_max -- an unconditional invariant. This does NOT mean C_max
# rises pairwise from dose to dose: this schedule's intervals grow (21 days,
# then 35 days), so more washout occurs before dose 3 than before dose 2,
# and dose 3's C_max can legitimately dip slightly below dose 2's while
# staying above dose 1's. Either way the effect size is SUB-1% -- not a
# substantial accumulation, do not overstate it.
# =============================================================================

PARAMS = {
    "K12_PER_H": 0.30,
    "K21_PER_H": 0.05,
    "K_ELIM_PER_H": 0.08,
    "DOSE_UG": 10.0,
    "DOSE_DAYS": (0, 21, 56),
    "STUDY_END_DAY": 84,      # spec's "0-84 day course"
    "N_EVAL_PER_DAY": 24,     # hourly resolution
    "RTOL": 1e-10,
    "ATOL": 1e-12,
    "METHOD": "Radau",
    # --- depot (IM absorption) revision, added per STEP_5C_DEPOT_FIX_PROMPT.md ---
    # ka = k_IM, Section III.B line 44: CL_depot = k_IM * LNP_dose, k_IM = 0.15-0.25 /h.
    # The IV-bolus model above is kept and run unchanged (it is the faithful reading of
    # Section V.B as literally written); the depot model adds the missing absorption
    # compartment IM administration actually requires (Section V.A line 179).
    "KA_VALUES": (0.15, 0.20, 0.25),
    "BIOAVAILABILITY_OUT_OF_DEPOT": 1.0,   # all administered mass eventually reaches A1; no loss route in the depot itself
}


def odes_ivbolus(t, y, k12, k21, k_elim):
    a1, a2 = y
    return [-k12 * a1 - k_elim * a1 + k21 * a2, k12 * a1 - k21 * a2]


def simulate_continuous_ivbolus():
    """IV-BOLUS model -- the faithful reading of Section V.B as literally written
    (dose added directly to A1, no absorption compartment). Kept unchanged as the
    baseline/comparison run; the depot model below is the fix.

    One continuous integration across the whole study window, with a dose
    (an instantaneous +DOSE_UG jump to A1) injected as an event at each dose
    day. State carries over exactly (not reset) across dose boundaries --
    this is what makes it a single continuous run, not three independent
    ones. Returns both the full plotting curve (with duplicate timestamps at
    each jump, for an honest vertical-line plot) AND a list of clean
    per-dose segments (no duplicate boundary times) used for C_max/T_max/AUC,
    since slicing the plotting curve by time is ambiguous exactly at a jump."""
    k12, k21, k_elim = PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"]
    dose = PARAMS["DOSE_UG"]
    dose_hours = sorted(d * 24 for d in PARAMS["DOSE_DAYS"])
    end_hour = PARAMS["STUDY_END_DAY"] * 24
    boundaries = dose_hours[1:] + [end_hour]   # each dose's segment runs to the NEXT dose (or study end)

    all_t, all_a1, all_a2 = [], [], []
    segments = []   # one entry per dose: clean (t, a1) from just after that dose to just before the next
    y0 = [0.0, 0.0]
    t_prev = 0.0
    for i, dose_hour in enumerate(dose_hours):
        # advance from t_prev to dose_hour with whatever state currently holds (0 for dose 1,
        # nonzero residual for later doses) -- this is the tail of the PREVIOUS dose's segment
        if dose_hour > t_prev:
            n_pts = max(2, int((dose_hour - t_prev) * PARAMS["N_EVAL_PER_DAY"] / 24))
            t_eval = np.linspace(t_prev, dose_hour, n_pts, endpoint=False)
            sol = solve_ivp(odes_ivbolus, (t_prev, dose_hour), y0, args=(k12, k21, k_elim),
                             t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
            if not sol.success:
                raise RuntimeError(sol.message)
            all_t.extend(sol.t.tolist())
            all_a1.extend(sol.y[0].tolist())
            all_a2.extend(sol.y[1].tolist())
            y0 = [sol.y[0][-1] if len(sol.y[0]) else y0[0], sol.y[1][-1] if len(sol.y[1]) else y0[1]]

        # DOSE EVENT: instantaneous jump, A2 untouched, A1 += dose.
        all_t.append(dose_hour)
        all_a1.append(y0[0])          # pre-dose value (trough, if not the first dose) -- plotting only
        all_a2.append(y0[1])
        y0 = [y0[0] + dose, y0[1]]
        all_t.append(dose_hour)
        all_a1.append(y0[0])          # post-dose value (new peak-forming state) -- plotting only
        all_a2.append(y0[1])

        # This dose's OWN clean segment: solve from just after this dose to just before the next
        # (or study end). Using this dose's post-jump y0 as the sole initial condition for the
        # segment removes all boundary-time ambiguity from the metrics calculation.
        seg_end = boundaries[i]
        n_pts = max(2, int((seg_end - dose_hour) * PARAMS["N_EVAL_PER_DAY"] / 24))
        t_eval = np.linspace(dose_hour, seg_end, n_pts)
        seg_sol = solve_ivp(odes_ivbolus, (dose_hour, seg_end), y0, args=(k12, k21, k_elim),
                             t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
        if not seg_sol.success:
            raise RuntimeError(seg_sol.message)
        segments.append({"dose_number": i + 1, "dose_hour": dose_hour,
                          "t": seg_sol.t, "a1": seg_sol.y[0], "a2": seg_sol.y[1]})

        t_prev = dose_hour

    # BUG FIX: for every dose except the last, this loop's "advance to the next dose" step
    # (top of the next iteration) appends that same time range to the plotting curve again --
    # redundant but harmless. The LAST dose's own segment (computed just above, correct and
    # already used for its C_max/T_max/AUC/C_min) was never appended to the plotting arrays,
    # because there is no next iteration to do it. Without this, the full-curve CSV silently
    # ended at the second-to-last dose day instead of the true end of study -- caught only when
    # Step 5D tried to read Ag(t) past that point. Append it explicitly here.
    last_seg = segments[-1]
    all_t.extend(last_seg["t"].tolist())
    all_a1.extend(last_seg["a1"].tolist())
    all_a2.extend(last_seg["a2"].tolist())

    return np.array(all_t), np.array(all_a1), np.array(all_a2), segments


def per_dose_metrics_ivbolus(segments, vd):
    """C_max, T_max, AUC (to the next dose or end of study), C_min (trough
    just before the NEXT dose) per dose, computed from each dose's own clean
    segment (see simulate_continuous) -- no boundary-time ambiguity."""
    rows = []
    cumulative_auc = 0.0
    for seg in segments:
        t_seg, a1_seg = seg["t"], seg["a1"]
        conc = a1_seg / vd
        c_max = float(np.max(conc))
        t_max = float(t_seg[int(np.argmax(conc))] - seg["dose_hour"])
        auc = float(np.trapz(conc, t_seg))
        c_min = float(conc[-1])   # trough right before the next dose (or end of study for the last)
        cumulative_auc += auc
        rows.append({
            "Dose_number": seg["dose_number"],
            "Dose_day": seg["dose_hour"] / 24,
            "C_max_ug_per_L": round(c_max, 5),
            "T_max_h_post_dose": round(t_max, 3),
            "AUC_this_interval_ug_h_per_L": round(auc, 5),
            "C_min_trough_ug_per_L": round(c_min, 5),
            "Cumulative_AUC_ug_h_per_L": round(cumulative_auc, 5),
        })
    return rows


# =============================================================================
# DEPOT (IM ABSORPTION) MODEL -- added per STEP_5C_DEPOT_FIX_PROMPT.md.
#
#   dD/dt  = -ka*D
#   dA1/dt =  ka*D - k12*A1 - k_elim*A1 + k21*A2
#   dA2/dt =  k12*A1 - k21*A2
#   dE/dt  =  k_elim*A1      (cumulative eliminated mass -- a 4th state purely
#                              to make the mass-balance gate an independent
#                              numerical check, not an assumed algebraic identity)
#
# ka = k_IM (Section III.B line 44, 0.15-0.25 /h). D(0)=Dose, A1(0)=A2(0)=E(0)=0.
# k12/k21/k_elim unchanged from Step 5B/the IV-bolus model above.
# =============================================================================
def odes_depot(t, y, k12, k21, k_elim, ka):
    d, a1, a2, e = y
    dd = -ka * d
    da1 = ka * d - k12 * a1 - k_elim * a1 + k21 * a2
    da2 = k12 * a1 - k21 * a2
    de = k_elim * a1
    return [dd, da1, da2, de]


def verify_depot_single_dose():
    """Single 10ug dose, each ka, checked against the prompt's expected table
    before the full multi-dose schedule is trusted. Returns the per-ka
    (T_max, C_max/dose) so run() can assert against the expected values."""
    k12, k21, k_elim = PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"]
    dose = PARAMS["DOSE_UG"]
    results = {}
    for ka in PARAMS["KA_VALUES"]:
        t_eval = np.linspace(0, 200, 20001)
        sol = solve_ivp(odes_depot, (0, 200), [dose, 0.0, 0.0, 0.0], args=(k12, k21, k_elim, ka),
                         t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
        if not sol.success:
            raise RuntimeError(sol.message)
        a1 = sol.y[1]
        idx = int(np.argmax(a1))
        results[ka] = {"t_max": float(sol.t[idx]), "c_max_over_dose": float(a1[idx] / dose)}
    return results


def simulate_continuous_depot(ka):
    """Same continuous-with-dose-events structure as simulate_continuous_ivbolus,
    extended to the 4-state depot system. Returns the full plotting curve and
    clean per-dose segments (D, A1, A2, E each)."""
    k12, k21, k_elim = PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"]
    dose = PARAMS["DOSE_UG"]
    dose_hours = sorted(d * 24 for d in PARAMS["DOSE_DAYS"])
    end_hour = PARAMS["STUDY_END_DAY"] * 24
    boundaries = dose_hours[1:] + [end_hour]

    all_t, all_d, all_a1, all_a2, all_e = [], [], [], [], []
    segments = []
    y0 = [0.0, 0.0, 0.0, 0.0]
    t_prev = 0.0
    for i, dose_hour in enumerate(dose_hours):
        if dose_hour > t_prev:
            n_pts = max(2, int((dose_hour - t_prev) * PARAMS["N_EVAL_PER_DAY"] / 24))
            t_eval = np.linspace(t_prev, dose_hour, n_pts, endpoint=False)
            sol = solve_ivp(odes_depot, (t_prev, dose_hour), y0, args=(k12, k21, k_elim, ka),
                             t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
            if not sol.success:
                raise RuntimeError(sol.message)
            all_t.extend(sol.t.tolist())
            all_d.extend(sol.y[0].tolist())
            all_a1.extend(sol.y[1].tolist())
            all_a2.extend(sol.y[2].tolist())
            all_e.extend(sol.y[3].tolist())
            y0 = [sol.y[0][-1], sol.y[1][-1], sol.y[2][-1], sol.y[3][-1]] if len(sol.y[0]) else y0

        # DOSE EVENT: instantaneous jump to D only (drug re-enters the depot, not A1 directly --
        # this is the whole point of the fix: an IM injection re-loads the injection site, it
        # does not appear in the bloodstream instantaneously).
        all_t.append(dose_hour); all_d.append(y0[0]); all_a1.append(y0[1]); all_a2.append(y0[2]); all_e.append(y0[3])
        y0 = [y0[0] + dose, y0[1], y0[2], y0[3]]
        all_t.append(dose_hour); all_d.append(y0[0]); all_a1.append(y0[1]); all_a2.append(y0[2]); all_e.append(y0[3])

        seg_end = boundaries[i]
        n_pts = max(2, int((seg_end - dose_hour) * PARAMS["N_EVAL_PER_DAY"] / 24))
        t_eval = np.linspace(dose_hour, seg_end, n_pts)
        seg_sol = solve_ivp(odes_depot, (dose_hour, seg_end), y0, args=(k12, k21, k_elim, ka),
                             t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
        if not seg_sol.success:
            raise RuntimeError(seg_sol.message)
        segments.append({"dose_number": i + 1, "dose_hour": dose_hour, "t": seg_sol.t,
                          "d": seg_sol.y[0], "a1": seg_sol.y[1], "a2": seg_sol.y[2], "e": seg_sol.y[3]})
        t_prev = dose_hour

    # Same bug fix as simulate_continuous_ivbolus: the last dose's own segment was computed
    # (and used correctly for metrics) but never appended to the plotting-curve arrays, since
    # there is no next loop iteration to do it. Append it explicitly.
    last_seg = segments[-1]
    all_t.extend(last_seg["t"].tolist())
    all_d.extend(last_seg["d"].tolist())
    all_a1.extend(last_seg["a1"].tolist())
    all_a2.extend(last_seg["a2"].tolist())
    all_e.extend(last_seg["e"].tolist())

    return (np.array(all_t), np.array(all_d), np.array(all_a1), np.array(all_a2), np.array(all_e), segments)


def _refine_peak_depot(seg, k12, k21, k_elim, ka, vd, window_h=48.0):
    """The per-dose segment's own output grid is hourly (fine enough for AUC/
    trough, and for the slow IV-bolus curve where T_max=0 exactly), but too
    coarse to pinpoint an absorption peak that occurs within ~4-5h -- the
    single-dose validation used a 0.01h grid to hit the expected table, and
    the same resolution is needed here for T_max to be trustworthy. Re-solves
    a short dense window from the segment's own starting state (not stored to
    the main output CSV -- this is metric-refinement only)."""
    y0 = [seg["d"][0], seg["a1"][0], seg["a2"][0], seg["e"][0]]
    t_end = min(window_h, seg["t"][-1] - seg["dose_hour"])
    t_eval = np.linspace(0.0, t_end, 20001)
    sol = solve_ivp(odes_depot, (0.0, t_end), y0, args=(k12, k21, k_elim, ka),
                     t_eval=t_eval, method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"])
    if not sol.success:
        raise RuntimeError(sol.message)
    conc = sol.y[1] / vd
    idx = int(np.argmax(conc))
    return float(conc[idx]), float(sol.t[idx])


def per_dose_metrics_depot(segments, vd, k12, k21, k_elim, ka):
    rows = []
    cumulative_auc = 0.0
    for seg in segments:
        t_seg, a1_seg = seg["t"], seg["a1"]
        conc = a1_seg / vd
        c_max_refined, t_max_refined = _refine_peak_depot(seg, k12, k21, k_elim, ka, vd)
        auc = float(np.trapz(conc, t_seg))
        c_min = float(conc[-1])
        cumulative_auc += auc
        rows.append({
            "Dose_number": seg["dose_number"],
            "Dose_day": seg["dose_hour"] / 24,
            "C_max_ug_per_L": round(c_max_refined, 5),
            "T_max_h_post_dose": round(t_max_refined, 3),
            "AUC_this_interval_ug_h_per_L": round(auc, 5),
            "C_min_trough_ug_per_L": round(c_min, 5),
            "Cumulative_AUC_ug_h_per_L": round(cumulative_auc, 5),
            "D_end_of_interval_ug": round(float(seg["d"][-1]), 8),
        })
    return rows


def run():
    common.print_banner("PHASE 5C -- MULTI-DOSE SCHEDULE (Section V.C, spec lines 207-220)")

    ok_pre = common.verify_rule1("Step 5C pre-check")
    if not ok_pre:
        print("[ERROR] Rule 1 dirty before Step 5C even ran -- refusing to proceed.")
        sys.exit(2)

    ledger = common.Ledger()
    vd = 7.0   # midpoint of Step 5B's 6-8 L range, used consistently for concentration conversion

    print("\n" + "-" * 90)
    print(f"Simulating ONE continuous course: doses at days {PARAMS['DOSE_DAYS']}, "
          f"{PARAMS['DOSE_UG']} ug each, study 0-{PARAMS['STUDY_END_DAY']} days, V_d={vd} L (midpoint)")
    t, a1, a2, segments = simulate_continuous_ivbolus()

    neg_a1 = int(np.sum(a1 < -1e-9))
    neg_a2 = int(np.sum(a2 < -1e-9))
    print(f"  Negative A1 samples: {neg_a1}, negative A2 samples: {neg_a2} (must both be 0)")
    assert neg_a1 == 0 and neg_a2 == 0, "Negative compartment mass"

    rows = per_dose_metrics_ivbolus(segments, vd)

    print("\n" + "-" * 90)
    print("PER-DOSE METRICS (this model, corrected for 2-compartment washout per Step 5B)")
    for r in rows:
        print(f"  Dose {r['Dose_number']} (day {r['Dose_day']:.0f}): C_max={r['C_max_ug_per_L']:.4f} ug/L, "
              f"T_max={r['T_max_h_post_dose']:.2f}h post-dose, AUC(interval)={r['AUC_this_interval_ug_h_per_L']:.3f}, "
              f"C_min(trough)={r['C_min_trough_ug_per_L']:.5f} ug/L, cum.AUC={r['Cumulative_AUC_ug_h_per_L']:.3f}")

    print("\n" + "-" * 90)
    print("GATE -- every dose's C_max must be >= the first (baseline, zero-residual) dose's C_max")
    c_maxes = [r["C_max_ug_per_L"] for r in rows]
    # NOTE (corrected from an earlier, over-strong version of this gate): pre-dose residual is
    # ALWAYS >= 0 for this system (nonneg rate constants, nonneg compartments), so C_max_i =
    # (Dose + residual_i)/Vd >= Dose/Vd = C_max_1, for every i. That is the true, unconditional
    # invariant. It does NOT require C_max to be pairwise non-decreasing between consecutive
    # doses -- that stronger claim only holds if inter-dose intervals are non-increasing. Here
    # the gap grows (21 days, then 35 days), so MORE washout occurs before dose 3 than before
    # dose 2, and residual_3 < residual_2 is legitimate, not a bug (confirmed analytically below).
    ge_baseline = all(c >= c_maxes[0] - 1e-9 for c in c_maxes)
    print(f"  C_max sequence: {[f'{c:.5f}' for c in c_maxes]} (dose intervals: "
          f"{[PARAMS['DOSE_DAYS'][i+1]-PARAMS['DOSE_DAYS'][i] for i in range(len(PARAMS['DOSE_DAYS'])-1)]} days)")
    print(f"  All >= dose-1 baseline ({c_maxes[0]:.5f}): {ge_baseline}")
    assert ge_baseline, "A later dose's C_max fell below the first dose's baseline -- impossible for this LTI system, model or code error"
    pairwise_nondecreasing = all(c_maxes[i + 1] >= c_maxes[i] - 1e-9 for i in range(len(c_maxes) - 1))
    print(f"  Pairwise non-decreasing (consecutive doses): {pairwise_nondecreasing} -- NOT required when "
          "inter-dose intervals grow; reported for information only, not asserted as a gate.")
    print(f"  [SUCCESS] Every dose's C_max is >= the baseline ({c_maxes[0]:.5f} ug/L) -- the correct, "
          "unconditional invariant for this model.")
    for i, r in enumerate(rows):
        ledger.N(f"5C_ivbolus_dose{i+1}_Cmax_ug_L", r["C_max_ug_per_L"], "Phase5C_MultiDose_IVbolus")
        ledger.N(f"5C_ivbolus_dose{i+1}_Tmax_h", r["T_max_h_post_dose"], "Phase5C_MultiDose_IVbolus")
        ledger.N(f"5C_ivbolus_dose{i+1}_AUC_interval", r["AUC_this_interval_ug_h_per_L"], "Phase5C_MultiDose_IVbolus")
        ledger.N(f"5C_ivbolus_dose{i+1}_Cmin_trough_ug_L", r["C_min_trough_ug_per_L"], "Phase5C_MultiDose_IVbolus")
        ledger.N(f"5C_ivbolus_dose{i+1}_cumulative_AUC", r["Cumulative_AUC_ug_h_per_L"], "Phase5C_MultiDose_IVbolus")

    print("\n" + "-" * 90)
    print("COMPARISON WITH THE SPEC'S TABLE (lines 211-215) -- expected disagreement, logged not reproduced")
    spec_table = [
        {"dose": 1, "day": 0, "Cmax": 4.2, "Tmax": 6, "AUC": 19.2, "Cmin": 0.0},
        {"dose": 2, "day": 21, "Cmax": 3.8, "Tmax": 5, "AUC": 17.6, "Cmin": 0.15},
        {"dose": 3, "day": 56, "Cmax": 3.5, "Tmax": 4, "AUC": 16.2, "Cmin": 0.08},
    ]
    for spec_row, our_row in zip(spec_table, rows):
        print(f"  Dose {spec_row['dose']}: spec C_max={spec_row['Cmax']}, ours={our_row['C_max_ug_per_L']:.4f}; "
              f"spec C_min={spec_row['Cmin']}, ours={our_row['C_min_trough_ug_per_L']:.5f}")
    print("  Precise statement of the contradiction (corrected from an earlier, over-strong version): "
          "pre-dose residual is ALWAYS >= 0 in this model, so EVERY dose's C_max must be >= the first "
          "dose's C_max -- that is the unconditional invariant, true for any dosing interval. The "
          "spec's own dose-3 value (3.5) is LOWER than its own dose-1 value (4.2); that specific "
          "comparison is impossible. A same-magnitude or even decreasing step from dose 2 to dose 3 is "
          "NOT automatically impossible on its own -- our own model shows a small dose2->dose3 decrease "
          "too, because the spec's own schedule has a longer gap between doses 2 and 3 (35 days) than "
          "between doses 1 and 2 (21 days), giving more time to wash out. The decisive, unconditional "
          "problem is dose 3's C_max (3.5) falling BELOW dose 1's C_max (4.2) -- that specific "
          "relationship cannot happen regardless of interval spacing, since pre-dose-3 residual cannot "
          "be negative. We do NOT reverse-engineer parameters to reproduce the spec's table; we report "
          "our own values and log the contradiction.")
    print(f"  Magnitude, per the 5B follow-up: our own residuals are non-zero but SUB-1% of a single dose "
          f"(day-21 trough {rows[0]['C_min_trough_ug_per_L']:.5f} ug/L vs dose-1 C_max "
          f"{rows[0]['C_max_ug_per_L']:.5f} ug/L = {100*rows[0]['C_min_trough_ug_per_L']/rows[0]['C_max_ug_per_L']:.3f}% "
          "-- real but small, not a substantial accumulation. Our dose2->dose3 C_max actually decreases "
          f"slightly too ({c_maxes[1]:.5f} -> {c_maxes[2]:.5f}), which is legitimate given the longer gap, "
          f"and both remain >= the dose-1 baseline ({c_maxes[0]:.5f}).")

    cmax_growth_pct = 100 * (c_maxes[-1] - c_maxes[0]) / c_maxes[0]
    ledger.N("5C_ivbolus_Cmax_growth_dose1_to_dose3_pct", round(cmax_growth_pct, 4), "Phase5C_MultiDose_IVbolus")
    print(f"  Total C_max growth, dose 1 -> dose 3: {cmax_growth_pct:.4f}% (small, consistent with the "
          "sub-1% residual finding from Step 5B's corrected estimate).")

    out_dir = common.step_output_dir("C")
    ts = common.timestamp()

    curve_path = os.path.join(out_dir, f"Phase5C_FullCurve_IVbolus_{ts}.csv")
    with open(curve_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Time_h", "Time_day", "A1_ug", "A2_ug", "Concentration_ug_per_L"])
        for ti, a1i, a2i in zip(t, a1, a2):
            w.writerow([round(float(ti), 4), round(float(ti) / 24, 5), round(float(a1i), 8),
                        round(float(a2i), 8), round(float(a1i) / vd, 6)])
    print(f"\n  Full curve written: {curve_path} ({len(t)} points, includes dose-event jump rows)")

    metrics_path = os.path.join(out_dir, f"Phase5C_DoseMetrics_IVbolus_{ts}.csv")
    with open(metrics_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Dose metrics written: {metrics_path}")

    contradiction_path = os.path.join(out_dir, f"Phase5C_SpecContradiction_IVbolus_{ts}.csv")
    with open(contradiction_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Dose", "Day", "Spec_Cmax", "Our_Cmax", "Spec_Cmin", "Our_Cmin", "Spec_AUC", "Our_AUC"])
        for spec_row, our_row in zip(spec_table, rows):
            w.writerow([spec_row["dose"], spec_row["day"], spec_row["Cmax"], our_row["C_max_ug_per_L"],
                        spec_row["Cmin"], our_row["C_min_trough_ug_per_L"], spec_row["AUC"],
                        our_row["AUC_this_interval_ug_h_per_L"]])
    print(f"  Spec contradiction table written: {contradiction_path}")

    ledger_path = os.path.join(out_dir, f"Phase5C_Ledger_IVbolus_{ts}.csv")
    ledger.write(ledger_path)
    print(f"  Ledger written: {ledger_path} ({len(ledger.rows)} numbers)")

    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write("# Phase 5C -- Methodology Note\n\n")
        f.write("**Spec:** Section V.C, `REVISED_METHODOLOGY_MpoxHIV.md` lines 207-220.\n\n")
        f.write("**Model:** same 2-compartment ODE as Step 5B, integrated as ONE continuous run "
                "0-84 days with instantaneous +10ug dose events at days 0, 21, 56 -- not three "
                "independent single-dose runs. State (A1, A2) carries over exactly across each dose "
                f"boundary. Solver: `scipy.integrate.solve_ivp`, method=`{PARAMS['METHOD']}`, "
                f"rtol={PARAMS['RTOL']}, atol={PARAMS['ATOL']}, hourly output resolution "
                f"({PARAMS['N_EVAL_PER_DAY']}/day). Concentration = A1/V_d with V_d=7L (midpoint of "
                "Step 5B's 6-8L range).\n\n")
        f.write("**Consistent with the Step 5B correction, with the invariant stated precisely:** "
                "pre-dose residual is always >=0 in this model, so EVERY dose's C_max is >= dose 1's "
                f"C_max ({c_maxes[0]:.5f} ug/L) -- an unconditional invariant, verified as a gate here. "
                "This does NOT require C_max to rise pairwise from one dose to the next: this schedule's "
                "gaps grow (21 days, then 35 days), so more washout occurs before dose 3 than before "
                f"dose 2, and our own model shows dose 3's C_max ({c_maxes[2]:.5f}) dip slightly below "
                f"dose 2's ({c_maxes[1]:.5f}) while staying above dose 1's -- legitimate model behaviour, "
                f"not a bug. The residual is SUB-1% of a single dose throughout "
                f"({100*rows[0]['C_min_trough_ug_per_L']/rows[0]['C_max_ug_per_L']:.3f}% at the day-21 "
                f"trough) -- real but small, not a substantial accumulation.\n\n")
        f.write("**EXPECTED DISAGREEMENT with the spec's table (lines 211-215), logged not reproduced:** "
                "the spec's dose-3 C_max (3.5) is LOWER than its own dose-1 C_max (4.2) -- impossible "
                "under the unconditional invariant above, regardless of interval spacing or residual "
                "size, since pre-dose-3 residual cannot be negative. A same-magnitude or slightly lower "
                "dose-2-to-dose-3 step is not on its own impossible (we reproduce a small instance of it "
                "ourselves, for the reason above); the decisive problem is falling below the dose-1 "
                "baseline. We did not reverse-engineer parameters to reproduce the spec's table; our own "
                f"values are reported and the contradiction is logged in `Phase5C_SpecContradiction_IVbolus_{ts}.csv`.\n\n")
        f.write("**No escalation.** This is one of the brief's explicitly pre-identified non-triggers "
                "(\"The C_max-decreasing table (5C). Already identified as internally inconsistent.\").\n")
    print(f"  Methodology note written: {note_path}")

    # =========================================================================
    # DEPOT (IM ABSORPTION) REVISION -- STEP_5C_DEPOT_FIX_PROMPT.md
    # =========================================================================
    common.print_banner("PHASE 5C REVISION -- DEPOT (IM ABSORPTION) MODEL")

    print("\n" + "-" * 90)
    print("VALIDATION -- single 10ug dose, each ka, checked against the fix prompt's expected table "
          "BEFORE trusting the full multi-dose schedule")
    expected = {0.15: (4.57, 0.2264), 0.20: (3.92, 0.2685), 0.25: (3.48, 0.3040)}
    single_dose = verify_depot_single_dose()
    for ka in PARAMS["KA_VALUES"]:
        got_tmax, got_ratio = single_dose[ka]["t_max"], single_dose[ka]["c_max_over_dose"]
        exp_tmax, exp_ratio = expected[ka]
        print(f"  ka={ka}: T_max={got_tmax:.4f}h (expected {exp_tmax}h), "
              f"C_max/dose={got_ratio:.4f} (expected {exp_ratio})")
        assert abs(got_tmax - exp_tmax) < 0.1, f"ka={ka}: T_max does not match expected table to 0.1h"
        assert abs(got_ratio - exp_ratio) < 0.002, f"ka={ka}: C_max/dose does not match expected table"
        ledger.N(f"5C_depot_ka{ka}_singledose_Tmax_h", round(got_tmax, 4), "Phase5C (validation vs fix prompt)")
        ledger.N(f"5C_depot_ka{ka}_singledose_Cmax_over_dose", round(got_ratio, 4), "Phase5C (validation vs fix prompt)")
    print("  [SUCCESS] All three ka values reproduce the expected single-dose table to within 0.1h / 0.002.")

    depot_rows_all = []
    contradiction_depot_rows = []
    depot_full_curves = {}
    for ka in PARAMS["KA_VALUES"]:
        print("\n" + "-" * 90)
        print(f"FULL 5C SCHEDULE, DEPOT MODEL, ka={ka} /h")
        t_d, d_arr, a1_d, a2_d, e_arr, seg_d = simulate_continuous_depot(ka)
        depot_full_curves[ka] = (t_d, d_arr, a1_d, a2_d, e_arr)

        neg_d = int(np.sum(d_arr < -1e-9)); neg_a1d = int(np.sum(a1_d < -1e-9))
        neg_a2d = int(np.sum(a2_d < -1e-9)); neg_e = int(np.sum(e_arr < -1e-9))
        print(f"  Negative samples -- D:{neg_d} A1:{neg_a1d} A2:{neg_a2d} E:{neg_e} (all must be 0)")
        assert neg_d == 0 and neg_a1d == 0 and neg_a2d == 0 and neg_e == 0, f"ka={ka}: negative compartment mass"

        # Mass balance: checked on each dose's CLEAN segment (seg_d), not the plotting curve --
        # the plotting curve has duplicate timestamps at each dose event (pre-jump and post-jump
        # both recorded at the same t, for an honest vertical line), so a naive t<=dose_hour
        # comparison miscounts by exactly one dose's worth at those instants. Within segment i
        # (which starts right after dose i), administered-so-far is constant: i * DOSE_UG.
        mass_err = 0.0
        for seg in seg_d:
            expected_administered = seg["dose_number"] * PARAMS["DOSE_UG"]
            total_seg = seg["d"] + seg["a1"] + seg["a2"] + seg["e"]
            seg_err = float(np.max(np.abs(total_seg - expected_administered)))
            mass_err = max(mass_err, seg_err)
        print(f"  Mass balance: max|D+A1+A2+E - administered| over all segments = {mass_err:.2e} ug "
              f"(dose scale: {PARAMS['DOSE_UG']} ug)")
        assert mass_err < 1e-4, f"ka={ka}: mass balance violated by {mass_err:.2e} ug -- a leak exists"
        print("  [SUCCESS] Mass balance holds -- no leak.")

        d_end = float(seg_d[-1]["d"][-1])
        print(f"  Depot at end of study: D({PARAMS['STUDY_END_DAY']}d) = {d_end:.6e} ug (should be ~0)")
        assert d_end < PARAMS["DOSE_UG"] * 1e-4, f"ka={ka}: depot has not emptied by end of study"

        rows_d = per_dose_metrics_depot(seg_d, vd, PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"], ka)
        for r in rows_d:
            print(f"  Dose {r['Dose_number']} (day {r['Dose_day']:.0f}): C_max={r['C_max_ug_per_L']:.5f} ug/L, "
                  f"T_max={r['T_max_h_post_dose']:.2f}h post-dose, AUC={r['AUC_this_interval_ug_h_per_L']:.3f}, "
                  f"C_min={r['C_min_trough_ug_per_L']:.5f} ug/L, cum.AUC={r['Cumulative_AUC_ug_h_per_L']:.3f}, "
                  f"D_end={r['D_end_of_interval_ug']:.4f}")
            r_with_ka = dict(r); r_with_ka["Ka_per_h"] = ka
            depot_rows_all.append(r_with_ka)

        tmax_all_positive = all(r["T_max_h_post_dose"] > 0 for r in rows_d)
        print(f"  T_max > 0 for every dose: {tmax_all_positive} (this is the whole point of the fix)")
        assert tmax_all_positive, f"ka={ka}: a dose had T_max=0 -- the depot fix did not take effect"

        cmax_d = [r["C_max_ug_per_L"] for r in rows_d]
        ge_baseline_d = all(c >= cmax_d[0] - 1e-9 for c in cmax_d)
        print(f"  Every dose's C_max >= dose-1's C_max ({cmax_d[0]:.5f}): {ge_baseline_d} "
              "(our corrected invariant, not pairwise monotonicity -- Opus's original guidance on "
              "pairwise monotonicity does not apply to unequal dose spacing)")
        assert ge_baseline_d, f"ka={ka}: a later dose's C_max fell below dose 1's"

        cum_auc_depot = rows_d[-1]["Cumulative_AUC_ug_h_per_L"]
        cum_auc_ivbolus = rows[-1]["Cumulative_AUC_ug_h_per_L"]
        auc_pct_diff = 100 * abs(cum_auc_depot - cum_auc_ivbolus) / cum_auc_ivbolus
        print(f"  Cumulative AUC: depot={cum_auc_depot:.3f}, IV-bolus={cum_auc_ivbolus:.3f}, "
              f"difference={auc_pct_diff:.4f}% (gate: <1%, since bioavailability out of the depot is 1.0 "
              "-- absorption changes the curve's SHAPE, not its AREA)")
        assert auc_pct_diff < 1.0, f"ka={ka}: cumulative AUC differs from IV-bolus by >1% -- mass leak or bioavailability != 1"
        print(f"  [SUCCESS] ka={ka}: all gates passed.")

        for i, r in enumerate(rows_d):
            ledger.N(f"5C_depot_ka{ka}_dose{i+1}_Cmax_ug_L", r["C_max_ug_per_L"], "Phase5C_MultiDose_IMdepot")
            ledger.N(f"5C_depot_ka{ka}_dose{i+1}_Tmax_h", r["T_max_h_post_dose"], "Phase5C_MultiDose_IMdepot")
            ledger.N(f"5C_depot_ka{ka}_dose{i+1}_AUC_interval", r["AUC_this_interval_ug_h_per_L"], "Phase5C_MultiDose_IMdepot")
            ledger.N(f"5C_depot_ka{ka}_dose{i+1}_Cmin_trough_ug_L", r["C_min_trough_ug_per_L"], "Phase5C_MultiDose_IMdepot")
        ledger.N(f"5C_depot_ka{ka}_cumulative_AUC", round(cum_auc_depot, 5), "Phase5C_MultiDose_IMdepot")
        ledger.N(f"5C_depot_ka{ka}_AUC_pct_diff_vs_ivbolus", round(auc_pct_diff, 5), "Phase5C_MultiDose_IMdepot")

        for spec_row, our_row in zip(spec_table, rows_d):
            contradiction_depot_rows.append({
                "Ka_per_h": ka, "Dose": spec_row["dose"], "Day": spec_row["day"],
                "Spec_Cmax": spec_row["Cmax"], "Our_Cmax": our_row["C_max_ug_per_L"],
                "Spec_Cmin": spec_row["Cmin"], "Our_Cmin": our_row["C_min_trough_ug_per_L"],
                "Spec_AUC": spec_row["AUC"], "Our_AUC": our_row["AUC_this_interval_ug_h_per_L"],
                "Spec_Tmax": spec_row["Tmax"], "Our_Tmax": our_row["T_max_h_post_dose"],
            })

    print("\n" + "-" * 90)
    print("Depot vs IV-bolus, all ka: C_max is now attained hours after dosing, not at t=0 -- "
          "T_max lands in the spec's own stated 4-8h range (line 69), unlike the IV-bolus run.")

    depot_metrics_path = os.path.join(out_dir, f"Phase5C_DoseMetrics_IMdepot_{ts}.csv")
    with open(depot_metrics_path, "w", newline="") as f:
        fieldnames = ["Ka_per_h"] + [k for k in depot_rows_all[0].keys() if k != "Ka_per_h"]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(depot_rows_all)
    print(f"\n  Depot dose metrics written: {depot_metrics_path} ({len(depot_rows_all)} rows, "
          f"{len(PARAMS['KA_VALUES'])} ka values x 3 doses)")

    depot_curve_path = os.path.join(out_dir, f"Phase5C_FullCurve_IMdepot_{ts}.csv")
    with open(depot_curve_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ka_per_h", "Time_h", "Time_day", "D_ug", "A1_ug", "A2_ug", "Eliminated_ug", "Concentration_ug_per_L"])
        for ka in PARAMS["KA_VALUES"]:
            t_d, d_arr, a1_d, a2_d, e_arr = depot_full_curves[ka]
            for ti, di, a1i, a2i, ei in zip(t_d, d_arr, a1_d, a2_d, e_arr):
                w.writerow([ka, round(float(ti), 4), round(float(ti) / 24, 5), round(float(di), 8),
                            round(float(a1i), 8), round(float(a2i), 8), round(float(ei), 8), round(float(a1i) / vd, 6)])
    print(f"  Depot full curve written: {depot_curve_path}")

    depot_contradiction_path = os.path.join(out_dir, f"Phase5C_SpecContradiction_IMdepot_{ts}.csv")
    with open(depot_contradiction_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(contradiction_depot_rows[0].keys()))
        w.writeheader()
        w.writerows(contradiction_depot_rows)
    print(f"  Depot spec-contradiction table written: {depot_contradiction_path}")

    depot_ledger_path = os.path.join(out_dir, f"Phase5C_Ledger_IMdepot_{ts}.csv")
    depot_ledger_only = common.Ledger()
    depot_ledger_only.rows = [r for r in ledger.rows if "depot" in r["Key"]]
    depot_ledger_only.write(depot_ledger_path)
    print(f"  Depot ledger written: {depot_ledger_path} ({len(depot_ledger_only.rows)} numbers)")

    with open(note_path, "a") as f:
        f.write("\n---\n\n# Phase 5C -- Depot (IM Absorption) Revision\n\n")
        f.write("**Source:** `STEP_5C_DEPOT_FIX_PROMPT.md`. The IV-bolus model above is the faithful "
                "reading of Section V.B as literally written; it gives T_max=0h for every dose, which "
                "is correct FOR THAT MODEL but contradicts Section III.C's stated T_max of 4-8h (line "
                "69) and Section V.C's own table (6/5/4h) -- because Section V.B never defines an "
                "absorption route for the IM administration Section V.A specifies (line 179).\n\n")
        f.write("**Fix:** added a depot compartment using `k_IM` from Section III.B line 44 "
                "(`CL_depot = k_IM * LNP_dose`, k_IM=0.15-0.25/h): `dD/dt=-ka*D`, "
                "`dA1/dt=ka*D-k12*A1-k_elim*A1+k21*A2`, `dA2/dt=k12*A1-k21*A2`, D(0)=Dose, "
                "A1(0)=A2(0)=0. A 4th state `E` (cumulative eliminated mass, `dE/dt=k_elim*A1`) was "
                "added purely so mass balance is an independent numerical check, not an assumed "
                "algebraic identity.\n\n")
        f.write("**Validated against the fix prompt's expected single-dose table before trusting the "
                "full schedule** (all three ka, T_max within 0.1h, C_max/dose within 0.002):\n\n")
        for ka in PARAMS["KA_VALUES"]:
            f.write(f"- ka={ka}/h: T_max={single_dose[ka]['t_max']:.4f}h "
                    f"(expected {expected[ka][0]}h), C_max/dose={single_dose[ka]['c_max_over_dose']:.4f} "
                    f"(expected {expected[ka][1]})\n")
        f.write("\n**All gates passed for all three ka values:** T_max>0 for every dose; depot empties "
                "by end of study; mass balance (D+A1+A2+E = administered dose) holds to <1e-4 ug; no "
                "negative compartments; every dose's C_max >= dose 1's C_max (the corrected invariant, "
                "not pairwise monotonicity, per the fix prompt's explicit correction to the original "
                "Opus guidance); cumulative AUC within 1% of the IV-bolus run (bioavailability out of "
                "the depot = 1.0, so absorption reshapes the curve without changing its area).\n\n")
        f.write("**Manuscript corrections, 4 items (1-2 new, 3-4 carried from the IV-bolus run):**\n\n")
        f.write("1. Section V.B models an IV bolus, not the IM route the study specifies (T_max=0h vs "
                "the spec's own stated 4-8h and table values); adding the depot with k_IM from Section "
                f"III.B line 44 yields T_max {min(expected[k][0] for k in expected):.2f}-"
                f"{max(expected[k][0] for k in expected):.2f}h, consistent with the stated range. "
                "Recommend the model be respecified with three compartments.\n")
        f.write("2. The spec's C_max (4.2 ug/L at 10ug dose) implies V_d=2.38L, inconsistent with "
                "Section V.B line 203's stated V_d=6-8L (which gives 1.25-1.67 ug/L, IV-bolus; our "
                f"depot value at V_d=7L spans "
                f"{min(r['C_max_ug_per_L'] for r in depot_rows_all if r['Dose_number']==1):.4f}-"
                f"{max(r['C_max_ug_per_L'] for r in depot_rows_all if r['Dose_number']==1):.4f} ug/L "
                "across ka). A ~3x inconsistency internal to the spec.\n")
        f.write("3. The spec's dose-3 C_max (3.5) is unconditionally impossible: it is below the "
                "spec's own dose-1 C_max (4.2), and pre-dose-3 residual cannot be negative in either "
                "model, regardless of dose spacing.\n")
        f.write("4. Accumulation is negligible in both models, not substantial -- do not describe this "
                "schedule as producing meaningful accumulation.\n\n")
        f.write("**Both model variants are kept and clearly labelled** (`*_IVbolus_*`, `*_IMdepot_*`) "
                "-- the manuscript needs both to show why the change was made.\n")
    print(f"  Depot section appended to: {note_path}")

    ok_post = common.verify_rule1("Step 5C post-check")
    if not ok_post:
        print("[ERROR] Rule 1 violated during Step 5C.")
        sys.exit(2)

    print("\n[SUCCESS] Step 5C complete.")
    return ledger.rows


if __name__ == "__main__":
    run()
