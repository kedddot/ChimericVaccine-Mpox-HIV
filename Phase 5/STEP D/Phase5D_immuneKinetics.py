import os
import sys
import csv
import math

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_COMMON = os.path.normpath(os.path.join(_THIS_DIR, "..", "_common"))
if _COMMON not in sys.path:
    sys.path.insert(0, _COMMON)

import phase5_common as common
import numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr

# =============================================================================
# STEP 5D -- immune response kinetics (Section V.D, spec lines 222-258)
#
# Ag(t) comes from Step 5C's IM-DEPOT curve (never the IV-bolus one -- that
# was the whole point of the depot fix: antigen rises over ~4h, not
# instantly). Primary ka=0.20/h; sensitivity ka=0.15 and 0.25 reported
# alongside.
#
# ARBITRARY UNITS THROUGHOUT. No term here is calibrated to an assay. Only
# fold-change, time-to-peak and rank order are reported as findings; the
# spec's absolute pg/mL and log2 figures (lines 236-258) cannot be derived
# from these equations and are logged as a manuscript correction, not
# reproduced.
# =============================================================================

PARAMS = {
    "AG_SOURCE": "Phase5C_FullCurve_IMdepot",   # never the IV-bolus curve
    "KA_PRIMARY": 0.20,
    "KA_SENSITIVITY": (0.15, 0.25),
    # 180 days, NOT 84 (5C's dosing-schedule window) -- see run(): Ab_titer had not peaked by
    # day 84 for ANY (ka, scenario, b/g) combination when first tried (still rising at the end
    # of the window for all of them, end/peak=1.0000 exactly). This is not a bug: B_mem's decay
    # (delta, ~90-day half-life) is deliberately slow relative to the 84-day dosing schedule, and
    # the spec's own narrative (line 258) puts a "Day 180" persistence checkpoint, implying titres
    # are expected to still be near-plateau through day ~84 and only decline by day 180. The
    # window was extended to match the spec's own horizon so the model can actually be observed
    # to peak and decline, rather than forcing an artificial early peak into the equations.
    "STUDY_END_DAY": 180,
    "N_EVAL_PER_DAY": 24,
    "RTOL": 1e-10,
    "ATOL": 1e-12,
    "METHOD": "LSODA",   # Radau struggled with the interpolated-forcing RHS's discontinuous
                          # derivative at dose-event kinks; LSODA (variable-order, stiff-aware)
                          # handles a piecewise-linear-forced linear system robustly
    # T-cell equation (spec line 227): dT_eff/dt = a*Ag(t) - b*T_eff - g*T_eff*T_reg
    "B_RANGE_PER_DAY": (0.01, 0.02),     # T-cell death/exhaustion rate (spec)
    "G_RANGE_PER_DAY": (0.001, 0.005),   # T-reg suppression rate (spec)
    # Spec ambiguity (documented in METHODOLOGY_NOTE.md): g's own stated units already carry a
    # "(T_reg/T_total)" factor, and no T_reg(t) trajectory is defined anywhere in the spec. We
    # fold g*T_reg into a single effective suppression rate -(b+g)*T_eff, i.e. treat b and g as
    # two additive per-day rate constants rather than inventing an unstated T_reg(t) state.
    "BG_GRID": ((0.01, 0.001), (0.015, 0.003), (0.02, 0.005)),   # (b, g) low/mid/high combos

    # GC_output coupling (spec: "peaks 10-14 days post-antigen", coupling undefined -- our
    # choice, documented and numerically verified against this project's own Ag(t) curve, not
    # asserted blind): a two-stage transit (Erlang-2) delay filter on Ag(t):
    #   dG1/dt = Ag(t) - G1/tau ; dGC_output/dt = G1/tau - GC_output/tau
    # tau=7 days was found by sweeping tau against the actual Ag(t) curve (ka=0.20) and checking
    # where the PRIME-induced GC_output peak falls before boost 1 (day 21): tau=6->day 10.2,
    # tau=7->day 11.4, tau=8->day 12.6 -- all three land in the spec's 10-14d window; 7 (the
    # window's approximate midpoint) is used as the primary value.
    "TAU_GC_DAYS": 7.0,

    # B-cell/antibody cascade (spec lines 242-251): dB_mem/dt = s*GC_output - d*B_mem;
    # dPlasma_cell/dt = f*B_mem - r*Plasma_cell; dAb_titer/dt = p*Plasma_cell - z*Ab_titer.
    # The spec supplies a NUMBER for none of s, d, f, r, p -- only z is derivable (from the
    # stated IgG/IgM half-lives). Two of the missing four (s, p) are pure linear GAIN constants:
    # for this linear cascade, rescaling s or p rescales B_mem/Ab_titer's AMPLITUDE only, not
    # their timing, peak day, or the FOLD-CHANGE ratios between doses (which are gain-invariant
    # for a linear system) -- exactly the quantities this step is allowed to report (section 4).
    # They are set to 1.0 (arbitrary units) without loss of generality for any reported metric.
    # f (phi) is also a pure gain on the B_mem->Plasma_cell transfer for the same reason and is
    # likewise set to 1.0. d (delta, B_mem decay) and r (rho, Plasma_cell decay) DO set timing/
    # shape and have no derivation path in the spec; they are assigned standard immunology
    # literature-order-of-magnitude defaults (memory B-cell persistence ~90 days -> d=ln2/90;
    # plasma cell lifespan ~10 days -> r=ln2/10), explicitly flagged as illustrative, not derived
    # from this project's own data or the spec, and are documented as a manuscript-correction
    # item, not silently presented as if given by the spec.
    "SIGMA": 1.0, "PHI": 1.0, "PSI": 1.0,
    "DELTA_PER_DAY": math.log(2) / 90.0,   # memory B-cell persistence ~90 days (illustrative)
    "RHO_PER_DAY": math.log(2) / 10.0,     # plasma cell lifespan ~10 days (illustrative)
    "IGG_HALFLIFE_DAYS": 21.0,
    "IGM_HALFLIFE_DAYS": 5.0,

    # Alpha scenarios (Section 3 of the fix prompt) -- Phase IV ledger values, n=10 MHC-I
    # construct epitopes. BigMHC-IM: probability, higher=more immunogenic. PRIME %Rank: a rank,
    # lower=better. OPPOSITE DIRECTIONS -- direction-checked against Phase4B EL rank below
    # before use (both already verified in Phase IV Step 3; re-cited here, not re-derived).
    "BIGMHC_MEAN": 0.3502, "BIGMHC_HIV": 0.382, "BIGMHC_MPOX": 0.318,
    "PRIME_MEAN": 0.642, "PRIME_HIV": 0.278, "PRIME_MPOX": 0.177,
    "ALPHA_UNIFORM": 1.0,   # arbitrary units baseline -- constant alpha, as the spec implies
}


# ------------------------------------------------------------------ Ag(t) loader
def load_ag_curve(ka):
    path_candidates = [p for p in os.listdir(os.path.join(common.PROJECT_ROOT, "Step_Outputs", "Phase5", "StepC"))
                        if p.startswith("Phase5C_FullCurve_IMdepot_") and p.endswith(".csv")]
    if not path_candidates:
        raise FileNotFoundError("No Phase5C_FullCurve_IMdepot_*.csv found -- run Step 5C's depot revision first.")
    path_candidates.sort()
    path = os.path.join(common.PROJECT_ROOT, "Step_Outputs", "Phase5", "StepC", path_candidates[-1])
    t_days, conc = [], []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if abs(float(r["Ka_per_h"]) - ka) < 1e-9:
                t_days.append(float(r["Time_day"]))
                conc.append(float(r["Concentration_ug_per_L"]))
    pairs = sorted(zip(t_days, conc), key=lambda x: x[0])
    return np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs]), os.path.basename(path)


def make_ag_interpolator(t_arr, c_arr):
    def ag_of_t(t):
        return np.interp(t, t_arr, c_arr)
    return ag_of_t


# ------------------------------------------------------------------ combined ODE system
# State vector y = [T_eff, G1, GC_output, B_mem, Plasma_cell, Ab_titer]
def odes_full(t, y, ag_of_t, alpha, b, g, tau, sigma, delta, phi, rho, psi, zeta):
    t_eff, g1, gc, b_mem, plasma, ab = y
    ag = ag_of_t(t)
    d_t_eff = alpha * ag - (b + g) * t_eff
    d_g1 = ag - g1 / tau
    d_gc = g1 / tau - gc / tau
    d_b_mem = sigma * gc - delta * b_mem
    d_plasma = phi * b_mem - rho * plasma
    d_ab = psi * plasma - zeta * ab
    return [d_t_eff, d_g1, d_gc, d_b_mem, d_plasma, d_ab]


def run_scenario(ag_of_t, alpha, b, g, zeta, end_day=None):
    p = PARAMS
    end_day = end_day or p["STUDY_END_DAY"]
    n_pts = int(end_day * p["N_EVAL_PER_DAY"])
    t_eval = np.linspace(0, end_day, n_pts + 1)
    y0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    sol = solve_ivp(odes_full, (0, end_day), y0,
                     args=(ag_of_t, alpha, b, g, p["TAU_GC_DAYS"], p["SIGMA"], p["DELTA_PER_DAY"],
                           p["PHI"], p["RHO_PER_DAY"], p["PSI"], zeta),
                     t_eval=t_eval, method=p["METHOD"], rtol=p["RTOL"], atol=p["ATOL"])
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.t, sol.y


def peak_time_and_value(t, y):
    idx = int(np.argmax(y))
    return float(t[idx]), float(y[idx])


def zeta_from_halflife(t_half_days):
    return math.log(2) / t_half_days


def run():
    common.print_banner("PHASE 5D -- IMMUNE RESPONSE KINETICS (Section V.D, spec lines 222-258)")

    ok_pre = common.verify_rule1("Step 5D pre-check")
    if not ok_pre:
        print("[ERROR] Rule 1 dirty before Step 5D even ran -- refusing to proceed.")
        sys.exit(2)

    ledger = common.Ledger()
    p = PARAMS
    out_dir = common.step_output_dir("D")
    ts = common.timestamp()

    print("\n" + "-" * 90)
    print("DIRECTION CHECK -- BigMHC-IM and PRIME %Rank, re-cited from Phase IV Step 3 (not re-derived)")
    print("  BigMHC-IM best score vs Phase4B EL rank: Pearson r=-0.2984, p=0.4024, n=10 -- "
          "expected negative (higher BigMHC-IM = more immunogenic = should track a LOWER, better, "
          "EL rank); right sign, not distinguishable from zero at n=10.")
    print("  PRIME %Rank vs Phase4B EL rank: Pearson r=+0.9159, p=0.0002, n=10 -- expected positive "
          "(both are rank-like, lower=better, so they should co-vary positively); meaningful.")
    print("  Both pass their direction check. BigMHC-IM (higher=better) and PRIME %Rank (lower=better) "
          "run in OPPOSITE directions -- confirmed before use below.")
    ledger.N("5D_direction_check_bigmhc_r", -0.2984, "Phase4C_DirectionChecks (re-cited)")
    ledger.N("5D_direction_check_prime_r", 0.9159, "Phase4C_DirectionChecks (re-cited)")

    print("\n" + "-" * 90)
    print(f"LOADING Ag(t) from the IM-DEPOT curve (never IV-bolus) -- primary ka={p['KA_PRIMARY']}, "
          f"sensitivity ka={p['KA_SENSITIVITY']}")
    ag_curves = {}
    for ka in (p["KA_PRIMARY"],) + p["KA_SENSITIVITY"]:
        t_arr, c_arr, src_file = load_ag_curve(ka)
        ag_curves[ka] = (make_ag_interpolator(t_arr, c_arr), src_file)
        print(f"  ka={ka}: {len(t_arr)} points loaded from {src_file}")

    zeta_igg = zeta_from_halflife(p["IGG_HALFLIFE_DAYS"])
    zeta_igm = zeta_from_halflife(p["IGM_HALFLIFE_DAYS"])
    print(f"\n  zeta (IgG, t_half={p['IGG_HALFLIFE_DAYS']}d) = {zeta_igg:.5f} /day")
    print(f"  zeta (IgM, t_half={p['IGM_HALFLIFE_DAYS']}d) = {zeta_igm:.5f} /day")
    ledger.N("5D_zeta_IgG_per_day", round(zeta_igg, 5), "Phase5D (derived from stated half-life)")
    ledger.N("5D_zeta_IgM_per_day", round(zeta_igm, 5), "Phase5D (derived from stated half-life)")

    print("\n" + "-" * 90)
    print("ALPHA SCENARIOS (three, all reported -- no winner picked)")
    alpha_scenarios = {
        "1_uniform": {"HIV": p["ALPHA_UNIFORM"], "Mpox": p["ALPHA_UNIFORM"]},
        "2_bigmhc_weighted": {"HIV": p["BIGMHC_MEAN"], "Mpox": p["BIGMHC_MEAN"]},
        "3_per_pathogen": {"HIV": p["BIGMHC_HIV"], "Mpox": p["BIGMHC_MPOX"]},
    }
    for name, vals in alpha_scenarios.items():
        print(f"  {name}: HIV alpha={vals['HIV']}, Mpox alpha={vals['Mpox']}")
        ledger.N(f"5D_alpha_{name}_HIV", vals["HIV"], "STEP_5D_PROMPT.md / Phase4G ledger")
        ledger.N(f"5D_alpha_{name}_Mpox", vals["Mpox"], "STEP_5D_PROMPT.md / Phase4G ledger")
    print("  Scenario 2 uses the BigMHC-IM CONSTRUCT mean (0.3502) applied uniformly (no HIV/Mpox "
          "split) -- it tests whether weighting alpha by overall immunogenicity strength (vs the "
          "arbitrary uniform 1.0) changes T-cell kinetics shape, not the pathogen split. Scenario 3 "
          "is the one that tests the pathogen split, using BigMHC-IM's own HIV/Mpox means (0.382/0.318) "
          "since BigMHC-IM is the primary rung of the Phase IV TCR ladder and has the more direct "
          "immunogenicity interpretation (a probability); PRIME's HIV/Mpox split (0.278/0.177, lower="
          "better) is reported as a sensitivity check on scenario 3, not a fourth scenario.")

    # SENSITIVITY CHECK, caught by the direction discipline the prompt explicitly warns about:
    # does using PRIME instead of BigMHC-IM for the per-pathogen ratio agree on which pathogen
    # scores higher? PRIME is lower=better, so it is inverted (1/PRIME_rank) to a higher=better
    # quantity, matching BigMHC-IM's convention, before comparing.
    bigmhc_ratio = p["BIGMHC_MPOX"] / p["BIGMHC_HIV"]
    prime_ratio_inverted = (1.0 / p["PRIME_MPOX"]) / (1.0 / p["PRIME_HIV"])   # = PRIME_HIV/PRIME_MPOX
    print(f"  Cross-check: Mpox/HIV ratio via BigMHC-IM = {bigmhc_ratio:.4f} (<1: Mpox trails HIV) "
          f"vs via inverted PRIME %Rank = {prime_ratio_inverted:.4f} (>1: Mpox EXCEEDS HIV).")
    print("  *** The two tools DISAGREE on direction. *** BigMHC-IM says Mpox is less immunogenic "
          "than HIV; PRIME's own per-pathogen means say the opposite (Mpox's raw %Rank, 0.177, is "
          "numerically LOWER/better than HIV's, 0.278). This is exactly the kind of reversal the "
          "brief warns a sign error could produce silently -- here it is a genuine tool disagreement, "
          "not a coding error (both ratios were independently hand-verified). Scenario 3 uses "
          "BigMHC-IM per the fix prompt's explicit instruction (the primary TCR-ladder rung, a direct "
          "immunogenicity probability); this disagreement is reported as an explicit caveat on "
          "scenario 3, not resolved by picking whichever tool agrees with the narrative.")
    ledger.N("5D_mpox_hiv_ratio_bigmhc", round(bigmhc_ratio, 4), "Phase4G ledger (bm_hiv, bm_mpx)")
    ledger.N("5D_mpox_hiv_ratio_prime_inverted", round(prime_ratio_inverted, 4), "Phase4G ledger (pr_hiv, pr_mpx)")

    print("\n" + "-" * 90)
    print(f"SIMULATING -- {len(alpha_scenarios)} alpha scenarios x 2 pathogens x {len(p['BG_GRID'])} "
          f"(b,g) combos x {1+len(p['KA_SENSITIVITY'])} ka values")
    timecourse_rows = []
    summary_rows = []
    igg_ab_curves = {}   # keyed (scenario, pathogen, ka, bg_label) -> (t, T_eff, GC, B_mem, Plasma, Ab)

    for ka in (p["KA_PRIMARY"],) + p["KA_SENSITIVITY"]:
        ag_of_t, src_file = ag_curves[ka]
        is_primary_ka = (ka == p["KA_PRIMARY"])
        for scen_name, vals in alpha_scenarios.items():
            for pathogen in ("HIV", "Mpox"):
                alpha = vals[pathogen]
                for bg_idx, (b, g) in enumerate(p["BG_GRID"]):
                    bg_label = ("low", "mid", "high")[bg_idx]
                    t, y = run_scenario(ag_of_t, alpha, b, g, zeta_igg)
                    t_eff, g1, gc, b_mem, plasma, ab = y

                    neg = int(np.sum(y < -1e-9))
                    if neg > 0:
                        raise AssertionError(f"Negative state value: ka={ka} scen={scen_name} "
                                              f"pathogen={pathogen} bg={bg_label} -- {neg} samples")

                    t_peak_ag, v_peak_ag = peak_time_and_value(t, ag_of_t(t))
                    t_peak_teff, v_peak_teff = peak_time_and_value(t, t_eff)
                    t_peak_gc, v_peak_gc = peak_time_and_value(t[t <= 21], gc[t <= 21])
                    t_peak_ab, v_peak_ab = peak_time_and_value(t, ab)

                    # end-of-study return-toward-baseline check (nothing diverges)
                    end_over_peak = {
                        "T_eff": t_eff[-1] / max(v_peak_teff, 1e-12),
                        "Ab_titer": ab[-1] / max(v_peak_ab, 1e-12),
                    }

                    summary_rows.append({
                        "Ka_per_h": ka, "Alpha_scenario": scen_name, "Pathogen": pathogen,
                        "b_per_day": b, "g_per_day": g, "bg_label": bg_label,
                        "Ag_peak_day": round(t_peak_ag, 3), "Teff_peak_day": round(t_peak_teff, 3),
                        "GC_peak_day_before_boost1": round(t_peak_gc, 3),
                        "Ab_peak_day": round(t_peak_ab, 3),
                        "Teff_peak_value": round(v_peak_teff, 6),
                        "GC_peak_value": round(v_peak_gc, 6),
                        "Ab_peak_value": round(v_peak_ab, 6),
                        "Teff_end_over_peak": round(end_over_peak["T_eff"], 5),
                        "Ab_end_over_peak": round(end_over_peak["Ab_titer"], 5),
                        "Ordering_OK": bool(t_peak_ag <= t_peak_teff <= t_peak_ab),
                    })

                    if is_primary_ka and bg_label == "mid":
                        igg_ab_curves[(scen_name, pathogen)] = (t, t_eff, gc, b_mem, plasma, ab)
                        # store a decimated time-course (daily) for the primary scenario set only --
                        # the full grid's summary (above) already captures every (ka, scenario, bg) combo
                        for i in range(0, len(t), p["N_EVAL_PER_DAY"]):
                            timecourse_rows.append({
                                "Ka_per_h": ka, "Alpha_scenario": scen_name, "Pathogen": pathogen,
                                "bg_label": bg_label, "Time_day": round(float(t[i]), 3),
                                "Ag_ug_per_L": round(float(ag_of_t(t[i])), 6),
                                "T_eff": round(float(t_eff[i]), 6), "GC_output": round(float(gc[i]), 6),
                                "B_mem": round(float(b_mem[i]), 6), "Plasma_cell": round(float(plasma[i]), 6),
                                "Ab_titer": round(float(ab[i]), 6),
                            })
    print(f"  {len(summary_rows)} (ka, scenario, pathogen, b/g) combinations simulated.")

    print("\n" + "-" * 90)
    print("GATES")
    n_ok_order = sum(1 for r in summary_rows if r["Ordering_OK"])
    print(f"  Ag peak <= T_eff peak <= Ab peak: {n_ok_order}/{len(summary_rows)} combinations")
    assert n_ok_order == len(summary_rows), "Ordering violated in at least one combination -- coupling is wrong"
    print("  [SUCCESS] Ordering holds everywhere: antigen peak, then T-cell peak, then antibody peak.")

    gc_days = [r["GC_peak_day_before_boost1"] for r in summary_rows]
    n_in_window = sum(1 for d in gc_days if 10.0 <= d <= 14.0)
    print(f"  GC_output (prime-induced, pre-boost1) peak in the spec's 10-14 day window: "
          f"{n_in_window}/{len(gc_days)} (range observed: {min(gc_days):.2f}-{max(gc_days):.2f} days)")
    assert n_in_window == len(gc_days), "GC_output peak fell outside the spec's 10-14 day window in at least one combination"
    print("  [SUCCESS] GC_output peak lands in the 10-14 day window for every combination (tau=7 days).")

    max_end_over_peak = max(max(r["Teff_end_over_peak"], r["Ab_end_over_peak"]) for r in summary_rows)
    print(f"  Return-toward-baseline check: max(end-of-study / peak) across T_eff and Ab_titer = "
          f"{max_end_over_peak:.4f} (should be well below 1.0 -- nothing diverges)")
    assert max_end_over_peak < 0.95, "A state variable did not decay from its peak by end of study -- possible divergence"
    print("  [SUCCESS] All states return toward baseline; nothing diverges.")

    ledger.N("5D_n_combinations_simulated", len(summary_rows), "Phase5D_Summary")
    ledger.N("5D_gc_peak_day_min", round(min(gc_days), 3), "Phase5D_Summary")
    ledger.N("5D_gc_peak_day_max", round(max(gc_days), 3), "Phase5D_Summary")
    ledger.N("5D_tau_gc_days", p["TAU_GC_DAYS"], "Phase5D (numerically verified against Ag(t))")

    print("\n" + "-" * 90)
    print("SCENARIO 3 (per-pathogen alpha, BigMHC-IM): the interesting comparison, primary ka/bg")
    for scen_name in alpha_scenarios:
        hiv_row = next(r for r in summary_rows if r["Ka_per_h"] == p["KA_PRIMARY"] and r["Alpha_scenario"] == scen_name
                        and r["Pathogen"] == "HIV" and r["bg_label"] == "mid")
        mpx_row = next(r for r in summary_rows if r["Ka_per_h"] == p["KA_PRIMARY"] and r["Alpha_scenario"] == scen_name
                        and r["Pathogen"] == "Mpox" and r["bg_label"] == "mid")
        ratio = mpx_row["Teff_peak_value"] / hiv_row["Teff_peak_value"] if hiv_row["Teff_peak_value"] else float("nan")
        print(f"  {scen_name}: T_eff peak HIV={hiv_row['Teff_peak_value']:.4f} at day {hiv_row['Teff_peak_day']:.2f}, "
              f"Mpox={mpx_row['Teff_peak_value']:.4f} at day {mpx_row['Teff_peak_day']:.2f} "
              f"(Mpox/HIV peak ratio = {ratio:.4f})")
        ledger.N(f"5D_{scen_name}_Teff_peak_HIV", hiv_row["Teff_peak_value"], "Phase5D_Summary")
        ledger.N(f"5D_{scen_name}_Teff_peak_Mpox", mpx_row["Teff_peak_value"], "Phase5D_Summary")
        ledger.N(f"5D_{scen_name}_Mpox_HIV_ratio", round(ratio, 4), "Phase5D_Summary")
    scen3_ratio = next(r for r in summary_rows if r["Ka_per_h"] == p["KA_PRIMARY"] and r["Alpha_scenario"] == "3_per_pathogen"
                        and r["Pathogen"] == "Mpox" and r["bg_label"] == "mid")["Teff_peak_value"] / \
                  next(r for r in summary_rows if r["Ka_per_h"] == p["KA_PRIMARY"] and r["Alpha_scenario"] == "3_per_pathogen"
                        and r["Pathogen"] == "HIV" and r["bg_label"] == "mid")["Teff_peak_value"]
    print(f"  Scenario 3 magnitude: Mpox T-cell peak is {scen3_ratio*100:.1f}% of HIV's (a "
          f"{(1-scen3_ratio)*100:.1f} percentage-point gap) -- this equals the input alpha ratio "
          f"({bigmhc_ratio:.4f}) exactly, because the T-cell equation is linear in alpha and Ag(t) "
          "is identical for both pathogens (same construct, same PK) -- expected, not a new result.")
    ledger.N("5D_scenario3_Mpox_pct_of_HIV_Teff_peak", round(scen3_ratio * 100, 2), "Phase5D_Summary")

    print("\n" + "-" * 90)
    print("SENSITIVITY -- ka=0.15 vs 0.20 vs 0.25 (scenario 1, uniform alpha, mid b/g)")
    for ka in (p["KA_PRIMARY"],) + p["KA_SENSITIVITY"]:
        row = next(r for r in summary_rows if r["Ka_per_h"] == ka and r["Alpha_scenario"] == "1_uniform"
                   and r["Pathogen"] == "HIV" and r["bg_label"] == "mid")
        print(f"  ka={ka}: T_eff peak day={row['Teff_peak_day']:.2f}, value={row['Teff_peak_value']:.4f}; "
              f"Ab peak day={row['Ab_peak_day']:.2f}")

    print("\n" + "-" * 90)
    print("(b, g) SWEEP -- effect on T_eff peak timing/magnitude (scenario 1, primary ka, HIV)")
    for bg_idx, (b, g) in enumerate(p["BG_GRID"]):
        bg_label = ("low", "mid", "high")[bg_idx]
        row = next(r for r in summary_rows if r["Ka_per_h"] == p["KA_PRIMARY"] and r["Alpha_scenario"] == "1_uniform"
                   and r["Pathogen"] == "HIV" and r["bg_label"] == bg_label)
        print(f"  b={b}, g={g} ({bg_label}): T_eff peak day={row['Teff_peak_day']:.2f}, "
              f"value={row['Teff_peak_value']:.4f}")

    print("\n" + "-" * 90)
    print("DO NOT OVER-CLAIM -- absolute units check")
    print("  The spec's lines 236-258 state absolute IFN-gamma (50-100, 300-600, 400-700 pg/mL) and "
          "antibody titres (4-6, 8-10, 9-11 Log2) at named study days. No term in this model maps "
          "T_eff or Ab_titer (arbitrary units) onto either assay scale -- neither is reproduced. Only "
          "fold-change, time-to-peak, and rank order are reported below and in the CSVs.")
    print("  Fold-change (T_eff peak, scenario 1 uniform, primary ka/bg): dose response is encoded "
          "directly in Ag(t), already reported via the 5C dose metrics; this step's own reportable "
          "fold-change is scenario-to-scenario and pathogen-to-pathogen (see scenario 3 above), not a "
          "dose fold-change (Ag(t) already IS the multi-dose curve, T_eff responds continuously to it).")

    # ---------------------------------------------------------------- outputs
    summary_path = os.path.join(out_dir, f"Phase5D_Summary_{ts}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\n  Summary written: {summary_path} ({len(summary_rows)} rows: every ka x scenario x "
          "pathogen x b/g combination)")

    timecourse_path = os.path.join(out_dir, f"Phase5D_TimeCourse_{ts}.csv")
    with open(timecourse_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(timecourse_rows[0].keys()))
        w.writeheader()
        w.writerows(timecourse_rows)
    print(f"  Time-course written: {timecourse_path} ({len(timecourse_rows)} rows: daily resolution, "
          "primary bg='mid', all ka x scenario x pathogen)")

    # -------------------------------------------------------------- figures
    figures_written = []
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(11, 8))
        for pathogen, color in (("HIV", "tab:blue"), ("Mpox", "tab:red")):
            t, t_eff, gc, b_mem, plasma, ab = igg_ab_curves[("3_per_pathogen", pathogen)]
            axes[0, 0].plot(t, ag_curves[p["KA_PRIMARY"]][0](t), label=f"Ag(t) [{pathogen} curve identical]" if pathogen == "HIV" else None, color="gray", alpha=0.6)
            axes[0, 0].set_title("Ag(t), primary ka=0.20")
            axes[0, 1].plot(t, t_eff, label=pathogen, color=color)
            axes[0, 1].set_title("T_eff -- scenario 3 (per-pathogen alpha)")
            axes[1, 0].plot(t, gc, label=pathogen, color=color)
            axes[1, 0].set_title("GC_output")
            axes[1, 1].plot(t, ab, label=pathogen, color=color)
            axes[1, 1].set_title("Ab_titer (arbitrary units)")
        for ax in axes.flat:
            ax.set_xlabel("day"); ax.legend(fontsize=8)
        fig.tight_layout()
        fig_path = os.path.join(out_dir, f"Phase5D_Scenario3_Curves_{ts}.png")
        fig.savefig(fig_path, dpi=110)
        plt.close(fig)
        figures_written.append(fig_path)
        print(f"  Figure written: {fig_path}")
    except Exception as e:
        print(f"  [WARNING] matplotlib figure generation failed ({e}) -- CSV-only output for this step.")

    ledger_path = os.path.join(out_dir, f"Phase5D_Ledger_{ts}.csv")
    ledger.write(ledger_path)
    print(f"  Ledger written: {ledger_path} ({len(ledger.rows)} numbers)")

    corrections_path = os.path.join(out_dir, "MANUSCRIPT_CORRECTIONS_5D.md")
    with open(corrections_path, "w") as f:
        f.write("# Manuscript corrections -- Section V.D (Step 5D)\n\n")
        f.write("**1. Absolute IFN-gamma (pg/mL) and antibody titre (Log2) figures at spec lines "
                "236-258 cannot be derived from the spec's own equations.** The T-cell and B-cell/"
                "antibody ODEs (lines 227, 242-251) are unitless -- no term maps model output onto "
                "an assay scale. This step's T_eff and Ab_titer are dimensionless model states; only "
                "fold-change, time-to-peak, and rank order are reported as findings.\n\n")
        f.write("**2. Spec parameters with no supplied numeric value: sigma, delta, phi, rho "
                "(only zeta is derivable, from the stated IgG/IgM half-lives).** Sigma and psi (p) "
                "are pure linear gain constants for this cascade -- they rescale amplitude only, not "
                "timing, peak day, or fold-change ratios (the reportable quantities) -- and were set "
                "to 1.0 without loss of generality. Phi is a gain on the B_mem->Plasma_cell transfer "
                "for the same reason and was also set to 1.0. Delta and rho DO set timing/shape and "
                f"have no derivation path in the spec; they were assigned illustrative literature-order "
                f"defaults (delta=ln2/90d, memory B-cell persistence; rho=ln2/10d, plasma-cell "
                "lifespan), explicitly NOT derived from this project's own data, and are flagged here "
                "rather than presented as spec-supplied.\n\n")
        f.write("**3. The GC_output coupling to Ag(t) is undefined in the spec ('peaks 10-14 days "
                "post-antigen', no equation given).** A two-stage transit (Erlang-2) delay filter "
                "was used, `dG1/dt=Ag(t)-G1/tau`, `dGC/dt=G1/tau-GC/tau`. tau=7 days was chosen by "
                "sweeping tau against this project's own Ag(t) curve (ka=0.20) and verifying the "
                "resulting prime-induced peak lands in the spec's stated window (tau=6->day 10.2, "
                "tau=7->day 11.4, tau=8->day 12.6); 7 was taken as the window's approximate midpoint, "
                "not fitted to an external target.\n\n")
        f.write("**4. BigMHC-IM and PRIME %Rank disagree on which pathogen is more immunogenic, "
                "compared via their own per-pathogen construct means.** BigMHC-IM (higher=better): "
                f"HIV={p['BIGMHC_HIV']}, Mpox={p['BIGMHC_MPOX']} -> Mpox trails HIV (ratio "
                f"{bigmhc_ratio:.4f}). PRIME %Rank (lower=better): HIV={p['PRIME_HIV']}, "
                f"Mpox={p['PRIME_MPOX']} -- Mpox's raw rank is numerically LOWER (better) than HIV's, "
                f"i.e. PRIME says Mpox EXCEEDS HIV (inverted-ratio {prime_ratio_inverted:.4f}). "
                "Scenario 3 uses BigMHC-IM per this step's brief (the primary TCR-ladder rung, a "
                "direct immunogenicity probability); the disagreement with PRIME is reported as an "
                "explicit limitation on scenario 3's direction, not silently resolved.\n\n")
        f.write(f"**5. GC_output peak window:** {min(gc_days):.2f}-{max(gc_days):.2f} days across all "
                "combinations, inside the spec's stated 10-14 day range by construction (tau chosen "
                "and verified for this).\n")
    print(f"  Manuscript corrections written: {corrections_path}")

    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write("# Phase 5D -- Methodology Note\n\n")
        f.write("**Spec:** Section V.D, `REVISED_METHODOLOGY_MpoxHIV.md` lines 222-258.\n\n")
        f.write("**Ag(t) source:** `Phase5C_FullCurve_IMdepot_*.csv` (never the IV-bolus curve -- "
                "the depot fix's whole purpose was a realistic ~4h absorption delay, and an "
                "instantaneous-spike-driven T-cell model would give a materially different early "
                f"response). Concentration_ug_per_L column used as Ag(t), linearly interpolated "
                f"(`numpy.interp`) at arbitrary t. Primary ka={p['KA_PRIMARY']}/h; sensitivity "
                f"ka={p['KA_SENSITIVITY']}.\n\n")
        f.write(f"**Solver:** `scipy.integrate.solve_ivp`, method=`{p['METHOD']}`, rtol={p['RTOL']}, "
                f"atol={p['ATOL']}. Radau (used in Steps 5B/5C) struggled with this interpolated-"
                "forcing right-hand side, whose derivative has kinks at each dose event and at the "
                "underlying curve's own hourly sample points; LSODA (variable-order, automatically "
                "switches stiff/non-stiff) handled it robustly. This is a solver-choice deviation "
                "from Steps 5B/C, made because the RHS itself changed character (interpolated "
                "external forcing vs a pure closed-form linear system), not a loosening of accuracy "
                "standards -- rtol/atol are unchanged from 5B/5C.\n\n")
        f.write("**State vector:** [T_eff, G1, GC_output, B_mem, Plasma_cell, Ab_titer], all "
                "arbitrary units, all start at 0.\n\n")
        f.write("**T-cell equation, spec ambiguity resolved:** the spec writes "
                "`dT_eff/dt = a*Ag(t) - b*T_eff - g*T_eff*T_reg` where g's own stated units already "
                "carry a `(T_reg/T_total)` factor, and no T_reg(t) trajectory is defined anywhere in "
                "the spec. Rather than invent an unstated T_reg(t) state, b and g were folded into a "
                "single effective per-day suppression rate, `-(b+g)*T_eff`; both are still swept "
                "across their full stated ranges (b: 0.01-0.02/day, g: 0.001-0.005/day), as a 3-point "
                "(low/mid/high) grid.\n\n")
        f.write("**GC_output coupling and B-cell/antibody parameters:** see manuscript-corrections "
                "items 2-3 above (same file, `MANUSCRIPT_CORRECTIONS_5D.md`) for the full reasoning; "
                "not repeated here.\n\n")
        f.write("**Alpha scenarios:** three, per the brief, none picked as a winner. Scenario 2 "
                "(\"BigMHC-weighted\") uses the construct-wide BigMHC-IM mean (0.3502) uniformly "
                "(tests whether weighting by overall immunogenicity strength vs. an arbitrary "
                "uniform 1.0 changes T-cell kinetics shape). Scenario 3 (\"per-pathogen\") uses "
                f"BigMHC-IM's own HIV/Mpox means ({p['BIGMHC_HIV']}/{p['BIGMHC_MPOX']}) -- BigMHC-IM "
                "was chosen over PRIME for this split because it is the TCR ladder's primary rung and "
                "a direct immunogenicity probability, not because it gave an expected-looking answer "
                "(the PRIME-based cross-check, which disagrees in direction, is reported plainly, not "
                "hidden -- see manuscript corrections item 4).\n\n")
        f.write("**Direction checks (re-cited from Phase IV Step 3, not re-derived):** BigMHC-IM "
                "best score vs Phase4B EL rank, Pearson r=-0.2984, p=0.4024, n=10 (expected negative "
                "sign, correct, not significant at n=10). PRIME %Rank vs Phase4B EL rank, Pearson "
                "r=+0.9159, p=0.0002, n=10 (expected positive sign, meaningful). Both tools' scores "
                "were used only after confirming these signs match their stated conventions.\n\n")
        f.write("**Gates, all passed:** no negative state values anywhere; antigen peak <= T_eff "
                "peak <= antibody peak in every one of the "
                f"{len(summary_rows)} (ka, scenario, pathogen, b/g) combinations simulated; "
                f"GC_output's prime-induced peak (before boost 1, day 21) falls in "
                f"[{min(gc_days):.2f}, {max(gc_days):.2f}] days, inside the spec's stated 10-14 day "
                "window for every combination; T_eff and Ab_titer both decay well below their peak "
                "values by end of study (no divergence).\n\n")
        f.write("**Determinism:** no stochastic elements anywhere in this step (deterministic ODEs, "
                "deterministic solver, no RNG seed needed). Re-running reproduces output exactly.\n\n")
        f.write(f"**Figures:** {'written (' + ', '.join(os.path.basename(f) for f in figures_written) + ')' if figures_written else 'CSV-only -- matplotlib figure generation failed, see run log for the exception'}.\n")
    print(f"  Methodology note written: {note_path}")

    ok_post = common.verify_rule1("Step 5D post-check")
    if not ok_post:
        print("[ERROR] Rule 1 violated during Step 5D.")
        sys.exit(2)

    print("\n[SUCCESS] Step 5D complete. Not starting 5E per the brief.")
    return ledger.rows


if __name__ == "__main__":
    run()
