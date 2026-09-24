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

# =============================================================================
# STEP 5B -- two-compartment open PK model (spec lines 183-205)
#
# dA1/dt = -k12*A1 - k_elim*A1 + k21*A2
# dA2/dt =  k12*A1 - k21*A2
#
# A1 = central compartment (plasma), A2 = peripheral (lymphoid/liver/RES).
# Solved numerically via solve_ivp; derived PK parameters (t_half, V_d, CL,
# AUC) are checked against both the analytic solution and the spec's stated
# values (Rule 4: check sign/range/direction before using a result).
# =============================================================================

PARAMS = {
    "K12_PER_H": 0.30,      # central -> peripheral transfer rate
    "K21_PER_H": 0.05,      # peripheral -> central return rate
    "K_ELIM_PER_H": 0.08,   # elimination rate (hepatic RES uptake, mRNA degradation)
    "DOSE_UG": 10.0,        # reference dose for this step
    "VD_RANGE_L": (6.0, 8.0),   # spec's stated V_d range for CL/AUC checks
    "T_SPAN_H": (0.0, 1500.0),  # ~62.5 days -- the SYSTEM'S slow eigenvalue (peripheral-compartment
                                # return, k21=0.05/h) gives a ~73h decay half-life, far longer than
                                # k_elim's 8.66h alone; 1500h is long enough for total mass to clear
                                # the 1e-4*dose floor (checked below, not assumed)
    "N_EVAL": 6001,             # 0.25 h resolution for the trapezoid AUC over the longer window
    "RTOL": 1e-10,
    "ATOL": 1e-12,
    "METHOD": "Radau",      # stiff-capable implicit method; this system is not stiff at these
                             # rate constants, but Radau's accuracy at tight tol costs little here
}


def odes(t, y, k12, k21, k_elim):
    a1, a2 = y
    da1 = -k12 * a1 - k_elim * a1 + k21 * a2
    da2 = k12 * a1 - k21 * a2
    return [da1, da2]


def solve(dose_ug):
    k12, k21, k_elim = PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"]
    t_eval = np.linspace(*PARAMS["T_SPAN_H"], PARAMS["N_EVAL"])
    sol = solve_ivp(
        odes, PARAMS["T_SPAN_H"], y0=[dose_ug, 0.0],
        args=(k12, k21, k_elim), t_eval=t_eval,
        method=PARAMS["METHOD"], rtol=PARAMS["RTOL"], atol=PARAMS["ATOL"],
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"solve_ivp failed: {sol.message}")
    return sol.t, sol.y[0], sol.y[1]


def run():
    common.print_banner("PHASE 5B -- TWO-COMPARTMENT PK MODEL (Section V.B, spec lines 183-205)")

    ok_pre = common.verify_rule1("Step 5B pre-check")
    if not ok_pre:
        print("[ERROR] Rule 1 dirty before Step 5B even ran -- refusing to proceed.")
        sys.exit(2)

    ledger = common.Ledger()
    k12, k21, k_elim = PARAMS["K12_PER_H"], PARAMS["K21_PER_H"], PARAMS["K_ELIM_PER_H"]
    dose = PARAMS["DOSE_UG"]

    print("\n" + "-" * 90)
    print(f"PARAMETERS: k12={k12} /h, k21={k21} /h, k_elim={k_elim} /h, dose={dose} ug, "
          f"solver={PARAMS['METHOD']} rtol={PARAMS['RTOL']} atol={PARAMS['ATOL']}")

    t_half_micro = math.log(2) / k_elim
    print(f"\nMICRO-constant ln(2)/k_elim = {t_half_micro:.4f} h (spec line 202 labels this 'the' half-life, 8.7 h)")
    assert abs(t_half_micro - 8.6643) < 0.01, "micro-constant formula check failed"
    ledger.N("5B_t_half_MICRO_ln2_over_kelim_h", round(t_half_micro, 4), "Phase5B (analytic; NOT the terminal half-life)")

    print("\nCORRECTION (HALT 5B, trigger 2): ln(2)/k_elim is a MICRO-constant. It equals the terminal "
          "half-life only in a one-compartment model. This is a two-compartment system -- the observable "
          "macro-constants (alpha, beta) come from the eigenvalues of the system matrix.")
    M = np.array([[-(k12 + k_elim), k21], [k12, -k21]])
    eig = np.linalg.eigvals(M)
    alpha = float(max(abs(eig)))   # fast (distribution) macro-constant
    beta = float(min(abs(eig)))    # slow (terminal) macro-constant
    t_half_alpha = math.log(2) / alpha
    t_half_beta = math.log(2) / beta
    print(f"  alpha (distribution) = {alpha:.6f} /h -> distribution t_half = {t_half_alpha:.4f} h")
    print(f"  beta  (terminal)     = {beta:.6f} /h -> terminal t_half     = {t_half_beta:.4f} h")
    assert abs(alpha - 0.420487) < 1e-5 and abs(beta - 0.009513) < 1e-5, "macro-constants do not match expected eigenvalues"
    ledger.N("5B_alpha_per_h", round(alpha, 6), "Phase5B (eigenvalues of system matrix)")
    ledger.N("5B_beta_per_h", round(beta, 6), "Phase5B (eigenvalues of system matrix)")
    ledger.N("5B_t_half_distribution_h", round(t_half_alpha, 4), "Phase5B (analytic, alpha)")
    ledger.N("5B_t_half_terminal_h", round(t_half_beta, 4), "Phase5B (analytic, beta)")
    print(f"  [SUCCESS] macro-constants computed: alpha={alpha:.6f}/h (t_half={t_half_alpha:.3f}h), "
          f"beta={beta:.6f}/h (t_half={t_half_beta:.3f}h). The TERMINAL half-life is {t_half_beta:.2f}h, "
          f"not the {t_half_micro:.2f}h micro-constant the spec quotes.")

    # kept as an alias for the summary table / downstream code that expects "t_half" -- now points at
    # the correct terminal value, not the micro-constant
    t_half = t_half_beta

    print("\n" + "-" * 90)
    print("SOLVING the ODE system numerically")
    t, a1, a2 = solve(dose)

    neg_a1 = int(np.sum(a1 < -1e-9))
    neg_a2 = int(np.sum(a2 < -1e-9))
    print(f"  Negative A1 samples: {neg_a1}, negative A2 samples: {neg_a2} (must both be 0)")
    assert neg_a1 == 0 and neg_a2 == 0, "Negative compartment mass -- solver or model error"

    total = a1 + a2
    mass_drift = float(np.max(np.abs(total - dose)) / dose)
    # NOTE: this is an OPEN system (k_elim term), so total mass is NOT conserved --
    # it decays toward 0. "Mass balance" here means: does the analytic elimination-only
    # decay envelope bound the numeric total from above, and does total -> 0 as t -> inf.
    print(f"  Total mass A1+A2 at t=0: {total[0]:.6f} ug (dose: {dose} ug)")
    print(f"  Total mass A1+A2 at t={t[-1]:.0f}h: {total[-1]:.6e} ug (should be ~0, open system)")
    assert abs(total[0] - dose) < 1e-6, "Initial mass does not equal dose"
    assert total[-1] < dose * 1e-4, "System has not decayed sufficiently by end of window -- extend T_SPAN_H"
    print("  [SUCCESS] Initial mass = dose; system decays toward 0 as expected for an open 2-compartment model.")

    print("\n" + "-" * 90)
    print("GATE -- empirical terminal slope of the solved A1(t) curve vs analytic beta")
    # Terminal phase: fit log(A1) vs t on the LAST portion of the curve, where the fast (alpha)
    # component has decayed away and only the slow (beta) mode remains. Using the last 20% of the
    # window (t >= 1200h here) is well past 1/alpha (~2.4h) so the fast mode's contribution is
    # negligible there.
    tail_mask = t >= (PARAMS["T_SPAN_H"][1] * 0.8)
    t_tail, a1_tail = t[tail_mask], a1[tail_mask]
    valid = a1_tail > 0
    slope, _ = np.polyfit(t_tail[valid], np.log(a1_tail[valid]), 1)
    empirical_beta = -slope
    beta_err_pct = abs(empirical_beta - beta) / beta * 100
    print(f"  Empirical terminal slope (linear fit of ln(A1) over t={t_tail[0]:.0f}-{t_tail[-1]:.0f}h): "
          f"{empirical_beta:.6f} /h")
    print(f"  Analytic beta: {beta:.6f} /h -- relative error: {beta_err_pct:.4f}%")
    assert beta_err_pct < 1.0, "Empirical terminal slope does not match analytic beta to within 1% -- model or fit error"
    print(f"  [SUCCESS] The curve's own terminal slope confirms beta={beta:.6f}/h (t_half={t_half_beta:.2f}h) "
          f"as the real terminal decay rate, not the {k_elim} micro-constant.")
    ledger.N("5B_empirical_terminal_slope_per_h", round(empirical_beta, 6), "Phase5B_PKCurve (fit)")
    ledger.N("5B_empirical_vs_analytic_beta_pct_err", round(beta_err_pct, 4), "Phase5B_PKCurve (fit)")

    print("\n" + "-" * 90)
    print("PK PARAMETER DERIVATION")
    for vd_lo, vd_hi in [PARAMS["VD_RANGE_L"]]:
        cl_lo = k_elim * vd_lo
        cl_hi = k_elim * vd_hi
        print(f"  CL = k_elim * V_d: at V_d={vd_lo}-{vd_hi} L -> CL = {cl_lo:.4f}-{cl_hi:.4f} L/h (spec: 0.48-0.64 L/h)")
        assert abs(cl_lo - 0.48) < 0.001 and abs(cl_hi - 0.64) < 0.001, "CL range does not match spec"
        ledger.N("5B_CL_lo_L_per_h", round(cl_lo, 4), "Phase5B (analytic)")
        ledger.N("5B_CL_hi_L_per_h", round(cl_hi, 4), "Phase5B (analytic)")

        auc_analytic_lo = dose / cl_hi   # larger CL -> smaller AUC
        auc_analytic_hi = dose / cl_lo
        print(f"  AUC = Dose/CL: at CL={cl_lo:.4f}-{cl_hi:.4f} L/h, dose={dose} ug -> "
              f"AUC = {auc_analytic_lo:.4f}-{auc_analytic_hi:.4f} ug*h/L (spec: 15.6-20.8)")
        assert abs(auc_analytic_lo - 15.625) < 0.01 and abs(auc_analytic_hi - 20.8333) < 0.01, \
            "AUC analytic range does not match spec"
        ledger.N("5B_AUC_analytic_lo", round(auc_analytic_lo, 4), "Phase5B (analytic)")
        ledger.N("5B_AUC_analytic_hi", round(auc_analytic_hi, 4), "Phase5B (analytic)")

    print("\n" + "-" * 90)
    print("NUMERICAL AUC (trapezoid over the solved A1 curve) vs analytic Dose/CL")
    auc_numeric_a1 = float(np.trapz(a1, t))
    # The "AUC" the spec means is plasma (central-compartment) exposure. CL, as computed above,
    # is defined via k_elim*V_d, an independent quantity from V_d, not directly from this curve's
    # AUC (V_d is a separate stated range, not solved for from A1). To check the solver against
    # the closed-form single-dose-IV-bolus two-compartment identity for compartment 1:
    #   AUC_A1(0-inf) = Dose / k_elim   (standard result: total elimination-only clearance of the
    #   *amount* in compartment 1, since only k_elim removes mass from the whole system and A1
    #   returns to 0; this holds regardless of k12/k21, by mass-balance integration of dA1/dt+dA2/dt).
    auc_analytic_identity = dose / k_elim
    print(f"  Numeric AUC(0-{t[-1]:.0f}h) via trapezoid on A1(t): {auc_numeric_a1:.4f} ug*h/L-equivalent "
          f"(units: A1 is in ug, so this is ug*h; treated as ug*h/L with implicit unit V=1L central "
          f"compartment, matching the spec's own dose-in-ug/AUC-in-ug*h/L convention)")
    print(f"  Analytic identity AUC = Dose/k_elim = {auc_analytic_identity:.4f} ug*h "
          f"(exact for the FULL 0-infinity integral of A1, any k12/k21)")
    rel_err = abs(auc_numeric_a1 - auc_analytic_identity) / auc_analytic_identity
    print(f"  Relative error (numeric truncated-window trapezoid vs analytic 0-inf identity): {rel_err*100:.4f}%")
    assert rel_err < 0.01, "Numerical AUC does not match analytic identity to within 1% -- solver tolerance too loose"
    print("  [SUCCESS] Numerical AUC matches the analytic 0-infinity identity to within 1% "
          f"(solved 0-{PARAMS['T_SPAN_H'][1]:.0f}h is long enough to capture effectively all of the AUC).")
    ledger.N("5B_AUC_numeric_ug_h", round(auc_numeric_a1, 4), "Phase5B_PKCurve")
    ledger.N("5B_AUC_analytic_identity_ug_h", round(auc_analytic_identity, 4), "Phase5B (analytic)")
    ledger.N("5B_AUC_numeric_vs_analytic_pct_err", round(rel_err * 100, 5), "Phase5B (analytic)")

    print("\n" + "-" * 90)
    print("Reconciling the two AUC statements: the spec's Dose/CL AUC (15.6-20.8) assumes CL is an "
          "INDEPENDENT clinical estimate (k_elim * an assumed apparent V_d of 6-8 L), a common PK "
          "convention where V_d is fit from data, not this model's own state variables. This model's "
          "own compartment-1 AUC identity (Dose/k_elim = {:.1f} ug*h) is a DIFFERENT quantity -- it is "
          "not divided by an apparent V_d, so its units are ug*h (mass*time), not ug*h/L "
          "(concentration*time). The spec's 15.6-20.8 ug*h/L figure implicitly assumes C = A1/V_d with "
          "V_d = 6-8 L; applying that here: AUC_conc = (Dose/k_elim)/V_d.".format(auc_analytic_identity))
    for vd in PARAMS["VD_RANGE_L"]:
        auc_conc = auc_analytic_identity / vd
        print(f"    V_d={vd} L -> AUC_conc = {auc_conc:.4f} ug*h/L")
    auc_conc_lo = auc_analytic_identity / PARAMS["VD_RANGE_L"][1]
    auc_conc_hi = auc_analytic_identity / PARAMS["VD_RANGE_L"][0]
    print(f"  This reproduces the spec's 15.6-20.8 ug*h/L range exactly (same Dose/(k_elim*V_d) = Dose/CL "
          f"identity used above): {auc_conc_lo:.4f}-{auc_conc_hi:.4f} ug*h/L.")
    assert abs(auc_conc_lo - auc_analytic_lo) < 1e-9 and abs(auc_conc_hi - auc_analytic_hi) < 1e-9, \
        "AUC-via-concentration does not match AUC-via-Dose/CL -- these must be algebraically identical"
    print("  [SUCCESS] Confirmed algebraically identical to the Dose/CL calculation above -- one formula, "
          "consistently applied.")

    out_dir = common.step_output_dir("B")
    ts = common.timestamp()

    curve_path = os.path.join(out_dir, f"Phase5B_PKCurve_{ts}.csv")
    with open(curve_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Time_h", "A1_ug", "A2_ug", "Total_ug"])
        for ti, a1i, a2i in zip(t, a1, a2):
            w.writerow([round(float(ti), 4), round(float(a1i), 8), round(float(a2i), 8), round(float(a1i + a2i), 8)])
    print(f"\n  PK curve written: {curve_path} ({len(t)} time points)")

    summary_path = os.path.join(out_dir, f"Phase5B_Summary_{ts}.csv")
    with open(summary_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Metric", "Value", "Unit", "Spec_value", "Match"])
        w.writerow(["k12", k12, "/h", 0.30, "YES"])
        w.writerow(["k21", k21, "/h", 0.05, "YES"])
        w.writerow(["k_elim", k_elim, "/h", 0.08, "YES"])
        w.writerow(["t_half_MICRO_ln2_over_kelim", round(t_half_micro, 4), "h", "8.7",
                    "matches spec's number, but this is NOT the terminal half-life (see below)"])
        w.writerow(["alpha_macro_constant", round(alpha, 6), "/h", "n/a (not in spec)", "n/a"])
        w.writerow(["beta_macro_constant", round(beta, 6), "/h", "n/a (not in spec)", "n/a"])
        w.writerow(["t_half_distribution_alpha", round(t_half_alpha, 4), "h", "n/a", "n/a"])
        w.writerow(["t_half_TERMINAL_beta", round(t_half_beta, 4), "h", "8.7 (spec, WRONG)",
                    "NO -- spec's 8.7h is the micro-constant, not this terminal value"])
        w.writerow(["empirical_terminal_slope", round(empirical_beta, 6), "/h", "n/a",
                    f"matches analytic beta to {beta_err_pct:.3f}%"])
        w.writerow(["CL_lo", round(cl_lo, 4), "L/h", 0.48, "YES"])
        w.writerow(["CL_hi", round(cl_hi, 4), "L/h", 0.64, "YES"])
        w.writerow(["AUC_lo", round(auc_analytic_lo, 4), "ug*h/L", 15.6, "YES"])
        w.writerow(["AUC_hi", round(auc_analytic_hi, 4), "ug*h/L", 20.8, "YES"])
        w.writerow(["AUC_numeric_vs_analytic_pct_err", round(rel_err * 100, 5), "%", "<1", "YES"])
        w.writerow(["Dose_reference", dose, "ug", 10, "YES"])
    print(f"  Summary written: {summary_path}")

    ledger_path = os.path.join(out_dir, f"Phase5B_Ledger_{ts}.csv")
    ledger.write(ledger_path)
    print(f"  Ledger written: {ledger_path} ({len(ledger.rows)} numbers)")

    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write("# Phase 5B -- Methodology Note\n\n")
        f.write("**Spec:** Section V.B, `REVISED_METHODOLOGY_MpoxHIV.md` lines 183-205.\n\n")
        f.write("**Model:** two-compartment open PK, `dA1/dt = -k12*A1 - k_elim*A1 + k21*A2`, "
                "`dA2/dt = k12*A1 - k21*A2`. Solved via `scipy.integrate.solve_ivp`, method=`{}`, "
                "rtol={}, atol={}, over t=0-{}h at {} points.\n\n".format(
                    PARAMS["METHOD"], PARAMS["RTOL"], PARAMS["ATOL"], PARAMS["T_SPAN_H"][1], PARAMS["N_EVAL"]))
        f.write("**CORRECTION (post-review, HALT 5B / trigger 2): the spec's half-life is a category "
                "error.** `ln(2)/k_elim = 8.66h` is a MICRO-constant. It equals the terminal half-life "
                "only in a one-compartment model; this is a two-compartment system, so the *observable* "
                "half-lives come from the eigenvalues of the 2x2 system matrix, not from k_elim alone.\n\n")
        f.write(f"- alpha (distribution, fast) = {alpha:.6f} /h -> t_half = {t_half_alpha:.3f} h\n")
        f.write(f"- beta (terminal, slow) = {beta:.6f} /h -> t_half = {t_half_beta:.2f} h\n\n")
        f.write(f"The empirical terminal slope of the solved A1(t) curve ({empirical_beta:.6f} /h, fit over "
                f"the last 20% of the {PARAMS['T_SPAN_H'][1]:.0f}h window) confirms beta to "
                f"{beta_err_pct:.3f}% -- the curve was always right, only the reported single-number "
                "summary was wrong. Both half-lives, both macro-constants, and the micro-constant are now "
                "reported, each explicitly labelled, and the empirical-vs-analytic beta check is a "
                "permanent verification gate.\n\n")
        f.write(f"**CL and AUC are unaffected** -- both are defined via k_elim directly "
                "(`CL = k_elim*V_d`, `AUC = Dose/CL`), which is correct terminology in either a one- or "
                "two-compartment model; they still reproduce the spec exactly (0.48-0.64 L/h, "
                "15.6-20.8 ug*h/L).\n\n")
        f.write("**Consequence flagged for Step 5C:** with beta=0.0095/h (t_half_terminal~73h), the "
                "prime-dose peripheral reservoir has NOT washed out by day 21 (504h, ~6.9 terminal "
                "half-lives puts ~1% of the dose remaining, not zero). Step 5C must therefore show "
                "non-zero troughs and a NON-DECREASING C_max across doses, not the spec's table showing "
                "C_max falling with zero troughs.\n\n")
        f.write("**AUC derivation, made explicit (the spec shows no derivation for this step):** two "
                "equivalent routes were checked against each other. (1) `AUC = Dose/CL` where "
                "`CL = k_elim*V_d`, treating V_d (6-8 L) as an independently stated apparent volume "
                "of distribution, the standard clinical-PK convention. (2) This model's own compartment-1 "
                "mass-integral identity, `AUC_mass(0-inf) = Dose/k_elim` (exact for any k12/k21, since only "
                "k_elim removes mass from the closed A1+A2 system), converted to a concentration-AUC by "
                "dividing by the same V_d. Both routes give the identical 15.6-20.8 ug*h/L range -- "
                "confirmed algebraically in the run log, not just numerically.\n\n")
        f.write(f"**Numerical vs analytic check:** trapezoid integration of the solved A1(t) curve over "
                f"0-{PARAMS['T_SPAN_H'][1]:.0f}h gives {auc_numeric_a1:.4f} ug*h, matching the analytic "
                f"Dose/k_elim identity ({auc_analytic_identity:.4f} ug*h) to {rel_err*100:.4f}% -- well under "
                "the 1% gate. The solve window (240h, ~28 half-lives) is long enough that the truncated "
                "trapezoid integral is effectively the full 0-infinity integral.\n\n")
        f.write("**Checks passed:** no negative compartment masses at any solved time point; initial "
                "A1+A2 equals the dose exactly; A1+A2 decays toward 0 (open system, not conserved) rather "
                "than plateauing or growing.\n\n")
        f.write("**Escalation: trigger 2 fired and was resolved by direct ruling (HALT 5B), not by an "
                "independent Opus escalation cycle** -- the reviewer identified the half-life category "
                "error directly, with the exact macro-constants, before this step's own report was "
                "written. CL and AUC remain within spec rounding; only the half-life claim was wrong.\n")
    print(f"  Methodology note written: {note_path}")

    ok_post = common.verify_rule1("Step 5B post-check")
    if not ok_post:
        print("[ERROR] Rule 1 violated during Step 5B.")
        sys.exit(2)

    print("\n[SUCCESS] Step 5B complete.")
    return ledger.rows


if __name__ == "__main__":
    run()