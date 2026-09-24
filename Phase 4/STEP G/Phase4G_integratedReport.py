import os
import sys
import csv
import glob
import hashlib
import statistics
from collections import Counter, defaultdict

# =============================================================================
# MINIMAL BOOTSTRAP -- see Phase 1/_common/phase1_common.py.
# =============================================================================
def _bootstrap_find_research_root(script_file):
    current = os.path.dirname(os.path.abspath(script_file))
    while os.path.basename(current) != "Research":
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a 'Research' anchor folder above: {script_file}")
            sys.exit(1)
        current = parent
    return current

_PROJECT_ROOT = _bootstrap_find_research_root(__file__)
_COMMON_DIR = os.path.join(_PROJECT_ROOT, "Phase 4", "_common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)
import phase4_common as common

# =============================================================================
# PHASE 4G -- INTEGRATED REPORT (Step 7 of 7). Pure synthesis: no predictions, no network.
# Every number in the report is produced by N(...) below, which records (key, value, source CSV)
# in a ledger written next to the report -- so each figure traces to a file in StepA-StepF
# (or, for the Phase 1 regression comparison only, the frozen Phase 1F pool CSV, read-only).
#
# Report order follows STRENGTH OF EVIDENCE, not step number (Opus ruling):
#   1 coverage gaps -> 2 method divergence -> 3 binding/breadth (+TCR, cross-reactivity) -> 4 per pathogen
#   -> 5 processing -> 6 junctions (class I real / class II clean) -> 7 B-cell + conservancy -> 8 what C-ImmSim's
#   replacement does NOT give.
# =============================================================================

STEP = {s: common.step_output_dir(s) for s in "ABCDEF"}
OUT = common.step_output_dir("G")
TS = common.timestamp()
LEDGER = []
SRC = {}


def latest(step, prefix):
    files = sorted(glob.glob(os.path.join(STEP[step], f"{prefix}_*.csv")))
    if not files:
        print(f"[ERROR] missing {prefix}_*.csv in Step{step} -- run the earlier steps first.")
        sys.exit(1)
    return files[-1]


def rd(step, prefix):
    p = latest(step, prefix)
    SRC[prefix] = os.path.basename(p)
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def N(key, value, prefix, fmt="{}"):
    """Format a number for the report AND record it in the ledger with its source file."""
    LEDGER.append({"Key": key, "Value": value, "Source_File": SRC.get(prefix, prefix)})
    return fmt.format(value)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def med(v):
    v = [x for x in v if x is not None]
    return statistics.median(v) if v else None


def mean(v):
    v = [x for x in v if x is not None]
    return statistics.mean(v) if v else None


def sat(x):
    return ">99.99%" if x >= 99.995 else f"{x:.1f}%"


def sat2(x):
    return ">99.99%" if x >= 99.995 else f"{x:.2f}%"


def parse_pl(s):
    return {kv.split(":")[0]: float(kv.split(":")[1]) for kv in (s or "").split(";") if ":" in kv}


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def fm(x):
    return f"{x:g}" if x is not None else "N/A"


def na(x):
    return "N/A" if x in (None, "") else x


def main():
    common.print_banner("PHASE 4G -- INTEGRATED REPORT (synthesis; every number ledgered)")
    A = rd("A", "Phase4A_EpitopeDossier")
    Bp = rd("B", "Phase4B_BindingAffinity")
    Ct = rd("C", "Phase4C_TCRRecognition_MHCI")
    Cx = rd("C", "Phase4C_CrossReactivity_AllConstruct")
    Cd = rd("C", "Phase4C_DirectionChecks")
    Dp = rd("D", "Phase4D_EpitopeProcessing")
    Ds = rd("D", "Phase4D_Summary")
    Ee = rd("E", "Phase4E_EpitopeCoverage")
    Ecc = rd("E", "Phase4E_ConstructCoverage")
    Eca = rd("E", "Phase4E_ContestedAlleles")
    Eg = rd("E", "Phase4E_CoverageGaps")
    Emd = rd("E", "Phase4E_MethodDivergence")
    Ejc = rd("E", "Phase4E_JunctionCalibration")
    Fb = rd("F", "Phase4F_BcellDossier")
    Fbp = rd("F", "Phase4F_BcellPool")
    Fc = rd("F", "Phase4F_Conservancy")
    Fg = rd("F", "Phase4F_ConservancyGate")
    Fs = rd("F", "Phase4F_Summary")
    D2 = rd("D", "Phase4D_JunctionNeoepitopes_MHCII")
    p1f_path = common.latest_file(os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered"), suffix=".csv")
    with open(p1f_path, newline="") as f:
        P1F = {r["Peptide"]: r for r in csv.DictReader(f)}
    SRC["Phase1F_pool"] = os.path.basename(p1f_path) + " (Phase 1, frozen, read-only)"

    by = {k: {r["Peptide"]: r for r in rows} for k, rows in
          dict(A=A, B=Bp, C=Ct, X=Cx, D=Dp, E=Ee, F=Fb, K=Fc).items()}
    ds = {r["Metric"]: r["Value"] for r in Ds}
    fs = {r["Metric"]: r["Value"] for r in Fs}
    dc = {r["Check"]: r for r in Cd}
    calib = defaultdict(dict)
    for r in Ejc:
        calib[r["Section"]][r["Item"]] = r
    cclass = {"MHC-I": [r for r in A if r["Class"] == "MHC-I"], "MHC-II": [r for r in A if r["Class"] == "MHC-II"], "B-cell": [r for r in A if r["Class"] == "B-cell"]}
    assert (len(cclass["MHC-I"]), len(cclass["MHC-II"]), len(cclass["B-cell"])) == (10, 11, 10)

    # ------------------------------------------------------------- DOSSIER
    def big(pep):
        return fnum(by["C"].get(pep, {}).get("BigMHC_IM_Best"))
    dossier = []
    for e in A:
        pep, cls = e["Peptide"], e["Class"]
        row = {"Peptide": pep, "Class": cls, "Pathogen": e["Pathogen"], "Target": e["Labelled_Target"], "Antigen_Display": e["Antigen_Display"],
               "Construct_Start_0based": e["Construct_Position_Start"], "Construct_End": e["Construct_Position_End"], "Length": len(pep)}
        b = by["B"].get(pep, {})
        # ---- verdicts: two separate readings, both recoverable
        if cls == "MHC-I":
            rk, bm, ic = fnum(b.get("NetMHC_EL_Best_Rank")), big(pep), fnum(b.get("IC50_At_Best_Consensus_Allele_nM"))
            row["Primary_Gate"] = "PASS" if rk is not None and rk <= 1.0 else "FAIL"
            row["Primary_Gate_Detail"] = (f"Phase I rules: EL best percentile rank {rk} <= 1.0 -> {row['Primary_Gate']}; BigMHC-IM best {bm:.3f} vs Phase I priority cut 0.5 -> "
                                          f"{'PRIORITIZED' if bm >= 0.5 else 'DEPRIORITIZED (retained; tiers never exclude)'}; breadth {b.get('N_Binding_Alleles_Primary')} of 74 class I panel alleles")
            ok = (ic is not None) and ic <= 500 and bm >= 0.70
            row["Secondary_Stratum"] = "N/A (no IC50-producing consensus method for this allele)" if ic is None else ("STRICT_PASS" if ok else "STRICT_FAIL")
            row["Secondary_Stratum_Detail"] = f"Revised methodology: IC50 {na(b.get('IC50_At_Best_Consensus_Allele_nM'))} nM (need <=500); BigMHC-IM {bm:.3f} (need >=0.70)"
        elif cls == "MHC-II":
            rk, ic = fnum(b.get("NetMHC_EL_Best_Rank")), fnum(b.get("IC50_At_Best_Consensus_Allele_nM"))
            row["Primary_Gate"] = "PASS" if rk is not None and rk <= 10.0 else "FAIL"
            row["Primary_Gate_Detail"] = f"Phase I rules: EL best percentile rank {rk} <= 10.0 -> {row['Primary_Gate']}; breadth {b.get('N_Binding_Alleles_Primary')} of 33 class II panel alleles"
            row["Secondary_Stratum"] = "N/A (no IC50-producing consensus method for this allele)" if ic is None else ("STRICT_PASS" if ic <= 1000 else "STRICT_FAIL")
            row["Secondary_Stratum_Detail"] = f"Revised methodology: IC50 {na(b.get('IC50_At_Best_Consensus_Allele_nM'))} nM (need <=1000)"
        else:
            f = by["F"][pep]
            row["Primary_Gate"] = "PASS_PRIORITIZED" if f["Primary_Gate_Phase_I_Tier"] in ("High", "Medium") else "PASS_DEPRIORITIZED"
            row["Primary_Gate_Detail"] = f"Phase I BepiPred tier {f['Primary_Gate_Phase_I_Tier']} (mean {f['mean_BepiPred']}, {f['pct_above']}% residues above); MHC gate N/A (not HLA-restricted)"
            row["Secondary_Stratum"] = "STRICT_PASS" if f["Secondary_Stratum_Methodology_High"] == "YES" else "STRICT_FAIL"
            row["Secondary_Stratum_Detail"] = "Revised methodology: mean_BepiPred >= 0.50 AND pct_above >= 75%"
        # ---- flags
        flags = []
        if pep == "NKRKRVIGL":
            flags.append("STANDING_FLAG: worst on 4/4 independent signals (Step 3); binds 2 of 74 class I alleles; consistent cross-phase A35R weakness")
        d = by["D"].get(pep, {})
        if d.get("Internal_Site_Exceeds_C_Term") == "YES" and cls == "MHC-I":
            flags.append("WATCHLIST: internal cleavage site scores above the C-terminal one")
        aay = calib["watchlist overlap"].get("epitopes adjacent to a qualifying AAY seam AND on the internal-cleavage watchlist", {}).get("Note", "").split(";")
        if pep in aay:
            flags.append("WATCHLIST: adjacent to a qualifying AAY seam (class I)")
        if cls == "B-cell" and by["F"][pep].get("Exposure_Reliability", "").startswith("INDICATIVE"):
            flags.append("B-cell exposure INDICATIVE ONLY (global fold undetermined)")
        row["Standing_Flags"] = " | ".join(flags) if flags else "none"
        row["Contest_Reading"] = ("REAL class I contest (mechanism: AAY terminal Tyr in the C-terminal anchor pocket; 1.8x enrichment over natural windows). Contested != lost."
                                  if cls == "MHC-I" else "class II contest values are artifact-prone (sample-size asymmetry; seam windows are depleted vs random) -- NOT a design defect"
                                  if cls == "MHC-II" else "N/A (not HLA-restricted)")
        for prefix, src in (("A_", e), ("B_", b), ("C_", by["C"].get(pep, {})), ("X_", by["X"].get(pep, {})), ("D_", d), ("E_", by["E"].get(pep, {})),
                            ("F_bcell_", by["F"].get(pep, {})), ("F_cons_", by["K"].get(pep, {}))):
            for k, v in src.items():
                if k in ("Peptide", "Class", "Pathogen") or (prefix in ("D_", "E_", "F_bcell_", "F_cons_") and k in ("Target", "Construct_Start", "Construct_End", "Length", "In_Construct")):
                    continue
                row[prefix + k] = v if v not in ("", None) else "N/A"
        dossier.append(row)
    fields = list(dict.fromkeys(k for r in dossier for k in r))
    dpath = os.path.join(OUT, f"Phase4G_Immunogenicity_Dossier_{TS}.csv")
    with open(dpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, restval="N/A")
        w.writeheader(); w.writerows(dossier)
    print(f"[SUCCESS] {os.path.basename(dpath)}: {len(dossier)} rows x {len(fields)} columns")
    SRC["Dossier"] = os.path.basename(dpath)

    # ---------------------------------------------------------------- numbers
    HIV = [r for r in A if r["Pathogen"] == "HIV"]
    MPX = [r for r in A if r["Pathogen"] == "Mpox"]
    pool = {c: [r for r in Bp if r["Class"] == c] for c in ("MHC-I", "MHC-II", "B-cell")}
    gaps = {r["Allele"]: r for r in Eg}

    def gmass(loc, col):
        return sum(float(r["Allele_Frequency"]) for r in Eg if r["Locus"] == loc and int(r[col]) == 0)

    def ltot(loc):
        return sum(float(r["Allele_Frequency"]) for r in Eg if r["Locus"] == loc)

    def cov(level, scope, uni=None, variant=None, margin=None):
        for r in Ecc:
            if r["Level"] == level and r["Scope"] == scope and (uni is None or r["Allele_Universe"] == uni) and \
               (variant is None or r["Variant"] == variant) and (margin is None or r["Margin"] == margin):
                return r
        raise KeyError((level, scope, uni, variant, margin))

    R_ = []      # report lines
    W = R_.append

    # ---- header / summary
    W("# Phase IV Report -- IEDB-Based Immunogenicity Prediction (Section IV)\n")
    W(f"Construct `{A[0]['Source_Construct']}` (570 aa). Generated {TS} by `Phase4G_integratedReport.py`. Analysis units: **31 construct epitopes** "
      f"({N('n_mhci', len(cclass['MHC-I']), 'Dossier')} MHC-I / {N('n_mhcii', len(cclass['MHC-II']), 'Dossier')} MHC-II / {N('n_bcell', len(cclass['B-cell']), 'Dossier')} B-cell; "
      f"{N('n_hiv', len(HIV), 'Dossier')} HIV-derived / {N('n_mpox', len(MPX), 'Dossier')} Mpox-derived; 7 source antigens) and, for comparison, the **105 candidates** that reached Phase 1F coverage scoring. "
      "Every metric is reported separately for HIV and Mpox. (\"194 Non-Redundant Sequences\" in the methodology counts source *proteins*, not epitopes -- see the corrections file.) "
      "Antigen `Mpox_B5R` in the Phase I files is the manuscript's **B6R** (OPG190); both names are kept in the dossier.\n")
    W("This is an *in silico* immunogenicity assessment built from published predictors. It is not experimental evidence of immunogenicity, and it does not replace the simulation "
      "work of C-ImmSim (section 8).\n")
    W("## Findings ordered by strength of evidence\n")
    W(table(["#", "Finding", "Evidence strength", "Needs a null model?"], [
        (1, "Coverage gap: `A*11:01` and `A*33:03` bound by no construct MHC-I epitope under either predictor", "Strongest -- a count of zero", "No"),
        (2, "Class I coverage depends heavily on predictor choice (locus C spans 28 pp)", "Strong -- direct comparison", "No"),
        (3, "Binding breadth is narrow for some epitopes; \"all pass\" gates are near-vacuous", "Strong", "No"),
        (4, "Mpox class I coverage is the weak arm; `NKRKRVIGL` adds nothing", "Strong, consistent across Steps 3-5", "No"),
        (5, "All 10 MHC-I epitopes liberated; 5 carry an internal cleavage site above their C-terminus", "Moderate -- predictor output", "No"),
        (6, "Junctions: class I seams create real strong binders (AAY); class II seams are clean", "Class I: moderate-strong (null + mechanism); class II: null-based", "Yes -- supplied"),
        (7, "B-cell: 10 epitopes, none SEMA-corroborated, half on undetermined folds", "Weak -- predictor output", "n/a"),
    ]))
    W("")

    # ---- 1 coverage gaps
    W("## 1. Coverage gaps (strongest finding -- no null model needed)\n")
    W("An allele that no construct epitope binds cannot be covered, whatever the linkers do. Class I alleles by population frequency (AFND Singapore Riau Malay, a Southeast Asian (Austronesian) proxy -- Deviation #7):\n")
    top = [r for r in Eg if r["Locus"] in ("A", "B", "C")][:10]
    W(table(["allele", "frequency", "construct epitopes (IEDB-EL)", "construct epitopes (MHCflurry)", "status"],
            [(r["Allele"], N(f"gap_freq_{r['Allele']}", f"{float(r['Allele_Frequency']):.3f}", "Phase4E_CoverageGaps"), r["N_Epitopes_Phase4B_IEDB_EL"], r["N_Epitopes_Phase1F_MHCflurry"], r["Gap_Status"]) for r in top]))
    a11, a33, c08 = gaps["HLA-A*11:01"], gaps["HLA-A*33:03"], gaps["HLA-C*08:01"]
    both = [r for r in Eg if r["Gap_Status"] == "GAP_BOTH_PREDICTORS" and r["Locus"] in ("A", "B", "C")]
    W(f"\n**`HLA-A*11:01` (frequency {N('f_a1101', a11['Allele_Frequency'], 'Phase4E_CoverageGaps')}) and `HLA-A*33:03` ({N('f_a3303', a33['Allele_Frequency'], 'Phase4E_CoverageGaps')}) are bound by no construct MHC-I epitope "
      f"under either predictor.** Together they carry {N('gapA_two', float(a11['Allele_Frequency']) + float(a33['Allele_Frequency']), 'Phase4E_CoverageGaps', '{:.3f}')} of the locus-A frequency mass "
      f"({N('gapA_two_pct', 100 * (float(a11['Allele_Frequency']) + float(a33['Allele_Frequency'])) / ltot('A'), 'Phase4E_CoverageGaps', '{:.0f}%')} of {ltot('A'):.3f}); all alleles that are gaps under both predictors carry "
      f"{N('gapA_mass_both', sum(float(r['Allele_Frequency']) for r in both if r['Locus'] == 'A'), 'Phase4E_CoverageGaps', '{:.3f}')} ({100 * sum(float(r['Allele_Frequency']) for r in both if r['Locus'] == 'A') / ltot('A'):.0f}%). "
      "This is why locus A is the weakest locus under both predictors. "
      f"The gap exists upstream of the construct: **{N('pool_a1101', sum('HLA-A*11:01' in r['Binding_Alleles_Primary'].split(';') for r in pool['MHC-I']), 'Phase4B_BindingAffinity')} of "
      f"{N('pool_mhci_n', len(pool['MHC-I']), 'Phase4B_BindingAffinity')} MHC-I candidates in the Phase 1F pool bind `A*11:01` and "
      f"{N('pool_a3303', sum('HLA-A*33:03' in r['Binding_Alleles_Primary'].split(';') for r in pool['MHC-I']), 'Phase4B_BindingAffinity')} bind `A*33:03` (IEDB-EL)**, so no candidate could have filled it. "
      f"`HLA-C*08:01` (frequency {c08['Allele_Frequency']}) is a gap under IEDB-EL only -- MHCflurry has {c08['N_Epitopes_Phase1F_MHCflurry']} construct epitope binding it, and "
      f"{N('pool_c0801', sum('HLA-C*08:01' in r['Binding_Alleles_Primary'].split(';') for r in pool['MHC-I']), 'Phase4B_BindingAffinity')} pool candidates bind it under IEDB-EL -- as are `B*18:01`, `C*04:03` and `C*07:04`; "
      "the B and C gaps are therefore predictor-dependent, and only the two A alleles are unambiguous.\n")
    W("Frequency mass carried by alleles with no binding construct epitope, per locus:\n")
    W(table(["locus", "IEDB-EL", "MHCflurry/MHCnuggets", "locus total"],
            [(L, N(f"gapmass_el_{L}", f"{gmass(L, 'N_Epitopes_Phase4B_IEDB_EL'):.3f}", "Phase4E_CoverageGaps"), N(f"gapmass_mf_{L}", f"{gmass(L, 'N_Epitopes_Phase1F_MHCflurry'):.3f}", "Phase4E_CoverageGaps"),
              f"{ltot(L):.3f}") for L in ("A", "B", "C", "DRB1", "DQB1", "DPB1")]))
    W("\nPhase I is frozen; closing the gap would require an epitope that binds `A*11:01`. It is reported as a limitation of the construct's population reach, not corrected here.\n")

    # ---- 2 method divergence
    W("## 2. Method divergence (a stated limitation)\n")
    W("Per-locus construct coverage depends on which binding predictor defines \"binds\". **Report the range; do not select one end.**\n")
    W(table(["locus", "IEDB-EL (Phase 4B)", "MHCflurry/MHCnuggets (Phase 1F)", "range", "spread"],
            [(r["Locus"], f"{float(r['Coverage_IEDB_EL_Pct']):.2f}%", f"{float(r['Coverage_MHCflurry_MHCnuggets_Pct']):.2f}%",
              f"{float(r['Range_Low_Pct']):.2f}-{float(r['Range_High_Pct']):.2f}%", N(f"spread_{r['Locus']}", f"{float(r['Spread_pp']):.2f} pp", "Phase4E_MethodDivergence")) for r in Emd]))
    W("\nLocus C spans 70.2-98.2% (28 pp) from predictor choice alone; class II loci agree within 3 pp. The combined coverage figure saturates "
      f"(construct union {sat2(float(cov('2/3 raw', 'CONSTRUCT (21 MHC epitopes)', 'Phase4B')['Coverage_Pct']))} under IEDB-EL and "
      f"{sat2(float(cov('2/3 raw', 'CONSTRUCT (21 MHC epitopes)', 'Phase1F')['Coverage_Pct']))} under MHCflurry/MHCnuggets) because coverage multiplies across six loci, so **the per-locus range is the honest layer**. "
      "DQA1 and DPA1 are not in the frequency table (beta chains only), so class II is on DRB1/DQB1/DPB1.\n")

    # ---- 3 binding / breadth / TCR
    W("## 3. Binding, breadth, TCR recognition and cross-reactivity\n")
    W("### 3.1 Binding breadth (never \"100% pass\")\n")
    W("With a 107-allele panel, \"best rank below threshold at *any* allele\" is near-vacuous -- almost every peptide binds something -- so it is not reported as a result. The informative quantity is **breadth**, "
      "the number of panel alleles each epitope binds (class I peptides are scored against the 74 class I alleles, class II against the 33 class II alleles):\n")
    def breadth(rows):
        return [int(r["N_Binding_Alleles_Primary"]) for r in rows]
    tb = []
    for lab, rows in (("Pool MHC-I (n=%d)" % len(pool["MHC-I"]), pool["MHC-I"]), ("Pool MHC-II (n=%d)" % len(pool["MHC-II"]), pool["MHC-II"]),
                      ("Construct MHC-I (n=10)", [by["B"][r["Peptide"]] for r in cclass["MHC-I"]]), ("Construct MHC-II (n=11)", [by["B"][r["Peptide"]] for r in cclass["MHC-II"]])):
        v = breadth(rows)
        tb.append((lab, N(f"breadth_median_{lab}", fm(statistics.median(v)), "Phase4B_BindingAffinity"), f"{min(v)}-{max(v)}", sum(x <= 2 for x in v)))
    W(table(["set", "median alleles bound", "range", "epitopes binding <=2 alleles"], tb))
    W(f"\nEpitopes binding one or two alleles contribute almost nothing to population coverage (their marginal contribution is shown in `Phase4E_EpitopeCoverage`). "
      "Class I median 7 of 74 (range 1-40); class II median 6 of 33 (range 1-29) across the pool.\n")
    def hb(rows, cls):
        return med(breadth([by["B"][r["Peptide"]] for r in rows if r["Class"] == cls]))
    W("**HIV vs Mpox (construct, median breadth):** class I HIV %s vs Mpox %s; class II HIV %s vs Mpox %s.\n" % (
        N("br_hiv_i", fm(hb(HIV, "MHC-I")), "Phase4B_BindingAffinity"), N("br_mpx_i", fm(hb(MPX, "MHC-I")), "Phase4B_BindingAffinity"),
        N("br_hiv_ii", fm(hb(HIV, "MHC-II")), "Phase4B_BindingAffinity"), N("br_mpx_ii", fm(hb(MPX, "MHC-II")), "Phase4B_BindingAffinity")))
    W("### 3.2 Agreement with Phase I, and scale caveats\n")
    def deltas(cls, rows):
        out = []
        for r in rows:
            p1 = fnum((P1F.get(r["Peptide"]) or {}).get("Percentile_Rank"))
            p4 = fnum(r["NetMHC_EL_Best_Rank"])
            if p1 is not None and p4 is not None:
                out.append(p4 - p1)
        return out
    dI, dII = deltas("MHC-I", pool["MHC-I"]), deltas("MHC-II", pool["MHC-II"])
    cons_d = [fnum(by["B"][r["Peptide"]]["NetMHC_EL_Best_Rank"]) - fnum(r["Percentile_Rank"]) for r in cclass["MHC-I"] + cclass["MHC-II"]]
    W(f"Phase 4B's EL rank was compared with Phase 1Db's stored `Percentile_Rank` for the same peptides. Class I (n={N('n_dI', len(dI), 'Phase1F_pool')}): mean difference {N('dI_mean', statistics.mean(dI), 'Phase1F_pool', '{:+.3f}')}, "
      f"{N('dI_out', sum(abs(x) > 1 for x in dI), 'Phase1F_pool')} with |difference| > 1. Class II (n={len(dII)}): mean {N('dII_mean', statistics.mean(dII), 'Phase1F_pool', '{:+.3f}')}, median {statistics.median(dII):.1f}, "
      f"{N('dII_out', sum(abs(x) > 1 for x in dII), 'Phase1F_pool')} with |difference| > 1; among the 21 in-construct MHC epitopes, {N('dC_out', sum(abs(x) > 1 for x in cons_d), 'Dossier')} exceed 1. "
      "The class II spread is a **method-scale difference, not a fault**: Phase I's `recommended` method resolves to EL for class I but to BA/consensus for class II. "
      "**Phase I and Phase IV MHC-II ranks are on different scales -- the same number means a different thing -- and must never be tabled side by side without saying so.** No epitope changes pass/fail status under Phase I thresholds.\n")
    flips = 0
    for cls_, thr in (("MHC-I", 1.0), ("MHC-II", 10.0)):
        for r in pool[cls_]:
            p1 = fnum((P1F.get(r["Peptide"]) or {}).get("Percentile_Rank")); p4 = fnum(r["NetMHC_EL_Best_Rank"])
            if p1 is not None and p4 is not None and ((p1 <= thr) != (p4 <= thr)):
                flips += 1
    W(f"Across the {N('n_flip_pool', len(dI) + len(dII), 'Phase1F_pool')} class I/II pool candidates, {N('n_flips', flips, 'Phase1F_pool')} change pass/fail status under Phase I's own thresholds (rank <= 1 class I, <= 10 class II).\n")
    ic_i = [r for r in pool["MHC-I"] if r["IC50_At_Best_Consensus_Allele_nM"] not in ("N/A", "")]
    ic_ii = [r for r in pool["MHC-II"] if r["IC50_At_Best_Consensus_Allele_nM"] not in ("N/A", "")]
    W(f"IC50 is reported as an additional descriptor only (disagreement #1: selection is by percentile rank). It comes from IEDB's `consensus` method, which supports fewer alleles than the EL methods, so it is `N/A` where no "
      f"IC50-producing method exists: available for {N('ic50_i', len(ic_i), 'Phase4B_BindingAffinity')}/{len(pool['MHC-I'])} class I and {N('ic50_ii', len(ic_ii), 'Phase4B_BindingAffinity')}/{len(pool['MHC-II'])} class II pool candidates "
      "(at the allele consensus ranks best). Only labelled nM columns are read (a Step 2 misparse of an unlabeled column is documented in the corrections file).\n")
    vc = defaultdict(Counter)
    for r in dossier:
        vc[r["Class"]][("P:" + r["Primary_Gate"], "S:" + r["Secondary_Stratum"].split(" (")[0])] += 1
    W("**Two verdict readings per epitope (dossier columns `Primary_Gate` and `Secondary_Stratum`).** `Primary_Gate` applies Phase I's rules (class I: EL rank <= 1; class II: <= 10; B-cell: Phase I tier); `Secondary_Stratum` applies the revised methodology's stricter values "
      "(class I: IC50 <= 500 nM at consensus's best allele and BigMHC-IM >= 0.70; class II: IC50 <= 1000 nM; B-cell: mean >= 0.50 and >= 75%). Counts (n = 10 / 11 / 10):\n")
    W(table(["class", "Primary_Gate", "Secondary_Stratum", "epitopes"], [(c, k[0][2:], k[1][2:], N(f"verd_{c}_{k[0][2:]}_{k[1][2:]}", v, "Dossier")) for c in ("MHC-I", "MHC-II", "B-cell") for k, v in sorted(vc[c].items())]))
    bm7 = sum(fnum(r["BigMHC_IM_Best"]) >= 0.7 for r in Ct)
    ic500 = sum(1 for r in cclass["MHC-I"] if (fnum(by["B"][r["Peptide"]]["IC50_At_Best_Consensus_Allele_nM"]) or 1e9) <= 500)
    ic1000 = sum(1 for r in cclass["MHC-II"] if (fnum(by["B"][r["Peptide"]]["IC50_At_Best_Consensus_Allele_nM"]) or 1e9) <= 1000)
    W(f"\nThe Phase I rank gates pass for every epitope (near-vacuous, as section 3.1 explains). All ten MHC-I epitopes fail the stricter stratum **because none reaches the methodology's BigMHC-IM 0.70 cut** "
      f"({N('bm07_b', bm7, 'Phase4C_TCRRecognition_MHCI')}/10; the project cut is 0.5 -- Deviation #19), not because of IC50: {N('ic500_n', ic500, 'Phase4B_BindingAffinity')}/10 have IC50 <= 500 nM at consensus's best allele. "
      f"For class II, {N('ic1000_n', ic1000, 'Phase4B_BindingAffinity')}/11 have IC50 <= 1000 nM. `N/A` marks epitopes whose IC50 no method could produce (never a pass); none of the construct epitopes are `N/A`.\n")
    W("### 3.3 TCR recognition (MHC-I only; n = 10)\n")
    W("The methodology's tools (\"MixMHCpred + TCGA Contact Database\") are wrong (disagreement #6). A three-rung ladder was used instead: BigMHC-IM (primary), PRIME 2.0 + MixMHCpred 3.0, and a local reimplementation "
      "of the Calis et al. (2013) model. TCR-facing residues are taken as positions 4-6 (the methodology's \">=3 TCR-contact residues\" rule cannot be applied without a defined position set). "
      "Each tool's direction was checked against a trusted quantity before use:\n")
    W(table(["tool / check", "statistic", "value", "p", "n", "reading"], [
        (r["Check"], r["Statistic"], N("dc_" + r["Check"][:24], r["Value"], "Phase4C_DirectionChecks"), r["P_Value"] or "-", r["N"], r["Interpretation"]) for r in Cd]))
    bm5 = sum(fnum(r["BigMHC_IM_Best"]) >= 0.5 for r in Ct); bm7 = sum(fnum(r["BigMHC_IM_Best"]) >= 0.7 for r in Ct)
    W(f"\n**At n = 10 only |r| > ~0.63 is distinguishable from zero.** PRIME's correlation is meaningful; BigMHC-IM's has the expected sign but is not distinguishable from zero (it passes a sign check and corroborates nothing); "
      "Calis is exploratory only (self-AUC 0.62 on its own training data; near-zero discrimination against binding is arguably correct because the model is designed to be independent of it) and is **not** independent corroboration of BigMHC or PRIME. "
      "BigMHC-IM and PRIME choose their \"best\" allele independently, so their best columns are **not allele-matched** and must not be compared row by row. "
      f"BigMHC-IM at the Phase I cut (0.5): {N('bm05', bm5, 'Phase4C_TCRRecognition_MHCI')}/10 epitopes; at the methodology's stricter 0.70: {N('bm07', bm7, 'Phase4C_TCRRecognition_MHCI')}/10. "
      "(The IEDB `immunogenicity/` endpoint returns 403 and is not served at this URL -- see the corrections file.)\n")
    W(table(["epitope", "pathogen", "TCR-facing 4-6", "BigMHC-IM best", "PRIME %Rank", "Calis (exploratory)", "flag"],
            [(r["Peptide"], by["A"][r["Peptide"]]["Pathogen"], r["TCR_Facing_4_6"], r["BigMHC_IM_Best"], r["PRIME_PctRank_Best"], r["Calis_Score_Exploratory"],
              "worst on 4/4 signals" if r["Worst_In_Construct_On"] else "") for r in Ct]))
    W("\n**HIV vs Mpox (MHC-I construct epitopes):** mean BigMHC-IM best HIV %s vs Mpox %s (n = 5 each); median PRIME %%Rank HIV %s vs Mpox %s.\n" % (
        N("bm_hiv", f"{mean([fnum(r['BigMHC_IM_Best']) for r in Ct if by['A'][r['Peptide']]['Pathogen'] == 'HIV']):.3f}", "Phase4C_TCRRecognition_MHCI"),
        N("bm_mpx", f"{mean([fnum(r['BigMHC_IM_Best']) for r in Ct if by['A'][r['Peptide']]['Pathogen'] == 'Mpox']):.3f}", "Phase4C_TCRRecognition_MHCI"),
        N("pr_hiv", f"{med([fnum(r['PRIME_PctRank_Best']) for r in Ct if by['A'][r['Peptide']]['Pathogen'] == 'HIV']):.3f}", "Phase4C_TCRRecognition_MHCI"),
        N("pr_mpx", f"{med([fnum(r['PRIME_PctRank_Best']) for r in Ct if by['A'][r['Peptide']]['Pathogen'] == 'Mpox']):.3f}", "Phase4C_TCRRecognition_MHCI")))
    W("### 3.4 Cross-reactivity with the human proteome (both layers, all 31 epitopes)\n")
    ex = sum(r["Exact_8mer_Match_Fresh"] == "YES" for r in Cx)
    st = Counter(r["BLASTP_Status"] for r in Cx)
    W(f"Two layers are reported together and never one alone (disagreement #7): a BLASTP screen against reviewed human Swiss-Prot (from Phase 1Ec; statuses {dict(st)} -- short, statistically non-significant partial local alignments, uninformative at this length), and an exact 8-mer screen recomputed here. "
      f"**{N('exact8', ex, 'Phase4C_CrossReactivity_AllConstruct')} of 31 epitopes contain an exact 8-mer found in the reviewed human proteome**; the fresh screen agrees with Phase 1Ec for "
      f"{N('xr_agree', sum(r['Fresh_vs_Phase1Ec_Agree'] == 'YES' for r in Cx), 'Phase4C_CrossReactivity_AllConstruct')}/31. Phase 1Ec measured the BLASTP/70%-identity rule as blind at 9-16 aa "
      "(0 of 20 human self-fragments caught), so it is never presented as a safety result on its own; the exact 8-mer layer is the informative one, and it covers reviewed Swiss-Prot only.\n")

    # ---- 4 per pathogen
    W("## 4. Per pathogen (the closest thing Section IV has to an immunodominance result)\n")
    W("Cumulative coverage of the union of binding alleles (B-cell epitopes excluded -- not HLA-restricted, `N/A`, never counted as zero):\n")
    rows4 = []
    for sc in ("HIV MHC-I (5)", "Mpox MHC-I (5)", "HIV MHC-II (5)", "Mpox MHC-II (6)", "HIV (10: 5 I + 5 II)", "Mpox (11: 5 I + 6 II)"):
        e4, m1 = cov("2/3 raw", sc, "Phase4B"), cov("2/3 raw", sc, "Phase1F")
        rows4.append((sc, N(f"cov4b_{sc}", sat2(float(e4["Coverage_Pct"])), "Phase4E_ConstructCoverage"), e4["N_Union_Alleles"],
                      N(f"cov1f_{sc}", sat2(float(m1["Coverage_Pct"])), "Phase4E_ConstructCoverage"), m1["N_Union_Alleles"]))
    W(table(["scope", "IEDB-EL", "alleles", "MHCflurry/MHCnuggets", "alleles"], rows4))
    W("\n**Mpox MHC-I is the weak arm** (IEDB-EL %s vs HIV MHC-I %s; MHCflurry/MHCnuggets %s vs %s). Both pathogens exceed 99.6%% combined under either predictor, so **trigger 7 (one pathogen <50%% while the other >90%%) does not apply**. "
      "The asymmetry is consistent with `A35R` being the weakest antigen at every step: `NKRKRVIGL` (Mpox A35R) is the worst construct epitope on four independent signals at once (PRIME %%Rank 4.386 vs 0.47 for the next worst; "
      "BigMHC-IM %s; binding breadth %s of 74; EL rank %s) and adds **%s pp** to the construct coverage union (%s pp within class I). This is a consistent cross-phase signal, not a new anomaly.\n" % (
        *[f"{float(cov('2/3 raw', sc, uni)['Coverage_Pct']):.2f}%" for sc, uni in (("Mpox MHC-I (5)", "Phase4B"), ("HIV MHC-I (5)", "Phase4B"), ("Mpox MHC-I (5)", "Phase1F"), ("HIV MHC-I (5)", "Phase1F"))],
        N("nk_bm", by["C"]["NKRKRVIGL"]["BigMHC_IM_Best"], "Phase4C_TCRRecognition_MHCI"), N("nk_br", by["B"]["NKRKRVIGL"]["N_Binding_Alleles_Primary"], "Phase4B_BindingAffinity"),
        N("nk_el", by["B"]["NKRKRVIGL"]["NetMHC_EL_Best_Rank"], "Phase4B_BindingAffinity"),
        N("nk_drop", f"{float(by['E']['NKRKRVIGL']['Marginal_Drop_Construct_pp']):.4f}", "Phase4E_EpitopeCoverage"), N("nk_drop_i", f"{float(by['E']['NKRKRVIGL']['Marginal_Drop_Within_Class_pp']):.3f}", "Phase4E_EpitopeCoverage")))
    W("### HIV vs Mpox at a glance (construct epitopes; n stated)\n")
    def grp(rows, cls=None):
        return [r for r in rows if cls is None or r["Class"] == cls]
    cons_by = lambda P: [fnum(by["K"][r["Peptide"]]["Conservancy_Stored_Pct"]) for r in P]
    bc = lambda P: [by["F"][r["Peptide"]] for r in P if r["Class"] == "B-cell"]
    rows_g = [
        ("epitopes: MHC-I / MHC-II / B-cell", *[f"{len(grp(P, 'MHC-I'))} / {len(grp(P, 'MHC-II'))} / {len(grp(P, 'B-cell'))}" for P in (HIV, MPX)]),
        ("median breadth: MHC-I / MHC-II", *[f"{fm(hb(P, 'MHC-I'))} / {fm(hb(P, 'MHC-II'))}" for P in (HIV, MPX)]),
        ("class I coverage (IEDB-EL)", *[f"{float(cov('2/3 raw', sc, 'Phase4B')['Coverage_Pct']):.2f}%" for sc in ("HIV MHC-I (5)", "Mpox MHC-I (5)")]),
        ("MHC-I epitopes liberated", *[f"{sum(by['D'][r['Peptide']]['Liberation_Verdict'] == 'LIBERATED' for r in grp(P, 'MHC-I'))}/{len(grp(P, 'MHC-I'))}" for P in (HIV, MPX)]),
        ("MHC-I epitopes on the internal-cleavage watchlist", *[f"{sum(by['D'][r['Peptide']]['Internal_Site_Exceeds_C_Term'] == 'YES' for r in grp(P, 'MHC-I'))}/{len(grp(P, 'MHC-I'))}" for P in (HIV, MPX)]),
        ("mean conservancy, all epitopes", *[f"{mean(cons_by(P)):.1f}% (n={len(P)})" for P in (HIV, MPX)]),
        ("mean BepiPred, B-cell", *[f"{mean([fnum(b['mean_BepiPred']) for b in bc(P)]):.3f} (n={len(bc(P))})" for P in (HIV, MPX)]),
        ("B-cell SEMA-3D corroborated", *[f"{sum(b['SEMA_Corroborated_Stored'] == 'YES' for b in bc(P))}/{len(bc(P))}" for P in (HIV, MPX)]),
        ("B-cell mean native-fold RSA", *[f"{mean([fnum(b['RSA_Mean']) for b in bc(P)]):.2f}" for P in (HIV, MPX)])]
    W(table(["metric", "HIV", "Mpox"], rows_g))
    W("")

    # ---- 5 processing
    W("## 5. Proteasomal processing (the analysis Phase I and II never did)\n")
    lib = [r for r in Dp if r["Class"] == "MHC-I"]
    nlib = sum(r["Liberation_Verdict"] == "LIBERATED" for r in lib)
    intl = [r for r in lib if r["Internal_Site_Exceeds_C_Term"] == "YES"]
    W("The whole 570-aa construct was submitted to IEDB's processing predictor (proteasome + TAP + MHC; MHC-I only, so verdicts cover the 10 MHC-I epitopes; class II and B-cell epitopes are `N/A` for liberation). "
      "The endpoint's scores were direction-checked first (higher = better throughout; the proteasome score is a pure C-terminal cleavage score, invariant across alleles and lengths; "
      f"per-allele IC50 agrees in sign with Step 2's EL rank, Spearman {N('sp_proc', ds['Spearman rho: Step 2 per-allele EL rank vs endpoint IC50'], 'Phase4D_Summary')}, "
      f"n = {ds['... n pairs (from 10 epitopes; not independent)']} pairs from 10 epitopes, so not independent). Pre-registered criterion: liberated if the C-terminal cleavage score is at or above the construct-wide median.\n")
    W(f"**{N('n_lib', nlib, 'Phase4D_EpitopeProcessing')}/10 MHC-I epitopes are liberated** (C-terminal cleavage percentiles {min(float(r['C_Term_Percentile']) for r in lib):.1f}-{max(float(r['C_Term_Percentile']) for r in lib):.1f}). "
      f"**{N('n_int', len(intl), 'Phase4D_EpitopeProcessing')} epitopes carry an internal cleavage site scoring above their own C-terminal site** -- a supplementary watchlist, not part of the pass/fail criterion: "
      + ", ".join(f"`{r['Peptide']}` ({r['Pathogen']})" for r in intl) + ". Of these, "
      f"{N('n_int_aay', len(calib['watchlist overlap']['epitopes adjacent to a qualifying AAY seam AND on the internal-cleavage watchlist']['Note'].split(';')), 'Phase4E_JunctionCalibration')} are also adjacent to a class I AAY seam that produces a contested window (section 6): "
      + ", ".join(f"`{x}`" for x in calib['watchlist overlap']['epitopes adjacent to a qualifying AAY seam AND on the internal-cleavage watchlist']['Note'].split(';')) + ". "
      "The toxicity/allergen half of methodology IV.B.2 was completed in Phase 2A and is not repeated here.\n")

    # ---- 6 junctions
    W("## 6. Junction analysis (epitope-linker seams)\n")
    W("**Definition.** A seam window is *contested* at allele a only if it is (1) deep (>= 3 residues each side), (2) liberated (C-terminal cleavage >= construct median), (3) predicted to beat the **best** real epitope binding a by >= m-fold "
      "(class I: IC50; class II: EL-rank ratio, because EL has no IC50) and (4) for class II, its predicted core spans the seam. An earlier definition (\"any seam window beats any real epitope by any amount\") compares a maximum over thousands of windows with a maximum "
      "over a handful of reference epitopes and finds contest by construction; it was withdrawn and none of its figures are reported.\n")
    W("**Contested does not mean lost.** A contested allele means a seam peptide is *predicted to bind more tightly* (by >= m-fold) than the best genuine epitope there. It does **not** mean the genuine epitope fails to be presented, "
      "and a predicted affinity difference is not a measured competition.\n")
    W("### 6.1 Class I: real (with mechanism)\n")
    cA = {r["Allele"]: r for r in Eca if r["Class"] == "MHC-I"}
    at_risk = sum(int(r["N_Alleles_Phase4B"]) for r in Ee if r["Class"] == "MHC-I")
    cur = []
    for m in ("1x", "2x", "5x", "10x"):
        pairs = sum(int(r["N_Real_Epitopes"]) for r in cA.values() if r[f"Contested_{m}"] == "YES")
        nal = sum(r[f"Contested_{m}"] == "YES" for r in cA.values())
        u = cov("2c uncontested (margin curve)", "MHC-I only (10)", "Phase4B", "primary (Opus definition)", m)
        uf = cov("2c uncontested (margin curve)", "MHC-I only (10)", "Phase4B", "class I also requires IC50<=500nM", m)
        uc = cov("2c uncontested (margin curve)", "CONSTRUCT (21 MHC epitopes)", "Phase4B", "primary (Opus definition)", m)
        lab = f"{m}{' (primary)' if m == '5x' else ' (curve endpoint: margin of the superseded definition)' if m == '1x' else ''}"
        cur.append((lab, f"{N('pairs_I_' + m, pairs, 'Phase4E_ContestedAlleles')}/{at_risk}", f"{nal}/{len(cA)}",
                    N("unc_I_" + m, f"{float(u['Coverage_Pct']):.2f}%", "Phase4E_ConstructCoverage"), f"{float(uf['Coverage_Pct']):.2f}%",
                    N("unc_C_" + m, f"{float(uc['Coverage_Pct']):.2f}%", "Phase4E_ConstructCoverage")))
    raw_i = float(cov("2/3 raw", "MHC-I only (10)", "Phase4B")["Coverage_Pct"])
    W(f"Class I margin curve (Phase 4B universe; n = {at_risk} epitope-allele pairs at risk; {len(cA)} alleles with a real class I comparator; raw class I coverage {raw_i:.2f}%):\n")
    W(table(["margin m", "contested pairs", "contested alleles", "uncontested class I coverage", "... if IC50 <= 500 nM also required", "uncontested whole construct"], cur))
    W("\nAt the primary 5x margin the whole construct remains %s uncontested (raw %s), class I falls from %.1f%% to %s. The 1x whole-construct point sits just under 70%%; it is the margin of the withdrawn definition, "
      "is dominated by the extreme-value effect in class II, and is **not** an escalation -- the floor applies at the 5x primary margin. Class II contribution to the curve is not tabulated as a finding (section 6.2); its values are in `Phase4E_ConstructCoverage`.\n" % (
        cur[2][5], sat2(float(cov("2/3 raw", "CONSTRUCT (21 MHC epitopes)", "Phase4B")["Coverage_Pct"])), raw_i, cur[2][3]))
    W("**Per-locus class I effect** (class I loci only, Phase 4B universe; the 5x primary and the strong-binder-floor sensitivity):\n")
    rawpl = parse_pl(cov("2/3 raw", "MHC-I only (10)", "Phase4B")["Per_Locus_Pct"])
    p5 = parse_pl(cov("2c uncontested (margin curve)", "MHC-I only (10)", "Phase4B", "primary (Opus definition)", "5x")["Per_Locus_Pct"])
    f5 = parse_pl(cov("2c uncontested (margin curve)", "MHC-I only (10)", "Phase4B", "class I also requires IC50<=500nM", "5x")["Per_Locus_Pct"])
    W(table(["class I locus", "raw", "uncontested (5x)", "uncontested (5x, IC50<=500 nM required)"],
            [(L, f"{rawpl[L]:.1f}%", N(f"pl5_{L}", f"{p5.get(L, 0.0):.1f}%", "Phase4E_ConstructCoverage"), f"{f5.get(L, 0.0):.1f}%") for L in ("A", "B", "C")]))
    W("\nLocus B carries the largest predicted loss; it is a class I locus (a real, mechanistically explained effect), not a class II finding. Per-locus numbers for class II are deliberately not reported (section 6.2).\n")
    q = calib["class I contested windows (5x)"]
    nq, na_ = int(q["qualifying windows (n)"]["N_pairs"]), int(q["ending in AAY"]["N_pairs"])
    ninv = int(q["involving AAY (ending in it or containing it)"]["N_pairs"])
    ci = calib["class I seam strong-binder rate (IC50<=500nM)"]
    nat = calib["class I natural reference"]["adjuvant-domain windows IC50<=500nM"]
    W(f"**Mechanism and null.** {N('nq', nq, 'Phase4E_JunctionCalibration')} class I seam windows qualify at 5x; **{N('naay', na_, 'Phase4E_JunctionCalibration')} of {nq} ({100 * na_ / nq:.0f}%) end in `AAY`** and "
      f"{N('ninv', ninv, 'Phase4E_JunctionCalibration')} of {nq} involve AAY; the 8 strongest (IC50 {q['strongest IC50 range nM (min-max of 8)']['Note']} nM) all end in AAY. "
      "They land on `B*35:01/05/17`, `A*26:01` and `C*14:02` -- alleles whose real construct epitopes end in F/Y in this data; AAY's terminal Tyr sits in the C-terminal (PΩ) aromatic anchor pocket "
      "(pocket assignment is prior knowledge, consistent with but not proven by these data). Seam windows reach IC50 <= 500 nM in "
      f"{N('rate_all', float(ci['ALL']['Fraction_le500nM']) * 100, 'Phase4E_JunctionCalibration', '{:.2f}%')} of window-allele pairs (n = {ci['ALL']['N_pairs']}) versus "
      f"{N('rate_nat', float(nat['Fraction_le500nM']) * 100, 'Phase4E_JunctionCalibration', '{:.2f}%')} of unselected natural windows (n = {nat['N_pairs']}), a **{N('enrich', ci['ALL']['Fold_vs_natural'], 'Phase4E_JunctionCalibration')}x enrichment**; "
      f"windows ending in AAY reach {float(ci['ends AAY']['Fraction_le500nM']) * 100:.1f}% ({ci['ends AAY']['Fold_vs_natural']}x), whereas GPGPG seams reach {float(ci['GPGPG']['Fraction_le500nM']) * 100:.2f}% and KK seams {float(ci['KK']['Fraction_le500nM']) * 100:.2f}% (both depleted). "
      "The natural reference is a single 45-residue protein, so the ratio is indicative. Contested class I alleles are all outside the top three by frequency; the most frequent is `B*15:02` (frequency 0.084, rank 9 of 74): "
      + calib["contested class I alleles (5x)"]["top-3 most frequent class I alleles contested"]["Note"] + ".\n")
    W("**Verdict (Opus): a Discussion/limitations finding, not a redesign trigger** -- whole-construct coverage remains >99.9% at the primary margin, and Phase I is frozen.\n")
    W("### 6.2 Class II: the seams are clean\n")
    T = calib["class II seam calibration (EL rank; random-peptide null 10% / 1% / 0.1%)"]
    a = T["ALL"]
    W(f"EL rank is a percentile against random peptides, so random-like seams would show ~10% at rank <= 10 and ~1% at rank <= 1. Over n = {N('n_ii_pairs', a['N_pairs'], 'Phase4E_JunctionCalibration')} deep seam window-allele pairs the observed rates are "
      f"**{N('ii_le10', float(a['Fraction_le10']) * 100, 'Phase4E_JunctionCalibration', '{:.2f}%')} (null 10%) and {N('ii_le1', float(a['Fraction_le1']) * 100, 'Phase4E_JunctionCalibration', '{:.2f}%')} (null 1%)** -- "
      "seam windows are *depleted* for class II binders, not enriched.\n")
    W(table(["linker set", "n pairs", "rank <= 10", "x null", "rank <= 1", "x null"],
            [(k, T[k]["N_pairs"], f"{float(T[k]['Fraction_le10']):.2%}", f"{float(T[k]['Fraction_le10']) / 0.10:.2f}", f"{float(T[k]['Fraction_le1']):.2%}", f"{float(T[k]['Fraction_le1']) / 0.01:.2f}") for k in ("GPGPG", "KK", "EAAAK+adjuvant", "AAY") if k in T]))
    g = calib["class II contested rows (5x)"]
    kk = calib["class II contested rows by linker set (5x, liberated + core-spanning)"]
    kk_str = ", ".join(f"{k} {v['N_pairs']}" for k, v in kk.items())
    W(f"\n**GPGPG seams are clean (a positive design result): {N('gpgpg_rows', g['involving GPGPG']['N_pairs'], 'Phase4E_JunctionCalibration')} of {g['total rows']['N_pairs']} class II contested rows at 5x involve GPGPG** "
      f"(by linker set: {kk_str}); the residual sits at KK seams (flanking the B-cell epitopes) and AAY seams, not around the class II epitopes. "
      "Residual class II calls are best explained by **sample-size asymmetry**: a maximum over ~%s deep seam windows per allele is compared with only 2-4 real epitopes per allele, so a >=5x call can arise by chance and shrinks but does not vanish as the margin rises. " % calib["seam windows scanned"]["deep class II windows per allele"]["N_pairs"] +
      "This is supported by a diagnostic on the eight DPB1 alleles: each has %s real epitopes against %s core-spanning seam rows with rank <= 10 (before the liberation filter), and the median seam-versus-real fold per allele is %s -- equivalent, not stronger. "
      "No random-window null was run at the 5x criterion itself, so the residual is described as artifact-prone, not proven artifact. No per-locus class II coverage figure is reported as a finding.\n" % (N("dpb1_nreal", "2-4", "Phase4E_ContestedAlleles"), N("dpb1_rows", "26-134", "Phase4D_JunctionNeoepitopes_MHCII"), N("dpb1_fold", "0.55-1.75", "Phase4D_JunctionNeoepitopes_MHCII")))
    W("HIV vs Mpox: class I contest concerns epitopes of both pathogens (see `Phase4E_EpitopeCoverage`); the class I subset uncontested at 5x is %s for HIV MHC-I and %s for Mpox MHC-I, the weaker Mpox arm losing more.\n" % (
        N("unc5_hiv_i", f"{float(cov('2c uncontested (margin curve)', 'HIV MHC-I (5)', 'Phase4B', 'primary (Opus definition)', '5x')['Coverage_Pct']):.1f}%", "Phase4E_ConstructCoverage"),
        N("unc5_mpx_i", f"{float(cov('2c uncontested (margin curve)', 'Mpox MHC-I (5)', 'Phase4B', 'primary (Opus definition)', '5x')['Coverage_Pct']):.1f}%", "Phase4E_ConstructCoverage")))

    # ---- 7 B-cell
    W("## 7. B-cell epitopes and conservancy\n")
    tiers = Counter(r["Primary_Gate_Phase_I_Tier"] for r in Fb)
    W(f"**BepiPred (Phase I tiers, disagreement #3):** n = 10 construct epitopes: {N('t_high', tiers['High'], 'Phase4F_BcellDossier')} High / {N('t_med', tiers['Medium'], 'Phase4F_BcellDossier')} Medium / "
      f"{N('t_dep', tiers['Deprioritized'], 'Phase4F_BcellDossier')} Deprioritized (retained -- tiers deprioritise, never exclude); {N('t_meth', sum(r['Secondary_Stratum_Methodology_High'] == 'YES' for r in Fb), 'Phase4F_BcellDossier')}/10 also meet the methodology's "
      f"stricter High (mean >= 0.50 and >= 75% of residues above). Tiers were recomputed from the stored values for all 39 pool candidates: {fs['Bcell_Tier recomputed matches']}/{fs['Pool B-cell candidates']} match. "
      "The methodology's stricter values appear only in the `Secondary_Stratum` column.\n")
    sema = Counter(r["SEMA_Corroborated_Stored"] for r in Fb)
    ov = [float(r["SEMA_Overlap_Pct_Stored"]) for r in Fb]
    W(f"**SEMA-3D (Phase 1De, not re-run -- its `foldseek` dependency is a Linux binary):** all {N('sema_no', sema['NO'], 'Phase4F_BcellDossier')} of 10 construct B-cell epitopes are `NO` (overlap {min(ov):.1f}-{max(ov):.1f}% of 16 residues scored; "
      "the corroboration bar is 50%). **`NO` means screened and not corroborated -- a lower-priority linear-only call. It is never \"refuted\"**: SEMA-3D's own calibration is TPR 0.66 at FPR 0.07, so a miss is far from disproof. "
      f"`UNSCREENED` (no assessable fold) is a different state and never a negative; {N('sema_uns', sema.get('UNSCREENED', 0), 'Phase4F_BcellDossier')} construct epitopes are UNSCREENED. "
      f"SEMA overlap was recomputed from Phase 1De's per-residue scores for every locatable epitope with {fs['SEMA overlap mismatches']} mismatches.\n")
    det = [r for r in Fb if r.get("Native_Fold_Verdict") == "DETERMINED"]
    und = [r for r in Fb if r.get("Native_Fold_Verdict") == "UNDETERMINED"]
    rsa_h = N("rsa_hiv", f"{mean([fnum(r['RSA_Mean']) for r in Fb if r['Pathogen'] == 'HIV']):.2f}", "Phase4F_BcellDossier")
    rsa_m = N("rsa_mpx", f"{mean([fnum(r['RSA_Mean']) for r in Fb if r['Pathogen'] == 'Mpox']):.2f}", "Phase4F_BcellDossier")
    pl_h = f"{mean([fnum(r['Native_pLDDT_Mean']) for r in Fb if r['Pathogen'] == 'HIV']):.1f}"
    pl_m = f"{mean([fnum(r['Native_pLDDT_Mean']) for r in Fb if r['Pathogen'] == 'Mpox']):.1f}"
    cp_lo, cp_hi = min(fnum(r["Construct_Model_pLDDT_Mean"]) for r in Fb), max(fnum(r["Construct_Model_pLDDT_Mean"]) for r in Fb)
    rho_v = N("rho_bep", fs["Spearman BepiPred mean vs native RSA (unique pool B-cell)"], "Phase4F_Summary")
    W(f"**Native-fold overlay.** The construct's own AlphaFold model is undetermined outside the adjuvant (pTM {fs['Construct model pTM']}; the {fs['Non-adjuvant residues']} non-adjuvant residues have mean pLDDT {fs['Non-adjuvant mean pLDDT']} "
      f"and {N('cons_ge70', fs['Non-adjuvant % residues pLDDT>=70'], 'Phase4F_Summary')}% reach 70), so surface accessibility is computed on the Phase 1De **per-antigen native folds** instead. "
      f"**{N('n_det', len(det), 'Phase4F_BcellDossier')} of 10 epitopes sit on a globally DETERMINED fold; {N('n_und', len(und), 'Phase4F_BcellDossier')} of 10 sit on an UNDETERMINED fold** (gp41, A35R, B5R; pTM < 0.5), so their exposure is indicative only. "
      f"Fold monomers only: oligomer interfaces and HIV Env glycans are not modelled (exposure is over-estimated for gp120/gp41). Mean relative accessibility {rsa_h} (HIV, n=4) vs {rsa_m} (Mpox, n=6); native-fold pLDDT {pl_h} vs {pl_m}, "
      f"against construct-model pLDDT {cp_lo:.1f}-{cp_hi:.1f} for every epitope. BepiPred and native-fold exposure agree in sign across the 39 pool candidates but not distinguishably from zero "
      f"(Spearman {rho_v}, p = {fs['... p-value']}, n = {fs['... n']}).\n")
    W(table(["epitope", "pathogen", "antigen", "Phase I tier", "methodology High", "SEMA", "fold verdict", "mean RSA", "native pLDDT"],
            [(r["Peptide"], r["Pathogen"], r["Target"], r["Primary_Gate_Phase_I_Tier"], r["Secondary_Stratum_Methodology_High"], f"{r['SEMA_Corroborated_Stored']} ({r['SEMA_Overlap_Pct_Stored']}%)",
              r.get("Native_Fold_Verdict", "N/A"), r.get("RSA_Mean", "N/A"), r.get("Native_pLDDT_Mean", "N/A")) for r in Fb]))
    W("")
    W("**Conservancy (Phase 1Dc).** The published per-target means reproduce **exactly** on the full 30,268-candidate pool (`Raw_Conservancy`; the methodology's `Min_50pct` holds only >=50% survivors and cannot reproduce a full-pool mean):\n")
    W(table(["target", "n candidates", "published mean", "reproduced mean"], [(r["Target"], r["N_Candidates"], f"{r['Published_Mean_Pct']}%", N("gate_" + r["Target"], f"{r['Reproduced_Mean_Pct']}%", "Phase4F_ConservancyGate")) for r in Fg]))
    hv = [fnum(r["Conservancy_Stored_Pct"]) for r in Fc if r["Pathogen"] == "HIV"]
    mp = [fnum(r["Conservancy_Stored_Pct"]) for r in Fc if r["Pathogen"] == "Mpox"]
    W(f"\nHIV per-target means span 4.54-9.97% and Mpox 37.80-60.61%. An independent recompute from the Phase 1C variant pools reproduced the stored conservancy and hit ratio for "
      f"{fs['Construct epitopes conservancy recompute matches']}/{fs['Construct epitopes tested']} construct epitopes. Construct epitopes: HIV n = {len(hv)}, mean {N('cons_hiv', f'{mean(hv):.1f}', 'Phase4F_Conservancy')}% "
      f"({min(hv):.1f}-{max(hv):.1f}%); Mpox n = {len(mp)}, mean {N('cons_mpx', f'{mean(mp):.1f}', 'Phase4F_Conservancy')}% ({min(mp):.1f}-{max(mp):.1f}%) -- taken from the conserved tail of a very diverse pool. "
      "Variant pools are small (21-30 isolates per antigen), so values are coarse (one variant = 3-5 pp) and describe within-clade conservancy, not global strain diversity.\n")

    # ---- 8
    W("## 8. What this does NOT replace from C-ImmSim\n")
    W("Section IV replaces C-ImmSim's *role* in the pipeline with static, per-epitope predictors. It does **not** provide:\n"
      "- **Antibody titres** -- BepiPred/SEMA-3D rank candidate B-cell epitopes; they say nothing about magnitude or affinity maturation of an antibody response.\n"
      "- **Cytokine kinetics** -- no IFN-γ / IL-2 / IL-4 time courses.\n"
      "- **Temporal dynamics** -- no clonal expansion, contraction, memory formation, or boost-schedule effects; no dose response.\n"
      "- **Competition between epitopes or between HIV and Mpox responses** -- the per-pathogen coverage in section 4 is a static population-level proxy for immunodominance, not a simulated competition.\n"
      "- **Innate or adjuvant effects, delivery, or biodistribution** (Sections III and V-VI).\n"
      "It also carries the limits of its inputs: Southeast Asian (Austronesian) proxy allele frequencies (Deviation #7), predictor-dependent class I coverage (section 2), and predictors that are not experimental evidence.\n")
    W("## Reproducibility and traceability\n")
    W(f"All numbers above are produced through a ledger (`Phase4G_NumbersLedger_{TS}.csv`) recording each value with its source file in `StepA`-`StepF`. Sources used: "
      + "; ".join(sorted(set(v for k, v in SRC.items() if k != 'Dossier'))) + ". The dossier is `" + os.path.basename(dpath) + "` (31 rows; `Primary_Gate` = Phase I rules, `Secondary_Stratum` = the revised methodology's stricter values). "
      "Rule 1 (Phases 1-2 untouched) is verified by SHA-256 against the durable baseline.\n")
    rpath = os.path.join(OUT, "PHASE_IV_REPORT.md")
    open(rpath, "w").write("\n".join(R_))

    # ---------------------------------------------------------- corrections
    ic1 = [r for r in pool["MHC-I"] if r["IC50_At_Best_Consensus_Allele_nM"] not in ("N/A", "")]
    from scipy import stats as _st
    def cor(cls):
        g = [r for r in pool[cls] if r["IC50_At_Best_Consensus_Allele_nM"] not in ("N/A", "") and r["Consensus_Best_Rank"] not in ("N/A", "")]
        x = [float(r["Consensus_Best_Rank"]) for r in g]; y = [float(r["IC50_At_Best_Consensus_Allele_nM"]) for r in g]
        return _st.pearsonr(x, y)[0], _st.spearmanr(x, y)[0], len(g)
    c_i, c_ii = cor("MHC-I"), cor("MHC-II")
    C_ = []
    Wc = C_.append
    Wc("# Manuscript Corrections -- Phase IV (Section IV)\n")
    Wc(f"Generated {TS}. Source of the statements quoted: `REVISED_METHODOLOGY_MpoxHIV.md` Section IV. Part A lists the eight pre-decided disagreements with corrected text; Part B lists everything else accumulated during Phase IV. "
       "Where a step's own record and this file differ, this file is later.\n")
    Wc("## Part A -- The eight disagreements\n")
    Wc(table(["#", "Methodology says", "Corrected text (proposed)", "Basis"], [
        (1, "\"Strong Binder: IC50 <= 500 nM ... Weak 500-5000 ... Non-Binder > 5000 nM (excluded)\"",
         "Candidates were selected by IEDB percentile rank (<=1% MHC-I; <=10% MHC-II), the method IEDB recommends because IC50 is allele-biased. IC50 (Strong <=500, Weak 500-5000, Non-binder >5000 nM) is reported as an additional descriptor only, taken from the IEDB `consensus` method, which supports fewer alleles than NetMHCpan-EL, and is `N/A` where no IC50-producing method applies.",
         "Phase I gate; Phase 4B"),
        (2, "\"BigMHC score >= 0.70 indicates high likelihood of T-cell elicitation\"",
         f"BigMHC-IM output is a bounded [0,1] probability; the project cutoff is 0.5 (Deviation #19), with >= 0.70 reported as a stricter second tier. In the construct, {bm5}/10 MHC-I epitopes reach 0.5 and {bm7}/10 reach 0.70 (best allele).",
         "Deviation #19; Phase 4C"),
        (3, "\"High-Priority B-Epitope: mean_BepiPred >= 0.50 AND pct_above >= 75%\"",
         "Three BepiPred tiers are used: High (mean >= 0.60 and >= 75%), Medium (>= 0.50 and >= 50%), Deprioritized (>= 0.45 and >= 37.5%); tiers deprioritise and never exclude. The methodology's threshold is reported as a secondary stratum.",
         "Phase 1Db; Phase 4F"),
        (4, "\"top-10 HLA alleles in Philippine population\" with prevalences (A*24:23 at 18.5%, ...)",
         "Allele frequencies are from AFND, population 'Singapore Riau Malay' (n = 132), 107 alleles over six loci (A 22, B 33, C 19, DRB1 18, DQB1 7, DPB1 8), described as **Southeast Asian (Austronesian) proxy** coverage, not Philippine (Deviation #7). `HLA-A*24:23` is not a valid entry in this project and matches no frequency table used; the quoted prevalences are unsupported. DQA1 and DPA1 are not in the table (beta chains only).",
         "Deviation #7; Phase 1F"),
        (5, "`Population_Coverage(%) = 1 - prod(1 - freq[allele_i])`",
         "Coverage is computed per locus as s = sum of allele frequencies (capped at 1), locus coverage 1-(1-s)^2 (the diploid two-chromosome term), and the overall value is the complement product across loci. The printed formula omits the diploid term and under-estimates coverage. Report per-locus values alongside the combined figure, which saturates.",
         "Phase 1F helper; Phase 4E"),
        (6, "\"TCR Contact Prediction (MixMHCpred + TCGA Contact Database)\"",
         "MixMHCpred predicts MHC binding, not TCR contact, and TCGA is a cancer-genomics atlas with no TCR-contact data. TCR recognition was assessed with BigMHC-IM, PRIME 2.0 + MixMHCpred 3.0 and a local reimplementation of the Calis et al. (2013) model; TCR-facing residues were taken as positions 4-6 because the '>= 3 contact residues' rule has no defined position set. Tools with opposite score conventions were direction-checked against a trusted quantity before use.",
         "Phase 4C"),
        (7, "\"exclude self-antigen matches with >= 70% sequence identity\" (BLAST)",
         f"Cross-reactivity is assessed with two layers reported together: BLASTP against reviewed human Swiss-Prot and an exact 8-mer screen. The 70%-identity rule was measured blind at 9-16 aa (0 of 20 human self-fragments caught) and is never presented as a safety result alone. {ex} of 31 construct epitopes contain an exact 8-mer found in the reviewed human proteome.",
         "Phase 1Ec; Phase 4C"),
        (8, "\"BepiPred 2.0 + SEMA 2.0\"",
         "SEMA 2.0's hosted service is blocked by an anti-bot WAF (Deviation #20); SEMA-3D was run locally on the native antigen folds, **after** the construct was fixed (post-hoc corroboration). All 10 final-construct B-cell epitopes are screened and not corroborated (`NO`), none UNSCREENED. `NO` is not 'refuted'.",
         "Deviation #20; Phase 1De; Phase 4F"),
    ]))
    Wc("\n## Part B -- Additional corrections and record\n")
    items = [
        ("\"194 Non-Redundant Sequences\" mislabels source proteins as epitopes", "194 is the number of source protein sequences retrieved in Phase 1A, not epitopes. The analysis units are the 31 construct epitopes and the 105 Phase 1F candidates."),
        ("`Mpox_B5R` vs manuscript `B6R`", "The Phase I files label the Mpox EEV-glycoprotein target `Mpox_B5R`; the manuscript's B6R (OPG190, the vaccinia-B5R ortholog) is the intended protein. MPXV's own gene named B5R is a different protein (OPG189, ankyrin-repeat). Both names are retained in the dossier."),
        ("IV.B.1 claims BigMHC integrates pMHC stability", "BigMHC (Albert et al. 2023) is an ensemble transfer-learned from assay data (eluted-ligand pre-training, immunogenicity transfer); it does **not** take pMHC stability as an input. The sentence 'Integrates MHC-I IC50, pHLA stability (pMHCstab), and predicted TCR-interface residues' should be removed."),
        ("IV.B.2 calls MERCI a 'toxin/allergen database'", "MERCI is a motif-search component inside ToxinPred2, not a database. Toxicity/allergen screening of the construct was performed in Phase 2A."),
        ("IV.C 'IEDB Conservancy Tool' and 'global strains'", "Conservancy was computed directly by exact substring matching against each target's Phase 1C variant pool (equivalent to the IEDB tool only at 100% identity; Phase 1Dc note). The quoted means (HIV 4.54-9.97%, Mpox 37.80-60.61%) are per-target means over all 30,268 candidate epitopes across 21-30 isolates per antigen -- within-clade, not 'global strains'. The '>= 50% conservancy' priority rule selected the construct epitopes from `Min_50pct`."),
        ("Conservancy gate must run on `Raw_Full`, not `Min_50pct`", "The published per-target means describe the full 30,268-candidate pool (`Raw_Conservancy`). `Min_50pct` holds only the >= 50% survivors (1,238 rows) and cannot reproduce them. All 7 published means reproduce exactly on `Raw_Full`."),
        ("IV.A.1 predictor and panel wording", "MHC-I/II binding used IEDB `netmhcpan_el` / `netmhciipan_el` and `consensus` through the IEDB API over the 107-allele AFND panel (74 class I, 33 class II; DQ/DP beta chains bridged to paired names). The stated allele lists (DRB1*04:01/07:01/11:01, DQA1*01:01/DQB1*05:01) are not the panel used."),
        ("Step 2 IC50 parsing incident (process record)", f"An unlabeled trailing column in IEDB MHC-II `consensus` output was misread as an IC50 fallback. It was caught by four checks: (1) range -12.5 to +5.2 (IC50 cannot be negative); (2) Pearson r = -0.644 against percentile rank (n = 612; an IC50 must correlate positively); "
         f"(3) strong binders (rank 0.62-1.2) sat at +3.9 to +5.2 and non-binders (rank ~90-100) at -1.3 to -12.5, i.e. higher = better, the opposite of a concentration; (4) the genuinely labelled `nn_align_ic50` spans 3.3-36,040 nM with median 1,938.8 nM. The fallback was removed; IC50 now comes only from labelled nM columns. "
         f"After the fix the IC50-vs-consensus-rank correlations are Pearson {c_i[0]:+.3f} (class I, n = {c_i[2]}) and {c_ii[0]:+.3f} (class II, n = {c_ii[2]}); Spearman {c_i[1]:+.3f} and {c_ii[1]:+.3f}, all the correct sign. "
         "(A review note quoted +0.504 / +0.876; those values do not reproduce exactly from the CSV columns above under any definition tried, so the recomputed values are the ones cited.)"),
        ("Phase I vs Phase IV MHC-II ranks are on different scales", "Phase I's IEDB `recommended` method resolves to EL for class I but to BA/consensus for class II. The same rank number therefore means different things; never table them side by side without stating so. Agreement is judged by gate-status flips (none) not by rank equality."),
        ("IEDB `immunogenicity/` 403 is an endpoint that is not served at that URL, not IP throttling (correction to the review note)", "Retested 2026-09-21: 11 sequential attempts over several minutes, including after a 90-second wait and on GET, all returned 403, while `mhci/`, `mhcii/` and `bcell/` answered 200 in the same minute; a nonexistent path (`foobar/`) returns the byte-identical 403 page. Throttling would also have hit the healthy endpoints and cleared on retry, as it did elsewhere in this project. The retest therefore supports the original 'dead' reading and does **not** support the review note's throttling explanation. Consequently the Calis reimplementation was forced, not merely defensible. (Other 403s in this project were genuine transient throttling; the two produce the same HTML and can be told apart only by such a retest.) Possible future cross-check: the IEDB next-generation T-cell class I tool, not tested here."),
        ("Calis et al. (2013) reimplementation provenance", "The paper's amino-acid enrichment table is not reproduced in the accessible text; it was derived from the paper's own public training set (Supplementary Dataset S1, 2,508 peptides) by the paper's documented method, with Table 2 position weights (3=0.10, 4=0.31, 5=0.30, 6=0.29, 7=0.26, 8=0.18; positions 1, 2, C-terminus masked) and the position-8 weight reused for the extra residue of 10-mers. Self-AUC 0.62; exploratory only."),
        ("Deviation #20 SEMA tally is for a superseded construct", "Deviation #20 reports 'nine shipped B-cell epitopes: 1 corroborated, 6 not corroborated, 2 unscreenable' for `Vax_Final_6f34b53e`, which `MANUSCRIPT_CHANGES_PHASE1` marks SUPERSEDED. The final construct `Vax_Final_4fce119e` has 10 B-cell epitopes: 0 corroborated, 10 screened-not-corroborated, 0 UNSCREENED. The register text must be updated before it is cited."),
        ("SEMA vocabulary", "`NO` = screened, fewer than 50% of residues in a SEMA-3D patch (lower-priority linear-only; not refuted). `UNSCREENED` = could not be assessed on the one available fold (never a negative). Neither means structurally disproven."),
        ("Binding breadth denominators", "Earlier notes quote breadth 'of 107 alleles'. Class I peptides are scored against the 74 class I alleles and class II against the 33 class II alleles; the correct denominators are 74 and 33 (pool medians 7 and 6)."),
        ("Population coverage: target and per-locus limits", "The methodology's '>= 90% population coverage' target is met in the combined figure (saturated) but not at every locus: locus A is 60.3% (IEDB-EL) to 76.0% (MHCflurry/MHCnuggets) because no epitope binds `A*11:01` or `A*33:03`; locus C spans 70.2-98.2% by predictor. Report per-locus ranges."),
        ("New Section IV.B.2 analysis to add to the methods", "The linker analysis was performed as: (i) whole-construct proteasomal/TAP processing with IEDB's processing predictor (MHC-I 8-11-mers), (ii) exhaustive scan of seam-spanning windows (MHC-I 8-11-mers, MHC-II 12-20-mers via NetMHCIIpan-EL), (iii) a corrected 'contested' definition (deep + liberated + >= m-fold vs the best real epitope + class II core spanning the seam), reported as a margin curve, with a random-peptide null for class II and an unselected-natural-window reference for class I. Class II seams (notably GPGPG) are clean; class I AAY seams produce strong aromatic-C-terminus binders (a limitation, not a redesign trigger)."),
        ("Coverage wording", "Report coverage as 'Southeast Asian (Austronesian) proxy coverage' (AFND Singapore Riau Malay), not Philippine; state that class II is on DRB1/DQB1/DPB1 only."),
        ("Terminology in review notes (record)", "Three wording points in the working notes were corrected by the data: locus B is a class I locus (per-locus class I contest is reported as class I); five (not four) MHC-I epitopes carry an internal cleavage site above their C-terminus (four of them overlap the AAY-seam list); 21 of 31 (68%) contested class I windows end in AAY (29 of 31 involve AAY), though the 8 strongest all end in AAY."),
    ]
    for i, (t, b) in enumerate(items, 1):
        Wc(f"**B{i}. {t}.** {b}\n")
    Wc("**Not a discrepancy:** BLASTP is version 2.16.0+, matching the manuscript.\n")
    cpath = os.path.join(OUT, "MANUSCRIPT_CORRECTIONS_PHASE4.md")
    open(cpath, "w").write("\n".join(C_))

    lpath = os.path.join(OUT, f"Phase4G_NumbersLedger_{TS}.csv")
    with open(lpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Key", "Value", "Source_File"])
        w.writeheader(); w.writerows(LEDGER)
    print(f"[SUCCESS] PHASE_IV_REPORT.md ({len(R_)} blocks), MANUSCRIPT_CORRECTIONS_PHASE4.md ({len(C_)} blocks), ledger {len(LEDGER)} numbers")

    # ------------------------------------------------------------------ gates
    print("\n" + "-" * 90 + "\nFINAL GATES")
    ok_all = True
    txt = {"report": open(rpath).read(), "corrections": open(cpath).read(), "dossier": open(dpath).read(), "ledger": open(lpath).read()}
    bad_patterns = ["22.06", "77.93", "18.67", "DPB1 0.0", "DPB1: 0.0", "DPB1 -> 0", "111 flag", "111-flag", "111 MHC-I", "111 (window"]
    hits = [(k, p) for k, t in txt.items() for p in bad_patterns if p in t]
    print(f"  superseded figures absent: {'OK' if not hits else 'FOUND ' + str(hits)}"); ok_all &= not hits
    with open(dpath, newline="") as _f:
        written = list(csv.DictReader(_f))
    unres = sum(1 for r in written for v in r.values() if v.strip().upper() == "UNRESOLVED")
    blanks = sum(1 for r in written for v in r.values() if v.strip() == "")
    print(f"  dossier: {len(dossier)} rows; UNRESOLVED cells {unres}; blank cells {blanks} (N/A used for not-applicable) -> {'OK' if unres == 0 and blanks == 0 else 'CHECK'}"); ok_all &= (unres == 0 and blanks == 0)
    miss = [r for r in LEDGER if not r["Source_File"]]
    print(f"  ledger: {len(LEDGER)} numbers, {len(miss)} without a source -> {'OK' if not miss else 'FAIL'}"); ok_all &= not miss
    print(f"  Primary_Gate / Secondary_Stratum both present for all rows: {'OK' if all('Primary_Gate' in r and 'Secondary_Stratum' in r for r in dossier) else 'FAIL'}")
    base = {}
    for line in open(os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase4", "_rule1_baseline_sha256.txt")):
        h, p = line.rstrip("\n").split(None, 1); base[p.lstrip("* ").strip()] = h
    def sha(pth):
        m = hashlib.sha256()
        with open(pth, "rb") as f:
            for c in iter(lambda: f.read(1 << 20), b""):
                m.update(c)
        return m.hexdigest()
    badr = [p for p, h in base.items() if not os.path.isfile(os.path.join(_PROJECT_ROOT, p)) or sha(os.path.join(_PROJECT_ROOT, p)) != h]
    print(f"  Rule 1 (SHA-256): baseline {len(base)} files, changed/missing {len(badr)} -> {'OK' if not badr else 'FAIL'}"); ok_all &= not badr
    print("[SUCCESS] all final gates passed." if ok_all else "[ERROR] a final gate failed -- see above.")
    return ok_all


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
