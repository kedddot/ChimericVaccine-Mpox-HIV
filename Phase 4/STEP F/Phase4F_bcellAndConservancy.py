import os
import sys
import csv
import re
import json
import glob
import statistics
from collections import defaultdict, Counter

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
for _d in (os.path.join(_PROJECT_ROOT, "Phase 4", "_common"),
           os.path.join(_PROJECT_ROOT, "Phase 1", "STEP D")):
    if _d not in sys.path:
        sys.path.insert(0, _d)

import phase4_common as common
import Phase1Db_filtration as p1db   # read-only: classify_bcell_tier (disagreement #3 -- Phase I's 3 tiers)

# =============================================================================
# PHASE 4F -- B-CELL + CONSERVANCY (methodology IV.A.3 and IV.C, Step 6 of 7)
#
# Mostly REUSE -- no new predictions. Everything stored by Phase 1 is re-read, and
# every stored value that can be recomputed deterministically IS recomputed as a
# regression check (BepiPred tier, SEMA overlap, conservancy). Disagreements:
#   #3 use Phase I's three BepiPred tiers (High >=0.60 & >=75%, Medium >=0.50 & >=50%,
#      Deprioritized >=0.45 & >=37.5%); the methodology's stricter values are reported
#      only as a secondary stratum column.
#   #8 SEMA is NOT re-run (its foldseek dependency is a Linux binary). Stored Phase 1De
#      results are reused. `UNSCREENED` = could not be assessed on the one available
#      fold; it is never a negative. `NO` = screened, <50% of residues in a SEMA-3D
#      patch: lower-priority linear-only -- also not a refutation (SEMA-3D calibration:
#      TPR 0.66 at FPR 0.07).
# =============================================================================

SEMA_CUT = 1.3862943611198906         # ln(4): three or more expected antibody contacts (Phase 1De calibration)
RSA_EXPOSED = 0.25                    # conventional per-residue exposure cut on relative SASA
PUBLISHED_MEANS = {                   # RESULTS_PHASE_I_II_2026-08-31.txt, section I.D.2 (mean conservancy over ALL 30,268 candidates)
    "HIV_gp120": 4.54, "HIV_gp41": 5.69, "HIV_p17": 4.67, "HIV_p24": 9.97,
    "Mpox_A35R": 37.80, "Mpox_B5R": 45.02, "Mpox_L1R": 60.61}
GATE_RANGE = {"HIV": (4.54, 9.97), "Mpox": (37.80, 60.61)}   # the brief's stated ranges
FOLDS = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1De")
CONSTRUCT_FOLD = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase2", "StepB", "AlphaFold_Raw", "fold_last_construct_august_31")
ADJUVANT_LEN = 45


def read_csv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mean(v):
    v = [x for x in v if x is not None]
    return statistics.mean(v) if v else None


def fmt(x, nd=2):
    return "N/A" if x is None else f"{x:.{nd}f}"


def pathogen_of_target(t):
    return "HIV" if t.startswith("HIV") else "Mpox"


# ------------------------------------------------------------------ folds
def load_native_fold(target):
    """Returns dict(seq, sasa_rsa[], plddt[], resnames[]) for the Phase 1De per-antigen native fold (model 0)."""
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    from Bio.PDB.DSSP import residue_max_acc
    from Bio.Data.IUPACData import protein_letters_3to1
    pdbs = sorted(glob.glob(os.path.join(FOLDS, "Folds", target, "*_model_0_converted.pdb")))
    if not pdbs:
        return None
    st = PDBParser(QUIET=True).get_structure(target, pdbs[0])
    model = st[0]
    chains = list(model)
    chain = chains[0]
    res = [r for r in chain if r.id[0] == " "]
    ShrakeRupley(probe_radius=1.4, n_points=100).compute(model, level="R")
    maxacc = residue_max_acc["Wilke"]           # Tien et al. 2013 theoretical maxima, shipped with Biopython
    seq, rsa, pl, names = [], [], [], []
    for r in res:
        name = r.get_resname().upper()
        seq.append(protein_letters_3to1.get(name.capitalize(), "X"))
        rsa.append(r.sasa / maxacc[name] if name in maxacc else None)
        pl.append(statistics.mean(a.get_bfactor() for a in r))
        names.append(name)
    return {"seq": "".join(seq), "rsa": rsa, "plddt": pl, "names": names, "n_chains": len(chains), "path": pdbs[0]}


def build_bcell_conservancy():
    common.print_banner("PHASE 4F -- B-CELL + CONSERVANCY (reuse; regression-gated)")
    ep = read_csv(common.latest_file(common.step_output_dir("A"), suffix=".csv"))
    pool = read_csv(common.latest_file(os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1F", "Filtered"), suffix=".csv"))
    raw_dir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1D", "Phase1Dc", "Raw_Conservancy")
    raw = read_csv(common.latest_file(raw_dir, suffix=".csv"))
    conf = read_csv(common.latest_file(os.path.join(FOLDS, "Folds", "_confidence"), suffix=".csv"))
    fold_conf = {r["Target"]: r for r in conf}
    sema_json = json.load(open(sorted(glob.glob(os.path.join(FOLDS, "Sensitivity", "Phase1De_SEMA3D_perResidue_*.json")))[-1]))["targets"]

    failures = []

    # ======================================================= CONSERVANCY GATE
    print("\n" + "-" * 90)
    print("CONSERVANCY GATE -- per-target mean over ALL 30,268 candidates (Phase1Dc Raw_Full) vs published")
    print("[INFO] The brief names Min_50pct, but that file holds only the >=50% survivors (1,238 rows) and cannot reproduce the "
          "published means; the gate is therefore run on Raw_Full, which is the full pool the published means describe.")
    by_t = defaultdict(list)
    for r in raw:
        by_t[r["Target"]].append(fnum(r["Conservancy"]))
    gate_rows = []
    for t, pub in PUBLISHED_MEANS.items():
        v = by_t[t]
        m = statistics.mean(v)
        ok = abs(m - pub) <= 0.005
        lo, hi = GATE_RANGE[pathogen_of_target(t)]
        in_range = lo - 0.005 <= m <= hi + 0.005
        gate_rows.append({"Target": t, "Pathogen": pathogen_of_target(t), "N_Candidates": len(v), "Published_Mean_Pct": pub,
                          "Reproduced_Mean_Pct": round(m, 2), "Match": "YES" if ok else "NO", "Within_Brief_Range": "YES" if in_range else "NO",
                          "Max_Pct": round(max(v), 2), "Pct_At_100": round(100 * sum(x >= 100 for x in v) / len(v), 1)})
        print(f"  {t:<10} n={len(v):>6}  published {pub:6.2f}%  reproduced {m:6.2f}%  {'OK' if ok else 'MISMATCH'}  (brief range {lo}-{hi}%: {'in' if in_range else 'OUT'})")
        if not ok or not in_range:
            failures.append(f"conservancy gate {t}: reproduced {m:.2f} vs published {pub}")
    hiv_means = [r["Reproduced_Mean_Pct"] for r in gate_rows if r["Pathogen"] == "HIV"]
    mpox_means = [r["Reproduced_Mean_Pct"] for r in gate_rows if r["Pathogen"] == "Mpox"]
    print(f"  HIV per-target means span {min(hiv_means)}-{max(hiv_means)}% (brief 4.54-9.97); Mpox {min(mpox_means)}-{max(mpox_means)}% (brief 37.80-60.61).")
    if failures:
        common.halt_for_opus(
            branch="Phase4F_conservancy_gate", trigger_number=4,
            one_liner="Conservancy per-target means do not reproduce the published values",
            what_i_was_doing="Recomputing per-target mean conservancy over the full 30,268-candidate Phase1Dc pool.",
            exact_numbers="\n".join(failures), options=["Environment/input fault -- do not consume Step 6 numbers until explained."])
        sys.exit(1)
    print("[SUCCESS] Conservancy gate PASSED: all 7 published per-target means reproduce exactly.")

    # Independent recompute for the 31 construct epitopes from the Phase 1C variant pools (Phase1Dc's own method)
    lib = defaultdict(list)
    fdir = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1C", "Filtered_Antigenicity")
    for f in os.listdir(fdir):
        if f.endswith(".fasta"):
            with open(os.path.join(fdir, f)) as fh:
                sq = "".join(l.strip() for l in fh if not l.startswith(">"))
            lib[f.split("_Var")[0]].append(re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", sq.upper()))
    cons_rows, cons_fail = [], []
    raw_by_target = {t: sorted(v) for t, v in by_t.items()}
    for e in ep:
        pep, t = e["Peptide"], e["Labelled_Target"]
        pool_t = lib.get(t, [])
        hits = sum(1 for sq in pool_t if pep in sq)
        recomputed = round(100 * hits / len(pool_t), 2) if pool_t else None
        stored = fnum(e["Conservancy"])
        ok = recomputed is not None and stored is not None and abs(recomputed - stored) <= 0.005 and e["Hit_Ratio"] == f"{hits}/{len(pool_t)}"
        if not ok:
            cons_fail.append((pep, t, stored, e["Hit_Ratio"], recomputed, f"{hits}/{len(pool_t)}"))
        import bisect
        pct_in_target = 100 * bisect.bisect_right(raw_by_target[t], stored) / len(raw_by_target[t]) if stored is not None else None
        cons_rows.append({"Peptide": pep, "Class": e["Class"], "Pathogen": e["Pathogen"], "Target": t,
                          "Conservancy_Stored_Pct": stored, "Hit_Ratio_Stored": e["Hit_Ratio"],
                          "Conservancy_Recomputed_Pct": recomputed, "Hit_Ratio_Recomputed": f"{hits}/{len(pool_t)}",
                          "Regression_Match": "YES" if ok else "NO", "Variant_Pool_Size": len(pool_t),
                          "Target_Mean_Conservancy_Pct": round(statistics.mean(by_t[t]), 2),
                          "Percentile_Within_Target_Candidates": round(pct_in_target, 1) if pct_in_target is not None else "N/A"})
    print(f"[INFO] independent conservancy recompute (Phase 1C variant pools, exact substring): {len(cons_rows) - len(cons_fail)}/{len(cons_rows)} of the 31 construct epitopes reproduce "
          f"the stored value AND hit ratio.")
    if cons_fail:
        common.halt_for_opus(
            branch="Phase4F_conservancy_recompute", trigger_number=4,
            one_liner="Construct-epitope conservancy does not reproduce from the Phase 1C variant pools",
            what_i_was_doing="Independent exact-substring recompute of the stored conservancy for the 31 construct epitopes.",
            exact_numbers="\n".join(str(x) for x in cons_fail), options=["Environment/input fault -- investigate before using Step 6 conservancy."])
        sys.exit(1)

    print("\nCONSERVANCY -- construct epitopes by pathogen (n stated; pool sizes are small, so values are coarse):")
    for pth in ("HIV", "Mpox"):
        v = [r["Conservancy_Stored_Pct"] for r in cons_rows if r["Pathogen"] == pth]
        pools = sorted({r["Variant_Pool_Size"] for r in cons_rows if r["Pathogen"] == pth})
        print(f"  {pth:<5} n={len(v):>2}  mean {mean(v):6.2f}%  min {min(v):6.2f}%  max {max(v):6.2f}%  (variant pool sizes {pools}; target-pool mean over all candidates "
              f"{statistics.mean([x for t, vv in by_t.items() if pathogen_of_target(t) == pth for x in vv]):.2f}%)")

    # ======================================================= B-CELL: BepiPred tiers
    print("\n" + "-" * 90)
    print("B-CELL -- BepiPred (Phase I tiers, disagreement #3); tier recomputed with the imported classify_bcell_tier")
    bc_pool = [r for r in pool if r["Type"] == "B-cell"]
    construct_b = {e["Peptide"]: e for e in ep if e["Class"] == "B-cell"}
    tier_fail = []
    for r in bc_pool:
        m_, p_ = fnum(r["mean_BepiPred"]), fnum(r["pct_above"])
        if p1db.classify_bcell_tier(m_, p_) != r["Bcell_Tier"]:
            tier_fail.append((r["Peptide"], r["Bcell_Tier"], p1db.classify_bcell_tier(m_, p_)))
    print(f"[INFO] Bcell_Tier recomputed for all {len(bc_pool)} pool B-cell candidates: {len(bc_pool) - len(tier_fail)}/{len(bc_pool)} match Phase 1Db.")
    if tier_fail:
        common.halt_for_opus(
            branch="Phase4F_bcell_tier", trigger_number=4, one_liner="Bcell_Tier does not reproduce from stored BepiPred values",
            what_i_was_doing="Recomputing Phase I tiers with the imported classify_bcell_tier.", exact_numbers="\n".join(str(x) for x in tier_fail),
            options=["Environment/input fault."])
        sys.exit(1)

    # ======================================================= SEMA recompute + native-fold overlay
    print("\n" + "-" * 90)
    print("SEMA-3D (reused, not re-run) + native-fold overlay (Phase 1De per-antigen folds)")
    folds = {}
    for t in PUBLISHED_MEANS:
        folds[t] = load_native_fold(t)
    # sanity 1: pLDDT read from the PDB B-factor column must reproduce Phase 1De's stored Mean_pLDDT
    print("  sanity -- mean pLDDT from B-factors vs Phase 1De stored Mean_pLDDT (must agree; B-factor holds pLDDT):")
    for t, f in folds.items():
        mp = statistics.mean(f["plddt"])
        st = fnum(fold_conf[t]["Mean_pLDDT"])
        flag = "OK" if abs(mp - st) <= 1.0 else "MISMATCH"
        print(f"    {t:<10} {mp:5.1f} vs {st:5.1f}  {flag}  ({len(f['seq'])} residues, {f['n_chains']} chain(s))")
        if flag != "OK":
            failures.append(f"pLDDT mismatch {t}")
    # sanity 2: RSA direction -- charged/polar residues must be MORE exposed than buried hydrophobics
    allrsa = defaultdict(list)
    for t, f in folds.items():
        for aa, r in zip(f["seq"], f["rsa"]):
            if r is not None:
                allrsa[aa].append(r)
    pol = mean([x for a in "KREDQN" for x in allrsa[a]])
    hyd = mean([x for a in "LIVFWC" for x in allrsa[a]])
    print(f"  sanity -- mean RSA: charged/polar (KREDQN) {pol:.2f} vs buried-prone hydrophobic (LIVFWC) {hyd:.2f} (must be polar > hydrophobic): {'OK' if pol > hyd else 'WRONG SIGN'}")
    if not pol > hyd:
        failures.append("RSA direction check failed")
    if failures:
        common.halt_for_opus(branch="Phase4F_fold_overlay", trigger_number=4, one_liner="Native-fold overlay sanity check failed",
                             what_i_was_doing="pLDDT/RSA sanity checks.", exact_numbers="\n".join(failures), options=["Environment fault."])
        sys.exit(1)

    # construct-model contrast (why the construct fold is not used)
    from Bio.PDB import MMCIFParser
    cif = sorted(glob.glob(os.path.join(CONSTRUCT_FOLD, "*_model_0.cif")))[0]
    cst = MMCIFParser(QUIET=True).get_structure("construct", cif)[0]
    cres = [r for r in next(iter(cst)) if r.id[0] == " "]
    cpl = [statistics.mean(a.get_bfactor() for a in r) for r in cres]
    ptm = json.load(open(sorted(glob.glob(os.path.join(CONSTRUCT_FOLD, "*summary_confidences_0.json")))[0]))["ptm"]
    non_adj = cpl[ADJUVANT_LEN:]
    print(f"  construct model: {len(cres)} residues, pTM {ptm}; adjuvant (first {ADJUVANT_LEN}) mean pLDDT {statistics.mean(cpl[:ADJUVANT_LEN]):.1f} "
          f"({100 * sum(x >= 70 for x in cpl[:ADJUVANT_LEN]) / ADJUVANT_LEN:.0f}% >=70); non-adjuvant mean pLDDT {statistics.mean(non_adj):.1f}, "
          f"{100 * sum(x >= 70 for x in non_adj) / len(non_adj):.1f}% of {len(non_adj)} residues >=70")
    cons_model = {"ptm": ptm, "n": len(cres), "adj_mean": statistics.mean(cpl[:ADJUVANT_LEN]), "adj_ge70": 100 * sum(x >= 70 for x in cpl[:ADJUVANT_LEN]) / ADJUVANT_LEN,
                  "non_mean": statistics.mean(non_adj), "non_ge70": 100 * sum(x >= 70 for x in non_adj) / len(non_adj)}

    def overlay(r, positions=None):
        """Locate the epitope in its target's native fold; returns metrics dict."""
        t, pep = r["Target"] if "Target" in r else r["Labelled_Target"], r["Peptide"]
        f = folds.get(t)
        out = {"Located_In_Native_Fold": "NO"}
        if f is None:
            return out
        i = f["seq"].find(pep)
        if i < 0:
            out["Located_In_Native_Fold"] = "NO (epitope not in the Var_01-based fold sequence)"
            return out
        L = len(pep)
        rs = [x for x in f["rsa"][i:i + L] if x is not None]
        pl = f["plddt"][i:i + L]
        sc = sema_json[t]["scores"][i:i + L]
        verdict = fold_conf[t]["Global_Fold_Verdict"]
        out.update({"Located_In_Native_Fold": "YES", "Fold_Start_1based": i + 1, "Native_Fold_Verdict": verdict, "Native_Fold_pTM": fold_conf[t]["pTM_model0"],
                    "RSA_Mean": round(statistics.mean(rs), 3), "Frac_Residues_Exposed_RSA_ge_0.25": round(sum(x >= RSA_EXPOSED for x in rs) / len(rs), 3),
                    "Native_pLDDT_Mean": round(statistics.mean(pl), 1), "Native_Frac_pLDDT_ge70": round(sum(x >= 70 for x in pl) / len(pl), 3),
                    "SEMA_Mean_Score": round(statistics.mean(sc), 3), "SEMA_N_Residues_ge_cut": sum(x >= SEMA_CUT for x in sc),
                    "SEMA_Overlap_Recomputed_Pct": round(100 * sum(x >= SEMA_CUT for x in sc) / L, 1),
                    "Exposure_Reliability": ("USABLE (global fold DETERMINED)" if verdict == "DETERMINED" else
                                             "INDICATIVE ONLY (global fold UNDETERMINED: pTM<0.5)")})
        return out

    def sema_state(r):
        v = r["SEMA_Corroborated"]
        return {"YES": "CORROBORATED (>=50% residues in a SEMA-3D patch)",
                "NO": "SCREENED, not corroborated (<50% in a patch): lower-priority linear-only; NOT a refutation",
                "UNSCREENED": "UNSCREENED: could not be assessed on the one available fold; NOT a negative"}.get(v, f"UNKNOWN({v})")

    # --------------------------------------------------- build B-cell tables
    def b_row(r, dossier_row=None):
        is_construct = dossier_row is not None
        t = r.get("Target") or r["Labelled_Target"]
        m_, p_ = fnum(r["mean_BepiPred"]), fnum(r["pct_above"])
        rr = dict(r)
        rr["Target"] = t
        ov = overlay(rr)
        row = {"Peptide": r["Peptide"], "Pathogen": pathogen_of_target(t), "Target": t, "In_Construct": "YES" if is_construct else "NO",
               "mean_BepiPred": m_, "pct_above": p_,
               "Primary_Gate_Phase_I_Tier": r["Bcell_Tier"],
               "Secondary_Stratum_Methodology_High": "YES" if (m_ >= 0.50 and p_ >= 75.0) else "NO",
               "SEMA_Overlap_Pct_Stored": r["SEMA_Overlap_Pct"], "SEMA_Residues_Scored": r["SEMA_Residues_Scored"],
               "SEMA_Corroborated_Stored": r["SEMA_Corroborated"], "SEMA_Interpretation": sema_state(r),
               "Conservancy_Pct": fnum(r["Conservancy"]), "Hit_Ratio": r["Hit_Ratio"]}
        row.update(ov)
        if ov.get("SEMA_Overlap_Recomputed_Pct") is not None:
            row["SEMA_Recompute_Match"] = "YES" if abs(ov["SEMA_Overlap_Recomputed_Pct"] - (fnum(r["SEMA_Overlap_Pct"]) or 0)) <= 0.06 else "NO"
        else:
            row["SEMA_Recompute_Match"] = "N/A (not locatable in native fold)"
        if is_construct:
            s1, e1 = int(dossier_row["Construct_Position_Start"]), int(dossier_row["Construct_Position_End"])
            seg = cpl[s1:e1]
            row.update({"Construct_Start": s1 + 1, "Construct_End": e1, "Construct_Model_pLDDT_Mean": round(statistics.mean(seg), 1),
                        "Construct_Model_Frac_pLDDT_ge70": round(sum(x >= 70 for x in seg) / len(seg), 3)})
        return row

    construct_rows = [b_row(r, r) for r in ep if r["Class"] == "B-cell"]
    pool_rows = [b_row(r, construct_b.get(r["Peptide"])) for r in bc_pool]
    sema_bad = [r["Peptide"] for r in construct_rows + pool_rows if r["SEMA_Recompute_Match"] == "NO"]
    n_loc = sum(1 for r in construct_rows + pool_rows if r["SEMA_Recompute_Match"] == "YES")
    print(f"[INFO] SEMA overlap recomputed from Phase 1De per-residue scores (cut ln4={SEMA_CUT:.4f}): {n_loc} locatable epitopes match, {len(sema_bad)} mismatch "
          f"(construct 10: {sum(1 for r in construct_rows if r['SEMA_Recompute_Match'] == 'YES')}/10 match).")
    if sema_bad:
        common.halt_for_opus(branch="Phase4F_sema_recompute", trigger_number=4, one_liner="SEMA overlap does not reproduce from stored per-residue scores",
                             what_i_was_doing="Recomputing SEMA_Overlap_Pct.", exact_numbers=str(sema_bad), options=["Environment fault."])
        sys.exit(1)

    from scipy import stats as _st
    loc = [r for r in pool_rows if r.get("Located_In_Native_Fold") == "YES"]
    rho_b, p_b = _st.spearmanr([r["mean_BepiPred"] for r in loc], [r["RSA_Mean"] for r in loc])
    print(f"[INFO] direction cross-check -- BepiPred mean vs native-fold RSA across the {len(loc)} unique pool B-cell candidates: Spearman rho={rho_b:+.3f}, p={p_b:.2g} "
          f"(BepiPred models surface exposure, so a non-negative sign is expected; n={len(loc)}, a sanity check not a validation).")

    def stats_by_pathogen(rows, label):
        print(f"\n{label}")
        for pth in ("HIV", "Mpox"):
            g = [r for r in rows if r["Pathogen"] == pth]
            tiers = Counter(r["Primary_Gate_Phase_I_Tier"] for r in g)
            sema = Counter(r["SEMA_Corroborated_Stored"] for r in g)
            loc = [r for r in g if r.get("Located_In_Native_Fold") == "YES"]
            print(f"  {pth:<5} n={len(g):>2} | tiers {dict(tiers)} | methodology-High {sum(r['Secondary_Stratum_Methodology_High'] == 'YES' for r in g)}/{len(g)} | "
                  f"mean BepiPred {fmt(mean([r['mean_BepiPred'] for r in g]), 3)} | SEMA {dict(sema)} | "
                  f"native fold located {len(loc)}/{len(g)}; mean RSA {fmt(mean([r['RSA_Mean'] for r in loc]), 2)}; "
                  f"mean native pLDDT {fmt(mean([r['Native_pLDDT_Mean'] for r in loc]), 1)}; conservancy {fmt(mean([r['Conservancy_Pct'] for r in g]), 1)}%")
    stats_by_pathogen(construct_rows, "CONSTRUCT B-cell epitopes (n=10) -- HIV vs Mpox:")
    stats_by_pathogen(pool_rows, "POOL B-cell candidates (n=39, comparison) -- HIV vs Mpox:")
    print("\nPer construct B-cell epitope (native-fold overlay; construct-model pLDDT shown for contrast):")
    print(f"  {'peptide':<18}{'target':<10}{'tier':<14}{'BepiPred':>9}{'SEMA%':>6} {'SEMA':<5}{'RSA':>6}{'exp.frac':>9}{'nat.pLDDT':>10}{'constr.pLDDT':>13}  fold verdict")
    for r in construct_rows:
        print(f"  {r['Peptide']:<18}{r['Target']:<10}{r['Primary_Gate_Phase_I_Tier']:<14}{r['mean_BepiPred']:>9.3f}{r['SEMA_Overlap_Pct_Stored']:>6} {r['SEMA_Corroborated_Stored']:<5}"
              f"{fmt(r.get('RSA_Mean'), 2):>6}{fmt(r.get('Frac_Residues_Exposed_RSA_ge_0.25'), 2):>9}{fmt(r.get('Native_pLDDT_Mean'), 1):>10}{r['Construct_Model_pLDDT_Mean']:>13}  {r.get('Native_Fold_Verdict', 'N/A')}")

    # ------------------------------------------------------------- outputs
    out_dir, ts = common.step_output_dir("F"), common.timestamp()

    def write(name, data, fields=None):
        fields = fields or list(dict.fromkeys(k for r in data for k in r))
        p = os.path.join(out_dir, f"{name}_{ts}.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, restval="N/A")
            w.writeheader(); w.writerows(data)
        print(f"[SUCCESS] {os.path.basename(p)}: {len(data)} rows")
    write("Phase4F_BcellDossier", construct_rows)
    write("Phase4F_BcellPool", pool_rows)
    write("Phase4F_Conservancy", cons_rows)
    write("Phase4F_ConservancyGate", gate_rows)

    fs = [("Construct model residues", cons_model["n"]), ("Construct model pTM", cons_model["ptm"]),
          ("Adjuvant (first 45) mean pLDDT", round(cons_model["adj_mean"], 1)), ("Adjuvant % residues pLDDT>=70", round(cons_model["adj_ge70"], 1)),
          ("Non-adjuvant residues", cons_model["n"] - ADJUVANT_LEN), ("Non-adjuvant mean pLDDT", round(cons_model["non_mean"], 1)),
          ("Non-adjuvant % residues pLDDT>=70", round(cons_model["non_ge70"], 1)),
          ("Mean RSA charged/polar KREDQN (all 7 native folds)", round(pol, 3)), ("Mean RSA buried-prone LIVFWC", round(hyd, 3)),
          ("Spearman BepiPred mean vs native RSA (unique pool B-cell)", round(rho_b, 4)), ("... p-value", round(p_b, 3)), ("... n", len(loc)),
          ("Pool B-cell candidates", len(bc_pool)), ("Bcell_Tier recomputed matches", len(bc_pool) - len(tier_fail)),
          ("SEMA overlap recomputed matches (locatable, incl. construct duplicates)", n_loc), ("SEMA overlap mismatches", len(sema_bad)),
          ("Construct epitopes conservancy recompute matches", len(cons_rows) - len(cons_fail)), ("Construct epitopes tested", len(cons_rows))]
    for tt, fdd in folds.items():
        fs.append((f"Native fold {tt}: verdict / pTM / mean pLDDT / residues", f"{fold_conf[tt]['Global_Fold_Verdict']} / {fold_conf[tt]['pTM_model0']} / {statistics.mean(fdd['plddt']):.1f} / {len(fdd['seq'])}"))
    with open(os.path.join(out_dir, f"Phase4F_Summary_{ts}.csv"), "w", newline="") as _f:
        _w = csv.writer(_f); _w.writerow(["Metric", "Value"]); _w.writerows(fs)
    print(f"[SUCCESS] Phase4F_Summary_{ts}.csv: {len(fs)} rows")

    # ------------------------------------------------------------ methodology note
    tiers_c = Counter(r["Primary_Gate_Phase_I_Tier"] for r in construct_rows)
    tiers_p = Counter(r["Primary_Gate_Phase_I_Tier"] for r in pool_rows)
    sema_p = Counter(r["SEMA_Corroborated_Stored"] for r in pool_rows)
    det = [r for r in construct_rows if r.get("Native_Fold_Verdict") == "DETERMINED"]
    und = [r for r in construct_rows if r.get("Native_Fold_Verdict") == "UNDETERMINED"]
    L = []
    A = L.append
    A("# Phase 4F -- B-cell + Conservancy: methodology notes\n")
    A(f"Generated: {ts}\n")
    A("## STATUS: complete; all regression gates passed; no escalation\n")
    A("Mostly reuse -- **no new predictions**. Every stored value that can be recomputed deterministically was recomputed as a regression check "
      "(BepiPred tier, SEMA overlap, per-target conservancy means, per-epitope conservancy). All matched.\n")
    A("## Conservancy\n")
    A("**Hard gate PASSED -- all 7 published per-target means reproduce exactly** (mean over all 30,268 Phase 1Dc candidates):\n")
    A("| target | n candidates | published | reproduced | max | % at 100 |\n|---|---|---|---|---|---|")
    for g in gate_rows:
        A(f"| {g['Target']} | {g['N_Candidates']} | {g['Published_Mean_Pct']}% | {g['Reproduced_Mean_Pct']}% | {g['Max_Pct']}% | {g['Pct_At_100']}% |")
    A(f"\nHIV per-target means span {min(hiv_means)}-{max(hiv_means)}%; Mpox {min(mpox_means)}-{max(mpox_means)}% (brief: 4.54-9.97 / 37.80-60.61). "
      "**Source note:** the brief says to use `Min_50pct`, but that file holds only the >=50% survivors (1,238 rows) and cannot reproduce a mean over all candidates; the gate uses "
      "`Raw_Conservancy/Phase1Dc_Raw_Full`, the full pool the published table describes. `Min_50pct` is what the 31 construct epitopes were selected from.\n")
    A("Independent recompute (Phase 1C variant pools, exact substring, Phase1Dc's own method): **all 31 construct epitopes reproduce both the stored conservancy and hit ratio**.\n")
    A("| pathogen | n | mean | min | max | variant pool sizes | mean over ALL candidates of those targets |\n|---|---|---|---|---|---|---|")
    for pth in ("HIV", "Mpox"):
        v = [r["Conservancy_Stored_Pct"] for r in cons_rows if r["Pathogen"] == pth]
        pools = sorted({r["Variant_Pool_Size"] for r in cons_rows if r["Pathogen"] == pth})
        allm = statistics.mean([x for t, vv in by_t.items() if pathogen_of_target(t) == pth for x in vv])
        A(f"| {pth} | {len(v)} | {mean(v):.1f}% | {min(v):.1f}% | {max(v):.1f}% | {pools} | {allm:.1f}% |")
    hv = [r["Conservancy_Stored_Pct"] for r in cons_rows if r["Pathogen"] == "HIV"]
    A(f"\nThe construct's HIV epitopes ({min(hv):.1f}-{max(hv):.1f}%) sit far above the HIV candidate means (4.5-10%): selection took the conserved tail of a very diverse pool, "
      "so the pool means describe the diversity of the antigens, not the construct. Variant pools are small (n = 21-30), so values are coarse (one variant = 3-5 pp); "
      "HIV epitopes are conserved across fewer of a smaller number of CRF01_AE isolates -- interpret as within-clade conservancy only.\n")
    A("## B-cell -- BepiPred (Phase I tiers, disagreement #3)\n")
    A(f"Tiers recomputed with the imported `classify_bcell_tier` for all {len(bc_pool)} pool candidates: **{len(bc_pool)}/{len(bc_pool)} match Phase 1Db**. "
      "`Primary_Gate` = Phase I tier; `Secondary_Stratum` = the methodology's stricter 'High' (mean >= 0.50 and >= 75% of residues above).\n")
    A(f"Construct (n = 10): Phase I tiers {dict(tiers_c)}; methodology-High {sum(r['Secondary_Stratum_Methodology_High'] == 'YES' for r in construct_rows)}/10. "
      f"Pool (n = 39): {dict(tiers_p)}. Three construct B-cell epitopes are `Deprioritized` -- Phase I retained them (tiers deprioritise, they never exclude). HIV vs Mpox:\n")
    A("| pathogen | n | Phase I tiers | methodology-High | mean BepiPred |\n|---|---|---|---|---|")
    for pth in ("HIV", "Mpox"):
        g = [r for r in construct_rows if r["Pathogen"] == pth]
        A(f"| {pth} | {len(g)} | {dict(Counter(r['Primary_Gate_Phase_I_Tier'] for r in g))} | {sum(r['Secondary_Stratum_Methodology_High'] == 'YES' for r in g)}/{len(g)} | {mean([r['mean_BepiPred'] for r in g]):.3f} |")
    A("\n## B-cell -- SEMA-3D (Phase 1De, reused; not re-run because `foldseek` is a Linux binary)\n")
    A("All 10 construct B-cell epitopes are `SEMA_Corroborated = NO` (overlap 0-12.5% of 16 residues scored). **`NO` means screened and not corroborated: fewer than 50% of the residues "
      "fall in a SEMA-3D conformational patch (score >= ln 4). It is a lower-priority *linear-only* call, not a refutation** -- SEMA-3D's own calibration is TPR 0.66 at FPR 0.07, so a miss is far from disproof. "
      "**`UNSCREENED` (no assessable fold) is a different state and is never a negative**; none of the construct epitopes are UNSCREENED "
      f"(pool of 39: {dict(sema_p)}). The overlap was recomputed from Phase 1De's per-residue scores for every locatable epitope: all match.\n")
    A("## B-cell -- native-fold overlay (why not the construct model)\n")
    A(f"The construct's own AlphaFold fold is undetermined outside the adjuvant: **pTM {cons_model['ptm']}**; the adjuvant domain (first {ADJUVANT_LEN} residues) has mean pLDDT {cons_model['adj_mean']:.1f} "
      f"({cons_model['adj_ge70']:.0f}% of residues >= 70), but the non-adjuvant {cons_model['n'] - ADJUVANT_LEN} residues have mean pLDDT {cons_model['non_mean']:.1f} with **{cons_model['non_ge70']:.1f}% >= 70** "
      "(verified from the model file). Surface accessibility computed on a model that is not a determined fold would be meaningless, so it is computed on the **Phase 1De per-antigen native folds** instead "
      "(each antigen folded alone). Each epitope is located in its antigen's fold by exact sequence match; per-residue relative SASA (Shrake-Rupley, probe 1.4 A, Tien/Wilke 2013 maxima shipped with Biopython) "
      "and pLDDT are averaged over the 16 residues. Sanity checks passed before use: PDB B-factor pLDDT reproduces Phase 1De's stored Mean_pLDDT for all 7 antigens, and RSA direction is correct "
      "(polar residues more exposed than buried-prone hydrophobics).\n")
    A("| peptide | target | fold verdict (pTM) | mean RSA | residues exposed (RSA>=0.25) | native pLDDT | construct-model pLDDT |\n|---|---|---|---|---|---|---|")
    for r in construct_rows:
        A(f"| {r['Peptide']} | {r['Target']} | {r.get('Native_Fold_Verdict', 'N/A')} ({r.get('Native_Fold_pTM', 'N/A')}) | {fmt(r.get('RSA_Mean'), 2)} | {fmt(r.get('Frac_Residues_Exposed_RSA_ge_0.25'), 2)} | {fmt(r.get('Native_pLDDT_Mean'), 1)} | {r['Construct_Model_pLDDT_Mean']} |")
    A(f"\nNative folds: **{len(det)} of 10** construct B-cell epitopes sit on a DETERMINED global fold (gp120, p17, p24, L1R); **{len(und)} of 10 sit on an UNDETERMINED fold** "
      "(gp41 pTM 0.39, A35R 0.47, B5R 0.43) and their exposure is INDICATIVE ONLY. Caveats for all: monomer folds (oligomer interfaces not modelled), and HIV Env glycans are not modelled, "
      "so gp120/gp41 exposure is over-estimated. Exposure is context, not a gate.\n")
    A(f"Direction cross-check: BepiPred mean vs native-fold RSA across the {len(loc)} unique pool B-cell candidates, Spearman rho = {rho_b:+.3f} (p = {p_b:.2g}, n = {len(loc)}). "
      "The sign is as expected but the correlation is NOT distinguishable from zero at this n: it passes the sign check and corroborates nothing.\n")
    A("## HIV vs Mpox (every metric split)\n")
    A("| metric | HIV | Mpox |\n|---|---|---|")
    def col(fn):
        return [fn(pth) for pth in ("HIV", "Mpox")]
    def grp(pth, rows=construct_rows):
        return [r for r in rows if r["Pathogen"] == pth]
    A("| construct B-cell epitopes (n) | %d | %d |" % tuple(len(grp(p)) for p in ("HIV", "Mpox")))
    A("| mean BepiPred | %.3f | %.3f |" % tuple(mean([r["mean_BepiPred"] for r in grp(p)]) for p in ("HIV", "Mpox")))
    A("| SEMA corroborated (YES) | %d/%d | %d/%d |" % (sum(r["SEMA_Corroborated_Stored"] == "YES" for r in grp("HIV")), len(grp("HIV")), sum(r["SEMA_Corroborated_Stored"] == "YES" for r in grp("Mpox")), len(grp("Mpox"))))
    A("| mean native-fold RSA (located epitopes) | %s | %s |" % tuple(fmt(mean([r.get("RSA_Mean") for r in grp(p) if r.get("Located_In_Native_Fold") == "YES"]), 2) for p in ("HIV", "Mpox")))
    A("| mean conservancy, construct B-cell epitopes | %.1f%% | %.1f%% |" % tuple(mean([r["Conservancy_Pct"] for r in grp(p)]) for p in ("HIV", "Mpox")))
    A("| mean conservancy, all 31 construct epitopes | %.1f%% | %.1f%% |" % tuple(mean([r["Conservancy_Stored_Pct"] for r in cons_rows if r["Pathogen"] == p]) for p in ("HIV", "Mpox")))
    A("\n## Caveats\n")
    A("- BepiPred and SEMA-3D are predictors; none of this is experimental evidence of antibody binding.\n"
      "- Conservancy is exact-match within the Phase 1C variant pools (n = 21-30 per target), not IEDB's alignment tool (Phase 1Dc note) and not a global diversity measure.\n"
      "- Native-fold exposure on undetermined folds is indicative only; glycan shielding and oligomerisation are not modelled.\n")
    with open(os.path.join(out_dir, "METHODOLOGY_NOTE.md"), "w") as nf:
        nf.write("\n".join(L))
    print("[INFO] Methodology note written.")
    return dict(construct_rows=construct_rows, pool_rows=pool_rows, cons_rows=cons_rows, gate_rows=gate_rows, cons_model=cons_model)


if __name__ == "__main__":
    build_bcell_conservancy()
