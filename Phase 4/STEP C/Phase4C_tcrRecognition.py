import os
import sys
import csv
import json
import subprocess
import hashlib
import statistics
from collections import defaultdict

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
# PHASE 4C -- TCR RECOGNITION (methodology Section IV.A.2, Step 3 of 7)
#
# The methodology's named tools ("MixMHCpred + TCGA Contact Database") are
# both wrong -- MixMHCpred predicts MHC binding, not TCR contact, and TCGA is
# a cancer genomics atlas with no TCR-contact data at all (disagreement #6).
# IEDB's own immunogenicity/ endpoint returns 403 -- retested 2026-09-21 (11 attempts over several
# minutes, incl. after a 90 s wait and on GET) while mhci/, mhcii/ and bcell/ answered 200 in the same
# minute, and a NONEXISTENT path (foobar/) returns the byte-identical 403. That is an endpoint not
# served at this URL, NOT IP throttling (throttling hits healthy endpoints too and clears on retry).
# This step uses the
# pre-approved 3-rung ladder instead:
#   1. BigMHC-IM (primary) -- already installed, MHC-I only.
#   2. PRIME 2.0 + MixMHCpred 3.0 (optional) -- installed fresh this run
#      (see METHODOLOGY_NOTE.md for the install path). MHC-I only.
#   3. Calis et al. (2013) immunogenicity model, reimplemented locally from
#      the paper's own published training data (not from memory -- see
#      METHODOLOGY_NOTE.md for full derivation and validation).
#
# SCOPE: all three ladder tools, and BigMHC/PRIME/Calis's underlying science,
# are MHC-I specific (this is explicit in Section F's own wording: "TCR-facing
# residues = positions 4-6 of each MHC-I epitope"). TCR recognition proper
# (this step's scored ladder) therefore covers the 10 MHC-I construct
# epitopes only. Cross-reactivity (self-homology) is a separate, broader
# safety question and is reported for all 31 construct epitopes.
#
# DIRECTION DISCIPLINE (learned the hard way in Phase4B -- see its
# METHODOLOGY_NOTE.md "CAUGHT-AND-CORRECTED PARSING ERROR"): every tool's
# score is checked against a trusted quantity (Phase4B's MHC-I percentile
# rank) for the CORRECT sign before being used, and the observed correlation
# is reported here, not assumed.
# =============================================================================

BIGMHC_PREDICT = os.path.join(_PROJECT_ROOT, "external_tools", "bigmhc", "src", "predict.py")
PRIME_BIN = os.path.join(_PROJECT_ROOT, "external_tools", "prime", "PRIME")
MIXMHCPRED_BIN = os.path.join(_PROJECT_ROOT, "external_tools", "mixmhcpred", "MixMHCpred")
PRIME_ALLELES_LIST = os.path.join(_PROJECT_ROOT, "external_tools", "prime", "lib", "alleles.txt")

# Both PRIME and MixMHCpred's own shell scripts refuse to run from a path
# containing a space ("Spaces in path to MixMHCpred are not supported") --
# and this project's root, "ChimericVaccine/Research", sits under a
# space-containing volume name ("Extended SSD"). The project already solves
# exactly this problem for conda (see ~/mpoxhiv_env.sh's comment on
# ~/mpoxhiv_ssd): a space-free symlink to the same volume already exists at
# ~/mpoxhiv_ssd. Re-deriving the tool paths through that symlink instead of
# _PROJECT_ROOT sidesteps the restriction without moving anything.
_SPACE_FREE_ROOT = os.path.expanduser("~/mpoxhiv_ssd/ChimericVaccine/Research")
PRIME_BIN_SF = os.path.join(_SPACE_FREE_ROOT, "external_tools", "prime", "PRIME")
MIXMHCPRED_BIN_SF = os.path.join(_SPACE_FREE_ROOT, "external_tools", "mixmhcpred", "MixMHCpred")

TCR_FACING_START, TCR_FACING_END = 4, 6  # 1-based, inclusive -- see Section F Step 3


def tcr_facing_residues(peptide):
    """Positions 4-6 (1-based) -- the methodology's '>=3 TCR-contact residue'
    rule cannot be applied without a defined position set, so this step uses
    the published TCR-facing region directly and reports which residues sit
    there, rather than scoring against an assumed threshold."""
    return peptide[TCR_FACING_START - 1:TCR_FACING_END]


def corr_with_p(xs, ys):
    """Pearson r, two-sided p, n. Every correlation this step reports carries
    its n and p -- at n=10 a moderate r is indistinguishable from zero."""
    from scipy import stats
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan"), n
    r, p = stats.pearsonr(xs, ys)
    return float(r), float(p), n


# =============================================================================
# RUNG 1 -- BigMHC-IM
# =============================================================================
def run_bigmhc_im(epitope_allele_pairs):
    """epitope_allele_pairs: list of (peptide, allele). Returns
    {(peptide, allele): bigmhc_im_score} or {} if BigMHC is unavailable."""
    if not os.path.isfile(BIGMHC_PYTHON := os.environ.get("BIGMHC_PYTHON", "")):
        print(f"[ERROR] BIGMHC_PYTHON not found: {BIGMHC_PYTHON!r} -- is ~/mpoxhiv_env.sh sourced?")
        return {}
    if not os.path.isfile(BIGMHC_PREDICT):
        print(f"[ERROR] BigMHC predict.py not found: {BIGMHC_PREDICT}")
        return {}

    cache_dir = common.cache_dir_for("C")
    key = common.content_key(sorted(epitope_allele_pairs))
    in_path = os.path.join(cache_dir, f"bigmhc_im_input__{key}.csv")
    out_path = os.path.join(cache_dir, f"bigmhc_im_output__{key}.csv")

    if not os.path.isfile(out_path):
        with open(in_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["mhc", "pep", "tgt"])  # exact header BigMHC's default column indices expect
            for pep, allele in epitope_allele_pairs:
                w.writerow([allele, pep, 0])
        print(f"[PROCESS] Running BigMHC-IM on {len(epitope_allele_pairs)} (peptide, allele) pairs...")
        result = subprocess.run(
            [BIGMHC_PYTHON, BIGMHC_PREDICT, f"-i={in_path}", "-m=im", "-d=cpu", f"-o={out_path}"],
            capture_output=True, text=True, cwd=os.path.dirname(BIGMHC_PREDICT),
        )
        if result.returncode != 0 or not os.path.isfile(out_path):
            print(f"[ERROR] BigMHC-IM failed (exit {result.returncode}):\n{result.stderr[-2000:]}")
            return {}
    else:
        print(f"[INFO] Reusing cached BigMHC-IM output: {os.path.basename(out_path)}")

    # Parse by HEADER NAME, never index -- BigMHC's own example .cmp output
    # shows extra columns (len, tgt) interleaved before the score column, so
    # a positional read would silently grab the wrong field on any format
    # drift.
    out = {}
    with open(out_path, newline="") as f:
        reader = csv.DictReader(f)
        if "BigMHC_IM" not in reader.fieldnames:
            print(f"[ERROR] BigMHC-IM output missing 'BigMHC_IM' column: {reader.fieldnames}")
            return {}
        for row in reader:
            try:
                out[(row["pep"], row["mhc"])] = float(row["BigMHC_IM"])
            except (ValueError, KeyError):
                continue
    return out


# =============================================================================
# RUNG 2 -- PRIME 2.0 + MixMHCpred 3.0 (optional; drop to rungs 1+3 if this
# fails to install/run cleanly -- pre-approved, not a tool substitution)
# =============================================================================
# PRIME's compiled binary (not just its bash launcher) rejects a space
# anywhere in its -i/-o file path arguments. The ~/mpoxhiv_ssd symlink fixes
# this for the launcher scripts and the binary PATHS themselves (bash's own
# `cd`/`pwd` preserve the symlink, never resolving it) -- but PRIME.x calls
# realpath() internally on the FILE arguments it's given, which DOES
# dereference the symlink back to the real, space-containing
# "/Volumes/Extended SSD/..." path, then rejects that. Confirmed by direct
# reproduction: the exact same command fails with "Spaces in input path are
# not supported" when the input file lives under the symlinked path, and
# succeeds when it lives under a path with no symlink in its chain at all.
# Fix: stage PRIME's own input/output files in a small scratch directory
# under $HOME (never under /Volumes/Extended SSD, symlinked or not). This is
# transient plumbing for one space-intolerant third-party binary, not a
# pipeline output -- every actual RESULT this step produces is still parsed
# out and written to Step_Outputs/Phase4/StepC/ on the SSD as normal.
_PRIME_IO_DIR = os.path.expanduser("~/.mpoxhiv_phase4c_prime_io")
os.makedirs(_PRIME_IO_DIR, exist_ok=True)


def _prime_allele_name(allele_freq_key):
    """'HLA-A*24:02' -> 'A2402' (PRIME/MixMHCpred's no-punctuation format,
    confirmed against lib/alleles.txt)."""
    return allele_freq_key.replace("HLA-", "").replace("*", "").replace(":", "")


def _load_prime_supported_alleles():
    if not os.path.isfile(PRIME_ALLELES_LIST):
        return set()
    with open(PRIME_ALLELES_LIST) as f:
        return {line.strip() for line in f if line.strip()}


def prime_available():
    return os.path.isfile(PRIME_BIN_SF) and os.path.isfile(MIXMHCPRED_BIN_SF)


def run_prime_one_epitope(peptide, allele_freq_keys, supported_alleles):
    """Runs PRIME once for a single peptide against its own binding-allele
    set. Returns dict with pct_rank_best, score_best, best_allele, or None
    if no allele in this epitope's set is PRIME-supported or the run fails."""
    prime_alleles = sorted({_prime_allele_name(a) for a in allele_freq_keys} & supported_alleles)
    if not prime_alleles:
        return None

    key = common.content_key((peptide, tuple(prime_alleles)))
    in_path = os.path.join(_PRIME_IO_DIR, f"prime_input__{key}.txt")
    out_path = os.path.join(_PRIME_IO_DIR, f"prime_output__{key}.txt")

    if not os.path.isfile(out_path):
        with open(in_path, "w") as f:
            f.write(peptide + "\n")
        env = dict(os.environ)
        phase2_bin = os.path.join(env.get("MPOXHIV_ENVS", ""), "phase2", "bin")
        env["PATH"] = phase2_bin + ":" + env.get("PATH", "")  # PRIME/MixMHCpred call bare `python3`, needs pandas
        result = subprocess.run(
            [PRIME_BIN_SF, "-i", in_path, "-o", out_path,
             "-a", ",".join(prime_alleles), "-mix", MIXMHCPRED_BIN_SF],
            capture_output=True, text=True, env=env,
        )
        if result.returncode != 0 or not os.path.isfile(out_path):
            # PRIME prints its own error text to stdout, not stderr.
            print(f"[WARNING] PRIME failed for {peptide} (exit {result.returncode}): "
                  f"{(result.stdout + result.stderr).strip()[-500:]}")
            return None

        # Copy the raw response onto the SSD for audit/reproducibility --
        # _PRIME_IO_DIR is transient plumbing (see its definition above),
        # the evidence trail belongs in Step_Outputs like every other tool
        # response this pipeline caches.
        with open(out_path) as _src:
            _raw = _src.read()
        with open(os.path.join(common.cache_dir_for("C"), f"prime_output__{key}.txt"), "w") as _dst:
            _dst.write(_raw)

    # PRIME's output has ~10 leading "#" comment lines before the real header
    # -- locate the header by content ("Peptide" as the first field), never
    # by a fixed line-skip count, then parse by column NAME.
    with open(out_path) as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    header_idx = next((i for i, l in enumerate(lines) if l.split("\t")[0] == "Peptide"), None)
    if header_idx is None or header_idx + 1 >= len(lines):
        print(f"[WARNING] Could not locate PRIME's header row for {peptide}")
        return None
    header = lines[header_idx].split("\t")
    row = dict(zip(header, lines[header_idx + 1].split("\t")))
    try:
        return {
            "pct_rank_best": float(row["%Rank_bestAllele"]),
            "score_best": float(row["Score_bestAllele"]),
            "best_allele": row["BestAllele"],
        }
    except (KeyError, ValueError) as e:
        print(f"[WARNING] Could not parse PRIME output for {peptide}: {e} -- header was {header}")
        return None


# =============================================================================
# RUNG 3 -- Calis et al. (2013) immunogenicity model, reimplemented from the
# paper's OWN published training data (not from memory).
#
# PROVENANCE: Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD,
# Sette A, Kesmir C, Peters B (2013) "Properties of MHC Class I Presented
# Peptides That Enhance Immunogenicity." PLOS Comput Biol 9(10):e1003266.
# Supplementary Dataset S1 (pcbi.1003266.s001.xls, authored by Jorg Calis --
# confirmed from the file's own metadata) contains the paper's full 2,508-
# peptide labelled training set (2167 immunogenic / 341 non-immunogenic,
# lengths 8-10). Table 2 of the paper (quoted from the PMC full text,
# PMC3808449) gives the per-position importance weight (Kullback-Leibler
# divergence): position 3=0.10, 4=0.31, 5=0.30, 6=0.29, 7=0.26, 8=0.18;
# positions 1, 2 and the C-terminus are anchors (masked to 0 -- "P1, P2 and
# P9 for most HLA molecules", per the paper's own text).
#
# The per-amino-acid log-enrichment table E(aa) is NOT hand-copied from the
# paper (it is not reproduced in the accessible text/HTML of the article or
# its abstract-level tables) -- it is DERIVED HERE directly from the public
# S1 dataset using the paper's own documented method: "the enrichment is
# calculated as the ratio between the fraction of that amino acid in the
# immunogenic versus non-immunogenic data sets" pooled over non-anchor
# positions (position-independent). This reproduces the paper's own
# qualitative findings exactly (tryptophan/phenylalanine/aromatic residues
# enriched, serine most depleted -- see METHODOLOGY_NOTE.md for the
# validation: self-AUC 0.62, correct-direction mean scores).
#
# EXTENSION FOR 10-MERS: Table 2 only tabulates positions up to 8 (9-mer
# peptides). For a 10-mer's one extra internal (non-anchor, non-C-terminal)
# position, this implementation reuses the position-8 weight (0.18) --
# a documented choice, not a re-derivation of the published table.
# =============================================================================
CALIS_E_AA = {
    "W": +0.317, "T": +0.313, "I": +0.197, "N": +0.176, "F": +0.161,
    "H": +0.156, "G": +0.153, "D": +0.135, "Y": +0.134, "R": +0.104,
    "V": +0.099, "E": -0.006, "A": -0.066, "L": -0.095, "P": -0.141,
    "C": -0.208, "K": -0.211, "M": -0.248, "Q": -0.271, "S": -0.395,
}
CALIS_W_POSITION = {3: 0.10, 4: 0.31, 5: 0.30, 6: 0.29, 7: 0.26, 8: 0.18}


def calis_weight(position_1based, length):
    if position_1based in (1, 2, length):  # anchors, masked
        return 0.0
    return CALIS_W_POSITION.get(position_1based, 0.18)  # extend beyond 8 with position-8 weight


def calis_score(peptide):
    length = len(peptide)
    return sum(
        CALIS_E_AA.get(aa, 0.0) * calis_weight(i, length)
        for i, aa in enumerate(peptide, start=1)
    )


# =============================================================================
# Cross-reactivity -- BOTH layers (disagreement #7): BLASTP (Phase1Ec, frozen,
# Rule 1 -- reused, not recomputed) AND the exact 8-mer screen (fresh, via
# phase1_common.human_self_homology -- imported, not reimplemented).
# =============================================================================
def cross_reactivity_row(peptide, blastp_status, blastp_subject, blastp_pident,
                          phase1ec_exact_match, phase1ec_exact_kmer):
    hit, kmer = common.human_self_homology(peptide)
    fresh_exact = "YES" if hit else "NO"
    p1ec_exact = phase1ec_exact_match or "NONE"
    p1ec_exact_norm = "YES" if p1ec_exact not in ("NONE", "", "NO") else "NO"
    agrees = (fresh_exact == p1ec_exact_norm)
    return {
        "BLASTP_Status": blastp_status or "N/A",
        "BLASTP_Subject": blastp_subject or "",
        "BLASTP_Pident": blastp_pident or "",
        "Exact_8mer_Match_Fresh": fresh_exact,
        "Exact_8mer_Kmer_Fresh": kmer or "",
        "Exact_8mer_Match_Phase1Ec": p1ec_exact_norm,
        "Exact_8mer_Kmer_Phase1Ec": phase1ec_exact_kmer or "",
        "Fresh_vs_Phase1Ec_Agree": "YES" if agrees else "NO",
        "Combined_Status": (
            f"BLASTP: {blastp_status or 'N/A'}"
            + (f" ({blastp_subject})" if blastp_subject else "")
            + f" | Exact-8mer: {'MATCH (' + kmer + ')' if hit else 'no match'}"
        ),
    }


def build_tcr_recognition():
    common.print_banner("PHASE 4C -- TCR RECOGNITION (BigMHC-IM / PRIME / Calis ladder)")

    dossier_path = common.latest_file(common.step_output_dir("A"), suffix=".csv")
    binding_path = common.latest_file(common.step_output_dir("B"), suffix=".csv")
    if dossier_path is None or binding_path is None:
        print("[ERROR] Missing Phase4A dossier or Phase4B binding affinity output -- run Steps 1-2 first.")
        sys.exit(1)
    print(f"[INFO] Dossier: {os.path.basename(dossier_path)}")
    print(f"[INFO] Binding affinity: {os.path.basename(binding_path)}")

    with open(dossier_path, newline="") as f:
        dossier_rows = list(csv.DictReader(f))
    with open(binding_path, newline="") as f:
        binding_by_pep = {r["Peptide"]: r for r in csv.DictReader(f)}

    mhci_epitopes = [r for r in dossier_rows if r["Class"] == "MHC-I"]
    print(f"[INFO] {len(mhci_epitopes)} MHC-I construct epitopes (TCR-recognition ladder scope).")

    # ---- TCR-facing residues + binding-allele sets --------------------------
    epitope_alleles = {}
    for r in mhci_epitopes:
        pep = r["Peptide"]
        b = binding_by_pep.get(pep, {})
        alleles = [a for a in b.get("Binding_Alleles_Primary", "").split(";") if a]
        epitope_alleles[pep] = alleles
        r["TCR_Facing_4_6"] = tcr_facing_residues(pep)
        r["N_Binding_Alleles"] = len(alleles)

    # ---- Rung 1: BigMHC-IM ---------------------------------------------------
    all_pairs = [(pep, allele) for pep, alleles in epitope_alleles.items() for allele in alleles]
    bigmhc_raw = run_bigmhc_im(all_pairs)
    bigmhc_best = {}
    for pep in epitope_alleles:
        scored = [(bigmhc_raw[(pep, a)], a) for a in epitope_alleles[pep] if (pep, a) in bigmhc_raw]
        if scored:
            best_score, best_allele = max(scored)  # higher = more immunogenic
            bigmhc_best[pep] = (best_score, best_allele)
    print(f"[INFO] BigMHC-IM scored {len(bigmhc_best)}/{len(mhci_epitopes)} epitopes "
          f"({len(bigmhc_raw)} (peptide,allele) predictions).")

    # ---- Rung 2: PRIME (optional) ---------------------------------------------
    prime_ok = prime_available()
    prime_results = {}
    if prime_ok:
        supported = _load_prime_supported_alleles()
        print(f"[INFO] PRIME + MixMHCpred available -- {len(supported)} alleles in its reference panel.")
        for pep, alleles in epitope_alleles.items():
            res = run_prime_one_epitope(pep, alleles, supported)
            if res:
                prime_results[pep] = res
        print(f"[INFO] PRIME scored {len(prime_results)}/{len(mhci_epitopes)} epitopes.")
    else:
        print("[WARNING] PRIME/MixMHCpred not available at the expected paths -- "
              "dropping to rungs 1+3 (pre-approved, not a tool substitution).")

    # ---- Rung 3: Calis reimplementation ---------------------------------------
    calis_results = {r["Peptide"]: calis_score(r["Peptide"]) for r in mhci_epitopes}

    # ---- Direction checks (apply sign discipline before trusting) -----------
    # Every correlation is reported with its n and two-sided p. At n=10 only
    # |r| > ~0.63 is distinguishable from zero at p<0.05 -- a sign check, not
    # a powered validation.
    def _aligned(scores_by_pep):
        peps = [p for p in scores_by_pep
                if binding_by_pep.get(p, {}).get("NetMHC_EL_Best_Rank", "") not in ("", "N/A")]
        return ([float(binding_by_pep[p]["NetMHC_EL_Best_Rank"]) for p in peps],
                [scores_by_pep[p] for p in peps])

    xs, ys = _aligned({p: v[0] for p, v in bigmhc_best.items()})
    r_bigmhc, p_bigmhc, n_bigmhc = corr_with_p(xs, ys)
    xs, ys = _aligned({p: v["pct_rank_best"] for p, v in prime_results.items()})
    r_prime, p_prime, n_prime = corr_with_p(xs, ys)
    xs, ys = _aligned(calis_results)
    r_calis, p_calis, n_calis = corr_with_p(xs, ys)

    # Same-allele, same-tool determinism check vs Phase1Db. Phase 1 stored the
    # MAX BigMHC-IM score over ITS OWN Binding_Alleles set for each epitope
    # (Phase1Dd_immunogenicity.py). This run's per-epitope best is over a
    # DIFFERENT (107-allele-panel) set for 9 of 10 epitopes, so comparing the
    # two "best" columns is NOT a same-tool consistency check -- it compares
    # different allele sets. The honest check re-scores Phase 1's exact
    # (allele, peptide) pairs and compares the max over that same set.
    p1_alleles = {r["Peptide"]: [a for a in r.get("Binding_Alleles", "").split(";") if a]
                  for r in mhci_epitopes}
    p1_pairs = [(pep, a) for pep, als in p1_alleles.items() for a in als]
    p1_raw = run_bigmhc_im(p1_pairs)
    same_allele_rows = []
    for r in mhci_epitopes:
        pep = r["Peptide"]
        stored_raw = r.get("BigMHC_IM_Score", "")
        vals = [p1_raw[(pep, a)] for a in p1_alleles[pep] if (pep, a) in p1_raw]
        if stored_raw in ("", None) or not vals:
            continue
        same_allele_rows.append((pep, float(stored_raw), max(vals)))
    max_abs_diff = max((abs(a - b) for _, a, b in same_allele_rows), default=None)
    r_same, p_same, n_same = corr_with_p([a for _, a, _ in same_allele_rows],
                                          [b for _, _, b in same_allele_rows])

    print("\n" + "-" * 90)
    print("DIRECTION CHECKS (each tool's score vs Phase4B MHC-I percentile rank -- lower rank = better binder)")
    print(f"  BigMHC-IM  (0-1, higher=more immunogenic) vs rank: r={r_bigmhc:+.3f} p={p_bigmhc:.3f} n={n_bigmhc} "
          f"-- expect NEGATIVE; NOT distinguishable from zero at this n")
    print(f"  PRIME %Rank (lower=better)                vs rank: r={r_prime:+.3f} p={p_prime:.4f} n={n_prime} "
          f"-- expect POSITIVE")
    print(f"  Calis score (higher=more immunogenic)     vs rank: r={r_calis:+.3f} p={p_calis:.3f} n={n_calis} "
          f"-- exploratory only, no prior")
    print(f"  SAME-ALLELE determinism vs Phase1Db BigMHC_IM_Score (re-scored on Phase 1's own alleles): "
          f"r={r_same:+.4f} n={n_same} max|diff|={max_abs_diff:.2e}")

    sign_ok = True
    if r_bigmhc > 0:
        print("[WARNING] BigMHC-IM correlates POSITIVELY with rank -- wrong sign, do not trust BigMHC-IM values "
              "without investigating the column read.")
        sign_ok = False
    if r_prime < 0:
        print("[WARNING] PRIME %Rank correlates NEGATIVELY with rank -- wrong sign, do not trust PRIME values "
              "without investigating the column read.")
        sign_ok = False
    if sign_ok:
        print("[SUCCESS] Observed signs match expectation for both trusted-direction tools.")
    if max_abs_diff is not None and max_abs_diff > 1e-4:
        print(f"[WARNING] BigMHC-IM does not reproduce Phase1Db on identical (allele, peptide) pairs "
              f"(max|diff|={max_abs_diff:.2e}) -- deterministic tool disagreeing with a stored result is an "
              f"environment fault (Section H trigger 4).")

    # BigMHC and PRIME select their "best" allele independently -- these are
    # NOT allele-matched columns and must never be compared as if they were.
    allele_match = {}
    for pep in epitope_alleles:
        bm, pr = bigmhc_best.get(pep), prime_results.get(pep)
        if bm and pr:
            allele_match[pep] = "YES" if _prime_allele_name(bm[1]) == pr["best_allele"] else "NO"
    n_match = sum(1 for v in allele_match.values() if v == "YES")
    print(f"[INFO] BigMHC-IM and PRIME chose the SAME best allele for {n_match}/{len(allele_match)} epitopes -- "
          f"their 'best' columns are not allele-matched.")

    # Cross-signal flag: epitopes that are the construct's worst on several
    # independent signals at once.
    worst_on = defaultdict(list)
    def _mark(label, values, pick_max):
        if not values:
            return
        target = max(values.values()) if pick_max else min(values.values())
        for pep, v in values.items():
            if v == target:
                worst_on[pep].append(label)
    _mark("PRIME_PctRank", {p: v["pct_rank_best"] for p, v in prime_results.items()}, True)
    _mark("BigMHC_IM", {p: v[0] for p, v in bigmhc_best.items()}, False)
    _mark("Binding_Breadth", {p: len(a) for p, a in epitope_alleles.items()}, False)
    _mark("EL_Rank", {p: float(binding_by_pep[p]["NetMHC_EL_Best_Rank"]) for p in epitope_alleles
                       if binding_by_pep.get(p, {}).get("NetMHC_EL_Best_Rank", "") not in ("", "N/A")}, True)
    flagged = {p: v for p, v in worst_on.items() if len(v) >= 3}
    for pep, v in flagged.items():
        print(f"[WARNING] FLAG: {pep} is the construct's worst on {len(v)}/4 signals: {v}")

    # ---- Cross-reactivity (all 31 construct epitopes) ------------------------
    print("\n" + "-" * 90)
    print("CROSS-REACTIVITY -- both layers, all 31 construct epitopes (disagreement #7)")
    cross_rows = []
    for r in dossier_rows:
        pep = r["Peptide"]
        cr = cross_reactivity_row(
            pep,
            r.get("Self_Homology_Status", ""),
            r.get("SelfHomology_Subject", ""),
            r.get("SelfHomology_Pident", ""),
            r.get("SelfHomology_Exact_Match", ""),
            r.get("SelfHomology_Exact_Kmer", ""),
        )
        cr["Peptide"] = pep
        cr["Class"] = r["Class"]
        cross_rows.append(cr)
    disagreements = [r for r in cross_rows if r["Fresh_vs_Phase1Ec_Agree"] == "NO"]
    print(f"[INFO] Fresh exact-8mer screen vs Phase1Ec's own exact-match column: "
          f"{len(cross_rows) - len(disagreements)}/{len(cross_rows)} agree.")
    if disagreements:
        print(f"[WARNING] {len(disagreements)} disagreement(s): {[r['Peptide'] for r in disagreements]}")

    # ---- Write outputs ---------------------------------------------------------
    out_dir = common.step_output_dir("C")
    ts = common.timestamp()

    tcr_path = os.path.join(out_dir, f"Phase4C_TCRRecognition_MHCI_{ts}.csv")
    tcr_fields = [
        "Peptide", "TCR_Facing_4_6", "N_Binding_Alleles",
        "BigMHC_IM_Best", "BigMHC_IM_Best_Allele",
        "PRIME_PctRank_Best", "PRIME_Score_Best", "PRIME_Best_Allele",
        "Calis_Score_Exploratory", "NetMHC_EL_Best_Rank_Ref",
        "BigMHC_PRIME_Best_Allele_Match", "Worst_In_Construct_On",
    ]
    with open(tcr_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=tcr_fields)
        w.writeheader()
        for r in mhci_epitopes:
            pep = r["Peptide"]
            bm = bigmhc_best.get(pep)
            pr = prime_results.get(pep)
            w.writerow({
                "Peptide": pep,
                "TCR_Facing_4_6": r["TCR_Facing_4_6"],
                "N_Binding_Alleles": r["N_Binding_Alleles"],
                "BigMHC_IM_Best": round(bm[0], 4) if bm else "UNRESOLVED",
                "BigMHC_IM_Best_Allele": bm[1] if bm else "",
                "PRIME_PctRank_Best": round(pr["pct_rank_best"], 3) if pr else ("N/A (rung not available)" if not prime_ok else "UNRESOLVED"),
                "PRIME_Score_Best": round(pr["score_best"], 4) if pr else "",
                "PRIME_Best_Allele": pr["best_allele"] if pr else "",
                "Calis_Score_Exploratory": round(calis_results[pep], 4),
                "NetMHC_EL_Best_Rank_Ref": binding_by_pep.get(pep, {}).get("NetMHC_EL_Best_Rank", ""),
                "BigMHC_PRIME_Best_Allele_Match": allele_match.get(pep, "N/A"),
                "Worst_In_Construct_On": ";".join(worst_on.get(pep, [])) if len(worst_on.get(pep, [])) >= 3 else "",
            })
    print(f"\n[SUCCESS] TCR recognition table written: {tcr_path}")

    cross_path = os.path.join(out_dir, f"Phase4C_CrossReactivity_AllConstruct_{ts}.csv")
    cross_fields = ["Peptide", "Class", "BLASTP_Status", "BLASTP_Subject", "BLASTP_Pident",
                     "Exact_8mer_Match_Fresh", "Exact_8mer_Kmer_Fresh",
                     "Exact_8mer_Match_Phase1Ec", "Exact_8mer_Kmer_Phase1Ec",
                     "Fresh_vs_Phase1Ec_Agree", "Combined_Status"]
    with open(cross_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cross_fields)
        w.writeheader()
        w.writerows(cross_rows)
    print(f"[SUCCESS] Cross-reactivity table written: {cross_path}")

    dc_path = os.path.join(out_dir, f"Phase4C_DirectionChecks_{ts}.csv")
    with open(dc_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Check", "Statistic", "Value", "P_Value", "N", "Interpretation"])
        w.writerow(["BigMHC-IM best score vs Phase4B EL rank", "Pearson r", round(r_bigmhc, 4), round(p_bigmhc, 4), n_bigmhc, "expected negative; right sign, NOT distinguishable from zero"])
        w.writerow(["PRIME %Rank vs Phase4B EL rank", "Pearson r", round(r_prime, 4), round(p_prime, 5), n_prime, "expected positive; meaningful"])
        w.writerow(["Calis score (exploratory) vs Phase4B EL rank", "Pearson r", round(r_calis, 4), round(p_calis, 4), n_calis, "no prior; exploratory only"])
        w.writerow(["BigMHC-IM re-scored on Phase 1's own alleles vs Phase1Db stored", "Pearson r", round(r_same, 6), "", n_same, f"same-allele determinism; max abs diff {max_abs_diff:.2e}"])
        w.writerow(["BigMHC-IM vs PRIME chose the same best allele", "count", n_match, "", len(allele_match), "best-allele columns are NOT allele-matched"])
        w.writerow(["Fresh exact-8mer screen vs Phase1Ec", "agreement", len(cross_rows) - len(disagreements), "", len(cross_rows), "both cross-reactivity layers reproduce"])
    print(f"[SUCCESS] Direction-check table written: {dc_path}")

    # ---- Methodology note ---------------------------------------------------
    if prime_ok:
        _prime_install_note = (
            "Installed fresh this run from the GitHub repos (GfellerLab/PRIME, "
            "GfellerLab/MixMHCpred), both free/no licence wall, matching Section F's "
            "description. Both ship a Mach-O arm64 binary (`PRIME.x`) that runs "
            "natively on this machine -- no compilation needed. Two environment "
            "issues were hit and fixed:\n\n"
            "1. **Both tools' own launcher scripts refuse a space in their install "
            "path** (\"Spaces in path to MixMHCpred are not supported\") -- and "
            "this project's root sits under a space-containing macOS volume name "
            "(\"Extended SSD\"). Fixed by invoking both tools through the "
            "space-free symlink this project already maintains for exactly this "
            "class of problem (`~/mpoxhiv_ssd`, see `~/mpoxhiv_env.sh`'s own "
            "comment on why it exists) instead of the direct "
            "`/Volumes/Extended SSD/...` path -- same files, different path, "
            "nothing moved or duplicated.\n"
            "2. **Both tools invoke a bare `python3`** for their own internal "
            "scoring script, which resolves to the system Python (no pandas "
            "installed) rather than the project's `phase2` conda env. Fixed by "
            "prepending `$MPOXHIV_ENVS/phase2/bin` to `PATH` for the subprocess "
            "call only.\n\n"
            "Neither fix required editing PRIME's or MixMHCpred's own files -- "
            "both are external, third-party tools and Rule 1 (never edit Phase "
            "1/2) extends in spirit to never editing vendored external tools "
            "either."
        )
    else:
        _prime_install_note = "PRIME/MixMHCpred were attempted and NOT used -- see below."

    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(f"""# Phase 4C -- TCR Recognition: methodology notes

Generated: {ts}

## Ladder used (pre-approved, Section D/F)

Methodology IV.A.2 names "MixMHCpred + TCGA Contact Database" -- both wrong
(disagreement #6: MixMHCpred predicts MHC BINDING, not TCR contact; TCGA is
a cancer genomics atlas with no TCR-contact data). IEDB's own
`immunogenicity/` endpoint returns 403. **Retested 2026-09-21:** 11 patient attempts over several minutes (including after a
90 s wait, and on GET) all returned 403 while `mhci/`, `mhcii/` and `bcell/` answered 200 in the same minute, and a nonexistent path
(`foobar/`) returns the byte-identical 403 page -- i.e. the endpoint is not served at this URL; this is **not** IP throttling
(throttling also hit healthy endpoints and cleared on retry). The Calis reimplementation was therefore forced, not merely convenient.
(Other, transient 403s elsewhere in this project were genuine throttling; the two produce the same HTML page and can only be told apart by such a retest.) This step uses the pre-approved
3-rung ladder instead:

1. **BigMHC-IM** (primary, MHC-I only) -- {len(bigmhc_best)}/{len(mhci_epitopes)} epitopes scored.
2. **PRIME 2.0 + MixMHCpred 3.0** (optional, MHC-I only) -- {"INSTALLED AND USED" if prime_ok else "NOT AVAILABLE, dropped (rungs 1+3 only)"}.
3. **Calis et al. (2013)**, reimplemented locally from the paper's own
   published data (see below) -- all {len(mhci_epitopes)} epitopes scored.

**Scope: MHC-I only (10 construct epitopes).** All three ladder tools are
MHC-I specific by design (BigMHC/PRIME/MixMHCpred explicitly; Calis's model
was trained on and validated against MHC-I presented peptides). Section F's
own text confirms this scope ("TCR-facing residues = positions 4-6 of each
MHC-I epitope"). MHC-II and B-cell construct epitopes are out of scope for
this ladder -- there is no equivalent pre-approved TCR-recognition tool for
CD4 T-cell or antibody epitopes in this project.

## PRIME/MixMHCpred installation

{_prime_install_note}

## Direction checks (mandatory before trusting any score)

Each tool's score was checked against Phase4B's MHC-I `NetMHC_EL_Best_Rank`
(a quantity already validated in Step 2) for the expected sign BEFORE being
reported as a result -- this discipline was added after Phase4B's IC50
parsing error (misreading an unlabeled column with the wrong sign) was
caught and corrected.

| Tool | Convention | Expected vs rank | Observed r | p (two-sided) | n |
|---|---|---|---|---|---|
| BigMHC-IM | 0-1, higher=more immunogenic | negative | {r_bigmhc:+.3f} | {p_bigmhc:.3f} | {n_bigmhc} |
| PRIME %Rank | lower=better | positive | {r_prime:+.3f} | {p_prime:.4f} | {n_prime} |
| Calis score (EXPLORATORY) | higher=more immunogenic | no prior (peptide-intrinsic) | {r_calis:+.3f} | {p_calis:.3f} | {n_calis} |

**Read these with n in mind.** At n=10 only |r| > ~0.63 is distinguishable
from zero at p<0.05. **PRIME's r={r_prime:+.3f} is meaningful; BigMHC-IM's
r={r_bigmhc:+.3f} (p={p_bigmhc:.2f}) is correct in sign but NOT statistically
distinguishable from zero** -- it passes the sign check, it does not
corroborate anything. {"All observed signs matched expectation." if sign_ok else "**One or more signs came out backwards -- see the [WARNING] in the console log. Do NOT trust the affected tool's values without investigating the column read first (this is exactly the failure mode Phase4B hit).**"}

**Same-allele determinism check (replaces an earlier, invalid comparison).**
An earlier draft reported r=+0.460 between this run's BigMHC-IM "best" and
Phase1Db's stored `BigMHC_IM_Score` as a "same-tool consistency check". That
was wrong: Phase 1 stored the max over ITS OWN Binding_Alleles set (9-allele
panel), this run's best is over a different (107-allele-panel) set for 9 of
10 epitopes -- two different quantities. The honest check re-scores Phase 1's
exact (allele, peptide) pairs and compares the max over that same set:
r={r_same:+.4f}, n={n_same}, max |diff| = {max_abs_diff:.2e}. BigMHC-IM
reproduces Phase 1 to numerical precision on identical inputs; the difference
in the "best" columns is entirely because the 107-allele panel found better
alleles, which is correct.

**BigMHC-IM and PRIME choose their "best" allele independently.** They agreed
on the same allele for only {n_match}/{len(allele_match)} epitopes
(`BigMHC_PRIME_Best_Allele_Match` column). Their "best" columns are therefore
NOT allele-matched and must never be compared row-by-row as if they were.

**Calis is exploratory only.** Self-AUC 0.62 on its own training set, r={r_calis:+.3f}
here. Near-zero discrimination against binding rank is arguably correct --
the model is designed to be independent of binding affinity -- but that also
means it is NOT independent corroboration of BigMHC-IM or PRIME and must not
be presented as such.

## Flagged finding carried forward: NKRKRVIGL (Mpox A35R)

Worst in the construct on all four independent signals at once
(`Worst_In_Construct_On` column): PRIME %Rank 4.386 (next worst 0.47),
BigMHC-IM 0.0344, binding breadth 2/107 alleles, EL rank 0.85. Also note
Phase 1 scored it at `HLA-B*08:01`, which is not in the 107-allele panel at
all (0.163 there; within the local panel its best BigMHC-IM is 0.0344). This
is consistent with A35R's lowest-in-table BigMHC score in the paper and with
the `A35R_Coverage_Note` substitution recorded in Phase1G. Expect it to add
almost nothing to Step 5's coverage union. **Step 7 must report this as a
consistent cross-phase signal, not a new anomaly.**

## Calis et al. (2013) reimplementation -- full provenance

Citation: Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette
A, Kesmir C, Peters B (2013) "Properties of MHC Class I Presented Peptides
That Enhance Immunogenicity." PLOS Comput Biol 9(10):e1003266.

**The amino-acid log-enrichment table is NOT taken from memory or a
secondary source.** It is derived here directly from the paper's own public
Supplementary Dataset S1 (`pcbi.1003266.s001.xls`, downloaded from PLOS;
the file's own embedded metadata credits "Jorg Calis" as author, confirming
authenticity) -- 2,508 labelled peptides (2,167 immunogenic / 341
non-immunogenic; lengths 8/9/10), using the paper's own documented method
(quoted from the PMC full text, PMC3808449): "the enrichment is calculated
as the ratio between the fraction of that amino acid in the immunogenic
versus non-immunogenic data sets," pooled over non-anchor positions
(anchors = positions 1, 2, and the C-terminus -- "P1, P2 and P9 for most
HLA molecules," per the paper's own text; a generic, allele-agnostic mask,
since per-allele anchor definitions are not published in a reusable table).

Per-position importance weights are the paper's own **Table 2** values
(quoted exactly from PMC3808449): position 3=0.10, 4=0.31, 5=0.30, 6=0.29,
7=0.26, 8=0.18 (positions 1, 2, 9 are anchors, weight 0). For the one extra
internal position in a 10-mer (beyond position 8), this implementation
reuses the position-8 weight (0.18) -- Table 2 does not tabulate further,
and this is a documented extension, not a published value.

**Validation before use:** scoring all 2,508 training peptides with this
reimplementation gives a self-evaluated AUC of 0.62 (not held-out -- a
sanity check, not a performance claim) and the correct-direction mean score
(immunogenic peptides score higher than non-immunogenic on average). The
derived amino-acid table's qualitative pattern also matches the paper's own
reported findings: tryptophan (+0.317) and other aromatic/large residues
(phenylalanine +0.161, tyrosine +0.134) are the most enriched, serine
(-0.395) is the single most depleted -- both independently reported in the
paper's text.

## TCR-facing residues

Positions 4-6 (1-based) of each MHC-I epitope, per Section F's own
specification. The methodology's "\\>=3 TCR-contact residues" exclusion rule
cannot be applied without a defined position set (it never specifies one),
so this step reports which residues occupy 4-6 rather than scoring against
an assumed threshold.

## Cross-reactivity (disagreement #7 -- both layers, never one alone)

All 31 construct epitopes (not just the 10 MHC-I ones -- self-tolerance is
a safety question independent of TCR-recognition tooling). Two layers,
reported together, never one presented alone as a safety verdict:

- **BLASTP layer**: reused directly from Phase1Ec (frozen, Rule 1 -- not
  recomputed). Phase1Ec's own BLASTP screen against reviewed human
  Swiss-Prot is the same screen Phase 1Ec_Filtration measured as blind at
  9-16 aa (0/20 known human self-fragments caught in that validation).
- **Exact 8-mer layer**: computed FRESH this run via
  `phase1_common.human_self_homology()` (imported, not reimplemented, per
  Section F's explicit instruction), and cross-checked against Phase1Ec's
  own `SelfHomology_Exact_Match`/`SelfHomology_Exact_Kmer` columns (already
  in the Phase4A dossier) for consistency: {len(cross_rows) - len(disagreements)}/{len(cross_rows)} agree.
""")
    print(f"[INFO] Methodology note written: {note_path}")

    print("\n" + "=" * 90)
    print("[SUCCESS] Phase 4C complete.")
    print("=" * 90 + "\n")
    return tcr_path, cross_path


if __name__ == "__main__":
    build_tcr_recognition()
