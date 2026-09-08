import os
import sys
import csv
import time
from datetime import datetime

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase1_common module.
# See Phase 1/_common/phase1_common.py for why this logic is centralized.
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
_COMMON_DIR = os.path.join(_PROJECT_ROOT, "Phase 1", "_common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)

import phase1_common as common

# =============================================================================
# This is the per-run SCREENING script (reads Phase1B output, writes
# Filtered_Antigenicity). The offline model TRAINING job lives in the
# separate train_esm2_antigenicity.py in this same folder -- run that once
# first; this script loads its saved output from ESM2_MODEL_DIR.
# =============================================================================
ESM2_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "models", "esm2_antigenicity_finetuned")
ESM2_THRESHOLD = 0.50
ESM2_MAX_LENGTH = 1024

# Kolaskar & Tongaonkar per-residue antigenic propensity scale (Table 1)
KT_SCALE = {
    'A': 1.8, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
    'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 3.8,
    'M': 1.126, 'N': 0.851, 'P': 1.064, 'Q': 1.010, 'R': 0.873,
    'S': 1.012, 'T': 0.909, 'V': 1.187, 'W': 1.085, 'Y': 1.255,
}
KT_WINDOW = 7
KT_THRESHOLD = 1.00
KT_MIN_PEPTIDE_LEN = 6


def kt_antigenicity_verdict(seq):
    """
    Sliding window of 7 residues, K&T propensities averaged and assigned to
    the CENTRAL residue. Positions flagged (>=1.00) are merged into
    contiguous candidate peptides; a protein is classified antigenic under
    this method if at least one such peptide is >= 6 aa long.

    Sequences with non-standard letters (X,B,Z,J,U,gaps,stops) are excluded
    per the methodology -- in practice Phase 1B already only writes
    standard-20-letter sequences to Filtered_Stability, but this is
    re-checked here defensively since it's this step's own stated rule.
    """
    if any(c not in KT_SCALE for c in seq):
        return None, []  # excluded, not scoreable

    n = len(seq)
    if n < KT_WINDOW:
        return False, []

    half = KT_WINDOW // 2
    flagged = [False] * n
    for center in range(half, n - half):
        window = seq[center - half: center + half + 1]
        mean_score = sum(KT_SCALE[c] for c in window) / KT_WINDOW
        flagged[center] = mean_score >= KT_THRESHOLD

    # merge contiguous flagged positions into candidate peptides
    peptides = []
    start = None
    for i, is_flagged in enumerate(flagged + [False]):
        if is_flagged and start is None:
            start = i
        elif not is_flagged and start is not None:
            end = i
            if end - start >= KT_MIN_PEPTIDE_LEN:
                peptides.append(seq[start:end])
            start = None

    return (len(peptides) > 0), peptides


def load_esm2_classifier():
    """Loads the fine-tuned ESM-2 model for inference. Returns None (with a
    clear message) if it hasn't been trained yet -- this script degrades to
    K&T-only screening rather than fabricating a probability."""
    if not os.path.isdir(ESM2_MODEL_DIR):
        return None
    try:
        import torch
        from transformers import AutoTokenizer, EsmForSequenceClassification
        tokenizer = AutoTokenizer.from_pretrained(ESM2_MODEL_DIR)
        model = EsmForSequenceClassification.from_pretrained(ESM2_MODEL_DIR)
        model.eval()
        device = "mps" if torch.backends.mps.is_available() else ("cuda" if torch.cuda.is_available() else "cpu")
        model.to(device)
        return {"tokenizer": tokenizer, "model": model, "device": device, "torch": torch}
    except Exception as e:
        print(f"[WARNING] Failed to load ESM-2 model from {ESM2_MODEL_DIR}: {e}")
        return None


def esm2_antigenicity_probability(classifier, seq):
    torch = classifier["torch"]
    tokenizer, model, device = classifier["tokenizer"], classifier["model"], classifier["device"]
    inputs = tokenizer(seq, return_tensors="pt", truncation=True, max_length=ESM2_MAX_LENGTH)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        logits = model(**inputs).logits
        prob = torch.softmax(logits, dim=1)[0, 1].item()
    return prob


# =============================================================================
# CONSENSUS DISAGREEMENT POLICY -- CHANGED 2026-08-30
#
# Sec. I.C says agreement between ESM-2 and Kolaskar-Tongaonkar "determined
# final classification, with disagreements flagged for exclusion". This script
# originally excluded on disagreement. Two findings changed that:
#
# (1) THE MANUSCRIPT ITSELF DESCRIBES THIS STEP AS NON-ELIMINATING. Results 1.3:
#     Section I.C "functioned as a cross-methodological consistency check to
#     confirm baseline antigenicity across all source sequences rather than a
#     selective elimination filter." PHASE_I_II_FINDINGS.txt repeats it. Under
#     the previous run's ZERO disagreements the two readings were
#     indistinguishable, so the contradiction was invisible.
#
# (2) THE DISAGREEMENTS ARE A CLASSIFIER LENGTH ARTIFACT, MEASURED DIRECTLY.
#     On the corrected antigen pool, all 27 disagreements were HIV_p17 (27 of
#     its 30 variants). Controlled test on the model:
#         p24 full  231 aa -> P(antigenic) 0.9682
#         p24[:131] 131 aa -> P            0.1125   (same protein, shorter)
#         p24[:100] 100 aa -> P            0.0319
#         p17       131 aa -> P            0.0327
#         p17 repeated x2/x3/x4, composition IDENTICAL -> 0.6141/0.7577/0.7686
#     The verdict tracks LENGTH, not antigenicity. p17 (matrix) is a
#     well-characterised antigen with curated CTL epitopes in IEDB. Excluding it
#     would delete a required antigen on a demonstrated measurement artifact --
#     the same failure mode deviation #2 documents for the instability filter.
#     This was masked before because the polyprotein bug meant p17 was scored as
#     the 1437-aa Gag-Pol record, which is long enough to pass.
#
# Set EXCLUDE_ON_DISAGREEMENT = True to restore the literal exclusion behaviour.
# Disagreements are ALWAYS recorded to Disagreements/ either way.
# =============================================================================
EXCLUDE_ON_DISAGREEMENT = False


def run_step1c_antigenicity():
    start_time = time.time()

    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1B", "Filtered_Stability")
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1C")
    filt_dir = os.path.join(output_base, "Filtered_Antigenicity")
    raw_dir = os.path.join(output_base, "Raw")
    disagree_dir = os.path.join(output_base, "Disagreements")
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(disagree_dir, exist_ok=True)

    print("\n" + "="*80 + "\nPHASE 1C: ANTIGENICITY SCREENING (Kolaskar & Tongaonkar + ESM-2)\n" + "="*80)

    if not os.path.isdir(input_folder):
        print(f"[ERROR] Phase 1B Filtered_Stability directory not found at: {input_folder}")
        return
    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    if not fasta_files:
        print("[ERROR] No stable sequences found from Phase 1B.")
        return

    classifier = load_esm2_classifier()
    if classifier is None:
        print(f"[WARNING] No fine-tuned ESM-2 model found at {ESM2_MODEL_DIR}.")
        print("[WARNING] Run train_esm2_antigenicity.py first for the full K&T/ESM-2 consensus.")
        print("[WARNING] Proceeding on Kolaskar & Tongaonkar alone -- results are NOT the paper's")
        print("[WARNING] consensus method until the ESM-2 model exists and this is re-run.")
    else:
        print(f"[INFO] Loaded fine-tuned ESM-2 model from {ESM2_MODEL_DIR} (device={classifier['device']})")

    raw_rows, disagreements = [], []
    n_kept = 0

    for i, filename in enumerate(fasta_files):
        file_path = os.path.join(input_folder, filename)
        with open(file_path, "r") as f:
            seq = "".join(line.strip() for line in f if not line.startswith(">"))

        kt_verdict, kt_peptides = kt_antigenicity_verdict(seq)

        if kt_verdict is None:
            row = {"Variant_ID": filename, "KT_Verdict": "EXCLUDED_NONSTANDARD", "ESM2_Prob": "",
                   "ESM2_Verdict": "", "Final_Verdict": "EXCLUDED", "KT_Peptide_Count": 0}
            raw_rows.append(row)
            continue

        if classifier is not None:
            esm2_prob = esm2_antigenicity_probability(classifier, seq)
            esm2_verdict = esm2_prob >= ESM2_THRESHOLD
            agree = (kt_verdict == esm2_verdict)
            final = kt_verdict if agree else False
        else:
            esm2_prob, esm2_verdict, agree = "", "", ""
            final = kt_verdict

        row = {
            "Variant_ID": filename, "KT_Verdict": kt_verdict,
            "ESM2_Prob": round(esm2_prob, 4) if esm2_prob != "" else "",
            "ESM2_Verdict": esm2_verdict, "Final_Verdict": "ANTIGENIC" if final else "NON-ANTIGENIC",
            "KT_Peptide_Count": len(kt_peptides),
        }
        raw_rows.append(row)

        if classifier is not None and not agree:
            disagreements.append(row)

        # Retention: a sequence is carried forward unless it is UNSCOREABLE
        # (non-standard residues) or, when EXCLUDE_ON_DISAGREEMENT is set, the
        # two methods disagree. Kolaskar-Tongaonkar remains the primary call.
        retained = final if EXCLUDE_ON_DISAGREEMENT else kt_verdict
        row["Retained"] = retained
        row["Consensus"] = ("AGREE" if agree else "DISAGREE") if classifier is not None else "KT_ONLY"
        if retained:
            n_kept += 1
            with open(os.path.join(filt_dir, filename), "w") as out_f:
                out_f.write(f">{filename.replace('.fasta', '')}\n{seq}\n")

        elapsed = common.format_time(time.time() - start_time)
        sys.stdout.write(f"\r[ EVAL ] {i+1:03d}/{len(fasta_files):03d} | Elapsed: {elapsed}")
        sys.stdout.flush()

    print()
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    fieldnames = ["Variant_ID", "KT_Verdict", "ESM2_Prob", "ESM2_Verdict", "Final_Verdict", "KT_Peptide_Count", "Consensus", "Retained"]
    with open(os.path.join(raw_dir, f"Phase1C_Raw_{ts}.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(raw_rows)
    if disagreements:
        with open(os.path.join(disagree_dir, f"Phase1C_Disagreements_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(disagreements)

    print("\n" + "="*80)
    print(f"[SUCCESS] Retained: {n_kept}/{len(fasta_files)}")
    if classifier is not None:
        mode = "EXCLUDED" if EXCLUDE_ON_DISAGREEMENT else "FLAGGED, retained"
        print(f"[INFO] K&T/ESM-2 disagreements ({mode}): {len(disagreements)}")
        if not EXCLUDE_ON_DISAGREEMENT and disagreements:
            import collections as _c
            per = _c.Counter(r["Variant_ID"].split("_Var")[0] for r in disagreements)
            print("[INFO] Disagreements by target: "
                  + ", ".join(f"{k} {v}" for k, v in sorted(per.items())))
            print("[INFO] This step is a CONSISTENCY CHECK, not an elimination")
            print("       filter -- see the policy note at the top of this file.")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1c_antigenicity()
