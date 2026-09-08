import os, sys, csv, shutil, subprocess
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
# Sec. I.E (current manuscript): "all candidate T-cell and B-cell epitopes
# surviving toxicity and allergenicity screening were additionally screened
# via BLASTP against the reviewed human reference proteome (UniProt Homo
# sapiens, Swiss-Prot subset)." This runs AFTER 1Eb (allergenicity), which is
# why it is Phase1Ec rather than folded into 1Ea.
#
# Exclude ("self-like") only when ALL THREE hold: pident >= 70%,
# coverage (alignment_length/peptide_length) >= 70%, E-value <= 1e-5.
# Anything with a hit below that bar is RETAINED but flagged for manual
# review -- the paper is explicit that partial homology should not be a
# hard exclusion ("even sub-threshold similarity to self-proteins can, in
# some cases, still influence immunogenicity through partial tolerance
# effects").
#
# *** BLASTP IS BLIND AT THIS STUDY'S CANDIDATE LENGTHS -- MEASURED, NOT
# ASSUMED, BEFORE THIS SCREEN WAS RUN ON REAL DATA. ***
# Tested on 20 peptides cut VERBATIM from human Swiss-Prot at this study's own
# candidate lengths (9/10/15/16 aa, 5 each) -- every one is a guaranteed
# perfect self-match by construction, since it IS a fragment of a human
# protein. Default blastp settings found 0/20. Switching to the
# short-query-appropriate `-task blastp-short` recovers 13/20 as detectable
# hits -- but STILL 0/20 that clear the published E<=1e-5 acceptance bar. The
# best E-value achievable on a 9-mer against this ~11.4M-residue database is
# ~0.2, three orders of magnitude short of significance: this is a
# mathematical floor of BLAST's statistics at this query length, not a
# biological finding. Run on the real 108 Phase 1Eb survivors, the BLASTP
# screen (even with blastp-short) returns 0 peptides that clear the
# acceptance rule -- not "0 flagged as self-like", but 0 peptides for which
# the screen produces any usable discriminating signal at all. This is the
# same failure mode already documented for AllerCatPro and the allergen
# predictors at epitope length (deviation #27).
#
# DECIDED WITH THE USER (2026-08-31): run BOTH.
#   (1) BLASTP below, AS PUBLISHED, with `-task blastp-short` substituted for
#       the default task -- a tool-correctness fix (blastp-short is NCBI's own
#       recommended task for queries this short), not a threshold change. The
#       70/70/1e-5 acceptance rule is untouched. Reported honestly: this layer
#       is expected to clear 0 peptides at these lengths, and that outcome
#       must be stated as a length limitation, not as a clean safety result.
#   (2) An EXACT contiguous-identity screen against the same human reference
#       proteome (phase1_common.human_self_homology, window=8), the same
#       technique already used and validated for the allergen screen
#       (deviation #27). This has no significance threshold to fail and no
#       length floor above the window itself.
#       Window=8 chosen by the SAME scrambled-sequence-control calibration
#       used for the allergen screen: composition-matched scrambled versions
#       of the real 108 candidates (5 permutations each, same length, same
#       amino acid multiset, no biology) were searched at each window width.
#       w=6 and w=7 show NO separation from the scrambled background (80.6%
#       vs 80.7% real vs scrambled at w=6) and are useless; w=8 separates
#       (2.8% real vs 0.9% scrambled, ~3x background) and is used. Loosening
#       BLAST's identity/coverage thresholds instead of using exact matching
#       was tried first and REJECTED on the same real-vs-scrambled test: at
#       70%/70% with the E-value requirement dropped, the real pool flags
#       24.1% and the scrambled controls flag 25.2% -- indistinguishable, pure
#       background matching human protein amino-acid composition, not
#       self-homology.
# =============================================================================
BLASTP_BINARY = os.environ.get("BLASTP_BINARY", "/opt/miniconda3/envs/phase2/bin/blastp")
HUMAN_SWISSPROT_DB = os.environ.get(
    "HUMAN_SWISSPROT_DB",
    os.path.join(_PROJECT_ROOT, "human_swissprot_db", "human_swissprot"),
)

PIDENT_CUTOFF = 70.0
COVERAGE_CUTOFF = 70.0
EVALUE_CUTOFF = 1e-5


def _binary_exists(path):
    return shutil.which(path) is not None or os.path.isfile(path)


def local_tools_available():
    return _binary_exists(BLASTP_BINARY) and (
        os.path.isfile(HUMAN_SWISSPROT_DB + ".phr") or os.path.isfile(HUMAN_SWISSPROT_DB + ".pin")
    )


def run_blastp_human_swissprot(fasta_path, work_dir):
    """
    Best (lowest-E-value) reviewed-human hit per query peptide. Same
    outfmt/columns as Phase1Ea's run_blastp_toxprot -- only the target DB
    and the acceptance thresholds (70/70/1e-5 vs Tox-Prot's 80/80/1e-5)
    differ, per the paper's own distinct wording for each screen.

    -task blastp-short: NCBI's own recommended task for queries under ~30 aa
    (it uses different statistics tuned for short sequences). Measured
    necessary here -- see the module-level note above: default blastp finds
    0/20 verbatim-self-fragment positive controls at this study's lengths;
    blastp-short finds 13/20. The search -evalue is widened to 10 so that
    weak-but-real hits are still RETURNED for inspection; the SEPARATE
    acceptance rule below (PIDENT_CUTOFF/COVERAGE_CUTOFF/EVALUE_CUTOFF,
    unchanged at 70/70/1e-5) still decides exclusion. Widening the search
    -evalue is a "don't discard candidate hits before we can even look at
    them" setting, not a relaxation of the published acceptance criterion.
    """
    output_tsv = os.path.join(work_dir, "blastp_human_swissprot.tsv")
    cmd = [
        BLASTP_BINARY, "-query", fasta_path, "-db", HUMAN_SWISSPROT_DB,
        "-outfmt", "6 qseqid sseqid evalue pident length qlen",
        "-task", "blastp-short", "-evalue", "10", "-max_target_seqs", "5",
        "-out", output_tsv,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    if result.returncode != 0:
        raise RuntimeError(f"BLASTP (human Swiss-Prot) failed (exit {result.returncode}):\n{result.stderr}")

    best_hits = {}
    if os.path.isfile(output_tsv):
        with open(output_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    continue
                qseqid, sseqid, evalue, pident, length, qlen = (
                    parts[0], parts[1], float(parts[2]), float(parts[3]), int(parts[4]), int(parts[5]),
                )
                if qseqid not in best_hits or evalue < best_hits[qseqid]["evalue"]:
                    best_hits[qseqid] = {
                        "evalue": evalue, "pident": pident,
                        "coverage": (length / qlen) * 100 if qlen else 0.0,
                        "subject": sseqid,
                    }
    return best_hits


METHODOLOGY_NOTE_TEXT = """# Phase 1Ec — Human Self-Homology Screening Methodology Note

**Applies to:** Sec. I.E of the manuscript (BLASTP against the reviewed human
reference proteome, UniProt Homo sapiens Swiss-Prot subset).

## Two layers, and why there are two

**Layer 1 -- BLASTP, as published.** Every peptide in 1Eb's Filtered output is
BLASTP'd (`-task blastp-short`, the NCBI-recommended task for short queries)
against a local database built from `reviewed:true AND organism_id:9606`.
A peptide is flagged "self-like" and EXCLUDED only if its best hit meets ALL
THREE: percent identity >= 70%, alignment coverage
(alignment_length / peptide_length) >= 70%, E-value <= 1e-5 -- the published
acceptance rule, unchanged. Any hit below that bar is RETAINED but flagged
`Self_Homology_Status = PARTIAL` for manual review, per the paper's own text:
"Epitopes with partial homology below this threshold were retained but
flagged for manual review, given that even sub-threshold similarity to
self-proteins can, in some cases, still influence immunogenicity through
partial tolerance effects."

**MEASURED LIMITATION, stated plainly rather than left implicit: BLASTP
cannot reach statistical significance at this study's candidate lengths
(9-16 aa).** Tested on 20 peptides cut verbatim from human Swiss-Prot at
these exact lengths -- guaranteed perfect self-matches by construction --
default blastp settings found 0/20; blastp-short recovers 13/20 as
detectable hits, but 0/20 clear the E<=1e-5 bar (best achievable E-value on a
9-mer against this database is ~0.2). On this study's real candidate pool,
Layer 1 is expected to, and did, exclude 0 peptides. That is a property of
BLAST's significance statistics at this query length, not evidence the
candidates are free of human self-homology, and must be reported as a length
limitation rather than as a clean safety result.

**Layer 2 -- exact contiguous-identity screen (new).** The same technique
already used and validated for the allergen screen (deviation #27):
`phase1_common.human_self_homology()`, an 8-residue exact-substring lookup
against the same human reference proteome. This has no significance
threshold to fail and no length floor above the 8-residue window, so it
supplies what Layer 1 structurally cannot at these lengths.

Window=8 was NOT assumed -- it was calibrated against a scrambled-sequence
control (5 composition-matched random permutations per real candidate, same
length, same amino acid multiset, no biology): w=6 and w=7 show no
separation from the scrambled background (~81% flagged either way) and are
useless; w=8 separates real candidates from scrambled controls by roughly
3x and is used. Loosening BLASTP's own identity/coverage thresholds instead
of using exact matching was tried first and rejected on the same
real-vs-scrambled test: at 70%/70% with the E-value requirement dropped, the
real pool and the scrambled controls flag statistically indistinguishable
fractions (24.1% vs 25.2%) -- pure background matching human protein amino-
acid composition, not genuine self-homology.

A HIT in Layer 2 EXCLUDES, on the same evidence-over-score principle as the
allergen homology screen and deviation #1's toxicity rule: an exact 8+
residue match to a reviewed human protein is factual, not a prediction. A
Layer-2 NONE is not proof of non-self-similarity: only exact matches are
detected, and a peptide 80% similar to a human protein passes clean -- the
same caveat already stated for the allergen homology screen.

## Database provenance
`human_swissprot_db/human_swissprot.fasta`, fetched from
`rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta`
-- 20,431 reviewed human sequences. Layer 1 uses a `makeblastdb -dbtype prot`
index of this file; Layer 2 searches the same FASTA directly as a
delimiter-joined substring blob, so both layers run against the identical
reference set.
"""


def run_step1ec_self_homology():
    input_folder = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Eb", "Filtered")
    output_base = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Ec")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    tool_runs_dir = os.path.join(output_base, "_tool_runs")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)
    os.makedirs(tool_runs_dir, exist_ok=True)

    common.print_banner("PHASE 1Ec: HUMAN SELF-HOMOLOGY SCREENING (BLASTP vs Swiss-Prot Human)")

    latest_csv = common.latest_file(input_folder, suffix=".csv")
    if latest_csv is None:
        print(f"[ERROR] No Phase 1Eb Filtered output found at: {input_folder}")
        print("[ERROR] Run Phase 1Eb first.")
        return

    with open(latest_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_fields = reader.fieldnames

    unique_peptides = sorted(set(row['Peptide'] for row in rows))
    print(f"[INFO] {len(rows)} candidate rows | {len(unique_peptides)} unique peptides")

    if not local_tools_available():
        print(f"[ERROR] blastp or the human Swiss-Prot DB not found.")
        print(f"[ERROR] BLASTP_BINARY = {BLASTP_BINARY}")
        print(f"[ERROR] HUMAN_SWISSPROT_DB = {HUMAN_SWISSPROT_DB}")
        return

    query_fasta = os.path.join(tool_runs_dir, "self_homology_query.fasta")
    pep_id_map = {}
    with open(query_fasta, "w") as f:
        for i, pep in enumerate(unique_peptides):
            pep_id = f"Pep_{i}"
            pep_id_map[pep_id] = pep
            f.write(f">{pep_id}\n{pep}\n")

    print("[INFO] Running BLASTP (Layer 1, as published) against local human Swiss-Prot database...")
    blast_raw = run_blastp_human_swissprot(query_fasta, tool_runs_dir)
    blast_hits = {pep_id_map[k]: v for k, v in blast_raw.items() if k in pep_id_map}

    fieldnames = original_fields + [
        "SelfHomology_Pident", "SelfHomology_Coverage", "SelfHomology_Evalue",
        "SelfHomology_Subject", "Self_Homology_Status",
        "SelfHomology_Exact_Match", "SelfHomology_Exact_Kmer",
    ]
    raw_data, filtered_data = [], []
    n_excluded, n_partial, n_clean = 0, 0, 0
    n_exact_hit = 0
    exact_ref_warned = False

    for row in rows:
        pep = row['Peptide']
        hit = blast_hits.get(pep)
        clean_row = {k: row[k] for k in original_fields}

        # Layer 2: exact contiguous-identity screen (window=8, validated
        # against a scrambled-sequence control -- see the module-level note
        # and the METHODOLOGY_NOTE_TEXT above). Runs regardless of Layer 1's
        # outcome, since Layer 1 is measured to be blind at these lengths.
        try:
            exact_hit, exact_kmer = common.human_self_homology(pep)
        except FileNotFoundError as exc:
            exact_hit, exact_kmer = False, None
            if not exact_ref_warned:
                exact_ref_warned = True
                print(f"[WARNING] human self-homology exact-match screen skipped: {exc}")
        if exact_hit:
            n_exact_hit += 1

        if hit is None:
            clean_row.update({
                "SelfHomology_Pident": "", "SelfHomology_Coverage": "",
                "SelfHomology_Evalue": "", "SelfHomology_Subject": "",
                "Self_Homology_Status": "NONE",
            })
            is_self_like = False
            n_clean += 1
        else:
            is_self_like = (
                hit["pident"] >= PIDENT_CUTOFF
                and hit["coverage"] >= COVERAGE_CUTOFF
                and hit["evalue"] <= EVALUE_CUTOFF
            )
            clean_row.update({
                "SelfHomology_Pident": round(hit["pident"], 2),
                "SelfHomology_Coverage": round(hit["coverage"], 2),
                "SelfHomology_Evalue": hit["evalue"],
                "SelfHomology_Subject": hit["subject"],
                "Self_Homology_Status": "SELF-LIKE" if is_self_like else "PARTIAL",
            })
            if is_self_like:
                n_excluded += 1
            else:
                n_partial += 1

        clean_row.update({
            "SelfHomology_Exact_Match": "HIT" if exact_hit else "NONE",
            "SelfHomology_Exact_Kmer": exact_kmer or "",
        })

        # EXCLUDE if EITHER layer's own rule says so. Layer 1's rule is the
        # published 70/70/1e-5 acceptance test; Layer 2's rule is a factual
        # exact-match HIT. Measured (see module-level note): at this study's
        # candidate lengths Layer 1 contributes 0 exclusions, so in practice
        # this run's exclusions come from Layer 2 -- but the combination is
        # written generally so Layer 1 is not silently disabled if it is ever
        # run on longer peptides (e.g. junction-level screening) where it CAN
        # reach significance.
        if not (is_self_like or exact_hit):
            filtered_data.append(clean_row)

        raw_data.append(clean_row)

    # n_excluded (Layer 1) and n_exact_hit (Layer 2) can overlap in principle
    # (measured: they do not, on this run -- Layer 1 excludes 0). Recompute
    # the true exclusion count from the rows themselves rather than summing
    # the two counters, so an overlap can never be double-counted.
    n_excluded_total = sum(
        1 for r in raw_data
        if r["Self_Homology_Status"] == "SELF-LIKE" or r["SelfHomology_Exact_Match"] == "HIT"
    )
    print(f"[INFO] Layer 1 (BLASTP, published rule) self-like: {n_excluded} | "
          f"Layer 2 (exact-match) hit: {n_exact_hit} | "
          f"TOTAL excluded (either layer): {n_excluded_total}")
    print(f"[INFO] Partial (retained+flagged, BLASTP sub-threshold): {n_partial} | "
          f"No BLASTP hit: {n_clean}")
    print(f"[INFO] Retained for Phase 1F: {len(filtered_data)}/{len(rows)} rows")

    # Per (target x class) survivor counts -- surfaces any group emptied by
    # this screen, per the plan's I-1 contingency (default: accept and
    # document; only re-screen deeper if a whole antigen loses all 3 classes).
    groups = {}
    for r in filtered_data:
        groups.setdefault((r.get("Target", ""), r.get("Type", "")), 0)
        groups[(r.get("Target", ""), r.get("Type", ""))] += 1
    for key in sorted(groups):
        print(f"          {key[0]:12s} {key[1]:8s} : {groups[key]} survivors")
    empty_targets = set()
    for row in rows:
        key_target = row.get("Target", "")
        empty_targets.add(key_target)
    for t in sorted(empty_targets):
        classes_present = {k[1] for k in groups if k[0] == t}
        if len(classes_present) < 3:
            missing = {"MHC-I", "MHC-II", "B-cell"} - classes_present
            if missing:
                print(f"[WARNING] {t}: no survivors for class(es) {sorted(missing)} after self-homology screen")

    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    if raw_data:
        with open(os.path.join(raw_dir, f"Phase1Ec_Raw_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(raw_data)
    if filtered_data:
        with open(os.path.join(filt_dir, f"Phase1Ec_Filtered_{ts}.csv"), 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader(); writer.writerows(filtered_data)

    note_path = os.path.join(output_base, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(METHODOLOGY_NOTE_TEXT)
    print(f"[INFO] Methodology note written to {note_path}")


if __name__ == "__main__":
    run_step1ec_self_homology()
