import os
import sys
import csv
import re
from collections import Counter

# =============================================================================
# MINIMAL BOOTSTRAP -- locates Phase 4's shared module. Copied verbatim from
# the Phase 1 bootstrap pattern (see Phase 1/_common/phase1_common.py).
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
# PHASE 4A -- EPITOPE DOSSIER (methodology Section IV, Step 1 of 7)
#
# No network calls. This step only JOINS existing Phase 1 outputs into one
# table, one row per construct epitope -- it computes nothing new. Every
# join is by PEPTIDE SEQUENCE, never by row number (row order is not stable
# across these files and was never intended to be).
#
# Section IV of REVISED_METHODOLOGY_MpoxHIV.md says "194 Non-Redundant
# Sequences". That is a MISLABEL: 194 is the number of source PROTEIN
# sequences retrieved in Phase 1A, not epitopes. The real analysis unit is
# the 31 epitopes in the final construct (this step) and, for comparison,
# the 105 candidates that survived to Phase 1F coverage scoring (used from
# Phase 4B onward). Logged in MANUSCRIPT_CORRECTIONS_PHASE4.md (Phase 4G).
# =============================================================================

CONSTRUCT_SOURCES = {
    "construct": ("Step_Outputs", "Phase1", "Phase1G"),
    "binding_bcell": ("Step_Outputs", "Phase1", "Phase1D", "Phase1Db"),
    "conservancy": ("Step_Outputs", "Phase1", "Phase1D", "Phase1Dc",
                     "Filtered_Benchmarks", "Min_50pct"),
    "coverage_pool": ("Step_Outputs", "Phase1", "Phase1F", "Filtered"),
    "self_homology": ("Step_Outputs", "Phase1", "Phase1E", "Phase1Ec", "Filtered"),
}

# Epitope_Provenance writes the Mpox EEV-glycoprotein target as "Mpox_B5R",
# but the manuscript (and its correct ortholog identity, see Phase 1A's
# antigen-identity check) calls it B6R (OPG190). Both names are kept on every
# row so neither reading of the manuscript is lost; logged as a correction.
ANTIGEN_DISPLAY_NAME = {
    "Mpox_B5R": "Mpox_B5R (manuscript: B6R / OPG190)",
}


def _find_input(key, suffix=".csv"):
    folder = os.path.join(_PROJECT_ROOT, *CONSTRUCT_SOURCES[key])
    path = common.latest_file(folder, suffix=suffix)
    if path is None:
        print(f"[ERROR] No {suffix} file found under: {folder}")
    return path


def _load_by_peptide(path, wanted_fields):
    """Returns {peptide: {field: value, ...}} restricted to wanted_fields."""
    out = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        missing_cols = [c for c in wanted_fields if c not in reader.fieldnames]
        if missing_cols:
            print(f"[WARNING] {os.path.basename(path)} is missing expected columns: {missing_cols}")
        for row in reader:
            pep = row["Peptide"]
            out[pep] = {c: row.get(c, "") for c in wanted_fields}
    return out


def _parse_construct_row(construct_path):
    """
    Reads Epitope_Provenance ("PEPTIDE:Target;...") and Boundary_Map
    ("...;MHC-I:PEPTIDE[start-end];...") from the Phase 1G construct CSV.
    Returns (construct_id, sequence, provenance, klass, positions) where
    provenance/klass/positions are all keyed by peptide.
    """
    with open(construct_path, newline="") as f:
        row = next(csv.DictReader(f))

    provenance = {}
    for item in row["Epitope_Provenance"].split(";"):
        item = item.strip()
        if ":" in item:
            pep, target = item.rsplit(":", 1)
            provenance[pep.strip()] = target.strip()

    klass = {}
    positions = {}
    for seg in row["Boundary_Map"].split(";"):
        seg = seg.strip()
        for cls in ("MHC-I", "MHC-II", "B-cell"):
            prefix = cls + ":"
            if seg.startswith(prefix):
                rest = seg[len(prefix):]
                m = re.match(r"^(.*)\[(\d+)-(\d+)\]$", rest)
                if not m:
                    print(f"[WARNING] Could not parse Boundary_Map segment: {seg!r}")
                    break
                pep = m.group(1)
                klass[pep] = cls
                positions[pep] = (int(m.group(2)), int(m.group(3)))
                break

    return row["Construct_ID"], row["Sequence"], provenance, klass, positions


def build_dossier():
    common.print_banner("PHASE 4A -- EPITOPE DOSSIER (join Phase 1 outputs, no new predictions)")

    construct_path = _find_input("construct")
    if construct_path is None:
        print("[ERROR] Cannot proceed without the Phase 1G final construct.")
        sys.exit(1)
    print(f"[INFO] Construct source: {os.path.basename(construct_path)}")

    construct_id, sequence, provenance, klass, positions = _parse_construct_row(construct_path)
    if construct_id != common.CONSTRUCT_ID:
        print(f"[WARNING] Construct_ID in file ({construct_id}) does not match "
              f"expected ({common.CONSTRUCT_ID}) -- proceeding with the file's own ID.")

    peptides = set(provenance.keys())
    print(f"[INFO] {len(peptides)} epitopes read from {construct_id} "
          f"({len(sequence)} aa full construct).")

    # ---- Load the four Phase 1 join sources -----------------------------
    db_path = _find_input("binding_bcell")
    dc_path = _find_input("conservancy")
    f1f_path = _find_input("coverage_pool")
    ec_path = _find_input("self_homology")
    for label, p in (("Phase1Db", db_path), ("Phase1Dc", dc_path),
                      ("Phase1F", f1f_path), ("Phase1Ec", ec_path)):
        if p is None:
            print(f"[ERROR] Cannot proceed without {label}.")
            sys.exit(1)
        print(f"[INFO] {label} source: {os.path.basename(p)}")

    db = _load_by_peptide(db_path, [
        "Target", "Variant", "Type", "Length", "GRAVY",
        "Percentile_Rank", "mean_BepiPred", "pct_above", "Bcell_Tier",
        "Binding_Alleles", "BigMHC_IM_Score", "Immunogenicity_Priority",
        "SEMA_Overlap_Pct", "SEMA_Corroborated", "SEMA_Residues_Scored",
    ])
    dc = _load_by_peptide(dc_path, ["Conservancy", "Hit_Ratio"])
    f1f = _load_by_peptide(f1f_path, [
        "Overall_Coverage_Pct", "Coverage_Status", "Coverage_Priority",
    ])
    ec = _load_by_peptide(ec_path, [
        "SelfHomology_Pident", "SelfHomology_Coverage", "SelfHomology_Evalue",
        "SelfHomology_Subject", "Self_Homology_Status",
        "SelfHomology_Exact_Match", "SelfHomology_Exact_Kmer",
    ])

    # ---- Join ------------------------------------------------------------
    rows = []
    unjoined = {"Phase1Db": [], "Phase1Dc": [], "Phase1F": [], "Phase1Ec": []}
    for pep in sorted(peptides, key=lambda p: (klass.get(p, ""), p)):
        target = provenance[pep]
        pathogen = common.pathogen_of(target)
        cls = klass.get(pep, "")
        start, end = positions.get(pep, ("", ""))

        row = {
            "Peptide": pep,
            "Class": cls,
            "Pathogen": pathogen,
            "Labelled_Target": target,
            "Antigen_Display": ANTIGEN_DISPLAY_NAME.get(target, target),
            "Construct_Position_Start": start,
            "Construct_Position_End": end,
            "Source_Construct": construct_id,
        }

        if pep in db:
            row.update(db[pep])
        else:
            unjoined["Phase1Db"].append(pep)
        if pep in dc:
            row.update(dc[pep])
        else:
            unjoined["Phase1Dc"].append(pep)
        if pep in f1f:
            row.update(f1f[pep])
        else:
            unjoined["Phase1F"].append(pep)
        if pep in ec:
            row.update(ec[pep])
        else:
            unjoined["Phase1Ec"].append(pep)

        rows.append(row)

    # ---- Cross-check: Phase1F carries its own copy of self-homology status;
    # Phase1Ec is the source of record per the task spec. They should agree
    # for every construct epitope since both derive from the same BLASTP run.
    f1f_full = _load_by_peptide(f1f_path, ["Self_Homology_Status"])
    disagreements = []
    for pep in peptides:
        if pep in f1f_full and pep in ec:
            a, b = f1f_full[pep]["Self_Homology_Status"], ec[pep]["Self_Homology_Status"]
            if a != b:
                disagreements.append((pep, a, b))
    if disagreements:
        print(f"[WARNING] {len(disagreements)} epitope(s) have disagreeing Self_Homology_Status "
              f"between Phase1F and Phase1Ec: {disagreements}")
    else:
        print("[INFO] Self_Homology_Status agrees between Phase1F and Phase1Ec for all construct epitopes.")

    # ---- Assertions (Section E Step 1 checklist) --------------------------
    failed = []

    if len(rows) != common.EXPECTED_TOTAL_EPITOPES:
        failed.append(f"row count = {len(rows)}, expected {common.EXPECTED_TOTAL_EPITOPES}")

    class_counts = Counter(r["Class"] for r in rows)
    if dict(class_counts) != common.EXPECTED_BY_CLASS:
        failed.append(f"class counts = {dict(class_counts)}, expected {common.EXPECTED_BY_CLASS}")

    pathogen_counts = Counter(r["Pathogen"] for r in rows)
    if pathogen_counts.get("HIV", 0) != common.EXPECTED_HIV_COUNT or \
       pathogen_counts.get("Mpox", 0) != common.EXPECTED_MPOX_COUNT:
        failed.append(f"pathogen counts = {dict(pathogen_counts)}, "
                       f"expected HIV={common.EXPECTED_HIV_COUNT} Mpox={common.EXPECTED_MPOX_COUNT}")

    n_antigens = len(set(r["Labelled_Target"] for r in rows))
    if n_antigens != common.EXPECTED_SOURCE_ANTIGENS:
        failed.append(f"distinct source antigens = {n_antigens}, expected {common.EXPECTED_SOURCE_ANTIGENS}")

    any_unjoined = {k: v for k, v in unjoined.items() if v}
    if any_unjoined:
        failed.append(f"epitopes failed to join: {any_unjoined}")

    print("\n" + "-" * 90)
    print("ASSERTION SUMMARY")
    print(f"  Total epitopes:        {len(rows)} (expected {common.EXPECTED_TOTAL_EPITOPES})")
    print(f"  By class:              {dict(class_counts)} (expected {common.EXPECTED_BY_CLASS})")
    print(f"  By pathogen:           {dict(pathogen_counts)} "
          f"(expected HIV={common.EXPECTED_HIV_COUNT}, Mpox={common.EXPECTED_MPOX_COUNT})")
    print(f"  Distinct source antigens: {n_antigens} (expected {common.EXPECTED_SOURCE_ANTIGENS})")
    print(f"  Join failures:          {sum(len(v) for v in unjoined.values())}")

    if failed:
        details = "\n".join(f"  - {f}" for f in failed)
        print("\n[ERROR] One or more Step 1 assertions failed:")
        print(details)
        common.halt_for_opus(
            branch="Phase4A_epitopeDossier",
            trigger_number=2,
            one_liner="Epitope dossier assertion failure (join or count mismatch)",
            what_i_was_doing="Joining the 31 Vax_Final_4fce119e construct epitopes against "
                              "Phase1Db/Phase1Dc/Phase1F/Phase1Ec by peptide sequence.",
            exact_numbers=details,
            options=[
                "Re-check the Phase 1G Epitope_Provenance / Boundary_Map parsing for a format change.",
                "Re-check whether the Phase 1 source CSVs picked up by latest_file() are the intended ones "
                "(multiple timestamped files could exist under a folder).",
                "If a genuine count mismatch, this changes what '31 epitopes / 10-11-10 / 14-17' means "
                "for the rest of Phase IV and needs Opus judgment before Steps B-G proceed.",
            ],
        )
        sys.exit(1)

    print("\n[SUCCESS] All Step 1 assertions passed.")

    # ---- Write output ------------------------------------------------------
    out_dir = os.path.join(common.step_output_dir("A"))
    ts = common.timestamp()
    out_path = os.path.join(out_dir, f"Phase4A_EpitopeDossier_{ts}.csv")
    fieldnames = [
        "Peptide", "Class", "Pathogen", "Labelled_Target", "Antigen_Display",
        "Construct_Position_Start", "Construct_Position_End", "Source_Construct",
        "Target", "Variant", "Type", "Length", "GRAVY",
        "Percentile_Rank", "mean_BepiPred", "pct_above", "Bcell_Tier",
        "Binding_Alleles", "BigMHC_IM_Score", "Immunogenicity_Priority",
        "SEMA_Overlap_Pct", "SEMA_Corroborated", "SEMA_Residues_Scored",
        "Conservancy", "Hit_Ratio",
        "Overall_Coverage_Pct", "Coverage_Status", "Coverage_Priority",
        "SelfHomology_Pident", "SelfHomology_Coverage", "SelfHomology_Evalue",
        "SelfHomology_Subject", "Self_Homology_Status",
        "SelfHomology_Exact_Match", "SelfHomology_Exact_Kmer",
    ]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\n[SUCCESS] Dossier written: {out_path}")

    # ---- Methodology note ---------------------------------------------------
    note_path = os.path.join(out_dir, "METHODOLOGY_NOTE.md")
    with open(note_path, "w") as f:
        f.write(f"""# Phase 4A -- Epitope Dossier: methodology notes

Generated: {ts}

## Deviations from a literal reading of Section IV

1. **"194 Non-Redundant Sequences" is a mislabel.** Section IV's own text calls
   194 the epitope count. 194 is actually the number of source PROTEIN
   sequences retrieved in Phase 1A (across all HIV/Mpox targets and variants),
   not the number of epitopes analysed. Phase IV's real analysis units are the
   **31 epitopes** in the final construct {construct_id} (this dossier) and,
   from Phase 4B onward, the **105 candidates** that survived to Phase 1F
   coverage scoring. See MANUSCRIPT_CORRECTIONS_PHASE4.md for the corrected text.

2. **Mpox_B5R / B6R naming.** `Epitope_Provenance` in the Phase 1G construct
   CSV labels the Mpox EEV-glycoprotein epitopes' target as `Mpox_B5R`. The
   manuscript and Phase 1A's antigen-identity check both call the correct
   ortholog **B6R** (OPG190; VACV B5R ortholog -- MPXV's own gene literally
   named "B5R" is a different protein, OPG189, an ankyrin repeat protein; see
   Phase 1A_provenanceCorrection.py's ANTIGEN_REFERENCES). Both names are kept
   on every row (`Labelled_Target` = as-recorded, `Antigen_Display` = with the
   manuscript name) so neither reading is lost.

## Join method

All joins are by **peptide sequence** (never row number) against:
- Phase1G (`Epitope_Provenance`, `Boundary_Map`, `Sequence`) -- construct membership, class, position
- Phase1Db Elite -- percentile rank, BepiPred, BigMHC-IM, SEMA overlap
- Phase1Dc Min_50pct -- conservancy
- Phase1F Elite (Filtered) -- overall population coverage
- Phase1Ec Filtered -- self-homology status (cross-checked against Phase1F's
  own copy of the same columns; any disagreement is logged as a [WARNING] above)

## Assertions enforced

31 total rows; 10 MHC-I / 11 MHC-II / 10 B-cell; 14 HIV-derived / 17
Mpox-derived; 7 distinct source antigens; zero unjoined epitopes across all
four source files. A failure on any of these halts the step and escalates
per Section H trigger #2 (join/count failure), rather than silently
proceeding with partial data.
""")
    print(f"[INFO] Methodology note written: {note_path}")

    print("\n" + "=" * 90)
    print(f"[INFO] {len(rows)} epitopes | {dict(class_counts)} | {dict(pathogen_counts)} | "
          f"{n_antigens} source antigens")
    print("=" * 90 + "\n")

    return out_path


if __name__ == "__main__":
    build_dossier()
