"""
phase1_common.py

Shared utilities for the Phase I epitope-selection pipeline (Steps 1A-1G).

WHY THIS FILE EXISTS:
Every Phase 1 script computed its project root as
os.path.join(script_dir, "..", "..", "..") -- three fixed hops. From
Research/Phase 1/STEP X/ that resolves to the folder ABOVE Research
itself, so every step's output was silently written outside the repo.
Phase 2 solved the same problem with an anchor-folder walk instead of a
fixed hop count (see Phase 2/_common/phase2_common.py); this file ports
that fix to Phase 1 so both phases resolve the same way.
"""

import os
import sys

RESEARCH_ANCHOR = "Research"


def resolve_project_root(script_file):
    """
    Walk upward from a script's own location until a folder literally
    named "Research" is found. Anchoring on a named folder (rather than
    a fixed hop count) keeps working regardless of how deep a given
    STEP script sits -- which is exactly what broke the fixed "../../.."
    version of this logic.
    """
    script_dir = os.path.dirname(os.path.abspath(script_file))
    current = script_dir
    while os.path.basename(current) != RESEARCH_ANCHOR:
        parent = os.path.dirname(current)
        if parent == current:
            print(f"\n[FATAL ERROR] Could not locate a '{RESEARCH_ANCHOR}' anchor folder above: {script_dir}")
            sys.exit(1)
        current = parent
    return current


def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"


def print_banner(text, width=90):
    print("\n" + "=" * width)
    print(f"{text:^{width}}")
    print("=" * width)


def latest_file(folder, suffix=".csv"):
    """
    Returns the path to the most recently created file with the given
    suffix in `folder`, or None if the folder is missing or empty.
    Callers must check for None -- this never raises on a missing dir.
    """
    if not os.path.isdir(folder):
        return None
    candidates = [f for f in os.listdir(folder) if f.endswith(suffix)]
    if not candidates:
        return None
    candidates.sort(key=lambda f: os.path.getctime(os.path.join(folder, f)))
    return os.path.join(folder, candidates[-1])


# =============================================================================
# HIV SUBUNIT RANGES -- one definition, read by 1A's provenance check and 1De.
#
# WHY THIS EXISTS: Phase 1A defines its four HIV targets as four Entrez QUERY
# STRINGS ("gp120 AND ...", "gp41 AND ...", "p17 AND ...", "p24 AND ..."), and
# stores whatever record NCBI returns whole, labelled with the query name. For
# CRF01_AE, NCBI returns POLYPROTEINS: HIV_gp120_Var_01 and HIV_gp41_Var_01 are
# the same 856-aa Env (YES72107.1), and HIV_p17_Var_01 and HIV_p24_Var_01 are
# the same 1437-aa Gag-Pol (YES72110.1). Every epitope sliced from them
# therefore inherited its query's label regardless of where in the polyprotein
# it actually lies -- see the Phase 1A provenance-correction report.
#
# PROVENANCE OF THE NUMBERS: the GenBank records carry NO mat_peptide features,
# so boundaries come from a global BLOSUM62 alignment of each Var_01 sequence
# against the annotated HXB2 reference, and the reference's own UniProt CHAIN
# features mapped through that alignment:
#   Env      YES72107.1 (856 aa)  vs  P04578 ENV_HV1H2 (856 aa)
#            CHAIN 33..511  "Surface protein gp120"      -> Var_01  32..504
#            CHAIN 512..856 "Transmembrane protein gp41" -> Var_01 505..856
#   Gag-Pol  YES72110.1 (1437 aa) vs  P04585 POL_HV1H2 (1435 aa)
#            CHAIN 2..132   "Matrix protein p17"         -> Var_01   2..135
#            CHAIN 133..363 "Capsid protein p24"         -> Var_01 136..366
# Independently corroborated by sequence landmarks in the stored Var_01 files:
# the Env furin site + fusion peptide "...RVVERPKR | AVGIGAMIFGF" puts gp41 at
# 505, and the Gag MA/CA junction "...SQNY | PIVQ" puts p24 at 136, ending at
# "...KARVL" (366). Both methods agree exactly.
#
# 1-BASED AND INCLUSIVE, matching UniProt/GenBank convention. Slice with
# seq[start-1:end].
# =============================================================================
HIV_SUBUNIT_RANGES = {
    "HIV_gp120": {"parent": "Env",     "accession": "YES72107.1", "start": 32,  "end": 504},
    "HIV_gp41":  {"parent": "Env",     "accession": "YES72107.1", "start": 505, "end": 856},
    "HIV_p17":   {"parent": "Gag-Pol", "accession": "YES72110.1", "start": 2,   "end": 135},
    "HIV_p24":   {"parent": "Gag-Pol", "accession": "YES72110.1", "start": 136, "end": 366},
}

# Targets whose Var_01 file IS the mature protein -- no slicing needed.
SINGLE_PROTEIN_TARGETS = ("Mpox_A35R", "Mpox_B5R", "Mpox_L1R")


def subunit_of(target, position_1based):
    """
    Given a 1-based offset into a target's Var_01 parent record, returns the
    name of the HIV subunit that offset actually falls in, or None if the
    target is not one of the polyprotein-derived HIV four (Mpox targets are
    already mature proteins, so their label is always correct).

    Returns "Gag_downstream" for Gag-Pol offsets past p24 (p2/p7/p6/pol) --
    those are real regions, just not one of this study's four HIV antigens.
    """
    if target in SINGLE_PROTEIN_TARGETS:
        return target
    info = HIV_SUBUNIT_RANGES.get(target)
    if info is None:
        return None
    siblings = [(n, d) for n, d in HIV_SUBUNIT_RANGES.items() if d["parent"] == info["parent"]]
    for name, d in siblings:
        if d["start"] <= position_1based <= d["end"]:
            return name
    if info["parent"] == "Gag-Pol" and position_1based > 366:
        return "Gag_downstream"
    return None


# =============================================================================
# LANDMARK-BASED SUBUNIT BOUNDARIES FOR NON-Var_01 RECORDS.
#
# HIV_SUBUNIT_RANGES above is exact but applies ONLY to the Var_01 records it
# was aligned against. Epitopes selected from other variants live in records of
# different length (Env variants here run 854-868 aa), so those fixed offsets
# do not transfer. Rather than align every variant, locate the SAME cleavage
# landmarks the alignment confirmed on Var_01:
#
#   Env      gp120 | gp41 at the host-furin site immediately followed by the
#            gp41 fusion peptide: "...RVVERPKR | AVGIGAMIFGF". On Var_01 this
#            puts gp41 at 505 -- identical to the alignment-derived boundary.
#   Gag-Pol  p17 | p24 at the MA/CA junction "...SQNY | PIVQ", and p24's C
#            terminus at "...KARVL". On Var_01 these give 136 and 366 --
#            again identical to the alignment-derived boundaries.
#
# Both landmarks are the actual protease/furin recognition sites, so they track
# indels correctly in a way fixed offsets cannot. Returns None when a landmark
# is absent, and callers must treat that as UNRESOLVED rather than guessing.
# =============================================================================
def subunit_boundaries_by_landmark(seq, parent):
    """
    Returns {subunit_name: (start, end)} 1-based inclusive, or None if the
    defining landmark is absent. Used ONLY to corroborate the alignment-derived
    boundaries -- never as the primary source.

    CORRECTED 2026-08-30. The previous implementation had two defects, both
    found by testing it against live CRF01_AE records rather than assuming:

    (1) The Gag landmark was "SQNYPIVQ", the subtype-B / HXB2 MA-CA junction.
        CRF01_AE -- THIS STUDY'S OWN TARGET SUBTYPE -- does not carry it: real
        records read "...VSHNYPIVQ..." (S-H-N-Y, not S-Q-N-Y). Measured on 6
        live CRF01_AE Gag records: SQNYPIVQ absent from 6/6, so the function
        returned None for the ENTIRE pool and p17/p24 could never be resolved.
        Now anchored on "PIVQ" alone -- the capsid N-terminus. Pro-1 of CA is
        required to form the buried salt bridge with Asp-51 that folds the
        mature capsid beta-hairpin, so it is conserved where the upstream
        matrix residues are not. Verified present in 6/6 of the same records.

    (2) The Env branch returned HIV_gp120 as (1, fp) -- starting at residue 1,
        which folds the ~32-residue signal peptide into gp120. UniProt P04578
        annotates SIGNAL 1..32 and CHAIN 33..511 for gp120, so the mature
        surface protein starts at 33, not 1. This contradicted
        HIV_SUBUNIT_RANGES in this same file (start 32). The signal peptide is
        cleaved co-translationally and is not part of the mature antigen, so
        including it would have put 32 residues that never appear on the virion
        into the k-mer window pool. Signal-peptide cleavage is not reliably
        detectable from a short motif, so the gp120 START is left to the
        alignment mapping; this function reports only the gp120/gp41 boundary
        it CAN see, and returns start=None to force the caller to use the
        alignment rather than silently defaulting to 1.
    """
    if parent == "Env":
        fp = seq.find("AVGIG")           # gp41 fusion peptide = first residue of gp41
        if fp == -1:
            return None
        # start of gp120 is the signal-peptide cleavage site -- alignment only.
        return {"HIV_gp120": (None, fp), "HIV_gp41": (fp + 1, len(seq))}

    if parent in ("Gag", "Gag-Pol"):
        j = seq.find("PIVQ")             # capsid N-terminus; P is CA residue 1
        if j == -1:
            return None
        p24_start = j + 1                # 1-based position of that P
        ca_end = seq.find("KARVL", j)    # capsid C-terminus; L is the last residue
        if ca_end == -1:
            return None
        p24_end = ca_end + 5
        # p17 starts at 2 (initiator Met removed), per UniProt CHAIN 2..132.
        return {"HIV_p17": (2, p24_start - 1), "HIV_p24": (p24_start, p24_end)}

    return None


# =============================================================================
# ANTIGEN IDENTITY GATE -- ADDED 2026-08-30
#
# WHY THIS EXISTS. The construct Vax_Final_6f34b53e shipped with 12 of its 32
# epitopes drawn from proteins the study never intended to target, because
# Phase 1A selected antigens by literal gene-name STRING and stored whatever
# NCBI returned. Three separate ways that failed:
#
#   * MPXV has its OWN gene called B5R that is unrelated to vaccinia's B5R.
#     Querying "B5R" returned OPG189, an intracellular ankyrin-repeat protein,
#     100% identity to the wrong thing and 25.2% to OPG190/B6R, the actual EEV
#     glycoprotein the study wanted. The query was answered correctly; the
#     question was wrong.
#   * NCBI gives OPG095 and OPG053 the IDENTICAL product name "IMV membrane
#     protein L1R". Verified: XYR87238.1 (250 aa, /gene="OPG095") and
#     XYR87197.1 (212 aa, /gene="OPG053") share Protein Name AND Title. No
#     name-based query can separate them -- only /gene differs.
#   * The four HIV "antigens" were two polyproteins, so an epitope's Target
#     recorded which QUERY fetched its parent, not which protein it sits in.
#
# THE GATE. Every candidate is aligned against a PANEL holding the intended
# reference AND the known decoys. A record is accepted only if the intended
# reference is its BEST match and clears the floor. A wrong protein therefore
# surfaces as a measured result rather than being assumed away -- and it is
# rejected before it can enter the pipeline, not audited afterwards.
#
# Kept here rather than in Phase1A so retrieval and the standalone provenance
# audit share ONE definition and cannot drift apart.
# =============================================================================

# Identity floors are per-target because the comparisons are not equivalent.
# Mpox records are the same species as their reference, so genuine variants sit
# near-identical and 90% is a wide margin that still rejects the ~25% decoys
# decisively. The HIV references are subtype-B HXB2 while the study targets
# CRF01_AE, so a genuine CRF01_AE record is only ~70-80% identical to its own
# reference -- the observed minimum on the previous pool was 73.4%. Setting the
# HIV floor at Mpox levels would reject the study's entire target subtype.
ANTIGEN_REFERENCES = {
    "Mpox_A35R": {
        "intended": "A0A7H0DND2",
        "gene": "OPG161",
        "desc": "OPG161 -- VACV A33R ortholog, EEV envelope glycoprotein",
        "decoys": {},
        "floor": 90.0,
        "length_band": (170, 195),
    },
    "Mpox_L1R": {
        "intended": "M1LBP0",
        "gene": "OPG095",
        "desc": "OPG095 / M1R -- VACV L1R ortholog, IMV membrane protein",
        # OPG053 is the gene NCBI labels with the identical product name.
        "decoys": {"A0A7H0DN26": "OPG053 (EFC-associated protein -- same NCBI product name)"},
        "floor": 90.0,
        "length_band": (240, 262),
    },
    "Mpox_B5R": {
        "intended": "P0DTN2",
        "gene": "OPG190",
        "desc": "OPG190 / B6R -- VACV B5R ortholog, EEV envelope glycoprotein",
        # OPG189 is MPXV's own unrelated gene that is literally named B5R.
        "decoys": {"A0A7H0DNF1": "OPG189 (ankyrin repeat protein -- MPXV's own gene named B5R)"},
        "floor": 90.0,
        "length_band": (300, 332),
    },
    # HIV parents. gp120/gp41 are sliced from Env; p17/p24 from Gag.
    "HIV_Env": {
        "intended": "P04578",
        "gene": None,
        "desc": "Env gp160 precursor (HXB2 reference)",
        "decoys": {},
        "floor": 65.0,
        "length_band": (820, 900),
    },
    "HIV_Gag": {
        # P04591 (Gag, 500 aa) NOT P04585 (Gag-Pol, 1435 aa): the records this
        # query actually returns are ~493-aa Gag proteins, and aligning 493 aa
        # against a 1435-aa Gag-Pol invites spurious placement. Both carry the
        # identical p17/p24 CHAIN coordinates (2..132 / 133..363), verified
        # against UniProt on 2026-08-30, so nothing is lost by using the
        # shorter, better-matched reference.
        "intended": "P04591",
        "gene": None,
        "desc": "Gag polyprotein (HXB2 reference)",
        "decoys": {},
        "floor": 65.0,
        "length_band": (470, 520),
    },
}

# UniProt CHAIN coordinates on the REFERENCE, mapped onto each query record
# through a global alignment. Confirmed live against UniProt 2026-08-30:
#   P04578  SIGNAL 1..32 | gp120 CHAIN 33..511 | gp41 CHAIN 512..856
#   P04591  p17 CHAIN 2..132 | p24 CHAIN 133..363
HIV_SUBUNIT_CHAINS = {
    "HIV_Env": {"HIV_gp120": (33, 511), "HIV_gp41": (512, 856)},
    "HIV_Gag": {"HIV_p17": (2, 132), "HIV_p24": (133, 363)},
}

# Canonical mature-subunit lengths, used as an independent sanity band on the
# slices the alignment produces. A slice outside its band is reported, never
# silently kept.
HIV_SUBUNIT_LENGTH_BAND = {
    "HIV_gp120": (455, 510),
    "HIV_gp41":  (330, 365),
    "HIV_p17":   (125, 140),
    "HIV_p24":   (225, 240),
}

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def sanitise_sequence(seq):
    return "".join(seq.split()).upper().replace("*", "")


def has_nonstandard(seq):
    """True if the sequence carries anything outside the 20 standard residues."""
    return any(ch not in STANDARD_AA for ch in seq)


def fetch_uniprot_fasta(accession, cache_dir):
    """Downloads a UniProt reference once and caches it on disk."""
    import urllib.request
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{accession}.fasta")
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            return sanitise_sequence("".join(l.strip() for l in fh if not l.startswith(">")))
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    with urllib.request.urlopen(url, timeout=60) as resp:
        text = resp.read().decode()
    if not text.startswith(">"):
        raise RuntimeError(f"UniProt returned no FASTA for {accession}")
    tmp = path + ".part"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)
    return sanitise_sequence("".join(l.strip() for l in text.splitlines() if not l.startswith(">")))


_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        a = Align.PairwiseAligner(mode="global", open_gap_score=-11, extend_gap_score=-1)
        a.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
        _ALIGNER = a
    return _ALIGNER


def global_identity(query, reference):
    """Percent identity over the SHORTER sequence, BLOSUM62 global alignment.

    Same definition Phase1A_provenanceCorrection.py used, kept identical so the
    audit and the gate cannot disagree about what "identity" means.
    """
    aln = _aligner().align(query, reference)[0]
    a, b = str(aln[0]), str(aln[1])
    matches = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    return 100.0 * matches / min(len(query), len(reference))


def best_reference_match(seq, target, cache_dir):
    """
    Scores a sequence against the intended reference AND every known decoy.

    Returns (best_accession, best_pct, intended_pct, [(acc, pct), ...]).
    Reporting the whole panel rather than only the intended score is the point:
    it is what turns "this is B5R because we asked for B5R" into "this is
    OPG189 at 100.0% and OPG190 at 25.2%, so it is not the antigen we wanted".
    """
    spec = ANTIGEN_REFERENCES[target]
    panel = {spec["intended"]: spec["desc"]}
    panel.update(spec["decoys"])
    scored = []
    for acc in panel:
        ref = fetch_uniprot_fasta(acc, cache_dir)
        scored.append((acc, global_identity(seq, ref)))
    scored.sort(key=lambda t: t[1], reverse=True)
    best_acc, best_pct = scored[0]
    intended_pct = next(p for a, p in scored if a == spec["intended"])
    return best_acc, best_pct, intended_pct, scored


def identity_verdict(seq, target, cache_dir):
    """
    The gate itself. Returns (accepted: bool, reason: str, detail: dict).
    A record must clear ALL of: standard residues only, length band, intended
    reference is the best match, and identity above the per-target floor.
    """
    spec = ANTIGEN_REFERENCES[target]
    detail = {"length": len(seq)}
    if has_nonstandard(seq):
        bad = sorted({c for c in seq if c not in STANDARD_AA})
        detail["nonstandard"] = "".join(bad)
        return False, f"NONSTANDARD_RESIDUES({''.join(bad)})", detail
    lo, hi = spec["length_band"]
    if not (lo <= len(seq) <= hi):
        return False, f"LENGTH_OUT_OF_BAND({len(seq)} not in {lo}-{hi})", detail
    best_acc, best_pct, intended_pct, scored = best_reference_match(seq, target, cache_dir)
    detail.update({"best_acc": best_acc, "best_pct": round(best_pct, 2),
                   "intended_acc": spec["intended"], "intended_pct": round(intended_pct, 2),
                   "panel": "; ".join(f"{a}:{p:.1f}%" for a, p in scored)})
    if best_acc != spec["intended"]:
        return False, f"WRONG_PROTEIN(best={best_acc} @{best_pct:.1f}%)", detail
    if intended_pct < spec["floor"]:
        return False, f"BELOW_FLOOR({intended_pct:.1f}% < {spec['floor']}%)", detail
    return True, "ACCEPTED", detail


# ---- diversity-preserving selection -----------------------------------------
def _kmer_set(seq, k=4):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def _jaccard_distance(a, b):
    inter = len(a & b)
    union = len(a | b)
    return 1.0 - (inter / union) if union else 0.0


def diversity_select(seqs, n):
    """
    Farthest-first traversal over distinct sequences.

    WHY NOT JUST TAKE THE FIRST N: querying by /gene returns ~8,800 MPXV records
    per antigen, overwhelmingly near-identical 2022-outbreak genomes deposited
    in batches. Taking a contiguous block would hand Phase 1Dc a pool with
    almost no variation, and every epitope would score ~100% conservancy -- a
    measurement of NCBI's deposit pattern, not of cross-strain conservation.

    Seeded on the most frequently deposited exact sequence (the consensus
    strain, which becomes Var_01), then repeatedly adds whichever remaining
    sequence is FARTHEST from everything already chosen. Distance is Jaccard
    over 4-mer sets: length-robust, no all-pairs alignment needed.
    """
    from collections import Counter
    counts = Counter(seqs)
    distinct = list(counts)
    if len(distinct) <= n:
        return sorted(distinct, key=lambda s: (-counts[s], s))
    profiles = {s: _kmer_set(s) for s in distinct}
    seed = max(distinct, key=lambda s: (counts[s], -len(s)))
    chosen = [seed]
    mind = {s: _jaccard_distance(profiles[s], profiles[seed]) for s in distinct if s != seed}
    while len(chosen) < n and mind:
        nxt = max(mind, key=lambda s: (mind[s], counts[s]))
        chosen.append(nxt)
        del mind[nxt]
        for s in list(mind):
            d = _jaccard_distance(profiles[s], profiles[nxt])
            if d < mind[s]:
                mind[s] = d
    return chosen


def map_reference_span(query_seq, reference_seq, ref_start, ref_end):
    """
    Maps a 1-based inclusive span on the REFERENCE onto the QUERY through a
    global BLOSUM62 alignment. Returns (q_start, q_end) 1-based inclusive, or
    None if neither endpoint lands in an aligned block.

    This is how HIV mature-subunit boundaries are derived PER RECORD. Fixed
    offsets cannot be reused: the previous run pinned them to YES72107.1
    (856 aa) and YES72110.1 (1437 aa), and the records NCBI returns today are
    860-aa Env and 493-aa Gag, so those offsets no longer refer to anything.
    """
    aln = _aligner().align(query_seq, reference_seq)[0]
    q_blocks, r_blocks = aln.aligned
    pairs = list(zip(q_blocks, r_blocks))

    def to_query(ref_pos_0, prefer_forward):
        # exact hit inside an aligned block
        for (qs, qe), (rs, re) in pairs:
            if rs <= ref_pos_0 < re:
                return qs + (ref_pos_0 - rs)
        # reference position falls in a gap -- take the nearest aligned anchor
        if prefer_forward:
            best = None
            for (qs, qe), (rs, re) in pairs:
                if rs >= ref_pos_0 and (best is None or rs < best[1]):
                    best = (qs, rs)
            return best[0] if best else None
        best = None
        for (qs, qe), (rs, re) in pairs:
            if re <= ref_pos_0 and (best is None or re > best[1]):
                best = (qe - 1, re)
        return best[0] if best else None

    qs0 = to_query(ref_start - 1, True)
    qe0 = to_query(ref_end - 1, False)
    if qs0 is None or qe0 is None or qe0 < qs0:
        return None
    return qs0 + 1, qe0 + 1


# =============================================================================
# ALLERGEN HOMOLOGY SCREEN -- ADDED 2026-08-31
#
# WHY THIS EXISTS. Sec. I.E screens every candidate peptide for allergenicity
# with AllerTOP and AllergenFP. Both are whole-protein predictors, and measured
# on this project's own candidates they cannot operate at the lengths involved:
#
#   AllergenFP  refuses anything below 16 aa outright (deviation #6), so it can
#               never see the MHC-I 9/10-mers or the MHC-II 15-mers.
#   AllerCatPro 2.0  measured with a positive control on 2026-08-31: the SAME
#               allergen (Bet v 1) was detected at 160, 16 and 15 aa and became
#               INVISIBLE at 10 and 9 aa, returning "No significant hit (E-value
#               threshold 0.001)". A short query cannot reach BLAST significance;
#               that is a statistical floor, not a biological verdict.
#   AlgPred 2.0 measured on 1,600 labelled peptides: AUC 0.611, barely above
#               random, and 78.9% false positives on human peptides at its
#               default threshold.
#
# So a "no allergen" call on a 9-mer from any of those tools means "cannot be
# assessed", not "safe". Reporting it as safety would be reporting a number the
# method was not designed to produce -- the same error deviation #1 corrected
# for raw ToxinPred scores.
#
# WHAT THIS DOES INSTEAD. It applies the FAO/WHO Codex Alimentarius criterion:
# does the peptide share a contiguous identical stretch with a KNOWN allergen?
# That is a factual lookup rather than a prediction, and it has NO minimum
# length -- it works on a 9-mer where every predictor above fails.
#
# WINDOW CHOICE. Codex specifies >=6 contiguous identical residues. That rule is
# documented as over-sensitive (Silvanovich et al. 2006, Toxicol Sci; Hileman et
# al. 2002) and this project reproduced that independently: at w=6 it flags 6.1%
# of ordinary HUMAN peptides. Default is therefore w=8, where the human
# false-positive rate falls to 0.4% while positive controls are still detected
# 100% of the time at every length from 9 to 16 aa. The window is a parameter,
# and w=6 remains available for anyone wanting the literal Codex rule.
#
# WHAT A RESULT MEANS -- BOTH DIRECTIONS MATTER:
#   HIT  = positive evidence of allergen homology. Actionable.
#   NONE = no match in THIS reference set. NOT proof of safety: the set holds
#          1,020 reviewed allergens, novel allergens are in no database, and
#          only exact matches are detected (an 80%-similar peptide passes).
# =============================================================================

ALLERGEN_MATCH_WINDOW = 8
_ALLERGEN_BLOB = None


def load_allergen_reference(project_root=None):
    """Loads the curated allergen set as one delimiter-joined blob for
    substring search. Cached after first call."""
    global _ALLERGEN_BLOB
    if _ALLERGEN_BLOB is not None:
        return _ALLERGEN_BLOB
    if project_root is None:
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
    # Prefer the merged reference (UniProtKB KW-0020 + the AllerTOP v2 allergen
    # training set, 2,769 unique records / 665,125 residues) and fall back to
    # the UniProt-only set. A larger reference makes a NEGATIVE stronger, since
    # "no match" is only as informative as the set searched -- this is the
    # upgrade allergen_db/README.txt flagged. See allergen_db/README.txt for the
    # measured false-positive cost of the larger set.
    path = os.path.join(project_root, "allergen_db",
                        "allergen_reference_merged.fasta")
    if not os.path.isfile(path):
        path = os.path.join(project_root, "allergen_db",
                            "uniprot_reviewed_allergens.fasta")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Allergen reference not found at {path}. See allergen_db/README.txt "
            "for the UniProt query that regenerates it.")
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    # "#" cannot occur in a sequence, so joining on it prevents a match that
    # spans the boundary between two different allergens.
    _ALLERGEN_BLOB = "#".join(seqs)
    return _ALLERGEN_BLOB


def allergen_homology(peptide, window=ALLERGEN_MATCH_WINDOW, project_root=None):
    """
    FAO/WHO-style contiguous-identity check against known allergens.

    Returns (hit: bool, matched_kmer: str|None). A peptide shorter than the
    window cannot be assessed and returns (False, None) -- callers must not read
    that as a negative result.
    """
    pep = sanitise_sequence(peptide)
    if len(pep) < window:
        return False, None
    blob = load_allergen_reference(project_root)
    for i in range(len(pep) - window + 1):
        kmer = pep[i:i + window]
        if kmer in blob:
            return True, kmer
    return False, None


# =============================================================================
# HUMAN SELF-HOMOLOGY SCREEN -- exact contiguous-identity check.
# Added 2026-08-31, deviation TBD (Phase 1Ec restructure).
#
# WHY THIS EXISTS ALONGSIDE THE BLASTP SCREEN, NOT INSTEAD OF IT. Sec. I.E of
# the manuscript specifies BLASTP against the reviewed human proteome with an
# explicit 70% identity / 70% coverage / E<=1e-5 acceptance rule. That rule
# stays and Phase 1Ec still runs it, because it is the cited method. But it was
# TESTED on 20 peptides cut verbatim from human Swiss-Prot -- guaranteed
# perfect self-matches by construction, at this study's own candidate lengths
# (9/10/15/16 aa, 5 each) -- and returned 0/20 hits at default blastp settings,
# rising to 13/20 with -task blastp-short, but STILL 0/20 that clear the
# E<=1e-5 bar: the best E-value achievable on a 9-mer against this database is
# ~0.2, three orders of magnitude short of significance. On the real 108
# Phase 1Eb survivors the BLASTP screen (even with blastp-short) returns 0
# hits meeting the paper's threshold -- not "0 excluded", 0 peptides for which
# the screen produces ANY usable signal. A short peptide cannot reach BLAST
# statistical significance regardless of how self-like it truly is; this is
# the identical failure mode already documented for AllerCatPro at epitope
# length (deviation #27) and for the allergen predictors generally.
#
# So this screen supplies what BLASTP cannot at these lengths: a factual
# lookup for an identical contiguous stretch, which has no minimum length and
# no significance threshold to fail. It is the same technique already
# validated for the allergen screen (allergen_homology above), applied to the
# human reference proteome instead of the allergen reference.
#
# WINDOW CALIBRATED AGAINST A SCRAMBLED-SEQUENCE CONTROL, NOT ASSUMED.
# Composition-matched scrambled controls (5 random permutations per real
# candidate peptide, same amino acid multiset, same length, no biology) were
# BLASTP'd and substring-matched alongside the real candidates:
#     window 6   real 80.6%   scrambled 80.7%   -- NO SEPARATION, useless
#     window 7   real 14.8%   scrambled 15.0%   -- NO SEPARATION, useless
#     window 8   real  2.8%   scrambled  0.9%   -- separates, ~3x background
#     window 9+  real  0.0%   scrambled ~0.0%   -- no positives either way
# The same BLASTP-relaxed-threshold experiment (identity+coverage only, E-value
# dropped) was tried first and REJECTED: at 70/70 the real pool flagged 24.1%
# and the scrambled controls flagged 25.2% -- statistically indistinguishable,
# i.e. pure background matching human protein composition, not self-homology.
# Loosening BLAST's thresholds to get ANY signal at these lengths manufactures
# noise, not evidence; window=8 exact matching is what actually separates
# signal from background here, exactly as it did for the allergen screen.
#
# WHAT A RESULT MEANS -- BOTH DIRECTIONS MATTER, SAME CAVEATS AS THE ALLERGEN
# SCREEN:
#   HIT  = positive evidence of an identical 8+ residue stretch shared with a
#          reviewed human protein. Actionable -- flag for exclusion or review.
#   NONE = no exact match in this reference set. NOT proof the peptide is
#          non-self: near-identical (not exact) matches are invisible to this
#          method, and this checks reviewed Swiss-Prot human proteins only.
# =============================================================================

HUMAN_SELF_MATCH_WINDOW = 8
_HUMAN_SELF_BLOB = None


def load_human_reference(project_root=None):
    """Loads the reviewed human Swiss-Prot proteome as one delimiter-joined
    blob for substring search. Cached after first call. Reuses the same FASTA
    Phase 1Ec's BLASTP screen is built from (human_swissprot_db/), so both
    screens run against the identical reference set."""
    global _HUMAN_SELF_BLOB
    if _HUMAN_SELF_BLOB is not None:
        return _HUMAN_SELF_BLOB
    if project_root is None:
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
    path = os.path.join(project_root, "human_swissprot_db", "human_swissprot.fasta")
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Human reference proteome not found at {path}. This is the same "
            "FASTA human_swissprot_db/ was built from for the BLASTP self-"
            "homology screen -- see that folder's provenance notes.")
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    # "#" cannot occur in a sequence, so joining on it prevents a match that
    # spans the boundary between two different proteins.
    _HUMAN_SELF_BLOB = "#".join(seqs)
    return _HUMAN_SELF_BLOB


def human_self_homology(peptide, window=HUMAN_SELF_MATCH_WINDOW, project_root=None):
    """
    Exact contiguous-identity check against the reviewed human proteome --
    the human-reference counterpart to allergen_homology() above. Same
    method, same rationale: BLASTP cannot reach statistical significance on a
    9-16 aa query regardless of true self-similarity, so an identical-stretch
    lookup (no length floor, no significance threshold) is used instead.

    Returns (hit: bool, matched_kmer: str|None). A peptide shorter than the
    window cannot be assessed and returns (False, None) -- callers must not
    read that as a negative result.
    """
    pep = sanitise_sequence(peptide)
    if len(pep) < window:
        return False, None
    blob = load_human_reference(project_root)
    for i in range(len(pep) - window + 1):
        kmer = pep[i:i + window]
        if kmer in blob:
            return True, kmer
    return False, None
