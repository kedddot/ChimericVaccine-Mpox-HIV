"""
Phase1A_sequenceRetrival.py -- REWRITTEN 2026-08-30

WHAT WAS WRONG WITH THE PREVIOUS VERSION
----------------------------------------
It defined each antigen as an Entrez QUERY STRING built from a gene NAME:

    field_query = f'("{gene}"[Protein Name] OR "{gene}"[Title]) AND "Monkeypox virus"[Organism]'

and wrote whatever NCBI returned straight to disk, labelled with the query name.
A gene name is not a protein identity. There was no length check, no identity
check, no accession pinning -- `grep -i "validat|verify|identity"` over the old
file returned nothing. Every downstream step then trusted those labels forever.

The construct that resulted, Vax_Final_6f34b53e, carried 12 of its 32 epitopes
(161 of 566 residues) from proteins the study never intended to target:

  * Mpox_B5R was 100.0% OPG189, an intracellular ANKYRIN-REPEAT protein, and
    only 25.2% identical to OPG190/B6R, the EEV envelope glycoprotein actually
    wanted. Orthopoxvirus gene names are not conserved across species: vaccinia's
    B5R ortholog is B6R/OPG190 in MPXV, and MPXV has its OWN unrelated gene
    called B5R. The query returned exactly what it asked for. 0 of 25 records
    were the intended protein.
  * 16 of 30 Mpox_L1R records were OPG053, a different entry-fusion-complex
    protein only 25.9% identical to L1R. NCBI gives both genes the IDENTICAL
    product name -- verified, XYR87238.1 (/gene="OPG095") and XYR87197.1
    (/gene="OPG053") share both Protein Name and Title. NO name-based query can
    separate them. Only /gene differs.
  * The four HIV targets were really TWO polyproteins. gp120 and gp41 returned
    the same Env record; p17 and p24 the same Gag record. Every k-mer inherited
    the label of the QUERY that fetched its parent, not the subunit it lies in,
    so 9 of 14 HIV epitopes were mislabelled.

WHAT THIS VERSION DOES INSTEAD
------------------------------
1. QUERIES MPOX BY /gene (OPG number), never by protein name.
2. GATES EVERY RECORD ON IDENTITY BEFORE IT REACHES DISK. Each candidate is
   aligned against a panel holding the intended UniProt reference AND the known
   decoys; it is accepted only if the intended reference is its BEST match and
   clears a per-target floor. A wrong protein surfaces as a measured result and
   is rejected -- it is no longer possible for one to enter the pipeline.
3. SLICES THE HIV POLYPROTEINS TO MATURE SUBUNITS HERE, at retrieval, so one
   FASTA is one mature protein exactly as the Mpox targets already are. Every
   downstream label is then correct by construction rather than by assumption.
   Boundaries are derived PER RECORD by mapping UniProt CHAIN coordinates
   through a global alignment, and independently corroborated by cleavage
   landmarks.
4. SELECTS FOR DIVERSITY. Querying by /gene returns ~8,800 MPXV records per
   antigen instead of the 25-27 a name query found, so the binding constraint
   flipped from scarcity to redundancy. Taking the first 30 would hand Phase 1Dc
   a pool of near-identical 2022-outbreak genomes and every epitope would score
   ~100% conservancy -- a measurement of NCBI's deposit pattern, not of
   cross-strain conservation.
5. WRITES AN AUDIT ROW FOR EVERY CANDIDATE EXAMINED, accepted or rejected, with
   its identity against the whole reference panel and the reason for rejection.

The shared gate lives in phase1_common.py so this script and the standalone
Phase1A_provenanceCorrection.py audit cannot drift apart.
"""

import os
import sys
import csv
import time
import urllib.parse
import urllib.request
from datetime import datetime
from collections import defaultdict
import socket
from Bio import Entrez
from Bio import SeqIO

# Bound every network call -- see esearch_ids() for why this is not optional.
socket.setdefaulttimeout(120)

# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase1_common module.
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
# CONFIGURATION
# =============================================================================
Entrez.email = os.getenv("NCBI_EMAIL", "enzoleonor.3309@gmail.com")

VARIANTS_PER_TARGET = 30

# How many records to actually download and gate per source query. The full
# result set is far larger (~8,800 for each Mpox gene, ~12,000 for HIV Env);
# IDs are cheap to list but records are not, so the ID list is retrieved in
# full and then STRIDE-SAMPLED across its whole span. Sampling evenly rather
# than taking a contiguous block matters: NCBI returns IDs in accession order,
# and accession order is deposit order, so the first N records are typically
# one submission batch of near-identical genomes.
CANDIDATE_POOL_SIZE = 250
EFETCH_BATCH = 50

# Mpox: retrieved by /gene. The OPG identifier is the ONLY field that
# distinguishes OPG095 from OPG053, which share a product name.
MPOX_SOURCES = {
    "Mpox_A35R": "OPG161",
    "Mpox_L1R":  "OPG095",
    "Mpox_B5R":  "OPG190",
}

# HIV: retrieved as PARENTS, then sliced. Both subunits of a parent come from
# the same record set, which is honest about what NCBI actually holds -- the
# previous version issued four separate queries and pretended it had four
# independent antigens when two records were doing all four jobs.
HIV_SOURCES = {
    "HIV_Env": 'gp120 AND "HIV-1"[Organism] AND CRF01_AE',
    "HIV_Gag": 'p24 AND "HIV-1"[Organism] AND CRF01_AE',
}
HIV_SUBUNITS_OF = {"HIV_Env": ("HIV_gp120", "HIV_gp41"),
                   "HIV_Gag": ("HIV_p17", "HIV_p24")}

CLADE_PATTERN_HINTS = ("clade iia", "clade iib", "clade ii", "west african",
                       "genotype: wa", "genotype:wa")


def classify_clade(record):
    """Clade signal from the record's description and source qualifiers.

    NCBI has no taxonomy node for Clade IIa/IIb (Monkeypox virus is one taxid,
    10244), so clade is only recoverable from annotation text. Both confirmed
    and unconfirmed are kept and labelled, per deviation #4.
    """
    haystacks = [record.description.lower()]
    for feature in record.features:
        if feature.type == "source":
            for key in ("note", "isolate", "strain", "country", "collection_date"):
                for val in feature.qualifiers.get(key, []):
                    haystacks.append(str(val).lower())
    blob = " | ".join(haystacks)
    for hint in CLADE_PATTERN_HINTS[:3]:
        if hint in blob:
            return "Clade " + hint.replace("clade ", "").upper(), True
    for hint in CLADE_PATTERN_HINTS[3:]:
        if hint in blob:
            return "Clade II (legacy West African designation)", True
    return "Unconfirmed", False


def gene_qualifier(record):
    """The /gene value, which is the field that actually identifies the gene."""
    for feature in record.features:
        if feature.type in ("CDS", "gene", "Protein"):
            for key in ("gene", "locus_tag"):
                vals = feature.qualifiers.get(key, [])
                if vals:
                    return str(vals[0])
    return ""


# =============================================================================
# ENTREZ HELPERS
# =============================================================================
def esearch_ids(query, retmax=40000):
    """
    NOTE: Biopython's Entrez wrappers set NO socket timeout, so a stalled NCBI
    connection blocks forever with the process at 0% CPU. A first attempt at
    this run hung for 10 minutes that way. socket.setdefaulttimeout() above
    bounds every Entrez call, and each is retried.
    """
    for attempt in range(4):
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=retmax)
            result = Entrez.read(handle)
            handle.close()
            return result["IdList"], int(result["Count"])
        except Exception as exc:
            print(f"    [retry {attempt+1}/4] esearch failed: {exc}", flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"esearch failed after 4 attempts: {query}")


def stride_sample(ids, n):
    """Evenly sample n IDs across the full ordered list (deposit-order spread)."""
    if len(ids) <= n:
        return list(ids)
    step = len(ids) / float(n)
    return [ids[min(int(i * step), len(ids) - 1)] for i in range(n)]


def efetch_records(id_list, rettype="fasta"):
    """
    Fetches records in batches, with retries.

    WHY TWO PASSES (rettype matters for speed). A GenBank record carries the
    /gene qualifier and the source features needed for clade tagging, but it is
    an order of magnitude heavier than FASTA and NCBI serves it far more slowly
    -- a 100-record GenBank batch was taking minutes, which would have made the
    seven source queries an hour of pure waiting.

    Only the 30 records actually SELECTED per target need that annotation. So
    the candidate pool is screened in FASTA (identity is a sequence property and
    needs nothing else), and GenBank is fetched only for the final picks. Same
    checks, same evidence, a fraction of the transfer.
    """
    fmt = "genbank" if rettype == "gb" else "fasta"
    out = []
    for i in range(0, len(id_list), EFETCH_BATCH):
        chunk = id_list[i:i + EFETCH_BATCH]
        print(f"      efetch[{rettype}] {i+1}-{min(i+EFETCH_BATCH,len(id_list))}"
              f" of {len(id_list)}", flush=True)
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(chunk),
                                       rettype=rettype, retmode="text")
                out.extend(list(SeqIO.parse(handle, fmt)))
                handle.close()
                break
            except Exception as exc:
                print(f"      [retry {attempt+1}/3] efetch: {exc}", flush=True)
                if attempt == 2:
                    print(f"[WARN] efetch batch failed after 3 attempts: {exc}", flush=True)
                else:
                    time.sleep(3 * (attempt + 1))
        time.sleep(0.34)
    return out


def fetch_all_cached(query, label, cache_dir):
    """
    Fetches EVERY record for a query, caching to disk so a re-run is instant.

    WHY EXHAUSTIVE FOR MPOX. The first version stride-sampled 250 IDs evenly
    across the ~8,800 available and found only 2-8 DISTINCT protein sequences
    per gene. That was a sampling artifact, not biology: MPXV protein records
    are dominated by a clonal 2022-outbreak lineage, so an even sample lands
    almost entirely inside it, while the genuinely divergent isolates are rare
    and get missed. Fetching the whole set is the only way to see them. It is
    ~8,800 records per gene in batches of 500 -- around 18 requests, trivially
    within NCBI's limits -- and the result is cached, so this cost is paid once.

    HIV is deliberately NOT exhaustive: its 250-record sample already yields
    153-154 distinct sequences per subunit, so there is nothing to recover.
    """
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f"{label}_all.fasta")
    if os.path.isfile(cache) and os.path.getsize(cache) > 0:
        recs = list(SeqIO.parse(cache, "fasta"))
        print(f"    cache hit: {len(recs)} records from {os.path.basename(cache)}",
              flush=True)
        return recs
    ids, total = esearch_ids(query)
    print(f"    NCBI holds {total} records; fetching ALL of them "
          f"(cached for re-runs)", flush=True)
    recs = []
    BIG = 500
    for i in range(0, len(ids), BIG):
        chunk = ids[i:i + BIG]
        for attempt in range(4):
            try:
                handle = Entrez.efetch(db="protein", id=",".join(chunk),
                                       rettype="fasta", retmode="text")
                recs.extend(list(SeqIO.parse(handle, "fasta")))
                handle.close()
                break
            except Exception as exc:
                print(f"      [retry {attempt+1}/4] {exc}", flush=True)
                time.sleep(4 * (attempt + 1))
        if (i // BIG) % 4 == 0:
            print(f"      {min(i+BIG,len(ids))}/{len(ids)} fetched", flush=True)
        time.sleep(0.34)
    tmp = cache + ".part"
    with open(tmp, "w") as fh:
        for r in recs:
            fh.write(f">{r.id} {r.description}\n{str(r.seq)}\n")
    os.replace(tmp, cache)
    print(f"    cached {len(recs)} records -> {os.path.basename(cache)}", flush=True)
    return recs


def annotate_selected(records):
    """
    Second pass: re-fetch ONLY the selected records as GenBank so /gene and the
    source features are available for the independent gene check and the clade
    audit. Keyed by accession so the annotation is matched to the right record.
    """
    accs = [r.id for r in records]
    gb = efetch_records(accs, rettype="gb")
    return {g.id: g for g in gb}


# =============================================================================
# HIV SUBUNIT SLICING
# =============================================================================
def slice_hiv_subunits(seq, parent, ref_seq, audit):
    """
    Returns {subunit: (start, end, seq)} 1-based inclusive, or {} if the record
    cannot be resolved.

    PRIMARY method is the alignment: map the reference's UniProt CHAIN
    coordinates onto this record. Fixed offsets are NOT reused -- the previous
    run pinned them to YES72107.1 (856 aa) and YES72110.1 (1437 aa), and the
    records NCBI returns today are 860-aa Env and 493-aa Gag, so those numbers
    refer to nothing.

    CORROBORATION is the cleavage landmark. Where both resolve, they must agree
    within LANDMARK_TOLERANCE or the record is reported as a disagreement rather
    than one method being silently preferred.
    """
    LANDMARK_TOLERANCE = 3
    chains = common.HIV_SUBUNIT_CHAINS[parent]
    landmark_parent = "Env" if parent == "HIV_Env" else "Gag"
    landmarks = common.subunit_boundaries_by_landmark(seq, landmark_parent) or {}

    out = {}
    for subunit, (ref_s, ref_e) in chains.items():
        mapped = common.map_reference_span(seq, ref_seq, ref_s, ref_e)
        if mapped is None:
            audit.append(f"{subunit}:ALIGN_FAILED")
            continue
        qs, qe = mapped

        lm = landmarks.get(subunit)
        if lm:
            lm_s, lm_e = lm
            if lm_s is not None and abs(lm_s - qs) > LANDMARK_TOLERANCE:
                audit.append(f"{subunit}:START_DISAGREE(align={qs},landmark={lm_s})")
            if lm_e is not None and abs(lm_e - qe) > LANDMARK_TOLERANCE:
                audit.append(f"{subunit}:END_DISAGREE(align={qe},landmark={lm_e})")

        sub_seq = seq[qs - 1:qe]
        lo, hi = common.HIV_SUBUNIT_LENGTH_BAND[subunit]
        if not (lo <= len(sub_seq) <= hi):
            audit.append(f"{subunit}:LENGTH_OUT_OF_BAND({len(sub_seq)} not {lo}-{hi})")
            continue
        if common.has_nonstandard(sub_seq):
            audit.append(f"{subunit}:NONSTANDARD")
            continue
        out[subunit] = (qs, qe, sub_seq)
    return out


# =============================================================================
# MAIN
# =============================================================================
def run_retrieval():
    start_time = time.time()
    phase1a_path = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1A")
    ref_dir = os.path.join(phase1a_path, "Reference_Antigens")
    os.makedirs(ref_dir, exist_ok=True)

    common.print_banner("PHASE 1A: SEQUENCE RETRIEVAL WITH ANTIGEN IDENTITY GATE", 88)
    print(f"[INFO] Start                : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Output               : {phase1a_path}")
    print(f"[INFO] Variants per target  : {VARIANTS_PER_TARGET} (diversity-selected)")
    print(f"[INFO] Candidates gated     : up to {CANDIDATE_POOL_SIZE} per source query")
    print("-" * 88)

    # --- purge previous FASTAs so a re-run cannot mix two retrievals ---
    purged = 0
    if os.path.isdir(phase1a_path):
        for f in os.listdir(phase1a_path):
            if f.endswith(".fasta"):
                os.remove(os.path.join(phase1a_path, f))
                purged += 1
    print(f"[PROCESS] Purged {purged} prior FASTA files.")

    # --- reference panel, fetched fresh ---
    print("\n[STEP 1/4] Fetching UniProt reference panel...")
    needed = set()
    for spec in common.ANTIGEN_REFERENCES.values():
        needed.add(spec["intended"])
        needed.update(spec["decoys"].keys())
    for acc in sorted(needed):
        seq = common.fetch_uniprot_fasta(acc, ref_dir)
        print(f"    {acc:<12} {len(seq):>5} aa")

    audit_rows = []
    clade_rows = []
    written = defaultdict(int)

    def gate_pool(source_label, query, is_mpox, exhaustive=False):
        """Fetch -> identity-gate -> return accepted [(record, seq)]."""
        print(f"\n[STEP 2/4] {source_label}")
        print(f"    query : {query}")
        if exhaustive:
            cache_dir = os.path.join(phase1a_path, "_fetch_cache")
            records = fetch_all_cached(query, source_label, cache_dir)
        else:
            ids, total = esearch_ids(query)
            sampled = stride_sample(ids, CANDIDATE_POOL_SIZE)
            print(f"    NCBI holds {total} records; downloading {len(sampled)} "
                  f"stride-sampled across the full result set")
            records = efetch_records(sampled, rettype="fasta")
        print(f"    fetched {len(records)} records; applying identity gate...")

        expected_gene = None   # /gene is verified on the SELECTED records only

        # DEDUPLICATE BEFORE GATING. Identity is a property of the SEQUENCE
        # alone, so aligning 8,700 records that collapse to 23 distinct
        # sequences repeats the same alignment hundreds of times. Gating the
        # distinct set instead is ~300x faster and decides exactly the same
        # thing. Record multiplicity is retained so diversity_select can still
        # seed on the most frequently deposited sequence.
        from collections import Counter as _Counter
        seq_of = {}
        counts = _Counter()
        for rec in records:
            sq = common.sanitise_sequence(str(rec.seq))
            counts[sq] += 1
            seq_of.setdefault(sq, rec)
        distinct = list(counts)
        print(f"    {len(records)} records collapse to {len(distinct)} distinct "
              f"sequences; gating the distinct set", flush=True)

        accepted = []
        for _i, sq in enumerate(distinct, 1):
            rec = seq_of[sq]
            if _i % 25 == 0:
                print(f"      gated {_i}/{len(distinct)} distinct  "
                      f"accepted so far {len(accepted)}", flush=True)
            seq = sq
            gene = gene_qualifier(rec)
            ok, reason, detail = common.identity_verdict(seq, source_label, ref_dir)

            # /gene is an INDEPENDENT check, not a substitute for identity.
            # Where NCBI records it, it must agree; where it is absent the
            # alignment still decides, so a missing qualifier cannot let a
            # wrong protein through.
            gene_status = "ABSENT"
            if expected_gene and gene:
                gene_status = "MATCH" if gene.upper() == expected_gene.upper() else f"MISMATCH({gene})"
                if ok and gene_status.startswith("MISMATCH"):
                    ok, reason = False, f"GENE_QUALIFIER_{gene_status}"

            audit_rows.append({
                "Source": source_label, "Accession": rec.id, "Length": len(seq),
                "N_Records_With_This_Sequence": counts[sq],
                "Gene_Qualifier": gene or "", "Gene_Check": gene_status,
                "Best_Match": detail.get("best_acc", ""),
                "Best_Identity_Pct": detail.get("best_pct", ""),
                "Intended_Ref": detail.get("intended_acc", ""),
                "Intended_Identity_Pct": detail.get("intended_pct", ""),
                "Panel_Scores": detail.get("panel", ""),
                "Accepted": ok, "Reason": reason,
                "Description": rec.description[:120],
            })
            if ok:
                accepted.append((rec, seq, counts[sq]))

        rej = len(distinct) - len(accepted)
        print(f"    ACCEPTED {len(accepted)} / {len(distinct)} distinct "
              f"(rejected {rej})")
        if rej:
            from collections import Counter
            reasons = Counter(r["Reason"].split("(")[0] for r in audit_rows
                              if r["Source"] == source_label and r["Accepted"] is False)
            for why, n in reasons.most_common():
                print(f"        rejected {n:>4}  {why}")
        return accepted

    def emit(target, picked):
        """Write the selected variants as one FASTA each."""
        for i, (rec, seq) in enumerate(picked, start=1):
            fname = f"{target}_Var_{i:02d}_{rec.id}.fasta"
            with open(os.path.join(phase1a_path, fname), "w") as fh:
                fh.write(f">{target}_Var_{i:02d}_{rec.id} {rec.description[:100]}\n")
                for j in range(0, len(seq), 60):
                    fh.write(seq[j:j + 60] + "\n")
            written[target] += 1

    # ---------------- Mpox ----------------
    for target, opg in MPOX_SOURCES.items():
        query = f'{opg}[Gene Name] AND "Monkeypox virus"[Organism]'
        accepted = gate_pool(target, query, is_mpox=True, exhaustive=True)
        if not accepted:
            print(f"[ERROR] {target}: zero records survived the identity gate.")
            continue
        # Pass the FULL list including duplicates: diversity_select seeds on the
        # most frequently deposited exact sequence (the consensus strain, which
        # becomes Var_01). Deduplicating first would flatten every count to 1
        # and make that seed arbitrary.
        all_seqs = []
        by_seq = {}
        for rec, seq, mult in accepted:
            all_seqs.extend([seq] * mult)   # replicate so frequency is preserved
            by_seq.setdefault(seq, rec)
        # Over-select so a /gene rejection in the annotation pass below can be
        # backfilled without leaving the target short of its 30.
        overselect = min(VARIANTS_PER_TARGET + 12, len(by_seq))
        chosen_seqs = common.diversity_select(all_seqs, overselect)
        shortlist = [(by_seq[s], s) for s in chosen_seqs]
        print(f"    distinct sequences {len(by_seq)} of {len(all_seqs)} accepted; "
              f"diversity-selected {len(shortlist)} (over-selected for /gene check)")

        # SECOND, INDEPENDENT CHECK on the shortlist only: NCBI's own /gene
        # qualifier must agree with the OPG number we asked for. The alignment
        # already decided identity; this catches the separate failure mode where
        # a record is the right protein family but filed under another gene.
        print(f"    annotating {len(shortlist)} selected records (GenBank) for "
              f"/gene verification and clade tagging...", flush=True)
        ann = annotate_selected([r for r, _s in shortlist])
        picked, gene_rejected = [], 0
        for rec, seq in shortlist:
            if len(picked) >= VARIANTS_PER_TARGET:
                break
            g = ann.get(rec.id)
            gene = gene_qualifier(g) if g is not None else ""
            if gene and gene.upper() != opg.upper():
                gene_rejected += 1
                audit_rows.append({
                    "Source": target, "Accession": rec.id, "Length": len(seq),
                    "N_Records_With_This_Sequence": "",
                    "Gene_Qualifier": gene, "Gene_Check": f"MISMATCH({gene})",
                    "Best_Match": "", "Best_Identity_Pct": "",
                    "Intended_Ref": "", "Intended_Identity_Pct": "",
                    "Panel_Scores": "", "Accepted": False,
                    "Reason": f"GENE_QUALIFIER_MISMATCH(expected {opg}, got {gene})",
                    "Description": rec.description[:120]})
                continue
            picked.append((rec, seq))
            tag, confirmed = classify_clade(g) if g is not None else ("Unconfirmed", False)
            clade_rows.append({"Target": target, "Accession": rec.id,
                               "Clade_Tag": tag, "Confirmed": confirmed,
                               "Gene_Qualifier": gene or "ABSENT"})
        print(f"    /gene check: {len(picked)} kept, {gene_rejected} rejected on "
              f"gene mismatch")
        emit(target, picked)

    # ---------------- HIV ----------------
    for parent, query in HIV_SOURCES.items():
        accepted = gate_pool(parent, query, is_mpox=False)
        if not accepted:
            print(f"[ERROR] {parent}: zero records survived the identity gate.")
            continue
        ref_seq = common.fetch_uniprot_fasta(
            common.ANTIGEN_REFERENCES[parent]["intended"], ref_dir)

        print(f"[STEP 3/4] {parent}: slicing to mature subunits (alignment + landmark)")
        per_sub = {s: {} for s in HIV_SUBUNITS_OF[parent]}
        per_sub_all = {s: [] for s in HIV_SUBUNITS_OF[parent]}
        notes = defaultdict(list)
        for rec, seq, mult in accepted:
            issues = []
            subs = slice_hiv_subunits(seq, parent, ref_seq, issues)
            for sub, (qs, qe, sseq) in subs.items():
                per_sub[sub].setdefault(sseq, (rec, qs, qe))
                per_sub_all[sub].extend([sseq] * mult)
            for msg in issues:
                notes[msg.split(":")[0]].append(msg)
        for sub in HIV_SUBUNITS_OF[parent]:
            n_issue = len(notes.get(sub, []))
            print(f"    {sub:<12} distinct subunit sequences: {len(per_sub[sub]):>4}"
                  f"   records with issues: {n_issue}")
            if n_issue:
                from collections import Counter
                for msg, c in Counter(m.split("(")[0] for m in notes[sub]).most_common(3):
                    print(f"        {c:>4}  {msg}")

        for sub in HIV_SUBUNITS_OF[parent]:
            pool = per_sub[sub]
            if not pool:
                print(f"[ERROR] {sub}: no subunit sequences resolved.")
                continue
            chosen = common.diversity_select(per_sub_all[sub], VARIANTS_PER_TARGET)
            picked = [(pool[s][0], s) for s in chosen]
            emit(sub, picked)
            for s in chosen:
                rec, qs, qe = pool[s]
                audit_rows.append({
                    "Source": sub, "Accession": rec.id, "Length": len(s),
                    "N_Records_With_This_Sequence": "",
                    "Gene_Qualifier": "", "Gene_Check": "N/A",
                    "Best_Match": f"sliced from {parent}",
                    "Best_Identity_Pct": "", "Intended_Ref": "",
                    "Intended_Identity_Pct": "",
                    "Panel_Scores": f"subunit span {qs}-{qe} on parent",
                    "Accepted": True, "Reason": "SUBUNIT_SLICE",
                    "Description": rec.description[:120],
                })

    # ---------------- outputs ----------------
    print("\n[STEP 4/4] Writing audit files...")
    audit_path = os.path.join(phase1a_path, "Phase1A_IdentityGate.csv")
    fields = ["Source", "Accession", "Length", "N_Records_With_This_Sequence",
              "Gene_Qualifier", "Gene_Check",
              "Best_Match", "Best_Identity_Pct", "Intended_Ref",
              "Intended_Identity_Pct", "Panel_Scores", "Accepted", "Reason",
              "Description"]
    with open(audit_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(audit_rows)
    print(f"    identity gate audit : {audit_path}  ({len(audit_rows)} rows)")

    if clade_rows:
        clade_path = os.path.join(phase1a_path, "Phase1A_Clade_Audit.csv")
        with open(clade_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["Target", "Accession", "Clade_Tag", "Confirmed", "Gene_Qualifier"])
            w.writeheader()
            w.writerows(clade_rows)
        n_conf = sum(1 for r in clade_rows if r["Confirmed"])
        print(f"    clade audit         : {clade_path}  "
              f"({n_conf}/{len(clade_rows)} clade-confirmed)")

    common.print_banner("RETRIEVAL COMPLETE", 88)
    total = 0
    for t in list(MPOX_SOURCES) + ["HIV_gp120", "HIV_gp41", "HIV_p17", "HIV_p24"]:
        print(f"    {t:<12} {written.get(t, 0):>3} variants written")
        total += written.get(t, 0)
    print(f"\n[SUCCESS] {total} sequences written across 7 targets "
          f"in {common.format_time(time.time() - start_time)}")
    print("[INFO] Every sequence above passed an identity check against a UniProt")
    print("       reference BEFORE being written. See Phase1A_IdentityGate.csv.")
    print("=" * 88 + "\n")


if __name__ == "__main__":
    run_retrieval()
