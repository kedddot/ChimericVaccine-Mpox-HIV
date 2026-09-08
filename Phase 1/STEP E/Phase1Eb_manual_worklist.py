"""
Phase 1Eb -- MANUAL WEB-GATE HELPER
===================================
Generates the two worklists a human fills in at the AllerTOP / AllergenFP web
forms, and ingests them back into Manual_Allergenicity_Results.csv.

WHY THIS EXISTS. The alternative was hand-editing a 184-row CSV in a
spreadsheet while alt-tabbing to a browser. That is the single most
error-prone action in this whole pipeline: one skipped row silently shifts
every verdict below it onto the wrong peptide, the file stays perfectly
well-formed, and Phase 1Eb would filter on it without complaint. That is the
exact failure class this rebuild exists to remove.

SO THE WORKLIST PRINTS THE SEQUENCE ON EVERY LINE and the merge is keyed BY
SEQUENCE, never by line number. An off-by-one cannot corrupt the result; it can
only produce a mismatch, which is reported and refused.

USAGE
    python Phase1Eb_manual_worklist.py            # write the blank worklists
    python Phase1Eb_manual_worklist.py --ingest   # read them back, validate,
                                                  # write the results CSV

The web tools are behind an interactive Cloudflare challenge (verified
2026-08-31: 'cf-mitigated: challenge', HTTP 403 to any scripted request).
That is deliberate anti-bot protection on the authors' side and is NOT
circumvented here -- a person using a browser is the supported route, and the
same position was taken for SEMA in deviation #20.
"""

import os, sys, csv, re

# A data line, and ONLY a data line, looks like:
#     12  PEPTIDESEQUENCE   (16 aa, Pep_11)  = N
# Parsing by regex rather than by "does the line contain '='" is deliberate:
# the instructions above the table talk about the "=" sign, and a looser rule
# tried to read those sentences as peptide answers. Caught in testing, and the
# lesson is the general one -- never let free text and data share a namespace.
_DATA_LINE = re.compile(
    r"^\s*(?P<idx>\d+)\s+(?P<pep>[A-Za-z]+)\s+\([^)]*\)\s*=\s*(?P<ans>.*?)\s*$")

# Everything above this marker is instructions and is never parsed.
_SENTINEL = "--- ANSWERS BELOW THIS LINE, ONE PER PEPTIDE ---"

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
EB_DIR = os.path.join(_PROJECT_ROOT, "Step_Outputs", "Phase1", "Phase1E", "Phase1Eb")

ALLERTOP_WORKLIST = os.path.join(EB_DIR, "WORKLIST_1_AllerTOP.txt")
ALLERGENFP_WORKLIST = os.path.join(EB_DIR, "WORKLIST_2_AllergenFP.txt")
MANUAL_CSV = os.path.join(EB_DIR, "Manual_Allergenicity_Results.csv")

# What a person may type. Deliberately permissive about FORM (case, the full
# word, the letter) and completely strict about MEANING -- anything not on this
# list is refused rather than guessed at.
VERDICTS = {
    "A": "ALLERGEN", "ALLERGEN": "ALLERGEN", "Y": "ALLERGEN", "YES": "ALLERGEN",
    "N": "NON-ALLERGEN", "NON-ALLERGEN": "NON-ALLERGEN", "NONALLERGEN": "NON-ALLERGEN",
    "NO": "NON-ALLERGEN", "PROBABLE NON-ALLERGEN": "NON-ALLERGEN",
    "PROBABLE ALLERGEN": "ALLERGEN",
    # AllerTOP v2.1's result page phrases its answer as the sequence
    # "IS a probable allergen" / "is NOT a probable allergen", and the operator
    # transcribes that wording directly. Added after it appeared in the real
    # worklist rather than guessed at in advance.
    "IS": "ALLERGEN", "NOT": "NON-ALLERGEN",
    "IS ALLERGEN": "ALLERGEN", "NOT ALLERGEN": "NON-ALLERGEN",
    "IS A PROBABLE ALLERGEN": "ALLERGEN",
    "IS NOT A PROBABLE ALLERGEN": "NON-ALLERGEN",
}

HEADER = """\
================================================================================
{title}
================================================================================
TOOL : {url}
JOB  : {n} peptides{extra}

HOW TO FILL THIS IN
  1. Open the tool in a normal browser tab.
  2. Paste the peptide (the middle column) into the sequence box and submit.
  3. Type the answer after the "=" on that line:
         A  = allergen        (also: IS, ALLERGEN, Y, YES,
                               "is a probable allergen")
         N  = non-allergen    (also: NOT, NON-ALLERGEN, NO,
                               "is not a probable allergen")
  4. Save this file. Leave anything you have not done yet blank.
  5. Run:  python Phase1Eb_manual_worklist.py --ingest

YOU CAN STOP AND RESUME AT ANY POINT. Blank lines are simply "not done yet".
Re-running --ingest merges whatever is filled in so far and tells you what is
left. Nothing is lost by quitting halfway.

DO NOT reorder, delete or retype the peptide column. Answers are matched to
peptides BY SEQUENCE, so a mistyped sequence is reported as an error rather
than silently attached to the wrong peptide -- but it still has to be fixed.
================================================================================
{sentinel}
"""


def _load_peptides():
    """Reads the exact peptide set 1Eb is waiting on, from the query FASTAs it
    already wrote -- not from a re-derivation, so the worklist cannot drift
    from what was actually submitted."""
    def readfa(path):
        out, name = [], None
        if not os.path.isfile(path):
            print(f"[ERROR] Missing {os.path.basename(path)}. Run Phase1Eb_allergenicity.py first.")
            sys.exit(1)
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif line:
                out.append((name, line))
        return out
    return (readfa(os.path.join(EB_DIR, "Allergenicity_Query.fasta")),
            readfa(os.path.join(EB_DIR, "AllergenFP_Query_16aa.fasta")))


def _parse_worklist(path):
    """Returns {peptide: verdict} plus a list of problems. Tolerates partial
    completion; refuses to interpret anything it does not recognise."""
    answers, problems, blank = {}, [], 0
    if not os.path.isfile(path):
        return answers, [f"{os.path.basename(path)} not found"], 0

    lines = open(path).read().splitlines()
    try:
        start = next(i for i, l in enumerate(lines) if _SENTINEL in l) + 1
    except StopIteration:
        return answers, [f"{os.path.basename(path)}: the "
                         f"'ANSWERS BELOW THIS LINE' marker is missing -- the file "
                         f"was edited or truncated. Regenerate it by running this "
                         f"script with no arguments."], 0

    seen_idx = set()
    for offset, raw in enumerate(lines[start:]):
        lineno = start + offset + 1
        if not raw.strip():
            continue
        m = _DATA_LINE.match(raw)
        if not m:
            problems.append(f"line {lineno}: not a readable worklist row -- "
                            f"{raw.strip()[:60]!r}")
            continue
        idx, pep, ans = m.group("idx"), m.group("pep").upper(), m.group("ans").upper()
        if idx in seen_idx:
            problems.append(f"line {lineno}: row number {idx} appears twice")
        seen_idx.add(idx)
        if not ans:
            blank += 1
            continue
        if ans not in VERDICTS:
            problems.append(f"line {lineno}: cannot interpret answer "
                            f"{m.group('ans')!r} for {pep} -- use A or N")
            continue
        if pep in answers and answers[pep] != VERDICTS[ans]:
            problems.append(f"line {lineno}: {pep} answered twice, and the two "
                            f"answers disagree ({answers[pep]} vs {VERDICTS[ans]})")
            continue
        answers[pep] = VERDICTS[ans]
    return answers, problems, blank


def write_worklists():
    at, fp = _load_peptides()
    with open(ALLERTOP_WORKLIST, "w") as f:
        f.write(HEADER.format(
            title="WORKLIST 1 of 2  --  AllerTOP v2.1  --  ALL 184 PEPTIDES",
            url="https://ddg-pharmfac.net/AllerTOP/",
            n=len(at),
            extra=", 9 to 16 amino acids",
            sentinel=_SENTINEL))
        for i, (name, pep) in enumerate(at, 1):
            f.write(f"{i:3d}  {pep:<18s}  ({len(pep):2d} aa, {name})  = \n")
    with open(ALLERGENFP_WORKLIST, "w") as f:
        f.write(HEADER.format(
            title="WORKLIST 2 of 2  --  AllergenFP v1.0  --  THE 57 16-MERS ONLY",
            url="https://ddg-pharmfac.net/AllergenFP/",
            n=len(fp),
            extra=", all exactly 16 aa\n"
                  "       The other 127 peptides are NOT in this list: AllergenFP\n"
                  "       cannot compute its fingerprint below 16 aa (deviation #6).\n"
                  "       They are already recorded as N/A and need no submission.",
            sentinel=_SENTINEL))
        for i, (name, pep) in enumerate(fp, 1):
            f.write(f"{i:3d}  {pep:<18s}  (16 aa, {name})  = \n")
    print(f"[OK] {ALLERTOP_WORKLIST}   {len(at)} peptides")
    print(f"[OK] {ALLERGENFP_WORKLIST}   {len(fp)} peptides")
    print("\nFill them in a browser, then re-run with --ingest.")


def ingest():
    at, fp = _load_peptides()
    at_set = {p for _, p in at}
    fp_set = {p for _, p in fp}

    at_ans, at_prob, at_blank = _parse_worklist(ALLERTOP_WORKLIST)
    fp_ans, fp_prob, fp_blank = _parse_worklist(ALLERGENFP_WORKLIST)

    # A sequence in the worklist that is NOT a real query peptide means the
    # peptide column was edited. Refuse rather than drop it quietly.
    for label, ans, valid in (("AllerTOP", at_ans, at_set), ("AllergenFP", fp_ans, fp_set)):
        for pep in ans:
            if pep not in valid:
                (at_prob if label == "AllerTOP" else fp_prob).append(
                    f"{pep} is not one of the {label} query peptides -- the "
                    f"peptide column was altered")

    problems = at_prob + fp_prob

    # A DELETED ROW cannot corrupt anything -- the merge is keyed by sequence,
    # so every answer still lands on its own peptide. But the deleted peptide
    # would then never be asked about, and the run would look merely
    # "incomplete" rather than "edited". Tested and confirmed: deleting a row
    # produced no misattribution and no warning. So the row count is checked
    # against the query set and a missing peptide is named outright.
    for label, path, expected in (("AllerTOP", ALLERTOP_WORKLIST, at),
                                  ("AllergenFP", ALLERGENFP_WORKLIST, fp)):
        if not os.path.isfile(path):
            continue
        present = set()
        lines = open(path).read().splitlines()
        try:
            start = next(i for i, l in enumerate(lines) if _SENTINEL in l) + 1
        except StopIteration:
            continue
        for raw in lines[start:]:
            m = _DATA_LINE.match(raw)
            if m:
                present.add(m.group("pep").upper())
        missing = [p for _, p in expected if p not in present]
        if missing:
            problems.append(
                f"{os.path.basename(path)}: {len(missing)} peptide row(s) are "
                f"GONE from the file, not merely unanswered -- the worklist was "
                f"edited. First few: {', '.join(missing[:5])}. Regenerate the "
                f"worklist (run with no arguments) and re-enter, or paste the "
                f"missing rows back.")

    if problems:
        print("\n[REFUSED] The worklists could not be read cleanly. Nothing was written.\n")
        for p in problems[:40]:
            print("   -", p)
        if len(problems) > 40:
            print(f"   ... and {len(problems)-40} more")
        print("\nFix the lines above and re-run --ingest.")
        sys.exit(1)

    rows = []
    for _, pep in at:
        a = at_ans.get(pep, "")
        if len(pep) < 16:
            fpv = "N/A"           # tool genuinely cannot score it -- a resolved state
        else:
            fpv = fp_ans.get(pep, "")
        rows.append([pep, a, fpv])

    with open(MANUAL_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Peptide", "AllerTOP_Result", "AllergenFP_Result"])
        w.writerows(rows)

    done_at = sum(1 for r in rows if r[1])
    done_fp = sum(1 for r in rows if len(r[0]) >= 16 and r[2] in ("ALLERGEN", "NON-ALLERGEN"))
    ready = sum(1 for r in rows if r[1] and r[2])

    print(f"\n[OK] Wrote {os.path.basename(MANUAL_CSV)}  ({len(rows)} peptides)")
    print(f"     AllerTOP   {done_at:3d}/{len(at):3d} answered   ({len(at)-done_at} left)")
    print(f"     AllergenFP {done_fp:3d}/{len(fp):3d} answered   ({len(fp)-done_fp} left)")
    print(f"     fully resolved and ready for the filtering matrix: {ready}/{len(rows)}")
    if ready == len(rows):
        print("\n[READY] All peptides answered. Run Phase1Eb_allergenicity.py to filter.")
    else:
        print(f"\n[PARTIAL] {len(rows)-ready} peptide(s) still unanswered. Phase 1Eb will")
        print("          mark those UNRESOLVED and route them to Needs_Review rather")
        print("          than guess. Fill in more and re-run --ingest any time.")


def write_chunks(size=10):
    """Emit the queries in pasteable blocks.

    If the web form accepts several sequences at once, a block goes in as one
    submission and 184 becomes ~19. If it only takes one, the blocks are still
    a convenient place to work from. Written both ways because which one is
    true can only be found out in a browser."""
    at, fp = _load_peptides()
    for tag, recs in (("AllerTOP", at), ("AllergenFP", fp)):
        path = os.path.join(EB_DIR, f"CHUNKS_{tag}_by{size}.txt")
        with open(path, "w") as f:
            f.write(f"{tag} -- {len(recs)} sequences in blocks of {size}\n")
            f.write("Paste ONE BLOCK at a time. If the form rejects a whole\n"
                    "block, it only takes one sequence at a time -- use the\n"
                    "WORKLIST file instead.\n")
            f.write("Record answers in the WORKLIST file either way; that is\n"
                    "what --ingest reads.\n\n")
            for start in range(0, len(recs), size):
                block = recs[start:start + size]
                f.write("=" * 60 + "\n")
                f.write(f"BLOCK {start//size + 1}  (worklist rows "
                        f"{start+1}-{start+len(block)})\n")
                f.write("-- plain, for a sequence box --\n")
                for _, pep in block:
                    f.write(pep + "\n")
                f.write("-- FASTA, for a file/upload box --\n")
                for name, pep in block:
                    f.write(f">{name}\n{pep}\n")
                f.write("\n")
        print(f"[OK] {os.path.basename(path)}  "
              f"{len(recs)} sequences in {(len(recs)+size-1)//size} blocks")


if __name__ == "__main__":
    if "--chunks" in sys.argv:
        n = 10
        for a in sys.argv:
            if a.startswith("--size="):
                n = int(a.split("=")[1])
        write_chunks(n)
    elif "--ingest" in sys.argv:
        ingest()
    else:
        write_worklists()
