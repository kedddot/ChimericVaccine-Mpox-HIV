"""
phase5_common.py

Shared utilities for Section V (Pharmacokinetics & Dose-Response, Steps A-E).

Phase V never edits Phase 1/Phase 2/Phase 4 outputs (Rule 1) -- this module only
READS from those trees and re-exports a few of their helper functions so
Phase 5 steps do not reimplement logic that already exists (latest_file,
print_banner, etc). Everything Phase-5 specific (timestamps, the parameter
ledger, the HALT-FOR-OPUS writer) lives here so all five steps behave
identically -- mirrors Phase 4's phase4_common.py exactly.
"""

import os
import sys
import csv
from datetime import datetime


# =============================================================================
# MINIMAL BOOTSTRAP -- copied verbatim from phase4_common.py / any Phase 1
# step script, so all phases resolve the project root identically regardless
# of how deep a script sits.
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


PROJECT_ROOT = _bootstrap_find_research_root(__file__)

_PHASE1_COMMON_DIR = os.path.join(PROJECT_ROOT, "Phase 1", "_common")
if _PHASE1_COMMON_DIR not in sys.path:
    sys.path.insert(0, _PHASE1_COMMON_DIR)

import phase1_common as _p1  # noqa: E402  (import after sys.path fixup, by design)

# Re-exported as-is -- Rule 1 forbids editing Phase 1, importing is safe
# because every Phase 1 script is __main__-guarded.
print_banner = _p1.print_banner
format_time = _p1.format_time
latest_file = _p1.latest_file

CONSTRUCT_ID = "Vax_Final_4fce119e"


def timestamp():
    """Phase 5's own timestamp format: no dashes, since steps pick inputs by
    sorting filenames and '-' sorts before '0' (matches Phase 4's convention)."""
    return datetime.now().strftime("%Y%m%d_%H%M")


def step_output_dir(step_letter):
    """e.g. step_output_dir('A') -> Step_Outputs/Phase5/StepA"""
    d = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase5", f"Step{step_letter}")
    os.makedirs(d, exist_ok=True)
    return d


def run_log_path(step_name):
    d = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase5", "_run_logs")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{step_name}_{timestamp()}.log")


def phase4_ledger_path():
    """The Step 7 (Phase4G) numbers ledger -- read-only source for Step 5D's
    alpha parameterisation. Fixed path, not latest_file()'d, because Phase IV
    is closed and this is its one and only final ledger."""
    return os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase4", "StepG",
                         "Phase4G_NumbersLedger_20260921_2022.csv")


def read_phase4_ledger():
    path = phase4_ledger_path()
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {r["Key"]: r["Value"] for r in rows}


# =============================================================================
# NUMBERS LEDGER -- Phase IV style. Every number this phase reports is
# recorded as (key, value, source_file) so every figure in the final report
# traces back to a CSV. One ledger instance per step; steps hand their rows
# to Step E, which merges them into the master ledger (matches Phase 4G).
# =============================================================================
class Ledger:
    def __init__(self):
        self.rows = []

    def N(self, key, value, source_file, fmt=None):
        """Record a number and return it (optionally formatted) for printing.
        `value` is stored as given; `fmt` is a format spec applied only to the
        returned/printed string, never to what's stored."""
        self.rows.append({"Key": key, "Value": value, "Source_File": source_file})
        if fmt is not None and value is not None:
            try:
                return format(value, fmt)
            except (ValueError, TypeError):
                return str(value)
        return value

    def write(self, path):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["Key", "Value", "Source_File"])
            w.writeheader()
            w.writerows(self.rows)
        return path

    def extend(self, other_rows):
        self.rows.extend(other_rows)


# =============================================================================
# ESCALATION -- mirrors phase4_common.halt_for_opus exactly (idempotent per
# branch: a re-run does not re-open a branch whose HALT section already
# contains "RESOLVED").
# =============================================================================
def halt_for_opus(branch, trigger_number, one_liner, what_i_was_doing, exact_numbers, options):
    path = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase5", "HALT-FOR-OPUS.md")
    if os.path.isfile(path):
        with open(path) as f:
            head, *sections = f.read().split("## Branch: ")
        prior = [x for x in sections if x.startswith(branch + "\n")]
        if prior and "RESOLVED" in prior[0]:
            print(f"[INFO] HALT branch '{branch}' was already RESOLVED by a ruling -- a re-run does not re-open it.")
            return path
        kept = [x for x in sections if not x.startswith(branch + "\n")]
        with open(path, "w") as f:
            f.write(head + "".join("## Branch: " + x for x in kept))
    exists = os.path.isfile(path)
    with open(path, "a") as f:
        if not exists:
            f.write("# Section V -- Halted Branches (escalated to Opus)\n\n")
        f.write(f"## Branch: {branch}\n\n")
        f.write(f"**Trigger fired:** #{trigger_number} -- {one_liner}\n\n")
        f.write(f"**What I was doing:**\n\n{what_i_was_doing}\n\n")
        f.write(f"**Exact numbers / output:**\n\n{exact_numbers}\n\n")
        f.write("**Options I see:**\n\n")
        for i, opt in enumerate(options, 1):
            f.write(f"{i}. {opt}\n")
        f.write(f"\n_Halted: {datetime.now().isoformat(timespec='seconds')}_\n\n")
        f.write("-" * 90 + "\n\n")
    print(f"[ACTION REQUIRED] Branch '{branch}' halted -- trigger #{trigger_number}: {one_liner}")
    print(f"[ACTION REQUIRED] Details written to: {path}")
    return path


# =============================================================================
# RULE 1 -- verify Phase 1/2/4 are untouched. Every Phase 5 step calls this
# before and after doing anything, and refuses to proceed if it fails.
# =============================================================================
def verify_rule1(where):
    import hashlib
    baseline = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase4", "_rule1_baseline_sha256.txt")
    changed, missing, total = [], [], 0
    with open(baseline) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or "  " not in line:
                continue
            sha, path = line.split("  ", 1)
            total += 1
            full = os.path.join(PROJECT_ROOT, path)
            if not os.path.isfile(full):
                missing.append(path)
                continue
            h = hashlib.sha256()
            with open(full, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            if h.hexdigest() != sha:
                changed.append(path)
    ok = not changed and not missing
    print(f"[RULE 1] {where}: baseline {total} files, changed {len(changed)}, missing {len(missing)} "
          f"-> {'CLEAN' if ok else 'DIRTY'}")
    if not ok:
        print(f"[ERROR] Rule 1 violated ({where}). Changed: {changed[:10]} Missing: {missing[:10]}")
    return ok
