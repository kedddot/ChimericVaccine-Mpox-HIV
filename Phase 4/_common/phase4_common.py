"""
phase4_common.py

Shared utilities for Phase IV (IEDB-based immunogenicity prediction, Steps A-G).

Phase IV never edits Phase 1/Phase 2 outputs (Rule 1) -- this module only
READS from those trees and re-exports a few of their helper functions so
Phase 4 steps do not reimplement logic that Phase 1 already got right
(latest_file, print_banner, human_self_homology, etc). Everything Phase-4
specific (timestamps, API response caching, the HALT-FOR-OPUS writer) lives
here so all seven steps behave identically.
"""

import os
import sys
import json
import hashlib
from datetime import datetime


# =============================================================================
# MINIMAL BOOTSTRAP -- locates the shared phase1_common module and the
# Research root. Copied verbatim from the Phase 1 bootstrap pattern (e.g.
# Phase 1/STEP A/Phase1A_provenanceCorrection.py) so all phases resolve the
# project root identically regardless of how deep a script sits.
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
human_self_homology = _p1.human_self_homology
allergen_homology = _p1.allergen_homology
global_identity = _p1.global_identity

CONSTRUCT_ID = "Vax_Final_4fce119e"
EXPECTED_TOTAL_EPITOPES = 31
EXPECTED_BY_CLASS = {"MHC-I": 10, "MHC-II": 11, "B-cell": 10}
EXPECTED_HIV_COUNT = 14
EXPECTED_MPOX_COUNT = 17
EXPECTED_SOURCE_ANTIGENS = 7


def pathogen_of(target):
    """'HIV_gp120' -> 'HIV', 'Mpox_A35R' -> 'Mpox'. Raises on anything else --
    a Phase IV dossier should never contain a target this can't classify."""
    if target.startswith("HIV"):
        return "HIV"
    if target.startswith("Mpox"):
        return "Mpox"
    raise ValueError(f"Unrecognised target prefix (not HIV/Mpox): {target!r}")


def timestamp():
    """Phase 4's own timestamp format: no dashes, since steps pick inputs by
    sorting filenames and '-' sorts before '0' (see Section E of the task)."""
    return datetime.now().strftime("%Y%m%d_%H%M")


def step_output_dir(step_letter):
    """e.g. step_output_dir('A') -> Step_Outputs/Phase4/StepA"""
    d = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase4", f"Step{step_letter}")
    os.makedirs(d, exist_ok=True)
    return d


def run_log_path(step_name):
    d = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase4", "_run_logs")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{step_name}_{timestamp()}.log")


# =============================================================================
# CONTENT-HASH RESPONSE CACHE -- mirrors the pattern in Phase 1Db / Phase 2A.
# Re-running a step must make zero network calls (Section G final check), so
# every IEDB/BigMHC/PRIME call in Steps B-D goes through this.
# =============================================================================
def cache_dir_for(step_letter):
    d = os.path.join(step_output_dir(step_letter), "_tool_runs")
    os.makedirs(d, exist_ok=True)
    return d


def content_key(payload):
    """Stable content hash for a request payload (dict, list, or str)."""
    if not isinstance(payload, str):
        payload = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def cache_get(step_letter, tool, payload):
    key = content_key(payload)
    path = os.path.join(cache_dir_for(step_letter), f"{tool}__{key}.json")
    if os.path.isfile(path):
        with open(path) as f:
            return json.load(f)
    return None


def cache_put(step_letter, tool, payload, response_obj):
    key = content_key(payload)
    path = os.path.join(cache_dir_for(step_letter), f"{tool}__{key}.json")
    with open(path, "w") as f:
        json.dump(response_obj, f, indent=2, default=str)
    return path


# =============================================================================
# IEDB POST WITH DISK CACHE, IPv4 FORCING, AND RETRY -- mirrors the pattern in
# Phase1Db_filtration.py's cached_post() and Phase2A_secondaryScreening.py's
# run_iedb_mhc_binding_batch() (IPv4 forcing avoids the ~30s IPv6 dead-end on
# this network; disk caching keyed by request content means a rerun makes
# zero network calls). A 403 from IEDB has been observed here to be
# transient throttling from calling too fast, not a rejected allele/method --
# confirmed live by immediately retrying the exact same request after a
# short delay and getting 200 -- so 403/5xx get retried with backoff before
# being treated as a real failure.
# =============================================================================
_IPV4_FORCED = False


def _force_ipv4():
    global _IPV4_FORCED
    if _IPV4_FORCED:
        return
    import socket
    import urllib3.util.connection as urllib3_cn
    urllib3_cn.allowed_gai_family = lambda: socket.AF_INET
    _IPV4_FORCED = True


def iedb_post_cached(step_letter, tool_name, url, payload, cache_key, max_retries=3):
    """
    Returns (text, was_cached). text is None if every attempt failed --
    callers must treat that as "no prediction for this allele/batch",
    never as a peptide-level UNRESOLVED/fail.
    """
    import time
    import requests

    _force_ipv4()

    digest = content_key(cache_key)
    path = os.path.join(cache_dir_for(step_letter), f"{tool_name}__{digest}.txt")
    if os.path.isfile(path):
        with open(path) as f:
            return f.read(), True

    text = None
    for attempt in range(max_retries):
        try:
            resp = requests.post(url, data=payload, timeout=180)
            if resp.status_code == 200:
                text = resp.text
                break
            print(f"[WARNING] IEDB {tool_name} HTTP {resp.status_code} "
                  f"(attempt {attempt + 1}/{max_retries}) -- allele={payload.get('allele')}")
        except Exception as e:
            print(f"[WARNING] IEDB {tool_name} request error (attempt {attempt + 1}/{max_retries}): {e}")
        if attempt < max_retries - 1:
            time.sleep(3.0 * (attempt + 1))

    if text is not None:
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            f.write(text)
        os.replace(tmp, path)
    else:
        print(f"[ERROR] IEDB {tool_name} failed after {max_retries} attempts -- "
              f"allele={payload.get('allele')} -- treated as no-prediction, NOT cached.")
    return text, False


# =============================================================================
# ESCALATION -- Section H. Each halting branch gets its own section appended
# to one shared file; a halt on one branch never touches another branch's
# section, matching "pause only that branch, keep working on anything
# independent."
# =============================================================================
def halt_for_opus(branch, trigger_number, one_liner, what_i_was_doing, exact_numbers, options):
    path = os.path.join(PROJECT_ROOT, "Step_Outputs", "Phase4", "HALT-FOR-OPUS.md")
    if os.path.isfile(path):
        # Idempotent per branch: a re-run replaces that branch's section instead of duplicating it.
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
            f.write("# Phase IV -- Halted Branches (escalated to Opus)\n\n")
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
