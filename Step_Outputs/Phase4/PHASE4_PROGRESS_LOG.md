# Phase IV Progress Log

**Purpose:** running record of what's been done, decided, and found, so a new
session (or a new person) can resume without re-deriving context. Update this
file at the end of every step — don't wait until Phase 4G.

Source instructions: `$R/PHASE_IV_SONNET_PROMPT.md` (the task brief). The 3
rules, the 8 pre-decided disagreements (Section C), and the escalation
triggers (Section H) are not re-litigated here — this log tracks *execution*
against that brief.

`$R` = `/Volumes/Extended SSD/ChimericVaccine/Research`

---

## Status at a glance

| Step | Script | Status | Output |
|---|---|---|---|
| 0 | toolchain fix | ✅ DONE | `Step_Outputs/Phase4/_env_check.md` |
| 1 (Phase4A) | epitope dossier | ✅ DONE | `StepA/Phase4A_EpitopeDossier_20260920_1855.csv` |
| 2 (Phase4B) | binding affinity | ✅ DONE (revised 19:44, then corrected 19:53 after a caught parsing error) | `StepB/Phase4B_BindingAffinity_20260920_1954.csv` |
| 3 (Phase4C) | TCR recognition | ✅ DONE | `StepC/Phase4C_TCRRecognition_MHCI_20260921_1704.csv`, `StepC/Phase4C_CrossReactivity_AllConstruct_20260921_1704.csv` |
| 4 (Phase4D) | proteasomal processing | ✅ DONE (trigger 3 resolved by Opus ruling) | `StepD/Phase4D_*_20260921_1844.csv` |
| 5 (Phase4E) | population coverage + junction weighting | ✅ DONE — RE-ISSUED under Opus's corrected definition; escalation resolved | `StepE/Phase4E_*_20260921_1934.csv` |
| 6 (Phase4F) | B-cell + conservancy | ✅ DONE (all gates passed) | `StepF/Phase4F_*_20260921_1933.csv` |
| 7 (Phase4G) | integrated report | ✅ DONE (2026-09-21; independently audited 2026-09-24) | `StepG/PHASE_IV_REPORT.md`, dossier, corrections, ledger |

**No open escalation.** HALT-FOR-OPUS.md holds three RESOLVED branches (Step 4 trigger 3; Step 5 contested-coverage; superseded figures there must not be reported). The Phase I/II-only Phase IV branch is complete; later Phase III integration is a separate task described in `Research/Handoff.md`.

**Rule 1 (never touch Phase 1/Phase 2) verified clean again on 2026-09-24**
(diffed all 1,810 files against the user-provided durable baseline
`_rule1_baseline_20260920_1850.tsv` — zero mtime changes, zero added/removed
files). Re-check this diff again at the end of Phase IV (Section G checklist).

---

## Step 0 — Toolchain fix (DONE)

**Problem:** conda envs + Phenix were moved into `Research (ARCHIVE)/`,
breaking every script's shebang.

**Fix:** two symlinks (no conda rebuild needed):
```
/Volumes/Extended SSD/MpoxHIV_deps    -> /Volumes/Extended SSD/Research (ARCHIVE)/MpoxHIV_deps
/Volumes/Extended SSD/phenix-2.2-6143 -> /Volumes/Extended SSD/Research (ARCHIVE)/phenix-2.2-6143
```
`~/mpoxhiv_ssd -> /Volumes/Extended SSD` already existed. Source
`~/mpoxhiv_env.sh` after creating the symlinks, before running anything.

**4 checks — all PASSED** (details in `Step_Outputs/Phase4/_env_check.md`):
1. `$MPOXHIV_ENVS/phase2/bin/python3 -c "import Bio,numpy,requests"` → Biopython 1.88, numpy 1.26.4
2. `$BIGMHC_PYTHON external_tools/bigmhc/src/predict.py --help` → full usage text
3. `$BLASTP_BINARY -version` → 2.16.0+ (matches manuscript — **not** a discrepancy, per user correction 2026-09-20)
4. `human_swissprot_db/human_swissprot.p*` → all 8 BLAST DB files present

**Rule-1 baseline:** the user supplied a durable baseline (not my first,
session-scratch attempt, which wouldn't have survived) at
`Step_Outputs/Phase4/_rule1_baseline_20260920_1850.tsv` (path/mtime/size,
1810 files) and `_rule1_baseline_sha256.txt` (content hashes). **Use these
for the final Section G check — do not regenerate them.**

---

## Step 1 — Phase4A epitope dossier (DONE)

**Scripts:** `Phase 4/_common/phase4_common.py`, `Phase 4/STEP A/Phase4A_epitopeDossier.py`
**Output:** `Step_Outputs/Phase4/StepA/Phase4A_EpitopeDossier_20260920_1855.csv` (31 rows + header)

**What it does:** no network calls. Joins the Phase 1G final construct
(`Vax_Final_4fce119e`) against four Phase 1 sources, by peptide sequence
(never row number):
- Phase1G `Phase1G_FinalConstruct_2026-08-31_1908.csv` → provenance, class, construct position, full sequence
- Phase1Db `Phase1Db_Elite_20260831_0501.csv` → Percentile_Rank, BepiPred, BigMHC-IM, SEMA overlap
- Phase1Dc `Phase1Dc_Min_50pct_2026-08-31_0501.csv` → conservancy
- Phase1F `Phase1F_Elite_Vaccine_Candidates_20260831_1742.csv` → Overall_Coverage_Pct, Coverage_Status
- Phase1Ec `Phase1Ec_Filtered_2026-08-31_1733.csv` → self-homology status

**All assertions passed:** 31 total; 10 MHC-I / 11 MHC-II / 10 B-cell; 14
HIV-derived / 17 Mpox-derived; 7 distinct source antigens (HIV_gp120,
HIV_gp41, HIV_p17, HIV_p24, Mpox_A35R, Mpox_B5R, Mpox_L1R); zero join
failures. Phase1F's and Phase1Ec's self-homology columns agree for every
epitope (bonus cross-check, not required by the brief).

**Key fact discovered:** the 31 construct epitopes are a **verified subset**
of the 105-candidate Phase 1F comparison pool (zero missing when checked).
This matters for Step 2 onward — score the 105-pool once, the 31 come along
for free, no duplicate IEDB calls needed.

**Corrections logged** (destined for `MANUSCRIPT_CORRECTIONS_PHASE4.md` at Step 7):
1. "194 Non-Redundant Sequences" (Section IV) is a mislabel — 194 is the
   count of source *proteins* from Phase 1A, not epitopes. Real units: 31
   construct epitopes / 105 comparison pool.
2. `Epitope_Provenance` labels the Mpox EEV-glycoprotein target `Mpox_B5R`;
   manuscript/correct ortholog identity is **B6R** (OPG190). MPXV's own gene
   literally named "B5R" is a different protein (OPG189, ankyrin repeat) —
   see Phase 1A's `ANTIGEN_REFERENCES`. Both names kept on every dossier row.

---

## Step 2 — Phase4B binding affinity (DONE)

**Script:** `Phase 4/STEP B/Phase4B_bindingAffinity.py`
**Output:** `Step_Outputs/Phase4/StepB/Phase4B_BindingAffinity_20260920_1954.csv`
(105 rows — the full comparison pool; 31 are construct members, flagged via
`In_Construct`) — this is the FINAL, corrected file. Earlier timestamped
CSVs from the same step (1927, 1930, 1944) are superseded; 1944 in
particular contains the IC50 parsing bug described below and should not be
used.
**Run history:** first pass finished 2026-09-20 19:27 (~430 IEDB calls,
background, ~18 min cold). A design gap (see "IC50 coverage gap" below) was
patched at 19:30 (1.9s, fully cached). A second, incorrect fix at 19:44
introduced a parsing bug that was caught and reverted at 19:53 (see "Second
finding" below) — that is the current, correct state. Every fix after the
initial cold run made zero new network calls.

**Scope:** scores the 105-candidate pool (not the 31 separately — see Step 1
finding above) against IEDB `mhci/` and `mhcii/`, methods `netmhcpan_el` and
`consensus`, using Phase 1F's full 107-allele `ALLELE_FREQ` table
(disagreement #4 — imported directly, not retyped).

**Results (revised 2026-09-20 19:44 after user review — see corrections below):**
- **Do NOT report "37/37 MHC-I pass, 29/29 MHC-II pass" as a result.**
  User correction: with a 107-allele panel, "best rank ≤ threshold across
  ANY allele" is near-vacuous — almost every peptide binds *something*.
  The informative metric is **breadth** (`N_Binding_Alleles_Primary`):
  MHC-I median 7 alleles bound (range 1–40, n=37), MHC-II median 6 (range
  1–29, n=29). **8 epitopes bind ≤2 alleles in their class-specific panels** (e.g. `LGLNKIVRM`,
  `NGLYYQGSCY`, `SSTTQYDHKESCNGL` each bind only 1/107) — these contribute
  almost nothing to population coverage despite technically "passing" the
  gate. Step 5 (Phase4E) will quantify this numerically. Pass/fail counts
  are now printed under a "context only, NOT a headline result" banner.
- Regression check vs Phase1Db, **corrected accounting**: **2 outliers among
  the 21 in-construct MHC epitopes** (`HNVWATHACVPTDPN` Δ=+5.2,
  `YKRWIILGLNKIVRM` Δ=−1.3, both MHC-II), **3 across the full 29-candidate
  MHC-II pool** (same two plus `GKLDAWEKIRLRPGG` Δ=−3.8, pool-only). Ruling:
  **benign, not an environment fault** — Phase1Db's `recommended` resolves
  to EL for MHC-I (full pool: mean Δ −0.058, **zero** |Δ|>1, n=37) but to
  BA/consensus for MHC-II (full pool: mean Δ −0.025, median 0 — no
  systematic bias — but range −3.8…+5.2, 3 of 29 with |Δ|>1). Same number,
  different scale, not a bug. **Zero gate-status flips** confirmed under
  Phase 1's own thresholds across the entire pool. **PASSED — no
  escalation.** Full reasoning recorded in `StepB/METHODOLOGY_NOTE.md`.
  **Standing instruction for Step 7:** never table Phase I's MHC-II rank
  (BA/consensus) beside Phase IV's MHC-II rank (EL) without stating they're
  different scales.
- B-cell supplementary MHC-II scan: all 39/39 pool B-cell epitopes (10 of
  them construct members) scored cleanly at length=16; got a usable
  supplementary IC50 for most (blanks eliminated — see IC50 fix below).

**Finding — IC50 coverage gap, found and fixed post-run (not a STOP
trigger, logged for transparency):** IEDB's `consensus` method supports a
much narrower allele set than `netmhcpan_el`/`netmhciipan_el` — measured
directly from cached responses: only ~31/74 (42%) MHC-I alleles and ~18/33
(55%) MHC-II alleles return a valid `consensus` result; the rest come back
`Invalid allele name X found` from IEDB itself (confirmed not a formatting
bug — the identical allele string succeeds under `netmhcpan_el` in the same
run). Consequence: the original design (IC50 read only at the primary
gate's best `netmhcpan_el`/`netmhciipan_el` allele) left
`IC50_At_Primary_Allele_nM` blank for 17/37 MHC-I and 18/29 MHC-II epitopes,
even though `consensus` *did* return a usable IC50 at some other allele.
**Fix:** added a second column pair, `IC50_At_Best_Consensus_Allele_nM` /
`IC50_Band_Best_Consensus` (IC50 at whichever allele `consensus` itself
ranked best, independent of the primary gate's allele). After the fix, IC50
coverage went from 20/37→37/37 (MHC-I) and 11/29→26/29 (MHC-II). The
original primary-gate-linked column is kept alongside it — never removed,
just supplemented. **The primary gate itself was never affected** — it's
always `netmhcpan_el`/`netmhciipan_el` at the full 107-allele panel,
regardless of `consensus` coverage.

**Two user clarifications applied (2026-09-20), both about the 10 B-cell
epitopes:**
1. B-cell epitopes have no percentile rank by design (not HLA-restricted).
   `Passes_Primary_Gate` = `"N/A (not HLA-restricted, by design)"` — never
   `UNRESOLVED`, never a fail. Matches Phase 1F's own `NOT_APPLICABLE`
   coverage status for these epitopes.
2. B-cell epitopes are 16-mers, long enough for NetMHCIIpan to score
   directly (confirmed live before coding). Captured in a **separate**
   supplementary column pair — `Bcell_MHCII_Rank_Supplementary` /
   `_Allele` / `Bcell_MHCII_IC50_Supplementary_nM` — informational only
   (T-help signal), never feeds `Binding_Alleles_Primary` or population
   coverage (Phase4E).

**Design decisions made while building this step (not pre-decided in the
brief — logged here for transparency, none are STOP-worthy):**

- **Gating thresholds:** kept Phase 1's own percentile-rank cutoffs
  (≤1.0 MHC-I, ≤10.0 MHC-II — see `Phase1Db_filtration.py`), per
  disagreement #1 ("select by percentile rank", already decided).
- **IC50 source:** `netmhcpan_el`/`netmhciipan_el` return NO ic50 column at
  all (elution-score methods only — verified live). IC50 always comes from
  the paired `consensus` call: `ann_ic50` for MHC-I, `nn_align_ic50` for
  MHC-II, read at whichever allele gave the *primary-gate* best rank. Bands:
  Strong ≤500nM / Weak 500–5000nM / Non-binder >5000nM, reported as EXTRA
  columns only (per Section E Step 2 spec) — never gates.
- **DQ/DP allele-name bridging:** `ALLELE_FREQ`'s DQB1/DPB1 entries are
  single-chain; IEDB's `mhcii/` needs paired alpha/beta notation.
  Reused+inverted `Phase1Db_filtration.py`'s existing
  `MHCII_ALLELE_TO_SINGLE_CHAIN` (covers 12 of 15 non-DRB1 alleles), added 5
  more pairings for the alleles outside Phase1Db's smaller panel (DPA1*01:03
  default for 3 DPB1 alleles, matching Phase1Db's own majority convention;
  documented DQA1*01:02/DQA1*01:03 haplotype linkage for 2 DQB1 alleles).
  **All 5 verified live against IEDB (200 OK) before use.**
- **Transient 403s:** observed IEDB returning HTTP 403 on rapid-fire
  requests that succeeded immediately on retry with a few seconds' delay —
  confirmed this is throttling, not a rejected allele/method (same request,
  same allele, retried → 200). `iedb_post_cached()` in `phase4_common.py`
  retries up to 3x with backoff (3s, 6s) before treating a call as failed.
- **Local 5-allele profile** (A*24:02, B*15:02, B*40:01, DRB1*15:01,
  DRB1*12:02): all 5 are members of the 107-allele panel already being
  scanned, so this "secondary view" required zero extra API calls — just a
  filtered read of the same results.
- **Regression check design:** Phase1Db's own allele panel (9 MHC-I / 21
  MHC-II) overlaps almost completely with `ALLELE_FREQ` (verified: 8/9
  MHC-I, ~30/33 MHC-II alleles in common) — so Phase4B's best rank across
  the larger 107-allele panel should be ≤ or close to Phase1Db's stored
  rank. Checked on all 21 non-B-cell construct epitopes (≥20 required by
  Section G). An epitope is flagged an "outlier" only if Phase4B's rank is
  *meaningfully worse* despite the broader panel (`ours > theirs + max(2, theirs)`).
  If >20% of checked epitopes are outliers → STOP trigger 4 (environment
  fault, not a scientific finding).

**Unmappable/failed alleles:** the 5 extra DQ/DP pairings added for
`netmhcpan_el`/`netmhciipan_el` all worked (pre-verified live before coding).
No `ALLELE_FREQ` allele was unmappable for the PRIMARY gate tools. A couple
of transient HTTP 500/403s occurred mid-run (`HLA-B*27:04`, `HLA-C*07:01`,
`HLA-B*58:01`, a few DP/DRB1 alleles under `consensus`) — all succeeded on
automatic retry (`iedb_post_cached`'s 3x backoff). No allele was lost to
transient failure; the `consensus`-coverage gap above is a real IEDB
allele-support limit, not a transient failure.

**Second finding, INITIALLY MISDIAGNOSED then corrected (2026-09-20
19:44→19:53) — do not repeat the 19:44 version of this fix.** User flagged 3
pool-only MHC-II candidates (`KNKRKRVIGLCIRIS`, `NKRKRVIGLCIRISM`,
`KRKRVIGLCIRISMV`) with a populated `Consensus_Best_Rank` but blank
`IC50_At_Best_Consensus_Allele_nM`. Root cause: **every MHC-II `consensus`
response row carries 4 more tab-separated fields than its own header names**
(verified: 1224/1224 cached rows) — an unlabeled trailing block, mostly `-`,
that holds the only non-dash values for a few rare alleles (e.g.
`HLA-DRB1*13:01`) where the three NAMED methods are all `-`.

**19:44 fix (WRONG, reverted at 19:53):** read that trailing block's 2nd
field as a 4th `nn_align_ic50` fallback and reported it in nM. Looked
plausible (its rank sub-value matched the row's own `percentile_rank`
exactly) but was never checked against the actual value distribution.

**User halted Step 3 and caught it** with four proofs, independently
re-verified against the raw cache before accepting:
1. Column range is **-12.5 to +5.2** (n=612) — IC50 cannot be negative.
2. Pearson r = **-0.644** against `percentile_rank` (n=612) — a genuine
   IC50 must correlate *positively* with rank; this runs backwards.
3. Strong-rank rows (0.62–1.2) sit at **+3.9 to +5.2**; weak-rank rows
   (90–98) sit at **-1.3 to -3.85** — a "higher-is-better" score, not a
   concentration that would run the other way.
4. The genuinely-labeled `nn_align_ic50`, read at its own correct header
   position in the same cache, spans **3.3–36040.7 nM, median 1938.8 nM** —
   normal IC50 magnitudes, nothing like the bad column.

**Consequence while the bug was live:** the Strong<=500nM band applied to a
[-12.5, +5.2] column marked essentially every row `Strong`, including rows
whose real (named) methods were all non-binders — `IC50_Band_Best_Consensus`
was corrupted for every row that hit the fallback, not just the 3 that
happened to surface it.

**19:53 fix (correct, current state):** the trailing-column fallback is
removed entirely — unlabeled columns are never read as typed nM values
anywhere in this step. IC50 now comes only from the two genuinely
nM-labeled columns per class (`ann_ic50`→`smm_ic50` for MHC-I,
`nn_align_ic50`→`smm_align_ic50` for MHC-II, via a `genuine_ic50()`
helper); where both are `-`, the value is honestly `"N/A"`. The 3 flagged
peptides now correctly show `IC50_At_Best_Consensus_Allele_nM = "N/A"` (no
IC50-producing method ran for `HLA-DRB1*13:01` on those peptides) — not a
fabricated 4.7nM. The unlabeled block itself is NOT exposed as its own
column (its identity is unconfirmed; not worth the risk for 3/105 rows).

**Corrected, honest IC50 coverage (`IC50_At_Best_Consensus_Allele_nM`):**
MHC-I 37/37 (100%), MHC-II 26/29 (~90%, 3 honestly N/A) — lower than the
buggy run's numbers and the correct ones. All 4 proofs and the full
before/after are recorded in `StepB/METHODOLOGY_NOTE.md`'s "CAUGHT-AND-
CORRECTED PARSING ERROR" section.

**Lesson for future steps:** never read a column by position alone without
checking it against the header's own column count, and never trust a typed
numeric value (nM, percent, etc.) without sanity-checking its range and
correlation direction against a known-good reference column in the same
response.

**Rule 1 re-verified after both the 19:44 and 19:53 fixes:** diffed against
the same durable baseline — zero mtime/file changes under Phase 1/Phase 2
both times. All fixes required zero new IEDB calls (each re-run <0.35s,
fully served from cache).

**Reminder for next session:** `HALT-FOR-OPUS.md` does not exist — no
escalation triggers fired in this step. Phase4B's `Binding_Alleles_Primary`
column (netmhcpan_el/netmhciipan_el gate) is what Step 5 (Phase4E, coverage)
should consume — not `Consensus_Best_Rank` or either IC50 column, which are
reporting-only. **Never present "N pass the gate" as a finding** — use
`N_Binding_Alleles_Primary` breadth instead. **Never table Phase I vs
Phase IV MHC-II ranks together without noting the BA/consensus-vs-EL scale
difference** (Step 7 standing instruction, see above).

---

## Step 3 — Phase4C TCR recognition (DONE)

**Script:** `Phase 4/STEP C/Phase4C_tcrRecognition.py`
**Outputs:**
- `Step_Outputs/Phase4/StepC/Phase4C_TCRRecognition_MHCI_20260921_1704.csv` (10 rows — the MHC-I construct epitopes)
- `Step_Outputs/Phase4/StepC/Phase4C_CrossReactivity_AllConstruct_20260921_1704.csv` (31 rows — all construct epitopes)

**Scope decision:** the TCR-recognition ladder (BigMHC-IM, PRIME/MixMHCpred,
Calis) is **MHC-I only (10 epitopes)** — all three tools are MHC-I specific
by design, and Section F's own wording ("TCR-facing residues = positions
4-6 of each MHC-I epitope") confirms this scope. Cross-reactivity
(self-homology) is broader — applied to all 31 construct epitopes, since
self-tolerance is a safety question independent of TCR-recognition tooling.

**Ladder result — all 3 rungs used, nothing dropped:**
1. **BigMHC-IM** (primary): 10/10 epitopes scored (70 peptide×allele pairs,
   using each epitope's own Phase4B `Binding_Alleles_Primary` set).
2. **PRIME 2.0 + MixMHCpred 3.0** (optional): installed fresh this run from
   GitHub (GfellerLab/PRIME, GfellerLab/MixMHCpred) — both free, no licence
   wall. 10/10 epitopes scored. See install friction below.
3. **Calis et al. (2013)**: reimplemented from the paper's own public
   training data (not from memory). All 10 epitopes scored.

**Direction/sign checks (mandatory before trusting any score) — n and p
stated for every correlation (corrected 2026-09-21 after review):**

| Tool | Convention | Expected vs MHC-I rank | r | p | n |
|---|---|---|---|---|---|
| BigMHC-IM | 0-1, higher=more immunogenic | negative | -0.298 | 0.402 | 10 |
| PRIME %Rank | lower=better | positive | **+0.916** | 0.0002 | 10 |
| Calis (EXPLORATORY) | higher=more immunogenic | no prior | +0.097 | 0.790 | 10 |

- At n=10 only |r| > ~0.63 is distinguishable from zero. **PRIME's +0.916 is
  meaningful. BigMHC-IM's -0.298 has the right sign but is NOT statistically
  distinguishable from zero** — it passes the sign check, it corroborates
  nothing.
- **Earlier "BigMHC-IM vs Phase1Db r=+0.460 (same-tool check)" was WRONG and
  is withdrawn.** Phase 1 stored the max over ITS OWN allele set; this
  run's best is over the 107-allele-panel set — for 9/10 epitopes a
  different set, so the two columns are different quantities. Replaced by a
  true same-allele check: re-scored Phase 1's exact (allele, peptide) pairs
  → **r=+1.0000, n=10, max |diff| = 3.97e-07**. BigMHC-IM reproduces
  Phase 1 to numerical precision; the "best" columns differ only because the
  107-allele panel found better alleles (correct). (Reviewer independently
  confirmed the one shared-allele case, QTSVFSATVY/HLA-B*15:02, is identical
  to 4dp.)
- **Calis is exploratory only** (self-AUC 0.62, r=+0.097). Near-zero
  discrimination vs binding rank is arguably correct — the model is designed
  to be independent of binding affinity — but that also means it is NOT
  independent corroboration of BigMHC/PRIME and must never be presented as
  such. Column renamed `Calis_Score_Exploratory`.
- **BigMHC-IM and PRIME pick their "best" allele independently** — same
  allele for only 3/10 epitopes (`BigMHC_PRIME_Best_Allele_Match` column).
  **Never compare their "best" columns as if allele-matched; say so in the
  Step 7 report.**

**STANDING FLAG for Steps 5 and 7 — `NKRKRVIGL` (Mpox A35R):** worst in the
construct on ALL FOUR independent signals at once (`Worst_In_Construct_On`
column, computed not hand-picked): PRIME %Rank 4.386 (next worst 0.47),
BigMHC-IM 0.0344, binding breadth 2/74 class I alleles, EL rank 0.85. Phase 1
scored it at `HLA-B*08:01`, which is not in the 107-allele panel at all
(0.163 there). Ties to A35R's lowest-in-table BigMHC score in the paper and
to the `A35R_Coverage_Note` substitution in Phase1G. **Expect it to add
almost nothing to Step 5's coverage union. Step 7 must report this as a
consistent cross-phase signal, NOT a new anomaly.**

**PRIME/MixMHCpred install — worked, but hit two environment problems along
the way (both fixed, neither required editing the vendored tools):**
1. Both tools' own launcher scripts refuse a space anywhere in their
   install/working path ("Spaces in path... are not supported") — this
   project's root sits under the space-containing "Extended SSD" volume
   name. The project's existing `~/mpoxhiv_ssd` space-free symlink fixes
   this for the bash launchers and binary paths themselves. **But PRIME's
   *compiled* binary additionally calls `realpath()` on its `-i`/`-o` FILE
   arguments internally, which dereferences that same symlink straight back
   to the real, space-containing path** — so the symlink trick that worked
   for Step 0's conda envs does NOT fully solve it here. Fixed by staging
   PRIME's own input/output files in a small scratch dir under `$HOME`
   (`~/.mpoxhiv_phase4c_prime_io/`, no symlink anywhere in that path at
   all) instead. This is transient plumbing only — every actual result is
   still parsed out and written to `Step_Outputs/Phase4/StepC/` on the SSD
   as normal; the raw PRIME response is also copied into the SSD's
   `_tool_runs/` cache for the audit trail.
2. Both tools invoke a bare `python3` internally (needs pandas) which
   resolves to system Python, not the project's `phase2` conda env. Fixed
   by prepending `$MPOXHIV_ENVS/phase2/bin` to `PATH` for the subprocess
   call only.

**Calis et al. (2013) reimplementation — full provenance, not fabricated
from memory.** Initial web searches for the published coefficient table
came up empty (IEDB's help pages don't document it; the PLOS article's
supplementary table S1 turned out, on inspection, to be the paper's full
2,508-peptide labelled TRAINING DATASET, not a coefficient table). Given
that, the amino-acid log-enrichment table was **derived directly from that
public dataset** using the paper's own documented method (quoted from PMC
full text): enrichment = ratio of an amino acid's frequency in immunogenic
vs non-immunogenic peptides, pooled over non-anchor positions (anchors =
positions 1, 2, and the C-terminus, generic/allele-agnostic). Per-position
importance weights are the paper's own quoted Table 2 values (position
3=0.10, 4=0.31, 5=0.30, 6=0.29, 7=0.26, 8=0.18; positions 1/2/9 are
anchors). **Validated before use:** self-scored AUC 0.62 on the 2,508
training peptides (not held-out — a sanity check only), and the derived
table's own qualitative pattern independently matches the paper's reported
findings (tryptophan/aromatic residues most enriched, serine most
depleted). Full derivation, quotes, and validation numbers are in
`StepC/METHODOLOGY_NOTE.md`.

**Cross-reactivity (disagreement #7, both layers, all 31 construct
epitopes):** BLASTP layer reused directly from Phase1Ec (frozen, Rule 1 —
not recomputed). Exact 8-mer layer computed fresh via
`phase1_common.human_self_homology()` (imported per Section F's explicit
instruction, not reimplemented) and cross-checked against Phase1Ec's own
stored exact-match columns: **31/31 agree.**

**Rule 1 re-verified after this step:** zero changes under Phase 1/Phase 2.
A second run of the script completed in ~1s using cached BigMHC-IM and
PRIME responses (no re-invocation of either tool).

---

## Step 4 — Phase4D proteasomal processing (DONE; trigger 3 RESOLVED by Opus ruling — see Step 5)

**Script:** `Phase 4/STEP D/Phase4D_processing.py`
**Outputs (`Step_Outputs/Phase4/StepD/`, timestamp 20260921_1844):** `Phase4D_EpitopeProcessing` (31 rows), `Phase4D_CleavageProfile` (563 end-positions),
`Phase4D_JunctionNeoepitopes_MHCI` (126 rows), `Phase4D_JunctionNeoepitopes_MHCII` (4022 rows), `METHODOLOGY_NOTE.md`. Escalation: `Step_Outputs/Phase4/HALT-FOR-OPUS.md`.

**Endpoint semantics (verified live before use):** `total = proteasome + tap + mhc`; `mhc_score = -log10(IC50 nM)`; higher = better for every column.
`proteasome_score` is a pure C-terminal cleavage score, one value per end position (0 deviations across 165,641 comparisons over all alleles/lengths) — so the
whole construct has one cleavage profile; tap/processing are allele-independent. Direction checks: residues by cleavage score W>L>F>Y>…>P>G (matches known
proteasome specificity); endpoint IC50 vs Step 2 per-allele EL rank, Spearman rho=+0.468, p=4.5e-5, n=70 pairs from 10 epitopes (NOT independent) — correct sign.
(My earlier "IC50<=500nM" check was an unvalidated prior and was removed.) No `tap=` parameter sent.

**Scope decision:** processing endpoint is MHC-I only -> liberation verdicts for the 10 MHC-I epitopes; MHC-II/B-cell get boundary cleavage scores as supplementary,
verdict N/A. MHC-II junctions (12-20 aa) scored with NetMHCIIpan-EL; the 11 real MHC-II epitopes scored in the same calls. Tox/allergen half of IV.B.2 = Phase 2A (referenced, not re-run).

**Pre-registered definitions (fixed before results):** NOT liberated = C-terminal cleavage percentile <50. Junction window = spans >=1 segment boundary; "deep" = >=3
residues each side of every boundary (fewer = shifted copy of a real epitope, reported but excluded from the trigger). Neoepitope flag = deep AND strong binder
(IC50<=500nM / EL rank<=10) AND out-scores >=1 real epitope binding the SAME allele.

**Results (complete run: 296/296 MHC-I tables, 297/297 MHC-II calls, 0 failures):**
- Liberation: **10/10 MHC-I epitopes liberated**; that half of trigger 3 did NOT fire. Secondary observation (not part of the trigger): **5 epitopes
  (GPKEPFRDY, HHFNCRGEF, HWTTYMDTF, RFALNPGLL, VYSTCTVPTM) have an internal cleavage site scoring higher than their own C-terminal one** — possible destruction risk.
- MHC-I junctions: 1475 windows cross a boundary (526 deep). **126 (window,allele) flags on 32 deep windows; 103 beat ALL real epitopes at that allele.**
  26/32 windows involve the AAY linker (GPGPG 4, EAAAK 1, KK 1). Mechanistic observation only: AAY ends in Tyr, a classic MHC-I C-terminal anchor and favoured
  cleavage residue, so epitope-end+AAY windows (e.g. TYMDTFAAY @ C*14:02, IC50 4.7 nM; CTVPTMAAY @ A*26:01, 6.2 nM; EPFRDYAAY / FSATVYAAY @ B*35:xx) are strong and well-cleaved.
- MHC-II junctions: 4465 windows (12-20 aa) x 33 alleles -> **4022 flags on 1117 windows; 3033 have the binding core spanning the junction; 1757 beat ALL real
  epitopes at that allele.** (e.g. PQDLAAYNKRKR @ DRB1*13:01, rank 0.01, core LAAYNKRKR.)
- **Trigger 3 fired** (junction neoepitopes out-scoring real epitopes). Per escalation rules only the *interpretation* is paused; all tables are final.
  Options recorded for Opus: (1) accept as inherent to any concatenation of binders, report severity tiers; (2) report as a design finding/limitation; (3) redesign
  linkers (changes the frozen construct — Opus/PI decision).
- Caveats: "out-scores a real epitope" is a weak bar where an allele has few/weak real comparators (see `Real_Epitopes_At_Allele`, `Weakest_Real_Total`, `Best_Real_Rank`
  columns); counts are per (window, allele), not frequency-weighted — Step 5 will show how much population weight they carry. EL rank (MHC-II) and processing
  IC50/total (MHC-I) are different scales; Phase I MHC-II ranks (BA/consensus) are not comparable to EL ranks (Step 7 standing instruction).

**Run-integrity incident (caught and corrected; do not repeat):** the first 4-worker run was throttled by IEDB (HTTP 403): 64/296 MHC-I tables (all HLA-C alleles +
a few others) and ~220 MHC-II calls failed silently while the script analysed the partial data (it reported "60/74 alleles usable", 111 flags, and a spurious
"15/61" consistency figure). All of it was discarded; direct tests showed the "missing" alleles ARE supported. Fix: 2 workers, 6 retries, and the script now
**refuses to analyse any incomplete scan** (exit code 2). Complete-data numbers differ (126 vs 111 MHC-I flags). Lesson: parallelism against IEDB needs a completeness
guard, and "unusable allele" must be verified, not inferred from a failed call. Speed check: server-bound (~19 s per 500-peptide MHC-II call); 2 workers is the
safe ceiling; skipping shallow windows would save ~27% but drop informational counts, so it was not done.
**Also:** `halt_for_opus()` is now idempotent per branch; `PHASE4_PROGRESS_LOG`/memory files updated; Rule 1 re-verified clean (0 changes).

---

## Step 5 — Phase4E population coverage + frequency-weighted junction analysis (DONE; RE-ISSUED; no open escalation)

**Script:** `Phase 4/STEP E/Phase4E_populationCoverage.py` (pure local; runs with the network disabled).
**Outputs (`StepE/`, 20260921_1934):** `Phase4E_EpitopeCoverage` (31), `Phase4E_ConstructCoverage` (margin curve + levels), `Phase4E_ContestedAlleles`, `Phase4E_CoverageGaps` (107 alleles), `Phase4E_MethodDivergence` (6 loci), `METHODOLOGY_NOTE.md`.

**History (do not repeat):** the first version used the original "contested" definition (any seam window beats any real epitope by any amount). I ran an unrequested null calibration that showed it is an
extreme-value comparison (2,491 class II seam windows per allele vs 11 reference epitopes) and reported both. **Opus ruled it a specification flaw**, discarded its figures (they must not be reported),
and issued a corrected definition. Everything below is the corrected re-issue.

**Two allele universes, never mixed:** Phase 1F (`Binding_Alleles_Recomputed`, MHCflurry/MHCnuggets) only for the level-1 regression gate; Phase 4B (IEDB EL) for everything else. All arithmetic via Phase 1F's imported helpers;
the methodology's formula is implemented nowhere (grep clean). **Gates passed:** Level-1 **31/31** reproduce Phase 1F (B-cell `NOT_APPLICABLE`, never 0.0); loci grouped separately {A,B,C,DRB1,DQB1,DPB1} — **DQA1/DPA1 do not
exist in ALLELE_FREQ** (beta chains only); contested ≤ epitope coverage (21/21); uncontested ≤ raw (every scope/margin); zero network; Rule 1 SHA-256 1,810/1,810.

**1. Coverage gaps (Step 7 must lead with this).** `A*11:01` (f 0.177) and `A*33:03` (0.109): no construct MHC-I epitope binds them under EITHER predictor — unambiguous, no null needed; explains locus A 60.3% (IEDB-EL) / 76.0% (MHCflurry).
`C*08:01` (0.187), `B*18:01`, `C*04:03`, `C*07:04` are gaps under IEDB-EL only (predictor-dependent). **Correction to the ruling:** it called C*08:01 a gap without qualification; MHCflurry has TMGAASITL binding it.
Frequency mass with no binding epitope (IEDB-EL / MHCflurry): A .446/.306, B .224/.000, C .491/.080.
**2. Method divergence (report the range, never the favourable end):** per-locus IEDB-EL vs MHCflurry/nuggets: A 60.3/76.0, B 89.4/99.0, **C 70.2/98.2 (28 pp)**, DRB1 99.5/99.5, DQB1 96.0/95.3, DPB1 87.0/84.1.
**3. Levels 2–3:** combined saturates >99.99% (per-locus/per-class is the informative layer). Mpox MHC-I 83.9% vs HIV MHC-I 97.5% (IEDB-EL) — the asymmetry, consistent with A35R weakest throughout; no trigger 7 (both >99.98% combined).
`NKRKRVIGL`: own coverage 4.72%, marginal drop 0.0000 pp (construct) — standing flag, consistent not new.

**4. Corrected junction analysis.** Contested at allele a iff: deep + liberated + beats the BEST real epitope at a by ≥m-fold (MHC-I IC50; MHC-II EL-rank ratio, EL has no IC50) + (class II) core spans the seam. Curve, n stated
(pairs at risk: MHC-I 70, MHC-II 81; alleles with a real comparator: I 43, II 33):

| m | contested pairs I / II | contested alleles I / II | uncontested: construct / MHC-I / MHC-II |
|---|---|---|---|
| 1× | 50/70, 66/81 | 32 / 28 | 69.4 / 59.9 / 23.6 (extreme-value regime — curve only, not a finding) |
| 2× | 40/70, 60/81 | 26 / 25 | 83.9 / 68.3 / 49.3 |
| **5× (primary)** | 28/70, 43/81 | 16 / 18 | **99.88 / 88.7 / 98.9** |
| 10× | 24/70, 24/81 | 14 / 11 | 99.98 / 89.1 / 99.8 |

Sensitivity (class I also requires IC50 ≤500 nM): uncontested MHC-I 76.5 / 81.1 / 90.0 / 90.2%. **Escalation: none** — whole-construct uncontested 99.88% (gap 0.12 pp) at the primary margin; trigger 7/4 not fired.
**Flagged for Opus:** the whole-construct **1× point is 69.36%, just under the ~70% floor**; not treated as an escalation because the ruling defines contest at ≥5× (1× is the extreme-value regime) — please confirm.
**Per locus at 5× (the honest layer; combined hides it because loci multiply):** A 60.3→57.8, **B 89.4→32.6, DPB1 87.0→0.0** (every DPB1 allele a real epitope binds is contested via KK/AAY seams), C 70.2→60.3, DRB1 99.5→88.2, DQB1 96.0→91.0.
HIV/Mpox: pathogen unions 99.0 / 99.4% uncontested; MHC-I subsets HIV 97.5→80.9, Mpox 83.9→62.0 (the weak arm gets weaker).

**Verified the ruling's claims against the data:** ✓ all 8 strongest class I windows end in AAY (IC50 4.7–13.3 nM; B*35:01/05/17, A*26:01, C*14:02); ✓ no top-3-frequency allele contested (most frequent contested: B*15:02, f .084, rank 9/74);
✓ real epitopes at those alleles end in F/Y (aromatic-PΩ reading consistent, not proven by this data). **Corrections:** overall 21/31 (68%) qualifying windows end in AAY (29/31 involve AAY, 8 contain it without ending in it, 2 involve only EAAAK/GPGPG);
watchlist∩AAY-seam epitopes = **four** (GPKEPFRDY, HWTTYMDTF, RFALNPGLL, VYSTCTVPTM), not three.
**GPGPG class II seams have no 5× contested rows; residual KK/AAY calls remain:** EL-rank calibration n = 82,203 → rank ≤10 7.17% (null 10%, ×0.72), rank ≤1 0.76% (null 1%). By linker (×null at rank ≤10): GPGPG 0.59, KK 0.65, EAAAK+adjuvant 0.49, AAY 1.14 (only mild enrichment). **0 of 193** class II contested rows (m=5) involve GPGPG
(108 KK, 83 AAY, 2 EAAAK/adjuvant) — the residual sits at KK seams (flanking B-cell epitopes) and AAY seams, not around class II epitopes. Caveat: no random-window null was run at the 5× criterion itself, so no claim on whether that residual exceeds chance (cost 1.1 pp).
**Class I:** IC50 ≤500 nM in 0.75% of seam pairs vs 0.42% of unselected natural adjuvant windows (×1.81); windows ending in AAY 6.12% (×14.7); GPGPG 0.06%, KK 0.09% (depleted). Mechanism: AAY's terminal Tyr in the MHC-I C-terminal pocket. **Discussion/limitations finding; no redesign** (Opus).

---

## Step 6 — Phase4F B-cell + conservancy (DONE; all regression gates passed)

**Script:** `Phase 4/STEP F/Phase4F_bcellAndConservancy.py` (no network, no new predictions). **Outputs (`StepF/`, 20260921_1933):** `Phase4F_BcellDossier` (10), `Phase4F_BcellPool` (39), `Phase4F_Conservancy` (31), `Phase4F_ConservancyGate` (7), `METHODOLOGY_NOTE.md`.

**Conservancy gate — reproduced EXACTLY:** HIV gp120 4.54, gp41 5.69, p17 4.67, p24 9.97; Mpox A35R 37.80, B5R 45.02, L1R 60.61 (means over all 30,268 candidates). **Source note:** the brief says `Min_50pct`, but that file holds only the ≥50% survivors
(1,238 rows) and cannot reproduce a full-pool mean — the gate runs on `Raw_Conservancy/Phase1Dc_Raw_Full`. Independent recompute for the 31 construct epitopes from the Phase 1C variant pools: **31/31 reproduce stored value AND hit ratio**.
Construct conservancy: HIV n=14 mean 75.5% (53.3–96.7%); Mpox n=17 mean 99.7% (95.7–100%) — the conserved tail of a diverse pool (candidate means 5.5% / 46.6%). Variant pools small (n = 21–30) → coarse, within-clade only.

**B-cell (10 construct; 39 pool comparison; HIV n=4 / Mpox n=6):** Bcell_Tier recomputed with the imported `classify_bcell_tier`: 39/39 match (Phase I tiers kept, disagreement #3). Construct: 1 High, 6 Medium, 3 Deprioritized (retained — tiers deprioritise, never exclude);
methodology-High (secondary stratum) 4/10 (HIV 1/4, Mpox 3/6). **SEMA-3D (not re-run):** all 10 `NO` (0–12.5% of 16 residues); overlap recomputed from Phase 1De per-residue scores for every locatable epitope (49/49 match). `NO` = screened, <50% in a patch — lower-priority
linear-only, **not** a refutation (SEMA-3D TPR 0.66 at FPR 0.07); `UNSCREENED` = no assessable fold, never a negative (none in the construct). **Native-fold overlay:** construct model verified unusable (pTM 0.17; adjuvant mean pLDDT 87.2, 98% ≥70; the 525 non-adjuvant residues mean pLDDT 27.2,
**0.0% ≥70**) so exposure uses Phase 1De per-antigen folds. Sanity: B-factor pLDDT reproduces Phase 1De Mean_pLDDT for 7/7; RSA polar 0.47 > hydrophobic 0.30 (correct direction). 5/10 epitopes on a DETERMINED fold (gp120, p17, p24, L1R), 5/10 on UNDETERMINED (gp41 pTM .39, A35R .47, B5R .43 — exposure indicative only).
Mean RSA HIV 0.37 / Mpox 0.34; native pLDDT 79.6 / 87.3 vs construct-model pLDDT 20.8–49.8 for every epitope. BepiPred vs RSA across the 39 unique pool candidates: Spearman +0.112, p = 0.5 — right sign, not distinguishable from zero (passes the sign check, corroborates nothing).
Caveats: monomer folds; HIV Env glycans not modelled (exposure over-estimated); predictors, not antibody-binding evidence.

**Rule 1: SHA-256 1,810/1,810 byte-identical after Step 6.**

---

## Step 7 — Phase4G integrated report (DONE — **PHASE IV COMPLETE**)

**Script:** `Phase 4/STEP G/Phase4G_integratedReport.py` (~630 lines; pure synthesis, no predictions, no network). Ran clean on the 5th iteration ("run5") after fixing 4 self-caught script
bugs along the way (invalid nested f-string; a stray `%` inside an f-string later `%`-formatted, `TypeError`; a blank-cell gate false positive from checking in-memory dicts instead of the
written CSV; an `UnboundLocalError` from a variable used before its later computation).

**Outputs (`StepG/`, timestamp `20260921_2022`):**
- `Phase4G_Immunogenicity_Dossier_20260921_2022.csv` — 31 rows × 153 columns. `Primary_Gate` (Phase I's own rank thresholds) and `Secondary_Stratum` (the revised methodology's stricter
  values) as two independent verdict columns per epitope, plus `Standing_Flags` and `Contest_Reading`.
- `PHASE_IV_REPORT.md` — 8 sections ordered by evidence strength (coverage gaps → method divergence → binding/TCR/cross-reactivity → per-pathogen → proteasomal processing → junction analysis
  → B-cell/conservancy → what C-ImmSim's replacement does not give).
- `MANUSCRIPT_CORRECTIONS_PHASE4.md` — Part A: the 8 pre-decided disagreements with corrected text and citations back to the deciding step. Part B: 18 further items (B1–B18): the
  194-sequences-vs-epitopes mislabel, B5R/B6R naming, BigMHC "stability" input claim, MERCI "database" wording, the Step 2 IC50 parsing incident (with all 4 proofs), the Phase I/IV MHC-II
  scale difference, the IEDB `immunogenicity/` 403 (dead endpoint, not throttling — corrects an earlier review note, with the retest evidence), the conservancy `Min_50pct`→`Raw_Full`
  file-pointer correction, and more.
- `Phase4G_NumbersLedger_20260921_2022.csv` — every number quoted in the report recorded as `(key, value, source_file)`; 155 rows.

**Self-review before close-out:** read the full generated `PHASE_IV_REPORT.md` and `MANUSCRIPT_CORRECTIONS_PHASE4.md` end to end (not just the script's own pass/fail gates) and found/fixed
~10 wording/accuracy issues: odd non-integer breadth-median formatting (added an `fm()` helper); a coverage-gap sentence that double-counted frequency mass; per-pathogen coverage percentages
that had been truncated with `[:5]` string slicing instead of `.2f` (also corrected a resulting ">99.9%" claim to the accurate ">99.6%"); a leftover internal planning-note sentence; imprecise
class-II "sample-size asymmetry" wording tightened with an explicit caveat that no null was run at the 5× criterion itself; a grammar fix; a gloss explaining BLASTP `PARTIAL` statuses are
non-significant at this length; a full-pool gate-status-flip check (0/66 flips, added as a corroborating check); a new "Two verdict readings per epitope" subsection after discovering all 10
MHC-I epitopes fail `Secondary_Stratum` specifically because none reaches BigMHC-IM's 0.70 cut (not because of IC50). None of these changed a numeric finding, only presentation/accuracy.

**Final gates — independently re-verified in this session, not just trusted from the script's internal check:**
- Rule 1: `shasum -a 256 -c Step_Outputs/Phase4/_rule1_baseline_sha256.txt` → **1810/1810 OK, 0 changed/missing.** (An earlier ad-hoc Python re-check misread the wrong column of the wrong
  baseline file — the `.tsv` holds mtime/size, not sha256 — and falsely flagged all 1809 data rows as "changed"; the actual sha256 file, checked properly, is clean.)
- Dossier: 31 rows, 0 cells containing `UNRESOLVED`, 0 blank cells, `Primary_Gate`/`Secondary_Stratum` present on every row (re-read from the written CSV with a fresh script, not the
  in-memory dicts).
- Ledger: 155/155 rows have a non-empty `Source_File`.
- No superseded figure (22.06 %, 77.93 pp, 18.67 pp, DPB1 0.0 %, the 111-flag count) appears in either report.

**Escalation status:** Step 7's only escalation condition (trigger 4, a regression gate failing) did **not** fire. All 3 prior escalations (Step 2 misparse; Step 4 trigger 3; Step 5 second
firing) remain resolved and were not re-opened. `HALT-FOR-OPUS.md` still holds exactly those two RESOLVED branches; no third branch exists.

**Documentation gates closed this step:** `HANDOFF_SECTIONS_III_V_VI.md` §0 marked **PHASE IV COMPLETE**; §7 gained findings 11–12 (TCR recognition significance/allele-matching caveat; the
Primary/Secondary verdict split); §9 gained trap 17 (immunogenicity/ endpoint is dead) and trap 1 was corrected twice on the IC50-vs-rank correlation: an initial pass wrongly reported
+0.504/+0.876 as non-reproducing; on independent recheck (`IC50_At_Best_Consensus_Allele_nM` vs `Consensus_Best_Rank`) **both pairs are real** — Pearson +0.487/+0.752 on raw IC50 and
+0.504/+0.876 on log10(IC50), same two columns, different transform (Spearman +0.508/+0.909 either way, transform-invariant) — see corrections file B8; §12 got the final changelog entry.

**Phase IV is complete: all 7 steps done, 0 open escalations, Rule 1 clean (1,810/1,810 files byte-identical throughout).**

---

## Open questions / things to double-check later

- None outstanding. Phase IV is closed. If Section III/V/VI work surfaces a Phase IV number that doesn't reproduce, or a reason to re-open a branch, log it here before acting on it.
