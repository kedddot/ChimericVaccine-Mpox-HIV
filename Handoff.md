# Research handoff: completed Phase IV and the next evidence boundary

**Repository:** `/Volumes/Extended SSD/ChimericVaccine/Research`  
**Construct in the completed branch:** `Vax_Final_4fce119e`, 570 amino acids  
**Consolidated:** 2026-09-24  
**Owner of this document:** future researcher integrating verified Phase III evidence with the completed Phase IV analysis

## 1. Current status and which methods control

The **revised-methodology Section IV** is an IEDB-centered, epitope-level analysis. Its Phase I/II-only branch is complete through Steps A–G. The appended “Proposed Methods (Additional)” Section IV on pages 40–42 of the latest manuscript PDF describes an older daughter-construct and repeated-docking/C-ImmSim plan; those daughter-construct comparisons were **not outputs of this completed branch**. When the standalone revised methodology differs from the frozen construct or the executed Phase IV pipeline, describe the executed method and record the difference in `PAPER_UPDATES_PHASE_IV.md`.

`PHASE_IV_RESULTS_DETAILED.md` is the **single current narrative of Phase IV results**. `PAPER_UPDATES_PHASE_IV.md` is the **single Section IV manuscript revision package**. The Step G report and Step G corrections document now point to these files; they are no longer competing narratives. The original numerical evidence remains in the Step A–G CSVs, especially the final 31-row dossier and 155-entry numbers ledger under `Step_Outputs/Phase4/StepG/`.

| Scheme | Relevant meaning | What exists |
|---|---|---|
| Manuscript old Section III | Docking, MD, mRNA and C-ImmSim | Reported in the PDF; this repo visibly holds Phase 3 Step A receptor/ligand preparation. Raw later-run provenance is unresolved. |
| Manuscript old Section IV | Daughter constructs and repeated Phase II/III comparisons | Appended proposal, not the executed revised Phase IV branch. |
| Revised Section IV | IEDB-centered immunogenicity analysis | Steps A–G and final CSV outputs complete, using Phase I/II inputs and no Phase III/docking input. |
| Revised Sections III, V, VI | Delivery, pharmacokinetics, validation roadmap | Separate scopes. `Phase 5/` exists but was not audited for this handoff. |

## 2. What was verified in this handoff

- The final construct dossier has **31 rows and 153 columns**: 10 MHC-I, 11 MHC-II, 10 B-cell epitopes; 14 HIV-derived and 17 Mpox-derived. No cell is blank or `UNRESOLVED`; both primary and strict-secondary verdict fields are present. The Phase 1F comparison table has **105 candidates** (37 MHC-I, 29 MHC-II, 39 B-cell).
- The Step G numbers ledger contains **155 entries with source fields**. A source field is a pointer, not a substitute for checking the named CSV. Key claims in the consolidated results were rechecked against Step B–F tables; the remote predictors were **not** rerun.
- All **1,810** files in the saved Phase I/II SHA-256 baseline matched on 2026-09-24. The Phase IV code has no Phase III path dependency. No Phase I, II or III source/output file was edited in this consolidation.
- The existing progress log and per-step `METHODOLOGY_NOTE.md` files record execution history, including rejected intermediate calculations. They are historical records, **not** alternate current results. The final CSVs and `PHASE_IV_RESULTS_DETAILED.md` control quantitative reporting.

## 3. Minimum findings the next researcher must preserve

1. The two reference-panel HLA-A alleles `A*11:01` (frequency 0.177) and `A*33:03` (0.109) have **no binding final-construct MHC-I epitope under either predictor**. Under Phase IV IEDB-EL, none of the 37 comparison-pool MHC-I peptides binds either allele. Phase I's nine-allele discovery panel omitted them; this is a plausible contributor, not a proven counterfactual cause.
2. Per-locus projected coverage depends on the predictor, particularly class I: A **60.31/75.99%**, B **89.37/98.96%**, C **70.19/98.18%** (Phase IV IEDB-EL / Phase I MHCflurry-MHCnuggets). DRB1 is 99.54/99.54%, DQB1 96.00/95.25%, DPB1 86.97/84.08%. The combined model value saturates and is not a measured population response. Frequencies are from AFND Singapore Riau Malay (n=132), a Southeast Asian proxy, **not a Philippine HLA panel**.
3. Under IEDB-EL, the HIV and Mpox MHC-I subsets have projected coverage **97.52% and 83.87%**, respectively (n=5 selected peptides per pathogen). `NKRKRVIGL` is a weak predicted Mpox A35R peptide in the saved metrics. This does not establish immune dominance, interference, or lack of protection.
4. All **10/10** selected MHC-I peptides meet the pipeline's **predicted** C-terminal processing rule; five have a stronger predicted internal cleavage site. At the primary 5× seam margin, **28 class-I epitope–allele pairs** are contested. The seam strong-binder rate is **293/38,924 = 0.75%** versus **45/10,804 = 0.42%** in the limited adjuvant-domain window reference (descriptive ratio 1.81×). Of 31 qualifying class-I windows, 29 involve AAY. A “contested” call does **not** mean measured loss of a pathogen epitope.
5. **0/193** class-II 5× contested rows involve GPGPG; the remaining rows are 108 KK, 83 AAY and 2 adjuvant. The full 5× criterion has no matched random-window null, so do not describe all class-II seams as risk-free. All ten selected B-cell epitopes were SEMA-3D screened without corroboration at the chosen threshold; `NO` is not structural refutation. Five map to globally undetermined native folds.

For all denominator, rank-scale, row-level and limitation details, use `PHASE_IV_RESULTS_DETAILED.md`. Do not replace the project's 0.5 BigMHC-IM selection convention with the standalone methodology's 0.70 strict stratum. All ten selected MHC-I peptides fail the **strict secondary** stratum because none reaches 0.70; four also fail its IC₅₀ ≤500 nM condition. Phase I MHC-II BA/consensus and Phase IV MHC-II EL ranks are different measures.

## 4. Phase III / molecular-docking provenance gate

The latest 42-page manuscript PDF reports old-scheme Section III docking, 100-ns molecular dynamics, mRNA design and C-ImmSim results. In this repository, `Phase 3/` and `Step_Outputs/Phase3/` show receptor/ligand preparation only. An older session note states that later work had not yet been implemented at that time, and a prior volume search found no raw HDOCK/GROMACS/C-ImmSim files. This is an **unresolved provenance conflict**. It does not prove that those analyses never happened, and it does not authenticate the PDF's reported numbers.

Before adding any Phase III result to Phase IV:

1. Obtain the exact raw runs and metadata: ligand/construct sequence and version; receptor structure and chains; preparation steps; software and version; input parameters; docking poses and ranking; MD topology, trajectories, energies and convergence checks; C-ImmSim inputs, seeds, schedules and outputs where relevant.
2. Check that every reported PDF table/figure can be reproduced from those files and that its ligand matches `Vax_Final_4fce119e` or the explicitly identified 45-residue adjuvant. Do not transfer results from a superseded construct without a new bounded analysis.
3. Verify structural scope. The Phase IV Step F summary records **whole-model pTM 0.17**, **adjuvant mean pLDDT 87.2 (97.8% ≥70)** and **non-adjuvant mean pLDDT 27.2 (0% ≥70 across 525 residues)**. Whole-construct docking against one arbitrary low-confidence conformation must not be presented as robust receptor binding. The Phase 3 Step A prepared receptor/ligand files are read-only references; inspect their identities before use.
4. Link each validated Phase III observable to a claim it can support. A predicted receptor contact cannot be converted into epitope-specific T-cell activation, an antibody titre, cytokine kinetics or protection efficacy without additional data and a validated model. Keep those evidence layers separate.
5. Create a **new integration addendum** with columns for each claim, raw-file path, input identity, method/version, validation check, source figure/table and the distinct Phase IV finding it complements. Preserve the current Phase IV-only results as the comparison baseline.

The prepared receptor files in `Step_Outputs/Phase3/StepA/` are `Phase3A_TLR2_6NIG_receptor.pdb`, `Phase3A_TLR4_8WTA_receptor.pdb` and `Phase3A_Vaccine_ligand.pdb`. The old handoff recorded 6NIG's chimeric non-human tail and bound agonist as preparation traps, and 8WTA's retained MD-2 lipid-A chains as a deliberate choice. **Verify these properties from the actual prepared structures and deposition records before any new docking interpretation**; do not treat this handoff as validation of a later docking run.

## 5. Files, safeguards and manuscript work

**Locked trees:** `Phase 1/`, `Phase 2/`, `Phase 3/`, `Step_Outputs/Phase1/`, `Step_Outputs/Phase2/`, `Step_Outputs/Phase3/`. Read them as needed; do not edit or rerun them as part of this continuation. Keep future integration artifacts in a new named location. `Phase 5/` is separate work; this handoff does not certify its status or values. The old broad handoff quoted an 8.7-hour terminal half-life from the standalone draft, while `FINDINGS.md` records a correction; **do not reuse that old number without independently checking the Phase V model and outputs**.

**Phase IV source map:**

- Code: `Phase 4/STEP A`–`STEP G`, with shared helpers in `Phase 4/_common/`.
- Row-level output: `Step_Outputs/Phase4/StepA`–`StepF` CSVs and methodology notes.
- Final selected-epitope dossier: `Step_Outputs/Phase4/StepG/Phase4G_Immunogenicity_Dossier_20260921_2022.csv`.
- Numeric source index: `Step_Outputs/Phase4/StepG/Phase4G_NumbersLedger_20260921_2022.csv`.
- Historical process record: `Step_Outputs/Phase4/PHASE4_PROGRESS_LOG.md` and `_run_logs/`; the intermediate `HALT-FOR-OPUS.md` contains superseded counts.
- Integrity baseline: `Step_Outputs/Phase4/_rule1_baseline_sha256.txt` (1,810 files).
- Current human-readable results and paper package: `PHASE_IV_RESULTS_DETAILED.md` and `PAPER_UPDATES_PHASE_IV.md`.

**Operational pitfalls from the recorded run:** parse IEDB values by named units, not column position; failed API calls must not be treated as “allele unsupported”; keep the Phase 1F and Phase 4B allele universes separate; do not interpret an extreme over many seam windows without a margin and an appropriate reference; treat `N/A` and `UNRESOLVED` differently. These are documented in the progress log. The prior toolchain used volume symlinks and several conda environments, but their current state must be checked before rerunning any Phase IV code; no rerun was required for this consolidation.

`Phase 4/STEP G/Phase4G_integratedReport.py` is a historical report generator. Running it would overwrite the two Step G redirect files with older prose, including claims superseded by the current results narrative. Review and update that generator before using it to regenerate a research report.

`PAPER_UPDATES_PHASE_IV.md` contains the replacement Section IV methods, manuscript-ready results/discussion and the full correction register. The reviewed 42-page paper is a local working PDF (SHA-256 `14d0b81a0950d4b8ef2ceb98c618003e9fa95336f413b1c0d4bc9d7b4fad74ef`) with no editable source found in the working research folder; it is excluded from the Phase IV GitHub publication set pending the Section III provenance review. The PDF itself was not rebuilt. The proposed old Section IV on pages 40–42 should be replaced in the author's editable manuscript. The old Section III, abstract and conclusion claims require the separate provenance audit described above before they are integrated with the revised Section IV.
