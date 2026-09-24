# Phase IV consolidated documents — quality check

**Checked:** 2026-09-24  
**Repository:** `/Volumes/Extended SSD/ChimericVaccine/Research`  
**Scope:** `PHASE_IV_RESULTS_DETAILED.md`, `Handoff.md`, `PAPER_UPDATES_PHASE_IV.md`, their three redirect documents, and `MANUSCRIPT_CORRECTIONS_MASTER.md`. The later `Phase 4/README.md` was checked for path and scope consistency before publication selection.

## Verdict

**PASS for consistency with the saved Phase IV evidence, within the checks below.** The documents distinguish prediction from biological observation and keep Phase III/docking results outside the completed Phase I/II-only branch. This is a document and saved-output audit, not a rerun of remote prediction services, an experimental validation, or an independent validation of every source-ledger value.

## Checks performed

| Check | Result |
|---|---|
| Final construct identity | The Phase 1G FASTA is 570 amino acids; its sequence MD5 prefix is `4fce119e`. All 31 dossier peptides match the FASTA at their recorded positions. |
| Final dossier and comparison pool | Dossier: 31 rows, 153 columns, no blank or `UNRESOLVED` cells; 10 MHC-I, 11 MHC-II, 10 B-cell. Phase 4B comparison pool: 105 rows, comprising 37/29/39 by class. |
| Selected-peptide appendix | All 31 appendix rows were compared field by field with the dossier for identity, binding, T-cell scores, processing, B-cell screening, fold assessment and conservancy where applicable. All matched at the displayed precision. |
| Coverage and binding gaps | Phase 4E per-locus values and Phase 4B calls support the reported A/B/C and class-II coverage values. No selected MHC-I peptide binds `A*11:01` or `A*33:03` under either saved predictor universe; 0/37 comparison peptides bind them under Phase IV IEDB-EL. |
| Secondary verdicts | The dossier records 10/10 MHC-I strict failures, 10/11 MHC-II strict passes and 4/10 B-cell strict passes. Phase 4C best-allele BigMHC-IM values reach 0.5 in 3/10 and 0.7 in 0/10. Four MHC-I peptides also exceed the supplementary 500 nM IC₅₀ cut. |
| Processing and seam analysis | Phase 4D records 10/10 selected MHC-I peptides meeting its C-terminal rule and five internal-cleavage watchlist members. Phase 4E records 28 contested class-I epitope–allele pairs at the 5× margin; 293/38,924 versus 45/10,804 strong-binding window–allele pairs (rounded rates 0.75% and 0.42%; descriptive ratio 1.81). Of 193 class-II contested rows, zero involve GPGPG. |
| B-cell and conservancy | Phase 4F records 1 High, 6 Medium and 3 Deprioritized B-cell tiers; all ten SEMA-3D calls are `NO`, and five map to globally undetermined native folds. Seven antigen mean conservancy values match at their reported precision; all 31 selected-peptide regression flags are `YES`. |
| Source trace and formatting | All 155 numbers-ledger rows name an existing source file, after resolving the seven Phase 1F rows' explanatory suffix. All 38 concrete relative paths checked in the seven documents exist. Markdown table column counts are consistent in the four main documents. |
| Locked data integrity | All 1,810 Phase I/II files in the saved SHA-256 baseline still match. Phase I, II and III source/output files were not edited during this quality check. |

## Corrections made during this check

- Replaced `B5R/B6R` shorthand with the frozen `Mpox_B5R` pipeline label and explicitly left the manuscript's B6R name pending a source-accession crosswalk.
- Changed the conservancy wording from “exactly” to “at reported precision.”
- Named the 45-residue adjuvant domain as the limited seam comparison reference, rather than implying a broad natural-window control.
- Expanded two shortened Step G evidence paths and warned in the handoff that the historical Step G generator would recreate superseded prose if rerun.

## Limits that must remain visible in a paper

1. Binding, processing, T-cell scores, SEMA calls and population coverage are computational predictions. The Phase IV files contain no measured immune response or protection result.
2. The AFND Singapore Riau Malay panel (n=132) is a Southeast Asian proxy, not a Philippine cohort. The project's projected-coverage formula is its own calculation and must not be described as a direct run of the IEDB Population Coverage Tool. [IEDB population-coverage help](https://tools.iedb.org/population/help/) describes the tool's own HLA genotypic-frequency input and coverage output.
3. `Mpox_B5R` versus the manuscript's B6R name requires an accession-level crosswalk before normalization. The saved values remain reproducible under their pipeline label.
4. The manuscript PDF reports docking, MD and C-ImmSim in its old Section III, but their underlying run provenance has not been established by this Phase IV audit. No Phase III result was incorporated into these documents.
5. The old Step G generator still contains superseded prose. The two redirect documents are current, but rerunning that generator would overwrite them until the generator is revised.

## GitHub preflight, without a push

- `origin` is `https://github.com/kedddot/ChimericVaccine-Mpox-HIV.git`; local `main` and remote `main` both point to `bc91c2c0485145eeb753c3ce810da38a4375887d` at this check.
- The consolidated Markdown documents, revised methodology, `Phase 4/` code and `Step_Outputs/Phase4/` are currently untracked. Phase IV outputs contained 1,106 files totaling about 27 MB at the check; one disposable `.DS_Store` was subsequently removed. No individual file above 90 MB was found in that tree.
- The worktree also contains unrelated untracked Phase V material. The pre-existing tracked deletion of `CHIMERIC VACCINE.pdf` was reviewed: the old 27-page draft has a placeholder abstract and reports 32 selected epitopes, whereas the final dossier has 31. The publication selection includes removing that obsolete PDF from the current GitHub branch; an identical copy remains in the parent project folder. The newer 42-page working PDF is excluded pending its old Section III raw-run provenance check. A broad `git add -A` would include material outside this Phase IV task. A GitHub commit should use an explicit file list and preserve the untouched locked Phase I/II/III trees.
- **No commit, pull request or push was made by this check.**
