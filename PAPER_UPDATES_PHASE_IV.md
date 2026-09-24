# Paper updates: revised Phase IV immunogenicity analysis

**Prepared and consolidated:** 2026-09-24  
**Manuscript reviewed:** `DSTF_MCS_TeamMPOXHIV_Multi-Epitope-Immunoinformatics-2.pdf` (42-page local working draft; SHA-256 `14d0b81a0950d4b8ef2ceb98c618003e9fa95336f413b1c0d4bc9d7b4fad74ef`; excluded from the Phase IV GitHub publication set pending the separate Section III provenance review)  
**Controlling methods:** `REVISED_METHODOLOGY_MpoxHIV.md` §IV, as implemented by `Phase 4/` when the two differ  
**Evidence boundary:** Phase I and II outputs only. No Phase III, molecular docking, MD, mRNA delivery, pharmacokinetic or C-ImmSim outputs are used below.

## Editorial placement and scope

The PDF's pages 40–42 contain an **older proposed** Section IV on daughter-construct splicing, repeat docking/MD and C-ImmSim comparisons. Replace that proposal with the current IEDB-centered Section IV below in the editable manuscript source. Add the Results and Discussion text below after the existing Phase I/II results, with new tables/figure references assigned during typesetting. The supplied PDF is a rendered artifact; no editable `.docx`, `.tex` or other paper source was found in the working research folder, so this file is the reviewable paper update rather than an altered PDF.

The existing PDF's old Section III and abstract/conclusion discuss docking, MD and C-ImmSim as completed. Their raw run provenance is unresolved within this repository; those claims require a separate evidence check. Do not use the proposed Section IV's splicing/independent simulations as if they were performed in the completed Phase IV branch. The completed Phase IV cannot substantiate antibody titres, cytokine time courses, clonal dynamics, receptor activation or protective efficacy.

## Replacement Section IV Methods — ready for manuscript editing

### IV. Advanced immunogenicity prediction

The final 570-amino-acid construct (`Vax_Final_4fce119e`) contained 31 selected epitopes: 10 MHC class I, 11 MHC class II and 10 linear B-cell epitopes. Fourteen epitopes originated from HIV antigens and 17 from Mpox antigens. The 105 candidates surviving Phase I population-coverage scoring (37 MHC-I, 29 MHC-II and 39 B-cell) served as a comparison pool. The 194 sequences mentioned elsewhere in the manuscript are source protein sequences, not an epitope pool. Construct membership, pathogen origin, source antigen, position and Phase I scores were joined by peptide sequence in the Phase IV dossier.

MHC-I and MHC-II candidates were scored through the IEDB MHC tools API using NetMHCpan-EL and NetMHCIIpan-EL percentile ranks, with consensus outputs supplying supplementary predicted IC₅₀ values where the allele was supported. Phase I selection thresholds remained the primary gate (MHC-I rank ≤1%; MHC-II rank ≤10%); IC₅₀ did not retroactively select or exclude frozen epitopes. Binding breadth was the number of alleles meeting the class-specific rank threshold, out of 74 class I or 33 class II alleles in the AFND Singapore Riau Malay panel. The population comprises 132 sampled individuals and is used as a Southeast Asian proxy, not as a Philippine population. Projected coverage was calculated from the union of binding alleles per locus using `1 − (1 − s)²`, where `s` is the capped sum of bound-allele frequencies at that locus, followed by complement multiplication across loci. Phase I and Phase IV predictor universes were calculated separately, and their per-locus differences were reported rather than pooled.

MHC-I T-cell recognition was assessed with BigMHC-IM, PRIME 2.0/MixMHCpred and an exploratory local implementation of the Calis model. The project's BigMHC-IM cutoff of 0.5 was retained; 0.70 was tabulated as a strict secondary stratum. These are model scores, not measured response probabilities. The proposed MixMHCpred/TCGA TCR-contact procedure in the standalone methodology was not executable as written: MixMHCpred predicts MHC presentation and TCGA is a cancer-genomics resource. A separate human cross-reactivity review combined Phase I BLASTP results against reviewed human Swiss-Prot with a fresh exact 8-mer screen. Neither computational layer establishes biological safety.

The complete construct was evaluated with IEDB's MHC-I processing predictor for 8–11-mers. A selected MHC-I epitope met the pipeline's liberation criterion when its C-terminal cleavage score was at or above the construct-wide median; stronger internal sites were recorded as a watchlist. Junction-spanning windows were screened across MHC-I 8–11-mers and MHC-II 12–20-mers. A seam window was called *contested* at a given allele only when it crossed at least three residues on both sides of the seam, passed the processing criterion, and had a predicted affinity or EL-rank advantage of at least `m`-fold over the best real epitope at that allele; a class II call additionally required its predicted core to span the seam. The primary margin was `m = 5`, with 1, 2 and 10 as sensitivity values. Class I window binding was compared with unselected windows from the 45-residue adjuvant domain; class II EL-rank prevalence was compared with its random-peptide percentile reference. A contested prediction was not treated as measured loss of a genuine epitope.

B-cell epitopes retained Phase I's High, Medium and Deprioritized BepiPred tiers. Existing SEMA-3D predictions on native antigen folds were used for conformational corroboration; SEMA was not rerun on the whole-construct model. Native-fold relative solvent accessibility was interpreted in light of each fold's global confidence and did not include oligomer interfaces or HIV Env glycans. Conservancy was recomputed by exact peptide matching in the Phase I variant pools. Full-pool antigen means used the 30,268-candidate `Raw_Conservancy` table; the 1,238-row `Min_50pct` table was used only to check selection after the ≥50% gate. Variant pools contained 21–30 isolates per antigen, so these estimates describe the sampled within-clade diversity.

### Why these methods differ from the standalone draft

The standalone §IV lists a ten-allele “Philippine” panel, IC₅₀ selection cuts, a 0.70 BigMHC threshold, one BepiPred high-priority cut, an incomplete population-coverage equation, an IEDB Conservancy Tool run, and a TCR-contact database combination that did not match the implemented pipeline. The frozen Phase I/II construct and completed Phase IV CSVs take precedence for reporting. The decision-by-decision correction register is consolidated below. The older Step G corrections file is now a pointer to this package. The IEDB documentation distinguishes EL rank from BA/IC₅₀ measures; keep these labels explicit in tables.

## Phase IV Results — ready for manuscript editing

The 31-epitope construct comprised 10 MHC-I, 11 MHC-II and 10 B-cell epitopes; the comparison pool comprised 37, 29 and 39, respectively. Median binding breadth was 7 of 74 alleles (range 1–40) for pool MHC-I candidates and 6 of 33 (range 1–29) for pool MHC-II candidates. Construct medians were 6 of 74 and 6 of 33. All 66 HLA-restricted pool candidates retained their Phase I gate status when rescored under the Phase IV analysis. Predicted IC₅₀ at the best consensus-supported allele was available for 37/37 MHC-I and 26/29 MHC-II candidates. The three missing class II values were recorded as unavailable, not interpreted as negative binders. Phase I MHC-II BA/consensus ranks and Phase IV EL ranks represent different prediction scales.

The clearest coverage gaps concerned `HLA-A*11:01` and `HLA-A*33:03`, with panel allele frequencies 0.177 and 0.109: no construct MHC-I epitope bound either allele under either predictor, and no candidate among the 37 MHC-I comparison peptides bound them **under Phase IV IEDB-EL**. Neither allele was included in Phase I's nine-allele discovery screen. This mismatch is a plausible contributor to the gap, although a counterfactual screen was not performed. `HLA-C*08:01` (frequency 0.187) was a gap only under IEDB-EL and is therefore predictor-dependent. Projected locus A coverage was 60.31% with IEDB-EL versus 75.99% with Phase I's MHCflurry/MHCnuggets universe; locus B was 89.37% versus 98.96%, and locus C 70.19% versus 98.18%. Class II loci differed less: DRB1 99.54% under both, DQB1 96.00% versus 95.25%, and DPB1 86.97% versus 84.08%. Combined multi-locus coverage saturated in the model and should be reported with these per-locus values, not alone as a broad population-protection claim.

For the HIV-derived MHC-I subset (n=5), IEDB-EL projected 97.52% coverage; the Mpox-derived MHC-I subset (n=5) reached 83.87%. The corresponding Phase I predictor values were 99.96% and 94.95%. The Mpox A35R epitope `NKRKRVIGL` bound only 2 of 74 class I alleles, had BigMHC-IM score 0.0344 and PRIME percentile rank 4.386, and added 0.0000 percentage points to the construct-wide allele-union coverage. Per-pathogen combined MHC-I/II coverage was 99.98% for HIV and >99.99% for Mpox under IEDB-EL, which does not establish balanced immune response or absence of immunodominance.

For the ten selected MHC-I epitopes, the BigMHC-IM best-allele score met the project 0.5 cutoff in 3/10 and the draft methodology's 0.70 cut in 0/10. The strict secondary stratum therefore failed all ten MHC-I epitopes, although the Phase I primary gate remained pass; four also failed the secondary IC₅₀ ≤500 nM condition. PRIME percentile rank correlated with the Phase IV EL rank (Pearson r=0.916, p=0.0002, n=10); BigMHC-IM's correlation was not distinguishable from zero (r=−0.298, p=0.402, n=10). The Calis implementation was exploratory (r=0.097, p=0.790, n=10). BigMHC-IM and PRIME selected the same best allele for only 3/10 peptides, so their best-allele scores cannot be compared as matched observations. No construct epitope contained an exact eight-residue match in the reviewed human Swiss-Prot set (0/31); that screen is limited by its reference database and length rule.

All ten selected MHC-I epitopes met the pipeline's C-terminal liberation criterion, but five (`GPKEPFRDY`, `HHFNCRGEF`, `HWTTYMDTF`, `RFALNPGLL` and `VYSTCTVPTM`) had a stronger predicted internal cleavage site. At the primary 5-fold seam margin, 28 class I epitope–allele pairs were contested and the modeled whole-construct uncontested coverage was 99.88%. Among 31 qualifying class I seam windows, 29 involved AAY and 21 ended in AAY. Strong-binding class I window–allele pairs occurred in 0.75% of deep seam predictions (n=38,924), versus 0.42% of unselected adjuvant-domain windows (n=10,804), a descriptive 1.81-fold rate ratio against this limited reference. This suggests a linker-associated prediction signal, especially at AAY junctions; it does not show that genuine pathogen epitopes are absent from the cell surface.

No GPGPG seam appeared in the 193 class II 5-fold contested rows. Across all 82,203 deep class II seam window–allele pairs, 7.17% had EL rank ≤10, below the 10% random-peptide percentile reference. However, KK and AAY seams accounted for 108 and 83 of the 193 contested rows, with two at the adjuvant junction. Because the complete 5-fold criterion was not calibrated against a matched random-window null, the remaining class II calls are unresolved predictions. This distinction should accompany the favorable GPGPG observation.

Among the ten B-cell epitopes, one was High, six Medium and three Deprioritized by the Phase I BepiPred tiers; four met the draft methodology's alternate high-priority cut. All ten had SEMA-3D status `NO` at the chosen 50% overlap threshold, meaning screened without conformational corroboration rather than disproven. Five were mapped to native folds judged globally undetermined, which limits surface-exposure interpretation. Across the full 30,268-candidate pool, mean conservancy reproduced the prior seven antigen values: 4.54–9.97% for HIV targets and 37.80–60.61% for Mpox targets. All 31 selected epitope conservancy values reproduced from the variant pools. Selected epitopes were more conserved than the full pool by construction (mean 75.5% HIV, 99.7% Mpox); the small within-clade reference sets do not support a global-strain claim.

### Suggested table for the manuscript

| Endpoint | IEDB-EL / Phase IV | Phase I predictor | Interpretation |
|---|---:|---:|---|
| Locus A coverage | 60.31% | 75.99% | A*11:01 and A*33:03 unbound in both |
| Locus B coverage | 89.37% | 98.96% | Predictor spread 9.59 pp |
| Locus C coverage | 70.19% | 98.18% | Largest spread, 27.99 pp |
| DRB1 coverage | 99.54% | 99.54% | Same projected value |
| DQB1 coverage | 96.00% | 95.25% | Different predictor universe |
| DPB1 coverage | 86.97% | 84.08% | Different predictor universe |
| HIV MHC-I subset | 97.52% | 99.96% | Five construct epitopes |
| Mpox MHC-I subset | 83.87% | 94.95% | Five construct epitopes |

*Table note:* model-derived coverage from AFND Singapore Riau Malay frequencies (n=132); allele universes and thresholds must be stated in the caption. These are projected binding coverage values, not observed vaccine efficacy or Philippine population coverage. Values are from `Step_Outputs/Phase4/StepE/Phase4E_MethodDivergence_20260921_1949.csv` and `Step_Outputs/Phase4/StepE/Phase4E_ConstructCoverage_20260921_1949.csv`.

## Discussion and limitations — ready for manuscript editing

This Phase IV analysis adds allele-specific and linker-processing resolution to the frozen construct. Its strongest finding is the lack of predicted MHC-I binding to two common reference-panel HLA-A alleles across both algorithms, alongside substantial class I predictor divergence. The narrow Phase I discovery panel is a plausible explanation, but this study did not test an alternative construct or estimate clinical response in a Philippine cohort. The AAY seam signal is mechanistically plausible but based on predicted affinity/processing and a limited natural-window reference; a contested call is not measured suppression of a pathogen epitope. GPGPG class II junctions generated no primary-margin contests, while other class II seams did, so the latter cannot be declared risk-free.

The static MHC and B-cell predictors do not reproduce C-ImmSim's temporal antibody, cytokine or memory-cell outputs. The SEMA-3D non-corroborations and uncertain native folds further limit B-cell interpretation. Experimental pMHC presentation, T-cell activation, antibody binding, human cross-reactivity and Philippine HLA sampling would be needed to test the computational predictions. Phase III docking or MD claims should be discussed alongside these results only after the manuscript's underlying runs and construct identity are verified; they cannot fill the missing immune-response evidence by themselves.

## Specific existing-paper changes

1. Replace pages 40–42's proposed daughter-construct Section IV. No such daughter constructs, repeat docking/MD, ten-seed C-ImmSim comparison or inferential t-test were part of the completed revised Phase IV.
2. In the abstract and conclusion, add only a concise statement that the final construct underwent *in silico* allele-specific immunogenicity and linker-seam assessment. Do not say Phase IV demonstrated immune protection, balanced HIV/Mpox responses or absence of immunodominance.
3. Change any “Philippine HLA coverage” attributed to this Phase IV calculation to “AFND Singapore Riau Malay Southeast Asian proxy coverage.” Give locus-level ranges and the two HLA-A gaps.
4. Use 31 final construct epitopes and 105 comparison candidates; reserve 194 for source protein sequences. Keep class-specific breadth denominators 74 and 33.
5. Report the final ten B-cell SEMA statuses (0 corroborated, 10 screened without corroboration, 0 unscreened). Older nine-epitope tallies belong to a superseded construct.
6. Do not describe BigMHC-IM score as a calibrated in-vivo response probability or as an input that integrates pMHC stability. Treat IC₅₀ and EL rank as different outputs.
7. Resolve the PDF's old Section III provenance separately. The PDF reports docking, MD and C-ImmSim, whereas this repository exposes receptor preparation only. Retain those claims only after their underlying runs have been audited.

## Source trail

- Single current results narrative: `PHASE_IV_RESULTS_DETAILED.md` (including the former Step G report’s distinct evidence and 31 selected-peptide summary rows)
- Final row-level data and number ledger: `Step_Outputs/Phase4/StepG/Phase4G_Immunogenicity_Dossier_20260921_2022.csv`; `Step_Outputs/Phase4/StepG/Phase4G_NumbersLedger_20260921_2022.csv`
- Consolidated Section IV correction register: the section below; the former Step G corrections file points here
- Phase IV methods: `Phase 4/STEP A` through `STEP G`; `Step_Outputs/Phase4/StepA` through `StepG` methodology notes and outputs
- External method guidance: [IEDB MHC-I help](https://tools.iedb.org/mhci/help/), [IEDB population-coverage help](https://tools.iedb.org/population/help/), [IEDB tools API](https://tools.iedb.org/main/tools-api/). These explain tool outputs; all numerical findings above come from the project's files.

## Consolidated Section IV correction register

This register merges the former Step G corrections document with the manuscript-ready update above. “Implemented” means the saved Phase IV pipeline or frozen Phase I data used the method; it does **not** imply independent experimental validation. Historical API and tool failures are dated observations, not claims about their present availability. The full execution record remains in the per-step methodology notes and progress log.

### A. Eight direct differences from the standalone revised methodology

| ID | Standalone §IV wording | Wording supported by the executed pipeline |
|---|---|---|
| A1 | Use IC₅₀ ≤500 nM (MHC-I) or ≤1000 nM (MHC-II) to select candidates | The frozen primary gates use EL percentile rank ≤1% for class I and ≤10% for class II. Consensus IC₅₀ is an additional descriptor and part of the separately labelled strict-secondary stratum. Its allele coverage is incomplete. |
| A2 | BigMHC-IM ≥0.70 indicates high in-vivo T-cell response likelihood | The project used a 0.5 score convention in Phase I and reports 0.70 only as a strict-secondary cut. In the final MHC-I cassette, 3/10 reach 0.5 and 0/10 reach 0.70 by the saved Step C best-allele scores. Neither value is a measured response probability. |
| A3 | One BepiPred high-priority rule: mean ≥0.50 and ≥75% above threshold | Phase I stored High (mean ≥0.60 and ≥75%), Medium (≥0.50 and ≥50%), and Deprioritized (≥0.45 and ≥37.5%) tiers. The alternate rule is reported separately; tiers did not remove frozen epitopes. |
| A4 | Ten “Philippine” HLA alleles, including `A*24:23` | Coverage uses the AFND Singapore Riau Malay frequency table: 107 entries (74 class I, 33 class II), n=132; a Southeast Asian proxy, not a Philippine sample. The listed ten alleles/frequencies were not the executed panel. |
| A5 | Coverage `1 − ∏(1 − freq[allele])` | Code sums bound-allele frequencies within each locus (capped at 1), applies diploid locus coverage `1 − (1 − s)²`, then multiplies uncovered fractions across loci. This is a model with frequency and cross-locus assumptions; report per-locus values. |
| A6 | “MixMHCpred + TCGA Contact Database” predicts TCR contacts | The executed MHC-I score ladder used BigMHC-IM, PRIME 2.0/MixMHCpred 3.0, and an exploratory Calis implementation. The draft’s named combination was not a validated TCR-contact procedure in this pipeline; positions 4–6 were used as a defined descriptive TCR-facing segment, not measured contacts. |
| A7 | A ≥70% BLASTP identity exclusion proves lack of human self-reactivity | The saved analysis reports Phase I BLASTP alongside a fresh exact 8-mer comparison to reviewed human Swiss-Prot. It found 0/31 exact 8-mer matches in that reference. Neither layer establishes biological safety, and the identity rule alone was shown in Phase I to miss short self-fragments. |
| A8 | SEMA 2.0 on the construct validates B-cell conformational epitopes | Phase IV reused Phase I's local SEMA-3D calls on native antigen folds after construct selection. All 10 final B-cell epitopes were screened and returned `NO` under the chosen overlap rule; none was unscreened. `NO` means no corroboration at that threshold, not refutation. |

### B. Further corrections and process qualifications

| ID | Item to correct or qualify | Evidence-constrained action |
|---|---|---|
| B1 | “194 non-redundant sequences” as an epitope count | There are 194 Phase 1A variant FASTA files representing source protein sequences. Phase IV’s analysis units are 31 selected epitopes and 105 Phase 1F comparison candidates. |
| B2 | `Mpox_B5R` versus manuscript `B6R` | Keep the frozen `Mpox_B5R` label in reproducibility tables and provide a source-accession crosswalk before replacing it with the manuscript’s antigen name. The two strings must not be silently treated as equivalent. |
| B3 | BigMHC integrates IC₅₀, pMHC stability and TCR-interface features | Remove this input list from the methods unless directly supported by the version of BigMHC run here. The executed Step C score should be described as a model output for peptide–HLA pairs. |
| B4 | MERCI described as a toxin/allergen “database” run at Phase IV | MERCI is part of the earlier ToxinPred2 motif-screening workflow; construct toxicity/allergenicity was evaluated in Phase II, not rerun as a new Phase IV result. |
| B5 | “IEDB Conservancy Tool” across “global strains” | Phase I/IV used exact peptide matching against seven local Phase 1C variant pools (21–30 isolates per antigen). Describe sampled within-clade conservancy, not a global estimate or an IEDB Conservancy Tool execution. |
| B6 | Reproducing seven antigen means from `Min_50pct` | Use `Raw_Conservancy` over the full 30,268-candidate pool. `Min_50pct` is the selected subset (1,238 rows) and changes the denominator. The seven stored means match the full-pool recomputation. |
| B7 | Fixed DR/DQ list and “top ten” allele methods | Report the actual NetMHCpan-EL/NetMHCIIpan-EL and consensus API calls over the 107-entry panel. DQA1 and DPA1 frequencies were not in this table. |
| B8 | IC₅₀ parsing incident | An unlabeled MHC-II consensus column was briefly read as nM; negative values and reversed rank correlation exposed the error. It was removed. Final IC₅₀ values use labelled nM columns only; 37/37 MHC-I and 26/29 MHC-II pool candidates have a value at their best consensus-supported allele. Preserve `N/A` for the remaining three. |
| B9 | Treating Phase I and Phase IV MHC-II ranks as identical | Phase I’s stored `recommended` BA/consensus rank and Phase IV’s EL rank are different scales. Show both only with method labels; the 0/66 gate-status-flip check is an implementation comparison, not proof of equivalent predictors. |
| B10 | IEDB `immunogenicity/` endpoint explanation | During the recorded **2026-09-21** run, repeated calls to that URL returned 403 while other IEDB endpoints answered. This explains why the Calis implementation was used then. Do not claim the endpoint is permanently unavailable or that current service status was checked in this consolidation. |
| B11 | Calis-model evidence strength | The local Calis implementation reconstructed an enrichment table from published training data; its training-set self-AUC of 0.62 is not an independent validation. No measured immunogenicity labels are available for the ten selected peptides. Keep the score exploratory. |
| B12 | Older SEMA tally on a nine-B-cell construct | That tally belongs to a superseded construct. The final 570-aa construct has ten B-cell epitopes: 0 corroborated, 10 screened without corroboration, 0 unscreened under the saved rule. |
| B13 | `SEMA NO` interpreted as structural refutation | Define `NO` as screened below the selected overlap threshold. It does not establish absence of a native antibody epitope. `UNSCREENED` is a distinct missing-assessment state. |
| B14 | Breadth reported “out of 107” for each class | Denominators are 74 class-I entries and 33 class-II entries; pool medians are 7 and 6 bound entries, respectively. |
| B15 | A ≥90% target claimed for every locus | Combined model coverage is saturated, but locus A is 60.31–75.99% and locus C 70.19–98.18% across predictors. State the per-locus values and model/proxy limitation. |
| B16 | Linker/seam analysis missing from methods | Include the pipeline's processing criterion, deep-junction and core requirements, 1/2/5/10× sensitivity curve, and reference limitations. Class-I AAY windows show a descriptive strong-binder rate ratio of 1.81 against a 45-aa adjuvant-domain reference; no GPGPG window appears in the 193 class-II 5× contested rows. Neither finding measures antigen presentation or immunodominance. |
| B17 | Phase IV coverage called “Philippine” | Say “projected coverage under AFND Singapore Riau Malay frequencies, used as a Southeast Asian proxy.” Do not generalize it to a measured Philippine population. |
| B18 | Review-note numerical mix-ups | Five, not four, selected MHC-I epitopes have a stronger predicted internal cleavage site; four of those overlap an AAY contested-seam list. Of 31 qualifying class-I windows, 29 involve AAY and 21 end in AAY. Locus B is class I. |

### Author checks before publication

- Resolve the manuscript's old Section III docking/MD/C-ImmSim raw-data provenance independently. The PDF alone does not validate those numerical results for a combined Phase IV analysis.
- Verify the `Mpox_B5R`/B6R source-accession crosswalk before normalizing antigen labels across tables and prose.
- Distinguish predicted EL rank, predicted IC₅₀, BigMHC-IM score, coverage projection, predicted processing, and measured biological response wherever they appear. None of the first five is a direct efficacy result.
- Typeset the replacement Section IV from this package in the author’s editable source; the current PDF was not modified by this consolidation.
