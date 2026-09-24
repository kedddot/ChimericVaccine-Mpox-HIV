# PHASE IV — DETAILED RESULTS
### IEDB-Based Immunogenicity Prediction · `Vax_Final_4fce119e`

**For:** researchers writing up Section IV, or building analyses on top of it.
**Compiled:** 2026-09-22 · **Consolidated and checked:** 2026-09-24 · **Status:** Phase IV complete for the Phase I/II-only branch (Steps A–G)

This is the **single current narrative of Phase IV results**. It combines the former Step G report and detailed summary. All numerical findings use Phase I/II inputs and Phase IV outputs; no Phase III or docking input was used. The 31-row, 153-column final dossier, 105-row comparison pool and 155-entry source ledger were re-read; no dossier cell was blank or `UNRESOLVED`. All 1,810 Phase I/II baseline files matched their saved SHA-256 hashes on 2026-09-24. This review checked saved outputs and calculations; it did not rerun remote prediction services or validate the models experimentally. The ledger records reported numbers and source filenames, but a filled source field alone is not independent validation of a value.

---

## 0. WHAT THIS ANALYSIS IS AND IS NOT

The revised Section IV uses IEDB and related predictors to assess **epitope-level binding and presentation proxies**. It does not simulate an immune system or directly measure immunogenicity. The older proposed daughter-construct Section IV appended to the PDF was not executed in this branch.

**It does not produce:** antibody titres, cytokine concentrations, cell-population counts, or any
time-course. Claims of that kind cannot be sourced from Phase IV.

**It reports:** predicted per-allele binding, model-based population coverage, MHC-I T-cell-recognition scores, processing scores, and an epitope–linker seam analysis. These quantities have different evidence strength and are not interchangeable.

---

## 1. SUBSTRATE

| Set | n | Composition |
|---|---|---|
| Construct epitopes | 31 | 10 MHC-I, 11 MHC-II, 10 B-cell |
| By pathogen | 31 | 14 HIV-derived, 17 Mpox-derived |
| Comparison pool | 105 | 37 MHC-I, 29 MHC-II, 39 B-cell |
| Source antigens | 7 | HIV gp120, gp41, p17, p24; Mpox A35R, `Mpox_B5R` (pipeline label), L1R. The manuscript's B6R name needs a source-accession crosswalk. |

**Allele panel:** 107 entries from the AFND *Singapore Riau Malay* dataset (n=132 individuals), used as a Southeast Asian proxy, **not** a Philippine panel: 74 class I and 33 class II. The population-coverage formula assumes a diploid locus and combines locus coverages by a complement product; it is a model projection, not observed response. Breadth denominators are 74 and 33. Class II includes DRB1, DQB1 and DPB1 beta-chain entries; DQA1 and DPA1 are absent.

---

## 1A. METHODS AND DECISION RULES USED FOR THESE RESULTS

- The final construct is frozen. Phase IV scored the 105 Phase 1F comparison candidates and identified the 31 members of the final construct by sequence. NetMHCpan-EL and NetMHCIIpan-EL percentile ranks supplied the primary class-specific binding calls (≤1% for MHC-I, ≤10% for MHC-II); IEDB consensus IC₅₀ values are supplementary and unavailable for some alleles. Phase I's selection gate and the revised methodology's stricter secondary stratum are both retained in the dossier, rather than retroactively excluding any epitope.
- Projected locus coverage was calculated as `1 − (1 − s)²`, with `s` the capped frequency sum of alleles bound by at least one epitope at that locus. The complement of uncovered fractions across loci gives the combined model value. Phase 1F MHCflurry/MHCnuggets and Phase 4B IEDB-EL allele sets were never mixed within one coverage calculation. This calculation assumes the frequency model used by the code and should not be described as a direct IEDB Population Coverage Tool run or as a Philippine cohort estimate.
- MHC-I recognition used BigMHC-IM, PRIME 2.0 with MixMHCpred, and an exploratory Calis-model implementation. Scores and best alleles are model outputs, not observed T-cell activation. BLASTP and an exact 8-mer comparison used the reviewed human Swiss-Prot reference; neither is a clinical safety test.
- IEDB processing predictions supplied the pipeline's C-terminal liberation rule: a selected MHC-I epitope qualifies when its terminal cleavage score is at least the construct-wide median. A separate internal-cleavage watchlist does not change that verdict. A seam window is *contested* only when it is deep (≥3 residues on each side), passes the processing criterion, and beats the best real epitope at the same allele by at least the stated margin; class II also requires a core spanning the seam. The primary margin is 5×, with 1×, 2× and 10× sensitivity points.
- Phase I BepiPred tiers and stored SEMA-3D calls were reused. Native antigen folds, rather than the low-confidence full-construct fold, supplied relative solvent accessibility. Conservancy was recomputed on the Phase I variant pools; the full-pool means use `Raw_Conservancy`, not the selected `Min_50pct` subset.

---

## 2. FINDING 1 — COVERAGE GAPS (strongest evidence; no null model required)

Two high-frequency class I alleles are bound by **zero** construct MHC-I epitopes under **both**
predictors:

| Allele | Frequency | Epitopes binding it (IEDB-EL) | In the 37-candidate pool |
|---|---|---|---|
| `HLA-A*11:01` | 0.177 | 0 | 0 |
| `HLA-A*33:03` | 0.109 | 0 | 0 |

Combined allelic mass of the two: **0.286** (35.0% of locus-A frequency mass).

Under the **Phase IV IEDB-EL** call, none of the 37 MHC-I comparison candidates binds either allele. The final construct also has zero binders under both predictor universes. The gap is present in the saved candidate pool; the nine-allele Phase I discovery panel omitted both alleles, a plausible contributor. The counterfactual effect of screening a different pool was not tested.

`HLA-C*08:01` (0.187) is a gap under IEDB-EL only; MHCflurry finds `TMGAASITL` binding it, and 5
pool candidates bind it. Report it as predictor-dependent, not as an unambiguous gap.

Within this reference panel, `A*11:01` is the highest-frequency **HLA-A** allele. `C*08:01` has a higher frequency across the full class-I panel (0.187), and its gap is predictor-dependent.

---

## 3. FINDING 2 — PREDICTOR DIVERGENCE (a limitation that bounds every coverage claim)

Per-locus coverage differs substantially by predictor:

| Locus | IEDB-EL | MHCflurry / MHCnuggets | Spread |
|---|---|---|---|
| A | 60.3% | 76.0% | **15.68 pp** |
| B | 89.4% | 99.0% | 9.59 pp |
| C | 70.2% | 98.2% | **27.99 pp** |
| DRB1 | 99.5% | 99.5% | 0.00 pp |
| DQB1 | 96.00% | 95.25% | 0.75 pp |
| DPB1 | 86.97% | 84.08% | 2.89 pp |

The class II locus values differ by at most 2.89 percentage points; class I differs by up to 27.99 points at locus C. The combined model value rounds to 100.0% in both universes and hides these differences. Report per-locus values and both predictor names.

---

## 4. FINDING 3 — BINDING AND BREADTH

A pass rate alone is weak evidence here: with a large allele panel, every saved MHC-I candidate meets the best-allele rank gate somewhere. Binding breadth and locus-specific coverage are more informative.

| Set | Median alleles bound | Range |
|---|---|---|
| Pool MHC-I (n=37) | 7 / 74 | 1–40 |
| Pool MHC-II (n=29) | 6 / 33 | 1–29 |
| Construct MHC-I (n=10) | 6 / 74 | 2–16 |
| Construct MHC-II (n=11) | 6 / 33 | 1–29 |

By pathogen: HIV MHC-I 6, Mpox MHC-I 6; HIV MHC-II 3, **Mpox MHC-II 7**.

**Agreement with Phase I:** mean Δ in percentile rank is −0.058 (class I, n=37, zero |Δ|>1) and
−0.025 (class II, n=29, three |Δ|>1). Medians exactly 0.000 in these saved outputs; the class-II methods are on different scales, so this is not proof of method equivalence. **Zero
gate-status flips** across 66 pool epitopes. The class II scatter is a scale difference: Phase I's
`recommended` resolves to EL for class I but BA/consensus for class II (correction B9).

**IC₅₀ availability:** 37/37 (class I) and 26/29 (class II) at the best consensus-supported allele.
The three unavailable values are recorded as `N/A`; they are not interpreted as nonbinding.

---

## 5. FINDING 4 — TCR RECOGNITION (MHC-I only, n = 10)

Three scoring approaches were used; they are **not independent experimental corroboration**. Each correlation uses only 10 selected MHC-I peptides, so report n and p as well as r.

| Tool | Direction check vs Phase4B EL rank | Interpretation |
|---|---|---|
| BigMHC-IM (primary) | r = −0.298 | expected sign; p=0.4024, not distinguishable from zero in this sample |
| PRIME 2.0 %Rank | r = +0.916 | expected sign; p=0.0002, but both are prediction scores |
| Calis et al. (2013) | r = +0.097 | exploratory; p=0.7899, training-set self-AUC 0.62 |

**Calis is exploratory only.** Its near-zero correlation with EL rank is not a validation of immunogenicity discrimination; this set has no measured immunogenicity outcomes. The local implementation reconstructed an amino-acid enrichment table from the published training set, so its training-set self-AUC is not an independent performance estimate.

**BigMHC reproducibility:** the saved Step C direction-check table reports r=1.0 for the same peptide–allele pairs relative to Phase 1Dd, with maximum absolute difference 3.97×10⁻⁷. This is an implementation check, not evidence of biological accuracy. Scores obtained from different best-allele searches cannot be interpreted as like-for-like biological comparisons.

PRIME and BigMHC select **different best alleles** for most peptides (3 agreements of 10). Never
compare their "best" columns as if allele-matched.

**Verdict split (Primary gate vs strict secondary stratum):**

| Class | PASS / STRICT_FAIL | PASS / STRICT_PASS |
|---|---|---|
| MHC-I | **10** | 0 |
| MHC-II | 1 | 10 |
| B-cell | 6 | 4 |

All 10 MHC-I epitopes fail the strict stratum because none reaches the proposed BigMHC-IM 0.70 cut; **four also fail the supplementary IC₅₀ ≤500 nM condition** (6/10 meet that IC₅₀ condition). Three reach the project BigMHC-IM score cut of 0.5. These are secondary labels on a frozen construct, not a retrospective candidate-selection step.

**Cross-reactivity (all 31 epitopes, both layers):** exact 8-mer matches against the reviewed human Swiss-Prot reference =
**0**. Fresh screen agrees with Phase 1Ec on **31/31**. Report both layers — the ≥70% identity
BLASTP rule alone is blind at 9–16 aa (0 of 20 verbatim human self-fragments caught).

---

## 6. FINDING 5 — PER-PATHOGEN COVERAGE (Section IV's immunodominance proxy)

| Subset | IEDB-EL (Phase4B) | MHCflurry (Phase1F) |
|---|---|---|
| HIV MHC-I (5) | 97.52% | 99.96% |
| **Mpox MHC-I (5)** | **83.87%** | 94.95% |
| HIV MHC-II (5) | 99.14% | >99.99% |
| Mpox MHC-II (6) | >99.99% | 92.79% |
| HIV combined (10) | 99.98% | >99.99% |
| Mpox combined (11) | >99.99% | 99.64% |

Both pathogen subsets have high combined model coverage, but that does not test immune competition. The informative class-I difference is
**Mpox MHC-I at 83.87% vs HIV MHC-I at 97.52%**, a ~14 pp deficit.

**The `NKRKRVIGL` thread.** Mpox A35R's MHC-I epitope is weakest on every independent measure:

| Measure | Value | Construct context |
|---|---|---|
| Binding breadth | 2 / 74 alleles | narrowest |
| BigMHC-IM | 0.0344 | lowest |
| PRIME %Rank | 4.386 | worst (next: 0.47) |
| NetMHC EL rank | 0.85 | worst |
| Contribution to construct coverage union | **0.0000 pp** | nil |
| Contribution to Mpox MHC-I coverage | 0.089 pp | negligible |

The saved Phase IV tables consistently rank `NKRKRVIGL` poorly among selected class-I peptides. Its zero marginal contribution means removal would not change the recorded allele-union coverage at reported precision; it does not imply zero chance of presentation or response.

---

## 7. FINDING 6 — PROTEASOMAL PROCESSING (new; never done in Phases I–II)

The whole 570-aa construct was submitted to IEDB's MHC-I processing predictor (proteasome + TAP + MHC) for 8–11-mers: **296/296 class-I tables**. Separate NetMHCIIpan-EL junction-binding requests returned **297/297 class-II calls**. Neither complete scan had failed calls. The class-II calls are binding predictions, not proteasomal processing results.

**Predicted liberation: 10/10 selected MHC-I epitopes meet the pipeline's C-terminal score criterion.** This is a model classification, not direct observation of proteasomal products.

**Internal-cleavage watchlist — 5 epitopes** carry an internal site scoring above their own
C-terminus, i.e. a possible destruction risk: `GPKEPFRDY`, `HHFNCRGEF`, `HWTTYMDTF`,
`RFALNPGLL`, `VYSTCTVPTM`. **Four of these also appear in the AAY-seam flag list** — they both
risk internal destruction and generate seam binders.

Sanity check: the proteasome score is a pure C-terminal cleavage score with zero deviation
across all alleles and lengths. Endpoint IC₅₀ and Step B EL ranks have Spearman ρ=+0.468 (p≈4.5×10⁻⁵ across 70 peptide–allele pairs from only 10 peptides); those pairs are **not independent biological observations**, so the p-value should not be used as a confirmatory test.

---

## 8. FINDING 7 — EPITOPE–LINKER SEAM ANALYSIS

### 8.1 Definition (read this before using any number below)

A seam window counts as **contested** at allele *a* only if it is (i) deep — ≥3 residues either
side of the seam, (ii) liberated, (iii) beats the best real epitope at *a* by **≥ m-fold**, and
(iv) for class II, its predicted binding core spans the junction.

**An earlier definition without the margin term produced badly inflated figures** (class II
"22.06% uncontested", a "77.93 pp gap"). Those are artifacts of comparing a maximum over
thousands of seam windows against a maximum over ~11 reference epitopes — near-certain by chance.
**They are superseded and must not be cited.**

Report as a margin curve, never a single number:

| Margin | Contested class I pairs | Uncontested MHC-I | Uncontested whole construct |
|---|---|---|---|
| 1× | 50 | 59.89% | 69.36% |
| 2× | 40 | 68.30% | 83.91% |
| **5× (primary)** | **28** | **88.70%** | **99.88%** |
| 10× | 24 | 89.11% | 99.98% |

At the 5× primary margin the whole-construct gap is **0.12 pp**.

### 8.2 Class II: GPGPG seams have no 5× contests

- **0 of 193** contested class II rows involve **GPGPG** — the MHC-II linker. The 193 split
  108 KK / 83 AAY / 2 adjuvant, i.e. they arise at the *other* cassettes' seams.
- Against an EL-rank random-peptide reference (n = 82,203 pairs): seam windows reach rank ≤10 at **7.17%**
  versus 10% expected, and rank ≤1 at 0.76%. The observed proportion is below the EL-rank percentile reference; the pairs are correlated and no inferential depletion test was performed.
- Residual class II 5× contests remain at KK, AAY and adjuvant seams. A separate null was not run for the full 5× contested criterion, so the residual cannot be declared an artifact. For the 8 DPB1 alleles, median seam-versus-real fold is 0.55–1.75× with 26–134 candidate windows versus 2–4 real epitopes per allele. Do not interpret the DPB1 0.0% derived coverage point as a measured presentation failure. Locus B 32.6% is a **class I** sensitivity result.

### 8.3 Class I: descriptive predicted enrichment with a plausible mechanism

The same IC₅₀ endpoint was applied to seam windows and unselected windows from the 45-aa adjuvant domain. This single-domain reference is not matched for composition or length:

| Measure | Value |
|---|---|
| Deep seam window–allele pairs with IC₅₀ ≤500 nM | 293/38,924 = 0.75% |
| Adjuvant-domain reference window–allele pairs with IC₅₀ ≤500 nM | 45/10,804 = 0.42% |
| **Enrichment** | **1.81×** |

Of 31 contested class I windows: **29 involve AAY**, and **21 (68%) end in AAY**. All 8 strongest
end in AAY.

**Proposed mechanism.** AAY terminates in tyrosine, and the highest-ranked AAY-ending windows also score well for processing and MHC-I binding. An aromatic C-terminal anchor is consistent with the flagged alleles; these predictions do not experimentally establish cleavage or presentation. Examples of flagged alleles with compatible C-terminal preferences are `B*35:01/05/17`, `A*26:01`, `C*14:02`.

Examples: `TYMDTFAAY` @ `C*14:02` (4.7 nM), `CTVPTMAAY` @ `A*26:01` (6.2 nM),
`EPFRDYAAY` @ `B*35:01` (6.7 nM).

**Per-locus class I at 5×:** A 57.8%, B 32.6%, C 60.3%.
**Per-pathogen uncontested MHC-I at 5×:** HIV 80.9%, **Mpox 62.0%**.

### 8.4 Required caveat

**Contested ≠ lost.** A 5-fold predicted affinity or rank advantage flags a possible competitor in this model. It does **not** show that a genuine epitope fails to be processed, presented or recognized. The method does not measure immunodominance or quantify prediction uncertainty.

The recorded 5× class-I contested alleles exclude the three highest-frequency class-I entries in this panel; this does not validate the contest model.

---

## 9. FINDING 8 — B-CELL EPITOPES AND CONSERVANCY

### 9.1 B-cell (n = 10)

| BepiPred tier | n |
|---|---|
| High | 1 |
| Medium | 6 |
| Deprioritized | 3 |

**SEMA-3D: all 10 return `NO`.** This means *screened without corroboration at the chosen threshold* — fewer than 50% of
residues fall in a SEMA-3D patch. It does **not** mean structurally refuted. Zero are
`UNSCREENED` (correction B13).

> Deviation #20's tally ("1 corroborated, 6 not, 2 unscreenable") describes the **superseded**
> construct `Vax_Final_6f34b53e`. The final construct is 0 / 10 / 0. The register must be updated
> before it is cited (correction B12).

**Surface exposure comes from Phase 1De per-antigen native folds, not the construct model.** The
construct model is unusable for this: pTM 0.17, and **none of its 525 non-adjuvant residues
reaches pLDDT 70**. Of the 10 epitopes, **5 sit on folds that are not globally determined** — their
exposure is indicative only. Mean RSA: HIV 0.37, Mpox 0.34. Correlation with BepiPred score is
weak (ρ = 0.112).

### 9.2 Conservancy — gate reproduces at reported precision

| Antigen | Mean conservancy |
|---|---|
| HIV gp120 | 4.54% |
| HIV gp41 | 5.69% |
| HIV p17 | 4.67% |
| HIV p24 | 9.97% |
| Mpox A35R | 37.80% |
| `Mpox_B5R` (pipeline label) | 45.02% |
| Mpox L1R | 60.61% |

All 7 published means reproduce **at the reported precision**, and all 31 construct epitopes reproduce their stored
conservancy and hit ratio from the Phase 1C variant pools.

> **The gate must run on `Raw_Full`, not `Min_50pct`.** The published means describe the full
> 30,268-candidate pool; `Min_50pct` holds only the 1,238 survivors and cannot reproduce them
> (correction B6).

Construct-epitope conservancy: HIV mean **75.5%** (n=14), Mpox mean **99.7%** (n=17). These values are conditioned on selection from the conserved tail and cannot be compared with unselected full-pool means as an unbiased pathogen difference.
Conservancy is measured within-clade against 21–30 isolates per antigen, **not** "global strains"
(correction B5).

---

## 10. HOW TO CITE THESE NUMBERS

**Traceability.** `Step_Outputs/Phase4/StepG/Phase4G_NumbersLedger_20260921_2022.csv` has 155 entries with source-file fields. This does not prove that every prose number appears as a separate ledger row; use the named Step A–F CSVs for independent checks. The per-epitope dossier is `Step_Outputs/Phase4/StepG/Phase4G_Immunogenicity_Dossier_20260921_2022.csv` — 31
rows, zero blanks, zero `UNRESOLVED`, with `Primary_Gate` and `Secondary_Stratum` as separate
columns so both readings stay recoverable.

**Reporting rules this project adopted:**
1. State **n** with every proportion and correlation.
2. Report **breadth and locus-specific coverage** alongside any gate status for a large allele panel.
3. Keep `N/A` (tool cannot score) distinct from `UNRESOLVED` (no answer yet). Never count
   `UNRESOLVED` as a pass.
4. Report per-locus **ranges** across predictors for class I coverage.
5. For any seam/junction claim, give the **margin and the null** — a bare "some window beats X"
   is an extreme-value artifact.

**Superseded figures — do not cite:** class II uncontested 22.06%; 77.93 pp class II gap;
18.67 pp combined gap; a DPB1 0.0% *presentation-failure interpretation*; the 111-flag and 126-flag pre-correction counts; Deviation #20's
SEMA tally.

---

## 11. WHAT PHASE IV DOES NOT ESTABLISH

- **No immune dynamics.** No antibody titres, cytokine levels, cell counts, or time-courses.
  Nothing here substitutes for C-ImmSim's temporal output.
- **No reliable whole-construct fold.** Whole-model pTM is 0.17; none of the 525 non-adjuvant residues reaches pLDDT 70.
- **No Phase III integration in this branch.** Docking, MD and C-ImmSim data were not used in Phase IV. The latest manuscript PDF reports such results, but its raw run provenance is not established by the Phase IV files. Validate it separately before a combined analysis (see `Handoff.md`).
- **No experimental confirmation.** Every figure is a prediction. Class I coverage additionally
  varies by up to ~28 pp with predictor choice.
- **Coverage is a proxy population.** AFND Singapore Riau Malay, Southeast Asian (Austronesian) —
  not a Philippine dataset.

---

## 12. THE THREE FINDINGS THAT SHOULD REACH THE PAPER

Ordered by strength of evidence:

1. **`A*11:01` (0.177) and `A*33:03` (0.109) are bound by no construct MHC-I epitope under either
   predictor**, and by no candidate in the 37-epitope pool under Phase IV IEDB-EL. No null model is needed to describe these saved binding calls. This is
   the most actionable finding if Phase I is ever reopened.
2. **Class I coverage is predictor-dependent by up to 27.99 pp** (locus C). This bounds the
   confidence of every class I claim in the paper.
3. **A predicted class I seam pattern** — 1.81× descriptive rate ratio against the 45-aa adjuvant-domain reference, with 29 of 31 qualifying windows involving AAY. This is a hypothesis for Discussion, not a measured immunodominance result.
   Belongs in Discussion/limitations, not as a redesign trigger. The GPGPG class II seams have zero 5× contested rows; other class II seams retain predicted contests and need cautious wording.

---

## 13. PROVENANCE

Phase I and Phase II remain byte-identical to the saved SHA-256 baseline across **1,810 files** as independently checked on 2026-09-24
(`Step_Outputs/Phase4/_rule1_baseline_sha256.txt`).

Three escalations were raised and resolved during execution:

| # | Trigger | Resolution |
|---|---|---|
| 1 | Step 2 — IC₅₀ misparse | Unlabeled column read as nM; caught by sign/range checks, corrected |
| 2 | Step 4 — seam peptides out-scoring epitopes | Reweighted by allele frequency in Step 5 |
| 3 | Step 5 — uncontested-coverage gap | Metric re-specified with margin + null; artifact removed |

The earlier parsing and metric incidents produced plausible-looking wrong intermediate numbers. The corrected outputs above supersede those intermediate figures. The handoff records the practical checks for future work.

---

## Appendix A. Selected epitope results from the final 31-row dossier

The following tables are a compact view of the saved dossier, not new predictions. “Breadth” counts alleles at the class-specific primary rank threshold. `N/A` is intentional for class-inapplicable metrics. The full 153-column dossier retains allele lists, positions and both verdict readings.

### A1. MHC-I selected peptides (n=10)

| Peptide | Origin | Target label | EL best rank (%) | Breadth /74 | BigMHC-IM best | PRIME best %Rank | Processing rule | Internal-cleavage watch | Conservancy (%) |
|---|---|---|---:|---:|---:|---:|---|---|---:|
| `APGSPTNLEF` | Mpox | Mpox_L1R | 0.11 | 14 | 0.0809 | 0.177 | LIBERATED | NO | 100.0 |
| `GPKEPFRDY` | HIV | HIV_p24 | 0.15 | 16 | 0.4149 | 0.193 | LIBERATED | YES | 96.67 |
| `HHFNCRGEF` | HIV | HIV_gp120 | 0.11 | 5 | 0.3846 | 0.47 | LIBERATED | YES | 83.33 |
| `HWTTYMDTF` | Mpox | Mpox_L1R | 0.12 | 6 | 0.54 | 0.018 | LIBERATED | YES | 100.0 |
| `NKRKRVIGL` | Mpox | Mpox_A35R | 0.85 | 2 | 0.0344 | 4.386 | LIBERATED | NO | 100.0 |
| `QTSVFSATVY` | Mpox | Mpox_A35R | 0.19 | 3 | 0.3187 | 0.112 | LIBERATED | NO | 100.0 |
| `RFALNPGLL` | HIV | HIV_p17 | 0.25 | 6 | 0.3604 | 0.282 | LIBERATED | YES | 53.33 |
| `SEGATPQDL` | HIV | HIV_p24 | 0.13 | 5 | 0.1838 | 0.216 | LIBERATED | NO | 86.67 |
| `TMGAASITL` | HIV | HIV_gp41 | 0.39 | 7 | 0.5664 | 0.278 | LIBERATED | NO | 80.0 |
| `VYSTCTVPTM` | Mpox | Mpox_B5R | 0.3 | 6 | 0.6178 | 0.284 | LIBERATED | YES | 100.0 |

### A2. MHC-II selected peptides (n=11)

| Peptide | Origin | Target label | EL best rank (%) | Breadth /33 | IC₅₀ at best consensus allele (nM) | Conservancy (%) |
|---|---|---|---:|---:|---:|---:|
| `EQEIESLEATYHIII` | Mpox | Mpox_B5R | 1.9 | 4 | 101.0 | 100.0 |
| `GIVQQQSNLLRAIEA` | HIV | HIV_gp41 | 1.3 | 3 | 31.5 | 83.33 |
| `GNPITKTTSDYQDSD` | Mpox | Mpox_A35R | 0.74 | 10 | 796.6 | 95.65 |
| `HNVWATHACVPTDPN` | HIV | HIV_gp120 | 7.9 | 2 | 204.8 | 63.33 |
| `IFGFLGAAGSTMGAA` | HIV | HIV_gp41 | 2.5 | 3 | 3.3 | 86.67 |
| `NDKIKLILANKENVH` | Mpox | Mpox_L1R | 0.15 | 29 | 34.9 | 100.0 |
| `QCTHGIKPVVSTQLL` | HIV | HIV_gp120 | 0.72 | 9 | 4.8 | 70.0 |
| `SNGLISGSTFSIGGV` | Mpox | Mpox_B5R | 2.7 | 6 | 56.2 | 100.0 |
| `SSTTQYDHKESCNGL` | Mpox | Mpox_A35R | 6.8 | 1 | 3617.6 | 100.0 |
| `TPEQKAYVPAMFTAA` | Mpox | Mpox_L1R | 0.59 | 8 | 201.6 | 100.0 |
| `YKRWIILGLNKIVRM` | HIV | HIV_p24 | 4.5 | 6 | 3.9 | 86.67 |

### A3. B-cell selected peptides (n=10)

| Peptide | Origin | Target label | BepiPred mean | Phase I tier | SEMA overlap (%) | Native fold | Native RSA mean | Conservancy (%) |
|---|---|---|---:|---|---:|---|---:|---:|
| `AATETYSGLTPEQKAY` | Mpox | Mpox_L1R | 0.47 | Deprioritized | 6.2 | DETERMINED | 0.306 | 100.0 |
| `ANASAQTKCDIEIGNF` | Mpox | Mpox_L1R | 0.506 | Medium | 0.0 | DETERMINED | 0.348 | 100.0 |
| `CQPLQLEHGSCQPVKE` | Mpox | Mpox_B5R | 0.549 | Medium | 0.0 | UNDETERMINED | 0.368 | 100.0 |
| `DHKESCNGLYYQGSCY` | Mpox | Mpox_A35R | 0.507 | Deprioritized | 0.0 | UNDETERMINED | 0.403 | 100.0 |
| `EDTWGSDGNPITKTTS` | Mpox | Mpox_A35R | 0.607 | High | 6.2 | UNDETERMINED | 0.258 | 100.0 |
| `EHGSCQPVKEKYSFGE` | Mpox | Mpox_B5R | 0.542 | Medium | 0.0 | UNDETERMINED | 0.349 | 100.0 |
| `EMMTACQGVGGPSHKA` | HIV | HIV_p24 | 0.521 | Medium | 12.5 | DETERMINED | 0.354 | 83.33 |
| `GGKLDAWEKIRLRPGG` | HIV | HIV_p17 | 0.568 | Medium | 12.5 | DETERMINED | 0.366 | 53.33 |
| `VWGIKQLQARVLAVER` | HIV | HIV_gp41 | 0.519 | Medium | 0.0 | UNDETERMINED | 0.46 | 76.67 |
| `WDQSLKPCVKLTPLCV` | HIV | HIV_gp120 | 0.488 | Deprioritized | 0.0 | DETERMINED | 0.281 | 53.33 |

**Scope of these tables:** EL rank, BigMHC-IM, PRIME, processing and native-fold values are predictions. The IC₅₀ column comes from a separate consensus-supported allele and must not be read as the affinity for the EL-best allele. SEMA overlap is corroboration screening, not experimental refutation. `Mpox_B5R` is the frozen pipeline label; the manuscript uses B6R for the intended target, pending a source-protein naming crosswalk.

*End of consolidated Phase IV results.*
