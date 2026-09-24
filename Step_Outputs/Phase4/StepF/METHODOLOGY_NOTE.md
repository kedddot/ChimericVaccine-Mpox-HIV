# Phase 4F -- B-cell + Conservancy: methodology notes

Generated: 20260921_1949

## STATUS: complete; all regression gates passed; no escalation

Mostly reuse -- **no new predictions**. Every stored value that can be recomputed deterministically was recomputed as a regression check (BepiPred tier, SEMA overlap, per-target conservancy means, per-epitope conservancy). All matched.

## Conservancy

**Hard gate PASSED -- all 7 published per-target means reproduce exactly** (mean over all 30,268 Phase 1Dc candidates):

| target | n candidates | published | reproduced | max | % at 100 |
|---|---|---|---|---|---|
| HIV_gp120 | 13251 | 4.54% | 4.54% | 90.0% | 0.0% |
| HIV_gp41 | 7783 | 5.69% | 5.69% | 96.67% | 0.0% |
| HIV_p17 | 3982 | 4.67% | 4.67% | 76.67% | 0.0% |
| HIV_p24 | 3141 | 9.97% | 9.97% | 96.67% | 0.0% |
| Mpox_A35R | 657 | 37.8% | 37.8% | 100.0% | 4.9% |
| Mpox_B5R | 931 | 45.02% | 45.02% | 100.0% | 10.4% |
| Mpox_L1R | 523 | 60.61% | 60.61% | 100.0% | 22.4% |

HIV per-target means span 4.54-9.97%; Mpox 37.8-60.61% (brief: 4.54-9.97 / 37.80-60.61). **Source note:** the brief says to use `Min_50pct`, but that file holds only the >=50% survivors (1,238 rows) and cannot reproduce a mean over all candidates; the gate uses `Raw_Conservancy/Phase1Dc_Raw_Full`, the full pool the published table describes. `Min_50pct` is what the 31 construct epitopes were selected from.

Independent recompute (Phase 1C variant pools, exact substring, Phase1Dc's own method): **all 31 construct epitopes reproduce both the stored conservancy and hit ratio**.

| pathogen | n | mean | min | max | variant pool sizes | mean over ALL candidates of those targets |
|---|---|---|---|---|---|---|
| HIV | 14 | 75.5% | 53.3% | 96.7% | [30] | 5.5% |
| Mpox | 17 | 99.7% | 95.7% | 100.0% | [21, 23, 30] | 46.6% |

The construct's HIV epitopes (53.3-96.7%) sit far above the HIV candidate means (4.5-10%): selection took the conserved tail of a very diverse pool, so the pool means describe the diversity of the antigens, not the construct. Variant pools are small (n = 21-30), so values are coarse (one variant = 3-5 pp); HIV epitopes are conserved across fewer of a smaller number of CRF01_AE isolates -- interpret as within-clade conservancy only.

## B-cell -- BepiPred (Phase I tiers, disagreement #3)

Tiers recomputed with the imported `classify_bcell_tier` for all 39 pool candidates: **39/39 match Phase 1Db**. `Primary_Gate` = Phase I tier; `Secondary_Stratum` = the methodology's stricter 'High' (mean >= 0.50 and >= 75% of residues above).

Construct (n = 10): Phase I tiers {'Deprioritized': 3, 'Medium': 6, 'High': 1}; methodology-High 4/10. Pool (n = 39): {'Deprioritized': 14, 'Medium': 24, 'High': 1}. Three construct B-cell epitopes are `Deprioritized` -- Phase I retained them (tiers deprioritise, they never exclude). HIV vs Mpox:

| pathogen | n | Phase I tiers | methodology-High | mean BepiPred |
|---|---|---|---|---|
| HIV | 4 | {'Medium': 3, 'Deprioritized': 1} | 1/4 | 0.524 |
| Mpox | 6 | {'Deprioritized': 2, 'Medium': 3, 'High': 1} | 3/6 | 0.530 |

## B-cell -- SEMA-3D (Phase 1De, reused; not re-run because `foldseek` is a Linux binary)

All 10 construct B-cell epitopes are `SEMA_Corroborated = NO` (overlap 0-12.5% of 16 residues scored). **`NO` means screened and not corroborated: fewer than 50% of the residues fall in a SEMA-3D conformational patch (score >= ln 4). It is a lower-priority *linear-only* call, not a refutation** -- SEMA-3D's own calibration is TPR 0.66 at FPR 0.07, so a miss is far from disproof. **`UNSCREENED` (no assessable fold) is a different state and is never a negative**; none of the construct epitopes are UNSCREENED (pool of 39: {'NO': 39}). The overlap was recomputed from Phase 1De's per-residue scores for every locatable epitope: all match.

## B-cell -- native-fold overlay (why not the construct model)

The construct's own AlphaFold fold is undetermined outside the adjuvant: **pTM 0.17**; the adjuvant domain (first 45 residues) has mean pLDDT 87.2 (98% of residues >= 70), but the non-adjuvant 525 residues have mean pLDDT 27.2 with **0.0% >= 70** (verified from the model file). Surface accessibility computed on a model that is not a determined fold would be meaningless, so it is computed on the **Phase 1De per-antigen native folds** instead (each antigen folded alone). Each epitope is located in its antigen's fold by exact sequence match; per-residue relative SASA (Shrake-Rupley, probe 1.4 A, Tien/Wilke 2013 maxima shipped with Biopython) and pLDDT are averaged over the 16 residues. Sanity checks passed before use: PDB B-factor pLDDT reproduces Phase 1De's stored Mean_pLDDT for all 7 antigens, and RSA direction is correct (polar residues more exposed than buried-prone hydrophobics).

| peptide | target | fold verdict (pTM) | mean RSA | residues exposed (RSA>=0.25) | native pLDDT | construct-model pLDDT |
|---|---|---|---|---|---|---|
| AATETYSGLTPEQKAY | Mpox_L1R | DETERMINED (0.68) | 0.31 | 0.56 | 92.5 | 49.8 |
| ANASAQTKCDIEIGNF | Mpox_L1R | DETERMINED (0.68) | 0.35 | 0.56 | 90.6 | 37.1 |
| CQPLQLEHGSCQPVKE | Mpox_B5R | UNDETERMINED (0.43) | 0.37 | 0.69 | 87.4 | 23.3 |
| DHKESCNGLYYQGSCY | Mpox_A35R | UNDETERMINED (0.47) | 0.40 | 0.69 | 80.5 | 20.8 |
| EDTWGSDGNPITKTTS | Mpox_A35R | UNDETERMINED (0.47) | 0.26 | 0.56 | 86.5 | 39.3 |
| EHGSCQPVKEKYSFGE | Mpox_B5R | UNDETERMINED (0.43) | 0.35 | 0.62 | 86.1 | 22.5 |
| EMMTACQGVGGPSHKA | HIV_p24 | DETERMINED (0.69) | 0.35 | 0.69 | 72.8 | 21.9 |
| GGKLDAWEKIRLRPGG | HIV_p17 | DETERMINED (0.74) | 0.37 | 0.69 | 89.5 | 40.1 |
| VWGIKQLQARVLAVER | HIV_gp41 | UNDETERMINED (0.39) | 0.46 | 0.94 | 70.7 | 41.1 |
| WDQSLKPCVKLTPLCV | HIV_gp120 | DETERMINED (0.89) | 0.28 | 0.50 | 85.4 | 31.8 |

Native folds: **5 of 10** construct B-cell epitopes sit on a DETERMINED global fold (gp120, p17, p24, L1R); **5 of 10 sit on an UNDETERMINED fold** (gp41 pTM 0.39, A35R 0.47, B5R 0.43) and their exposure is INDICATIVE ONLY. Caveats for all: monomer folds (oligomer interfaces not modelled), and HIV Env glycans are not modelled, so gp120/gp41 exposure is over-estimated. Exposure is context, not a gate.

Direction cross-check: BepiPred mean vs native-fold RSA across the 39 unique pool B-cell candidates, Spearman rho = +0.112 (p = 0.5, n = 39). The sign is as expected but the correlation is NOT distinguishable from zero at this n: it passes the sign check and corroborates nothing.

## HIV vs Mpox (every metric split)

| metric | HIV | Mpox |
|---|---|---|
| construct B-cell epitopes (n) | 4 | 6 |
| mean BepiPred | 0.524 | 0.530 |
| SEMA corroborated (YES) | 0/4 | 0/6 |
| mean native-fold RSA (located epitopes) | 0.37 | 0.34 |
| mean conservancy, construct B-cell epitopes | 66.7% | 100.0% |
| mean conservancy, all 31 construct epitopes | 75.5% | 99.7% |

## Caveats

- BepiPred and SEMA-3D are predictors; none of this is experimental evidence of antibody binding.
- Conservancy is exact-match within the Phase 1C variant pools (n = 21-30 per target), not IEDB's alignment tool (Phase 1Dc note) and not a global diversity measure.
- Native-fold exposure on undetermined folds is indicative only; glycan shielding and oligomerisation are not modelled.
