# Phase 4E -- Population Coverage + Frequency-Weighted Junction Analysis: methodology notes

Generated: 20260921_1949

## STATUS: complete; no open escalation

Opus's two rulings are applied: (1) the construct stays frozen; (2) **the original 'contested' definition was a specification flaw** (a maximum over 2491 seam windows vs a maximum over 11 reference epitopes finds contest by construction). Its outputs were discarded and are not reported. Everything below uses the corrected definition.

## Method

All coverage arithmetic uses Phase 1F's imported `compute_overall_coverage` / `compute_cumulative_coverage` (per-locus s = sum of allele frequencies capped at 1.0, diploid two-chromosome term, complement product across loci). The methodology's own formula omits the diploid term (disagreement #5) and is implemented nowhere in this code. Loci grouped separately: {'A': 22, 'B': 33, 'C': 19, 'DRB1': 18, 'DQB1': 7, 'DPB1': 8}. **DQA1 and DPA1 do not exist in `ALLELE_FREQ`** (beta chains only), so class II coverage is on DRB1/DQB1/DPB1 -- a property of the Phase 1F table.

**Two allele universes, never mixed.** Phase 1F's stored coverage used `Binding_Alleles_Recomputed` (MHCflurry/MHCnuggets); the Step 4 junction scan used `Binding_Alleles_Primary` (IEDB netmhcpan_el / netmhciipan_el at Phase I rank gates). The level-1 regression gate reproduces Phase 1F in ITS universe; everything else runs in the Phase 4B universe.

## Level 1 -- per-epitope (hard regression gate)

**31/31 reproduce Phase 1F `Overall_Coverage_Pct`** (21 MHC-restricted numeric values match the stored 2 dp; 10 B-cell epitopes are `NOT_APPLICABLE`, never 0.0). Trigger 4 did not fire. Per-epitope detail: `Phase4E_EpitopeCoverage_20260921_1949.csv`.

## 1. COVERAGE GAPS (report ahead of the junction analysis)

Alleles that **no real construct epitope binds** cannot be covered, whatever the seams do. This needs no null model. Class I, by frequency:

| allele | freq | epitopes (IEDB-EL) | epitopes (MHCflurry) | status |
|---|---|---|---|---|
| HLA-C*08:01 | 0.187 | 0 | 1 | GAP_IEDB_EL_ONLY |
| HLA-A*11:01 | 0.177 | 0 | 0 | GAP_BOTH_PREDICTORS |
| HLA-A*24:07 | 0.157 | 3 | 5 | COVERED_BOTH |
| HLA-C*04:01 | 0.122 | 1 | 5 | COVERED_BOTH |
| HLA-C*07:03 | 0.112 | 1 | 5 | COVERED_BOTH |
| HLA-A*33:03 | 0.109 | 0 | 0 | GAP_BOTH_PREDICTORS |
| HLA-B*18:01 | 0.099 | 0 | 2 | GAP_IEDB_EL_ONLY |
| HLA-C*04:03 | 0.098 | 0 | 1 | GAP_IEDB_EL_ONLY |
| HLA-B*15:02 | 0.084 | 3 | 4 | COVERED_BOTH |
| HLA-C*07:04 | 0.084 | 0 | 1 | GAP_IEDB_EL_ONLY |

**`HLA-A*11:01` (f = 0.177) and `HLA-A*33:03` (f = 0.109) are gaps under BOTH predictors** -- an unambiguous coverage gap, and the reason locus A sits at 60.3% (IEDB-EL) / 76.0% (MHCflurry). **Correction to the ruling:** `HLA-C*08:01` (f = 0.187) is a gap under IEDB-EL only (MHCflurry has `TMGAASITL` binding it), as are `B*18:01`, `C*04:03` and `C*07:04`; so the B and C gaps are predictor-dependent, not unambiguous.

Frequency mass with no binding epitope, per locus (of the locus total):

| locus | IEDB-EL | MHCflurry/MHCnuggets | locus total |
|---|---|---|---|
| A | 0.446 | 0.306 | 0.816 |
| B | 0.224 | 0.000 | 0.898 |
| C | 0.491 | 0.080 | 0.945 |
| DRB1 | 0.000 | 0.000 | 0.932 |
| DQB1 | 0.000 | 0.018 | 0.800 |
| DPB1 | 0.000 | 0.038 | 0.639 |

## 2. METHOD DIVERGENCE (a stated limitation -- report the range, never the favourable end)

| locus | IEDB-EL | MHCflurry/MHCnuggets | range | spread |
|---|---|---|---|---|
| A | 60.31% | 75.99% | 60.31-75.99% | 15.68 pp |
| B | 89.37% | 98.96% | 89.37-98.96% | 9.59 pp |
| C | 70.19% | 98.18% | 70.19-98.18% | 27.99 pp |
| DRB1 | 99.54% | 99.54% | 99.54-99.54% | 0.00 pp |
| DQB1 | 96.00% | 95.25% | 95.25-96.00% | 0.75 pp |
| DPB1 | 86.97% | 84.08% | 84.08-86.97% | 2.89 pp |

Locus C spans 70.2-98.2% (28 pp) from predictor choice alone; class II loci agree within 3 pp. Class I coverage therefore carries substantial method uncertainty.

## 3. Levels 2 and 3 -- construct and per pathogen (B-cell excluded, not counted as zero)

| scope | Phase 1F universe (n alleles) | Phase 4B universe (n alleles) |
|---|---|---|
| CONSTRUCT (21 MHC epitopes) | >99.99% (96) | >99.99% (76) |
| MHC-I only (10) | >99.99% (65) | 98.7425% (43) |
| MHC-II only (11) | >99.99% (31) | >99.99% (33) |
| HIV (10: 5 I + 5 II) | >99.99% (88) | 99.9786% (53) |
| Mpox (11: 5 I + 6 II) | 99.6361% (49) | >99.99% (55) |
| HIV MHC-I (5) | 99.9625% (57) | 97.5238% (37) |
| HIV MHC-II (5) | >99.99% (31) | 99.1362% (16) |
| Mpox MHC-I (5) | 94.9532% (35) | 83.8669% (23) |
| Mpox MHC-II (6) | 92.7900% (14) | >99.99% (32) |

Per-locus, whole construct -- Phase 4B universe: `A:60.31;B:89.37;C:70.19;DPB1:86.97;DQB1:96.00;DRB1:99.54`; Phase 1F universe: `A:75.99;B:98.96;C:98.18;DPB1:84.08;DQB1:95.25;DRB1:99.54`. The combined figure saturates (complement product across six loci), so per-locus and per-class numbers are the informative layer.

**Asymmetry:** Mpox MHC-I 83.9% vs HIV MHC-I 97.5% (IEDB-EL; MHCflurry: 95.0% vs 100.0%). Consistent with A35R being the weakest antigen throughout Steps 3-5. Not trigger 7: both pathogens exceed 99.9% combined; class II is >99% for both.

**Standing flag `NKRKRVIGL` (Mpox A35R):** binds 2/107 alleles; own coverage 4.72%; marginal drop when removed: 0.0000 pp (construct union), 0.0890 pp (within MHC-I). Consistent cross-phase signal, not a new anomaly.

## 4. Junction analysis -- corrected definition

A junction window is **contested at allele a** only if ALL hold: (1) deep (>=3 residues each side); (2) passes the liberation filter (C-terminal cleavage >= construct median; a proxy for class II); (3) beats the **best real epitope binding a** by >= m-fold (MHC-I: IC50 from the processing endpoint; MHC-II: EL-rank ratio, because EL has no IC50); (4) MHC-II only: the predicted core spans the junction. It follows that contest is an allele-level property. Primary margin m = 5; the whole curve is shown, with n.

(epitope,allele) pairs at risk: MHC-I n=70, MHC-II n=81; alleles with a real comparator: class I 43, class II 33.

| m | contested pairs I | contested pairs II | contested alleles I / II | seam windows I / rows II | pop. carrying >=1 contested allele (I / II) | uncontested: construct | MHC-I | MHC-II |
|---|---|---|---|---|---|---|---|---|
| 1x | 50/70 | 66/81 | 32 / 28 | 80 / 740 | 95.24% / 99.98% | 69.362% | 59.89% | 23.62% |
| 2x | 40/70 | 60/81 | 26 / 25 | 57 / 385 | 93.96% / 99.91% | 83.911% | 68.30% | 49.25% |
| 5x (primary) | 28/70 | 43/81 | 16 / 18 | 31 / 193 | 79.45% / 94.47% | 99.880% | 88.70% | 98.94% |
| 10x | 24/70 | 24/81 | 14 / 11 | 19 / 118 | 78.21% / 68.87% | 99.983% | 89.11% | 99.84% |

**Reading the curve.** At 1x the class II point is dominated by the extreme-value effect the null calibration exposes (below), so it is shown only to complete the curve, not as a finding. **The whole-construct 1x point (69.36%) is just under the ~70% floor; it was not treated as an escalation because the ruling defines contest at >=5x -- flagged here for Opus to confirm.** From 5x the class II contest largely disappears, while class I persists. **Caveat:** 18 of 33 class II alleles still have a qualifying seam window at 5x (KK and AAY seams; none at GPGPG). No random-window null was run at the 5x criterion specifically, so no claim is made about whether that residual exceeds chance; its coverage cost is small (see the primary table).

Sensitivity (class I additionally requires IC50 <= 500 nM, i.e. a genuine strong binder): uncontested MHC-I at 1/2/5/10x = 76.5% / 81.1% / 90.0% / 90.2%.

**Primary (m = 5x):**

| scope | raw | uncontested | gap pp | union alleles raw -> uncontested |
|---|---|---|---|---|
| CONSTRUCT (21 MHC epitopes) | >99.99% | 99.88% | 0.12 | 76 -> 42 |
| MHC-I only (10) | 98.7425% | 88.70% | 10.05 | 43 -> 27 |
| MHC-II only (11) | >99.99% | 98.94% | 1.06 | 33 -> 15 |
| HIV (10: 5 I + 5 II) | 99.9786% | 98.96% | 1.02 | 53 -> 31 |
| Mpox (11: 5 I + 6 II) | >99.99% | 99.43% | 0.57 | 55 -> 28 |
| HIV MHC-I (5) | 97.5238% | 80.91% | 16.62 | 37 -> 22 |
| Mpox MHC-I (5) | 83.8669% | 62.04% | 21.83 | 23 -> 14 |

Per-locus construct: raw `A:60.31;B:89.37;C:70.19;DPB1:86.97;DQB1:96.00;DRB1:99.54` -> uncontested `A:57.75;B:32.60;C:60.31;DQB1:91.00;DRB1:88.17`. **Cost of the junction problem in population terms: 0.12 pp for the whole construct; 10.0 pp within MHC-I; 1.1 pp within MHC-II.** Class II epitopes compensate for class I losses because they are presented by a disjoint set of alleles.

**Per-locus effect at 5x (the informative layer -- Phase 1F's own warning that the combined figure saturates):** A 60.3% -> 57.8%; B 89.4% -> 32.6%; C 70.2% -> 60.3%; DRB1 99.5% -> 88.2%; DQB1 96.0% -> 91.0%; DPB1 87.0% -> 0.0%. A locus with no uncontested allele reads 0.0% (here DPB1: every DPB1 allele a real epitope binds is contested at 5x, via KK/AAY seams). The whole-construct figure stays at 99.88% only because coverage multiplies across six loci; the per-locus losses at B and DPB1 are far larger than the combined or per-class numbers suggest. This does not change the thresholds in the ruling, which are defined on the whole construct, but it is the honest per-locus picture and is flagged for Opus.

### Per-epitope contested fraction at 5x (n = 21; ranked descending)

| peptide | class | pathogen | epitope cov | contested cov | fraction | contested / bound alleles |
|---|---|---|---|---|---|---|
| TPEQKAYVPAMFTAA | MHC-II | Mpox | 86.97% | 86.97% | 1.00 | 8/8 |
| NDKIKLILANKENVH | MHC-II | Mpox | 99.96% | 93.54% | 0.94 | 17/29 |
| YKRWIILGLNKIVRM | MHC-II | HIV | 41.07% | 37.45% | 0.91 | 5/6 |
| GPKEPFRDY | MHC-I | HIV | 70.13% | 59.06% | 0.84 | 10/16 |
| QTSVFSATVY | MHC-I | Mpox | 24.59% | 19.42% | 0.79 | 2/3 |
| APGSPTNLEF | MHC-I | Mpox | 65.08% | 47.73% | 0.73 | 7/14 |
| GNPITKTTSDYQDSD | MHC-II | Mpox | 83.33% | 54.19% | 0.65 | 2/10 |
| SEGATPQDL | MHC-I | HIV | 32.43% | 20.61% | 0.64 | 3/5 |
| SNGLISGSTFSIGGV | MHC-II | Mpox | 81.24% | 45.90% | 0.56 | 4/6 |
| IFGFLGAAGSTMGAA | MHC-II | HIV | 25.70% | 10.70% | 0.42 | 2/3 |
| QCTHGIKPVVSTQLL | MHC-II | HIV | 94.77% | 30.44% | 0.32 | 2/9 |
| TMGAASITL | MHC-I | HIV | 34.21% | 10.51% | 0.31 | 1/7 |
| NKRKRVIGL | MHC-I | Mpox | 4.72% | 1.00% | 0.21 | 1/2 |
| HWTTYMDTF | MHC-I | Mpox | 43.75% | 9.18% | 0.21 | 1/6 |
| RFALNPGLL | MHC-I | HIV | 43.75% | 9.18% | 0.21 | 1/6 |
| VYSTCTVPTM | MHC-I | Mpox | 43.75% | 9.18% | 0.21 | 1/6 |
| HHFNCRGEF | MHC-I | HIV | 45.27% | 7.84% | 0.17 | 1/5 |
| EQEIESLEATYHIII | MHC-II | Mpox | 60.31% | 3.57% | 0.06 | 2/4 |
| GIVQQQSNLLRAIEA | MHC-II | HIV | 32.25% | 1.79% | 0.06 | 1/3 |
| HNVWATHACVPTDPN | MHC-II | HIV | 60.69% | 0.00% | 0.00 | 0/2 |
| SSTTQYDHKESCNGL | MHC-II | Mpox | 8.99% | 0.00% | 0.00 | 0/1 |

Invariant checked for all 21: contested coverage <= epitope coverage; and uncontested <= raw for every scope and margin.

### Class I finding (real, reportable as a Discussion/limitations item)

At 5x, 31 class I seam windows are contested on 16 alleles; **21/31 (68%) end in `AAY`**, and the 8 strongest (IC50 4.7-13.3 nM) all do -- verified. They fall on B*35:01/05/17, A*26:01 and C*14:02, alleles whose real construct epitopes end in F/Y (data-side check: HLA-A*26:01:Y; HLA-B*35:01:FY; HLA-B*35:05:FY; HLA-B*35:17:FY; HLA-C*14:02:FLM); the aromatic-PΩ-pocket explanation is from prior knowledge of these alleles and is consistent with, but not proven by, this data. Of the remaining 10: 8 contain AAY without ending in it and 2 involve only EAAAK/GPGPG (so 29/31 windows involve AAY). Linker breakdown (class I strong-binder rate, exclusive linker sets, n = window-allele pairs): windows ending in AAY 6.12% (x14.7 the natural-window rate, n=2664); AAY set 2.87% (n=9324); GPGPG 0.06% (n=20350); KK 0.09% (n=7400); natural adjuvant windows 0.42% (n=10804). AAY's terminal Tyr is the mechanism; GPGPG and KK seams are depleted for strong class I binders.

Contested class I alleles are all outside the top three by frequency (none of C*08:01, A*11:01, A*24:07); the most frequent is B*15:02 (f = 0.084, rank 9 of 74). Contested alleles (frequency rank): HLA-B*15:02 (0.084, #9), HLA-B*35:05 (0.069, #13), HLA-B*15:13 (0.069, #12), HLA-B*44:03 (0.064, #15), HLA-B*13:01 (0.054, #18), HLA-C*14:02 (0.047, #22), HLA-B*15:21 (0.040, #26), HLA-B*40:06 (0.040, #27), HLA-C*12:02 (0.037, #28), HLA-B*35:01 (0.035, #29), HLA-A*26:01 (0.020, #34), HLA-B*18:03 (0.015, #38), HLA-B*35:02 (0.010, #47), HLA-B*42:02 (0.005, #57), HLA-B*35:17 (0.005, #54), HLA-B*40:02 (0.005, #56).

**Watchlist overlap (correction to the ruling):** epitopes adjacent to a qualifying AAY seam AND on the internal-cleavage watchlist = **4**: GPKEPFRDY, HWTTYMDTF, RFALNPGLL, VYSTCTVPTM (the ruling listed three; `RFALNPGLL` also qualifies at 5x). Supplementary observation, not an escalation.

### Class II: the seams are clean (positive result)

Null calibration: EL rank is a percentile against random peptides, so random-like seams would show ~10% at rank <= 10 and ~1% at rank <= 1. Over n = 82203 deep seam (window,allele) pairs: **7.17% (x0.72) and 0.76% (x0.76)** -- seam windows bind slightly WORSE than random.

| linker set | n pairs | rank<=10 | x null | rank<=1 | x null |
|---|---|---|---|---|---|
| GPGPG | 36234 | 5.86% | 0.59 | 0.34% | 0.34 |
| KK | 26730 | 6.49% | 0.65 | 1.04% | 1.04 |
| AAY | 15906 | 11.36% | 1.14 | 1.37% | 1.37 |
| EAAAK+adjuvant | 2475 | 4.89% | 0.49 | 0.20% | 0.20 |

**GPGPG seams are clean:** x0.59 at rank <= 10 and x0.34 at rank <= 1 (n = 36234), and 0 of the 193 class II contested rows at 5x involve GPGPG (class I: 1 of 31 windows). Contested class II rows by linker set at 5x: {'AAY': 83, 'KK': 108, 'EAAAK+adjuvant': 2} -- the residual class II contest sits at KK seams (flanking the B-cell epitopes) and AAY seams (flanking the class I epitopes), not around the class II epitopes' GPGPG seams. AAY is the only linker with mild class II enrichment (x1.14 at rank <= 10; windows overlap, so treat as indicative). At 5x, class II uncontested coverage is 98.9% (n = 43/81 pairs still contested on 18 alleles).

## Escalation (Opus: whole-construct uncontested <~70% or gap >~15 pp; trigger 7)

At the primary margin the whole construct is 99.88% uncontested (gap 0.12 pp). Fired: none. Trigger 7 (HIV vs Mpox combined): not fired. Trigger 4: not fired. **No open escalation.**

## Caveats

- Coverage is Southeast Asian (Austronesian) PROXY coverage from AFND Singapore Riau Malay (deviation #7), not Philippine.
- Hardy-Weinberg per locus; the union treats an allele as presented if ANY epitope presents it.
- Class II liberation filter is a proxy (proteasomal cleavage is not the class II pathway); class II margin uses EL-rank ratio (EL has no IC50).
- The class I margin uses the processing endpoint's IC50; real epitopes have wide IC50 (12.6-30,227 nM), so a 5x margin against a weak real epitope is not the same as a strong binder -- see the IC50<=500 nM sensitivity.
- Phase I MHC-II ranks (BA/consensus) and the EL ranks used here are different scales; never table them side by side.
