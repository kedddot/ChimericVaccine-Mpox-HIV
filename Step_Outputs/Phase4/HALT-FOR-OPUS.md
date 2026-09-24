# Phase IV -- Halted Branches (escalated to Opus)

## Branch: Phase4D_proteasomalProcessing (interpretation of liberation + neoepitope findings)

> **RESOLVED 2026-09-21 by Opus ruling:** construct stays frozen (Phase I is final); weight the junction flags by allele frequency instead (Step 5 Part 2). Superseded by the `Phase4E_contested_coverage` branch below, which is OPEN.

**Trigger fired:** #3 -- Step 4 found a junction neoepitope out-scoring a real epitope and/or an epitope predicted not liberated

**What I was doing:**

Scoring the full 570-aa construct with IEDB processing (MHC-I, 8-11-mers) and NetMHCIIpan-EL (MHC-II junctions, 12-20-mers) against pre-registered definitions (see Phase4D_processing.py header). Data tables are complete and written; the DESIGN interpretation (are these junctions a vaccine-design defect?) is paused.

**Exact numbers / output:**

- junction neoepitopes out-score real epitopes: MHC-I 126 (window,allele) rows / 32 windows (103 beat ALL real epitopes at that allele); MHC-II 4022 rows (3033 with core spanning the junction)

Top MHC-I flags:
  TYMDTFAAY @ HLA-C*14:02 total=2.0119 IC50=4.7nM segs=MHC-I:HWTTYMDTF|AAY beats_all=YES
  CTVPTMAAY @ HLA-A*26:01 total=1.9634 IC50=6.2nM segs=MHC-I:VYSTCTVPTM|AAY beats_all=YES
  EPFRDYAAY @ HLA-B*35:01 total=1.688 IC50=6.7nM segs=MHC-I:GPKEPFRDY|AAY beats_all=YES
  FSATVYAAY @ HLA-B*35:01 total=1.6498 IC50=7.0nM segs=MHC-I:QTSVFSATVY|AAY beats_all=YES
  FSATVYAAY @ HLA-B*35:05 total=1.4949 IC50=10.0nM segs=MHC-I:QTSVFSATVY|AAY beats_all=YES
  EPFRDYAAY @ HLA-B*35:17 total=1.4349 IC50=12.0nM segs=MHC-I:GPKEPFRDY|AAY beats_all=YES
  TCTVPTMAAY @ HLA-A*26:01 total=1.43 IC50=21.0nM segs=MHC-I:VYSTCTVPTM|AAY beats_all=YES
  FSATVYAAY @ HLA-B*35:17 total=1.4015 IC50=12.4nM segs=MHC-I:QTSVFSATVY|AAY beats_all=YES

Top MHC-II flags:
  PQDLAAYNKRKR @ HLA-DRB1*13:01 rank=0.01 core=LAAYNKRKR spans=YES beats_all=YES
  QDLAAYNKRKRV @ HLA-DRB1*13:01 rank=0.01 core=LAAYNKRKR spans=YES beats_all=YES
  EQKAYKKANASA @ HLA-DPB1*09:01 rank=0.01 core=KAYKKANAS spans=YES beats_all=YES
  PEQKAYKKANASA @ HLA-DPB1*09:01 rank=0.01 core=KAYKKANAS spans=YES beats_all=YES
  TPQDLAAYNKRKRV @ HLA-DRB1*13:01 rank=0.01 core=LAAYNKRKR spans=YES beats_all=YES
  PQDLAAYNKRKRVI @ HLA-DRB1*13:01 rank=0.01 core=LAAYNKRKR spans=YES beats_all=YES
  QARVLAVERKKGGK @ HLA-DRB1*11:04 rank=0.01 core=VLAVERKKG spans=YES beats_all=YES
  PEQKAYKKANASAQ @ HLA-DPB1*09:01 rank=0.01 core=KAYKKANAS spans=YES beats_all=YES

**Options I see:**

1. Accept as expected: any 570-aa concatenation of binders creates some seam binders; report counts and severity tiers, no redesign.
2. Treat beats-ALL-at-allele deep junctions with a spanning core as a design finding for the manuscript (Discussion/limitations).
3. Redesign the affected linker(s) (e.g. swap AAY/GPGPG at the flagged seams) and re-run Steps 4-5 -- this would change the frozen construct, so it is Opus/PI's call.

_Halted: 2026-09-21T18:44:31_

------------------------------------------------------------------------------------------

## Branch: Phase4E_contested_coverage

> **RESOLVED 2026-09-21 by Opus's second ruling:** the 'contested' definition was a specification flaw (extreme-value comparison); class II contest is an artifact, class I contest is real but a Discussion/limitations finding; no redesign; construct stays frozen. Step 5 re-issued under the corrected definition: no escalation condition met at the primary 5x margin. The figures quoted in this section are SUPERSEDED and must not be reported.

**Trigger fired:** #3 -- Step 5 escalation condition from the Opus ruling was met

**What I was doing:**

Frequency-weighted junction analysis (Part 2) and level 3 per-pathogen coverage.

**Exact numbers / output:**

- CONSTRUCT (21 MHC epitopes) [UNFILTERED] LITERAL (combined): uncontested 81.33% (raw 100.00%, gap 18.67 pp)
- MHC-I only (10) [LIBERATED] informative layer (per class): uncontested 81.93% (raw 98.74%, gap 16.81 pp)
- MHC-I only (10) [UNFILTERED] informative layer (per class): uncontested 80.99% (raw 98.74%, gap 17.75 pp)
- MHC-II only (11) [LIBERATED] informative layer (per class): uncontested 22.06% (raw 100.00%, gap 77.93 pp)
- MHC-II only (11) [UNFILTERED] informative layer (per class): uncontested 1.79% (raw 100.00%, gap 98.21 pp)

CALIBRATION (supplementary, no new network calls):
- MHC-II seam windows are NOT enriched over random peptides: rank<=10 7.17% (null 10%), rank<=1 0.76% (null 1%), n=82203 (window,allele) pairs; 2491 deep windows scanned per allele vs 11 reference epitopes -> the class II 'contest' is dominated by an extreme-value effect.
- MHC-I seam windows ARE mildly enriched: IC50<=500nM in 0.75% (n=38924) vs 0.42% (n=10804) of unselected natural adjuvant-domain windows (x1.81); reference is a single 45-aa protein, treat as indicative.
- Margin sensitivity (LIBERATED, uncontested coverage; construct / MHC-I / MHC-II): 1x: 85.9/81.9/22.1; 2x: 88.7/83.6/31.3; 5x: 99.7/84.2/97.9; 10x: 100.0/88.6/99.7 -> MHC-I contest persists at every margin; MHC-II contest evaporates by 5x.

**Options I see:**

1. Rule on the DEFINITION: the spec'd 'any junction out-ranks it' (margin 1x) is dominated by extreme-value chance for class II; adopt a margin (e.g. >=5x) or a null-calibrated criterion for MHC-II.
2. Treat the MHC-I finding (robust, mechanistic: AAY ends in Tyr, an MHC-I anchor) as the real design finding; decide between reporting it as a limitation and redesigning the AAY linker (changes the frozen construct -- Opus/PI only).
3. Note the literal combined figure is saturated by class II compensating; MHC-I-only coverage falls 98.7% -> ~82-89% depending on margin.

_Halted: 2026-09-21T19:03:33_

------------------------------------------------------------------------------------------

