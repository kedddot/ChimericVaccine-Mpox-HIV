# Phase 4D -- Proteasomal Processing: methodology notes

Generated: 20260921_1949

## STATUS: data complete; DESIGN INTERPRETATION PAUSED (escalation trigger 3 fired)

See `Step_Outputs/Phase4/HALT-FOR-OPUS.md`. Every table below is complete and final for the definitions stated;
whether the flagged seams are a vaccine-design defect is a judgment call reserved for Opus/PI. Nothing here is a design verdict.

## What the processing endpoint returns (verified live before use -- direction discipline from Step 2)

- `total_score = proteasome_score + tap_score + mhc_score` and `mhc_score = -log10(IC50 nM)`: exact identities
  (max residual 0.0001). **Higher = better for every column.**
- `proteasome_score` is a pure **C-terminal cleavage score**: identical for every window ending at the same position
  (max deviation across all alleles and lengths = 0.0000, 165641 comparisons), so it defines one cleavage profile
  for the whole construct. `tap_score`/`processing_score` are allele-independent (max deviation 0.0000); only `mhc_score` varies by allele.
- Direction checks: mean cleavage score by C-terminal residue, high to low = **WLFYMAIRDHCVSNTEQKPG** (W/L/F/Y lead; P/G trail --
  matches known proteasome specificity). Endpoint IC50 vs Step 2's per-allele EL rank for the same epitope-allele pairs:
  Spearman rho = **+0.468** (p = 4.5e-05, n = 70 pairs from 10 epitopes -- pairs are not independent), correct (positive) sign.
  (An earlier "<=500 nM" consistency check was an unvalidated prior of mine and has been removed.)
- Request shape: `method=netmhcpan`, whole 570-aa construct, one call per (allele, length 8-11); **no `tap=` parameter** (rejected).
  All 74 MHC-I alleles of Phase 1F's ALLELE_FREQ x 4 lengths = 296/296 tables.

## Scope

The endpoint is MHC-I only. Liberation verdicts cover the 10 MHC-I construct epitopes. MHC-II (endolysosomal pathway) and
B-cell (native-structure recognition) epitopes get boundary cleavage scores as supplementary context; verdict N/A.
MHC-II junction peptides (12-20 aa) are scored with NetMHCIIpan-EL (same method and allele bridging as Step 2), and the 11 real MHC-II
epitopes are scored in the SAME calls so comparisons are like-for-like. The toxicity/allergen half of methodology IV.B.2 was done in
Phase 2A (`Step_Outputs/Phase2/StepA/`) and is referenced, not re-run.

## Definitions fixed BEFORE results were seen

- NOT liberated: C-terminal proteasome_score below the construct-wide median end-position score (percentile < 50).
- Junction window: spans >=1 segment boundary. Deep: >= 3 residues each side of every crossed boundary (fewer = a shifted copy of a real epitope, reported but excluded from the trigger).
- Neoepitope flag: deep AND strong binder (IC50 <= 500 nM MHC-I; EL rank <= 10 MHC-II) AND out-scores >= 1 real construct epitope that binds the SAME allele.
  "Beats ALL" = out-scores every real epitope binding that allele. Comparators are only epitopes Step 2 called binders of that allele.

## Results

- **Liberation (MHC-I):** 10/10 liberated by the pre-registered criterion; not liberated: none.
  This part of trigger 3 did NOT fire. C-terminal cleavage percentiles: {'APGSPTNLEF': 96.3, 'GPKEPFRDY': 84.7, 'HHFNCRGEF': 84.9, 'HWTTYMDTF': 97.7, 'NKRKRVIGL': 98.2, 'QTSVFSATVY': 97.5, 'RFALNPGLL': 96.1, 'SEGATPQDL': 93.6, 'TMGAASITL': 99.8, 'VYSTCTVPTM': 89.5}.
- **Secondary observation (not part of the trigger):** 5 MHC-I epitopes have an internal cleavage site scoring higher than their own C-terminal one
  (GPKEPFRDY, HHFNCRGEF, HWTTYMDTF, RFALNPGLL, VYSTCTVPTM) -- a possible destruction risk that the pre-registered criterion does not test.
- **MHC-I junctions:** 1475 windows span a boundary, 526 deep. **126 (window, allele) flags on 32 distinct deep windows; 103 beat ALL real epitopes at that allele.**
  Shallow (shifted-copy) flags excluded: 138. Linkers involved (distinct flagged windows): {'AAY': 26, 'EAAAK': 1, 'GPGPG': 4, 'KK': 1}.
  Mechanistic observation only: the AAY linker ends in Tyr, a classic MHC-I C-terminal anchor and a favoured cleavage residue, so epitope-end+AAY windows are strong binders that are also efficiently cleaved.
- **MHC-II junctions:** 4465 windows (12-20 aa), 33 alleles. **4022 flags on 1117 windows; 3033 have the predicted binding core spanning the junction; 1757 beat ALL real epitopes at that allele.** Shallow flags excluded: 4052.

## Caveats for whoever interprets this

- "Out-scores a real epitope" is a weak bar where an allele has few or weak real comparators; the `Real_Epitopes_At_Allele`, `Weakest_Real_Total`/`Best_Real_Rank` columns let you judge each row.
- MHC-II core-spanning is inferred from NetMHCIIpan's reported core located in the window (first occurrence).
- EL rank (MHC-II) and processing IC50/total (MHC-I) are different scales; do not compare across classes. Phase I MHC-II ranks (BA/consensus scale) are not comparable to these EL ranks (Step 7 standing instruction).
- Weighted by allele frequency (Step 5's population coverage), most of these seam binders may matter far less than counts suggest; counts here are per (window, allele).

## Run integrity (what went wrong on the way, kept for the record)

A first 4-worker run was throttled by IEDB (HTTP 403): 64/296 MHC-I tables (every HLA-C allele and a few others) and ~220 MHC-II calls failed silently
and the script analysed the partial data. Every number from that run (e.g. 111 MHC-I flags on 60/74 alleles) was DISCARDED. Direct tests confirmed the
"missing" alleles are supported. The script now runs 2 workers with patient retries and **refuses to analyse any incomplete scan** (exit code 2).
All figures above come from the complete run (296/296 MHC-I tables, 297/297 MHC-II calls, zero failed calls).
