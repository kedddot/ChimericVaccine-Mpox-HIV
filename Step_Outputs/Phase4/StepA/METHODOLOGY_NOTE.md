# Phase 4A -- Epitope Dossier: methodology notes

Generated: 20260920_1855

## Deviations from a literal reading of Section IV

1. **"194 Non-Redundant Sequences" is a mislabel.** Section IV's own text calls
   194 the epitope count. 194 is actually the number of source PROTEIN
   sequences retrieved in Phase 1A (across all HIV/Mpox targets and variants),
   not the number of epitopes analysed. Phase IV's real analysis units are the
   **31 epitopes** in the final construct Vax_Final_4fce119e (this dossier) and,
   from Phase 4B onward, the **105 candidates** that survived to Phase 1F
   coverage scoring. See MANUSCRIPT_CORRECTIONS_PHASE4.md for the corrected text.

2. **Mpox_B5R / B6R naming.** `Epitope_Provenance` in the Phase 1G construct
   CSV labels the Mpox EEV-glycoprotein epitopes' target as `Mpox_B5R`. The
   manuscript and Phase 1A's antigen-identity check both call the correct
   ortholog **B6R** (OPG190; VACV B5R ortholog -- MPXV's own gene literally
   named "B5R" is a different protein, OPG189, an ankyrin repeat protein; see
   Phase 1A_provenanceCorrection.py's ANTIGEN_REFERENCES). Both names are kept
   on every row (`Labelled_Target` = as-recorded, `Antigen_Display` = with the
   manuscript name) so neither reading is lost.

## Join method

All joins are by **peptide sequence** (never row number) against:
- Phase1G (`Epitope_Provenance`, `Boundary_Map`, `Sequence`) -- construct membership, class, position
- Phase1Db Elite -- percentile rank, BepiPred, BigMHC-IM, SEMA overlap
- Phase1Dc Min_50pct -- conservancy
- Phase1F Elite (Filtered) -- overall population coverage
- Phase1Ec Filtered -- self-homology status (cross-checked against Phase1F's
  own copy of the same columns; any disagreement is logged as a [WARNING] above)

## Assertions enforced

31 total rows; 10 MHC-I / 11 MHC-II / 10 B-cell; 14 HIV-derived / 17
Mpox-derived; 7 distinct source antigens; zero unjoined epitopes across all
four source files. A failure on any of these halts the step and escalates
per Section H trigger #2 (join/count failure), rather than silently
proceeding with partial data.
