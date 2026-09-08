ALLERTOP v2 AND ALLERGENFP v1 TRAINING SETS
===========================================
Downloaded by the user from the tools' own websites, 2026-08-31.

  allertop_allergens.fasta        2,426   AllerTOP v2 positive class
  allertop_non_allergens.fasta    2,426   AllerTOP v2 negative class
  allergenfp_allergens.fasta      2,391   AllergenFP v1 positive class
  allergenfp_non_allergens.fasta  2,363   AllergenFP v1 negative class

The published AllerTOP v2 paper states 2,427 + 2,427; the distributed files
hold 2,426 each.

WHAT THESE ARE FOR IN THIS PROJECT. Two things, and NOT a third:

  1. PROOF OF THE ALLERGENFP 16-AA FLOOR (deviation #6).
     The AllergenFP sets are strict SUBSETS of the AllerTOP sets, and the
     dropped records are exactly the short ones:
         dropped from allergens      35 records, every one  6-15 aa
         dropped from non-allergens  63 records, every one  7-15 aa
         shortest record KEPT in either AllergenFP file: exactly 16 aa
         records <16 aa kept: 0 of 4,754
     The authors removed every sequence below 16 residues from their own
     training data. Deviation #6 recorded the 16-aa floor as OUR empirical
     observation because it is not stated in the AllergenFP paper. It is now
     demonstrable from the authors' own distributed dataset.

  2. MEASURING WHETHER ALLERTOP IS VALID AT EPITOPE LENGTHS (deviation #27).
     AllerTOP v2 is a k-nearest-neighbour classifier, so its training set is
     effectively its model, and the method can be reimplemented. See
     Phase 1/STEP E/_tool_validation/allertop_v2_reimplementation.py and the
     REBUILD_LOG entry of 2026-08-31 12:30 for the measurements.

  3. *** NOT USED TO GENERATE THE STUDY'S ALLERTOP VERDICTS. ***
     The reimplementation reproduces the PUBLISHED ACCURACY but not
     necessarily the SERVER'S INDIVIDUAL CALLS: the ACC lag is not stated in
     the paper, and lag 5 versus lag 7 flips 57 of this study's 184 candidates
     (31%) while both settings reproduce the published accuracy equally well.
     A per-peptide verdict from it is therefore not defensible as "AllerTOP
     v2.0 output". The cited tool's own server remains the source of record.
