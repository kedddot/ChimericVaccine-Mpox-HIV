ALLERGEN REFERENCE SET
======================
uniprot_reviewed_allergens.fasta

Source : UniProtKB, keyword KW-0020 ("Allergen"), reviewed (Swiss-Prot) only.
Query  : https://rest.uniprot.org/uniprotkb/stream?query=keyword:KW-0020%20AND%20reviewed:true&format=fasta
Fetched: 2026-08-31

WHY REVIEWED-ONLY: Swiss-Prot entries are manually curated. Adding unreviewed
TrEMBL entries would inflate the reference with unverified annotations and make
a homology hit less meaningful.

WHAT IT IS FOR: the FAO/WHO Codex allergenicity criterion asks whether a query
shares a contiguous identical stretch with a KNOWN allergen. That is a factual
lookup, not a prediction, and unlike similarity/ML methods it has no minimum
length -- it works on a 9-mer.

MEASURED BEHAVIOUR (2026-08-31, this reference set):
  positive control, peptides cut FROM allergens   9/10/15/16 aa -> 100% detected
  negative control, peptides cut from human proteome            -> 0.0-2.0%
  the study's own viral candidates (6,000 sampled)              -> 0.00% at w=8
  window sensitivity: w=6 flags 6.1% of HUMAN peptides (over-sensitive, the
  documented weakness of the classic 6-mer rule); w=8 flags 0.4%.

LIMITS TO STATE WHEREVER THIS IS REPORTED:
  1. A negative means "no match in THIS reference set", not "safe". It cannot
     detect a novel allergen absent from the database.
  2. Exact matches only -- a peptide 80% similar to an allergen passes clean.

POSSIBLE UPGRADE: AllergenOnline (FARRP, ~2,400 entries) and COMPARE (~2,600)
are larger curated sets, both freely downloadable after registration. A larger
reference makes a NEGATIVE stronger, since absence of a match is only as
informative as the set searched.
