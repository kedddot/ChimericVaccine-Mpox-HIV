# Phase 1Ec — Human Self-Homology Screening Methodology Note

**Applies to:** Sec. I.E of the manuscript (BLASTP against the reviewed human
reference proteome, UniProt Homo sapiens Swiss-Prot subset).

## Two layers, and why there are two

**Layer 1 -- BLASTP, as published.** Every peptide in 1Eb's Filtered output is
BLASTP'd (`-task blastp-short`, the NCBI-recommended task for short queries)
against a local database built from `reviewed:true AND organism_id:9606`.
A peptide is flagged "self-like" and EXCLUDED only if its best hit meets ALL
THREE: percent identity >= 70%, alignment coverage
(alignment_length / peptide_length) >= 70%, E-value <= 1e-5 -- the published
acceptance rule, unchanged. Any hit below that bar is RETAINED but flagged
`Self_Homology_Status = PARTIAL` for manual review, per the paper's own text:
"Epitopes with partial homology below this threshold were retained but
flagged for manual review, given that even sub-threshold similarity to
self-proteins can, in some cases, still influence immunogenicity through
partial tolerance effects."

**MEASURED LIMITATION, stated plainly rather than left implicit: BLASTP
cannot reach statistical significance at this study's candidate lengths
(9-16 aa).** Tested on 20 peptides cut verbatim from human Swiss-Prot at
these exact lengths -- guaranteed perfect self-matches by construction --
default blastp settings found 0/20; blastp-short recovers 13/20 as
detectable hits, but 0/20 clear the E<=1e-5 bar (best achievable E-value on a
9-mer against this database is ~0.2). On this study's real candidate pool,
Layer 1 is expected to, and did, exclude 0 peptides. That is a property of
BLAST's significance statistics at this query length, not evidence the
candidates are free of human self-homology, and must be reported as a length
limitation rather than as a clean safety result.

**Layer 2 -- exact contiguous-identity screen (new).** The same technique
already used and validated for the allergen screen (deviation #27):
`phase1_common.human_self_homology()`, an 8-residue exact-substring lookup
against the same human reference proteome. This has no significance
threshold to fail and no length floor above the 8-residue window, so it
supplies what Layer 1 structurally cannot at these lengths.

Window=8 was NOT assumed -- it was calibrated against a scrambled-sequence
control (5 composition-matched random permutations per real candidate, same
length, same amino acid multiset, no biology): w=6 and w=7 show no
separation from the scrambled background (~81% flagged either way) and are
useless; w=8 separates real candidates from scrambled controls by roughly
3x and is used. Loosening BLASTP's own identity/coverage thresholds instead
of using exact matching was tried first and rejected on the same
real-vs-scrambled test: at 70%/70% with the E-value requirement dropped, the
real pool and the scrambled controls flag statistically indistinguishable
fractions (24.1% vs 25.2%) -- pure background matching human protein amino-
acid composition, not genuine self-homology.

A HIT in Layer 2 EXCLUDES, on the same evidence-over-score principle as the
allergen homology screen and deviation #1's toxicity rule: an exact 8+
residue match to a reviewed human protein is factual, not a prediction. A
Layer-2 NONE is not proof of non-self-similarity: only exact matches are
detected, and a peptide 80% similar to a human protein passes clean -- the
same caveat already stated for the allergen homology screen.

## Database provenance
`human_swissprot_db/human_swissprot.fasta`, fetched from
`rest.uniprot.org/uniprotkb/stream?query=reviewed:true+AND+organism_id:9606&format=fasta`
-- 20,431 reviewed human sequences. Layer 1 uses a `makeblastdb -dbtype prot`
index of this file; Layer 2 searches the same FASTA directly as a
delimiter-joined substring blob, so both layers run against the identical
reference set.
