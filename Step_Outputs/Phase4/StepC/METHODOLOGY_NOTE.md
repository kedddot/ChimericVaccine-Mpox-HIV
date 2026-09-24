# Phase 4C -- TCR Recognition: methodology notes

Generated: 20260921_1949

## Ladder used (pre-approved, Section D/F)

Methodology IV.A.2 names "MixMHCpred + TCGA Contact Database" -- both wrong
(disagreement #6: MixMHCpred predicts MHC BINDING, not TCR contact; TCGA is
a cancer genomics atlas with no TCR-contact data). IEDB's own
`immunogenicity/` endpoint returns 403. **Retested 2026-09-21:** 11 patient attempts over several minutes (including after a
90 s wait, and on GET) all returned 403 while `mhci/`, `mhcii/` and `bcell/` answered 200 in the same minute, and a nonexistent path
(`foobar/`) returns the byte-identical 403 page -- i.e. the endpoint is not served at this URL; this is **not** IP throttling
(throttling also hit healthy endpoints and cleared on retry). The Calis reimplementation was therefore forced, not merely convenient.
(Other, transient 403s elsewhere in this project were genuine throttling; the two produce the same HTML page and can only be told apart by such a retest.) This step uses the pre-approved
3-rung ladder instead:

1. **BigMHC-IM** (primary, MHC-I only) -- 10/10 epitopes scored.
2. **PRIME 2.0 + MixMHCpred 3.0** (optional, MHC-I only) -- INSTALLED AND USED.
3. **Calis et al. (2013)**, reimplemented locally from the paper's own
   published data (see below) -- all 10 epitopes scored.

**Scope: MHC-I only (10 construct epitopes).** All three ladder tools are
MHC-I specific by design (BigMHC/PRIME/MixMHCpred explicitly; Calis's model
was trained on and validated against MHC-I presented peptides). Section F's
own text confirms this scope ("TCR-facing residues = positions 4-6 of each
MHC-I epitope"). MHC-II and B-cell construct epitopes are out of scope for
this ladder -- there is no equivalent pre-approved TCR-recognition tool for
CD4 T-cell or antibody epitopes in this project.

## PRIME/MixMHCpred installation

Installed fresh this run from the GitHub repos (GfellerLab/PRIME, GfellerLab/MixMHCpred), both free/no licence wall, matching Section F's description. Both ship a Mach-O arm64 binary (`PRIME.x`) that runs natively on this machine -- no compilation needed. Two environment issues were hit and fixed:

1. **Both tools' own launcher scripts refuse a space in their install path** ("Spaces in path to MixMHCpred are not supported") -- and this project's root sits under a space-containing macOS volume name ("Extended SSD"). Fixed by invoking both tools through the space-free symlink this project already maintains for exactly this class of problem (`~/mpoxhiv_ssd`, see `~/mpoxhiv_env.sh`'s own comment on why it exists) instead of the direct `/Volumes/Extended SSD/...` path -- same files, different path, nothing moved or duplicated.
2. **Both tools invoke a bare `python3`** for their own internal scoring script, which resolves to the system Python (no pandas installed) rather than the project's `phase2` conda env. Fixed by prepending `$MPOXHIV_ENVS/phase2/bin` to `PATH` for the subprocess call only.

Neither fix required editing PRIME's or MixMHCpred's own files -- both are external, third-party tools and Rule 1 (never edit Phase 1/2) extends in spirit to never editing vendored external tools either.

## Direction checks (mandatory before trusting any score)

Each tool's score was checked against Phase4B's MHC-I `NetMHC_EL_Best_Rank`
(a quantity already validated in Step 2) for the expected sign BEFORE being
reported as a result -- this discipline was added after Phase4B's IC50
parsing error (misreading an unlabeled column with the wrong sign) was
caught and corrected.

| Tool | Convention | Expected vs rank | Observed r | p (two-sided) | n |
|---|---|---|---|---|---|
| BigMHC-IM | 0-1, higher=more immunogenic | negative | -0.298 | 0.402 | 10 |
| PRIME %Rank | lower=better | positive | +0.916 | 0.0002 | 10 |
| Calis score (EXPLORATORY) | higher=more immunogenic | no prior (peptide-intrinsic) | +0.097 | 0.790 | 10 |

**Read these with n in mind.** At n=10 only |r| > ~0.63 is distinguishable
from zero at p<0.05. **PRIME's r=+0.916 is meaningful; BigMHC-IM's
r=-0.298 (p=0.40) is correct in sign but NOT statistically
distinguishable from zero** -- it passes the sign check, it does not
corroborate anything. All observed signs matched expectation.

**Same-allele determinism check (replaces an earlier, invalid comparison).**
An earlier draft reported r=+0.460 between this run's BigMHC-IM "best" and
Phase1Db's stored `BigMHC_IM_Score` as a "same-tool consistency check". That
was wrong: Phase 1 stored the max over ITS OWN Binding_Alleles set (9-allele
panel), this run's best is over a different (107-allele-panel) set for 9 of
10 epitopes -- two different quantities. The honest check re-scores Phase 1's
exact (allele, peptide) pairs and compares the max over that same set:
r=+1.0000, n=10, max |diff| = 3.97e-07. BigMHC-IM
reproduces Phase 1 to numerical precision on identical inputs; the difference
in the "best" columns is entirely because the 107-allele panel found better
alleles, which is correct.

**BigMHC-IM and PRIME choose their "best" allele independently.** They agreed
on the same allele for only 3/10 epitopes
(`BigMHC_PRIME_Best_Allele_Match` column). Their "best" columns are therefore
NOT allele-matched and must never be compared row-by-row as if they were.

**Calis is exploratory only.** Self-AUC 0.62 on its own training set, r=+0.097
here. Near-zero discrimination against binding rank is arguably correct --
the model is designed to be independent of binding affinity -- but that also
means it is NOT independent corroboration of BigMHC-IM or PRIME and must not
be presented as such.

## Flagged finding carried forward: NKRKRVIGL (Mpox A35R)

Worst in the construct on all four independent signals at once
(`Worst_In_Construct_On` column): PRIME %Rank 4.386 (next worst 0.47),
BigMHC-IM 0.0344, binding breadth 2/107 alleles, EL rank 0.85. Also note
Phase 1 scored it at `HLA-B*08:01`, which is not in the 107-allele panel at
all (0.163 there; within the local panel its best BigMHC-IM is 0.0344). This
is consistent with A35R's lowest-in-table BigMHC score in the paper and with
the `A35R_Coverage_Note` substitution recorded in Phase1G. Expect it to add
almost nothing to Step 5's coverage union. **Step 7 must report this as a
consistent cross-phase signal, not a new anomaly.**

## Calis et al. (2013) reimplementation -- full provenance

Citation: Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette
A, Kesmir C, Peters B (2013) "Properties of MHC Class I Presented Peptides
That Enhance Immunogenicity." PLOS Comput Biol 9(10):e1003266.

**The amino-acid log-enrichment table is NOT taken from memory or a
secondary source.** It is derived here directly from the paper's own public
Supplementary Dataset S1 (`pcbi.1003266.s001.xls`, downloaded from PLOS;
the file's own embedded metadata credits "Jorg Calis" as author, confirming
authenticity) -- 2,508 labelled peptides (2,167 immunogenic / 341
non-immunogenic; lengths 8/9/10), using the paper's own documented method
(quoted from the PMC full text, PMC3808449): "the enrichment is calculated
as the ratio between the fraction of that amino acid in the immunogenic
versus non-immunogenic data sets," pooled over non-anchor positions
(anchors = positions 1, 2, and the C-terminus -- "P1, P2 and P9 for most
HLA molecules," per the paper's own text; a generic, allele-agnostic mask,
since per-allele anchor definitions are not published in a reusable table).

Per-position importance weights are the paper's own **Table 2** values
(quoted exactly from PMC3808449): position 3=0.10, 4=0.31, 5=0.30, 6=0.29,
7=0.26, 8=0.18 (positions 1, 2, 9 are anchors, weight 0). For the one extra
internal position in a 10-mer (beyond position 8), this implementation
reuses the position-8 weight (0.18) -- Table 2 does not tabulate further,
and this is a documented extension, not a published value.

**Validation before use:** scoring all 2,508 training peptides with this
reimplementation gives a self-evaluated AUC of 0.62 (not held-out -- a
sanity check, not a performance claim) and the correct-direction mean score
(immunogenic peptides score higher than non-immunogenic on average). The
derived amino-acid table's qualitative pattern also matches the paper's own
reported findings: tryptophan (+0.317) and other aromatic/large residues
(phenylalanine +0.161, tyrosine +0.134) are the most enriched, serine
(-0.395) is the single most depleted -- both independently reported in the
paper's text.

## TCR-facing residues

Positions 4-6 (1-based) of each MHC-I epitope, per Section F's own
specification. The methodology's "\>=3 TCR-contact residues" exclusion rule
cannot be applied without a defined position set (it never specifies one),
so this step reports which residues occupy 4-6 rather than scoring against
an assumed threshold.

## Cross-reactivity (disagreement #7 -- both layers, never one alone)

All 31 construct epitopes (not just the 10 MHC-I ones -- self-tolerance is
a safety question independent of TCR-recognition tooling). Two layers,
reported together, never one presented alone as a safety verdict:

- **BLASTP layer**: reused directly from Phase1Ec (frozen, Rule 1 -- not
  recomputed). Phase1Ec's own BLASTP screen against reviewed human
  Swiss-Prot is the same screen Phase 1Ec_Filtration measured as blind at
  9-16 aa (0/20 known human self-fragments caught in that validation).
- **Exact 8-mer layer**: computed FRESH this run via
  `phase1_common.human_self_homology()` (imported, not reimplemented, per
  Section F's explicit instruction), and cross-checked against Phase1Ec's
  own `SelfHomology_Exact_Match`/`SelfHomology_Exact_Kmer` columns (already
  in the Phase4A dossier) for consistency: 31/31 agree.
