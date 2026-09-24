# Phase 4B -- MHC Binding Affinity: methodology notes

Generated: 20260920_1954

## What was scored

The 105-candidate Phase 1F pool ({'MHC-I': 37, 'MHC-II': 29, 'B-cell': 39}), scored once against IEDB's
`mhci/` (MHC-I, lengths 9/10) and `mhcii/` (MHC-II, length 15) endpoints, each
with both `netmhcpan_el` and `consensus`. The 31 construct epitopes are a
verified subset of this pool (Phase4A), so no epitope is scored twice.

## Allele panel (disagreement #4)

All 107 alleles from Phase 1F's `ALLELE_FREQ` (74 MHC-I + 33 MHC-II), imported
directly -- not retyped, not the methodology's "10 Philippine alleles" (that
table's own top allele, `A*24:23` at 18.5%, doesn't parse as a real allele
and matches nothing in this project; `ALLELE_FREQ`'s provenance is AFND
Singapore Riau Malay, a documented proxy -- see Phase1F_coverage.py header).

MHC-II's DQB1/DPB1 entries in `ALLELE_FREQ` are single-chain (beta only);
IEDB's `mhcii/` requires paired alpha/beta notation. Bridged using
Phase1Db_filtration.py's existing `MHCII_ALLELE_TO_SINGLE_CHAIN` mapping
(inverted) plus 5 additional pairings for alleles outside Phase1Db's smaller
21-allele panel, using the standard default DP alpha (DPA1*01:03) and
documented DQ haplotype linkage for the rest. All 5 were verified live
against IEDB (200 OK, valid scored response) before use. Any allele that
still fails to map or persistently errors is logged and excluded --
contributes 0 to Binding_Alleles for every epitope, same handling as
Phase1F's own MHCflurry no-percentile-calibration alleles.

## Gating (disagreement #1)

Percentile rank, at Phase 1's own thresholds: <=1.0 for MHC-I, <=10.0 for
MHC-II (`netmhcpan_el` / `netmhciipan_el`, matching what Phase1Db_filtration.py
already uses). IC50 (from the `consensus` call: `ann_ic50` for MHC-I,
`nn_align_ic50` for MHC-II, at whichever allele gave the primary-gate best
rank) is reported as an EXTRA column, banded Strong<=500nM / Weak
500-5000nM / Non-binder>5000nM, and never gates.

`netmhcpan_el`/`netmhciipan_el` (elution-score methods) carry no IC50 column
at all -- confirmed by a live test call before writing the parser. IC50
always comes from the paired `consensus` call.

**Finding: `consensus` supports a much narrower allele set than
`netmhcpan_el`/`netmhciipan_el`.** Measured directly from this run's cached
responses: `consensus` returned a valid scored result for only ~31/74
(~42%) of the MHC-I alleles and ~18/33 (~55%) of the MHC-II alleles in
`ALLELE_FREQ` -- every other allele came back `Invalid allele name X found`
from IEDB itself (consensus predates NetMHCpan/NetMHCIIpan's much broader
pan-specific HLA coverage; this is a known, documented IEDB limitation, not
a bug in this pipeline's allele-name formatting -- confirmed by the SAME
allele name succeeding under `netmhcpan_el`/`netmhciipan_el` in the same run).
Consequence: for 35
of the 66 MHC-I/MHC-II
pool epitopes, the primary gate's best `netmhcpan_el`/`netmhciipan_el`
allele has no `consensus` result, leaving `IC50_At_Primary_Allele_nM` = "N/A"
for that epitope even though `consensus` DID return a usable IC50 at some
other allele. `IC50_At_Best_Consensus_Allele_nM` / `IC50_Band_Best_Consensus`
report that instead (IC50 at whichever allele `consensus` itself ranked
best), so IC50 reporting isn't needlessly sparse. This never affects the
primary gate, which stays on `netmhcpan_el`/`netmhciipan_el` at the full
107-allele panel regardless of `consensus` coverage. Blank cells were
changed to explicit `"N/A"` (never truly blank) in both IC50 columns and in
the B-cell primary-gate columns, so nothing in this table can be misread as
`UNRESOLVED`.

### CAUGHT-AND-CORRECTED PARSING ERROR: an unlabeled column was misread as IC50

MHC-II `consensus` responses carry an extra, unlabeled 4-column block IEDB's
own header doesn't name (verified: 1224/1224 cached MHC-II `consensus` rows
have 24 tab-separated fields against a 20-column header). For most rows this
trailing block is all `-`. For a handful of rare alleles (e.g.
`HLA-DRB1*13:01`) the three NAMED methods (comblib/smm_align/nn_align) are
all `-` and this trailing block holds the only non-dash values in the row.

**An earlier version of this step read that trailing block's 2nd field
(row-relative index -3) as a 4th `nn_align_ic50` fallback and reported it in
nM.** This was wrong and has been reverted. It was caught by review, then
proven wrong on four independent grounds (all re-derived directly from the
cached responses, not taken on faith):

1. **Range is impossible for IC50.** The column spans **-12.5 to +5.2**
   across 612 real (non-dash) values in this cache. IC50 is a concentration
   in nM and cannot be negative.
2. **Wrong-signed correlation with rank.** Pearson r between this column and
   the row's own `percentile_rank` is **-0.644** (n=612). A genuine IC50
   must correlate *positively* with rank (lower IC50 = tighter binding =
   lower/better rank, so both should move together) -- a strong negative
   correlation means this column moves the *opposite* direction of an IC50.
3. **Direction check confirms it's a "higher-is-better" score, not a
   concentration.** Strong-rank rows (rank 0.62-1.2, i.e. the best binders)
   sit at **+3.9 to +5.2**; weak-rank rows (rank 90-98, i.e. non-binders)
   sit at **-1.3 to -3.85**. An IC50 would run the other way (low nM for
   strong binders, high nM for weak ones) -- this column is some kind of a
   log-odds or z-scored prediction score from an unidentified 4th method,
   not a concentration in any unit.
4. **The genuinely-labeled `nn_align_ic50` column, read at its own correct
   header position in the same cache, looks nothing like it:** 1020 real
   values spanning **3.3 to 36040.7 nM, median 1938.8 nM** -- normal IC50
   magnitudes, all positive, none near the bad column's -12.5..+5.2 range.

**Consequence of the bug while it was live:** the Strong<=500nM band applied
to a column ranging -12.5..+5.2 marked essentially every row `Strong`,
including rows whose real (header-correct) methods were all non-binders.
`IC50_Band_Best_Consensus` was corrupted for every row that hit this
fallback, not only the 3 that surfaced it (those 3 were simply the only
rows where the NAMED methods were all dash, making the wrong value visible
as a non-blank cell instead of silently overriding a real one).

**Fix applied:**
- The trailing-column fallback is removed entirely. Unlabeled columns are
  never read as typed (nM) values anywhere in this step.
- IC50 (`IC50_At_Primary_Allele_nM`, `IC50_At_Best_Consensus_Allele_nM`,
  `Bcell_MHCII_IC50_Supplementary_nM`) now comes ONLY from the two
  genuinely nM-labeled columns per class: `ann_ic50` then `smm_ic50` for
  MHC-I, `nn_align_ic50` then `smm_align_ic50` for MHC-II (`genuine_ic50()`
  helper). Where both are `-` for the chosen allele, the value is `"N/A"` --
  honestly "no IC50-producing method ran for this allele/peptide", not a
  fabricated number.
- `KNKRKRVIGLCIRIS`, `NKRKRVIGLCIRISM`, `KRKRVIGLCIRISMV` (the 3 that
  surfaced this) now correctly read `IC50_At_Best_Consensus_Allele_nM =
  "N/A"` -- their only real methods (comblib/smm_align/nn_align) all
  returned dash for `HLA-DRB1*13:01`, and this step no longer invents a
  value from the unlabeled block.
- The unlabeled block itself is NOT exposed as its own column. Its
  identity (which 4th method it represents) is unconfirmed, and adding an
  unvalidated score column carries the same risk this whole finding is
  about -- not worth it for 3-of-105 rows.
- **Corrected, honest IC50 coverage after the fix**
  (`IC50_At_Best_Consensus_Allele_nM`, the fuller of the two columns):
  **MHC-I 37/37** (100% -- every MHC-I epitope has a genuine `ann_ic50` or
  `smm_ic50` at consensus's own best-ranked allele), **MHC-II 26/29**
  (~90% -- 3 epitopes honestly `N/A`). This is lower than the number the
  bug produced and is the correct number.

## Binding breadth vs. primary-gate pass/fail (do not report pass/fail as
## the headline result)

With a 107-allele panel, "best rank <= threshold across ANY allele" is a
near-vacuous pass criterion -- 37/37 MHC-I and 29/29 MHC-II pool epitopes
pass, which says almost nothing on its own (with enough alleles tried,
nearly every peptide binds something). **The informative metric is breadth**
-- `N_Binding_Alleles_Primary`, already computed per epitope:
- MHC-I: median 7 alleles bound (of 107), range 1-40 (n=37)
- MHC-II: median 6 alleles bound (of 107), range 1-29 (n=29)

Some epitopes bind only 1/107 alleles and contribute almost nothing to
population coverage despite technically "passing" the gate -- Phase4E (Step
5, population coverage) surfaces this numerically per epitope. Pass/fail
counts are reported in the console log for context only, never as a result.

## B-cell epitopes (10 in the construct, 16-mers)

Not HLA-restricted by design. `Passes_Primary_Gate` = "N/A (not HLA-restricted,
by design)" -- never UNRESOLVED, never a fail, consistent with Phase 1F's own
`NOT_APPLICABLE` coverage status for these epitopes.

They ARE long enough for NetMHCIIpan to score directly (confirmed live: a
16-mer scores cleanly at length=16). Since a B-cell epitope that also binds
MHC-II supplies T-cell help and is worth reporting, each is additionally
scanned against the full MHC-II panel and recorded in
`Bcell_MHCII_Rank_Supplementary` / `_Allele` / `Bcell_MHCII_IC50_Supplementary_nM`.
This is purely informational -- it is never included in `Binding_Alleles_Primary`
and never feeds population coverage (Phase4E).

## Regression check vs Phase1Db -- RULING: benign method-scale difference, not an environment fault

21 of the 21 non-B-cell construct epitopes had a stored
`Percentile_Rank` in Phase1Db to compare against (>=20 required). Phase1Db's
own panel (9 MHC-I / 21 MHC-II alleles, via IEDB `method=recommended`)
overlaps almost completely with `ALLELE_FREQ` (8/9 MHC-I, ~30/33 MHC-II).

**Root cause of the spread, confirmed by comparing across the FULL pool
(37 MHC-I + 29 MHC-II candidates, not just the 21 construct epitopes):**
IEDB's `recommended` method resolves to a *different underlying tool per
class*:
- **MHC-I** -- resolves to EL, which is exactly what Phase4B calls directly
  (`netmhcpan_el`). Full pool: mean delta -0.058,
  median 0, range -0.68 to 0.78,
  **zero** epitopes with `|delta|>1` (n=37).
- **MHC-II** -- resolves to BA/consensus, a genuinely **different scale**
  from Phase4B's `netmhciipan_el`. Same number, different meaning -- this is
  not the same measurement re-run, it's two different published methods.
  Full pool: mean delta -0.025, median 0 (no
  systematic bias), but range -3.80 to
  5.20, with
  3 epitopes at `|delta|>1`
  (n=29).

**Among the 21 in-construct MHC epitopes specifically: 2 outliers**
(`HNVWATHACVPTDPN` delta=+5.2, `YKRWIILGLNKIVRM` delta=-1.3), both MHC-II.
**Across the full 29-candidate MHC-II pool: 3 outliers** (the above two plus
`GKLDAWEKIRLRPGG` delta=-3.8, which is pool-only, not a construct member).

**Why this is ruled benign, not escalated:** mean and median deltas are
both ~0 for MHC-II -- there is no systematic bias, only method-scale
variance concentrated in a few peptides. Confirmed **zero gate-status
flips** under Phase 1's own thresholds (<=1.0 MHC-I / <=10.0 MHC-II) across
the entire pool: every epitope keeps the identical pass/fail verdict under
both Phase1Db's `recommended` and Phase4B's `netmhcpan_el`/`netmhciipan_el`.
A deterministic tool disagreeing with a stored result normally means an
environment fault (Section H trigger 4) -- but here the "disagreement" is
between two DIFFERENT, correctly-functioning IEDB methods being compared
against each other, not the same method producing different answers. No
escalation.

**Standing note for Phase 4G (Step 7):** never table Phase 1's MHC-II
`Percentile_Rank` (BA/consensus-derived, via `recommended`) side by side
with Phase 4B's MHC-II rank (`netmhcpan_el`, elution-score) without stating
they are different scales. The same number from each means a different
thing -- presenting them as directly comparable would misrepresent the
Phase I vs Phase IV methodology difference this whole exercise exists to
surface.

## Local demographic secondary view

`Local5_Ranks` / `Local5_Binding_Alleles` report the 5-allele dominant local
profile (A*24:02, B*15:02, B*40:01, DRB1*15:01, DRB1*12:02) as a smaller,
more realistic complement to the broad 107-allele panel. This is drawn from
the SAME netmhcpan_el/netmhciipan_el scan already run against the full panel
(all 5 alleles are members of `ALLELE_FREQ`), so it required no extra API
calls. It does not gate and is reported for context only.
