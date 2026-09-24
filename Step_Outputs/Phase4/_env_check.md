# Phase IV — Step 0 Environment Check

Date: 2026-09-20

## Toolchain fix

Two symlinks created (conda envs and Phenix were moved into `Research (ARCHIVE)/`,
breaking every script's shebang line):

```
/Volumes/Extended SSD/MpoxHIV_deps    -> /Volumes/Extended SSD/Research (ARCHIVE)/MpoxHIV_deps
/Volumes/Extended SSD/phenix-2.2-6143 -> /Volumes/Extended SSD/Research (ARCHIVE)/phenix-2.2-6143
```

`~/mpoxhiv_ssd -> /Volumes/Extended SSD` already existed. Sourced `~/mpoxhiv_env.sh`
after creating the symlinks. Conda was NOT rebuilt.

## The 4 checks

| # | Check | Result |
|---|---|---|
| 1 | `$MPOXHIV_ENVS/phase2/bin/python3 -c "import Bio,numpy,requests"` | **PASS** — Biopython 1.88, numpy 1.26.4, requests imports cleanly |
| 2 | `$BIGMHC_PYTHON external_tools/bigmhc/src/predict.py --help` | **PASS** — full usage text printed (`-i/-m/-o/-s/-d/-v/-b/-c/-a/-p/-t/-j/-f/-z`) |
| 3 | `$BLASTP_BINARY -version` | **PASS** — `blastp: 2.16.0+` (build Mar 28 2025) |
| 4 | `ls human_swissprot_db/human_swissprot.p*` | **PASS** — 8 BLAST DB files present (`.pdb .phr .pin .pjs .pot .psq .ptf .pto`) |

All 4 checks pass. No Phenix/FoldX/APBS/Aggrescan3D/OpenMM/a3d/apbs_env/mhcpredict/sema_env
environments were touched or needed, per the prompt.

## Connectivity spot-check (not one of the 4 required checks, done for confidence)

`curl -4` GET to `https://tools-cluster-interface.iedb.org/tools_api/mhci/` returned
HTTP 405 (Method Not Allowed) — expected for a POST-only endpoint, confirms the host
is reachable over IPv4 without the ~60s IPv6 penalty.

## Rule 1 safety snapshot

Recorded mtimes for all 1,810 files currently under `Phase 1/`, `Phase 2/`,
`Step_Outputs/Phase1/`, `Step_Outputs/Phase2/` to a scratch file, to diff against
at the end of Phase IV and confirm zero modifications.

## Directory scaffolding created

```
Phase 4/_common/, STEP A.. STEP G/
Step_Outputs/Phase4/_run_logs/, StepA.. StepG/
```

No code has been written yet — per Section D, stopping here to report before
proceeding.
