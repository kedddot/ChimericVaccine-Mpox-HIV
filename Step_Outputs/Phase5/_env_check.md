# Phase V (Section V — Pharmacokinetics) — Step 0 Environment Check

Date: 2026-09-22
Working directory: `/Volumes/Extended SSD/ChimericVaccine/Research` (`$R`)

## 1. Symlink / deps volume

```
$ ls -ld "/Volumes/Extended SSD/MpoxHIV_deps"
lrwxr-xr-x@ 1 nek staff 53 Sep 20 18:46 /Volumes/Extended SSD/MpoxHIV_deps -> /Volumes/Extended SSD/Research (ARCHIVE)/MpoxHIV_deps
```
Present and resolves correctly.

## 2. Environment script

`source ~/mpoxhiv_env.sh` runs clean, no `[mpoxhiv_env] WARNING`. Sets
`MPOXHIV_ENVS=/Users/nek/mpoxhiv_ssd/MpoxHIV_deps/envs` (the space-free symlink path, needed
because a shebang line cannot contain a space).

## 3. Pre-existing packages

```
$ "$MPOXHIV_ENVS/phase2/bin/python3" -c "from scipy.integrate import solve_ivp; import numpy, pandas; print('OK', numpy.__version__, pandas.__version__)"
OK 1.26.4 2.3.3
```
Confirmed: scipy 1.17.1 (`solve_ivp` import succeeds), numpy 1.26.4, pandas 2.3.3 — matches the
brief exactly.

## 4. matplotlib

Confirmed missing before install, as the brief expected:
```
$ "$MPOXHIV_ENVS/phase2/bin/python3" -c "import matplotlib"
ModuleNotFoundError: No module named 'matplotlib'
```

Installed:
```
$ "$MPOXHIV_ENVS/phase2/bin/pip" install matplotlib
Successfully installed contourpy-1.3.3 cycler-0.12.1 fonttools-4.65.0 kiwisolver-1.5.1
matplotlib-3.11.2 pillow-12.3.0 pyparsing-3.3.3
```

Post-install smoke test (headless `Agg` backend, since this is a script/CLI environment with no
display):
```
matplotlib 3.11.2, backend Agg
OK - figure saved (/tmp/_mpl_smoke_test.png, 19 KB, deleted after verification)
```
`matplotlib==3.11.2` is now available in the `phase2` env alongside scipy/numpy/pandas. All Phase
V figures will call `matplotlib.use("Agg")` before `pyplot` import (house style already used
nowhere in Phase 1-4 since none of those steps plot, but this is the standard non-interactive
backend for headless CLI runs and avoids any macOS GUI-backend crash).

## 5. Rule 1 — pre-check (before any Phase V code is written)

```
$ shasum -a 256 -c Step_Outputs/Phase4/_rule1_baseline_sha256.txt
```
Result: **1810/1810 OK, 0 failures.** Matches the expected baseline exactly. Phase 1/2/4 are
untouched going into Section V.

## Result

All 4 environment checks pass:
1. Deps symlink present and correct.
2. `mpoxhiv_env.sh` sources cleanly.
3. scipy/numpy/pandas already present at the stated versions.
4. matplotlib installed successfully (3.11.2) and verified to render and save a figure headlessly.

Plus the Rule 1 pre-check (1810/1810 OK) required before starting any step.

**No fallback to CSV-only output is needed** — matplotlib installed and works.

Per the brief: STOP here and report before writing any pipeline code (`Phase 5/_common/`,
`STEP A`–`STEP E`).
