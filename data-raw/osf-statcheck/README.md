# OSF statcheck — in-the-wild GLMM parity reproduction

An "in the wild" check that `mixeff` reproduces a *published* `lme4::glmer`
analysis. Tracked by mote **bd-01KT3ZB649HHX4ZYMKTQ5A18XX**.

## Source

OSF node **538bc**, *"Journal Data Sharing Policies and Statistical Reporting
Inconsistencies in Psychology"*, **Study 3** (statcheck × Open-Practice badges):

- Analysis script: <https://osf.io/jz9r6/> → `171020MultilevelAnalysis.R`
- Data component: <https://osf.io/st2ex/> → `170329MergedDataStatcheckBadges.txt`

All 37 active models in the published script are binomial-logit GLMMs of the
form `Error`/`DecisionError ~ predictors + (1 | Source)` — a pure `glmer`-parity
surface.

## Data caveat (important)

The OSF-hosted merged file is a **lossy UTF-16 re-save** of the authors'
`write.csv` output: rows whose article title contains a comma lost their
field-level quoting (and the title was truncated). `reconstruct.R` recovers the
modeling columns with a **right-anchored parser** (the trailing 19 columns are
intact; the title truncation is deterministic, so the `(1 | Source)` grouping is
preserved). Fidelity is verified: `glmer(nAGQ=0)` on the reconstruction
reproduces the paper's reported `OpenPractice:Year` = 0.7958, Z = 1.825,
p = .0679 to 4 decimals.

## Files

| file | purpose |
| ---- | ------- |
| `reconstruct.R`     | download the OSF `.txt`, robustly parse it, write `inst/extdata/osf-statcheck-t2.csv` (Period-2 slice, anonymized `gid`). Needs network. |
| `parity-harness.R`  | fit each model with `glmer` (nAGQ=1) and `mixeff::glmm` (`joint_laplace`); print estimate/SE/logLik comparison. Runs offline on the committed fixture. |
| `../../inst/extdata/osf-statcheck-t2.csv` | committed Period-2 fixture (5279 rows, 426 groups, 380 `Error` events). |
| `../../tests/testthat/test-glmm-osf-statcheck-parity.R` | gated regression test. |

## Findings (2026-06-02, mixeff-rs pin d2145163)

- **Estimation parity holds.** `joint_laplace` ≈ `glmer` on well-conditioned
  (centered-`Year`) models: max |fixef diff| ≈ 3.4e-3, max |logLik diff| ≈ 1e-3.
- **Gap E (upstream bd-01KT3Z64AY45NHA5144G2ZBMSY).** On raw `Year` (large
  offset), `joint_laplace` returns `converged_interior` with a non-finite
  objective and a sub-optimal coefficient (interaction 0.7959 vs the
  offset-invariant MLE 0.8529/0.8544). `glmer` fails loud on the same design.
- **Gap B (upstream bd-01KT3Z64YE5QN7626PQRJSJJVA).** GLMM standard errors are
  uncertified (active-subspace Hessian unavailable); `summary()` gives
  `method=not_computed`, z/p = NA. Recorded as `unsupported` in the parity
  ledger.
- **Use `method="joint_laplace"`** for glmer parity — the default
  `pirls_profiled` is a distinct (profiled) approximation (already covered by
  `test-glmm-joint-laplace-parity.R`).
