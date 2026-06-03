# OSF "Willingness to wait" — in-the-wild GLMM parity reproduction

An "in the wild" check that `mixeff` reproduces a *published* `lme4::glmer`
analysis with **crossed random effects and correlated random slopes**. Tracked
by mote **bd-01KT3ZRCKWRZQFA4W7TXTGWAZ0**.

## Source

OSF node **[ftexh](https://osf.io/ftexh/)**, *"Willingness to wait study"*:

- Analysis script: <https://osf.io/ftexh/files/pv6u2> → `Analysis- Correlations, LME4.R`
- Study 1a trial data: <https://osf.io/download/u2zt4/> → `FullData_study1a.csv` (1427 rows)
- Study 1b (replication) trial data: <https://osf.io/download/3bxkt/> → `FullData_replication.csv` (1600 rows)

The script fits **9 binomial-logit GLMMs**. Two families:

| family | random effects | example formula |
| ------ | -------------- | --------------- |
| wait decisions (×4) | correlated slopes on **both** crossed groups | `wait_choice ~ Enjoyment_centered (+ arousal_centered) + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)` |
| comprehension (×5) | random intercepts / one slope | `score ~ Enjoyment_centered (+ ...) + (1 | ID) + (1 | Title)` |

`ID` is a participant id; `Title` is the (text) book title. In the raw data
**`ID` is stored as an integer and `Title` as a character string** — neither is
a factor. The comprehension `score` is `Q1_correct`/`Q2_correct` gathered long.

## Why this case matters

It exercises three surfaces the curated fixtures under-covered:

1. **Grouping coercion (fixed here).** lme4 silently `factor()`s its grouping
   variables; before this case `mixeff` refused an integer/numeric grouping
   column with `"grouping factor 'ID' not found or not categorical"`. Every one
   of the 9 models failed at construction. `lmm()`/`glmm()` now coerce
   non-categorical grouping columns to factors (announced via a suppressible
   notice, never silently) — see `R/fit-lmm.R::mm_coerce_grouping_factors()`.
2. **Correlated random slopes on crossed groups.** `joint_laplace` matches
   `glmer` to **~1e-4 on fixed effects and ~1e-6 on logLik** here (and the
   random-effect SDs to ~3 decimals). Strong evidence the engine handles the
   `(1 + x | g1) + (1 + x | g2)` surface correctly.
3. **High-baseline random-intercept shortfall (documented gap).** The
   comprehension models (high accuracy → intercept ≈ +1.3) reveal a *premature*
   convergence: `joint_laplace` (and `pirls_profiled`) report
   `converged_interior` after ~20 iterations at a point **0.01–0.05 logLik below
   the glmer MLE**, with the intercept biased toward 0 by ~0.02–0.04. More
   optimizer budget does not help. Upstream
   **bd-01KT40T6FGVXQQ9N50G2HM0ZZE** (distinct from Gap E
   bd-01KT3Z64AY45NHA5144G2ZBMSY, which is about *badly-scaled* designs).

## Findings (2026-06-02, mixeff-rs pin d2145163)

- **Estimation parity holds for the headline (correlated-slope) models.**
  `wait_*` models match `glmer` to `max|dFixef| < 5e-3`, `|dlogLik| < 5e-2`.
- **Random-intercept comp models land within ~0.04 of `glmer`** — still
  scientifically usable, conclusions preserved — but short of certified
  parity; tracked upstream as above.
- **GLMM Wald inference is uncertified** (upstream
  bd-01KT3Z64YE5QN7626PQRJSJJVA): `summary()` reports `method = "not_computed"`
  and withholds z / p. The published z/p values are therefore **not** yet
  reproducible through `mixeff`; point estimates and logLik are.

## Files

| file | purpose |
| ---- | ------- |
| `reconstruct.R`     | download the OSF CSVs and write the slimmed modeling fixtures `tests/fixtures/osf_willingness_to_wait_study1{a,b}.csv`. Needs network. |
| `parity-harness.R`  | fit each model with `glmer` and `mixeff::glmm` (`joint_laplace`) and print estimate / SE / logLik / RE-variance comparisons. Runs offline on the committed fixtures. |
| `../../tests/fixtures/osf_willingness_to_wait_study1a.csv` | committed Study-1a fixture (1427 rows, 8 modeling columns). |
| `../../tests/fixtures/osf_willingness_to_wait_study1b.csv` | committed Study-1b fixture (1600 rows). |
| `../../tests/testthat/test-glmm-osf-willingness-to-wait-parity.R` | gated regression test (full 9-model sweep behind `MIXEFF_RUN_SLOW_PARITY`). |
