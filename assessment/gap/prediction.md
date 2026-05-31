# Gap Report — Prediction & Residuals

**Family:** Prediction & residuals (`predict`, `fitted`, `residuals`, `simulate`, `bootMer`, simulate-based intervals)
**Date:** 2026-05-31
**Reference:** `assessment/survey/lme4-prediction.md` (lme4 2.0.1 / lmerTest 3.2.1)
**mixeff source:** `R/predict.R`, `R/simulate.R`, `R/methods-extract.R`, `NAMESPACE`
**Standard:** "everything lme4 does, faster, with clearer errors" (PRD §1). Non-goals from PRD §3 cited where relevant.

All statuses below were confirmed by running `library(mixeff); library(lme4)` side-by-side
on `sleepstudy` with `Reaction ~ Days + (Days | Subject)`, not by source-reading alone.

---

## Summary verdict

The core prediction surface that a typical LMM user reaches — `fitted()`, `predict()`
in-sample, `predict(newdata=, re.form=NULL)`, `predict(newdata=, re.form=NA)`,
`predict(allow.new.levels=TRUE/FALSE)`, `residuals(type="response")`,
`simulate(nsim=, seed=, re.form=)` — **works and matches lme4 within the stated
tolerances** (population path to ~1e-12; conditional newdata to ~3e-3; fitted/in-sample
to ~1e-2, consistent with the package's documented theta/fixef tolerances).

However there are **three confirmed silent-surgery bugs** where an unsupported argument
is swallowed by `...` and a *wrong but plausible* answer is returned with no error or
warning. These directly violate the project's "No silent surgery" contract (CLAUDE.md)
and the "clearer errors" mandate, and are the most serious findings in this family:

1. `predict(random.only=TRUE)` silently returns the full conditional prediction.
2. `residuals(scaled=TRUE)` silently returns unscaled residuals.
3. `simulate(newdata=...)` silently ignores `newdata` and simulates the training rows.

Beyond those, the remaining gaps are in-scope-missing (other residual types,
`na.action` propagation, offset-in-newdata, `bootMer`/parametric bootstrap) or
out-of-scope-by-design (conditional-path `se.fit`, `interval=`, `re.form` subset,
GLMM prediction).

---

## Capability table

| lme4 capability | mixeff status | classification | severity | evidence / repro |
|---|---|---|---|---|
| `fitted(m)` | works | works | — | `max abs diff fitted = 0.0114` (within theta-driven tolerance); names attached from `model_frame`. `R/predict.R:101`. |
| `predict(m)` in-sample conditional | works | works | — | `max diff vs lme4 = 0.0114`. Reuses `object$fitted`. `R/predict.R:78-84`. |
| `predict(m, newdata=, re.form=NULL)` conditional | works | works | — | mixeff `312.659,220.243` vs lme4 `312.662,220.244` (~3e-3). Via Rust `predict_new` FFI. `R/predict.R:188-257`. |
| `predict(m, newdata=, re.form=NA / ~0)` population | works | works | — | exact match `282.807,303.742`; `re.form=NA` in-sample diff `1.1e-12`. R-side `X %*% beta`. `R/predict.R:263-300`. |
| `predict(m, allow.new.levels=FALSE)` (default) | works | works | — | Unseen level raises `mm_inference_unavailable` with a clear message naming the level and the remedy. **Clearer than lme4's terse error.** Repro: `predict(m, newdata=data.frame(Days=3, Subject=factor("ZZZ")))`. |
| `predict(m, allow.new.levels=TRUE)` | works | works | — | mixeff `282.807` == lme4 `282.807` (zero-RE / population fallback). Routes to `NewReLevels::population`. `R/predict.R:189`. |
| `predict(m, type="response"/"link")` | works | works | — | Accepted; identical for Gaussian LMM (correct — link == response). `R/predict.R:40,45`. |
| `predict(m, re.form=~(1|grp))` RE subset | partial | out-of-scope-by-design | major | Raises `mm_inference_unavailable` (clean refusal, documented). `R/predict.R:66-76`, `mm_prediction_target()` returns `"unsupported"`. Not in PRD §6 predict arg list (only `re.form`/`NA`/`~0` semantics certified). A real lme4 user plotting one grouping factor's BLUPs cannot do it, but the refusal is structured, not silent. |
| **`predict(m, random.only=TRUE)`** | **partial (SILENT WRONG ANSWER)** | **in-scope-missing** | **blocker** | Arg is swallowed by `...`; **returns full conditional prediction**, not the RE-only contribution. Confirmed: `predict(m, random.only=TRUE)` head `253.65,273.32,292.99` == `predict(m)`; lme4 RE-only head is `2.26,11.46,20.66`. No error, no warning. Violates "No silent surgery." Not in PRD §6 list, but silently returning a wrong answer is a defect regardless of scope. |
| `predict(m, newparams=)` | in-scope-missing | in-scope-missing | minor | Arg absent (caught by `...`, silently ignored — counterfactual prediction not honored). Niche; lme4 itself routes through `setParams`. Not in PRD §6 list. |
| `predict(m, terms=)` | works (parity) | out-of-scope-by-design | cosmetic | lme4 2.0.1 itself throws "not yet implemented"; mixeff ignores it. No user-visible gap. |
| **`predict(m, se.fit=TRUE)` — population/`re.form=NA` path** | partial | partial | major | Returns `$se.fit` all `NA` with `mm_unavailable_reason="prediction_se_unavailable_phase_2"`. **But `vcov(m)` returns the 2x2 fixed-effect covariance** — the fixed-only/population SE (`sqrt(diag(X_new %*% vcov %*% t(X_new)))`) is computable from materials already on the R side and is exact for LMM. Returning NA here is more conservative than necessary. PRD §6 lists `se.fit=` as in-scope API. |
| `predict(m, se.fit=TRUE)` — conditional path | partial | out-of-scope-by-design | major | Same NA stub. The conditional SE needs the joint FE+RE/prediction-variance from Rust, which is not certified. PRD line 380-382 explicitly designs this as `$se.fit=NA` until Rust certifies. Out-of-scope-by-design for the conditional path. |
| `predict(m, na.action=na.pass)` NA rows in newdata | in-scope-missing | in-scope-missing | major | mixeff **hard-errors** on an NA row: `mm_data_error: numeric column 'Days' contains a non-finite value (NaN)`. lme4 returns `312.66, NA, 220.24` (NA-in NA-out, row count preserved). `na.action` arg absent from `predict.mm_lmm`. A common cross-validation/grid workflow breaks. Error is clear, but the capability is missing. |
| Offset in newdata (formula `offset()` terms) | in-scope-missing | in-scope-missing | minor | `mm_predict_fixed_only()` (`R/predict.R:263-300`) builds the FE design from `mm_fixed_formula` and does not evaluate `offset()` terms against newdata. PRD §3 line 42 lists *in-fit* offsets as a non-goal, but predict-time offset evaluation for an already-fit model is a distinct surface. Confirm-on-demand; not exercised here. |
| `predict(m, interval="confidence"/"prediction")` | partial | out-of-scope-by-design | minor | Raises `mm_inference_unavailable` (depends on SE). PRD line 382-383 explicitly: `interval != "none"` without SE raises rather than fabricating. Designed refusal. |
| GLMM `predict()` | in-scope-missing | out-of-scope-by-design | major | `predict.mm_glmm` always raises `mm_inference_unavailable` ("not certified by the current Rust contract"). `R/predict.R:126-132`. GLMM prediction SE/profile paths are PRD §3 v2 deferrals; the point estimate refusal is a structured diagnostic, not silent. Crate survey (`crate-8.md:56`) marks GLMM `predict(newdata=)` as surfaced-but-not-delivered. |
| `residuals(m, type="response")` | works | works | — | Default; `max diff vs lme4 = 0.0114`. `R/predict.R:109-114`. |
| `residuals(m, type="pearson")` | in-scope-missing | in-scope-missing | major | `match.arg` only allows `"response"`; `type="pearson"` throws `simpleError` from `match.arg`. For LMM this is just `response/sigma` and is trivially computable. PRD §6 line 229 lists `residuals(type=)`. |
| `residuals(m, type="deviance")` | in-scope-missing | in-scope-missing | major | Same `match.arg` rejection. This is the **lme4 default for GLMM**, so GLMM users expect it. (LMM deviance residual == response residual.) |
| `residuals(m, type="working")` | in-scope-missing | in-scope-missing | minor | Same rejection. For LMM equals response residuals; mainly matters for GLMM. |
| `residuals(m, type="partial")` | works (parity) | out-of-scope-by-design | cosmetic | lme4 2.0.1 itself throws "not implemented yet"; mixeff rejects via `match.arg`. No real-world gap. |
| **`residuals(m, scaled=TRUE)`** | **partial (SILENT WRONG ANSWER)** | **in-scope-missing** | **major** | `scaled` arg swallowed by `...`; **returns unscaled residuals**. Confirmed: `residuals(m, scaled=TRUE)` == `residuals(m)`; lme4 scaled head `-0.16, -0.571, -1.649`. No error/warning. Trivially `residuals/sigma`. Violates "No silent surgery." |
| `residuals` NA re-insertion via `naresid` | in-scope-missing | in-scope-missing | minor | No `na.action` tracking; residual vector length equals fitted-row count only. Edge case unless data had NA at fit time. |
| GLMM `residuals()` | partial | partial | major | `residuals.mm_glmm <- residuals.mm_lmm` (`R/predict.R:122`) — inherits the single-type "response" limitation; GLMM users routinely need `deviance`/`pearson`. |
| `simulate(m, nsim=, seed=)` | works | works | — | `simulate(m, nsim=2, seed=42)` returns 180x2 `data.frame` with `sim_1,sim_2`; seed restored via `mm_with_seed`. `R/simulate.R:57-84`. |
| `simulate(m, re.form=NULL)` new RE draws | works | works | — | `mm_simulate_random_mean` draws RE from reconstructed covariance. `R/simulate.R:95-119`. |
| `simulate(m, re.form=NA)` marginal | works | works | — | Uses `fixed_fitted` + residual noise. `R/simulate.R:87-92`. |
| **`simulate(m, newdata=...)`** | **partial (SILENT WRONG ANSWER)** | **in-scope-missing** | **major** | `newdata` swallowed by `...`; **simulates the training rows, ignoring newdata**. Confirmed: `simulate(m, newdata=d[1:5,])` returns 180 rows; lme4 returns 5. No error/warning. Violates "No silent surgery." |
| `simulate(m, use.u=TRUE)` | in-scope-missing | in-scope-missing | minor | Arg absent (silently ignored). Deprecated in lme4 in favour of `re.form`, so low priority — but should error if passed, not be swallowed. |
| `simulate(m, newparams=)` | in-scope-missing | in-scope-missing | minor | Absent / silently ignored. |
| `simulate(m, re.form=~(...))` subset | partial | out-of-scope-by-design | minor | Raises `mm_inference_unavailable` (clean). `R/simulate.R:66-73`. Same RE-subset deferral as predict. |
| `simulate(m, cond.sim=, cluster.rand=, weights=, offset=, family=, formula=)` | in-scope-missing | out-of-scope-by-design | minor | Not in signature; advanced/experimental simulate knobs. Most silently ignored via `...`. |
| GLMM `simulate()` | in-scope-missing | out-of-scope-by-design | major | No `simulate.mm_glmm` method (only `S3method(simulate, mm_lmm)` in NAMESPACE). GLMM simulate is part of PRD §10 Phase 4 but not yet present for GLMM objects; dispatches to default `stats::simulate` (likely errors or wrong). |
| `bootMer(m, FUN, type="parametric")` | in-scope-missing | in-scope-missing | major | No `bootMer` in mixeff namespace (`exists("bootMer", asNamespace("mixeff")) == FALSE`). PRD §10 Phase 4 (line 807) explicitly lists "parametric bootstrap" as in-scope. `refit()` + `simulate()` building blocks exist (`R/simulate.R:13-38`), so a `bootMer`-equivalent is feasible but unwritten. |
| `bootMer(type="semiparametric")` | in-scope-missing | out-of-scope-by-design | minor | lme4 calls this experimental/LMM-only; not in PRD. |
| Simulate-based prediction intervals | in-scope-missing | in-scope-missing | major | No CI infrastructure; `interval=` refuses and no `bootMer`. The standard lme4 PI workflow (`bootMer(FUN=predict)` + quantiles) is unavailable end-to-end. |
| `bootMer(parallel="multicore")` | in-scope-missing | out-of-scope-by-design | cosmetic | Depends on bootMer; not applicable until that exists. |
| `ranef(m)` BLUPs | works | works | — | (Covered in ranef family; used here in prediction workflows.) `ranef.mm_lmm`. |
| `ranef(m, condVar=TRUE)` | works | works | — | Implemented via Rust `cond_var()`; parity within 1e-3. (ranef family.) |
| `refit(m, newresp)` | works | works | — | Refits via `lmm()` on stored model frame; tags provenance. `R/simulate.R:13-38`. Stricter than lme4: rejects NA/length-mismatch with `mm_arg_error`. |

---

## Severity rationale

- **Blocker:** `random.only=TRUE` silent wrong answer. A documented lme4 argument returns a
  numerically wrong result with no diagnostic — the exact failure mode CLAUDE.md forbids.
- **Major (silent wrong answers):** `residuals(scaled=TRUE)` and `simulate(newdata=)` —
  same class of bug, slightly lower blast radius (scaling is a constant factor; simulate
  newdata is a less common path), but each silently betrays the audit-first contract.
- **Major (clean but missing):** `na.action`/NA-in-newdata, `pearson`/`deviance` residuals,
  GLMM residual types, `bootMer`/parametric bootstrap + prediction intervals, population-path
  `se.fit` (computable from `vcov`), GLMM predict/simulate. These are capabilities a real
  lme4 user hits; mixeff either errors cleanly or refuses.
- **Minor/cosmetic:** offset-in-newdata, `newparams`, `use.u`, experimental simulate knobs,
  `terms=`/`type="partial"` (parity with lme4's own non-implementation).

## Highest-leverage fixes

1. **Make swallowed args loud.** `predict.mm_lmm`, `residuals.mm_lmm`, `simulate.mm_lmm`
   should detect unsupported-but-recognized lme4 args (`random.only`, `scaled`, `newdata`
   for simulate, `newparams`, `use.u`, `na.action`) in `...` and raise a structured
   `mm_inference_unavailable`/`mm_arg_error` instead of silently dropping them. Removes all
   three silent-wrong-answer defects at once.
2. **Implement the cheap residual/scaled paths in R:** `pearson` = response/√Var, `scaled`
   = response/sigma, `deviance` (== response for LMM). No Rust work needed for LMM.
3. **Population-path `se.fit`** from the existing `vcov(m)` via the delta method (exact for LMM).
4. **`na.action` propagation** for NA rows in newdata (na.pass → NA prediction, row count preserved).
