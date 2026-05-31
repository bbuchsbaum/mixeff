# Error-message quality probe: newdata-missing-col

**Scenario:** Call `predict()` with a `newdata` data frame that is missing a column required by the
model formula. Tested on Gaussian LMM (`lmm()` / `lmer()`) with the `sleepstudy` dataset: model
formula is `Reaction ~ Days + (1 | Subject)`, `newdata` contains only `Subject` (omitting `Days`).

**Date:** 2026-05-31  
**mixeff version:** installed (current main branch)  
**lme4 version:** 2.0.1

---

## Verbatim messages

### lme4 — `predict(m_lme4, newdata = newdata_bad)` (re.form = NULL)

```
Error: object 'Days' not found
Class: simpleError, error, condition
```

### lme4 — `predict(m_lme4, newdata = newdata_bad, re.form = NA)` (population-level)

```
Error: object 'Days' not found
Class: simpleError, error, condition
```

### mixeff — `predict(m_mixeff, newdata = newdata_bad)` (re.form = NULL)

```
Error: `newdata` is missing variable(s) required by the model formula: Days.
Class: mm_data_error, mm_condition, rlang_error, error, condition
```

### mixeff — `predict(m_mixeff, newdata = newdata_bad, re.form = NA)` (population-level)

```
Error: `newdata` is missing variable(s) required by the model formula: Days.
Class: mm_data_error, mm_condition, rlang_error, error, condition
```

---

## GLMM paths

- **lme4 `glmer`** (cbpp, newdata missing `period`): `Error: object 'period' not found` —
  same generic R eval error as the LMM path.
- **mixeff `glmm`** with a `cbind()` response: fit itself fails with a clear
  `mm_formula_error` explaining that `cbind(...)` is not in the engine's stateless transform
  subset. This is expected and out-of-scope (PRD §3 non-goals / v2 deferred).
- **mixeff `glmm` predict** is unconditionally refused with `mm_inference_unavailable`:
  "GLMM prediction is not certified by the current Rust contract." This fires before any column
  check, so the missing-column path for GLMM is not exercised (by design; PRD §3).

---

## Assessment

| Criterion | lme4 | mixeff |
|---|---|---|
| Error fires (not silent) | yes | yes |
| Names the missing column | no — only "object 'Days' not found" | **yes** — "missing variable(s) required by the model formula: Days." |
| Points the user at `newdata` | no | **yes** — message explicitly says `` `newdata` is missing... `` |
| Typed condition class | no (plain `simpleError`) | **yes** (`mm_data_error`) — machine-catchable |
| Carries structured metadata | no | **yes** — `missing` field carries the character vector of absent names |
| Both `re.form` paths covered | yes (same generic error) | **yes** (same guard fires first) |

**Verdict: good.**

mixeff's message is materially clearer than lme4's on every axis:

1. lme4 bubbles R's internal evaluation error ("object 'Days' not found") with no indication
   that the source is `newdata`. Users familiar with `predict.lmerMod` eventually learn this
   means a missing column, but the message gives no direct guidance.

2. mixeff fires an R-level guard in `mm_predict_newdata()` (predict.R lines 164-176) *before*
   the Rust FFI is ever called. The message names `newdata`, names the missing variable(s),
   and is caught by the structured `mm_data_error` class with a `missing` attribute—so
   callers can inspect the exact absent names programmatically.

3. The same guard fires identically on both the conditional (`re.form = NULL`) and
   population-level (`re.form = NA`) paths, so the quality is uniform.

**Classification:** `works` — in-scope, already handled correctly. No code change needed.

---

## Source location

`R/predict.R`, function `mm_predict_newdata()`, lines 164–176:

```r
needed <- all.vars(stats::delete.response(stats::terms(fit$formula)))
missing_vars <- setdiff(needed, names(newdata))
if (length(missing_vars)) {
  mm_abort(
    message = sprintf(
      "`newdata` is missing variable(s) required by the model formula: %s.",
      paste(missing_vars, collapse = ", ")
    ),
    class = "mm_data_error",
    input = newdata,
    missing = missing_vars
  )
}
```
