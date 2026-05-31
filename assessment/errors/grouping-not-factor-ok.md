# Error-Message Quality Assessment: grouping-not-factor-ok

**Scenario:** Grouping variable in `(1 | subject)` is numeric/integer (not a factor or character).  
**Date:** 2026-05-31  
**mixeff version:** installed (main branch, post d030be6)  
**lme4 version:** 2.0.1

---

## Setup

```r
set.seed(42)
n_subjects <- 30; n_obs_per <- 5; N <- 150
df <- data.frame(
  y       = rnorm(N),
  x       = rnorm(N),
  subject = rep(seq_len(n_subjects), each = n_obs_per)  # integer, NOT factor
)
# class(df$subject) == "integer", is.factor(df$subject) == FALSE
```

---

## lme4 behavior

### `lmer()` with integer grouping (30 levels, 5 obs each)

**Result:** Fits silently, no warning, no error.

lme4 implicitly coerces the integer column to a factor at fit time (`class(fit@flist$subject) == "factor"`). The user gets a working model with no indication that a coercion happened.

### `lmer()` with continuous double grouping (150 unique values)

**Exact error:**
```
Error in lmer(y ~ x + (1 | subject), data = df_cont, REML = FALSE) :
  number of levels of each grouping factor must be < number of observations (problems: subject)
```

This error fires because levels == observations (a structural impossibility), not because the column type is wrong. The message names the problem (levels >= observations) but not the root cause (non-factor column).

### `glmer()` with integer grouping (30 levels)

**Result:** Fits silently, no warning, no error (same implicit coercion as `lmer`).

---

## mixeff behavior

### `lmm()` with integer grouping (30 levels)

**Error class:** `mm_fit_error, mm_condition, rlang_error, error, condition`

**Exact error message:**
```
failed to construct LMM: Invalid argument: Grouping factor 'subject' not found or not categorical
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct LMM: Invalid argument: Grouping factor 'subject' not found or not categorical
```

**Error fields:** `message`, `trace`, `parent`, `formula`, `rlang`, `call`  
- `err$formula`: `"y ~ x + (1 | subject)"`  
- No `hint`, no `fix`, no `column_type`, no actionable suggestion field.

### `glmm()` with integer grouping (binomial)

**Exact error message:**
```
failed to construct GLMM: Invalid argument: Grouping factor 'subject' not found or not categorical
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: Invalid argument: Grouping factor 'subject' not found or not categorical
```

Same class, same structure, same absence of actionable guidance.

### `lmm()` with continuous double grouping (150 unique values)

**Exact error message:**
```
failed to construct LMM: Invalid argument: Grouping factor 'subject' not found or not categorical
```

Identical message — the error gives no information distinguishing "you used an integer" from "you used a float" from "the column name is misspelled".

---

## Comparative assessment

| Dimension | lme4 | mixeff |
|---|---|---|
| Crashes / panics | No | No |
| Correct typed error class | N/A (base R condition) | Yes (`mm_fit_error`) |
| Names the offending variable | No | Yes (`subject`) |
| Explains the root cause | No (silently coerces integer) | Partially ("not categorical") |
| Actionable fix suggestion | No | No |
| Distinguishes integer vs float vs misspelling | No | No |
| Formula field on error object | No | Yes |

### What mixeff does better than lme4

- It **refuses** to silently coerce a numeric grouping variable to a factor, which is the audit-first design intent. lme4 silently fits with a coerced factor — a "silent surgery" that mixeff explicitly rejects.
- The error is **typed** (`mm_fit_error`) and carries the `formula` field, enabling programmatic handling.
- The error fires at **compile time** (before any numerical work), not mid-optimization.

### What mixeff does worse / needs improvement

The error message `"Grouping factor 'subject' not found or not categorical"` is **ambiguous and non-actionable**:

1. **"not found OR not categorical"** conflates two completely different problems — a misspelled variable name vs. a wrong column type. A user who spelled the variable name correctly will be confused by "not found."

2. **No fix suggestion.** The message does not say "use `factor(subject)` in your formula or coerce the column before fitting." A user seeing this for the first time has no path forward without reading documentation.

3. **No column-type diagnostic.** The message does not report the actual type observed (`integer`, `numeric`, `character`) vs. what was expected (`factor` or `character`). A message like:  
   `"Grouping factor 'subject' must be a factor or character column, but it is integer. Use factor(subject) in the formula or convert the column with as.factor() before fitting."`  
   would fully explain the problem and provide an actionable fix.

4. **The "Caused by" duplication** repeats the same sentence twice, adding noise without adding information.

5. **The R-side has enough information to provide a better message.** At the point where `compile_model()` / `lmm()` calls `mm_translate_data()`, it already knows: (a) the formula specifies `subject` as a grouping variable, and (b) `df$subject` is `integer`. The R wrapper could intercept this and emit a specific diagnostic before even calling into Rust.

---

## Verdict

**needs-work**

The error is not a crash, not a Rust panic, not a silent wrong answer, and it is typed correctly. It is materially better than lme4 in that it refuses silent coercion. However, the message conflates "not found" with "not categorical", provides no actionable fix suggestion, and omits the observed column type. For a package whose core promise is "clearer errors than lme4," this is an improvement opportunity: the R wrapper has all the information needed to emit a precise, self-correcting diagnostic before delegating to Rust.

### Recommended fix (R-side, in `compile_model()` or `lmm()`)

After extracting `vars` and before calling `mm_translate_data()`, inspect each grouping variable mentioned in the formula's random-effect terms. If any is present in `data` but is numeric/integer (not factor or character), raise an `mm_data_error` with a message such as:

```
Grouping variable `subject` must be a factor or character, but is integer.
Coerce it before fitting: lmm(y ~ x + (1 | factor(subject)), data = df)
or: df$subject <- factor(df$subject)
```

This keeps the "no silent surgery" contract while giving the user a precise, actionable error at the R layer rather than an ambiguous bridge error from Rust.

---

## Classification

**in-scope-missing** — the typed error fires correctly and is not a crash, but the message content needs improvement to meet the "clearer than lme4" bar. The R wrapper has sufficient information to emit an actionable diagnostic; the Rust-side message leaks through without enrichment.
