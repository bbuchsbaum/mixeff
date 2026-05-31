# Error-message quality probe: empty-data

**Scenario:** zero-row data frame (correct columns, correct types, 0 observations).

**Probe script:** `assessment/errors/probe-empty-data.R`

---

## Verbatim messages

### lme4::lmer (LMM)

```
type   : error
class  : simpleError, error, condition
message: 0 (non-NA) cases
```

### mixeff::lmm (LMM)

```
type   : error
class  : mm_bridge_error, mm_condition, rlang_error, error, condition
message: Matrix index out of bounds.
Caused by error in `doTryCatch()`:
! Matrix index out of bounds.
```

### lme4::glmer (GLMM, binomial)

```
type   : error
class  : simpleError, error, condition
message: Invalid grouping factor specification, subject
```

### mixeff::glmm (GLMM, binomial)

```
type   : success   <-- NO ERROR RAISED
class  : NA
message: NA

Post-fit inspection:
  class    : mm_glmm, mm_fit, mm_compiled
  fixef    : (Intercept) = 0, x = 0
  nobs     : 0
  logLik   : 0
  fitted() : length 0
```

---

## Assessment

### mixeff::lmm — verdict: **needs-work** (worse than lme4)

lme4 says `"0 (non-NA) cases"` — terse but immediately actionable: the caller
knows the data has no rows.

mixeff raises an `mm_bridge_error` with `"Matrix index out of bounds."`. This
is a raw Rust-engine panic that leaks implementation internals. It does not:

- name the actual problem (zero observations),
- tell the caller what to fix,
- use the typed `mm_data_error` class that `R/conditions.R` documents for
  "data shape / type problems".

The condition reaches the R surface as `mm_bridge_error` (the untagged fallback
class), which means the Rust engine never emitted a structured `mm_data_error:`
tag. The R wrapper (`fit-lmm.R`) has no pre-flight `nrow(data) == 0` guard
before delegating to Rust, so the panic propagates unfiltered.

**Classification:** `in-scope-missing` — the PRD audit-first contract requires
structured diagnostics for detectable data-shape problems; empty data is the
simplest such case and should be caught R-side with a typed `mm_data_error`
before the Rust call is ever made.

### mixeff::glmm — verdict: **bug** (silent wrong answer)

mixeff::glmm succeeds silently and returns a structurally valid `mm_glmm`
object with:
- all fixed-effect coefficients = 0 (fabricated),
- logLik = 0 (fabricated),
- nobs = 0.

This is a direct violation of the package's audit-first design principle:
*"every model reduction or refusal crosses the boundary as a structured
diagnostic … code that hides a transformation from the user is a bug, not a
convenience"* (CLAUDE.md). Returning a zero-coefficient fit on zero-row data
fabricates inference without any diagnostic.

lme4::glmer errors with `"Invalid grouping factor specification, subject"` —
also terse, but at least it refuses rather than fabricating.

**Classification:** `bug` — `glmm()` must detect `nrow(data) == 0` R-side
(after `compile_model`) and raise a typed `mm_data_error` before any Rust call.
The Rust engine evidently accepts empty input and returns a degenerate zero
solution, which the R wrapper surfaces as a success.

---

## Recommended fix

Add a guard in both `lmm()` and `glmm()` immediately after `compile_model()`
(where the model frame is available), before `mm_translate_data()`:

```r
if (nrow(spec$model_frame) == 0L) {
  mm_abort(
    message = paste0(
      "Data has 0 rows after applying the formula and any NA removal. ",
      "A model cannot be fit without observations."
    ),
    class   = "mm_data_error",
    input   = data
  )
}
```

This would:
1. Use the correct typed condition class (`mm_data_error`).
2. Name the exact problem ("0 rows") and what the caller must fix.
3. Prevent the Rust engine from ever seeing empty input (eliminating both the
   matrix-index panic in lmm and the silent zero-fit in glmm).
4. Produce a message that is materially clearer than both lme4 messages.
