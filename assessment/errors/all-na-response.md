# Error-message quality: response variable is entirely NA

**Scenario:** `all-na-response`
**Date probed:** 2026-05-31
**mixeff version:** installed from repo (main branch)
**lme4 version:** 2.0.1 / lmerTest 3.2.1 (R 4.5.2)

---

## Setup

```r
set.seed(42)
n <- 40
df <- data.frame(
  y       = rep(NA_real_, n),   # response is entirely NA
  x       = rnorm(n),
  subject = factor(rep(seq_len(10), each = 4))
)
```

---

## Verbatim messages

### lme4::lmer (gaussian, y = all NA)

```
Class:  simpleError, error, condition
Message:
0 (non-NA) cases
```

### lme4::glmer (gaussian/identity shortcut → delegates to lmer)

```
Class:  simpleError, error, condition
Message:
0 (non-NA) cases
```

### lme4::glmer (binomial or poisson, y = all NA)

With non-gaussian families, lme4 first strips NA rows from the model frame.
All 40 rows are dropped, leaving an empty data frame. Instead of reporting the
cause (all-NA response), it then fails during grouping-factor setup:

```
Class:  simpleError, error, condition
Message:
Invalid grouping factor specification, subject
```

This is a secondary, misleading error that does not mention NA values at all.

### mixeff::lmm (y = all NA)

```
Class:  mm_data_error, mm_condition, rlang_error, error, condition
Message:
Missing values in design variable(s): `y` (40 NA). mixeff requires complete
cases; pass na.omit(data) explicitly before fitting.
```

### mixeff::glmm (binomial, y = all NA)

```
Class:  mm_data_error, mm_condition, rlang_error, error, condition
Message:
Missing values in design variable(s): `y` (40 NA). mixeff requires complete
cases; pass na.omit(data) explicitly before fitting.
```

### mixeff::glmm (gaussian/identity)

Not directly comparable: mixeff correctly refuses `gaussian/identity` as outside
its certified GLMM contract (`mm_inference_unavailable`). The lmm() path is the
right one for gaussian data and is tested above.

---

## Analysis

| Axis | lme4::lmer | lme4::glmer (non-gauss) | mixeff::lmm / glmm |
|---|---|---|---|
| Error type | `simpleError` (untyped) | `simpleError` (untyped) | `mm_data_error` (typed) |
| Message text | `0 (non-NA) cases` | `Invalid grouping factor specification, subject` | `Missing values in design variable(s): \`y\` (40 NA). mixeff requires complete cases; pass na.omit(data) explicitly before fitting.` |
| Names the affected variable? | No | No | Yes — `` `y` `` |
| Reports NA count? | No | No | Yes — `(40 NA)` |
| Actionable remedy? | No | No | Yes — `pass na.omit(data) explicitly` |
| Typed / machine-catchable? | No | No | Yes — `mm_data_error` |
| Misleading secondary cause? | No | Yes — "Invalid grouping factor" masks the real cause | No |

**lme4::lmer** produces `0 (non-NA) cases`, which is marginally more accurate
than nothing, but:
- does not name which variable is all-NA,
- does not count how many NAs,
- gives no remedy, and
- is untyped so cannot be programmatically caught by class.

**lme4::glmer** (with non-gaussian families) is actively misleading: because it
silently drops all NA rows before checking grouping factors, the user sees
`Invalid grouping factor specification, subject` — a confusing secondary failure
that names the grouping variable rather than the all-NA response. A user
encountering this would likely investigate the `subject` column, not `y`.

**mixeff** (both `lmm` and `glmm`) catches the all-NA condition during the
pre-fit `mm_check_no_na()` validation in `compile_model()`, before any Rust
call. The error:
- names the specific column (`y`),
- reports the exact NA count (40 NA),
- provides a concrete remedy (`pass na.omit(data) explicitly before fitting`),
- is typed as `mm_data_error` (machine-catchable),
- does not fabricate or guess — it reports exactly what it measured.

---

## Verdict

**clearer-than-lme4** — mixeff's `mm_data_error` message is materially clearer
than lme4 in all compared variants. lme4's best case (`0 (non-NA) cases`) is
terse and untyped; its non-gaussian glmer case is actively misleading.

**Overall verdict: good** — no defect; this scenario is handled well and is a
concrete illustration of the "less inscrutable errors" promise. The pre-fit NA
check fires early (before any FFI boundary), names the variable and count, and
supplies a concrete remedy.

---

## Possible improvement (not a bug)

The message could distinguish between the *response* being all-NA versus a
*predictor* being all-NA, e.g.:

> "Response variable `y` is entirely NA (40/40 missing). No observations
> remain for fitting. Pass `na.omit(data)` explicitly before fitting."

This would make the message even more specific for the all-NA-response case
versus a partially-NA predictor (which hits the same code path with a different
count). This is a polish item, not a defect.
