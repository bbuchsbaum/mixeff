# Error-message quality probe: "one-row" scenario

**Date:** 2026-05-31  
**Scenario:** single-row (n=1) data frame passed to a mixed model fitting function.  
**Packages:** lme4 2.0.1 / lmerTest 3.2.1 vs mixeff (current main).

---

## Exact messages captured

### Case 1 — truly one row, lme4

```
lme4::lmer(y ~ x + (1 | subject), data = one_row)
```
```
[ERROR] grouping factors must have > 1 sampled level
```

### Case 1 — truly one row, mixeff::lmm

```
mixeff::lmm(y ~ x + (1 | subject), data = one_row,
            control = mm_control(verbose = -1))
```
```
[ERROR class=mm_fit_error, mm_condition, rlang_error, error, condition]
failed to fit LMM: Constant response: model fitting failed
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to fit LMM: Constant response: model fitting failed
```

### Case 2 — 2 rows, 1 group level, varied y (isolates grouping vs constant-response), lme4

```
lme4::lmer(y ~ x + (1 | subject), data = two_rows_one_group)
```
```
[ERROR] grouping factors must have > 1 sampled level
```

### Case 2 — 2 rows, 1 group level, varied y, mixeff::lmm

```
mixeff::lmm(y ~ x + (1 | subject), data = two_rows_one_group,
            control = mm_control(verbose = -1))
```
```
[ERROR class=mm_fit_error, mm_condition, rlang_error, error, condition]
failed to fit LMM: Fixed-effect design is rank-saturated: rank(X) = 2 and n = 2,
leaving zero residual degrees of freedom. Ordinary unpenalized LMM fitting is not
identifiable; use fewer fixed effects or an explicit penalized/MAP fixed-effect prior.
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to fit LMM: Fixed-effect design is rank-saturated: rank(X) = 2
and n = 2, leaving zero residual degrees of freedom. Ordinary unpenalized LMM fitting
is not identifiable; use fewer fixed effects or an explicit penalized/MAP fixed-effect
prior.
```

### Case 3 — one-row, binomial GLMM, lme4

```
lme4::glmer(y ~ x + (1 | subject), data = one_row_bin, family = binomial)
```
```
[ERROR] grouping factors must have > 1 sampled level
```

### Case 3 — one-row, binomial GLMM, mixeff::glmm

```
mixeff::glmm(y ~ x + (1 | subject), data = one_row_bin, family = binomial,
             control = mm_control(verbose = -1))
```
```
[ERROR class=mm_fit_error, mm_condition, rlang_error, error, condition]
failed to construct GLMM: Invalid argument: response is constant;
GLMM construction requires variation in the response
Caused by error in `doTryCatch()`:
! mm_fit_error: failed to construct GLMM: Invalid argument: response is constant;
GLMM construction requires variation in the response
```

---

## Analysis

### lme4 behavior

lme4 fires a **pre-fit guard** at the design-matrix stage and reports the
underlying structural problem in one terse, accurate sentence:
`"grouping factors must have > 1 sampled level"`. This catches the issue before
any numerical work begins, and the message names the structural cause precisely.

### mixeff behavior — bug found

**n=1, lmm():** mixeff does NOT fire a pre-fit guard for insufficient data / too
few group levels. Instead it passes the single-row frame into the Rust engine,
which detects that y has zero variance and reports `"Constant response: model
fitting failed"`. This is **misdiagnosed**: the actual problem is n=1 (or
equivalently 1 group level), not constant y. A user whose y column happened to
hold a constant would get the same message, masking the true cause. The condition
class (`mm_fit_error`) is correct and typed, but the message text is wrong.

**n=1, glmm():** Same issue — `"response is constant; GLMM construction requires
variation in the response"` — again misdiagnosing cause. A 1-row binary response
is always 0 or 1, which happens to have zero variance; the real fault is the
structural data deficiency.

**2 rows / 1 group, lmm():** Once y is varied, the Rust engine reaches a
different check and reports `"Fixed-effect design is rank-saturated"`. This is
_accurate_ for that specific data shape (intercept + x fills rank 2 with n=2
rows) but is a downstream symptom of the true problem (only 1 random-effect
group). lme4 catches the grouping-level constraint earlier and gives the more
direct message.

### Summary verdict

| Dimension | lme4 | mixeff |
|---|---|---|
| Pre-fit guard for 1 group level | Yes — direct, accurate | **No** — falls through to Rust engine |
| Message for n=1 | "grouping factors must have > 1 sampled level" | "Constant response: model fitting failed" (misdiagnosis) |
| Typed / catchable condition class | No (plain `simpleError`) | Yes (`mm_fit_error / mm_condition`) |
| Message clarity | Terse but correct | Typed and structured but **diagnoses the wrong root cause** |

The condition infrastructure (typed class, `mm_fit_error` inheritance, structured
`tryCatch`-ability) is better than lme4. The content of the message is worse: it
names a symptom observed by the Rust engine rather than the structural data
problem that R could and should detect before calling into Rust.

---

## Classification

**Verdict:** `bug` — mixeff is materially LESS clear than lme4 for this scenario.
The message misidentifies the cause (constant response / rank-saturated design
instead of "only 1 group level" or "fewer than 2 observations").

**Fix location:** `R/fit-lmm.R` (and `R/glmm.R`) — add a pre-fit guard after
`compile_model()` / `mm_translate_data()` that checks `nrow(data) >= 2` and that
each grouping factor has at least 2 distinct levels, raising an `mm_data_error`
with a message parallel to lme4's.

**Scope:** in-scope-missing — the "less inscrutable errors" promise (PRD §2)
requires accurate pre-fit guards for this basic structural check.
