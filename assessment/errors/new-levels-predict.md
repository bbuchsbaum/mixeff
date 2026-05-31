# Error-message quality probe: new-levels-predict

**Scenario:** `predict(newdata)` where `newdata` contains a grouping level
(`subject = "99"`) that was never seen during training, and
`allow.new.levels = FALSE` (the default).

**Date probed:** 2026-05-31  
**mixeff version:** installed from /Users/bbuchsbaum/code/mixeff  
**lme4 version:** 2.0.1  

---

## Probe script

`assessment/errors/probe-new-levels-predict.R`

Key setup:

```r
train   <- data.frame(y = ..., x = ..., subject = factor(rep(1:10, each=4)))
fit_lme4 <- lmer(y ~ x + (1 | subject), data = train)
fit_mm   <- lmm (y ~ x + (1 | subject), data = train, control = mm_control(verbose=-1))

newdata_new_level <- data.frame(x = 0.5, subject = factor("99"))

predict(fit_lme4, newdata = newdata_new_level, allow.new.levels = FALSE)
predict(fit_mm,   newdata = newdata_new_level, allow.new.levels = FALSE)
```

---

## Verbatim messages

### lme4

```
Condition class: simpleError, error, condition
Message        : new levels detected in newdata: 99
```

### mixeff

```
Condition class: mm_inference_unavailable, mm_condition, rlang_error, error, condition
Message        : predict_new failed: Invalid argument: New level '99' in grouping
                 factor 'subject'. Use NewReLevels::Population or ::Missing to allow this.
Caused by error in `doTryCatch()`:
! mm_inference_unavailable: predict_new failed: Invalid argument: New level '99'
  in grouping factor 'subject'. Use NewReLevels::Population or ::Missing to allow this.
```

---

## allow.new.levels = TRUE (positive control)

Both engines succeed and return numerically identical predictions:

```
lme4   allow.new.levels=TRUE: 2.0868
mixeff allow.new.levels=TRUE: 2.0868
```

---

## Assessment

### Information content comparison

| Dimension | lme4 | mixeff |
|---|---|---|
| Class | generic `simpleError` | typed `mm_inference_unavailable` |
| Names the new level | yes ("99") | yes ("'99'") |
| Names the grouping factor | no | yes ("grouping factor 'subject'") |
| Suggests the fix | no | yes (use allow.new.levels=TRUE / Population or Missing) |
| Machine-catchable typed condition | no | yes |

mixeff's message is **clearer than lme4's** on every dimension:

- lme4 says only `new levels detected in newdata: 99` — a terse, untyped
  `simpleError` with no indication of which grouping factor is affected and
  no suggestion of how to proceed.
- mixeff names both the offending level (`'99'`) and the grouping factor
  (`subject`), and explicitly tells the caller what to do: set
  `allow.new.levels = TRUE` (the R-side spelling) which maps to the
  `NewReLevels::Population` or `::Missing` policies documented in the error.
- The condition class `mm_inference_unavailable` is machine-catchable,
  allowing callers to distinguish "the engine refused to certify this" from
  argument errors or data-shape problems.

### Leakage of internal implementation detail

One minor rough edge: the suggestion reads `Use NewReLevels::Population or
::Missing to allow this.` — these are **Rust enum variant names**, not the
R-level argument. The R user needs `allow.new.levels = TRUE`, not
`NewReLevels::Population`. The message leaks the Rust API surface through the
R boundary, which is confusing for users who don't know about the Rust layer.

This is a **needs-work** issue, not a bug: the message is already more
informative than lme4's, but the actionable instruction should be
`allow.new.levels = TRUE` rather than a Rust enum name.

### Pre-probe bug found (separate from the new-levels scenario)

When `newdata` contains the response column with `NA_real_`, mixeff raises an
`mm_data_error` about a non-finite value in the `y` column **before** it ever
reaches the new-level check. This masks the new-level error entirely when a
caller naively passes `newdata` that includes the response column. lme4 is
tolerant of `y = NA` in `newdata` for predict and ignores it. This secondary
issue is `needs-work`; it is distinct from the new-level message quality and
is filed separately.

---

## Verdict

**Verdict:** needs-work  
**Quality vs lme4:** clearer-than-lme4 (more informative, typed, actionable)  
**Needs-work reason:** The fix suggestion leaks Rust enum names
(`NewReLevels::Population`) rather than the R-side argument
(`allow.new.levels = TRUE`). Recommend rewriting the Rust error string (or
the `mm_abort_from_bridge` post-processing step) to say:
`"set allow.new.levels = TRUE to use population-mean predictions for unseen levels."`.

A secondary data-validation bug (NA response column in newdata triggers
mm_data_error before the new-level check) should also be fixed so the
correct new-level error surfaces regardless of whether the response is
present in newdata.
