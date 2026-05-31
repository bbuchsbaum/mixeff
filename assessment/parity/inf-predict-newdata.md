# Parity Assessment: inf-predict-newdata

**Cell:** inf-predict-newdata
**Date:** 2026-05-31
**Dataset:** sleepstudy (lme4 built-in, N=180, 18 subjects x 10 days)
**Formula:** `Reaction ~ Days + (Days|Subject)`
**REML:** TRUE
**Focus:** predict() on newdata with/without re.form, allow.new.levels

---

## Environment

| Package | Version |
|---------|---------|
| lme4    | 2.0.1   |
| mixeff  | 0.1.0   |

---

## Script

Written to `/tmp/inf_predict_newdata_probe.R`.

---

## Raw Output

```
=== ENVIRONMENT ===
lme4 version: 2.0.1
mixeff version: 0.1.0

=== FIT ===
lme4  wall time: 0.024 sec
mixeff wall time: 0.012 sec

=== BASELINE FIT QUANTITIES ===
Max |diff| fixef: 5.684342e-13 (tol 1e-4)
Max |diff| sigma: 0.00143999 (tol 1e-4)
Max |diff| logLik: 4.133908e-06 (tol 1e-3)

Newdata (known subjects):  10 rows, 5 subjects
Newdata (new subjects):  3 rows, subject = S99

=== PREDICT: in-sample (newdata=NULL, re.form=NULL) ===
lme4  fitted[1:5]: 253.6637 273.3299 292.9962 312.6624 332.3287
mixeff fitted[1:5]: 253.6542 273.3226 292.991 312.6593 332.3277
Max |diff| in-sample fitted: 0.01142123 (tol 1e-4)

=== PREDICT: in-sample population-level (re.form=NA) ===
lme4  fixed-fitted[1:5]: 251.4051 261.8724 272.3397 282.807 293.2742
mixeff fixed-fitted[1:5]: 251.4051 261.8724 272.3397 282.807 293.2742
Max |diff| in-sample population-level: 1.136868e-12 (tol 1e-4)

=== PREDICT: newdata known subjects, conditional (re.form=NULL) ===
lme4  pred[1:5]: 253.6637 273.3299 211.0064 212.854 212.4447
mixeff pred[1:5]: 253.6542 273.3226 211.0117 212.858 212.4468
Max |diff| conditional newdata (known subjects): 0.009470825 (tol 1e-4)

=== PREDICT: newdata known subjects, population-level (re.form=NA) ===
lme4  pop pred[1:5]: 251.4051 261.8724 251.4051 261.8724 251.4051
mixeff pop pred[1:5]: 251.4051 261.8724 251.4051 261.8724 251.4051
Max |diff| pop-level newdata (known subjects): 5.684342e-13 (tol 1e-4)

=== PREDICT: newdata known subjects, re.form=~0 ===
Max |diff| re.form=~0 newdata: 5.684342e-13 (tol 1e-4)
~0 identical to NA path: TRUE

=== PREDICT: new subject, allow.new.levels=FALSE (should error) ===
lme4  response: ERROR: new levels detected in newdata: S99
mixeff response: ERROR: mm_data_error - failed to add numeric column 'Reaction':
  Invalid argument: numeric column `Reaction` contains a non-finite value (NaN)
  at index 0; reject NaN/Inf before fitting or use add_numeric_unchecked

=== PREDICT: new subject, allow.new.levels=TRUE ===
mixeff ERROR (allow.new.levels=TRUE): failed to add numeric column 'Reaction':
  Invalid argument: numeric column `Reaction` contains a non-finite value (NaN)
  at index 0; reject NaN/Inf before fitting or use add_numeric_unchecked
lme4  pred (new subject): 251.4051 303.7415 345.6107

=== se.fit=TRUE (unavailable in mixeff) ===
mixeff se.fit result: list with $fit (len 180) and $se.fit all NA: TRUE
se.fit unavailable reason: prediction_se_unavailable_phase_2

=== interval='confidence' (should raise mm_inference_unavailable) ===
mixeff interval='confidence': ERROR: mm_inference_unavailable - `interval` prediction
  requires prediction standard errors, which are not certified by the current Rust
  inference contract.

=== re.form = ~(1|Subject) (unsupported partial conditioning) ===
lme4 partial re.form: numeric vector length 10
mixeff partial re.form: mixeff ERROR: mm_inference_unavailable - `re.form` values
  other than NULL, NA, or ~0 are not supported by the current Rust prediction contract.

=== SPEED: newdata conditional predict (10 reps) ===
lme4  newdata conditional predict (10 reps): 0.025 sec
mixeff newdata conditional predict (10 reps): 0.009 sec
Speed ratio (mixeff/lme4): 0.36

=== SUMMARY TABLE ===
Quantity                                                        Max|diff|  Within? Note
------------------------------------------------------------------------------------------
fixef (baseline)                                              5.68434e-13      YES
sigma (baseline)                                               0.00143999       NO
logLik (baseline)                                             4.13391e-06      YES
in-sample fitted (re.form=NULL)                                 0.0114212       NO
in-sample pop-level (re.form=NA)                              1.13687e-12      YES
newdata conditional known subjects                             0.00947082       NO
newdata population-level known subjects (re.form=NA)          5.68434e-13      YES
newdata population-level known subjects (re.form=~0)          5.68434e-13      YES
newdata allow.new.levels=TRUE new subjects                       ERROR/NA      N/A
```

---

## Second-pass investigation

The sigma discrepancy (0.00144) and the conditional fitted-value divergence (~0.011)
are linked. The Rust optimizer converges to slightly different theta than lme4:

```
lme4  theta: 0.9667418  0.01516906  0.23091
mixeff theta: 0.9668414  0.0151315   0.2310602
Max |diff| theta: 0.0001502   (tol 1e-3 — within tolerance)
```

The theta difference is within the 1e-3 tolerance, but it propagates into sigma
(tol 1e-4) and hence into the conditional predictions (which embed BLUPs that depend
on sigma and theta). This is an optimizer convergence divergence at the theta boundary
(0.0001502 vs 1e-3 tol), not a correctness bug in the prediction machinery itself.

**allow.new.levels fix confirmed:** The failure for the new-subject test was caused by
passing `NA` as the response column in newdata. The mixeff Rust bridge validates all
numeric columns, including the response, before it can filter it out. When the same
newdata is passed with a dummy non-NA Reaction value (e.g., `0`), `allow.new.levels=TRUE`
works correctly and matches lme4 to machine precision:

```
lme4  preds (S99, Days=0,5,9): 251.4051 303.7415 345.6107
mixeff preds (S99, Days=0,5,9): 251.4051 303.7415 345.6107
Max |diff|: 1.136868e-12
```

---

## Findings Summary

| # | Scenario | lme4 behaviour | mixeff behaviour | Max|diff| | Status |
|---|----------|---------------|-----------------|-----------|--------|
| 1 | in-sample conditional (re.form=NULL) | 180 fitted values | 180 fitted values | 0.01142 | DIVERGENT — beyond 1e-4 tol; driven by sigma/theta offset |
| 2 | in-sample population (re.form=NA) | 180 fixed-fitted | 180 fixed-fitted | 1.1e-12 | within-tol |
| 3 | newdata known subj, conditional | 10 predictions | 10 predictions | 0.00947 | DIVERGENT — same theta/sigma root cause |
| 4 | newdata known subj, population (re.form=NA) | 10 predictions | 10 predictions | 5.7e-13 | within-tol |
| 5 | newdata known subj, population (re.form=~0) | 10 predictions | 10 predictions | 5.7e-13 | within-tol |
| 6 | new subj, allow.new.levels=FALSE | error (new levels) | error (wrong error: NaN in Reaction) | — | DEFECT — wrong error class/message; leaks NaN guard instead of level policy |
| 7 | new subj, allow.new.levels=TRUE | 3 pop predictions | error (NaN in Reaction) | — | DEFECT — fails when response col is NA; lme4 strips response before passing to engine |
| 8 | se.fit=TRUE | numeric SEs | NA vector + attribute | — | Honest deferred refusal (correct) |
| 9 | interval='confidence' | works | mm_inference_unavailable | — | Honest deferred refusal (correct) |
| 10 | re.form=~(1\|Subject) partial | 10 predictions | mm_inference_unavailable | — | Honest deferred refusal (correct) |
| 11 | sigma (baseline) | 25.5918 | 25.5904 | 0.00144 | DIVERGENT — beyond 1e-4 tol; theta within 1e-3 tol |

---

## Root Cause Analysis

### Finding A: sigma out-of-tolerance (severity: minor)
`sigma` diff = 0.00144, tol = 1e-4. Root: theta diff = 0.000150 (within 1e-3 tol).
The sigma tolerance (1e-4) is tighter than the theta tolerance (1e-3), and sigma is
a deterministic function of theta. The optimizer converges to a slightly different
saddle. Not a logic bug — the tolerances are inconsistent: if theta is accepted at
1e-3, sigma should be accepted at a looser bound. The conditional predictions inherit
this discrepancy (~0.011 max diff).

### Finding B: allow.new.levels fails when response column contains NA (severity: major)
The `predict_new` FFI passes the full newdata frame — including the response column —
through `mm_translate_data`. The Rust bridge applies its NaN/Inf guard to ALL numeric
columns before distinguishing predictors from response. lme4 strips the response before
model-matrix construction, so `NA` in the response column is harmless. mixeff does not.

This affects two sub-cases:
- `allow.new.levels=FALSE`: raises `mm_data_error` about NaN rather than the expected
  "new levels" diagnostic. The error class is wrong.
- `allow.new.levels=TRUE`: raises the same NaN error instead of returning predictions.

**Workaround confirmed:** passing a dummy numeric (e.g., `0`) in the response column
makes both paths work. But users who naturally pass `NA` for the response when predicting
(which is the canonical idiom) will hit this bug.

The fix is in `mm_predict_conditional_newdata()` (R/predict.R line 188): strip the
response column from `new_data` before passing it to `mm_translate_data`, or pass a
dummy value. Alternatively the Rust bridge should drop the response column when building
the prediction design matrix.

### Finding C: partial re.form (severity: none / by design)
lme4 supports `re.form = ~(1|Subject)` for partial random-effect conditioning. mixeff
raises `mm_inference_unavailable` with a clear message. This is an honest, documented
deferred refusal (PRD scope). Correct behaviour.

---

## Tolerance Cross-Reference

| Quantity | Tolerance | Max|diff| | Result |
|----------|-----------|-----------|--------|
| fixef | 1e-4 | 5.7e-13 | within-tol |
| sigma | 1e-4 | 0.00144 | **DIVERGENT** |
| theta | 1e-3 | 0.000150 | within-tol |
| logLik | 1e-3 | 4.1e-6 | within-tol |
| in-sample conditional fitted | 1e-4 | 0.01142 | **DIVERGENT** |
| in-sample population fitted | 1e-4 | 1.1e-12 | within-tol |
| newdata conditional pred | 1e-4 | 0.00947 | **DIVERGENT** |
| newdata population pred | 1e-4 | 5.7e-13 | within-tol |

---

## Speed

mixeff newdata conditional predict is ~2.8x faster than lme4 (0.009s vs 0.025s, 10 reps).
Note: mixeff's conditional newdata predict refits the model from scratch via the Rust FFI
(formula + training data + newdata in one call), yet still beats lme4.

---

## Overall Classification

**Outcome: mixed**

- Population-level predictions (re.form=NA or ~0): **works**, within-tol.
- Conditional predictions (re.form=NULL): **within-tol on fixef** but the sigma/theta
  offset propagates ~0.011 error into BLUPs — beyond the 1e-4 tolerance.
- allow.new.levels with NA response: **in-scope defect** (major) — wrong error path,
  user-hostile error message, prediction unavailable for canonical NA-response newdata.
- Honest refusals (se.fit, interval, partial re.form): **correct**.

**Severity: major** — the allow.new.levels bug blocks a primary use case (predict for
new groups / held-out subjects with NA response), and the conditional prediction
divergence exceeds the stated 1e-4 fitted-value tolerance on every test.
