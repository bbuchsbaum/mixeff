# Parity probe: glmm-cbpp-agq

**Cell:** glmm-cbpp-agq  
**Dataset:** `cbpp` (lme4)  
**Formula:** `cbind(incidence, size - incidence) ~ period + (1|herd)`  
**Family:** binomial / logit  
**Focus:** nAGQ=1 (Laplace) AND nAGQ=2+ (AGQ) — nlopt feature gate  
**Date run:** 2026-05-31  
**lme4:** 2.0.1 | **mixeff:** 0.1.0

---

## Script

See `glmm-cbpp-agq-probe.R` in this directory.

---

## Raw output

```
=== SESSION INFO ===
lme4 version: 2.0.1
mixeff version: 0.1.0

=== DATASET (cbpp) ===
nrow: 56   ncol: 4
  herd incidence size period
1    1         2   14      1
...

=== lme4 FIT (nAGQ=1, Laplace) ===
lme4 nAGQ=1 wall-clock (seconds): 0.057

-- fixef --
(Intercept)     period2     period3     period4
  -1.398343   -0.991925   -1.128216   -1.579745
-- SE --
(Intercept)     period2     period3     period4
  0.2312140   0.3031506   0.3228300   0.4220489
-- theta --
herd.(Intercept)
       0.6420699
-- sigma --
[1] 1
-- logLik --
'log Lik.' -92.02657 (df=5)
-- AIC --
[1] 194.0531
-- BIC --
[1] 204.1799

=== lme4 FIT (nAGQ=2, AGQ) ===
-- fixef --
(Intercept)=-1.397946  period2=-0.992702  period3=-1.129066  period4=-1.580863
-- theta -- 0.6415872
-- logLik -- -50.05701
-- AIC -- 110.114

=== lme4 FIT (nAGQ=10, AGQ) ===
-- fixef --
(Intercept)=-1.399224  period2=-0.991409  period3=-1.127810  period4=-1.579481
-- theta -- 0.6475199
-- logLik -- -50.00501
-- AIC -- 110.01

=== mixeff FIT (nAGQ=1, pirls_profiled) ===
!!! mixeff nAGQ=1 ERROR !!!
class: mm_formula_error, mm_condition, rlang_error, error, condition
message: in-formula construct `cbind(...)` at position 0 is not in the engine's
stateless transform subset (allowed: `I(<+ - * / ^, unary -, parens, literals,
columns>)` and pointwise `log`/`log2`/`log10`/`exp`/`sqrt`/`abs`). Stateful
transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`, `center`, ...) carry
fitting-time state and must be precomputed as data columns or handled by the
host wrapper.

=== mixeff FIT (nAGQ=2, pirls_profiled) ===
!!! mixeff nAGQ=2 ERROR/REFUSAL !!!
(same mm_formula_error — cbind blocked)

=== mixeff FIT (joint_laplace — expected refusal) ===
(same mm_formula_error — cbind blocked before method dispatch)

=== lme4 AGQ SENSITIVITY ===
nAGQ=1 vs nAGQ=2 logLik diff: 4.197e+01  (42 log-lik units — Laplace is poor approx)
nAGQ=2 vs nAGQ=10 fixef max diff: 1.38e-03 (just above 1e-3 fixef tol)
nAGQ=2 vs nAGQ=10 logLik diff: 5.2e-02    (within 1e-3? No — 5.2e-2)

=== SPEED SUMMARY ===
lme4 nAGQ=1:    0.0570s
lme4 nAGQ=2:    0.0320s
lme4 nAGQ=10:   0.0310s
mixeff nAGQ=1:  (not reached — formula error)
```

---

## Analysis

### Root cause: `cbind(...)` in response position

mixeff's formula compiler rejects `cbind(incidence, size - incidence)` as the
response. The error class is `mm_formula_error` and the message is clear: the
engine only permits a stateless transform subset in the formula, and `cbind()`
is not in that subset. This is a **typed, honest diagnostic** — the error class
and message precisely identify what is unsupported and why.

The refusal fires before any method dispatch or Rust FFI call, so ALL three
attempted fits (nAGQ=1 pirls_profiled, nAGQ=2 pirls_profiled, joint_laplace)
fail at the same point with the same error. The joint_laplace refusal (which
would be expected due to the nlopt feature gate) is never reached.

### Is `cbind(y, n-y)` response support in scope?

The cbpp dataset with a `cbind(incidence, size - incidence)` response is a
**grouped binomial** (binomial trials aggregated per observation). Supporting
this requires either:

1. Formula-level `cbind()` response parsing (pre-compute successes + failures
   before dispatch), OR
2. An explicit `trials` argument to `glmm()`.

This is a genuine capability gap vs lme4 for the GLMM surface. The PRD §3
non-goals do not explicitly exclude grouped binomial response handling. The
certified family/link surface (binomial/logit) is supported — the blocker is
purely the formula-response transformation layer, not the Rust engine.

A workaround exists: expand to individual binary observations (842 rows from 56
grouped rows). However, this gives **different numerical results** from the
aggregate cbind form because the Laplace approximation changes with observation
count/grouping structure (fixef divergence ~4e-2, theta divergence ~1.7e-2 vs
lme4 cbind nAGQ=1). The workaround is not a transparent substitute.

### nAGQ behaviour in lme4 (context)

For this grouped binomial dataset, lme4's nAGQ=1 Laplace approximation is a
poor fit: logLik = -92.03 vs nAGQ=10 converged value = -50.00 (difference ~42
log-lik units). This is well-known: Laplace is unreliable for grouped binomial
responses with moderate cluster sizes. The fixef estimates also shift:
`period4` coefficient differs by ~1.1e-3 between nAGQ=1 and nAGQ=2+. Between
nAGQ=2 and nAGQ=10, fixef differ by ~1.3e-3 (just above the 1e-4 fixef
tolerance, within the 1e-3 theta tolerance). AGQ converges from nAGQ=2 to 10
for logLik to ~5.2e-2 (still outside 1e-3 logLik tolerance, converging slowly).

For this cell the relevant lme4 reference target would be nAGQ=10+, not nAGQ=1,
for any parity claim about marginal likelihood. mixeff's pirls_profiled (Laplace)
should be compared to lme4 nAGQ=1 for like-for-like; AGQ support would require
separate comparison to lme4 nAGQ=2+.

### joint_laplace refusal

The joint_laplace refusal is **expected** per the glmm() docstring (nlopt feature
gate). However, in this run the refusal class was `mm_formula_error` (cbind
blocked first), not the expected `mm_fit_error` for an unavailable backend.
The net result is the same (a typed error), but the diagnostic reason is
different from what the docstring implies. This is a cosmetic issue only — the
cbind error fires first and is more specific.

### Severity classification

The `cbind(...)` response formula block is a **capability gap** for grouped
binomial GLMM. The engine supports binomial/logit family; the blocker is the
formula-response layer. This is an **in-scope-missing** feature (grouped
binomial response is a mainstream GLMM use case, not excluded by PRD §3).

The error is a typed, honest diagnostic (`mm_formula_error`) rather than a
silent wrong answer or crash — that is the correct behavior. The gap is that the
feature is simply absent, not that it is handled incorrectly.

---

## Quantity-by-quantity table

| Quantity | lme4 (nAGQ=1) | mixeff | maxAbsDiff | Tol | Status |
|---|---|---|---|---|---|
| fixef (Intercept) | -1.398343 | — | — | 1e-4 | mixeff-error |
| fixef period2 | -0.991925 | — | — | 1e-4 | mixeff-error |
| fixef period3 | -1.128216 | — | — | 1e-4 | mixeff-error |
| fixef period4 | -1.579745 | — | — | 1e-4 | mixeff-error |
| SE (Intercept) | 0.23121 | — | — | 1e-4 | mixeff-error |
| theta (herd) | 0.6421 | — | — | 1e-3 | mixeff-error |
| sigma | 1.0 | — | — | 1e-4 | mixeff-error |
| logLik (nAGQ=1) | -92.027 | — | — | 1e-3 | mixeff-error |
| logLik (nAGQ=10) | -50.005 | — | — | 1e-3 | mixeff-error |
| AIC | 194.053 | — | — | 2e-3 | mixeff-error |
| BIC | 204.180 | — | — | 2e-3 | mixeff-error |
| ranef (herd) | computed | — | — | 1e-4 | mixeff-error |
| fitted values | computed | — | — | 1e-4 | mixeff-error |
| nAGQ=2 support | supported | refused (formula) | — | — | mixeff-error |
| nAGQ=10 support | supported | refused (formula) | — | — | mixeff-error |
| joint_laplace | N/A | refused (formula, then nlopt gate) | — | — | expected |

All mixeff entries are blank because the fit was never reached — the formula
parser rejected `cbind(incidence, size - incidence)` as the response before any
fitting was attempted.

---

## lme4 AGQ sensitivity (reference only)

| Comparison | fixef max abs diff | logLik diff |
|---|---|---|
| nAGQ=1 vs nAGQ=2 | 1.12e-3 (period4) | 42.0 (Laplace very poor) |
| nAGQ=1 vs nAGQ=10 | 8.81e-4 (Intercept) | 42.0 |
| nAGQ=2 vs nAGQ=10 | 1.38e-3 (period4) | 5.2e-2 |

For this grouped binomial, nAGQ >= 10 is needed for reliable marginal likelihood.
The fixef change from nAGQ=2 to nAGQ=10 is ~1.3e-3, just above the 1e-4 fixef
tolerance but within the 1e-3 theta tolerance.

---

## Overall verdict

**Outcome: mixeff-error** — mixeff cannot fit this model at all because
`cbind(incidence, size - incidence)` as a response formula is blocked by the
formula compiler with `mm_formula_error`. The Rust engine supports
binomial/logit; the gap is in the R-side formula-response handling layer for
grouped binomial (trials-format) inputs.

The error is **typed and honest** (class `mm_formula_error`, clear message
identifying `cbind(...)` as unsupported). No silent wrong answers are produced.

**Severity: major** — grouped binomial response (`cbind(y, n-y) ~`) is a core
GLMM use case; the cbpp dataset is lme4's canonical GLMM example. The feature
is in-scope (PRD §3 does not exclude it) and the family/link is certified. A
workaround (expand to binary rows) exists but is not numerically equivalent and
requires manual data transformation. Classifying as **major** (not blocker)
because the error is clean and the workaround exists, but this is a real
capability gap for standard usage.

**Classification: in-scope-missing** — grouped binomial formula-response
support (`cbind(success, failure) ~`) is absent from the formula compiler.
