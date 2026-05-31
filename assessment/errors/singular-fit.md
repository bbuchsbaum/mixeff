# Error-message quality probe: singular-fit

**Scenario:** Models that converge to a singular/boundary/reduced-rank
covariance matrix — the canonical "model is too complex for the data" case.

**Probe script:** `assessment/errors/probe-singular-fit.R`
**Date:** 2026-05-31
**mixeff version:** installed (crate 0.1.0)
**lme4 version:** 2.0.1

---

## Three sub-scenarios tested

### Scenario A — random slope on all-identical x (zero slope variance)

This scenario turned out **not** to be singular: both lme4 and mixeff converge
in the interior with a small-but-positive slope variance (~0.013).  Neither
package raised a warning.  Both agree on the non-singular diagnosis.
**Not a probe of singular messaging; included for completeness.**

---

### Scenario B — two groups, random intercept collapses to boundary = 0

**lme4 output (verbatim):**
```
isSingular: TRUE
 Groups   Name        Std.Dev.
 subject  (Intercept) 0.000   
 Residual             1.143   
```
lme4 produces no warning or message at fit time. `isSingular()` returns `TRUE`
after the fact. The user must know to call `isSingular()` manually or check the
VarCorr table for a zero std dev. No explanation is given.

**mixeff output (verbatim):**
```
fit_status: converged_reduced_rank
is_singular: TRUE

Linear mixed model fit by REML
Formula: y ~ x + (1 | subject)
Fit status: converged_reduced_rank
...
Fitted covariance state:
The fitted covariance matrix is rank-deficient.
  r0: requested rank 1; fitted effective rank 0.
Use changes(fit) to see which dimension was unsupported.
Use random_options(spec, group = subject) to inspect lower-dimensional covariance choices.

VarCorr:
   group        name variance std_dev correlation       note
 subject (Intercept)        0       0             [boundary]
[boundary]: variance component is at the boundary of the parameter space.
```

**Comparison:** mixeff is clearly better here. It:
- Names the problem at the top: `Fit status: converged_reduced_rank`
- Explains it inline in `print()`: "The fitted covariance matrix is rank-deficient."
- Quantifies it: "requested rank 1; fitted effective rank 0."
- Flags the affected row in VarCorr: `[boundary]`
- Points to concrete remediation verbs: `changes(fit)` and `random_options()`

lme4 is silent at fit time. The user gets a zero std.dev. and must call
`isSingular()` separately to confirm.

---

### Scenario C — overparameterised random structure (n_obs <= n_random_effects)

**lme4 output (verbatim):**
```
Error: number of observations (=20) <= number of random effects (=30)
for term (1 + x + z | subject); the random-effects parameters and the
residual variance (or scale parameter) are probably unidentifiable
```
lme4 **throws a hard error** and refuses to return a fitted object.

**mixeff output (verbatim):**
```
fit_status: converged_reduced_rank
is_singular: TRUE

Linear mixed model fit by REML
Formula: y ~ x + z + (1 + x + z | subject)
Fit status: converged_reduced_rank
...
Fitted covariance state:
The fitted covariance matrix is rank-deficient.
  r0: requested rank 3; fitted effective rank 1.
Use changes(fit) to see which dimension was unsupported.
Use random_options(spec, group = subject) to inspect lower-dimensional covariance choices.

VarCorr:
   group        name variance  std_dev correlation note
 subject (Intercept) 0.131934 0.363227                 
 subject           x 0.334642 0.578483       +1.00     
 subject           z 0.331625 0.575869 +1.00 +1.00     
```

**Comparison — needs-work signal:**
mixeff does not crash or panic (good). It returns a model object with a clear
`converged_reduced_rank` status and explains the rank deficit (requested 3,
effective 1). The VarCorr table shows all correlations fused to ±1.00, which is
informative.

However there is a notable gap relative to lme4's behavior in this scenario:
lme4 refuses outright because the model is mechanically non-identifiable
(n_obs=20 ≤ n_random_effects=30). mixeff silently fits and returns coefficients
without flagging the identifiability problem by name. The user sees a
rank-deficient diagnosis but **not** the statement that the model is
"probably unidentifiable." For this class of overparameterised model, the
identifiability framing is more actionable than the rank-deficiency framing
alone.

Additionally, the `VarCorr` all-correlations-at-1.00 pattern is a strong
indicator of a non-identified covariance structure, but mixeff does not surface
a message connecting this pattern to identifiability.

---

## Classification summary

| Sub-scenario | Gap | Classification |
|---|---|---|
| A (non-singular interior) | None — both agree | works |
| B (boundary variance = 0) | mixeff clearly better than lme4 | works |
| C (overparameterised, n_obs ≤ n_re) | mixeff fits silently; missing identifiability language | needs-work |

---

## Assessment

**Overall verdict: needs-work**

For the most common singular-fit case (B: a variance component hitting zero),
mixeff's messaging is **materially better** than lme4: the problem is named,
quantified, flagged in the VarCorr table, and paired with remediation verbs.
No crashes, no panics, no silent wrong answers.

For the extreme overparameterised case (C: structural non-identifiability),
mixeff fits and reports a rank-deficient solution rather than refusing.
The rank-deficit language is accurate, but it does not name identifiability
as the root cause. lme4's error message ("probably unidentifiable") is more
precise here. This is a `needs-work` item, not a bug — the fit object is
internally consistent and the audit-first design makes the diagnosis
inspectable — but the messaging could be strengthened with an identifiability
note when `n_obs <= n_random_effects`.

**No panics, no segfaults, no stack-trace-only errors, no silent wrong answers.**
