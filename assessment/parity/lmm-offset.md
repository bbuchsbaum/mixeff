# Parity cell: lmm-offset

**Formula:** `y ~ x + offset(o) + (1|g)`
**Dataset:** simulated (n=100, 20 groups of 5; true offset coefficient = 1.0)
**Date:** 2026-05-31
**Probe script:** `assessment/parity/lmm-offset-probe.R`

---

## lme4/lmerTest result

| Quantity     | Value |
|---|---|
| converged    | TRUE |
| `(Intercept)` | 1.4061 |
| `x`          | 0.6855 |
| SE(Intercept) | 0.2140 |
| SE(x)        | 0.0543 |
| theta        | 1.8702 |
| sigma        | 0.4977 |
| logLik       | -102.904 |
| AIC          | 213.808 |
| BIC          | 224.229 |

lme4 correctly incorporates the offset: the intercept (1.41) and slope (0.69) match the true DGP values (1.5 and 0.7), with the offset absorbing its unit contribution.

---

## mixeff result

mixeff **refused** at formula parsing time with a structured diagnostic:

```
mm_formula_error: in-formula construct `offset(...)` at position 8 is not in
the engine's stateless transform subset (allowed: `I(<+ - * / ^, unary -,
parens, literals, columns>)` and pointwise `log`/`log2`/`log10`/`exp`/`sqrt`/
`abs`). Stateful transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`,
`center`, …) carry fitting-time state and must be precomputed as data columns
or handled by the host wrapper.
```

Error class chain: `mm_formula_error > mm_condition > rlang_error > error`.

No fit was attempted. The wrapper's formula manifest does not list `offset` in
`formula_features$transformations` (only `implicit_intercept`,
`nested_grouping_expansion`, `interaction_grouping`).

---

## Comparison table

| Quantity | lme4 | mixeff | |diff| | Verdict |
|---|---|---|---|---|
| fixef `(Intercept)` | 1.4061 | — | — | mixeff-error |
| fixef `x`           | 0.6855 | — | — | mixeff-error |
| SE(Intercept)       | 0.2140 | — | — | mixeff-error |
| SE(x)               | 0.0543 | — | — | mixeff-error |
| theta               | 1.8702 | — | — | mixeff-error |
| sigma               | 0.4977 | — | — | mixeff-error |
| logLik              | -102.904 | — | — | mixeff-error |
| AIC                 | 213.808  | — | — | mixeff-error |
| BIC                 | 224.229  | — | — | mixeff-error |
| convergence         | TRUE   | formula refused | — | mixeff-error |

---

## Diagnosis

**Classification:** `in-scope-missing`

`offset()` is a standard lme4/glm term used to fix a linear predictor
contribution (e.g. log-exposure in Poisson models, or any known linear shift).
It is not a "stateful transform" in the sense that poly/scale/bs are — its
semantics are: "extract the vector `o` from data and add it to the linear
predictor with coefficient fixed at 1." The engine's error message incorrectly
groups `offset()` with stateful transforms; the real issue is that the formula
parser does not recognise `offset(...)` as a special term and the Rust engine
has no offset-vector pathway.

**Workaround exists?** Partially. A user could precompute `y_adj = y - o` and
fit `y_adj ~ x + (1|g)`, which reproduces the same fixed-effect estimates but
with a shifted intercept. However this is not equivalent for models where the
offset has a coefficient != 1 (e.g. GLMM Poisson with `offset(log(n))`), and
it is not a transparent or documented workaround.

**Severity:** `major`

- offset() is required for any rate/exposure model (Poisson GLMM, negative
  binomial, any trial-count model)
- lme4 supports it transparently
- The error message is honest but the classification of `offset` as a
  "stateful transform" is misleading — it is not stateful in the same sense
  as poly/splines
- No `offset=` argument path exists on `lmm()` or `glmm()` as an alternative

**What needs to change:**

1. The R wrapper should detect `offset(...)` terms in the formula, extract the
   offset vector from the model frame, and either (a) pass it as a separate
   numeric vector to the Rust engine alongside the design matrix, or (b) apply
   it as a y-adjustment before fitting (only valid for Gaussian/LMM).
2. The Rust engine needs an offset-vector slot in its LMM/GLMM fit payload, or
   the wrapper must handle it transparently.
3. The formula manifest should expose an `offset` capability flag so downstream
   code can branch cleanly.

---

## Raw output (verbatim)

```
=== Dataset summary ===
n=100, n_groups=20
y: mean=1.902, sd=1.276
x: mean=0.033, sd=1.041
o: mean=0.474, sd=0.271

=== lme4 fit ===
lme4 converged: TRUE
fixef:
(Intercept)           x 
  1.4060612   0.6854671 
SE:
(Intercept)           x 
 0.21400377  0.05431564 
theta: 1.870203 
sigma: 0.4976892 
logLik: -102.9041 
AIC: 213.8083 
BIC: 224.229 
ranef(g) head: -1.52588 0.4165383 0.8608904 1.660357 -0.8187092 
fitted head: 1.680217 0.3065242 0.3281309 1.368532 0.4572615 

=== mixeff fit ===
mixeff ERROR: in-formula construct `offset(...)` at position 8 is not in the
engine's stateless transform subset (allowed: `I(<+ - * / ^, unary -, parens,
literals, columns>)` and pointwise `log`/`log2`/`log10`/`exp`/`sqrt`/`abs`).
Stateful transforms (`poly`, `scale`, `ns`, `bs`, `cut`, `factor`, `center`,
…) carry fitting-time state and must be precomputed as data columns or handled
by the host wrapper.
mixeff ERROR class: mm_formula_error, mm_condition, rlang_error, error, condition

=== Comparison ===
lme4 succeeded, mixeff FAILED/REFUSED
Error class: mm_formula_error, mm_condition, rlang_error, error, condition

=== VERDICT: mixeff-error ===

=== mixeff formula manifest ===
formula_features$operators:
[1] "+"  "-"  "*"  ":"  "/"  "&"  "|"  "||"
formula_features$transformations:
[1] "implicit_intercept"        "nested_grouping_expansion"
[3] "interaction_grouping"
```
