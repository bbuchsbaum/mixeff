# lme4 / lmerTest Surface Survey — "lmerTest surface" capability family

**Date:** 2026-05-31  
**Packages surveyed:** lme4 2.0.1, lmerTest 3.2.1  
**Purpose:** Exhaustive reference of every user-facing function, argument, and
behavior in the lmerTest surface family, against which mixeff parity is judged.
This document describes lme4/lmerTest only; mixeff assessment is separate.

---

## 1. Overview

lmerTest extends lme4's `lmerMod` objects into a richer class
(`lmerModLmerTest`) that carries pre-computed Satterthwaite machinery
(Jacobians of `vcov(beta)` w.r.t. variance parameters, asymptotic
variance-covariance of variance parameters). Every inference function in this
family relies on that auxiliary data.

The user workflow is:

1. Fit with `lmerTest::lmer()` (a thin wrapper that upgrades the object class).
2. Call `summary()`, `anova()`, `contest*()`, `ranova()`, `ls_means()`, or
   `step()` directly on the result.
3. Existing `lmerMod` fits from `lme4::lmer()` can be upgraded with
   `as_lmerModLmerTest()`.

---

## 2. Core model class and fitting

### 2.1 `lmerTest::lmer()`

```r
lmer(formula, data = NULL, REML = TRUE, control = lmerControl(),
     start = NULL, verbose = 0L, subset, weights, na.action,
     offset, contrasts = NULL, devFunOnly = FALSE)
```

- Signature is identical to `lme4::lmer()`.
- Returns an object of class `lmerModLmerTest` (extends `lmerMod`).
- The upgrade from `lmerMod` to `lmerModLmerTest` happens silently: it
  computes and caches `vcov_varpar` (asymptotic var-cov of (theta, sigma)),
  `Jac_list` (Jacobians of `vcov(beta)` w.r.t. variance parameters),
  `vcov_beta`, and `sigma`.
- Users rely on this wrapper because **all lmerTest inference methods
  require the extra slots**; calling them on a plain `lmerMod` object silently
  falls back or errors.

### 2.2 `lmerModLmerTest` class slots (beyond `lmerMod`)

| Slot | Type | Content |
|------|------|---------|
| `vcov_varpar` | numeric matrix | Asymptotic var-cov of (theta, sigma) |
| `Jac_list` | list of matrices | Gradient of `vcov(beta)` w.r.t. each variance parameter |
| `vcov_beta` | numeric matrix | Asymptotic var-cov of fixed-effect parameters (beta) |
| `sigma` | numeric scalar | Residual standard deviation |

### 2.3 `as_lmerModLmerTest()`

```r
as_lmerModLmerTest(model, tol = 1e-08)
```

- Coerces a plain `lme4::lmer()` result (`lmerMod`) to `lmerModLmerTest`.
- `tol`: tolerance for deciding whether eigenvalues of the Hessian are
  negative, zero, or positive (used in the Satterthwaite Jacobian computation).
- Allows users with existing `lme4` workflows to access lmerTest methods
  without re-fitting.
- Returns an `lmerModLmerTest` object.

### 2.4 `merModLmerTest` (legacy)

A legacy S4 class from earlier lmerTest versions. Retained for backward
compatibility; current versions use `lmerModLmerTest`.

---

## 3. `summary()` — coefficient table with ddf

```r
## S3 method for class 'lmerModLmerTest'
summary(object, ..., ddf = c("Satterthwaite", "Kenward-Roger", "lme4"))
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `object` | — | An `lmerModLmerTest` fit |
| `...` | — | Passed to `lme4::summary.merMod` |
| `ddf` | `"Satterthwaite"` | Method for denominator df and t-statistics |

### `ddf` choices

| Value | Behaviour |
|-------|-----------|
| `"Satterthwaite"` | Welch-Satterthwaite approximation to denominator df; t-test for each fixed-effect coefficient |
| `"Kenward-Roger"` | Kenward-Roger scaled F→t statistics and df (via `pbkrtest::KRmodcomp`); adjusts for small-sample bias in variance estimation |
| `"lme4"` | Passes through to `lme4::summary.merMod` with no p-values; identical to calling summary on an unextended `lmerMod` |

### Output

- Returns a `summary.lmerModLmerTest` object (inherits from `summary.merMod`).
- The coefficient table (extractable with `coef(summary(model))`) has columns:
  `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)`.
- The `df` column contains Satterthwaite or KR denominator degrees of freedom
  per coefficient (not a single value — each row may differ).
- `print`, `coef`, and other `summary.merMod` methods apply.

### User reliance

Users call `summary()` to obtain p-values for each fixed-effect coefficient.
Without lmerTest the `df` and `Pr(>|t|)` columns are absent. This is the
single most commonly used entry point.

---

## 4. `anova()` — F-table for fixed effects

```r
## S3 method for class 'lmerModLmerTest'
anova(object, ..., type = c("III", "II", "I", "3", "2", "1"),
      ddf = c("Satterthwaite", "Kenward-Roger", "lme4"))
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `object` | — | An `lmerModLmerTest` fit |
| `...` | — | Additional `lmer` or `lm` model objects for model comparison |
| `type` | `"III"` | Sum-of-squares type (SAS terminology) |
| `ddf` | `"Satterthwaite"` | Denominator df method |

### `type` choices

| Value | Behaviour |
|-------|-----------|
| `"I"` / `"1"` | Sequential (Type I) ANOVA — each term adjusted for only those preceding it |
| `"II"` / `"2"` | Hierarchical (Type II) — each term adjusted for all other terms at the same level |
| `"III"` / `"3"` | Orthogonal (Type III) — each term adjusted for all other terms; default; SAS-compatible |

### Model comparison mode

When `...` contains additional model objects, `anova()` performs a likelihood
ratio test between nested models (passes through to `lme4::anova.merMod` with
`refit = TRUE` by default). In this mode `type` and `ddf` are ignored.

`lme4::anova.merMod` signature for comparison:
```r
anova(object, ..., refit = TRUE, model.names = NULL)
```
- `refit`: refit REML models with ML before comparing (prevents comparing
  incompatible REML log-likelihoods).
- `model.names`: character vector of labels for the table rows.

### Output

An ANOVA table with columns: `Sum Sq`, `Mean Sq`, `NumDF`, `DenDF`,
`F value`, `Pr(>F)`.

The `hypotheses` attribute stores the contrast matrices used for each term
(accessible via `show_tests()`).

### User reliance

The primary tool for testing whether fixed-effect terms as a whole are
significant. The `type` argument matters when terms are not orthogonal (e.g.,
unbalanced designs, interaction models).

---

## 5. `drop1()` — marginal single-term deletion

```r
## S3 method for class 'lmerModLmerTest'
drop1(object, scope, ddf = c("Satterthwaite", "Kenward-Roger", "lme4"),
      force_get_contrasts = FALSE, ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `object` | — | An `lmerModLmerTest` fit |
| `scope` | — | Terms to consider dropping; defaults to all marginal terms |
| `ddf` | `"Satterthwaite"` | Denominator df method |
| `force_get_contrasts` | `FALSE` | Internal; force re-derivation of contrast matrices |

### Output

A `drop1` table listing each term with `Sum Sq`, `F value`, `Df`, `Pr(>F)`.
This is the single-term-deletion analogue of `anova(type="III")` but considers
only terms that are marginally removable (i.e., not marginal to a retained
interaction).

### User reliance

Used for step-by-step model simplification; also called internally by `step()`.

---

## 6. `ranova()` / `rand()` — ANOVA-like table for random effects

```r
ranova(model, reduce.terms = TRUE, ...)
rand(model, reduce.terms = TRUE, ...)  # alias
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | An `lmerMod` or `lmerModLmerTest` fit |
| `reduce.terms` | `TRUE` | Reduce terms (e.g., `(f1+f2|gr)` → `(f2|gr)`) rather than remove them outright |
| `...` | — | Currently ignored |

### Reduction rules

| Term form | `reduce.terms = TRUE` | `reduce.terms = FALSE` |
|-----------|----------------------|----------------------|
| `(f1 + f2 \| gr)` | compared to `(f2\|gr)` and `(f1\|gr)` separately | removed entirely (3-df test) |
| `(f1 \| gr)` | reduced to `(1\|gr)` | removed |
| `(1 \| gr)` | removed (no further reduction possible) | removed |
| `(0 + f1 \| gr)` | reduced to `(1\|gr)` | removed |
| `(1 \| gr1/gr2)` | auto-expanded to `(1\|gr2:gr1)` + `(1\|gr1)` | same |

Structured covariance terms (e.g., `diag(...)`, `cs(...)`) always use
`reduce.terms = FALSE` behaviour.

If the model is REML-fitted, the LRTs are REML-LRTs.

### Output

An ANOVA-like table (inherits from `anova` and `data.frame`) with columns:

| Column | Meaning |
|--------|---------|
| `npar` | Number of model parameters |
| `logLik` | Log-likelihood (REML or ML to match fit) |
| `AIC` | AIC = -2*(logLik - npar) |
| `LRT` | Likelihood ratio test statistic (chi-square distributed asymptotically) |
| `Df` | Degrees of freedom for the LRT |
| `Pr(>Chisq)` | p-value |

### User reliance

The standard tool for testing whether random-effect terms contribute
significantly (e.g., "is the random slope needed?"). LRT at the boundary of
the parameter space (variance → 0) produces conservative p-values; lmerTest
does not apply boundary corrections.

---

## 7. `contest()` — flexible contrast test dispatcher

```r
## S3 method for class 'lmerModLmerTest'
contest(model, L, rhs = 0, joint = TRUE, collect = TRUE,
        confint = TRUE, level = 0.95,
        check_estimability = FALSE,
        ddf = c("Satterthwaite", "Kenward-Roger", "lme4"), ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | `lmerModLmerTest` fit |
| `L` | — | Contrast vector, matrix, or list thereof; length/ncol must equal `length(fixef(model))` |
| `rhs` | `0` | Hypothesised value (scalar, recycled) |
| `joint` | `TRUE` | `TRUE`: multi-df F-test; `FALSE`: per-row t-tests |
| `collect` | `TRUE` | Collect list results into a single matrix |
| `confint` | `TRUE` | Include lower/upper CI columns (only when `joint = FALSE`) |
| `level` | `0.95` | Confidence level |
| `check_estimability` | `FALSE` | Check estimability of each contrast (only effective when `joint = FALSE`); needed for rank-deficient designs |
| `ddf` | `"Satterthwaite"` | Denominator df method |

### Behaviour

- `L` may be a list of vectors/matrices; each element is tested separately.
- When `joint = TRUE`, calls `contestMD()` internally.
- When `joint = FALSE`, calls `contest1D()` for each row.
- `check_estimability`: when `TRUE` and the design matrix is rank-deficient,
  `L` must include columns for aliased (dropped) coefficients; contrast is
  marked `NA` if not estimable.
- Returns a `data.frame` or a list of `data.frame`s.

### User reliance

The high-level entry point for any custom contrast: equivalences of factor
levels, polynomial trends, pairwise differences not covered by `ls_means()`.

---

## 8. `contest1D()` — 1-df contrast t-test

```r
## S3 method for class 'lmerModLmerTest'
contest1D(model, L, rhs = 0,
          ddf = c("Satterthwaite", "Kenward-Roger"),
          confint = FALSE, level = 0.95, ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | `lmerModLmerTest` fit |
| `L` | — | Numeric contrast vector, length == `length(fixef(model))` |
| `rhs` | `0` | Hypothesised value |
| `ddf` | `"Satterthwaite"` | Denominator df method |
| `confint` | `FALSE` | Include CI columns |
| `level` | `0.95` | Confidence level |

### Output

A `data.frame` with one row: `Estimate`, `Std. Error`, `t value`, `df`,
`Pr(>|t|)`, and optionally `lower`, `upper`.

### User reliance

Direct, low-overhead entry point for a single t-test of a linear combination
of fixed effects. Used to test individual contrasts (e.g., treatment A − B)
with Satterthwaite df rather than relying on `qt(0.975, Inf)`.

---

## 9. `contestMD()` — multi-df contrast F-test

```r
## S3 method for class 'lmerModLmerTest'
contestMD(model, L, rhs = 0,
          ddf = c("Satterthwaite", "Kenward-Roger"),
          eps = sqrt(.Machine$double.eps), ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | `lmerModLmerTest` fit |
| `L` | — | Contrast matrix; `ncol == length(fixef(model))`, `nrow >= 1` |
| `rhs` | `0` | Hypothesised values; scalar or length `nrow(L)` |
| `ddf` | `"Satterthwaite"` | Denominator df method |
| `eps` | `sqrt(.Machine$double.eps)` | Tolerance for determining positive eigenvalues (→ numerator df) |

### Output

A `data.frame` with one row: `Sum Sq`, `Mean Sq`, `F value`, `NumDF`,
`DenDF`, `Pr(>F)`.

`NumDF` is the row-rank of `L` (rank-deficient `L` is accepted). `DenDF` is
the Satterthwaite or KR denominator df for the multi-df test.

### Relationship to `calcSatterth()`

`calcSatterth(model, L)` is a lower-level helper that computes only the
Satterthwaite denominator df (not exported formally in all versions). It
returns a list with `ddf` (denominator df) and related quantities.
`contestMD()` calls this internally.

### User reliance

Used when testing a hypothesis involving multiple fixed-effect parameters
simultaneously (e.g., "is there any effect of either Days or Days²?"). Also
the internal engine for `anova()` rows.

---

## 10. `ls_means()` / `lsmeansLT()` / `difflsmeans()` — LS-means

### `ls_means()` / `lsmeansLT()`

```r
## S3 method for class 'lmerModLmerTest'
ls_means(model, which = NULL, level = 0.95,
         ddf = c("Satterthwaite", "Kenward-Roger"),
         pairwise = FALSE, ...)

lsmeansLT(model, ...)  # alias for ls_means
```

### `difflsmeans()`

```r
## S3 method for class 'lmerModLmerTest'
difflsmeans(model, which = NULL, level = 0.95,
            ddf = c("Satterthwaite", "Kenward-Roger"), ...)
```
Equivalent to `ls_means(..., pairwise = TRUE)`.

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | `lmerModLmerTest` fit |
| `which` | `NULL` | Character vector of factor names; `NULL` = all factors |
| `level` | `0.95` | Confidence level |
| `ddf` | `"Satterthwaite"` | Denominator df method |
| `pairwise` | `FALSE` | Compute pairwise differences of LS-means instead of means |

### Statistical details

- LS-means (SAS terminology; equivalent to "estimated marginal means") are
  predicted means for each level of each factor, averaged with **equal weight**
  over all levels of other factors (flat/unweighted average).
- Numeric covariates are held at their mean.
- Non-estimable contrasts appear as `NA` rows.
- CI and p-values use the t-distribution with Satterthwaite or KR df.

### Output

An object of class `c("ls_means", "data.frame")` with attributes including
`hypotheses` (the underlying contrast matrices). Has a custom print method.

### User reliance

Quick way to obtain group means with inferential intervals for all categorical
factors, analogous to SAS `LSMEANS` / `LSMEANS / PDIFF`. The lmerTest
implementation is superseded in flexibility by `emmeans` but is simpler for
basic cases.

---

## 11. `show_tests()` — inspect hypothesis matrices

```r
## Default method
show_tests(object, fractions = FALSE, names = TRUE, ...)

## S3 method for class 'anova'
show_tests(object, fractions = FALSE, names = TRUE, ...)

## S3 method for class 'ls_means'
show_tests(object, fractions = FALSE, names = TRUE, ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `object` | — | An `anova` table or `ls_means` table carrying a `hypotheses` attribute |
| `fractions` | `FALSE` | Display matrix entries as exact fractions rather than decimals |
| `names` | `TRUE` | Include row and column names of hypothesis matrices |

### Output

A list of hypothesis (contrast) matrices — one per ANOVA term or LS-mean
factor.

### User reliance

Transparency/audit tool: lets users verify exactly which linear combinations
of fixed-effect parameters are being tested in each row of an `anova()` or
`ls_means()` table. Essential for checking Type II/III SS and LS-mean
correctness.

---

## 12. `step()` — backward elimination of random and fixed effects

```r
## S3 method for class 'lmerModLmerTest'
step(object, ddf = c("Satterthwaite", "Kenward-Roger"),
     alpha.random = 0.1, alpha.fixed = 0.05,
     reduce.fixed = TRUE, reduce.random = TRUE,
     keep, ...)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `object` | — | An `lmerModLmerTest` fit |
| `ddf` | `"Satterthwaite"` | Denominator df method for fixed-effect tests |
| `alpha.random` | `0.1` | Significance threshold for dropping random-effect terms |
| `alpha.fixed` | `0.05` | Significance threshold for dropping fixed-effect terms |
| `reduce.fixed` | `TRUE` | Whether to attempt fixed-effect elimination |
| `reduce.random` | `TRUE` | Whether to attempt random-effect elimination |
| `keep` | — | Character vector of fixed-effect term labels to protect from elimination; marginal terms are also protected |

### Algorithm

1. Random-effect elimination using `ranova()` (with `reduce.terms = TRUE`):
   terms with p > `alpha.random` are dropped in order from least significant,
   one at a time.
2. Fixed-effect elimination using `drop1()`: terms with p > `alpha.fixed` are
   dropped in order from least significant, one at a time.
3. The process continues until no further terms are eliminable.

### Output

A `step_list` object (list) with two elements:
- `$random`: `ranova`-like elimination table with an `Eliminated` column
  (integer giving elimination order; `0` = retained).
- `$fixed`: `drop1`-like elimination table with `Eliminated` column.

The final reduced model is an attribute of the `step_list`, accessible via
`get_model(<step_result>)`.

### `get_model()`

```r
get_model(x, ...)
```
Extracts the final fitted model from a `step_list` object. Returns an
`lmerModLmerTest` fit.

### User reliance

Automatic model selection for mixed models, analogous to `PROC MIXED` stepwise
in SAS. Controversial statistically (sequential testing) but widely used in
practice for exploratory analyses.

---

## 13. `calcSatterth()` — low-level Satterthwaite df computation

```r
calcSatterth(model, L)
```

### Arguments

| Argument | Default | Meaning |
|----------|---------|---------|
| `model` | — | An `lmerModLmerTest` fit |
| `L` | — | Contrast matrix (same dimension as in `contestMD`) |

### Output

A list containing `ddf` (Satterthwaite denominator df) plus intermediate
quantities (numerator df, F value, sum of squares). This is the internal
computation engine exposed for programmatic use.

### User reliance

Rarely called directly by end users; used by package developers and those
building custom inference wrappers. Provides the raw denominator df without the
full test-statistic machinery of `contestMD()`.

---

## 14. lme4 base methods inherited / relied upon

The lmerTest surface functions build on these lme4 methods (users encounter
them directly or indirectly):

| Function | Signature snippet | Role |
|----------|-------------------|------|
| `fixef()` | `fixef(object, ...)` | Extract fixed-effect coefficient vector (beta) |
| `ranef()` | `ranef(object, condVar = FALSE, ...)` | Extract BLUPs for random effects |
| `VarCorr()` | `VarCorr(x, sigma = 1, ...)` | Variance-covariance of random effects |
| `coef()` | `coef(object, ...)` | Per-group coefficient estimates (fixef + ranef) |
| `vcov()` | `vcov(object, ...)` | Variance-covariance matrix of fixed-effect estimates |
| `sigma()` | `sigma(object, ...)` | Residual standard deviation |
| `logLik()` | `logLik(object, REML = NULL, ...)` | Log-likelihood (REML or ML) |
| `deviance()` | `deviance(object, REML = NULL, ...)` | Deviance |
| `REMLcrit()` | `REMLcrit(object)` | REML criterion value |
| `confint()` | `confint(object, parm, level = 0.95, method = c("profile","Wald","boot"), ...)` | Confidence intervals via profile, Wald, or parametric bootstrap |
| `predict()` | `predict(object, newdata = NULL, re.form = NULL, random.only = FALSE, type = "link", allow.new.levels = FALSE, se.fit = FALSE, ...)` | Predictions at observed or new data |
| `simulate()` | `simulate(object, nsim = 1, seed = NULL, use.u = FALSE, ...)` | Parametric simulation |
| `residuals()` | `residuals(object, type = "response", ...)` | Model residuals |
| `fitted()` | `fitted(object, ...)` | Fitted values |
| `nobs()` | `nobs(object, ...)` | Number of observations |
| `ngrps()` | `ngrps(object, ...)` | Number of levels per grouping factor |
| `isREML()` | `isREML(x, ...)` | Whether fit used REML |
| `isSingular()` | `isSingular(x, tol = 1e-4, ...)` | Whether the fit is singular (variance ≈ 0) |
| `getME()` | `getME(object, name, ...)` | Extract named model components (e.g., `"X"`, `"Z"`, `"theta"`, `"Zt"`, etc.) |
| `model.matrix()` | `model.matrix(object, type = "fixed", ...)` | Fixed or random design matrices |
| `terms()` | `terms(x, fixed.only = TRUE, ...)` | Model terms object |
| `formula()` | `formula(x, fixed.only = FALSE, ...)` | Model formula |
| `update()` | `update(object, formula., ...)` | Refit with modified formula or data |
| `drop1()` | `drop1(object, scope, ...)` | Single-term deletion (lme4 base, overridden by lmerTest) |
| `profile()` | `profile(fitted, which = NULL, alphamax = 0.01, ...)` | Likelihood profile for `confint(method="profile")` |
| `rePCA()` | `rePCA(x, ...)` | PCA of random-effect var-cov; diagnoses singularity |
| `refit()` | `refit(object, newresp, ...)` | Refit with a new response vector |
| `refitML()` | `refitML(x, ...)` | Refit a REML model with ML |
| `bootMer()` | `bootMer(x, FUN, nsim, ...)` | Parametric bootstrap for `confint(method="boot")` |

---

## 15. Datasets bundled with lmerTest (used in documentation)

| Dataset | Description |
|---------|-------------|
| `sleepstudy` | (from lme4) Reaction times for sleep deprivation study; 180 obs, 18 subjects |
| `cake` | (from lme4) Cake-baking experiment; angle ~ recipe * temp |
| `ham` | Sensory evaluation of cured ham; 5-way design with random consumer effects |
| `carrots` | Sensory data for carrots; mixed design |
| `TVbo` | TV picture quality ratings; mixed model benchmark dataset |

---

## 16. External dependency: `pbkrtest`

The Kenward-Roger path in `summary()`, `anova()`, `contest*()`, `ls_means()`,
and `step()` calls `pbkrtest::KRmodcomp` internally. If `pbkrtest` is not
installed, KR methods fail with a namespace error. lmerTest treats `pbkrtest`
as a soft dependency (uses `requireNamespace`).

---

## 17. Argument cross-reference: `ddf` values

Every inference function that takes `ddf` accepts these values (partial
matching allowed):

| `ddf` value | Description | Requires |
|-------------|-------------|---------|
| `"Satterthwaite"` | Welch-Satterthwaite approximation; default everywhere | Only lmerTest slots |
| `"Kenward-Roger"` | Small-sample bias-corrected scaled F/t; recommended for small n | `pbkrtest` package |
| `"lme4"` | No p-values; pass-through to base lme4 output | Nothing extra |

---

## 18. Object class hierarchy

```
lmerModLmerTest
  └── lmerMod
        └── merMod
              └── (R5 reference class)
```

`summary.lmerModLmerTest` inherits from `summary.merMod`; `print`, `coef`,
`confint` methods for `summary.merMod` all apply.

The `ls_means` result is `c("ls_means", "data.frame")` with a `hypotheses`
attribute.

A `step_list` object has a custom `print` method and carries the final model
as an attribute.

---

## 19. Function-to-output quick reference

| Function | Output class | Key columns / slots |
|----------|-------------|---------------------|
| `summary()` | `summary.lmerModLmerTest` | Coef table: `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)` |
| `anova()` (single model) | `anova` data.frame | `Sum Sq`, `Mean Sq`, `NumDF`, `DenDF`, `F value`, `Pr(>F)` |
| `anova()` (model comparison) | `anova` data.frame | `Df`, `AIC`, `BIC`, `logLik`, `deviance`, `Chisq`, `Chi Df`, `Pr(>Chisq)` |
| `drop1()` | `anova` data.frame | `Df`, `Sum Sq`, `F value`, `Pr(>F)` |
| `ranova()` | `anova` data.frame | `npar`, `logLik`, `AIC`, `LRT`, `Df`, `Pr(>Chisq)` |
| `contest()` (joint) | `data.frame` | `Sum Sq`, `Mean Sq`, `F value`, `NumDF`, `DenDF`, `Pr(>F)` |
| `contest()` (per-row) | `data.frame` | `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)`, optionally `lower`, `upper` |
| `contest1D()` | `data.frame` | `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)` |
| `contestMD()` | `data.frame` | `Sum Sq`, `Mean Sq`, `F value`, `NumDF`, `DenDF`, `Pr(>F)` |
| `ls_means()` | `c("ls_means","data.frame")` | `Estimate`, `Std. Error`, `df`, `t value`, `Pr(>|t|)`, `lower`, `upper` |
| `difflsmeans()` | `c("ls_means","data.frame")` | Same columns as `ls_means(pairwise=TRUE)` |
| `step()` | `step_list` | `$random` (ranova-like), `$fixed` (drop1-like), each with `Eliminated` column |
| `get_model()` | `lmerModLmerTest` | Final model after backward elimination |
| `show_tests()` | list of matrices | One contrast matrix per ANOVA term or LS-mean factor |
| `as_lmerModLmerTest()` | `lmerModLmerTest` | Upgraded lme4 model object |
| `calcSatterth()` | list | `ddf` (denominator df) and intermediate quantities |

---

*End of lme4 / lmerTest surface survey.*
