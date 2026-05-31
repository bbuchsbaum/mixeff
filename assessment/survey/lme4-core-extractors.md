# lme4 / lmerTest — Core Estimate Extractors: Reference Surface

Survey date: 2026-05-31  
Packages confirmed: lme4 2.0.1, lmerTest 3.2.1  
Purpose: exhaustive reference against which mixeff is judged. This document
describes only lme4/lmerTest; mixeff assessment is in a separate artefact.

---

## 1. `fixef(object, add.dropped = FALSE, noScale = NULL, ...)`

**What it does.** Returns the maximum-likelihood (ML) or REML estimates of the
fixed-effects parameters as a *named numeric vector*.

**Arguments users rely on.**

| Argument | Default | Why users use it |
|---|---|---|
| `object` | — | the fitted `merMod` |
| `add.dropped` | `FALSE` | rank-deficient models drop aliased columns; `TRUE` reinserts them as `NA` so the result is full-length and aligned with the original design matrix |
| `noScale` | `NULL` | when autoscaling was applied during model fit, `TRUE` returns parameters on the original (unscaled) scale |

**Return value.** Named numeric vector; names are coefficient labels
(`"(Intercept)"`, `"Days"`, etc.).  Length equals the number of non-dropped
fixed effects (or the full design-matrix width when `add.dropped = TRUE`).

**lmerTest interaction.** `fixef` itself is unchanged by lmerTest. However
`coef(summary(model))` gains `df`, `t value`, and `Pr(>|t|)` columns (Satterthwaite
or Kenward-Roger denominator df) after loading lmerTest.

**Live example (sleepstudy).**
```
(Intercept)        Days
  251.40510    10.46729
```

---

## 2. `ranef(object, condVar = TRUE, drop = FALSE, whichel = names(ans), postVar = FALSE, ...)`

**What it does.** Extracts the *conditional modes* (BLUPs) of the random effects
— for linear models these are also the conditional means.

**Arguments users rely on.**

| Argument | Default | Why users use it |
|---|---|---|
| `condVar` | `TRUE` | attach the `"postVar"` attribute — a `p × p × k` array of per-level conditional variance-covariance matrices; essential for caterpillar plots and uncertainty quantification |
| `drop` | `FALSE` | when `TRUE`, single-column data frames (intercept-only RE terms) are simplified to named numeric vectors |
| `whichel` | all groups | character vector selecting a subset of grouping factors |
| `postVar` | `FALSE` | deprecated synonym for `condVar`; exists for back-compat |

**Return value.** Object of class `"ranef.mer"` — a *named list of data frames*,
one per grouping factor.  Each data frame has rows = levels of the factor,
columns = random-effect terms.  When `condVar = TRUE`, each data frame carries
a `"postVar"` attribute: a `p × p × n` array (or list of such arrays when more
than one RE term maps to the same grouping factor).

**`as.data.frame` method on `ranef.mer`.**
Converts to a long-format data frame with columns:
`grpvar`, `term`, `grp`, `condval`, `condsd`.
`condsd` is derived from the diagonal of the `postVar` array.
This is the primary bridge to ggplot2 caterpillar plots.

**Visualization methods on `ranef.mer`.**
- `dotplot(ranef.mer)` — lattice caterpillar plot with CIs
- `qqmath(ranef.mer)` — Q-Q plot against normal quantiles

**Live example (`(Days|Subject)`).**
```
$Subject            (Intercept)       Days
          308       2.258551   9.198976
          309     -40.398738  -8.619685
...
```

---

## 3. `coef(object, ...)`  — `coef.merMod`

**What it does.** Returns *subject-specific (BLUP) coefficients* for each level
of each grouping factor: the sum of the fixed-effect estimate and the
corresponding conditional mode (i.e., `fixef + ranef`).

**Return value.** Object of class `"coef.mer"` — a named list of data frames,
one per grouping factor; columns correspond to fixed-effect terms that appear in
that group's random structure.  For `(Days|Subject)` the result is an 18 × 2
data frame with `(Intercept)` and `Days` columns.

**`plot.coef.mer`** method exists for visualisation but has no `as.data.frame`
equivalent in the class.

**Key identity.**
```
coef(fit)$Subject[i, j]  ==  fixef(fit)[j] + ranef(fit)$Subject[i, j]
```

**Live example (first three subjects).**
```
    (Intercept)      Days
308    253.6637 19.666262
309    211.0064  1.847605
310    212.4447  5.018429
```

---

## 4. `VarCorr(x, sigma = 1, ...)` and its methods

### 4a. `VarCorr.merMod`

**What it does.** Extracts estimated variances, standard deviations, and
correlations for all random-effects terms, plus the residual standard deviation
for LMMs.

**Arguments.**

| Argument | Default | Why users use it |
|---|---|---|
| `sigma` | `1` | scale multiplier applied to standard deviations (rarely changed by end users, but useful for custom reporting) |

**Return value.** Object of class `"VarCorr.merMod"` — a named list of matrices
(one per grouping factor).  Each matrix is the *variance-covariance matrix* for
that term; its `"stddev"` and `"correlation"` attributes carry SD and
correlation forms.  The scalar residual SD is stored as attribute `"sc"` on
the top-level list.  Additional attributes `"theta"`, `"rho"`, and `"profpar"`
expose the underlying parameterization.

### 4b. `print.VarCorr.merMod`

```r
print(vc, digits = max(3L, getOption("digits") - 2L),
          comp = "Std.Dev.", formatter = format, ...)
```

`comp` accepts any combination of `"Variance"` and `"Std.Dev."` to control
which scale is shown.

### 4c. `as.data.frame.VarCorr.merMod`

```r
as.data.frame(x, row.names = NULL, optional = FALSE,
              order = c("cov.last", "lower.tri"), ...)
```

Produces one row per variance or covariance parameter (plus residual):

| Column | Content |
|---|---|
| `grp` | grouping factor name or `"Residual"` |
| `var1` | first variable name (`NA` for Residual row) |
| `var2` | second variable name (`NA` for variance rows) |
| `vcov` | variance or covariance value |
| `sdcor` | standard deviation or correlation |

`order = "cov.last"` (default): variances first, covariances appended.
`order = "lower.tri"`: lower-triangle order matching the matrix layout.

This data frame is the standard input for downstream tools (broom, ggplot2,
custom tables).

---

## 5. `sigma(object, ...)`  — `sigma.merMod`

**What it does.** Extracts the *residual standard deviation* (dispersion
parameter) from the fitted model.

**Return value.** Scalar numeric.

**Family-specific behaviour.**

| Family | sigma value | Interpretation |
|---|---|---|
| Gaussian (LMM) | estimated | residual SD; `sqrt(MSE)` analogue |
| Gamma | `1/sqrt(shape)` | scale of the variance function |
| Inverse Gaussian | `1/sqrt(lambda)` | shape-based |
| Binomial, Poisson | `1` (fixed) | no free scale parameter |

**Notes.** `getME(object, "sigma")` is identical; `sigma()` is preferred.
The value appears as `"Residual Std.Dev."` in `print.VarCorr`.

---

## 6. `fitted(object, ...)` — `fitted.merMod`

**What it does.** Returns *conditional fitted values* — predictions at the
observed covariate values conditioned on the estimated conditional modes of the
random effects.

**Return value.** Named numeric vector; names are row indices of the original
data frame (character, e.g., `"1"`, `"2"`, …).  Length equals `nobs(object)`.

**Relationship to `predict`.**
`fitted(object)` is identical to `predict(object, re.form = NULL)` — both
return conditional (level-1 BLUP-based) fitted values.

---

## 7. `predict(object, newdata = NULL, newparams = NULL, re.form = NULL, random.only = FALSE, terms = NULL, type = c("link", "response"), allow.new.levels = FALSE, na.action = na.pass, se.fit = FALSE, ...)`

**What it does.** Generates predictions from the fitted model, optionally for
new data and/or with different random-effect structures.

**Arguments users rely on.**

| Argument | Default | Why users use it |
|---|---|---|
| `newdata` | `NULL` | predict for new observations |
| `re.form` | `NULL` (all REs) | `NA` or `~0` → marginal/population prediction; a formula → include named REs only |
| `random.only` | `FALSE` | ignore fixed effects; useful for extracting pure RE contribution |
| `type` | `"link"` | `"response"` back-transforms via inverse link (GLMM) |
| `allow.new.levels` | `FALSE` | `TRUE` allows previously unseen group levels (uses unconditional RE = 0) |
| `se.fit` | `FALSE` | experimental; returns list with `$fit` and `$se.fit` |
| `newparams` | `NULL` | override fitted parameters (theta and/or beta) for sensitivity analysis |

**Return value.** Named numeric vector of predicted values (length = `nrow(newdata)`
or `nobs(object)` when `newdata = NULL`).  When `se.fit = TRUE`, a two-element
list `$fit` + `$se.fit`.

**Note on SEs.** lme4 marks `se.fit` as experimental; `bootMer` is the
recommended path for prediction SEs because it incorporates uncertainty in
variance parameters.

---

## 8. `residuals(object, type = ..., scaled = FALSE, ...)` — `residuals.merMod`

**What it does.** Returns residuals from the fitted model.

### 8a. For LMMs (`lmerMod`)

Default type: `"response"` (matches `residuals.lm`).

| `type` | Formula | Note |
|---|---|---|
| `"response"` | `y - fitted(object)` | default; raw residuals on the response scale |
| `"pearson"` | `(y - fitted) / sqrt(var(y))` | scaled by variance function; for LMM equals response residuals (variance is constant) |
| `"deviance"` | same as `"response"` for LMMs | deviance residuals coincide with response for Gaussian |
| `"working"` | same as `"response"` for LMMs | PIRLS working residuals; equal to response for Gaussian |
| `"partial"` | **not yet implemented** | raises an error |

**`scaled` argument.** When `scaled = TRUE`, residuals are divided by the
residual standard deviation `sigma(object)`.  This replicates `lme` (nlme)
behaviour.  Note: lme4's Pearson residuals are *not* scaled by default, unlike
`residuals.lme`.

### 8b. For GLMMs (`glmerMod`)

Default type: `"deviance"` (matches `residuals.glm`).

| `type` | Content |
|---|---|
| `"deviance"` | signed square root of contribution to deviance (default) |
| `"pearson"` | `(y - mu) / sqrt(V(mu))` where `V` is the variance function |
| `"response"` | `y - mu` on the response (probability/count) scale |
| `"working"` | PIRLS working residuals from the final IWLS iteration |
| `"partial"` | **not yet implemented** |

**Return value.** Named numeric vector; names are row indices; length = `nobs`.

---

## 9. `vcov(object, ...)` — `vcov.merMod`

**What it does.** Returns the variance-covariance matrix of the *fixed-effects*
parameter estimates.

**Return value.** A `dpoMatrix` (positive-definite, from the Matrix package);
dimension `p × p` where `p` is the number of fixed-effect parameters.
Dimnames are the coefficient names.

**Derived quantities.** Standard errors are `sqrt(diag(vcov(object)))`.
Used internally by `confint(method="Wald")`, `emmeans`, and custom contrast
machinery.

**Note.** lme4 exposes only the fixed-effect vcov via this generic.  Covariance
of the variance-component parameters (theta) is not directly available via
`vcov`; it can be obtained via `profile` + `confint`.

---

## 10. Additional closely related extractors (not core but adjacent)

These appear in `methods(class = "merMod")` and are part of the same user
workflow:

| Function | What it returns |
|---|---|
| `logLik(object, REML = NULL)` | `"logLik"` object with `df` and `nobs` attributes |
| `deviance(object, REML = NULL)` | scalar deviance (absolute, minus twice log-lik) |
| `nobs(object)` | integer; number of observations |
| `ngrps(object)` | named integer vector; number of levels per grouping factor |
| `df.residual(object)` | residual degrees of freedom |
| `model.matrix(object, type = c("fixed","random","randomListRaw"), noScale = NULL)` | fixed-effects design matrix (type="fixed") or random-effects matrix |
| `model.frame(object, fixed.only = FALSE)` | the model frame slot |
| `formula(x, fixed.only = FALSE, random.only = FALSE)` | formula |
| `getME(object, name)` | low-level access to internal slots: `"X"`, `"Z"`, `"theta"`, `"beta"`, `"b"`, `"u"`, `"sigma"`, `"fixef"`, `"flist"`, `"L"`, `"Lambda"`, `"Lambdat"`, `"devcomp"`, etc. (38 named components + `"ALL"`) |
| `REMLcrit(object)` | REML criterion at optimum |
| `isREML(object)` | logical; was model fitted by REML? |
| `isSingular(object, tol)` | logical; boundary/singular fit? |

---

## 11. lmerTest-specific additions (beyond base lme4)

lmerTest adds Satterthwaite (default) or Kenward-Roger denominator-df
computation for fixed-effect inference.  The additions surface via:

| Entry point | What changes |
|---|---|
| `summary(lmerMod)` | `coef(summary(.))` gains `df`, `t value`, `Pr(>|t|)` columns |
| `anova(lmerMod)` | gains `DenDF`, `F value`, `Pr(>F)` columns |
| `contest(model, L, ...)` | new function; general linear contrast test with Satterthwaite df |
| `ranova(model)` | random-effects ANOVA table |

`fixef`, `ranef`, `VarCorr`, `sigma`, `fitted`, `residuals`, `vcov`, and
`coef` are **unchanged** by lmerTest; they remain lme4 methods.

---

## 12. Tolerance conventions for statistical equivalence

Per the mixeff PRD §10 parity standard, these are the accepted tolerances when
comparing mixeff output to lme4 reference values on parity datasets:

| Quantity | Tolerance |
|---|---|
| Fixed effects (`fixef`) | 1e-4 |
| Variance-component parameters (`theta`) | 1e-3 |
| Log-likelihood (`logLik`) | 1e-3 |
| Residual SD (`sigma`) | 1e-4 |

---

*This document is a survey artefact. It describes what lme4 provides; it does
not assess mixeff coverage. See the companion parity assessment for the gap
analysis.*
