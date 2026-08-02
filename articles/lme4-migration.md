# Migrating from lme4

``` r

library(mixeff)
```

`mixeff` fits lme4-style formulas with the familiar extractors and, on
its supported envelope, gives statistical answers tested against pinned
lme4 references within documented tolerances. It is not a literal
drop-in — you call
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) /
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) rather
than `lmer()` / `glmer()`, and requests outside the envelope are changed
or declined with a diagnostic, never silently. This vignette is the
verb-for-verb, argument-for-argument map; the generated tables below are
the envelope’s exact boundary. It is the migration entry in the
package’s six-vignette set: the support contract itself is stated in
[`vignette("mixeff", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/mixeff.md),
the fitting workflows in
[`vignette("lmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/lmm.md)
and
[`vignette("glmm", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md),
convergence and boundary semantics in
[`vignette("convergence", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/convergence.md),
and the inference contract in
[`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

## The two edits

An `lmer` script becomes an `lmm` script with two changes: the fitting
verb and the control object.

``` r

fit <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy_data(),
           control = mm_control(verbose = -1))
fixef(fit)
#> (Intercept)        Days 
#>   251.40510    10.46729
```

``` r

m <- lme4::lmer(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)
lme4::fixef(m)
#> (Intercept)        Days 
#>   251.40510    10.46729
```

(`sleepstudy_data()` above just returns
[`lme4::sleepstudy`](https://rdrr.io/pkg/lme4/man/sleepstudy.html) when
lme4 is installed; use
[`lme4::sleepstudy`](https://rdrr.io/pkg/lme4/man/sleepstudy.html)
directly in your own code.)

## Verb map

| lme4 | mixeff | Notes |
|----|----|----|
| `lmer(y ~ x + (x \| g), data)` | `lmm(y ~ x + (x \| g), data)` | same formula language, incl. `(x\|\|g)`, `(1\|g1/g2)`, crossed |
| `glmer(y ~ ..., family = binomial)` | `glmm(y ~ ..., family = binomial())` | pass a family **object** ([`binomial()`](https://rdrr.io/r/stats/family.html)), not a string |
| `lmerControl(...)` / `glmerControl(...)` | `mm_control(verbose=, max_feval=)` | optimizer/tolerance knobs are engine-chosen (see below) |
| `fixef`, `ranef`, `VarCorr`, `coef`, `sigma`, `vcov` | identical | same generics |
| `logLik`, `AIC`, `BIC`, `deviance`, `nobs`, `confint` | identical | `confint` supports Wald, profile, and bootstrap on LMMs; Wald/profile/bootstrap availability per GLMM fit type is the route matrix in [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md) |
| `predict`, `fitted`, `residuals`, `simulate`, `refit` | identical (LMM) | [`predict()`](https://rdrr.io/r/stats/predict.html) supports `re.form = NULL/NA`, `se.fit`, `interval`; [`simulate()`](https://rdrr.io/r/stats/simulate.html)/[`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) are LMM-only and refuse GLMM fits with a typed error |
| `update(fit, . ~ . - x)` | identical | formula edits, `REML=`, `weights=`, etc. |
| `anova(m1, m2)`, `drop1`, `getME`, `ngrps`, `isSingular` | identical | `isSingular()` is [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md) |
| `broom.mixed::tidy/glance/augment` | identical | registered for `mm_lmm`/`mm_glmm` |
| `emmeans::emmeans(fit, ~ x)` | identical | mixeff registers an emmeans basis on the same covariance payload, so estimates and SEs agree; mixeff’s status/reliability/reason audit columns do not survive the emmeans trip, and a default profiled GLMM refuses inference through emmeans too |
| `lmerTest` p-values in [`summary()`](https://rdrr.io/r/base/summary.html) | built in | Satterthwaite/Kenward-Roger native, no extra package |

## Argument map for `lmm()` / `glmm()`

| lme4 argument | mixeff | Notes |
|----|----|----|
| `REML` | `lmm(..., REML=)` | same |
| `weights` | `weights=` | LMM and GLMM |
| `offset` | `glmm(..., offset=)` | GLMM only; LMM in-fit offset is not yet supported |
| `subset` | `lmm(..., subset=)` | supported for [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) |
| `na.action` | `lmm(..., na.action=)` | **default refuses NA**; pass `na.action = na.omit` for lme4’s complete-case behavior |
| `contrasts` | partial | unordered factors use treatment coding, ordered factors `contr.poly` (both matching R/lme4 defaults); other codings are refused — recode the factor. Contrast attributes set on factor columns and `options(contrasts=)` are currently ignored — the fit uses treatment coding — a known gap (bd-01KZ050ADJ2Y7SS20YH9CTE2NS); recode factors explicitly when porting |
| numeric/character grouping (accepted silently) | coerced to factor | non-factor grouping columns are coerced to factor, announced via an `mm_grouping_coercion_notice`; silence it with `mm_control(verbose = -1)` |
| `family = "binomial"` | `family = binomial()` | string families are not accepted |
| `nAGQ` | `glmm(..., nAGQ=)` | `>1` on the profiled path |
| `control = lmerControl(optimizer=, optCtrl=)` | `mm_control(optimizer=, max_feval=, ...)` | the engine picks a default optimizer; [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md) can override it or cap the evaluation budget |
| `start` | `mm_control(start=)` | theta warm starts |

## Four differences to check when porting

**1. Coefficient names match lme4 exactly.** Since 0.2.0,
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`summary()`](https://rdrr.io/r/base/summary.html) tables,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) dimnames, and
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
weight names use lme4’s naming and column order (`"recipeB"`,
`"temperature.L"`, `"recipeB:temperature.L"`), so lme4 code keyed on
coefficient names carries over without renaming. (Earlier versions used
an engine encoding like `"recipe: B"`; if you wrote normalization shims
for those, delete them.)

**2. Grouped binomial responses.**
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
accepts the `cbind(successes, failures)` spelling like `glmer`:

``` r

glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
     lme4::cbpp, family = binomial())
```

**3. The default GLMM estimator is not glmer’s.**
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
defaults to a fast profiled (PIRLS) estimator whose coefficients do
**not** match `glmer()` exactly; it prints a notice saying so. For
glmer-equivalent (joint Laplace) estimates, ask for them:

``` r

glmm(y ~ x + (1 | g), data, family = binomial(),
     method = "joint_laplace")
```

**4. `||` with a factor means full decorrelation.** In mixeff, `||`
fixes *every* covariance in the block at zero — including the
covariances among a factor’s level contrasts (each treatment-coded
contrast gets an independent variance), matching `MixedModels.jl`’s
`zerocorr()`. These are mixeff’s own semantics, stated directly: lme4
documents its `||` as intended for numeric predictors, and as currently
implemented it leaves factor terms intact (the factor keeps its full
within-factor covariance block) rather than promising a factor behavior.
So `(1 + cond + x || subj)` with a factor `cond` fits a strictly larger
model in lme4 than in mixeff. The two therefore disagree on the
parameter count, hence on `df` and AIC — and on the optimum itself
whenever the fitted within-factor covariance is non-zero. mixeff
announces the situation at compile time with an info diagnostic
(`covariance_assumption`, reason `double_bar_factor_term`). To reproduce
lme4’s model family exactly, write the expansion explicitly and give the
factor its own correlated block:

``` r

# mixeff `||`: independent variances for every column, factor levels included
glmm(y ~ cond * x + (1 + cond + x || subj), data, family = binomial())

# lme4-equivalent family: the factor keeps its within-factor covariance block
glmm(y ~ cond * x + (1 | subj) + (0 + cond | subj) + (0 + x | subj),
     data, family = binomial())
```

## The supported envelope, generated

The tables below are rendered from the package’s machine-readable
support registry (`inst/support/model-support.json`, via
[`supported_models()`](https://bbuchsbaum.github.io/mixeff/reference/supported_models.md)),
so this vignette cannot drift from what the code enforces — the test
suite holds the fitting code, this registry, and the
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
documentation to the same envelope.

| family | links | estimators | notes |
|:---|:---|:---|:---|
| binomial | logit, probit, cloglog | pirls_profiled, joint_laplace | Bernoulli 0/1, cbind(successes, failures), and proportion + trial weights responses. |
| poisson | log, sqrt | pirls_profiled, joint_laplace |  |
| Gamma | log | pirls_profiled, joint_laplace |  |
| negative_binomial | log | pirls_profiled | NB2 via mm_negative_binomial(); theta estimated or fixed. Intervals via confint(method = “bootstrap”). |

Supported GLMM family/link/estimator cells. {.table}

| feature | lmm | glmm | notes |
|:---|:---|:---|:---|
| case weights | supported | supported | GLMM binomial weights double as trial counts for proportion responses. |
| offset | unsupported | supported | lmm() has no offset argument; glmm(offset =) is supported. |
| simulate | supported | refused | GLMM simulation is not implemented in this version. |
| refit | supported | refused | Refit GLMMs by calling glmm() on the modified data. |
| subset / random / custom na.action / custom contrasts (GLMM) | see notes | refused | Reserved arguments on glmm(); recode data before fitting. LMMs: default refuses NA unless na.action = na.omit; non-default contrast codings are refused at compile time. |
| marginal means verbs (mm_means / mm_comparisons / mm_grid / mm_predictions / test_effect) | supported | refused | GLMM marginal means go through the emmeans bridge on certified fits. |
| double-bar \|\| with factor terms | mixeff semantics | mixeff semantics | mixeff fixes every covariance in the block at zero, including among a factor’s level contrasts (matching MixedModels.jl zerocorr()); lme4’s \|\| leaves factor terms intact. Announced at compile time (covariance_assumption / double_bar_factor_term); the migration vignette gives the lme4-equivalent expansion. |

Feature support (see notes for mixeff-specific semantics). {.table}

Explicit non-goals: nonlinear mixed models; arbitrary GLM family objects
or custom links; adaptive Gaussian quadrature as a general estimator
contract (nAGQ \> 1 is a profiled-path sensitivity mode).

## What is `NA`-with-a-reason (and why)

Where an inference method is unavailable, `mixeff` returns `NA` with a
machine-readable reason or raises a typed condition, rather than
printing a number without comment:

| Situation | lme4 (+lmerTest) | mixeff |
|----|----|----|
| `NA` in a model variable | silently dropped | refused unless `na.action = na.omit` |
| Boundary (singular) fit | one-time warning | persistent `[boundary]` tag + effective rank |
| Satterthwaite df at a boundary | lmerTest prints df anyway | refused with a reason; use bootstrap |
| Prediction SE for an unseen grouping level | — | `NA` with the reason in `mm_reason` |
| GLMM [`confint()`](https://rdrr.io/r/stats/confint.html) | computed | Wald intervals on `joint_laplace` fits via `method = "asymptotic"` (lowercase; lme4’s `method = "Wald"` spelling is not accepted); parametric-bootstrap intervals via `method = "bootstrap"` on profiled-estimator fits (like `bootMer`; refused on `joint_laplace` fits at this engine pin); Wald routes refused on the default `pirls_profiled` fit |

Use `inference_options(fit)` to see, before you run anything, which
inference routes are available on a given fit and why.
