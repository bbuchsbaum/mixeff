# Migrating from lme4

``` r

library(mixeff)
```

`mixeff` aims to be *functionally equivalent* to `lme4`: the same
formula language, the same extractor surface, and statistical answers
that agree within documented tolerances. It is **not** a literal drop-in
— you call
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) /
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) rather
than `lmer()` / `glmer()`, results are not bit-exact, and the package is
audit-first (it reports or refuses rather than silently transforming a
model). This vignette is the verb-for-verb, argument-for-argument map.

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
| `logLik`, `AIC`, `BIC`, `deviance`, `nobs`, `confint` | identical | `confint` supports Wald, profile (LMM), bootstrap (LMM) |
| `predict`, `fitted`, `residuals`, `simulate`, `refit` | identical | [`predict()`](https://rdrr.io/r/stats/predict.html) supports `re.form = NULL/NA`, `se.fit`, `interval` (population) |
| `update(fit, . ~ . - x)` | identical | formula edits, `REML=`, `weights=`, etc. |
| `anova(m1, m2)`, `drop1`, `getME`, `ngrps`, `isSingular` | identical | `isSingular()` is [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md) |
| `broom.mixed::tidy/glance/augment` | identical | registered for `mm_lmm`/`mm_glmm` |
| `emmeans::emmeans(fit, ~ x)` | identical | mixeff registers an emmeans basis |
| `lmerTest` p-values in [`summary()`](https://rdrr.io/r/base/summary.html) | built in | Satterthwaite/Kenward-Roger native, no extra package |

## Argument map for `lmm()` / `glmm()`

| lme4 argument | mixeff | Notes |
|----|----|----|
| `REML` | `lmm(..., REML=)` | same |
| `weights` | `weights=` | LMM and GLMM |
| `offset` | `glmm(..., offset=)` | GLMM only; LMM in-fit offset is not yet supported |
| `subset` | `lmm(..., subset=)` | supported for [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) |
| `na.action` | `lmm(..., na.action=)` | **default refuses NA**; pass `na.action = na.omit` for lme4’s complete-case behaviour |
| `contrasts` | partial | unordered factors use treatment coding, ordered factors `contr.poly` (both matching R/lme4 defaults); other codings are refused — recode the factor |
| `family = "binomial"` | `family = binomial()` | string families are not accepted |
| `nAGQ` | `glmm(..., nAGQ=)` | `>1` on the profiled path |
| `control = lmerControl(optimizer=, optCtrl=)` | `mm_control(optimizer=, max_feval=, ...)` | the engine picks a default optimizer; [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md) can override it or cap the evaluation budget |
| `start` | `mm_control(start=)` | theta warm starts |

## Four things that will bite, and the fix

**1. Coefficient names match lme4 exactly.** Since 0.2.0,
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`summary()`](https://rdrr.io/r/base/summary.html) tables,
[`vcov()`](https://rdrr.io/r/stats/vcov.html) dimnames, and
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
weight names use lme4’s naming and column order (`"recipeB"`,
`"temperature.L"`, `"recipeB:temperature.L"`), so name-keyed lme4 code
is drop-in compatible. (Earlier versions used an engine encoding like
`"recipe: B"`; if you wrote normalisation shims for those, delete them.)

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
contrast gets an independent variance). lme4’s `||` does **not** split
factor terms: a factor keeps its full within-factor covariance block. So
`(1 + cond + x || subj)` with a factor `cond` fits a strictly larger
model in lme4 than in mixeff, and the two disagree on the parameter
count (hence `df`, AIC, and — when the fitted within-factor covariance
is non-zero — the optimum itself). mixeff announces the situation at
compile time with an info diagnostic (`covariance_assumption`, reason
`double_bar_factor_term`). To reproduce lme4’s model family exactly,
write the expansion explicitly and give the factor its own correlated
block:

``` r

# mixeff `||`: independent variances for every column, factor levels included
glmm(y ~ cond * x + (1 + cond + x || subj), data, family = binomial())

# lme4-equivalent family: the factor keeps its within-factor covariance block
glmm(y ~ cond * x + (1 | subj) + (0 + cond | subj) + (0 + x | subj),
     data, family = binomial())
```

## What is `NA`-with-a-reason (and why)

`mixeff` never fabricates inference it cannot certify. Where lme4 would
silently return a number (or silently drop data), mixeff returns `NA`
with a machine- readable reason or raises a typed condition:

| Situation | lme4 | mixeff |
|----|----|----|
| `NA` in a model variable | silently dropped | refused unless `na.action = na.omit` |
| Boundary (singular) fit | one-time warning | persistent `[boundary]` tag + effective rank |
| Satterthwaite df at a boundary | may print unreliable df | refused with a reason; use bootstrap |
| Conditional prediction SE | not provided | `NA` with reason (population SE *is* provided) |
| GLMM `confint(method="profile")` | computed | refused (only Wald is certified for GLMMs) |

Use `inference_options(fit)` to see, before you run anything, which
inference routes are available on a given fit and why. \`\`\`
