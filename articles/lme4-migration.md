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
| `lmer(y ~ x + (x \| g), data)` | `lmm(y ~ x + (x \| g), data)` | same formula language, incl. `(x\|\|g)`, `(1\|g1/g2)`, crossed; fixed-part [`factor()`](https://rdrr.io/r/base/factor.html), [`poly()`](https://rdrr.io/r/stats/poly.html), [`scale()`](https://rdrr.io/r/base/scale.html), `ns()`/`bs()`, [`cut()`](https://rdrr.io/r/base/cut.html), [`relevel()`](https://rdrr.io/r/stats/relevel.html), [`offset()`](https://rdrr.io/r/stats/offset.html), `^`, `%in%` and `(a + b)^2` work with lme4’s coefficient names (evaluated in R and sent to the engine as columns); such transforms inside a random-effect bar are refused |
| `glmer(y ~ ..., family = binomial)` | `glmm(y ~ ..., family = binomial())` | pass a family **object** ([`binomial()`](https://rdrr.io/r/stats/family.html)), not a string |
| `lmerControl(...)` / `glmerControl(...)` | `mm_control(optimizer=, max_feval=, verbose=, start=, optCtrl=, disp_method=, disp_dof_correction=, max_phi_iter=)` | `optCtrl` names `maxfun`/`maxeval`, `ftol_*`, `xtol_*`, `rhobeg` are mapped (others refused); no `check.conv.*` (no post-hoc convergence warnings): use [`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md) / [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md); lme4 2.1-0’s GLMM dispersion controls map one to one: `disp_method = "moment"`/`"old/buggy"`, `disp_dof_correction`, and `maxPhiIter` -\> `max_phi_iter` (defaults `"moment"`, `TRUE`, `100`, as lme4 2.1-0); they only affect Gamma, inverse-Gaussian and Gaussian non-identity-link fits (elsewhere a `mm_control_ignored_notice` says they were ignored; lme4 is silent) and carry over to [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)/[`update()`](https://rdrr.io/r/stats/update.html) |
| `fixef`, `vcov` | same generics, same values | except for Gamma, inverse-Gaussian and Gaussian non-identity-link GLMMs where [`vcov()`](https://rdrr.io/r/stats/vcov.html) comes from RX (`nAGQ = 0` or the Hessian fallback): lme4 2.1-0 multiplies by `sigma()^2` although its weights already carry `1/phi`; mixeff reports the unscaled RX covariance (same for `predict(se.fit = TRUE)`) |
| `sigma` | identical | 1 for binomial, Poisson and negative-binomial GLMMs (the NB theta is `getME(fit, "glmer.nb.theta")`); Gamma, inverse-Gaussian and Gaussian non-identity-link GLMMs report `sqrt(phi)` with `phi = deviance / (n - rank([X, Z]))`, as lme4 \>= 2.1-0 (older lme4 used a different estimate, so those values, and their [`logLik()`](https://rdrr.io/r/stats/logLik.html), theta and [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md), match lme4 \>= 2.1-0 only) |
| `ranef` | same generic | `condVar = TRUE` for LMMs and GLMMs (GLMM conditional variances from the Laplace approximation at the fitted modes, as lme4; for Gamma, inverse-Gaussian and Gaussian non-identity-link GLMMs lme4 2.1-0 multiplies them by an extra `sigma()^2`, mixeff does not) |
| `coef` | identical | every fixed-effect column per group (fixed-only columns repeated), lme4’s column order |
| `VarCorr` | identical shape | named list of covariance matrices with `stddev`/`correlation` attributes and `sc`/`useSc`, so `VarCorr(m)$Subject` is the matrix; prints like lme4 (plus `[boundary]` markers); mixeff’s long table is still `VarCorr(m)$table` (and `$residual_sd`); [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) matches lme4; GLMMs have `useSc = FALSE` (no Residual row, absolute SDs) with `sc = sigma()`, as lme4 \>= 2.1-0 |
| `logLik`, `AIC`, `BIC`, `nobs` | identical | `AIC(m1, m2)` / `BIC(m1, m2)` return stats’ `df`/`AIC` data frame |
| `deviance` | identical | LMM: as lme4 for ML fits (REML fits return the REML criterion without lme4’s deprecation warning); GLMM: sum of squared deviance residuals, as lme4 |
| `REMLcrit` | [`REMLcrit()`](https://bbuchsbaum.github.io/mixeff/reference/REMLcrit.md) | exported by mixeff (lme4’s is not a generic; whichever package is attached last masks the other, and mixeff’s forwards lme4 fits to lme4) |
| `family`, `weights`, `model.frame` | identical | [`weights()`](https://rdrr.io/r/stats/weights.html) returns ones for unweighted fits and `type = "working"` for GLMMs; [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) has lme4’s transformed columns and `terms` attribute (the raw variables the engine used are `fit$model_frame`) |
| `confint` | same generic, **different default** | Wald (`"asymptotic"`, lme4 spelling `"Wald"`) by default, where lme4 defaults to profile; returns fixed effects only. `"profile"` (LMMs) returns lme4’s rows (`.sig01`, …, `.sigma`, fixed effects; REML fits profiled on the ML deviance, as lme4), and accepts `threads =`; `"bootstrap"` (lme4 spelling `"boot"`); GLMM availability per fit type is the route matrix in [`vignette("inference", package = "mixeff")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md) |
| `predict`, `fitted` | identical | `re.form = NULL/NA/~0` or a partial formula such as `~(1 \| g)` (no `se.fit`/`interval` for partial ones), `allow.new.levels`, `na.action` (default `na.pass`), `se.fit`, `interval`; `newdata` is narrowed to the needed variables and coerced like the training frame; an `offset =` fit offset must be resupplied with `predict(offset =)` (lme4 silently uses 0) |
| `residuals` | identical | LMM default `"response"`, `"pearson"`/`"deviance"` are `sqrt(w) * (y - mu)`; GLMM default `"deviance"`, plus `"pearson"`, `"working"`, `"response"`; `scaled = TRUE` divides by [`sigma()`](https://rdrr.io/r/stats/sigma.html) once |
| `simulate`, `refit` | LMM and GLMM | lme4’s `re.form`/`use.u` meaning (default `NA`: new random effects; `NULL`: condition on the BLUPs) and RNG order, so seeded draws match lme4; LMM draws honour prior weights (lme4 ignores them); Gamma (shape `w/phi`), Gaussian non-identity-link and inverse-Gaussian fits draw random effects on the absolute scale reported by lme4 \>= 2.1-0, and Gamma and unweighted Gaussian draws match lme4 2.1-0’s; Gaussian draws honour prior weights (lme4 ignores them); inverse-Gaussian draws use shape `w/phi` (needs statmod), where lme4 2.1-0 uses `w/sigma` and so simulates variance `sigma*mu^3` instead of the fitted `phi*mu^3`; a factor binomial response simulates as 0/1; a partial `re.form` such as `~(1 \| g)` conditions on those terms and redraws the others, as lme4; `newdata` and `newparams` are refused |
| `plot(fit)`, `dotplot(ranef(fit, condVar = TRUE))`, `qqmath(ranef(...))` | base-graphics [`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html); lattice methods | `plot(fit)` draws Pearson residuals vs fitted with base graphics (no formula interface); `plot(ranef(fit, condVar = TRUE))` gives a caterpillar plot |
| `hatvalues`, `cooks.distance`, `influence` | same generics | LMM values match lme4; GLMM hat values use the IRLS working weights (as `glm` does), so they differ from lme4’s; [`influence()`](https://rdrr.io/r/stats/lm.influence.html) takes `groups=` and `ncores=` and records failed refits instead of stopping |
| [`performance::r2()`](https://easystats.github.io/performance/reference/r2.html), [`performance::icc()`](https://easystats.github.io/performance/reference/icc.html), [`insight::get_variance()`](https://easystats.github.io/insight/reference/get_variance.html) | same calls, or [`mm_r2()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md), [`mm_icc()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md), [`mm_variance_components()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md) | Nakagawa R2/ICC agree with performance on lme4 fits; `by_group = TRUE` is not supported |
| `update(fit, . ~ . - x)` | same generic | a formula that needs new columns re-evaluates the original `data`, as lme4 does; arguments [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) cannot carry over (e.g. `subset`, `contrasts`) are refused rather than dropped; the fit’s `na.action` carries over |
| `anova(m1, m2)`, `drop1` | LRT with automatic ML refit | `anova(m1, m2)` is lme4’s data frame (`npar AIC BIC logLik -2*log(L) Chisq Df Pr(>Chisq)`, as lme4 \>= 2.0); single-model [`anova()`](https://rdrr.io/r/stats/anova.html) is lmerTest’s (`Sum Sq Mean Sq NumDF DenDF F value Pr(>F)`, Type I/II/III, Satterthwaite/KR) and LMM-only; mixeff’s provenance table is `$table`; a multi-column term such as `poly(x, 2)` is tested as one multi-df row, like lmerTest, unless another term contains it (then per column) |
| `getME` | [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md) | LMM and GLMM; lme4’s components including `u`, `b`, `L`, `RX`, `RZX`, `devcomp`, `lower`, `Gp`, `Ztlist`, `Tlist`, `ST`, `glmer.nb.theta`, `"ALL"`; `devfun` is refused (the optimizer lives in the engine) |
| `isSingular` | `is_singular(fit, tol = 1e-4)` | lme4’s rule (a zero-bounded theta below `tol`), LMM and GLMM; [`lme4::isSingular()`](https://rdrr.io/pkg/lme4/man/isSingular.html) is not a generic, so call [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md) |
| `ngrps` | identical | interaction/nested groups are named as in lme4 (`a:b`; `b:a` for the nested term of `a/b`) with levels `x:y` |
| `broom.mixed::tidy/glance/augment` | supported | registered for `mm_lmm`/`mm_glmm` |
| `emmeans::emmeans(fit, ~ x)` | supported | estimates and SEs agree with lme4’s; mixeff’s status/reliability/reason audit columns do not survive the emmeans trip, and a profiled (`pirls_profiled`) GLMM refuses inference through emmeans |
| `lmerTest` p-values in [`summary()`](https://rdrr.io/r/base/summary.html) | built in | Satterthwaite/Kenward-Roger native, no extra package; `summary()$coefficients` is lmerTest’s matrix (per-row method labels in `summary()$coef_table`); lmerTest’s `ddf = "Satterthwaite"`/`"Kenward-Roger"`/`"lme4"` is accepted by [`summary()`](https://rdrr.io/r/base/summary.html) and [`anova()`](https://rdrr.io/r/stats/anova.html) |
| `bootMer(m, FUN, nsim)` | `mm_bootmer(fit, FUN, nsim)` | LMM; a `boot`-compatible result for [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html); `type = "semiparametric"` refused |
| `rePCA`, `lmList`, `allFit` | [`rePCA()`](https://bbuchsbaum.github.io/mixeff/reference/rePCA.md), [`mm_lmlist()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmlist.md), [`mm_allfit()`](https://bbuchsbaum.github.io/mixeff/reference/mm_allfit.md) | `mm_*` names avoid masking lme4’s non-generic functions |
| `pbkrtest::KRmodcomp(big, small)` | `compare(small, big, method = "kenward_roger")` | also `"satterthwaite"` (`SATmodcomp`); KR F is pbkrtest’s scaled `Ftest` row (unscaled `FtestU` in `$fixed_f$unscaled_statistic`) |
| `lmerTest::step(m)` | `step(fit)` | random terms tested with the boundary-corrected LRT; whole terms only |
| `nlmer`, `devFunOnly = TRUE` | not provided | nonlinear mixed models are a non-goal; the engine owns the objective |

## Argument map for `lmm()` / `glmm()`

| lme4 argument | mixeff | Notes |
|----|----|----|
| `REML` | `lmm(..., REML=)` | same |
| `weights` | `weights=` | a column of `data` or a vector, as in lme4 (LMM and GLMM) |
| `offset` | `offset=` | a column of `data`, an expression, or a vector (LMM and GLMM); [`offset()`](https://rdrr.io/r/stats/offset.html) formula terms are added to it. LMMs fit `y - offset` and add it back in [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`predict()`](https://rdrr.io/r/stats/predict.html) |
| `subset` | `subset=` | LMM and GLMM; formula transforms are evaluated before subsetting, as in [`model.frame()`](https://rdrr.io/r/stats/model.frame.html); unused factor levels are dropped afterwards |
| `na.action` | `na.action=` | LMM and GLMM; default `getOption("na.action")` (`na.omit`) as in lme4; only model variables, `weights` and `offset` count; dropped rows are announced (`mm_rows_dropped` message) and recorded in `na.action(fit)`; `na.exclude` pads [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)/[`predict()`](https://rdrr.io/r/stats/predict.html); `na.fail`/`na.pass` raise `mm_data_error` |
| `contrasts` | partial | unordered factors use treatment coding, ordered factors `contr.poly` (both matching R/lme4 defaults). Any other coding is refused with an `mm_arg_error` rather than silently replaced, whether it is requested through `contrasts=`, a contrast attached to the factor (`contrasts(f) <- ...`), or `options(contrasts=)`; use [`relevel()`](https://rdrr.io/r/stats/relevel.html) to change the reference level, or build coded columns as numeric predictors |
| character predictors | converted to factors | sorted levels, as lme4’s [`factor()`](https://rdrr.io/r/base/factor.html) conversion does; logical predictors are coded as factors with levels `FALSE`/`TRUE`, as [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) does (e.g. `lFALSE`, `lTRUE` without an intercept) |
| numeric/character grouping (accepted silently) | coerced to factor | non-factor grouping columns are coerced to factor, announced via an `mm_grouping_coercion_notice`; silence it with `mm_control(verbose = -1)`; an ordered grouping factor is used as an unordered one (grouping needs only its levels) |
| `family = "binomial"` | `family = binomial()` | string families are not accepted |
| `nAGQ` | `glmm(..., nAGQ=)` | `1` (default) is joint Laplace, as in `glmer`; `nAGQ = 0` selects the profiled fast path (`method = "pirls_profiled"`, lme4’s PIRLS-only estimate); `>1` runs on the profiled path only and, as in lme4, needs a single scalar random effect (selected with an `mm_estimator_notice` when `method` is not given) |
| `control = lmerControl(optimizer=, optCtrl=)` | `mm_control(optimizer=, max_feval=, optCtrl=, ...)` | the engine picks a default optimizer; [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md) can override it or cap the evaluation budget |
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

**3. The GLMM estimator.**
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
defaults to `method = "joint_laplace"`, glmer’s joint Laplace
(`nAGQ = 1`) estimator, with certified Wald inference. The fast profiled
(PIRLS) estimator — lme4’s `nAGQ = 0` — is opt-in; its coefficients do
**not** match `glmer(nAGQ = 1)` exactly. Negative-binomial families (no
joint route yet), `nAGQ > 1`, and `inference = "working_hessian"` need
the profiled path: without an explicit `method`,
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
selects it and prints an `mm_estimator_notice` saying so.

``` r

glmm(y ~ x + (1 | g), data, family = binomial())   # joint Laplace (glmer)
glmm(y ~ x + (1 | g), data, family = binomial(),
     method = "pirls_profiled")                      # fast path (nAGQ = 0)
```

**4. `||` expands exactly like lme4’s, factors included.** mixeff
expands `(1 + x + f || g)` the way lme4’s `expandDoubleVerts()` does:
`(1 | g) + (0 + x | g) + (0 + f | g)` — the numeric columns get
independent variances, while the factor `f` keeps its full within-factor
covariance block, coded like R’s
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) (so
`(0 + f + h | g)` is full rank, as in lme4). The parameter count, `df`,
AIC and the optimum therefore match lme4, and
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
lists the pieces as separate entries named `g`, `g.1`, … exactly as lme4
does (also in `as.data.frame(VarCorr(fit))`). To fix *every* covariance
at zero — including those among a factor’s columns (MixedModels.jl’s
`zerocorr()`) — write `diag(1 + f | g)`, which has no lme4 equivalent.
mixeff still notes a factor inside `||` at compile time with an info
diagnostic (`covariance_assumption`, reason `double_bar_factor_term`).

``` r

# lme4 semantics: independent intercept and x slope; the factor keeps its block
glmm(y ~ cond * x + (1 + cond + x || subj), data, family = binomial())
# same model, written out
glmm(y ~ cond * x + (1 | subj) + (0 + x | subj) + (0 + cond | subj),
     data, family = binomial())

# MixedModels.jl zerocorr semantics: every covariance fixed at zero
glmm(y ~ cond * x + diag(1 + cond + x | subj), data, family = binomial())
```

Comparing a `||` model with its correlated `|` counterpart adds only
correlation parameters; the likelihood-ratio test (and
[`test_random_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_random_effect.md)’s
boundary route) then uses the ordinary chi-square, as
[`anova()`](https://rdrr.io/r/stats/anova.html) on lme4 fits does.

**5. Random-effect terms are presented in lme4’s order.**
[`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
[`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
(`theta`, `Lambdat`, `Zt`, `cnms`, …) and the `.sig01`, … rows of
`confint(method = "profile")` use lme4’s `mkReTrms` term order (formula
order, reordered by decreasing number of levels only when a later term
has more levels than an earlier one). The engine fits the same model
with its own internal (MixedModels.jl) order; the translation is exact.

**6. Profile intervals use lme4’s scale and names, but Wald is the
default.** `confint(fit, method = "profile")` returns lme4’s rows
(`.sig01`, …, `.sigma`, then the fixed effects) on the
standard-deviation/correlation scale, numbered as lme4 \>= 2.0 numbers
them (term by term, each term’s standard deviations before its
correlations; lme4 1.1 interleaved them in lower-triangle order), and,
like lme4, profiles a REML fit on the ML deviance. lme4’s
[`confint()`](https://rdrr.io/r/stats/confint.html) defaults to the
profile; mixeff’s default stays the Wald interval (the profile refits
the model many times), so pass `method = "profile"` to reproduce lme4’s
default output. Bootstrap and profile computations accept `threads =`
(via `bootstrap_control(threads = )` for bootstraps); results are
identical for every thread count, and the default of 1 respects CRAN’s
two-thread policy.

**7. Rank-deficient fixed effects.** As in lme4, of two collinear
columns the earlier is kept;
[`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
omits the dropped coefficient (use `fixef(fit, add.dropped = TRUE)` for
an `NA` in its place), and [`vcov()`](https://rdrr.io/r/stats/vcov.html)
and [`summary()`](https://rdrr.io/r/base/summary.html) cover the
estimable coefficients only.

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
| Gamma | inverse, log | pirls_profiled, joint_laplace | inverse is R’s default Gamma link. |
| inverse.gaussian | inverse, log | pirls_profiled, joint_laplace | R’s default 1/mu^2 link is not available in the engine and is refused; pass link = “inverse” or “log”. |
| gaussian | log, inverse, sqrt | pirls_profiled, joint_laplace | Non-identity links only; gaussian/identity is an LMM and is refused with a pointer to lmm(). |
| negative_binomial | log | pirls_profiled | NB2 via mm_negative_binomial(); theta estimated or fixed. Intervals via confint(method = “bootstrap”). |

Supported GLMM family/link/estimator cells. {.table}

| feature | lmm | glmm | notes |
|:---|:---|:---|:---|
| case weights | supported | supported | GLMM binomial weights double as trial counts for proportion responses. |
| offset | supported | supported | offset = argument and offset() formula terms. LMMs fit y - offset and add the offset back to fitted values and predictions (lmer semantics). predict(newdata) evaluates formula offsets from newdata; an offset = argument offset must be resupplied via predict(offset =). |
| simulate | supported | supported | lme4 re.form/use.u semantics (default NA: new random effects; NULL: conditional on the BLUPs). GLMM draws are family-aware (binomial, Poisson, Gamma, negative binomial); partial re.form formulas, newdata and newparams are refused. |
| refit | supported | supported | Re-fits with the stored settings; binomial fits accept a (successes, failures) matrix, a two-level factor or 0/1 values. |
| subset / na.action / contrasts | supported | supported | Same data preparation for lmm() and glmm(): na.action defaults to getOption(“na.action”) (na.omit) as in lme4, applied to the model variables, weights and offset only; dropped rows are announced (typed mm_rows_dropped message) and recorded in na.action(fit); na.exclude also pads fitted/residuals/predict; na.fail and na.pass raise mm_data_error. Unused factor levels are dropped after subset/NA removal. contrasts = is honoured only when it names the engine coding (contr.treatment / contr.poly). |
| random = (GLMM) | not applicable | refused | Reserved argument on glmm(); write random terms in the formula. |
| stateful fixed-effect terms (factor(), poly(), scale(), ns()/bs(), cut(), relevel(), offset(), ^, %in%, parenthesised groups) | supported | supported | Evaluated in R as stats::model.frame() does and sent to the engine as synthetic columns; coefficient names match lme4 and predict(newdata) reuses the training basis. Such terms inside random-effect bars are refused by the engine. |
| marginal means verbs (mm_means / mm_comparisons / mm_grid / mm_predictions / test_effect) | supported | refused | GLMM marginal means go through the emmeans bridge on certified fits. |
| double-bar \|\| with factor terms | supported | supported | Expanded like lme4: (1 + x + f \|\| g) is (1 \| g) + (0 + x \| g) + (0 + f \| g), the factor keeping its own block (model.matrix() coding); VarCorr() lists g, g.1, … as lme4. diag(1 + f \| g) fixes every covariance at zero (MixedModels.jl zerocorr()). Noted at compile time (covariance_assumption / double_bar_factor_term). |

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
| `NA` in a model variable | silently dropped (`na.omit`) | dropped under the same `na.action`, announced with a typed `mm_rows_dropped` message and recorded in `na.action(fit)` |
| Boundary (singular) fit | one-time warning | persistent `[boundary]` tag + effective rank |
| Satterthwaite df at a boundary | lmerTest prints df anyway | refused with a reason; use bootstrap |
| Prediction SE for an unseen grouping level | — | `NA` with the reason in `mm_reason` |
| GLMM [`confint()`](https://rdrr.io/r/stats/confint.html) | computed | Wald intervals on `joint_laplace` fits via `method = "asymptotic"` (lowercase; lme4’s `method = "Wald"` spelling is not accepted); parametric-bootstrap intervals via `method = "bootstrap"` on profiled-estimator fits (like `bootMer`; refused on `joint_laplace` fits at this engine pin); Wald routes refused on `pirls_profiled` fits |

Use `inference_options(fit)` to see, before you run anything, which
inference routes are available on a given fit and why.
