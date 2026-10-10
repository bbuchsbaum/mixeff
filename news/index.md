# Changelog

## mixeff 0.2.0

### Bridge performance: compile once, keep the fitted model (mixeff-rs bridge-perf)

- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  compile the formula and audit the design once per fit. The engine
  builds the model from the spec compiled for the pre-fit explanation
  (`CompiledModelSpec`), so the second compile/audit pass and the second
  translation of the data to the engine are gone. The pre-fit
  explanation now shows the artifact the model is actually built from:
  for a formula with in-formula transforms the engine evaluates (such as
  `log(x)`), or a `||` term on a factor, it reports the evaluated and
  expanded design.
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  is unchanged.
- Fits now keep the fitted engine model as a native handle
  (`fit$rust_handle`, see
  [`fit_handle_alive()`](https://bbuchsbaum.github.io/mixeff/reference/fit_handle_alive.md)).
  Computations that need the engine model reuse it instead of rebuilding
  and refitting the model from the stored data:
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`anova()`](https://rdrr.io/r/stats/anova.html) tests (Satterthwaite
  and Kenward-Roger),
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
  [`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md),
  [`test_random_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_random_effect.md),
  `confint(method = "profile")`, bootstrap contrasts and intervals,
  [`parametric_bootstrap()`](https://bbuchsbaum.github.io/mixeff/reference/parametric_bootstrap.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html) on new data and
  with intervals or `se.fit` (LMMs and GLMMs), `ranef(condVar = TRUE)`,
  GLMM bootstrap intervals, and
  [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md).
  Results are identical to the cold refit. Profiles and convergence
  checks work on a copy, so the cached model never changes. These calls
  stay interruptible. On a 200,000-row crossed LMM,
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  drops from 0.39 s to 0.002 s,
  `predict(newdata, interval = "prediction")` from 0.70 s to 0.04 s and
  [`summary()`](https://rdrr.io/r/base/summary.html) by the 0.45 s
  refit. On a 3,000-row binomial GLMM, `predict(newdata, se.fit = TRUE)`
  drops from 0.44 s to 0.006 s.
- The handle is a process-local cache. It does not survive
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) /
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) or a new R process.
  Those fits, and fits where the handle does not match the requested
  refit (for example a profile under ML of a REML fit), fall back to the
  cold refit automatically, with identical results.
  [`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
  keeps a live handle and clears a dead one. Copies of a fit share the
  handle.
- Memory: the handle keeps the engine model, including a copy of the
  model data, alive for as long as the fit object, outside R’s heap. To
  fit without keeping handles, set
  `options(mixeff.keep_handle = FALSE)`.
  [`update()`](https://rdrr.io/r/stats/update.html) and
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)
  create new fits with their own handles.

### Profile memory (mixeff-rs [\#12](https://github.com/bbuchsbaum/mixeff/issues/12))

- `confint(method = "profile")` and
  [`profile()`](https://rdrr.io/r/stats/profile.html) no longer build a
  dense n x n covariance for the fixed-effect profiles. A crossed LMM
  with 60,000 rows previously aborted R trying to allocate 28.8 GB; it
  now profiles in seconds with memory linear in the data. Leverage,
  Cook’s distance and new-data prediction variance also stop densifying
  the factor’s off-diagonal blocks once per row.

### Engine follow-ups (mixeff-rs engine-followups)

- GLMMs with a free dispersion parameter (Gamma, inverse Gaussian, and
  Gaussian with a non-identity link) now follow lme4 2.1-0’s estimated
  dispersion handling, which the engine adopts. The working weights
  carry `1/phi`, and `phi` is profiled during the fit.
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) is `sqrt(phi)` with
  `phi = deviance / (n - rank([X, Z]))`, and theta is the absolute
  random-effect SD. [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html),
  [`sigma()`](https://rdrr.io/r/stats/sigma.html), theta and
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  equal `glmer()` from lme4 \>= 2.1-0. Fits of these families no longer
  match older lme4, and their estimates change from earlier mixeff
  versions. mixeff no longer shifts the logLik of Gamma and
  inverse-Gaussian fits by 1; the engine’s value already matches lme4
  2.1-0.

- For these families `weights(fit, "working")` is divided by `phi`, as
  in lme4 2.1-0. This also changes
  [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)’s
  `L`, `RZX`, `RX` and `devcomp`. GLMM `devcomp` dims gain lme4 2.1-0’s
  `dispProfile`, `maxPhiIter` and `qEff`.
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  of every GLMM has `useSc = FALSE`, with no Residual row and absolute
  SDs, and keeps `sc = sigma()`.
  [`rePCA()`](https://bbuchsbaum.github.io/mixeff/reference/rePCA.md) no
  longer rescales GLMM covariances by
  [`sigma()`](https://rdrr.io/r/stats/sigma.html). Hat values and Cook’s
  distances use the new weights; the values are unchanged.

- Known lme4 2.1-0 difference: for these families lme4 2.1-0 multiplies
  by `sigma()^2` a second time in three places: the RX-based
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) (`nAGQ = 0`),
  `predict(se.fit = TRUE)` and `ranef(condVar = TRUE)`. mixeff reports
  the engine’s unscaled covariance there.

- The Nakagawa distribution-specific variance of a Gamma log-link GLMM
  is now insight 1.5’s log-normal approximation `log1p(sigma^2)`, not
  `sigma^2`. [`family()`](https://rdrr.io/r/stats/family.html) now works
  for Gaussian non-identity-link GLMMs.

- [`mm_r2()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md),
  [`mm_icc()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
  and
  [`mm_variance_components()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
  now support Gaussian GLMMs with a non-identity link. As insight 1.5 /
  performance 0.18 do for `glmer()` fits, the distribution-specific
  variance is `sigma()^2`, and observation-level terms add no
  dispersion. This residual variance is on the response scale, so the R2
  and ICC depend on the units of the response. Inverse-Gaussian GLMMs
  are still refused, now with reason code
  `"r2_distribution_variance_undefined"`. insight 1.5 has no
  inverse-Gaussian case and falls back to
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) itself, which is not a
  variance and changes when the response is rescaled.

- [`simulate()`](https://rdrr.io/r/stats/simulate.html) (and so
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)-based
  bootstraps) now works for inverse Gaussian and Gaussian
  non-identity-link GLMMs. Previously these were refused with
  `"simulate_family_unavailable"`. Draws follow the fitted variance
  function with `phi = sigma()^2`, and new random effects are drawn on
  lme4 2.1-0’s absolute scale. Gamma draws and unweighted Gaussian draws
  match lme4 2.1-0’s
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draw for draw.
  The tests compare seeded draws with live `glmer()` or with stored lme4
  2.1-0 draws. Gaussian draws use prior weights
  (`sd = sigma / sqrt(w)`), which lme4 ignores. Inverse-Gaussian draws
  use
  [`statmod::rinvgauss()`](https://rdrr.io/pkg/statmod/man/invgauss.html)
  (statmod is now in Suggests) with shape `w / phi`. lme4 2.1-0 instead
  uses `w / sigma()`, so its draws have variance `sigma * mu^3` rather
  than the fitted `phi * mu^3`.

- lme4 2.1-0’s `glmerControl()` dispersion settings are now available as
  `mm_control(disp_method = c("moment", "old/buggy"), disp_dof_correction, max_phi_iter)`.
  lme4’s `maxPhiIter` becomes `max_phi_iter`. They are sent to the
  engine (`set_dispersion_method()`, `set_dispersion_dof_correction()`,
  `set_max_phi_iter()`). Invalid values raise `mm_arg_error`, and unset
  ones keep lme4 2.1-0’s defaults (`"moment"`, `TRUE`, `100`). They
  affect Gamma, inverse-Gaussian and Gaussian non-identity-link GLMMs.
  For other families and for
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) they
  have no effect, and the fit says so with an
  `mm_control_ignored_notice` (lme4 ignores them silently). Because the
  control is stored on the fit,
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md),
  [`update()`](https://rdrr.io/r/stats/update.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html)-based refits,
  bootstraps and the R2 null model reuse it. `getME(fit, "devcomp")`
  dims `dispProfile`, `maxPhiIter` and `qEff` reflect the settings.
  Under `"old/buggy"`, the working weights do not carry `1/phi`, and
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  / [`simulate()`](https://rdrr.io/r/stats/simulate.html) use theta
  itself as the random-effect SD, as lme4 2.1-0 does. Fits match
  `glmer(control = glmerControl(...))` from lme4 2.1-0. Under
  `"old/buggy"` they agree to about 2e-4 in the estimates and 3e-6
  relative in [`logLik()`](https://rdrr.io/r/stats/logLik.html), at the
  same lme4 deviance. The tests compare with stored lme4 2.1-0 values
  when the installed lme4 is older.

- `confint(method = "profile")` and
  [`profile()`](https://rdrr.io/r/stats/profile.html) no longer fail
  when one parameter’s profile is irregular. The row is kept with a
  `status` (`"non_monotone"`, `"not_bracketing"` or `"failed"`), a
  `reason_code` (`"non_monotone_profile"`, `"profile_not_bracketing"`,
  `"profile_failed"`) and the engine’s `reason` in
  `attr(ci, "mm_profile")$table`. Bounds the profile cannot determine
  are `NA`, and printing the intervals adds a note. Failed `.sigNN` rows
  now appear; previously they were dropped. lme4 instead warns and falls
  back to linear interpolation, which often gives `[-1, 1]` for a
  correlation. The profile CI JSON schema accepts the new fields, null
  bounds and the engine’s string notes.

- The multi-df Kenward-Roger F is now pbkrtest’s scaled `Ftest`
  statistic (`F = lambda * F_U` with KR denominator df). Previously it
  was the unscaled `FtestU`. This affects
  `compare(method = "kenward_roger")` and
  `anova(ddf = "Kenward-Roger")`, which matches lmerTest: the F value is
  scaled and `Mean Sq` comes from the unscaled F. `$fixed_f` gains
  `statistic_scale`, `unscaled_statistic` and `unscaled_p_value`;
  `f_scaling` is lambda.

- When a joint-Laplace GLMM’s Hessian is not positive definite, mixeff
  now reports Wald inference instead of withholding it, as glmer does.
  The fixed-effect covariance falls back to RX conditional on theta, as
  in glmer’s `vcov(use.hessian = FALSE)`, with reliability `"low"`.
  [`summary()`](https://rdrr.io/r/base/summary.html) explains the
  fallback and prints the engine’s notes,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) warns (class
  `mm_vcov_rx_fallback`), and
  [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
  labels the route `glmm_laplace_rx_conditional_on_theta_wald`. These
  standard errors ignore uncertainty in theta.

- A joint-Laplace fit that the engine labels `not_optimized` is now
  shown as “convergence not certified” when the optimizer stopped
  normally within its budget and the estimated objective gap is at most
  `1e-2` deviance units. The engine certifies only at `1e-6`. This is a
  warning-level label, like lme4’s convergence warnings, and does not
  refuse the fit.
  [`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
  still returns the engine’s label.

- [`inference_options()`](https://bbuchsbaum.github.io/mixeff/reference/inference_options.md)
  on joint-Laplace GLMMs no longer describes withheld Wald rows as
  “uncertified for the profiled estimator”. It also reports the
  parametric bootstrap as refused, matching
  [`confint()`](https://rdrr.io/r/stats/confint.html).

### Correctness fixes (pre-release audit)

- Character fixed-effect predictors (e.g. columns from
  [`read.csv()`](https://rdrr.io/r/utils/read.table.html)) are now
  converted to factors with sorted levels, as lme4 does. Previously the
  fit aborted with an `mm_schema_error` while building coefficient
  names.
- A numeric predictor whose name starts with a factor’s name (`group`
  and `group_size`) no longer aborts the fit.
- Contrast codings mixeff cannot honour are refused with an
  `mm_arg_error` instead of being silently replaced by treatment coding:
  a contrast attached to an unordered factor, a non-`contr.treatment`
  unordered entry in `options(contrasts =)`, and
  `contrasts = list(f = "contr.SAS")` (which uses the last level as
  reference; the engine always uses the first).
- [`update()`](https://rdrr.io/r/stats/update.html) re-evaluates the
  original `data` when the updated formula needs columns the stored
  model frame lacks (as lme4 does), and refuses arguments it cannot
  carry over (e.g. `subset`, `contrasts`) instead of dropping them. It
  carries the fit’s `na.action` over (and accepts a new one).
- The declared Rust toolchain minimum is now 1.85, the bundled engine’s
  actual requirement (some vendored crates use edition 2024).
- The lme4 migration vignette no longer describes
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`anova()`](https://rdrr.io/r/stats/anova.html),
  [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md),
  `isSingular()` and `emmeans` as identical to lme4; it lists the actual
  differences.

### lme4 parity and checklist completion

- Missing values are handled as in `lmer()`/`glmer()`:
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  honour `na.action`, whose default is now `getOption("na.action")`
  (`na.omit` unless changed) instead of refusing `NA`. Only the
  variables the model uses count (fixed and random terms, evaluated
  transforms, `weights`, `offset`), so an `NA` in an unused column never
  drops a row. Dropped rows are announced with a typed `mm_rows_dropped`
  message giving the number of rows and the variables with missing
  values (silence with `mm_control(verbose = -1)`), recorded as
  `na.action(fit)` and `attr(model.frame(fit), "na.action")`, and
  excluded from [`nobs()`](https://rdrr.io/r/stats/nobs.html).
  `na.exclude` pads
  [`residuals()`](https://rdrr.io/r/stats/residuals.html)/[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`predict()`](https://rdrr.io/r/stats/predict.html)
  for both fitters; `na.fail` and `na.pass` are refused with a typed
  `mm_data_error`.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md),
  [`update()`](https://rdrr.io/r/stats/update.html), bootstrap and
  influence measures use the fitted rows, and
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)/[`anova()`](https://rdrr.io/r/stats/anova.html)
  refuse models fitted to different numbers of rows
  (`reason_code = "different_nobs"`), as lme4’s
  [`anova()`](https://rdrr.io/r/stats/anova.html) does.

- Shapes follow lme4 \>= 2.0 and insight \>= 1.5:

  - `confint(method = "profile")` numbers `.sigNN` term by term, each
    term’s standard deviations before its correlations. lme4 1.1
    numbered them in lower-triangle order.
  - `anova(m1, m2)` names its `-2 * logLik` column `-2*log(L)`, where it
    was `deviance`.
  - `getME(m, "devcomp")$dims` gains `npar`.
  - `print(VarCorr(m))` uses
    [`reformulas::formatVC()`](https://rdrr.io/pkg/reformulas/man/formatVC.html)
    when reformulas is installed.
  - For GLMMs,
    [`mm_r2()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
    /
    [`mm_icc()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
    divide the binomial distribution-specific variance of a
    [`cbind()`](https://rdrr.io/r/base/cbind.html) response by the mean
    number of trials. For count models, the null-model mean is
    `exp(b0 + v0/2)`, where `b0` is the null model’s intercept and `v0`
    its random-effect variance.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html) follows lme4:
  LMM `type = "pearson"`/`"deviance"` are `sqrt(w) * (y - mu)` (no
  longer divided by sigma), and `scaled = TRUE` divides by
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) exactly once (it
  divided twice for Pearson). GLMM residuals default to `"deviance"` and
  support `"pearson"`, `"working"` and `"response"`, computed from the
  family, prior weights and binomial trials.

- [`sigma()`](https://rdrr.io/r/stats/sigma.html) of a negative-binomial
  GLMM is 1, as in lme4; theta is `getME(fit, "glmer.nb.theta")`. GLMM
  [`deviance()`](https://rdrr.io/r/stats/deviance.html) is the sum of
  squared deviance residuals (as lme4), no longer `-2 * logLik()`.

- [`coef()`](https://rdrr.io/r/stats/coef.html) returns every
  fixed-effect column per group in lme4’s order.

- [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  returns lme4’s structure (named list of covariance matrices with
  `stddev`/`correlation` attributes, `sc`/`useSc`), so
  `VarCorr(m)$Subject` is the matrix; it prints like lme4.
  `VarCorr(m)$table` and `$residual_sd` still work.

- `AIC(m1, m2)` / `BIC(m1, m2)` return stats’ `df`/`AIC` data frame
  instead of refusing.

- New
  [`REMLcrit()`](https://bbuchsbaum.github.io/mixeff/reference/REMLcrit.md),
  [`family()`](https://rdrr.io/r/stats/family.html) methods,
  `weights(type = "working")`;
  [`weights()`](https://rdrr.io/r/stats/weights.html) returns ones for
  unweighted fits;
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html) has lme4’s
  transformed columns and `terms` attribute (raw variables stay in
  `fit$model_frame`).

- [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
  covers lme4’s components (`u`, `b`, `L`, `RX`, `RZX`, `devcomp`,
  `lower`, `Gp`, `Tp`, `Lind`, `Ztlist`, `Tlist`, `ST`, `offset`,
  `weights`, `glmer.nb.theta`, `"ALL"`, …) for LMMs and GLMMs; `devfun`
  is refused with a typed error. `Lambda`/`Lambdat` (and
  `model.matrix(fit, type = "random")`) now follow the engine’s (lme4’s)
  term order; previously theta was assigned in formula order, so crossed
  designs such as Penicillin’s `(1|sample) + (1|plate)` gave a term
  another term’s theta.

- [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
  uses lme4’s rule (a zero-bounded theta below `tol`), honours `tol`,
  and works for GLMMs.

- `ranef(condVar = TRUE)` returns Laplace conditional variances for
  GLMMs instead of `NA`.

- [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`anova()`](https://rdrr.io/r/stats/anova.html) accept lmerTest’s
  `ddf = "Satterthwaite"`, `"Kenward-Roger"` or `"lme4"`.
  `summary()$coefficients` is lmerTest’s numeric matrix (per-row method
  labels moved to `summary()$coef_table`). `anova(m)` is lmerTest’s data
  frame (`Sum Sq Mean Sq NumDF DenDF F value Pr(>F)`) and
  `anova(m1, m2)` lme4’s
  (`npar AIC BIC logLik deviance Chisq Df Pr(>Chisq)`); mixeff’s
  provenance table stays in `$table`.

- [`confint()`](https://rdrr.io/r/stats/confint.html) accepts lme4’s
  `method = "Wald"` and `"boot"` spellings.

- Stateful and expanded fixed-effect formula terms now work in
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md):
  [`factor()`](https://rdrr.io/r/base/factor.html),
  [`relevel()`](https://rdrr.io/r/stats/relevel.html),
  [`cut()`](https://rdrr.io/r/base/cut.html),
  [`scale()`](https://rdrr.io/r/base/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html)/`bs()`,
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html), `I(x > 0)`,
  [`offset()`](https://rdrr.io/r/stats/offset.html), `^`, `%in%`, and
  parenthesised groups such as `(a + b)^2` and `a * (b + c)`.
  Coefficient names, [`terms()`](https://rdrr.io/r/stats/terms.html),
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  [`drop1()`](https://rdrr.io/r/stats/add1.html) and `emmeans` grids
  follow lme4, and `predict(newdata)` reuses the training basis.

- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) gains
  `offset =` (and [`offset()`](https://rdrr.io/r/stats/offset.html)
  formula terms); fitted values and predictions include the offset, as
  in lmer.

- [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  accepts `subset`, `na.action` and `contrasts` with the same data
  preparation as
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md). Its
  `na.action` default is now `NULL` (refuse `NA`, like
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)); the
  old `na.omit` default was never applied.

- `na.action = na.exclude` pads
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) and in-sample
  [`predict()`](https://rdrr.io/r/stats/predict.html) to the original
  rows; the na.action record is kept on the fit (`fit$na.action`,
  model-frame attribute) and used by the emmeans bridge.

- Unused factor levels (fixed and grouping) are dropped after
  `subset`/NA removal, so
  [`ngrps()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  and coefficient sets match lme4.

- `predict(newdata)` keeps only the variables it needs (an unused Date
  column no longer fails), coerces integer/character grouping and
  character/logical fixed columns like the training frame, returns `NA`
  for incomplete rows (new `na.action` argument, default `na.pass` as in
  lme4), accepts partial `re.form` formulas such as `~ (1 | g)`, and
  evaluates formula offsets from `newdata` (an `offset =` argument
  offset must be resupplied via `predict(offset =)`).

- Interaction and nested grouping factors are labelled as in lme4
  (`a:b`, and `b:a` for the nested term of `a/b`, with levels `x:y`) in
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  and
  [`ngrps()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md);
  GLMM `predict(newdata)` works for them.

- Logical predictors are coded as factors with levels `FALSE`/`TRUE`, as
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) does, so
  no-intercept and margin-free interaction models match lme4 (`lFALSE`,
  `lTRUE`).

- Ordered grouping factors are grouped as unordered factors (no dense
  `contr.poly` basis is built for them).

- [`simulate()`](https://rdrr.io/r/stats/simulate.html) now works for
  GLMMs (binomial, Poisson, Gamma, negative binomial) and
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)
  re-fits a GLMM to a new response. Binomial
  [`cbind()`](https://rdrr.io/r/base/cbind.html) fits simulate
  two-column count matrices and proportion fits simulate proportions, as
  in lme4.

- Breaking: [`simulate()`](https://rdrr.io/r/stats/simulate.html) now
  uses lme4’s `re.form` meaning. The default `re.form = NA` (or `~0`, or
  `use.u = FALSE`) draws new random effects; `re.form = NULL` (or
  `use.u = TRUE`) conditions on the fitted random effects. Previously
  `NULL` drew new random effects and `NA` left them out entirely. Seeded
  draws now reproduce lme4’s when the estimates agree.

- New diagnostics:
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a fit
  draws Pearson residuals against fitted values (also Q-Q and
  scale-location via `which =`),
  [`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html) works on fits, and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`qqnorm()`](https://rdrr.io/r/stats/qqnorm.html)
  on
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  output draw caterpillar and Q-Q plots with conditional-variance
  intervals. With lattice loaded, `dotplot()` and `qqmath()` work on
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  output as in lme4.

- New [`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html),
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  and [`influence()`](https://rdrr.io/r/stats/lm.influence.html) (case
  or group deletion refits, with
  [`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html),
  [`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html) and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  methods for the result). LMM values match lme4. GLMM hat values use
  the working weights, as `glm` does, so they differ from lme4’s.

- New
  [`mm_r2()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md),
  [`mm_icc()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
  and
  [`mm_variance_components()`](https://bbuchsbaum.github.io/mixeff/reference/mm_variance_components.md)
  compute Nakagawa marginal/conditional R2 and adjusted/unadjusted ICC.
  With performance and insight installed,
  [`performance::r2()`](https://easystats.github.io/performance/reference/r2.html),
  [`performance::icc()`](https://easystats.github.io/performance/reference/icc.html)
  and
  [`insight::get_variance()`](https://easystats.github.io/insight/reference/get_variance.html)
  also work on mixeff fits and agree with their values on the matching
  lme4 fit.

- **Breaking:**
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  defaults to `method = "joint_laplace"` (glmer’s `nAGQ = 1` estimator,
  certified Wald inference). `method = "pirls_profiled"` remains
  available explicitly, and `nAGQ = 0` selects it (lme4’s fast
  estimate). Without an explicit `method`, negative-binomial families,
  `nAGQ > 1`, and `inference = "working_hessian"` use the profiled path
  with an `mm_estimator_notice`; explicit `method = "joint_laplace"`
  requests for them are refused.

- LMM bootstraps (`bootstrap_control(seed = NULL)`, the default) now
  draw their engine seed from R’s RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes them
  reproducible.

- [`simulate.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/simulate.mm_lmm.md)’s
  `"seed"` attribute follows
  [`stats::simulate`](https://rdrr.io/r/stats/simulate.html) (the prior
  `.Random.seed`, or the seed with its RNG `kind`).

- [`drop1()`](https://rdrr.io/r/stats/add1.html), random-term LRTs, and
  random-structure candidates keep a no-intercept model no-intercept
  instead of re-adding the intercept.

- [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) on
  a transformed response (`log(y) ~ ...`) refits to the new response
  (previously it returned the original fit); it keeps the fit’s
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).

- Internal refits no longer resolve `weights` through the data mask (a
  data column named `fit`, `full`, or `object` broke them), and
  REML-to-ML refits in
  [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)/[`anova()`](https://rdrr.io/r/stats/anova.html)
  keep the user’s
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md).

- [`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`ngrps()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
  and
  [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md)
  forward nlme (`lme`, `lmList`, `gls`) and lme4 objects to the owning
  package’s generic, so attaching mixeff no longer breaks nlme fits.

- [`summary()`](https://rdrr.io/r/base/summary.html)/`inference_table(method = )`
  compute all coefficient rows in one engine call instead of one refit
  per coefficient. Engine refits keep the fit’s full
  [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
  (a user `start` is no longer rounded to 4 digits on the wire).

- [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits
  every family/link pair the engine supports:
  [`Gamma()`](https://rdrr.io/r/stats/family.html) with its default
  inverse link,
  [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html) with
  `"inverse"`/`"log"`, and
  [`gaussian()`](https://rdrr.io/r/stats/family.html) with
  `"log"`/`"inverse"`/`"sqrt"` (glmer parity tests). Binomial
  cauchit/log/identity, Poisson identity, and inverse.gaussian’s default
  `1/mu^2` link are refused with a typed condition naming the supported
  set.

- Fits, refits, bootstraps and profiles can be interrupted (Ctrl-C/Esc):
  the engine checks for a pending interrupt between optimizer
  evaluations without longjmp-ing through Rust frames, and stops with a
  typed `mm_interrupted` error.

- New
  [`mm_bootmer()`](https://bbuchsbaum.github.io/mixeff/reference/mm_bootmer.md)
  ([`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html)
  counterpart; a `boot`-compatible result for
  [`boot::boot.ci()`](https://rdrr.io/pkg/boot/man/boot.ci.html)),
  [`rePCA()`](https://bbuchsbaum.github.io/mixeff/reference/rePCA.md),
  [`mm_lmlist()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmlist.md)
  (`lmList()`),
  [`mm_allfit()`](https://bbuchsbaum.github.io/mixeff/reference/mm_allfit.md)
  (`allFit()`), and `mm_control(optCtrl = )` (lme4-style optimizer
  controls, unknown names refused).

- `compare(small, big, method = "kenward_roger")` (and
  `"satterthwaite"`) gives pbkrtest’s `KRmodcomp()`/`SATmodcomp()` F
  test for nested fixed effects;
  [`step()`](https://bbuchsbaum.github.io/mixeff/reference/step.md)
  gives lmerTest-style backward elimination for `mm_lmm` fits (mixeff
  now exports a
  [`step()`](https://bbuchsbaum.github.io/mixeff/reference/step.md)
  generic whose default is
  [`stats::step()`](https://rdrr.io/r/stats/step.html)).

- Fast examples no longer use `\dontrun{}`; the fit result’s
  per-observation vectors cross the bridge as R doubles instead of JSON
  and are stored once in the fit object.

- Single-model [`anova()`](https://rdrr.io/r/stats/anova.html) tests a
  multi-column term such as `poly(x, 2)` as one multi-df row (joint F
  with the same df method), matching lmerTest, when no other term
  contains it.

- Engine error messages containing `%` (e.g. a formula with `%in%`, or a
  column name with `%s`) reach R verbatim: extendr passed the message to
  `Rf_error()` as a printf format string, which garbled them and could
  read out of bounds.

- Engine re-pinned to mixeff-rs 2873312. `||` with a factor now expands
  like lme4 (`(1 + x + f || g)` is
  `(1 | g) + (0 + x | g) + (0 + f | g)`), so parameter counts, `df`, AIC
  and estimates match lme4 (they used to differ when a factor sat inside
  `||`); random-effect factor bases use
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) coding.
  `diag(1 + f | g)` gives MixedModels.jl `zerocorr()` semantics.

- [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  `as.data.frame(VarCorr())`,
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`ngrps()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  and
  [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
  (`theta`, `Lambdat`, `Lind`, `Zt`, `cnms`, …) present random-effect
  terms in lme4’s `mkReTrms` order, with `||` pieces as separate entries
  `g`, `g.1`, …; `getME(fit, "theta")` is permuted to lme4’s order.
  [`getME()`](https://bbuchsbaum.github.io/mixeff/reference/getME.md)
  handles factor random-effect bases. `as.data.frame(VarCorr())` lists
  covariance rows in lme4’s order.

- `confint(fit, method = "profile")` /
  [`profile()`](https://rdrr.io/r/stats/profile.html) return lme4’s rows
  (`.sig01`, …, `.sigma`, fixed effects) on the SD/correlation scale and
  profile REML fits on the ML deviance, as lme4 does, so fixed-effect
  profile intervals are now available for REML fits (the
  `profile_beta_unavailable_under_reml` refusal is gone). Theta-scale
  rows stay in `attr(ci, "mm_profile")$table`. The LMM default remains
  Wald.

- New `threads` argument for bootstraps and profiles
  (`bootstrap_control(threads = )`,
  [`parametric_bootstrap()`](https://bbuchsbaum.github.io/mixeff/reference/parametric_bootstrap.md),
  `compare(method = "bootstrap")`, GLMM
  `confint(method = "bootstrap", threads = )`, LMM
  [`confint()`](https://rdrr.io/r/stats/confint.html)/[`profile()`](https://rdrr.io/r/stats/profile.html)):
  results are identical for every thread count; default 1 (CRAN’s
  two-thread policy); workers never call into R and Ctrl-C still
  interrupts.

- Rank-deficient fixed effects keep the earlier of two collinear columns
  (R’s rule).
  [`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  omits dropped coefficients (`add.dropped = TRUE` gives `NA`), and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) /
  [`summary()`](https://rdrr.io/r/base/summary.html) cover the estimable
  coefficients, matching lme4. Fixed: a dropped coefficient’s missing
  standard error shifted every later standard error by one position.

- Gamma and inverse-Gaussian GLMMs:
  [`logLik()`](https://rdrr.io/r/stats/logLik.html),
  [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html) and anova tables report
  glmer’s values (the engine’s density is 1 higher because glmer
  includes the family `aic()`’s `+2`); estimates,
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) and SEs match glmer.
  Inverse-Gaussian fits now have
  [`sigma()`](https://rdrr.io/r/stats/sigma.html),
  [`family()`](https://rdrr.io/r/stats/family.html) and a residual
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  row. Negative-binomial theta now maximizes the GLMM likelihood like
  `glmer.nb()`.

- [`test_random_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_random_effect.md)
  labels a correlation-only comparison (e.g. `||` versus `|`) as an
  ordinary chi-square test; the 50:50 boundary mixture is used only when
  exactly one variance is added.

- The aphantasia reproduction (`MIXEFF_RUN_APHANTASIA=true`) now reaches
  strict lme4 parity on every case with lme4’s `||` expansion: primary,
  intact and combined run the default joint-Laplace estimator, and the
  parity-ledger exemptions for the primary DiD estimate/SE, intact AIC
  and combined fixed effects are retired. The
  `inference = "working_hessian"` documentation no longer says its SEs
  run ~11% below glmer’s (that gap came from the old `||` family; they
  now agree within 1% on this dataset).

### Compatibility

- Fit-summary parsing accepts the additive `mixedmodels.fit_summary`
  1.1.0 schema emitted by mixeff-rs 1.0.0-rc.4, while retaining 1.0.0
  support. Covariance provenance is preserved in the stored payload.
  Unknown schema versions and malformed payloads remain errors.
- The bundled `mixeff-rs` engine moves from rc.1 (`1f3f689`) to rc.5
  plus pre-release fixes (`accd4b1`). This brings the upstream fixes for
  joint Laplace GLMM convergence without NLopt (and Ctrl-C during its
  inner PIRLS), GLMM parametric-bootstrap replicates that silently used
  the fast estimator, colliding `(1|a:b)` interaction keys, Wald
  p-values underflowing to 0, and negative-binomial deviance precision;
  it is also 2-7x faster on vector-valued and crossed models. New in
  this pin:
  - `y ~ 0 + f` codes the first factor with one column per level, as R’s
    [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) does;
    previously the reference level was silently dropped (its mean forced
    to zero).
  - `predict(newdata =)` uses the fit’s design coding, fixing wrong
    predictions for non-marginal formulas such as `y ~ f / h`.
  - LMM prior weights are validated (length, finite, positive) and
    simulation / parametric bootstrap draws residuals with sd
    `sigma / sqrt(w)`; binomial GLMM responses above 1 are refused
    instead of producing NaN estimates.
- `mixedmodels.fixed_effect_inference_table` 1.2.0 (per-row covariance
  provenance) is the current schema; stored 1.1.0 tables still parse.
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) for LMMs scales
  residual noise by prior weights, and keeps a zero-variance
  random-effect term at zero instead of borrowing another term’s
  variance.

### Breaking: API-shape stabilization

- [`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md) is
  now the single audit verb, dispatching on both compiled specs and
  fits;
  [`audit_design()`](https://bbuchsbaum.github.io/mixeff/reference/audit_design.md)
  forwards with a deprecation warning and will be removed later.
- [`df_for_contrast()`](https://bbuchsbaum.github.io/mixeff/reference/df_for_contrast.md)
  and
  [`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md)
  now return `mm_*` objects with `$table` (plus `$df` on the former,
  `$sections` for `reporting_table(section = "all")`), matching every
  sibling analysis verb, instead of a bare classed vector / data frame.
- The spec-accepting inspection verbs
  ([`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md),
  [`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md),
  [`reproducibility()`](https://bbuchsbaum.github.io/mixeff/reference/reproducibility.md),
  [`random_blocks()`](https://bbuchsbaum.github.io/mixeff/reference/random_blocks.md),
  [`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md),
  [`reporting_table()`](https://bbuchsbaum.github.io/mixeff/reference/model_report.md),
  [`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md))
  name their first argument `object` (previously `fit`, which was
  misleading for specs). Fit-only inference verbs keep `fit`.
  [`model.frame.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  keeps `formula` — that name is imposed by the
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  generic.
- [`confint()`](https://rdrr.io/r/stats/confint.html) presents
  `"asymptotic"` as the canonical method name (the package-wide term for
  the closed-form Wald interval); `"wald"` remains an accepted synonym.
  Computation is unchanged.
- [`drop1()`](https://rdrr.io/r/stats/add1.html) now matches
  [`stats::drop1()`](https://rdrr.io/r/stats/add1.html) marginality
  semantics: by default, main effects participating in an interaction
  are not offered for dropping. An explicit non-marginal `scope` is
  still honoured and fits normally. The result gains `status` and
  `reason` columns.

### Breaking: lme4-identical coefficient names

- Fixed-effect coefficient names now match
  `lme4`/[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  exactly — `"recipeB"`, `"temperature.L"`, `"recipeB:temperature.L"` —
  in [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  column order, on every programmatic surface:
  [`fixef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md),
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) tables,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) dimnames,
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  and
  [`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
  weight names, `tidy()`, `emmeans`, and
  [`predict()`](https://rdrr.io/r/stats/predict.html). Previously mixeff
  used its engine encoding (`"recipe: B"`) with a different interaction
  column order, so linear combinations and coefficient lookups
  copy-pasted from `lme4` code silently misaligned. Code written against
  the old names must switch to the lme4 forms. The engine encoding still
  appears inside engine-rendered
  `explain()`/[`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md)
  prose.
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  column names are stripped to the lme4 form (`"modalityAudio"`);
  [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  printing is unchanged. Fits saved with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) by older versions
  lack the stored name map and should be re-fit.

### Engine

- Engine pin bumped to `4a2abb3`: hardens convergence and runtime
  contracts. Non-marginal designs (`y ~ b + a:b`) now fit and match
  `lme4` exactly (previously refused). Convergence labelling is more
  conservative — a fit that reaches a flat/boundary region without a
  certified stationary point is now reported `not_optimized` (previously
  sometimes `converged_reduced_rank`) on small maximal random-slope
  models; genuine reduced-rank optima (e.g.
  [`lme4::Dyestuff2`](https://rdrr.io/pkg/lme4/man/Dyestuff.html)) still
  report `converged_reduced_rank`. The prior pin, `ee0c717`, carried the
  audit-render wording batch (policy recommendations phrased as options,
  boundary-sentence deduplication, humanized summary-view jargon).
- Earlier in this cycle the pin moved to `3b6ec69` (one commit past
  v1.0.0-rc.1), which fixes the native crossed-LMM trust-region start.
  Crossed-design fits that route through the trust-region optimizer may
  land on very slightly different (better-started) optima.
- The bundled `mixeff-rs` engine is now pinned to its first tagged
  release, v1.0.0-rc.1 (`3332f3e`). The two response-batch diagnostic
  reasons new in this release (`sink_stopped`, `adaptive_refinement`)
  are registered in the R-side reason registry, so the coverage contract
  (every engine reason has an R-side entry) holds.

### Extractors

- [`VarCorr()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  correlations are now stored at full precision as numeric columns of
  `$table` (`correlation`, plus `correlation2`, … for groups with three
  or more random terms; `NA` where no pair exists). Previously the
  column was a 2-decimal display string, forcing callers to parse text
  and costing ~2% precision on well-determined correlations. Printing is
  unchanged: rounding to 2 decimals now happens only at display time.

### Negative-binomial GLMMs

- [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  fits NB2 negative-binomial models (log link) on the default profiled
  path. `family = mm_negative_binomial()` estimates the size parameter
  theta alongside the model (the
  [`lme4::glmer.nb()`](https://rdrr.io/pkg/lme4/man/glmer.nb.html)
  route); `family = MASS::negative.binomial(theta)` or
  `mm_negative_binomial(theta)` fits conditional on a fixed theta. The
  fitted/fixed theta is recorded as `fit$family$nb_theta`.
  `method = "joint_laplace"` is not yet wired for this family at the
  pinned engine and is refused with a typed error.

### Clearer user-facing text (UX parity pass vs lme4)

A 13-scenario side-by-side battery against lme4 (graded independently)
drove a cleanup of every surface where engine internals leaked into
user-facing text:

- GLMM summaries with withheld inference now explain plainly that the
  fast default cannot certify SEs/z/p and to re-fit with
  `method = "joint_laplace"`; the engine’s covariance-geometry warrant
  moved behind `print(summary(fit), verbose = TRUE)`.
- Unsupported-family errors list the supported families and point to
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) for the
  rest.
- The new-grouping-level prediction error describes the R-level remedies
  (`re.form = NA`, `allow.new.levels = TRUE`) instead of Rust API names,
  and bridge errors no longer print a duplicated “Caused by” chain.
- [`anova()`](https://rdrr.io/r/stats/anova.html) prints a compact
  lme4-shaped table; single-df terms display as the equivalent F
  statistic (matching `lmerTest`), and provenance/list columns stay in
  `$table`.
- Aliased (rank-deficient) coefficients display as `NA` with an explicit
  note, instead of a misleading `0`.
- [`print()`](https://rdrr.io/r/base/print.html) no longer emits the
  artifact/crate provenance line (available on `fit$schema`);
  [`confint()`](https://rdrr.io/r/stats/confint.html)’s internal
  certification label is translated at display.

### Fixes and runtime notices

- [`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md)
  on a fitted GLMM reported the compile-time Lambda template (1s and 0s)
  as `theta_value` instead of the fitted theta; it now splices the
  fitted values in by index (LMM fits were unaffected; pre-fit specs
  keep the honest template).
- Logical random slopes now carry lme4-style
  [`ranef()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lmm-methods.md)
  column names (`"xTRUE"`), consistent with the fixed-effect naming and
  the conditional variance arrays.
- `glmm(method = "joint_laplace")` emits an up-front runtime notice: the
  joint route optimizes to an engine-chosen budget inside a single
  silent native call and can take minutes on large data (cap with
  `mm_control(max_feval = )`). Summary notes for completed joint fits no
  longer imply an unusable fit when the engine’s convergence label is
  `not_assessed`/`not_optimized` (label reliability is tracked
  upstream); they point to
  [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md).
- Bootstrap-based inference with 200+ replicates announces its scale
  before the single silent native call.
- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  now emit a rescaling advisory when a continuous predictor is on a
  scale far from 1 (matching `lme4`’s “predictors on very different
  scales” guidance): such fits can converge poorly, and
  [`scale()`](https://rdrr.io/r/base/scale.html) is the cheap fix. A
  notice, not a refusal; suppress with `mm_control(verbose = -1)`.
- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  document that optimization is silent and non-interruptible within one
  native call, with bounded budgets.

### Profiling and scope

- New
  [`profile.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/profile.mm_lmm.md)
  method: returns an `mm_profile` object over the engine’s certified
  profile-likelihood payload (`$table` with one row per parameter; REML
  fixed effects carry an explicit `profile_beta_unavailable_under_reml`
  reason instead of being dropped).
  [`confint()`](https://rdrr.io/r/stats/confint.html) on the profile
  reproduces `confint(fit, method = "profile")`.
- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) now
  refuses multivariate `cbind(y1, y2)` responses with a plain error (fit
  each outcome separately); shared-theta multivariate models are
  deferred post-release.
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  continues to accept `cbind(successes, failures)` for binomial
  responses.

### Prediction

- Population-level predictions (`re.form = NA` or `~0`) no longer
  require the random-effect grouping columns in `newdata`, matching
  `predict(lmer/glmer, re.form = NA)`. Only the fixed-part variables are
  needed; conditional predictions (`re.form = NULL`) still require the
  full formula’s variables.

### Contrasts

- Ordered factors are now coded with orthonormal polynomial contrasts
  (`contr.poly`) at fit time, matching R/lme4 defaults, instead of
  treatment coding. Fixed effects, random-slope (`Z`) coding,
  `logLik`/`AIC`/`BIC`, and predictions now reach parity with `lme4` on
  ordered-factor models (e.g.
  [`lme4::cake`](https://rdrr.io/pkg/lme4/man/cake.html)). Coefficient
  *names* still use mixeff’s engine encoding (`temperature: .L`) pending
  the lme4-identical renaming layer. If the global ordered-contrast
  option is not `contr.poly`, or an ordered column carries an explicit
  non-poly `contrasts` attribute,
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  refuse with a typed `mm_arg_error` rather than silently diverge from
  the requested coding.
- Behavior change: a *one-level* ordered factor now errors loudly
  (polynomial contrasts require at least two levels) instead of silently
  degenerating to an empty/near-empty design.

### Diagnostic clarity

- New verb
  [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md):
  re-runs a fitted LMM under the engine’s bounded verification workflow
  (restart from the optimum, jittered restarts, opt-in
  alternate-optimizer consensus) and reports the engine’s verdict with
  per-run objective/theta/beta deltas. This is the check the audit
  surface already pointed to for uncertain optima; the verdict and
  wording are engine-owned. `consensus` defaults to `FALSE` because this
  vendored build compiles without the optional `nlopt` backend, whose
  absence the consensus pass would otherwise report as a spurious
  `fragile`.
- [`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
  now prints one plain-language sentence per recorded change
  (e.g. `Fitted covariance for (1 | s): requested rank 1, fitted rank 0 [reduced_rank].`)
  instead of dumping the raw stage table. The certificate-time rank
  statement is treated as the canonical record of a boundary event, so
  its design/covariance restatements are not repeated (they remain in
  `$table`). A fit whose optimizer stopped early now says so explicitly
  (`none: no structural change was made; the optimizer stopped early (fit status \`not_optimized\`).`) instead of showing a misleading`unchanged
  / formula display\` row.
- The automatic pre-fit
  [`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
  block emitted by
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  now travels on the message stream (a typed `mm_explanation_notice`
  condition) instead of stdout, so
  [`suppressMessages()`](https://rdrr.io/r/base/message.html) and
  knitr’s `message = FALSE` can quiet it. It remains on by default;
  `mm_control(verbose = -1)` still suppresses it entirely, and an
  explicit `print(explain_model(spec))` still writes to stdout.
- [`summary()`](https://rdrr.io/r/base/summary.html) on a GLMM now
  defaults to `tests = "coefficients"` (matching
  [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html)). When the
  fit method cannot certify fixed-effect inference (the default
  `pirls_profiled` estimator), the SE/z/p columns are still withheld —
  but a `Notes:` line now states why, and that engine-certified Wald
  inference is available from a `method = "joint_laplace"` fit.
  Previously a default
  [`summary()`](https://rdrr.io/r/base/summary.html) printed `NA`
  columns with no explanation.
- [`summary()`](https://rdrr.io/r/base/summary.html) on a fit whose
  optimizer stopped without certifying an optimum (e.g. fit status
  `not_optimized`) now repeats that state as a plain-language `Notes:`
  line directly under the coefficient tests, instead of relying on the
  header status line alone.
- Displayed p-values now render through
  [`format.pval()`](https://rdrr.io/r/base/format.pval.html): an
  underflowed p-value prints as `< 1e-16` instead of `0.000000e+00`.
  Stored values are unchanged.
- The singular-fit [`print()`](https://rdrr.io/r/base/print.html) footer
  only advertises `random_options(spec, group = ...)` when that call can
  actually run for the fit (a slope candidate exists); previously the
  printed hint could error on the very fit that printed it.

### lme4 functional-equivalence layer

- Grouping-variable coercion:
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) /
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  coerce a non-categorical grouping variable (e.g. an integer subject id
  or numeric item code) to a factor for the random-effects structure,
  matching `lme4`/`nlme`/`glmmTMB`. Previously such a column was
  rejected by the native fit with “grouping factor not categorical”. The
  coercion is announced via a suppressible notice (class
  `mm_grouping_coercion_notice`; silence with
  `mm_control(verbose = -1)`), never silent. Surfaced by an in-the-wild
  OSF `glmer` reproduction with crossed `(1 | ID) + (1 | Title)`
  effects.
- Certified GLMM Wald inference: when fit with
  `method = "joint_laplace"`,
  [`summary()`](https://rdrr.io/r/base/summary.html),
  `confint(method = "wald")`,
  [`contrast()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md),
  and `tidy()` now report engine-certified fixed-effect standard errors,
  Wald *z* statistics, and *p*-values that match
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) within
  tolerance. The default `method = "pirls_profiled"` path is not
  certified for fixed-effect inference, so all four surfaces withhold
  SE/*z*/*p* (returning `NA` with a reason and a `vcov_status` of
  `"unsupported"`) rather than fabricate them from the uncertified
  working Hessian — consistent with the package’s “no fake certainty”
  contract.
- Conditional prediction standard errors and intervals:
  [`predict()`](https://rdrr.io/r/stats/predict.html) for `mm_lmm` /
  `mm_glmm` with `re.form = NULL` now routes `se.fit` and `interval`
  through the engine’s prediction-variance payload, which includes the
  random-effect (BLUP) variance and the fixed/random covariance — a
  surface
  [`lme4::predict.merMod`](https://rdrr.io/pkg/lme4/man/predict.merMod.html)
  does not offer at all. LMMs get conditional `se.fit` plus
  `"confidence"` and `"prediction"` intervals; GLMMs get conditional
  `se.fit` and `"confidence"` intervals on the link or response scale
  (variance propagated through the link by the engine). The engine
  certifies these rows for `method = "joint_laplace"` fits and — via a
  post-fit profiled-optimum certificate — for default `pirls_profiled`
  fits, so the default estimator now reports conditional SEs too. Rows
  the engine does not certify are withheld, not fabricated: uncertified
  fits (e.g. singular fits, whose certificate is never issued) and
  unseen grouping levels under `allow.new.levels = TRUE` return `NA`
  with the engine’s reason in the `mm_reason` attribute. Population
  (`re.form = NA`) SEs/intervals are unchanged.
- GLMM prediction (future-observation) intervals:
  `predict(interval = "prediction")` now works for conditional,
  response-scale GLMM predictions. Bounds are quantiles of the plug-in
  predictive distribution (the family conditional distribution mixed
  over link-scale fitted-mean uncertainty via Gauss–Hermite quadrature),
  so they are integers for count families and support points for
  Bernoulli; the interval is at least as wide as the corresponding
  confidence interval. Typed refusals remain for link-scale requests
  (future observations are response-scale objects), population-level
  requests, and grouped binomial fits (the future trial count is not
  representable in `newdata`).
- `||` factor-term semantics documented and contract-tested: in mixeff,
  zero-correlation syntax fully decorrelates the block — a factor’s
  treatment-coded level contrasts get independent variances with no
  within-factor covariances (the principled reading, shared by
  `afex::mixed(expand_re = TRUE)`, `glmmTMB::diag()`, and
  `MixedModels.jl zerocorr()`). `lme4`’s `||` instead leaves factor
  terms intact with a full within-factor covariance block, so the same
  formula fits a larger (and over-parameterized) model there. Fits
  announce the situation with an info diagnostic
  (`covariance_assumption`, reason `double_bar_factor_term`) naming the
  correlated-block rewrite (`(0 + f | g)`); the lme4-migration and
  formula vignettes carry the recipe.
- Binomial response coercion:
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) with
  `family = binomial()` now accepts a logical response (coerced 0/1
  silently) or a two-level factor response (coerced with the second
  level as success, announced via a suppressible `mm_factor_coercion`
  message), matching [`stats::glm()`](https://rdrr.io/r/stats/glm.html)
  / [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html). A factor
  with any other number of levels aborts with a typed `mm_data_error`.
  Previously these responses surfaced as an opaque engine error.
- [`update()`](https://rdrr.io/r/stats/update.html) for `mm_lmm` /
  `mm_glmm`: formula edits (`. ~ . - x`, preserving random-effect bars
  and `||`), `REML`/`weights`/`family`/ `offset`/`method`/`control`
  overrides, new `data`, and `evaluate = FALSE`.
- `broom` / `broom.mixed` support: `tidy()`, `glance()`, and `augment()`
  methods for `mm_lmm` / `mm_glmm` (registered with `generics`).
- [`confint.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/confint.mm_glmm.md):
  asymptotic Wald intervals for GLMM fixed effects (refuses
  profile/bootstrap with a typed reason).
- [`predict.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/predict.mm_glmm.md):
  `type = "link"`/`"response"`, population and conditional
  (`re.form = NULL`/`NA`) predictions with `allow.new.levels`, replacing
  the previous refusal. Validated against the engine’s
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).
- GLMM fixed-effect inference:
  [`contrast.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/contrast.md)
  (Wald),
  [`drop1.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/drop1.mm_glmm.md)
  (refit LRT), and
  [`anova.mm_glmm()`](https://bbuchsbaum.github.io/mixeff/reference/anova.mm_glmm.md)
  (sequential LRT for nested models).
- GLMM estimator transparency: the native `method = "joint_laplace"`
  path is certified against
  [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) within
  tolerance, and
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) now
  emits an informational notice (class `mm_estimator_notice`) when the
  default `pirls_profiled` estimator is used, since its coefficients are
  not glmer-equivalent (use `method = "joint_laplace"` for parity).

## mixeff 0.1.0

First public release of `mixeff`, an audit-first R wrapper around the
`mixedmodels` Rust crate. The package is distributed via R-Universe at
[bbuchsbaum.r-universe.dev](https://bbuchsbaum.r-universe.dev); the
upstream `nlopt` feature-gate PR that lands CRAN distribution is tracked
separately and ships as 0.2.0.

### Phase 0 — bridge and contract foundation

- `rextendr`/`extendr_api` bridge with vendored upstream `mixedmodels`
  crate; CRAN-compatible build with `cargo vendor` + `vendor.tar.xz`
  reconstitution at `R CMD INSTALL` time.
- [`mm_parse_formula()`](https://bbuchsbaum.github.io/mixeff/reference/mm_parse_formula.md)
  — R/Rust formula round-trip primitive.
- [`mm_formula_manifest()`](https://bbuchsbaum.github.io/mixeff/reference/mm_formula_manifest.md)
  — capability discovery for the bridge.
- [`mm_json_negotiate()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_negotiate.md)
  and
  [`mm_json_known_schemas()`](https://bbuchsbaum.github.io/mixeff/reference/mm_json_known_schemas.md)
  — schema versioning gate; mismatched artifacts raise `mm_schema_error`
  rather than silently misparse.
- Interrupt FFI: `Ctrl-C` during a long Rust fit cleanly returns to R.
  *\[Correction, 0.2.0: this applied only to the interrupt bridge demo
  (`mm_interrupt_demo`), never to real fits — ordinary
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  optimization runs in a single non-interruptible native call. See the
  0.2.0 documentation.\]*
- Typed condition catalog (`mm_condition` base class):
  `mm_formula_error`, `mm_data_error`, `mm_schema_error`,
  `mm_design_refusal`, `mm_inference_unavailable`, `mm_fit_error`,
  `mm_not_identifiable`, `mm_fit_not_optimized`.

### Phase 1 — audit-first construction surface

- [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  — formula + data → semantic IR + design audit, no fitting.
- [`audit_design()`](https://bbuchsbaum.github.io/mixeff/reference/audit_design.md)
  — structured design audit; raises `mm_design_refusal` for
  non-identifiable terms before any optimization runs.
- [`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
  — auto-printed once by
  [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) /
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  before the fit. Translates each random term into named-argument form,
  prints the per-block English gloss (authored in Rust), and emits the
  mandatory `No random slopes were added.` sentinel for intercept-only
  random terms.
- [`random_options()`](https://bbuchsbaum.github.io/mixeff/reference/random_options.md)
  — opt-in *map* of nearby random-effect spellings for a grouping factor
  (punt, slope-only, split-uncorrelated, double-bar, full). No
  “recommended” column; no preference ordering.
- [`compare_covariance()`](https://bbuchsbaum.github.io/mixeff/reference/compare_covariance.md)
  — full / diagonal / scalar comparison per random term.
- [`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md),
  [`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md),
  [`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md),
  [`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md),
  [`roles()`](https://bbuchsbaum.github.io/mixeff/reference/roles.md),
  [`as_json()`](https://bbuchsbaum.github.io/mixeff/reference/as_json.md),
  [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md).
- [`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) —
  REML/ML linear mixed-model fit via the upstream Rust
  `LinearMixedModel` engine. Returns a serializable `mm_lmm` carrying
  the JSON artifact, parsed state, beta/theta/sigma/logLik/deviance,
  fitted values, residuals, random effects, and varcorr.
- lme4-style extractor surface for `mm_lmm`: `fixef`, `ranef`, `coef`,
  `VarCorr`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `nobs`,
  `formula`, `model.frame`, `df.residual`, `fitted`,
  `residuals(type="response")`, basic
  [`predict()`](https://rdrr.io/r/stats/predict.html).
- Pedagogical `DiagnosticCode` variants surfaced from Rust: `ScopeNote`,
  `SupportNote`, `SyntaxExpansion`, `CovarianceAssumption`,
  `StructuralRefusal`.
- Random-term card schema (`RandomTermCard`) shipped per random term
  with `term_id`, `original_fragment`, `canonical_fragment`, group,
  blocks, `implied_constraints`, and `design_support`.
- R9 “no advice creep” contract enforced by `test-no-advice.R`: the
  strings `"suggested starting model"`, `"we recommend"`,
  `"you should"`, `"try ... instead"`, `"drop the random slope"`, and
  `"same model, different font"` cannot appear in package output.
- Vignettes: `intro.Rmd`, `lmm-basics.Rmd`, `demystifying-formulas.Rmd`.

### Phase 2 — saveRDS round-trip and lazy extractors

- `saveRDS` / `readRDS` survives without a live Rust handle — the
  artifact is the source of truth.
- [`revive()`](https://bbuchsbaum.github.io/mixeff/reference/revive.md)
  — rebuilds the Rust handle from the durable artifact when a live cache
  is needed.
- [`fit_handle_alive()`](https://bbuchsbaum.github.io/mixeff/reference/fit_handle_alive.md),
  `getME(fit, name)` for `X`, `Z`, `theta`, `Lambda`, `cnms`, `flist`,
  `Gp`, `lower`, `devcomp`, `optinfo`.
- `model.matrix(type=)`, `vcov(type="fixed")`.
- [`random_blocks()`](https://bbuchsbaum.github.io/mixeff/reference/random_blocks.md)
  — per-block decomposition of the random-effects matrix.
- [`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
  — convergence status, iterations, objective trace, verification trace.
- [`inference_table()`](https://bbuchsbaum.github.io/mixeff/reference/inference_table.md)
  — per-coefficient method/status/reliability rows read from the Rust
  inference contract.
- [`reproducibility()`](https://bbuchsbaum.github.io/mixeff/reference/reproducibility.md)
  — Rust-authored reproducibility envelope (engine version, schema
  version, seed, optimizer fingerprint).
- [`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
  — boolean predicate over the optimizer certificate.
- Vignette: `saving-and-reviving.Rmd`.

### Phase 3 — LMM inference

- `contrast(fit, L, rhs, method)` — fixed-effect contrast front door.
  Methods: `"auto"`, `"satterthwaite"`, `"kenward_roger"`,
  `"bootstrap"`, `"asymptotic"`, `"none"`. Returns `method` / `status` /
  `reliability` / `reason` columns; never fabricates p-values where the
  engine cannot certify a method.
- `test_effect(fit, term, method)` — term-level hypothesis tests.
  Bootstrap and bootstrap-LRT methods backed by the upstream Rust
  bootstrap entry points; cluster bootstrap is recognized but documented
  as estimator-distribution only (no certified p-value in schema 1.0.0).
- `inference_table(fit, method)` — multi-row inference table.
- [`df_for_contrast()`](https://bbuchsbaum.github.io/mixeff/reference/df_for_contrast.md),
  [`estimability()`](https://bbuchsbaum.github.io/mixeff/reference/estimability.md)
  — placeholders that return `NA` with a stable reason until 0.2.0 wires
  the Rust certificates end-to-end.
- [`anova()`](https://rdrr.io/r/stats/anova.html) — single and
  multi-model.
- [`drop1.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/drop1.mm_lmm.md).
- `confint(method = "wald")` — Wald asymptotic interval flagged with
  status `"not_certified_by_rust_inference_contract"`.
- `confint(method = "bootstrap")` — full-model bootstrap intervals with
  percentile / basic selection and bootstrap metadata.
- [`bootstrap_control()`](https://bbuchsbaum.github.io/mixeff/reference/bootstrap_control.md)
  — control object for bootstrap-backed methods (replicate count, seed,
  failed-refit policy).
- Vignettes: `inference.Rmd`, `inference-where-lme4-says-no.Rmd`.

### Phase 4 — GLMM boundary and LMM lifecycle

- [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) —
  Phase 4 boundary. The upstream Rust bridge does not yet expose a GLMM
  fit primitive, so
  [`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
  validates the family/link request, compiles the model spec, and raises
  a typed `mm_fit_error` with the expected `family` / `link` / `nAGQ`
  metadata until the bridge primitive lands. (Real GLMM fitting is
  queued for 0.2.0; upstream FFI is available.)
- [`simulate.mm_lmm()`](https://bbuchsbaum.github.io/mixeff/reference/simulate.mm_lmm.md)
  — simulate from a fitted LMM using the durable artifact state.
- [`refit()`](https://bbuchsbaum.github.io/mixeff/reference/refit.md) —
  refit with a new response.
- [`compare()`](https://bbuchsbaum.github.io/mixeff/reference/compare.md)
  — model comparison with auditable validity status.
- Multi-model [`anova()`](https://rdrr.io/r/stats/anova.html) and
  [`drop1()`](https://rdrr.io/r/stats/add1.html) over `mm_lmm` objects.
- [`parametric_bootstrap()`](https://bbuchsbaum.github.io/mixeff/reference/parametric_bootstrap.md)
  — parametric bootstrap distribution for fixed-effect tests.
- Manifest capabilities for `simulate`/`inference` exposed via the
  bridge contract.
- Vignettes: `glmm.Rmd` (boundary walkthrough), `benchmarking.Rmd`,
  `reporting-lmms.Rmd`.

### Cross-cutting infrastructure

- [`mm_control()`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
  — flat named list mirroring `lmerControl`. Honored fields include
  `optimizer`, `optimizer_max_iter`, `optimizer_xtol_abs`,
  `optimizer_ftol_abs`, `reml`, `nAGQ`, `verify_convergence`,
  `parallel_threads`, `seed`, `verbose`, `thresholds`, `schema_version`,
  `bridge_timeout_s`.
- `mm_thresholds()` — design/identifiability thresholds (byte-equivalent
  to the upstream `compiler_contract_v0_prd.md` §8).
- Parity ledger (`inst/extdata/expected-mismatches.json`) — every
  divergence from `lme4` is classified (`expected_mismatch` /
  `upstream_bug` / `unsupported`) with bounds enforced by
  `tests/testthat/helper-parity-scoreboard.R`.
- Parity scoreboard (`test-parity-scoreboard.R`) — emits a structured
  artifact recording observed differences against tolerances on the
  classic `lme4` parity baseline.
- Speedup vs `lme4` on the included scaling benchmark
  (`benchmarks/lme4-scaling/`) ranges from ~2× (small balanced LMMs) to
  ~5× (correlated random slopes on ≥30 grouping levels).

### Non-goals (preserved)

- `mixeff` is not a drop-in `lme4` replacement.
- No bit-exact numerical reproduction of `lme4`.
- No model-selection or random-effects recommendation engine (no
  `recommend_model()`, `auto_random_effects()`, `fix_singularity()`,
  `make_it_converge()`).
- [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) /
  [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) are not
  masked on attach.
