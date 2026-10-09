# mixeff 0.2.0

## Correctness fixes (pre-release audit)

* Character fixed-effect predictors (e.g. columns from `read.csv()`) are now
  converted to factors with sorted levels, as lme4 does. Previously the fit
  aborted with an `mm_schema_error` while building coefficient names.
* A numeric predictor whose name starts with a factor's name (`group` and
  `group_size`) no longer aborts the fit.
* Contrast codings mixeff cannot honour are refused with an `mm_arg_error`
  instead of being silently replaced by treatment coding: a contrast attached
  to an unordered factor, a non-`contr.treatment` unordered entry in
  `options(contrasts =)`, and `contrasts = list(f = "contr.SAS")` (which uses
  the last level as reference; the engine always uses the first).
* `update()` re-evaluates the original `data` when the updated formula needs
  columns the stored model frame lacks (as lme4 does), and refuses arguments
  it cannot carry over (e.g. `subset`, `na.action`) instead of dropping them.
* The declared Rust toolchain minimum is now 1.85, the bundled engine's
  actual requirement (some vendored crates use edition 2024).
* The lme4 migration vignette no longer describes `coef()`, `VarCorr()`,
  `residuals()`, `anova()`, `getME()`, `isSingular()` and `emmeans` as
  identical to lme4; it lists the actual differences.

## lme4 parity and checklist completion

* Shapes follow lme4 >= 2.0 and insight >= 1.5:
  * `confint(method = "profile")` numbers `.sigNN` term by term, each
    term's standard deviations before its correlations. lme4 1.1 numbered
    them in lower-triangle order.
  * `anova(m1, m2)` names its `-2 * logLik` column `-2*log(L)`, where it
    was `deviance`.
  * `getME(m, "devcomp")$dims` gains `npar`.
  * `print(VarCorr(m))` uses `reformulas::formatVC()` when reformulas is
    installed.
  * For GLMMs, `mm_r2()` / `mm_icc()` divide the binomial
    distribution-specific variance of a `cbind()` response by the mean
    number of trials. For count models, the null-model mean is
    `exp(b0 + v0/2)`, where `b0` is the null model's intercept and `v0` its
    random-effect variance.

* `residuals()` follows lme4: LMM `type = "pearson"`/`"deviance"` are
  `sqrt(w) * (y - mu)` (no longer divided by sigma), and `scaled = TRUE`
  divides by `sigma()` exactly once (it divided twice for Pearson). GLMM
  residuals default to `"deviance"` and support `"pearson"`, `"working"` and
  `"response"`, computed from the family, prior weights and binomial trials.
* `sigma()` of a negative-binomial GLMM is 1, as in lme4; theta is
  `getME(fit, "glmer.nb.theta")`. GLMM `deviance()` is the sum of squared
  deviance residuals (as lme4), no longer `-2 * logLik()`.
* `coef()` returns every fixed-effect column per group in lme4's order.
* `VarCorr()` returns lme4's structure (named list of covariance matrices with
  `stddev`/`correlation` attributes, `sc`/`useSc`), so `VarCorr(m)$Subject`
  is the matrix; it prints like lme4. `VarCorr(m)$table` and `$residual_sd`
  still work.
* `AIC(m1, m2)` / `BIC(m1, m2)` return stats' `df`/`AIC` data frame instead of
  refusing.
* New `REMLcrit()`, `family()` methods, `weights(type = "working")`;
  `weights()` returns ones for unweighted fits; `model.frame()` has lme4's
  transformed columns and `terms` attribute (raw variables stay in
  `fit$model_frame`).
* `getME()` covers lme4's components (`u`, `b`, `L`, `RX`, `RZX`, `devcomp`,
  `lower`, `Gp`, `Tp`, `Lind`, `Ztlist`, `Tlist`, `ST`, `offset`, `weights`,
  `glmer.nb.theta`, `"ALL"`, ...) for LMMs and GLMMs; `devfun` is refused with
  a typed error. `Lambda`/`Lambdat` (and `model.matrix(fit, type = "random")`)
  now follow the engine's (lme4's) term order; previously theta was assigned
  in formula order, so crossed designs such as Penicillin's
  `(1|sample) + (1|plate)` gave a term another term's theta.
* `is_singular()` uses lme4's rule (a zero-bounded theta below `tol`), honours
  `tol`, and works for GLMMs.
* `ranef(condVar = TRUE)` returns Laplace conditional variances for GLMMs
  instead of `NA`.
* `summary()` and `anova()` accept lmerTest's `ddf = "Satterthwaite"`,
  `"Kenward-Roger"` or `"lme4"`. `summary()$coefficients` is lmerTest's
  numeric matrix (per-row method labels moved to `summary()$coef_table`).
  `anova(m)` is lmerTest's data frame (`Sum Sq Mean Sq NumDF DenDF F value
  Pr(>F)`) and `anova(m1, m2)` lme4's (`npar AIC BIC logLik deviance Chisq Df
  Pr(>Chisq)`); mixeff's provenance table stays in `$table`.
* `confint()` accepts lme4's `method = "Wald"` and `"boot"` spellings.
* Stateful and expanded fixed-effect formula terms now work in `lmm()` and
  `glmm()`: `factor()`, `relevel()`, `cut()`, `scale()`, `poly()`,
  `splines::ns()`/`bs()`, `as.numeric()`, `I(x > 0)`, `offset()`, `^`,
  `%in%`, and parenthesised groups such as `(a + b)^2` and `a * (b + c)`.
  Coefficient names, `terms()`, `model.matrix()`, `drop1()` and `emmeans`
  grids follow lme4, and `predict(newdata)` reuses the training basis.
* `lmm()` gains `offset =` (and `offset()` formula terms); fitted values and
  predictions include the offset, as in lmer.
* `glmm()` now accepts `subset`, `na.action` and `contrasts` with the same
  data preparation as `lmm()`. Its `na.action` default is now `NULL`
  (refuse `NA`, like `lmm()`); the old `na.omit` default was never applied.
* `na.action = na.exclude` pads `fitted()`, `residuals()` and in-sample
  `predict()` to the original rows; the na.action record is kept on the fit
  (`fit$na.action`, model-frame attribute) and used by the emmeans bridge.
* Unused factor levels (fixed and grouping) are dropped after `subset`/NA
  removal, so `ngrps()` and coefficient sets match lme4.
* `predict(newdata)` keeps only the variables it needs (an unused Date
  column no longer fails), coerces integer/character grouping and
  character/logical fixed columns like the training frame, returns `NA` for
  incomplete rows (new `na.action` argument, default `na.pass` as in lme4),
  accepts partial `re.form` formulas such as `~ (1 | g)`, and evaluates
  formula offsets from `newdata` (an `offset =` argument offset must be
  resupplied via `predict(offset =)`).
* Interaction and nested grouping factors are labelled as in lme4 (`a:b`,
  and `b:a` for the nested term of `a/b`, with levels `x:y`) in `ranef()`,
  `VarCorr()` and `ngrps()`; GLMM `predict(newdata)` works for them.
* Logical predictors are coded as factors with levels `FALSE`/`TRUE`, as
  `model.matrix()` does, so no-intercept and margin-free interaction models
  match lme4 (`lFALSE`, `lTRUE`).
* Ordered grouping factors are grouped as unordered factors (no dense
  `contr.poly` basis is built for them).
* `simulate()` now works for GLMMs (binomial, Poisson, Gamma, negative
  binomial) and `refit()` re-fits a GLMM to a new response. Binomial
  `cbind()` fits simulate two-column count matrices and proportion fits
  simulate proportions, as in lme4.
* Breaking: `simulate()` now uses lme4's `re.form` meaning. The default
  `re.form = NA` (or `~0`, or `use.u = FALSE`) draws new random effects;
  `re.form = NULL` (or `use.u = TRUE`) conditions on the fitted random
  effects. Previously `NULL` drew new random effects and `NA` left them out
  entirely. Seeded draws now reproduce lme4's when the estimates agree.
* New diagnostics: `plot()` on a fit draws Pearson residuals against fitted
  values (also Q-Q and scale-location via `which =`), `qqnorm()` works on
  fits, and `plot()`/`qqnorm()` on `ranef()` output draw caterpillar and Q-Q
  plots with conditional-variance intervals. With lattice loaded,
  `dotplot()` and `qqmath()` work on `ranef()` output as in lme4.
* New `hatvalues()`, `cooks.distance()` and `influence()` (case or group
  deletion refits, with `dfbeta()`, `dfbetas()` and `cooks.distance()`
  methods for the result). LMM values match lme4. GLMM hat values use the
  working weights, as `glm` does, so they differ from lme4's.
* New `mm_r2()`, `mm_icc()` and `mm_variance_components()` compute
  Nakagawa marginal/conditional R2 and adjusted/unadjusted ICC. With
  performance and insight installed, `performance::r2()`,
  `performance::icc()` and `insight::get_variance()` also work on mixeff
  fits and agree with their values on the matching lme4 fit.
* **Breaking:** `glmm()` now defaults to `method = "joint_laplace"` (glmer's
  `nAGQ = 1` estimator, certified Wald inference). `method = "pirls_profiled"`
  remains available explicitly, and `nAGQ = 0` selects it (lme4's fast
  estimate). Without an explicit `method`, negative-binomial families,
  `nAGQ > 1`, and `inference = "working_hessian"` use the profiled path with
  an `mm_estimator_notice`; explicit `method = "joint_laplace"` requests for
  them are refused.
* LMM bootstraps (`bootstrap_control(seed = NULL)`, the default) now draw
  their engine seed from R's RNG, so `set.seed()` makes them reproducible.
* `simulate.mm_lmm()`'s `"seed"` attribute follows `stats::simulate` (the
  prior `.Random.seed`, or the seed with its RNG `kind`).
* `drop1()`, random-term LRTs, and random-structure candidates keep a
  no-intercept model no-intercept instead of re-adding the intercept.
* `refit()` on a transformed response (`log(y) ~ ...`) refits to the new
  response (previously it returned the original fit); it keeps the fit's
  `mm_control()`.
* Internal refits no longer resolve `weights` through the data mask (a data
  column named `fit`, `full`, or `object` broke them), and REML-to-ML refits
  in `compare()`/`anova()` keep the user's `mm_control()`.
* `fixef()`, `ranef()`, `VarCorr()`, `ngrps()`, `getME()` and `refit()`
  forward nlme (`lme`, `lmList`, `gls`) and lme4 objects to the owning
  package's generic, so attaching mixeff no longer breaks nlme fits.
* `summary()`/`inference_table(method = )` compute all coefficient rows in
  one engine call instead of one refit per coefficient. Engine refits keep
  the fit's full `mm_control()` (a user `start` is no longer rounded to 4
  digits on the wire).
* `glmm()` fits every family/link pair the engine supports: `Gamma()` with
  its default inverse link, `inverse.gaussian()` with `"inverse"`/`"log"`,
  and `gaussian()` with `"log"`/`"inverse"`/`"sqrt"` (glmer parity tests).
  Binomial cauchit/log/identity, Poisson identity, and inverse.gaussian's
  default `1/mu^2` link are refused with a typed condition naming the
  supported set.
* Fits, refits, bootstraps and profiles can be interrupted (Ctrl-C/Esc):
  the engine checks for a pending interrupt between optimizer evaluations
  without longjmp-ing through Rust frames, and stops with a typed
  `mm_interrupted` error.
* New `mm_bootmer()` (`lme4::bootMer()` counterpart; a `boot`-compatible
  result for `boot::boot.ci()`), `rePCA()`, `mm_lmlist()` (`lmList()`),
  `mm_allfit()` (`allFit()`), and `mm_control(optCtrl = )` (lme4-style
  optimizer controls, unknown names refused).
* `compare(small, big, method = "kenward_roger")` (and `"satterthwaite"`)
  gives pbkrtest's `KRmodcomp()`/`SATmodcomp()` F test for nested fixed
  effects; `step()` gives lmerTest-style backward elimination for `mm_lmm`
  fits (mixeff now exports a `step()` generic whose default is
  `stats::step()`).
* Fast examples no longer use `\dontrun{}`; the fit result's
  per-observation vectors cross the bridge as R doubles instead of JSON and
  are stored once in the fit object.

* Single-model `anova()` tests a multi-column term such as `poly(x, 2)` as
  one multi-df row (joint F with the same df method), matching lmerTest,
  when no other term contains it.
* Engine error messages containing `%` (e.g. a formula with `%in%`, or a
  column name with `%s`) reach R verbatim: extendr passed the message to
  `Rf_error()` as a printf format string, which garbled them and could read
  out of bounds.
* Engine re-pinned to mixeff-rs 2873312. `||` with a factor now expands
  like lme4 (`(1 + x + f || g)` is `(1 | g) + (0 + x | g) + (0 + f | g)`),
  so parameter counts, `df`, AIC and estimates match lme4 (they used to
  differ when a factor sat inside `||`); random-effect factor bases use
  `model.matrix()` coding. `diag(1 + f | g)` gives MixedModels.jl
  `zerocorr()` semantics.
* `VarCorr()`, `as.data.frame(VarCorr())`, `ranef()`, `ngrps()` and
  `getME()` (`theta`, `Lambdat`, `Lind`, `Zt`, `cnms`, ...) present
  random-effect terms in lme4's `mkReTrms` order, with `||` pieces as
  separate entries `g`, `g.1`, ...; `getME(fit, "theta")` is permuted to
  lme4's order. `getME()` handles factor random-effect bases.
  `as.data.frame(VarCorr())` lists covariance rows in lme4's order.
* `confint(fit, method = "profile")` / `profile()` return lme4's rows
  (`.sig01`, ..., `.sigma`, fixed effects) on the SD/correlation scale and
  profile REML fits on the ML deviance, as lme4 does, so fixed-effect
  profile intervals are now available for REML fits (the
  `profile_beta_unavailable_under_reml` refusal is gone). Theta-scale rows
  stay in `attr(ci, "mm_profile")$table`. The LMM default remains Wald.
* New `threads` argument for bootstraps and profiles
  (`bootstrap_control(threads = )`, `parametric_bootstrap()`,
  `compare(method = "bootstrap")`, GLMM `confint(method = "bootstrap",
  threads = )`, LMM `confint()`/`profile()`): results are identical for
  every thread count; default 1 (CRAN's two-thread policy); workers never
  call into R and Ctrl-C still interrupts.
* Rank-deficient fixed effects keep the earlier of two collinear columns
  (R's rule). `fixef()` omits dropped coefficients (`add.dropped = TRUE`
  gives `NA`), and `vcov()` / `summary()` cover the estimable coefficients,
  matching lme4. Fixed: a dropped coefficient's missing standard error
  shifted every later standard error by one position.
* Gamma and inverse-Gaussian GLMMs: `logLik()`, `AIC()`, `BIC()` and
  anova tables report glmer's values (the engine's density is 1 higher
  because glmer includes the family `aic()`'s `+2`); estimates, `sigma()`
  and SEs match glmer. Inverse-Gaussian fits now have `sigma()`,
  `family()` and a residual `VarCorr()` row. Negative-binomial theta now
  maximizes the GLMM likelihood like `glmer.nb()`.
* `test_random_effect()` labels a correlation-only comparison (e.g. `||`
  versus `|`) as an ordinary chi-square test; the 50:50 boundary mixture
  is used only when exactly one variance is added.
* The aphantasia reproduction (`MIXEFF_RUN_APHANTASIA=true`) now reaches
  strict lme4 parity on every case with lme4's `||` expansion: primary,
  intact and combined run the default joint-Laplace estimator, and the
  parity-ledger exemptions for the primary DiD estimate/SE, intact AIC and
  combined fixed effects are retired. The `inference = "working_hessian"`
  documentation no longer says its SEs run ~11% below glmer's (that gap came
  from the old `||` family; they now agree within 1% on this dataset).

## Compatibility

* Fit-summary parsing accepts the additive `mixedmodels.fit_summary` 1.1.0
  schema emitted by mixeff-rs 1.0.0-rc.4, while retaining 1.0.0 support.
  Covariance provenance is preserved in the stored payload. Unknown schema
  versions and malformed payloads remain errors.
* The bundled `mixeff-rs` engine moves from rc.1 (`1f3f689`) to rc.5 plus
  pre-release fixes (`accd4b1`). This brings the upstream fixes for joint
  Laplace GLMM convergence without NLopt (and Ctrl-C during its inner PIRLS),
  GLMM parametric-bootstrap replicates that silently used the fast
  estimator, colliding `(1|a:b)` interaction keys, Wald p-values underflowing
  to 0, and negative-binomial deviance precision; it is also 2-7x faster on
  vector-valued and crossed models. New in this pin:
  * `y ~ 0 + f` codes the first factor with one column per level, as R's
    `model.matrix()` does; previously the reference level was silently
    dropped (its mean forced to zero).
  * `predict(newdata =)` uses the fit's design coding, fixing wrong
    predictions for non-marginal formulas such as `y ~ f / h`.
  * LMM prior weights are validated (length, finite, positive) and
    simulation / parametric bootstrap draws residuals with sd
    `sigma / sqrt(w)`; binomial GLMM responses above 1 are refused instead of
    producing NaN estimates.
* `mixedmodels.fixed_effect_inference_table` 1.2.0 (per-row covariance
  provenance) is the current schema; stored 1.1.0 tables still parse.
* `simulate()` for LMMs scales residual noise by prior weights, and keeps a
  zero-variance random-effect term at zero instead of borrowing another
  term's variance.

## Breaking: API-shape stabilization

* `audit()` is now the single audit verb, dispatching on both compiled specs
  and fits; `audit_design()` forwards with a deprecation warning and will be
  removed later.
* `df_for_contrast()` and `reporting_table()` now return `mm_*` objects with
  `$table` (plus `$df` on the former, `$sections` for
  `reporting_table(section = "all")`), matching every sibling analysis verb,
  instead of a bare classed vector / data frame.
* The spec-accepting inspection verbs (`changes()`, `parameterization()`,
  `reproducibility()`, `random_blocks()`, `optimizer_certificate()`,
  `reporting_table()`, `audit()`) name their first argument `object`
  (previously `fit`, which was misleading for specs). Fit-only inference verbs
  keep `fit`. `model.frame.mm_lmm()` keeps `formula` — that name is imposed by
  the `stats::model.frame()` generic.
* `confint()` presents `"asymptotic"` as the canonical method name (the
  package-wide term for the closed-form Wald interval); `"wald"` remains an
  accepted synonym. Computation is unchanged.
* `drop1()` now matches `stats::drop1()` marginality semantics: by default,
  main effects participating in an interaction are not offered for dropping.
  An explicit non-marginal `scope` is still honoured and fits normally. The
  result gains `status` and `reason` columns.

## Breaking: lme4-identical coefficient names

* Fixed-effect coefficient names now match `lme4`/`model.matrix()` exactly —
  `"recipeB"`, `"temperature.L"`, `"recipeB:temperature.L"` — in
  `model.matrix()` column order, on every programmatic surface: `fixef()`,
  `coef()`, `summary()` tables, `vcov()` dimnames, `confint()`, `contrast()`
  and `mm_lincomb()` weight names, `tidy()`, `emmeans`, and `predict()`.
  Previously mixeff used its engine encoding (`"recipe: B"`) with a different
  interaction column order, so linear combinations and coefficient lookups
  copy-pasted from `lme4` code silently misaligned. Code written against the
  old names must switch to the lme4 forms. The engine encoding still appears
  inside engine-rendered `explain()`/`audit()` prose. `ranef()` column names
  are stripped to the lme4 form (`"modalityAudio"`); `VarCorr()` printing is
  unchanged. Fits saved with `saveRDS()` by older versions lack the stored
  name map and should be re-fit.

## Engine

* Engine pin bumped to `4a2abb3`: hardens convergence and runtime contracts.
  Non-marginal designs (`y ~ b + a:b`) now fit and match `lme4` exactly
  (previously refused). Convergence labelling is more conservative — a fit
  that reaches a flat/boundary region without a certified stationary point is
  now reported `not_optimized` (previously sometimes `converged_reduced_rank`)
  on small maximal random-slope models; genuine reduced-rank optima (e.g.
  `lme4::Dyestuff2`) still report `converged_reduced_rank`. The prior pin,
  `ee0c717`, carried the audit-render wording batch (policy recommendations
  phrased as options, boundary-sentence deduplication, humanized summary-view
  jargon).
* Earlier in this cycle the pin moved to `3b6ec69` (one commit past
  v1.0.0-rc.1), which fixes the native crossed-LMM trust-region start.
  Crossed-design fits that route through the trust-region optimizer may land
  on very slightly different (better-started) optima.
* The bundled `mixeff-rs` engine is now pinned to its first tagged release,
  v1.0.0-rc.1 (`3332f3e`). The two response-batch diagnostic reasons new in
  this release (`sink_stopped`, `adaptive_refinement`) are registered in the
  R-side reason registry, so the coverage contract (every engine reason has
  an R-side entry) holds.

## Extractors

* `VarCorr()` correlations are now stored at full precision as numeric
  columns of `$table` (`correlation`, plus `correlation2`, ... for groups
  with three or more random terms; `NA` where no pair exists). Previously
  the column was a 2-decimal display string, forcing callers to parse text
  and costing ~2% precision on well-determined correlations. Printing is
  unchanged: rounding to 2 decimals now happens only at display time.

## Negative-binomial GLMMs

* `glmm()` now fits NB2 negative-binomial models (log link) on the default
  profiled path. `family = mm_negative_binomial()` estimates the size
  parameter theta alongside the model (the `lme4::glmer.nb()` route);
  `family = MASS::negative.binomial(theta)` or `mm_negative_binomial(theta)`
  fits conditional on a fixed theta. The fitted/fixed theta is recorded as
  `fit$family$nb_theta`. `method = "joint_laplace"` is not yet wired for this
  family at the pinned engine and is refused with a typed error.

## Clearer user-facing text (UX parity pass vs lme4)

A 13-scenario side-by-side battery against lme4 (graded independently) drove
a cleanup of every surface where engine internals leaked into user-facing
text:

* GLMM summaries with withheld inference now explain plainly that the fast
  default cannot certify SEs/z/p and to re-fit with
  `method = "joint_laplace"`; the engine's covariance-geometry warrant moved
  behind `print(summary(fit), verbose = TRUE)`.
* Unsupported-family errors list the supported families and point to
  `lme4::glmer()` for the rest.
* The new-grouping-level prediction error describes the R-level remedies
  (`re.form = NA`, `allow.new.levels = TRUE`) instead of Rust API names, and
  bridge errors no longer print a duplicated "Caused by" chain.
* `anova()` prints a compact lme4-shaped table; single-df terms display as
  the equivalent F statistic (matching `lmerTest`), and provenance/list
  columns stay in `$table`.
* Aliased (rank-deficient) coefficients display as `NA` with an explicit
  note, instead of a misleading `0`.
* `print()` no longer emits the artifact/crate provenance line (available on
  `fit$schema`); `confint()`'s internal certification label is translated at
  display.

## Fixes and runtime notices

* `parameterization()` on a fitted GLMM reported the compile-time Lambda
  template (1s and 0s) as `theta_value` instead of the fitted theta; it now
  splices the fitted values in by index (LMM fits were unaffected; pre-fit
  specs keep the honest template).
* Logical random slopes now carry lme4-style `ranef()` column names
  (`"xTRUE"`), consistent with the fixed-effect naming and the conditional
  variance arrays.
* `glmm(method = "joint_laplace")` emits an up-front runtime notice: the
  joint route optimizes to an engine-chosen budget inside a single silent
  native call and can take minutes on large data (cap with
  `mm_control(max_feval = )`). Summary notes for completed joint fits no
  longer imply an unusable fit when the engine's convergence label is
  `not_assessed`/`not_optimized` (label reliability is tracked upstream);
  they point to `verify_convergence()`.
* Bootstrap-based inference with 200+ replicates announces its scale before
  the single silent native call.
* `lmm()`/`glmm()` now emit a rescaling advisory when a continuous predictor
  is on a scale far from 1 (matching `lme4`'s "predictors on very different
  scales" guidance): such fits can converge poorly, and `scale()` is the
  cheap fix. A notice, not a refusal; suppress with
  `mm_control(verbose = -1)`.
* `lmm()`/`glmm()` document that optimization is silent and non-interruptible
  within one native call, with bounded budgets.

## Profiling and scope

* New `profile.mm_lmm()` method: returns an `mm_profile` object over the
  engine's certified profile-likelihood payload (`$table` with one row per
  parameter; REML fixed effects carry an explicit
  `profile_beta_unavailable_under_reml` reason instead of being dropped).
  `confint()` on the profile reproduces `confint(fit, method = "profile")`.
* `lmm()` now refuses multivariate `cbind(y1, y2)` responses with a plain
  error (fit each outcome separately); shared-theta multivariate models are
  deferred post-release. `glmm()` continues to accept
  `cbind(successes, failures)` for binomial responses.

## Prediction

* Population-level predictions (`re.form = NA` or `~0`) no longer require the
  random-effect grouping columns in `newdata`, matching
  `predict(lmer/glmer, re.form = NA)`. Only the fixed-part variables are
  needed; conditional predictions (`re.form = NULL`) still require the full
  formula's variables.

## Contrasts

* Ordered factors are now coded with orthonormal polynomial contrasts
  (`contr.poly`) at fit time, matching R/lme4 defaults, instead of treatment
  coding. Fixed effects, random-slope (`Z`) coding, `logLik`/`AIC`/`BIC`, and
  predictions now reach parity with `lme4` on ordered-factor models (e.g.
  `lme4::cake`). Coefficient *names* still use mixeff's engine encoding
  (`temperature: .L`) pending the lme4-identical renaming layer. If the global
  ordered-contrast option is not `contr.poly`, or an ordered column carries an
  explicit non-poly `contrasts` attribute, `lmm()`/`glmm()` refuse with a typed
  `mm_arg_error` rather than silently diverge from the requested coding.
* Behavior change: a *one-level* ordered factor now errors loudly (polynomial
  contrasts require at least two levels) instead of silently degenerating to an
  empty/near-empty design.

## Diagnostic clarity

* New verb `verify_convergence()`: re-runs a fitted LMM under the engine's
  bounded verification workflow (restart from the optimum, jittered restarts,
  opt-in alternate-optimizer consensus) and reports the engine's verdict with
  per-run objective/theta/beta deltas. This is the check the audit surface
  already pointed to for uncertain optima; the verdict and wording are
  engine-owned. `consensus` defaults to `FALSE` because this vendored build
  compiles without the optional `nlopt` backend, whose absence the consensus
  pass would otherwise report as a spurious `fragile`.
* `changes()` now prints one plain-language sentence per recorded change
  (e.g. `Fitted covariance for (1 | s): requested rank 1, fitted rank 0
  [reduced_rank].`) instead of dumping the raw stage table. The
  certificate-time rank statement is treated as the canonical record of a
  boundary event, so its design/covariance restatements are not repeated
  (they remain in `$table`). A fit whose optimizer stopped early now says so
  explicitly (`none: no structural change was made; the optimizer stopped
  early (fit status \`not_optimized\`).`) instead of showing a misleading
  `unchanged / formula display` row.
* The automatic pre-fit `explain_model()` block emitted by `lmm()`/`glmm()`
  now travels on the message stream (a typed `mm_explanation_notice`
  condition) instead of stdout, so `suppressMessages()` and knitr's
  `message = FALSE` can quiet it. It remains on by default;
  `mm_control(verbose = -1)` still suppresses it entirely, and an explicit
  `print(explain_model(spec))` still writes to stdout.
* `summary()` on a GLMM now defaults to `tests = "coefficients"` (matching
  `lme4::glmer`). When the fit method cannot certify fixed-effect inference
  (the default `pirls_profiled` estimator), the SE/z/p columns are still
  withheld — but a `Notes:` line now states why, and that engine-certified
  Wald inference is available from a `method = "joint_laplace"` fit.
  Previously a default `summary()` printed `NA` columns with no explanation.
* `summary()` on a fit whose optimizer stopped without certifying an optimum
  (e.g. fit status `not_optimized`) now repeats that state as a plain-language
  `Notes:` line directly under the coefficient tests, instead of relying on
  the header status line alone.
* Displayed p-values now render through `format.pval()`: an underflowed
  p-value prints as `< 1e-16` instead of `0.000000e+00`. Stored values are
  unchanged.
* The singular-fit `print()` footer only advertises
  `random_options(spec, group = ...)` when that call can actually run for the
  fit (a slope candidate exists); previously the printed hint could error on
  the very fit that printed it.

## lme4 functional-equivalence layer

* Grouping-variable coercion: `lmm()` / `glmm()` now coerce a non-categorical
  grouping variable (e.g. an integer subject id or numeric item code) to a
  factor for the random-effects structure, matching `lme4`/`nlme`/`glmmTMB`.
  Previously such a column was rejected by the native fit with
  "grouping factor not categorical". The coercion is announced via a
  suppressible notice (class `mm_grouping_coercion_notice`; silence with
  `mm_control(verbose = -1)`), never silent. Surfaced by an in-the-wild OSF
  `glmer` reproduction with crossed `(1 | ID) + (1 | Title)` effects.
* Certified GLMM Wald inference: when fit with `method = "joint_laplace"`,
  `summary()`, `confint(method = "wald")`, `contrast()`, and `tidy()` now report
  engine-certified fixed-effect standard errors, Wald *z* statistics, and
  *p*-values that match `lme4::glmer()` within tolerance. The default
  `method = "pirls_profiled"` path is not certified for fixed-effect inference,
  so all four surfaces withhold SE/*z*/*p* (returning `NA` with a reason and a
  `vcov_status` of `"unsupported"`) rather than fabricate them from the
  uncertified working Hessian — consistent with the package's "no fake
  certainty" contract.
* Conditional prediction standard errors and intervals: `predict()` for
  `mm_lmm` / `mm_glmm` with `re.form = NULL` now routes `se.fit` and
  `interval` through the engine's prediction-variance payload, which includes
  the random-effect (BLUP) variance and the fixed/random covariance — a
  surface `lme4::predict.merMod` does not offer at all. LMMs get conditional
  `se.fit` plus `"confidence"` and `"prediction"` intervals; GLMMs get
  conditional `se.fit` and `"confidence"` intervals on the link or response
  scale (variance propagated through the link by the engine). The engine
  certifies these rows for `method = "joint_laplace"` fits and — via a
  post-fit profiled-optimum certificate — for default `pirls_profiled` fits,
  so the default estimator now reports conditional SEs too. Rows the engine
  does not certify are withheld, not fabricated: uncertified fits (e.g.
  singular fits, whose certificate is never issued) and unseen grouping
  levels under `allow.new.levels = TRUE` return `NA` with the engine's reason
  in the `mm_reason` attribute. Population (`re.form = NA`) SEs/intervals are
  unchanged.
* GLMM prediction (future-observation) intervals:
  `predict(interval = "prediction")` now works for conditional,
  response-scale GLMM predictions. Bounds are quantiles of the plug-in
  predictive distribution (the family conditional distribution mixed over
  link-scale fitted-mean uncertainty via Gauss–Hermite quadrature), so they
  are integers for count families and support points for Bernoulli; the
  interval is at least as wide as the corresponding confidence interval.
  Typed refusals remain for link-scale requests (future observations are
  response-scale objects), population-level requests, and grouped binomial
  fits (the future trial count is not representable in `newdata`).
* `||` factor-term semantics documented and contract-tested: in mixeff,
  zero-correlation syntax fully decorrelates the block — a factor's
  treatment-coded level contrasts get independent variances with no
  within-factor covariances (the principled reading, shared by
  `afex::mixed(expand_re = TRUE)`, `glmmTMB::diag()`, and
  `MixedModels.jl zerocorr()`). `lme4`'s `||` instead leaves factor terms
  intact with a full within-factor covariance block, so the same formula
  fits a larger (and over-parameterized) model there. Fits announce the
  situation with an info diagnostic (`covariance_assumption`, reason
  `double_bar_factor_term`) naming the correlated-block rewrite
  (`(0 + f | g)`); the lme4-migration and formula vignettes carry the
  recipe.
* Binomial response coercion: `glmm()` with `family = binomial()` now
  accepts a logical response (coerced 0/1 silently) or a two-level factor
  response (coerced with the second level as success, announced via a
  suppressible `mm_factor_coercion` message), matching `stats::glm()` /
  `lme4::glmer()`. A factor with any other number of levels aborts with a
  typed `mm_data_error`. Previously these responses surfaced as an opaque
  engine error.
* `update()` for `mm_lmm` / `mm_glmm`: formula edits (`. ~ . - x`,
  preserving random-effect bars and `||`), `REML`/`weights`/`family`/
  `offset`/`method`/`control` overrides, new `data`, and `evaluate = FALSE`.
* `broom` / `broom.mixed` support: `tidy()`, `glance()`, and `augment()`
  methods for `mm_lmm` / `mm_glmm` (registered with `generics`).
* `confint.mm_glmm()`: asymptotic Wald intervals for GLMM fixed effects
  (refuses profile/bootstrap with a typed reason).
* `predict.mm_glmm()`: `type = "link"`/`"response"`, population and
  conditional (`re.form = NULL`/`NA`) predictions with `allow.new.levels`,
  replacing the previous refusal. Validated against the engine's `fitted()`.
* GLMM fixed-effect inference: `contrast.mm_glmm()` (Wald), `drop1.mm_glmm()`
  (refit LRT), and `anova.mm_glmm()` (sequential LRT for nested models).
* GLMM estimator transparency: the native `method = "joint_laplace"` path is
  certified against `lme4::glmer` within tolerance, and `glmm()` now emits an
  informational notice (class `mm_estimator_notice`) when the default
  `pirls_profiled` estimator is used, since its coefficients are not
  glmer-equivalent (use `method = "joint_laplace"` for parity).

# mixeff 0.1.0

First public release of `mixeff`, an audit-first R wrapper around the
`mixedmodels` Rust crate. The package is distributed via R-Universe at
[bbuchsbaum.r-universe.dev](https://bbuchsbaum.r-universe.dev); the
upstream `nlopt` feature-gate PR that lands CRAN distribution is
tracked separately and ships as 0.2.0.

## Phase 0 — bridge and contract foundation

* `rextendr`/`extendr_api` bridge with vendored upstream `mixedmodels`
  crate; CRAN-compatible build with `cargo vendor` + `vendor.tar.xz`
  reconstitution at `R CMD INSTALL` time.
* `mm_parse_formula()` — R/Rust formula round-trip primitive.
* `mm_formula_manifest()` — capability discovery for the bridge.
* `mm_json_negotiate()` and `mm_json_known_schemas()` — schema
  versioning gate; mismatched artifacts raise `mm_schema_error` rather
  than silently misparse.
* Interrupt FFI: `Ctrl-C` during a long Rust fit cleanly returns to R.
  *[Correction, 0.2.0: this applied only to the interrupt bridge demo
  (`mm_interrupt_demo`), never to real fits — ordinary `lmm()`/`glmm()`
  optimization runs in a single non-interruptible native call. See the
  0.2.0 documentation.]*
* Typed condition catalog (`mm_condition` base class):
  `mm_formula_error`, `mm_data_error`, `mm_schema_error`,
  `mm_design_refusal`, `mm_inference_unavailable`, `mm_fit_error`,
  `mm_not_identifiable`, `mm_fit_not_optimized`.

## Phase 1 — audit-first construction surface

* `compile_model()` — formula + data → semantic IR + design audit, no
  fitting.
* `audit_design()` — structured design audit; raises
  `mm_design_refusal` for non-identifiable terms before any
  optimization runs.
* `explain_model()` — auto-printed once by `lmm()` / `glmm()` before
  the fit. Translates each random term into named-argument form,
  prints the per-block English gloss (authored in Rust), and emits the
  mandatory `No random slopes were added.` sentinel for
  intercept-only random terms.
* `random_options()` — opt-in *map* of nearby random-effect spellings
  for a grouping factor (punt, slope-only, split-uncorrelated,
  double-bar, full). No "recommended" column; no preference ordering.
* `compare_covariance()` — full / diagonal / scalar comparison per
  random term.
* `changes()`, `diagnostics()`, `fit_status()`,
  `parameterization()`, `roles()`, `as_json()`, `is_singular()`.
* `lmm()` — REML/ML linear mixed-model fit via the upstream Rust
  `LinearMixedModel` engine. Returns a serializable `mm_lmm` carrying
  the JSON artifact, parsed state, beta/theta/sigma/logLik/deviance,
  fitted values, residuals, random effects, and varcorr.
* lme4-style extractor surface for `mm_lmm`: `fixef`, `ranef`, `coef`,
  `VarCorr`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `nobs`,
  `formula`, `model.frame`, `df.residual`, `fitted`,
  `residuals(type="response")`, basic `predict()`.
* Pedagogical `DiagnosticCode` variants surfaced from Rust:
  `ScopeNote`, `SupportNote`, `SyntaxExpansion`,
  `CovarianceAssumption`, `StructuralRefusal`.
* Random-term card schema (`RandomTermCard`) shipped per random term
  with `term_id`, `original_fragment`, `canonical_fragment`, group,
  blocks, `implied_constraints`, and `design_support`.
* R9 "no advice creep" contract enforced by `test-no-advice.R`:
  the strings `"suggested starting model"`, `"we recommend"`,
  `"you should"`, `"try ... instead"`, `"drop the random slope"`,
  and `"same model, different font"` cannot appear in package output.
* Vignettes: `intro.Rmd`, `lmm-basics.Rmd`, `demystifying-formulas.Rmd`.

## Phase 2 — saveRDS round-trip and lazy extractors

* `saveRDS` / `readRDS` survives without a live Rust handle — the
  artifact is the source of truth.
* `revive()` — rebuilds the Rust handle from the durable artifact when
  a live cache is needed.
* `fit_handle_alive()`, `getME(fit, name)` for `X`, `Z`, `theta`,
  `Lambda`, `cnms`, `flist`, `Gp`, `lower`, `devcomp`, `optinfo`.
* `model.matrix(type=)`, `vcov(type="fixed")`.
* `random_blocks()` — per-block decomposition of the random-effects
  matrix.
* `optimizer_certificate()` — convergence status, iterations,
  objective trace, verification trace.
* `inference_table()` — per-coefficient method/status/reliability rows
  read from the Rust inference contract.
* `reproducibility()` — Rust-authored reproducibility envelope (engine
  version, schema version, seed, optimizer fingerprint).
* `is_singular()` — boolean predicate over the optimizer certificate.
* Vignette: `saving-and-reviving.Rmd`.

## Phase 3 — LMM inference

* `contrast(fit, L, rhs, method)` — fixed-effect contrast front door.
  Methods: `"auto"`, `"satterthwaite"`, `"kenward_roger"`,
  `"bootstrap"`, `"asymptotic"`, `"none"`. Returns
  `method` / `status` / `reliability` / `reason` columns; never
  fabricates p-values where the engine cannot certify a method.
* `test_effect(fit, term, method)` — term-level hypothesis tests.
  Bootstrap and bootstrap-LRT methods backed by the upstream Rust
  bootstrap entry points; cluster bootstrap is recognized but
  documented as estimator-distribution only (no certified p-value in
  schema 1.0.0).
* `inference_table(fit, method)` — multi-row inference table.
* `df_for_contrast()`, `estimability()` — placeholders that return
  `NA` with a stable reason until 0.2.0 wires the Rust certificates
  end-to-end.
* `anova()` — single and multi-model.
* `drop1.mm_lmm()`.
* `confint(method = "wald")` — Wald asymptotic interval flagged with
  status `"not_certified_by_rust_inference_contract"`.
* `confint(method = "bootstrap")` — full-model bootstrap intervals
  with percentile / basic selection and bootstrap metadata.
* `bootstrap_control()` — control object for bootstrap-backed methods
  (replicate count, seed, failed-refit policy).
* Vignettes: `inference.Rmd`, `inference-where-lme4-says-no.Rmd`.

## Phase 4 — GLMM boundary and LMM lifecycle

* `glmm()` — Phase 4 boundary. The upstream Rust bridge does not yet
  expose a GLMM fit primitive, so `glmm()` validates the family/link
  request, compiles the model spec, and raises a typed `mm_fit_error`
  with the expected `family` / `link` / `nAGQ` metadata until the
  bridge primitive lands. (Real GLMM fitting is queued for 0.2.0;
  upstream FFI is available.)
* `simulate.mm_lmm()` — simulate from a fitted LMM using the durable
  artifact state.
* `refit()` — refit with a new response.
* `compare()` — model comparison with auditable validity status.
* Multi-model `anova()` and `drop1()` over `mm_lmm` objects.
* `parametric_bootstrap()` — parametric bootstrap distribution for
  fixed-effect tests.
* Manifest capabilities for `simulate`/`inference` exposed via the
  bridge contract.
* Vignettes: `glmm.Rmd` (boundary walkthrough),
  `benchmarking.Rmd`, `reporting-lmms.Rmd`.

## Cross-cutting infrastructure

* `mm_control()` — flat named list mirroring `lmerControl`. Honored
  fields include `optimizer`, `optimizer_max_iter`, `optimizer_xtol_abs`,
  `optimizer_ftol_abs`, `reml`, `nAGQ`, `verify_convergence`,
  `parallel_threads`, `seed`, `verbose`, `thresholds`,
  `schema_version`, `bridge_timeout_s`.
* `mm_thresholds()` — design/identifiability thresholds (byte-equivalent
  to the upstream `compiler_contract_v0_prd.md` §8).
* Parity ledger (`inst/extdata/expected-mismatches.json`) — every
  divergence from `lme4` is classified
  (`expected_mismatch` / `upstream_bug` / `unsupported`) with bounds
  enforced by `tests/testthat/helper-parity-scoreboard.R`.
* Parity scoreboard (`test-parity-scoreboard.R`) — emits a structured
  artifact recording observed differences against tolerances on the
  classic `lme4` parity baseline.
* Speedup vs `lme4` on the included scaling benchmark
  (`benchmarks/lme4-scaling/`) ranges from
  ~2× (small balanced LMMs) to ~5× (correlated random slopes on
  ≥30 grouping levels).

## Non-goals (preserved)

* `mixeff` is not a drop-in `lme4` replacement.
* No bit-exact numerical reproduction of `lme4`.
* No model-selection or random-effects recommendation engine (no
  `recommend_model()`, `auto_random_effects()`, `fix_singularity()`,
  `make_it_converge()`).
* `lme4::lmer` / `lme4::glmer` are not masked on attach.
