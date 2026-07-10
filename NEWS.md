# mixeff (development version)

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

* Engine pin bumped to `3b6ec69` (one commit past v1.0.0-rc.1): fixes the
  native crossed-LMM trust-region start. Crossed-design fits that route
  through the trust-region optimizer may land on very slightly different
  (better-started) optima.
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
  ([`benchmarks/lme4-scaling/`](benchmarks/lme4-scaling/)) ranges from
  ~2× (small balanced LMMs) to ~5× (correlated random slopes on
  ≥30 grouping levels).

## Non-goals (preserved)

* `mixeff` is not a drop-in `lme4` replacement.
* No bit-exact numerical reproduction of `lme4`.
* No model-selection or random-effects recommendation engine (no
  `recommend_model()`, `auto_random_effects()`, `fix_singularity()`,
  `make_it_converge()`).
* `lme4::lmer` / `lme4::glmer` are not masked on attach.
