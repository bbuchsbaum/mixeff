# Changelog

## mixeff (development version)

### Engine

- Engine pin bumped to `3b6ec69` (one commit past v1.0.0-rc.1): fixes
  the native crossed-LMM trust-region start. Crossed-design fits that
  route through the trust-region optimizer may land on very slightly
  different (better-started) optima.
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
  ([`benchmarks/lme4-scaling/`](https://bbuchsbaum.github.io/mixeff/news/benchmarks/lme4-scaling/))
  ranges from ~2× (small balanced LMMs) to ~5× (correlated random slopes
  on ≥30 grouping levels).

### Non-goals (preserved)

- `mixeff` is not a drop-in `lme4` replacement.
- No bit-exact numerical reproduction of `lme4`.
- No model-selection or random-effects recommendation engine (no
  `recommend_model()`, `auto_random_effects()`, `fix_singularity()`,
  `make_it_converge()`).
- [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) /
  [`lme4::glmer`](https://rdrr.io/pkg/lme4/man/glmer.html) are not
  masked on attach.
