# mixeff (development version)

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
