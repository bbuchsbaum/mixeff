# PRD — `mixeff`: an R wrapper for the `mixedmodels` Rust crate

## 1. Context

The `mixedmodels` Rust crate at `/Users/bbuchsbaum/code/rust/mixedmodels` is a port of Julia's
`MixedModels.jl` providing LMM and GLMM fitting, formula parsing, design auditing, a
ThetaMap parameterization, and a serializable `CompiledModelArtifact`. It currently has no
language bindings.

The R proposal (`/Users/bbuchsbaum/code/rust/mixedmodels/docs/r_layer_proposal.md`)
specifies that R should be a **client** to the Rust compiler: R captures user intent
(formula, mode, inference target), Rust owns model semantics, convergence, covariance
reductions, and inference availability. Every printed claim in R must trace back to a
versioned Rust JSON artifact (proposal §"Source of truth"). The two contracts that govern
the wire surface are `compiler_contract_v0_prd.md` and
`mixed_model_compiler_inference_contract.md`.

This plan establishes the R package — provisional name **`mixeff`**, matching the empty
target directory `/Users/bbuchsbaum/code/mixeff` — its build system, API surface, internal
JSON contract plumbing, phased roadmap, and CRAN-readiness posture.

## 2. Goals

- An audit-first R package that exposes `mixedmodels` to R users via lme4-style formulas
  while *not* pretending to be lme4.
- The R object is a serializable record of Rust artifacts; the external pointer is a cache,
  not the source of truth (proposal §"Persistence guarantees", lines 300–308).
- Honest inference: where Rust says a method is unavailable, R prints `NA` with the Rust
  reason; we never fabricate p-values (contract §"No fake certainty").
- A path to CRAN, with R-Universe as the development distribution channel.

## 3. Explicit non-goals

- Not a drop-in `lme4` replacement (proposal lines 22–24).
- No silent model surgery; reductions and refusals always cross the boundary (proposal
  lines 61–62).
- **No bit-exact lme4 numerical reproduction** — statistical equivalence within documented
  tolerances on parity datasets only. This is a new non-goal not in the PRD that the plan
  recommends adding.
- v0 will *not* ship: profile-likelihood CIs, GAM smooths, residual covariance structures
  (AR(1), spatial), Kenward-Roger, regularized model search, multivariate cross-outcome
  residual covariance, `I()`/`poly()`/spline helpers in formulas, in-fit offsets.
- We do *not* mask `lme4::lmer`/`glmer` on attach.
- **No model selection or random-effects recommendation engine.** No
  `recommend_model()`, `auto_random_effects()`, `fix_singularity()`, or
  `make_it_converge()` verbs. The demystification helpers describe consequences
  and design support (parameter counts, what each option assumes is zero, what
  the data can support); they never rank, select, or substitute models. The
  package adds a *random-effects guidance layer*, not a *model-selection layer*
  (see §9.5).
- v1 ships passive explanation only. An interactive `re_builder()` is deferred
  to v2.

## 4. Architecture overview

```
+-------------------------------+        JSON (versioned)         +---------------------------+
| R user                        |  <----------------------------> | mixedmodels (Rust crate)  |
|  lmm(y ~ x + (x|g), data)     |  model_spec  +  data_payload    |  formula -> SemanticModel |
|        |                      |  -------------------------->    |  audit_design()           |
|        v                      |                                 |  fit() (PLS / PIRLS)      |
|  mm_fit (S3)                  |  artifact + audit + diagnostics |  CompiledModelArtifact    |
|    $artifact_json (truth)     |  + theta_map + cert + inference |  ThetaMap (sum type)      |
|    $rust_handle (cache)       |  + reproducibility + state      |  OptimizerCertificate     |
|    $lazy_cache (X, Z, Lambda) |  <-------------------------     |                           |
|        |                      |                                 |                           |
|        v                      |  prediction / simulation /      |                           |
|  S3 methods, audit(), revive()|  refit / contrast (lazy)        |                           |
+-------------------------------+                                 +---------------------------+
                                          extendr_api bridge
```

R holds the durable, parsed JSON state. The Rust handle is opaque and may die at any time
(GC, fork, restart); `revive(fit)` rebuilds it from `artifact_json`.

## 5. The Rust side — what exists, what we need

### 5.1 Available today (verified by exploration)

- Single library crate, edition 2021. No FFI yet.
- Formula parsing: `parse_formula(&str)` → `Formula` (lme4-compatible: `*`, `:`, `/`, `&`,
  `(re|g)`, `(re||g)`, intercept handling).
- `LinearMixedModel::new(formula, data, weights)` then `.fit(reml)`,
  `.predict()`/`.predict_new()`, `.simulate()`, `.refit()`, `.theta()`/`.set_theta()`,
  `.objective_at()`, `.varcorr()`, `.audit_report()`, `set_compiler_policy()`.
- `GeneralizedLinearMixedModel::pirls()`, `.deviance(n_agq)`.
- `CompiledModelArtifact`, `ThetaMap` (Scalar | Diagonal | FullCholesky | Structured |
  ReducedRank), all `serde`-serializable.
- Convergence statuses: `ConvergedInterior` | `ConvergedBoundary` | `ConvergedReducedRank`
  | `ConvergedPenalised` | `NotIdentifiable` | `NotOptimized`.
- Custom `DataFrame { columns: IndexMap<String, Column>, n_rows }` input. R must
  materialize and copy.
- No `rust-version` declared.

### 5.2 Upstream changes the wrapper requires (preconditions)

These should be filed as upstream issues in `mixedmodels` *before* Phase 1 begins:

1. **Gate `nlopt` behind a Cargo feature.** `nlopt` 0.8.1 wraps a C library
   that requires CMake at build time; this is CRAN-incompatible. The change
   needed is narrow because the optimizer abstraction already exists upstream
   — `src/types/opt_summary.rs` defines `enum Optimizer { Cobyla,
   PatternSearch, NloptNewuoa, NloptBobyqa }`, and `cobyla = "1"` is already
   a non-optional dep. The required upstream patch:
   `[features] nlopt = ["dep:nlopt"]; default = []`, mark `nlopt` as
   `optional = true`, and put the existing `nlopt` call sites (~10
   functions, concentrated in `src/model/linear.rs` and
   `src/types/opt_summary.rs`) behind `#[cfg(feature = "nlopt")]`. The two
   selectors (`use_nlopt_bobyqa_small_theta_optimizer`,
   `use_large_theta_nlopt_optimizer`) return `false` when the feature is
   off; the existing dispatch in `fit()` then routes to `cobyla` /
   `pattern_search`. No new optimizer trait, no new `argmin` dep. `nlopt`
   stays available behind `--features nlopt` for R-Universe / dev builds
   where BOBYQA parity with lme4 is the goal. Full draft of the upstream
   issue (with line numbers, acceptance criteria, CI matrix) lives at
   `planning/upstream-nlopt-issue.md`.
2. **Pin `rust-version = "1.78"`** in `Cargo.toml` to match nalgebra 0.33 +
   extendr_api, with a small buffer.
3. **Public schema-version constant** — a `pub const SCHEMA_VERSION: &str` on each artifact
   type so `mm_json_negotiate()` can fail fast on mismatch.
4. **Stable diagnostic codes** — confirm the contract's stable codes
   (`FixedRandomRedundant`, `RandomSlopeUnsupported`, `CovarianceTooRich`,
   `NotIdentifiable`, etc.) are exposed as a closed `enum DiagnosticCode` with
   `Display`/`serde` impls, not free strings. Extend the enum with the
   pedagogical (non-error) variants documented in §9.7: `ScopeNote`,
   `SupportNote`, `SyntaxExpansion`, `CovarianceAssumption`,
   `StructuralRefusal`.
5. **Random term card per random-effect term** — `audit_design()` and the
   `CompiledModelArtifact` must return, for each random term, a structured card
   with `term_id`, `original_fragment`, `canonical_fragment`, `group`, per-block
   `{basis, intercept, slopes, covariance, theta_parameters, english}`,
   `implied_constraints` (e.g., zero-covariance between blocks, with `reason`),
   and `design_support` (`group_levels`, `min_rows_per_group`,
   `median_rows_per_group`, `within_group_variation` per slope candidate, plus
   a `status` value drawn from the existing closed set). R formats these into
   `explain_model()`, `parameterization()`, `random_options()`,
   `compare_covariance()`, and `changes()` without reconstructing semantics.
   Plain-English block phrasing and constraint reasons are authored in Rust
   (single source of truth for wording across language bindings). See §9.6 for
   schema.

If item 1 is blocked, **fallback**: ship v0.1 to R-Universe only; defer CRAN until the
optimizer abstraction lands. Surfaced as Decision A below.

## 6. Package skeleton (day-1 file tree)

```
mixeff/
  DESCRIPTION                NAMESPACE        LICENSE        LICENSE.note
  NEWS.md  README.Rmd  README.md  cran-comments.md
  .Rbuildignore  .gitignore
  R/
    mixeff-package.R   zzz.R   extendr-wrappers.R (generated)
    formula.R   data-translate.R   compile.R
    fit-lmm.R   fit-glmm.R   mm-control.R
    audit.R   demystify.R   methods-print.R   methods-summary.R   methods-extract.R
    predict.R   simulate.R   compare.R
    revive.R   json.R   conditions.R
  src/
    Makevars   Makevars.win   Makevars.ucrt
    entrypoint.c   mixeff-win.def
    rust/
      Cargo.toml   Cargo.lock   .cargo/config.toml
      src/lib.rs       src/handle.rs       src/json_bridge.rs
      src/data.rs      src/interrupt.rs
      vendor/                       # cargo vendor output, checked in
  inst/
    schemas/   fixtures/   extdata/   AUTHORS   CITATION
  tests/testthat/        # see §10
  tools/
    vendor-rust.R   msrv-check.R   schema-snapshot.R   update-rust-crate.R
  vignettes/
    intro.Rmd   formula-dsl.Rmd   demystifying-formulas.Rmd
    audit-first-workflow.Rmd   lme4-migration.Rmd
  data/   data-raw/   pkgdown/_pkgdown.yml
  .github/workflows/
    R-CMD-check.yaml   pkgdown.yaml   schema-snapshot.yaml
    vendor-sync.yaml   coverage.yaml   lintr.yaml
```

## 7. Build & packaging

- **Bridge: `rextendr` + `extendr_api`.** Standard, expected by reviewers, generates
  `Makevars{,.win,.ucrt}`, `entrypoint.c`, and the Windows `.def` correctly.
- **MSRV: Rust 1.78** (`SystemRequirements: Cargo (>= 1.78.0), GNU make` in DESCRIPTION;
  `rust-version = "1.78"` in `src/rust/Cargo.toml`).
- **Vendoring:** `tools/vendor-rust.R` runs `cargo vendor`, sets
  `[source.crates-io] replace-with = "vendored-sources"`, regenerates `LICENSE.note` from
  every transitive crate's license, and stores the SHA256 of `Cargo.lock` in
  `tools/cargo-lock.sha256`. CI guards drift.
- **Tarball slimming:** `.Rbuildignore` excludes `src/rust/target`, `src/rust/vendor/*/{tests,examples,benches,*.png,*.md}`
  (preserving LICENSE/COPYING/NOTICE). `[profile.release] strip = "symbols", lto = "thin"`.
  Target tarball ≤ 10 MB; if it overshoots, justify in `cran-comments.md` with the
  precedent of `prqlr`, `gifski`, `salso`.
- **Cross-platform:** macOS arm64 + x86_64, Ubuntu LTS, Windows UCRT (Rtools43+ MinGW —
  `x86_64-pc-windows-gnu`, **not** MSVC). Drop Windows i386.
- **CRAN posture:** `lme4` in `Suggests` (used only in vignettes and parity tests under
  `skip_if_not_installed`). Document `R_RUST_TOOLCHAIN` precedent (`prqlr`, `string`,
  `polars`).

## 8. R-facing API

### 8.1 S3 class hierarchy (S3 only — push back on S4)

| Class | Inherits | Created by |
| --- | --- | --- |
| `mm_fit` | — | internal root |
| `mm_lmm` | `mm_fit` | `lmm()` |
| `mm_glmm` | `mm_fit` | `glmm()` |
| `mm_spec` | — | `compile_model()` |
| `mm_audit` / `mm_diagnostics` / `mm_theta_map` / `mm_inference_table` / `mm_change_log` / `mm_optimizer_certificate` / `mm_reproducibility` / `mm_estimability` / `mm_random_blocks` | `mm_report` (where applicable) | the matching verb |

Typed conditions inherit `mm_condition`: `mm_not_identifiable`, `mm_design_refusal`,
`mm_fit_not_optimized`, `mm_inference_unavailable`, `mm_schema_error`, `mm_formula_error`.

**Naming-collision call:** namespace verbs and rely on `::` for `compare()`, `contrast()`,
`audit()`, `changes()`. The package help for `mixeff::contrast` opens with the lead
"Note: this is not `emmeans::contrast`." A single opt-in helper
`mixeff::use_mixeff_verbs()` exposes unprefixed names — explicit and scriptable, never
default.

### 8.2 User-visible functions (full v0+ list, grouped by lifecycle)

- **Construct/compile:** `lmm()`, `glmm()`, `compile_model()`, `explain_model()`,
  `audit_design()`, `re()`, `vc()`, `roles()`, `mm_control()`, `mm_thresholds()`.
- **Print/summary:** `print`, `summary`.
- **Extractors (lme4-generic compatible):** `fixef`, `ranef(condVar=)`, `coef`, `VarCorr`,
  `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `nobs`, `model.frame`,
  `model.matrix(type=)`, `vcov(type=)`, `residuals(type=)`, `fitted`, `formula`, `family`,
  `terms`, `update`.
- **Audit-first verbs:** `audit`, `changes`, `parameterization`, `diagnostics`,
  `fit_status`, `optimizer_certificate`, `reproducibility`, `inference_table`,
  `estimability`, `df_for_contrast`, `fit_handle_alive`, `revive`, `random_blocks`,
  `is_singular`.
- **Reporting:** `model_report(fit, sections = "all")`,
  `reporting_table(fit, section = "all")`. These implement the
  Davies-Meteyard reporting contract in `planning/reporting_contract.md`:
  `mixeff` assembles publication-oriented tables from Rust-owned artifact
  fields plus R-owned provenance, preserving method/status/reason fields and
  marking unavailable sections explicitly.
- **Demystification & guidance (no recommendations — see §9.5):**
  - `explain_model(spec_or_fit)` — auto-invoked once by `lmm()`/`glmm()` *before*
    fitting. Translates each random term along the three axes (scope,
    covariance, support): the original formula fragment, its
    named-argument paraphrase, a one-line English gloss, what the term
    *does* model and what it explicitly *does not* model
    (`not modeled:` line), and design notes (which slope variables vary within
    group). Always emits the sentinel `No random slopes were added.` when the
    user wrote only intercept terms — the no-silent-surgery contract made
    visible. Suppress via `mm_control(verbose = -1)`.
  - `compare_covariance(spec_or_fit)` — opt-in alternate view. Side-by-side
    table of `full` / `diag` / `scalar` for each random term: parameter count,
    observations per parameter, budget threshold, status, and a plain-English
    "assumes zero" column. No recommended row.
  - `random_options(spec_or_fit, group, slope = NULL)` — opt-in *map* of nearby
    random-effect structures for a grouping factor (and optionally one
    within-group slope variable). Lists the punt, slope-only,
    additive-uncorrelated (`(1 | g) + (0 + x | g)`), the `||` synonym
    annotated as "same model, different font," and the maximal `(x | g)`.
    Each row prints varying coefficients, covariance family, theta parameters,
    design status, and a plain-English meaning. **It is a map, not a ranking**:
    no recommended column, no implicit ordering by preference, rung 0 (the
    punt) is a first-class entry.
  - **Three kinds of help.** `explain_model()` and `audit_design()` separate
    structural impossibility (firm — strict mode refuses), low information
    budget (factual, non-moralizing — reports parameters / levels / threshold /
    status), and unmodeled-but-possible structure (a quiet single-line scope
    note, never escalated to a warning). See §9.5.4.
  - **Split-term explanation.** When the formula contains a split-block
    pattern such as `(1 | subject) + (0 + time | subject)`, `explain_model()`
    prints both blocks separately, calls out the implied zero covariance
    between blocks with its reason ("separate random-effect blocks"), and
    enumerates nearby spellings (`(1 + time | subject)`,
    `(1 + time || subject)`) with what each estimates or fixes. See §9.5.2.
- **Predict/simulate:** `predict(re.form=, allow.new.levels=, type=, interval=, target=,
  condition_on=, se.fit=)`, `simulate`, `refit`.
- **Inference & comparison:** `contrast(L, rhs=, method=)`, `test_effect`, `anova`
  (single + multi), `compare`, `drop1`, `confint`.
- **Internals (read-only compat):** `getME(fit, name)` for `X`, `Z`, `theta`, `Lambda`,
  `cnms`, `flist`, `Gp`, `lower`, `devcomp`, `optinfo`.
- **emmeans (Phase 5):** `recover_data.mm_fit`, `emm_basis.mm_fit`.
- **Serialization:** `as_json(spec_or_fit)`.

### 8.3 `mm_control()` / `mm_thresholds()` shape

Flat named list, mirroring `lmerControl`. Serializes 1:1 to the Rust `CompilerPolicy` JSON.

```r
mm_control(
  optimizer = "auto",            # "auto" | "nlopt_bobyqa" | "cobyla" | "lbfgsb"
  optimizer_max_iter = 10000L,
  optimizer_xtol_abs = 1e-8,
  optimizer_ftol_abs = 1e-10,
  reml = TRUE,                   # LMM only
  nAGQ = 1L,                     # GLMM only
  verify_convergence = "bounded",  # "none"|"bounded"|"jittered"|"alt_optimizer"
  parallel_threads = 0L,
  seed = NULL,
  verbose = 0L,
  thresholds = mm_thresholds(),
  schema_version = "current",
  bridge_timeout_s = 0
)
mm_thresholds(
  min_levels_random_intercept_fit         = 2L,
  min_levels_random_intercept_reliability = 5L,
  min_levels_full_cov_intercept_const     = 10L,
  min_levels_full_cov_per_param           = 5L,
  max_condition_number                    = 1e10,
  min_within_group_sd                     = 1e-8,
  max_basis_pairwise_abs_corr             = 0.999,
  min_observations_per_supported_level    = 2L,
  effective_rank_relative_tolerance       = 1e-6,
  effective_rank_absolute_tolerance       = 1e-10
)
```

These names are byte-equivalent to `compiler_contract_v0_prd.md` §8 thresholds (lines
451–463).

## 9. Internal contract plumbing

### 9.1 `mm_fit` storage (named list)

| Field | Type | Origin |
| --- | --- | --- |
| `$call`, `$formula`, `$random_spec`, `$frame_schema`, `$mode`, `$inference_request`, `$control` | R | captured |
| `$artifact_json` | `mm_artifact` (parsed list) | Rust artifact — **truth** |
| `$state_json`, `$audit_json`, `$diagnostics_json`, `$theta_map_json`, `$cert_json`, `$inference_json`, `$reproducibility_json` | typed parsed lists | Rust |
| `$beta`, `$theta`, `$sigma`, `$logLik`, `$deviance`, `$df_residual`, `$fit_status` | small numerics / strings | duplicated for cheap access |
| `$rust_handle` | `externalptr` or `NULL` | bridge cache |
| `$schema` | `list(schema_name, schema_version, crate_version, package_version)` | header |
| `$lazy_cache` | environment | mutable, holds `X`, `Z`, `Lambda`, block decomps |

### 9.2 JSON helper layer (`R/json.R`, internal)

Per-blob parsers with explicit eager/lazy stance:

| Helper | Coercion |
| --- | --- |
| `mm_json_parse_artifact` | aggressive |
| `mm_json_parse_state` | aggressive |
| `mm_json_parse_audit` | preserve raw |
| `mm_json_parse_diagnostics` | aggressive → data.frame-ready |
| `mm_json_parse_theta_map` | preserve tag, coerce arrays |
| `mm_json_parse_certificate` | preserve raw |
| `mm_json_parse_inference` | aggressive → data.frame-ready |
| `mm_json_parse_repro` | preserve raw |
| `mm_json_parse_manifest` | aggressive |
| `mm_json_parse_prediction` | aggressive |
| `mm_json_negotiate(header)` | validates schema version, raises `mm_schema_error` |

`as_json(x)` is the inverse for `saveRDS` augmentation and for the user-facing
`as_json(spec)`.

### 9.3 Lazy extractor pattern

```r
.mm_lazy <- function(fit, key, producer) {
  cache <- fit$lazy_cache
  if (!exists(key, envir = cache, inherits = FALSE)) {
    if (!fit_handle_alive(fit)) fit <- revive(fit)
    assign(key, producer(fit$rust_handle), envir = cache)
  }
  get(key, envir = cache, inherits = FALSE)
}
```

`getME(fit, "Z")` → cached sparse `dgCMatrix`; `revive(fit)` returns a new `mm_fit` with a
fresh empty cache; `refit()` invalidates the full cache.

### 9.4 Honesty defaults

- `summary` in `mode = "exploratory"` or `"predictive"` prints
  `df = NA`, `statistic = NA`, `p.value = NA`, `method = "unavailable"` and a footer
  pointing the user to `mode = "confirmatory"`.
- `summary` on `ConvergedBoundary` is a `cli` *info*, not a warning; affected `VarCorr`
  rows tagged `[boundary]`; p-values printed only if Rust certifies a boundary-aware
  inference path.
- `predict(se.fit = TRUE)` when SE unavailable: returns `$se.fit = NA`,
  `attr(., "mm_unavailable_reason") = <Rust code>`, single `cli_inform()` per session per
  fit. `interval != "none"` without SE raises `mm_inference_unavailable` rather than
  fabricating intervals.
- **Confirmatory mode contract.** `mode = "confirmatory"` means: fit the
  *requested* scientific model if the design supports it. The compiler may
  refuse, canonicalize, or reduce unsupported requested structure (and must
  cross the boundary when it does — see `changes()`). It **never** adds
  unrequested random slopes, grouping factors, or covariance parameters. The
  phrase "maximal feasible" is reserved for `mode = "exploratory"`, where any
  future automatic random-structure discovery would live; confirmatory objects
  always reflect what the user wrote, modulo refusals and reductions
  documented in `changes()`. Documentation must avoid the unqualified phrase
  "deterministic design-compiled / maximal feasible" without naming the mode
  it applies to.
- **Singularity is reported model state, not a failure.** When the fitted
  covariance is rank-deficient, `summary()` and `print()` describe the state
  factually (requested dimension vs effective rank) and point the user to
  `changes(fit)` and `random_options(spec, group = ...)`. The package never
  prints "Try (1 | subject) instead." or "Drop the random slope." See §9.5.6.

### 9.5 Random-effects guidance without model recommendation

The package adds a *guidance* layer for random-effects syntax, **not** a
*model-selection* layer. The product position is:

> "Here is what the model you wrote can express, what it cannot express, what
> assumptions are being fixed by syntax, and what your data appear capable of
> supporting."

That fits the project mission: audit-first, no fake certainty, no silent
surgery. lme4-style formula syntax is an *on-ramp*, not the source of truth;
the semantic IR (random term: `group`, `slopes`, `cov`, `intercept`) is.

**Three axes** organize every explanation:

| Axis | User-facing question | Example |
| --- | --- | --- |
| **Scope** | Which coefficients are allowed to vary by group? | subject baselines only? subject time slopes too? |
| **Covariance** | Are those varying coefficients allowed to be related, or is some covariance fixed at zero by syntax? | do high-baseline subjects also tend to have steeper time trends? |
| **Support** | Does this dataset have enough design information to estimate that? | enough subjects? enough within-subject time variation? enough observations per subject? |

This frame replaces a "full / diag / scalar" comparison with a syntax-to-meaning
mapping that scales to the patterns users actually write.

#### 9.5.1 The contract — six rules binding on the v1 API

1. `lmm()` and `compile_model()` never add unrequested random effects, grouping
   factors, or covariance parameters.
2. `explain_model()` always translates random-effect syntax into the named
   semantic form (`re(group=, intercept=, slopes=, cov=)`).
3. `explain_model()` includes scope notes: which coefficients vary, which are
   fixed common across groups, which covariances are estimated, which are fixed
   at zero by syntax.
4. `audit_design()` reports structural impossibility and low-support
   covariance structures *before* fitting.
5. `random_options()` shows nearby random-effect structures, parameter costs,
   and data-support facts. It is a *map*, not a ranking.
6. `changes()` explains requested → effective → fitted random structure
   *after* fitting, including boundary and reduced-rank outcomes.

#### 9.5.2 Layer 1 — `explain_model()` auto-print on `lmm()` / `glmm()`

For each random term, the auto-print emits the formula fragment, its
named-argument paraphrase, a one-line English gloss, what the term *does*
model, what it *does not* model, and design notes. Default is terse (one
budget line per term plus targeted notes); suppress via `mm_control(verbose =
-1)`.

**The common punt — `(1 | subject)`.** When the user writes a fixed effect
that varies within group but does not request a corresponding random slope,
the package does not warn and does not say "you should add random slopes." It
quietly states what the model *does* and what it *does not* model:

```
Random effects
subject:
  wrote:        (1 | subject)
  named form:   re(group = subject, intercept = TRUE, slopes = NULL, cov = "scalar")
  meaning:      subjects may differ in average outcome.
  not modeled:  the effects of time and condition are constrained to be the same
                across subjects.
Design note:
  time varies within subject and could be considered as a subject-level slope.
  condition does not vary within subject, so a subject-level condition slope is
  not supported by this data layout.
No random slopes were added. Use random_options(spec, group = subject) to inspect
what alternative subject terms would mean.
```

The single sentence **`No random slopes were added.`** is mandatory. It
preserves the no-silent-surgery contract: the user knows the package
considered the question and chose to honor the syntax exactly. The package
does not treat `(1 | g)` as a "safe default"; it treats it as a specific model
statement.

**The split-term pattern — `(1 | subject) + (0 + time | subject)`.** The
folklore expression for "uncorrelated random intercept and slope." The
package explains it directly rather than leaving the user to read the lme4
source:

```
subject has two separate random-effect blocks.
Block 1:
  wrote:        (1 | subject)
  named form:   re(group = subject, intercept = TRUE, slopes = NULL, cov = "scalar")
  meaning:      subjects may differ in average outcome.
Block 2:
  wrote:        (0 + time | subject)
  named form:   re(group = subject, intercept = FALSE, slopes = time, cov = "scalar")
  meaning:      subjects may differ in their time slope.
Relationship between blocks:
  The subject intercept and subject time slope are modeled as uncorrelated because
  they are in separate covariance blocks.
Nearby spellings:
  (1 + time | subject)   estimates intercept variance, time-slope variance,
                         and their covariance.
  (1 + time || subject)  estimates intercept variance and time-slope variance,
                         but fixes their covariance to zero.
```

This direct explanation of split-block syntax is a first-class responsibility
of `explain_model()`, not a vignette-only lesson.

#### 9.5.3 Layer 2 — `random_options(spec_or_fit, group = <g>)`

Opt-in, prefit, non-fitting. Shows the *nearby* model space for a grouping
factor: the punt, slope-only, additive-uncorrelated (with its `||` synonym
annotated as "same model, different font"), and the maximal `(x | g)`. Each
row prints varying coefficients, covariance family, theta parameters, design
status, and a plain-English meaning. There is **no recommended column** and no
ordering by preference. Rung 0 (the punt) is a first-class entry, not a
starting point.

```
Random-effect options for group: subject
Current model:
  (1 | subject)                                <- this is what you wrote
  Allows subject baselines to vary.
  Does not allow the time effect to vary by subject.
Nearby options:
  (1 | subject)
    varying coefficients: intercept
    covariance family:    scalar
    theta parameters:     1
    design status:        supported
    plain meaning:        subjects differ in baseline only
  (0 + time | subject)
    varying coefficients: time
    covariance family:    scalar
    theta parameters:     1
    design status:        supported
    plain meaning:        subjects differ in time slope only
  (1 | subject) + (0 + time | subject)
    varying coefficients: intercept, time
    covariance family:    diagonal via separate blocks
    theta parameters:     2
    design status:        supported
    plain meaning:        subjects differ in baseline and time slope;
                          intercept-slope covariance fixed to zero
  (1 + time | subject)
    varying coefficients: intercept, time
    covariance family:    full
    theta parameters:     3
    design status:        low support: 14 subjects, policy floor 15
    plain meaning:        subjects differ in baseline and time slope;
                          estimate whether these are associated
```

`compare_covariance(spec_or_fit)` is a thin alternate view of the same
artifact data — three rows per term (full / diag / scalar) — for users who
prefer the column-major layout. It is *not* the primary surface;
`random_options()` is.

#### 9.5.4 Three kinds of help — each in its own register

The package distinguishes three kinds of pedagogical help. Each is rendered in
a different tone so users learn what to do with it.

| Kind | What it is | Tone | Where it surfaces |
| --- | --- | --- | --- |
| **Structural impossibility** | A requested term cannot be estimated from this data layout (e.g., between-subject `condition` requested as a subject-level slope). | Firm, factual. Strict mode refuses to fit. | `audit_design()` → `mm_design_refusal` |
| **Low information budget** | A requested term is information-hungry relative to the observed grouping levels (e.g., `(1 + time | subject)` with 8 subjects). | Factual, non-moralizing. Reports parameter count, levels, threshold, status. Does **not** claim the model is wrong. | `audit_design()` warning, `parameterization()`, `random_options()` design status |
| **Unmodeled but possible structure** | A within-group fixed effect could in principle vary by group but no random slope was requested. | Quiet, educational. A single neutral note. **Never** escalated to a warning. | `explain_model()` `Design note` line |

Sample renderings:

*Structural impossibility:*
```
The requested subject-level slope for condition is not supported:
condition does not vary within subject.
A subject-specific condition effect cannot be estimated from this design.
No model was fit in strict mode.
```

*Low information budget:*
```
The term (1 + time | subject) requests a 2x2 full covariance matrix:
  variance(intercept)
  variance(time)
  covariance(intercept, time)
Free covariance parameters: 3
Grouping levels: 8 subjects
Policy status: below reliability floor for full covariance
This does not prove the model is wrong. It means the requested covariance
structure is information-hungry relative to the observed grouping levels.
```

*Unmodeled but possible structure:*
```
No subject-level random slope was requested for time.
time varies within subject, so such a term is structurally possible.
This model treats the time effect as common across subjects.
```

#### 9.5.5 Tone of design refusals — "possible repairs," not "suggested model"

When `audit_design()` emits a `mm_design_refusal`, the rendered footer uses the
heading **`Possible repairs, not applied automatically:`** and lists each
repair with the *reason* the design supports it. The phrase "suggested
starting model" is **never** used; the package does not suggest models.

```
Possible repairs, not applied automatically:
1. Treat season as fixed
   Reason: season has 3 observed levels and appears to be a condition being
   compared.
2. Remove the full random slope block for sites
   Reason: sites has 3 levels, but (1 + duration | sites) requests 3 covariance
   parameters.
3. If sites are sampled from a larger population and the real study has many
   sites, keep site random effects in the larger dataset:
     lmm(log(effect) ~ duration * season + (1 | site), data = dat)
```

Each repair is a *coherent* alternative with a *reason*, not a ranked
recommendation.

#### 9.5.6 Singularity is reported model state, not a shameful failure

When the fitted covariance is rank-deficient, the rendered message is
descriptive and points to the audit verbs, not to a folk fix:

```
The fitted covariance matrix is rank-deficient.
The requested subject block had dimension 2; the fitted block has effective rank 1.
Use changes(fit) to see which dimension was unsupported.
Use random_options(spec, group = subject) to inspect lower-dimensional covariance
choices.
```

The package never prints `Try (1 | subject) instead.` or `Drop the random
slope.` Singularity is a state to inspect, not a flaw to apologize for.

#### 9.5.7 v1 syntax coverage for `explain_model()`

Phase 1 must produce deterministic, prefit, named-form translations for at
least these formula patterns, each with a snapshot test:

- `(1 | subject)` (the punt)
- `(0 + time | subject)` (slope-only)
- `(1 + time | subject)` (correlated intercept + slope)
- `(1 + time || subject)` (uncorrelated, double-bar form)
- `(1 | subject) + (0 + time | subject)` (uncorrelated, split-block form;
  both forms must be explicitly annotated as "same model, different font")
- `(1 | a/b)` (`syntax_expansion`: expands to `(1 | a) + (1 | a:b)`)
- `(1 | subject:item)` (interaction grouping)
- `(1 | subject) + (1 | item)` (crossed)

An interactive `re_builder()` is explicitly **deferred to v2**. v1 ships
passive explanation only — most of the value comes from passive explanation,
not interactivity.

### 9.6 Random term card — the per-term artifact

Each random term in `audit_design()` and the `CompiledModelArtifact` returns a
structured *card*. R never reconstructs semantics; it formats the card for
`explain_model()`, `parameterization()`, `random_options()`,
`compare_covariance()`, and `changes()`. Schema sketch:

```json
{
  "term_id": "re_002",
  "original_fragment": "(1 | subject) + (0 + time | subject)",
  "canonical_fragment": "(1 | subject) + (0 + time | subject)",
  "group": "subject",
  "blocks": [
    {
      "basis": ["Intercept"],
      "intercept": true,
      "slopes": [],
      "covariance": "scalar",
      "theta_parameters": 1,
      "english": "subjects may differ in average outcome"
    },
    {
      "basis": ["time"],
      "intercept": false,
      "slopes": ["time"],
      "covariance": "scalar",
      "theta_parameters": 1,
      "english": "subjects may differ in their time slope"
    }
  ],
  "implied_constraints": [
    {
      "type": "zero_covariance",
      "between": ["Intercept", "time"],
      "reason": "separate random-effect blocks"
    }
  ],
  "design_support": {
    "group_levels": 42,
    "min_rows_per_group": 3,
    "median_rows_per_group": 8,
    "within_group_variation": { "time": "present" },
    "status": "supported"
  }
}
```

`english` per-block phrasing and `implied_constraints[].reason` are authored
in Rust and shipped over the wire — single source of truth across language
bindings.

### 9.7 Pedagogical diagnostic taxonomy (non-error notes)

The existing stable `DiagnosticCode` enum (`FixedRandomRedundant`,
`RandomSlopeUnsupported`, `CovarianceTooRich`, `NotIdentifiable`, ...) covers
errors and warnings. The wrapper requires a parallel **non-error** class for
random-effects pedagogy, surfaced by `audit_design()` and `explain_model()`:

| Code | Surfaces when | Severity |
| --- | --- | --- |
| `ScopeNote` | A within-group fixed effect has no corresponding random slope. | informational |
| `SupportNote` | A requested term is below the reliability floor for its covariance family but is not refused. | informational |
| `SyntaxExpansion` | Syntactic sugar expands to a longer canonical form, e.g., `(1 | a/b)` → `(1 | a) + (1 | a:b)`. | informational |
| `CovarianceAssumption` | A formula choice fixes a covariance to zero by syntax (`||`, split blocks). | informational |
| `StructuralRefusal` | A requested random slope cannot be estimated because the slope variable does not vary within group. | error in strict mode |

These live as variants of the extended `DiagnosticCode` enum on the Rust side
(see §5.2 item 4) and serialize identically to the existing closed enum. R
formats; R does not author wording.

### 9.8 `roles()` — declared design roles, with observed fallback

`roles()` is the user's optional channel for declaring how each variable
participates in the design. v1 accepts simple strings; the conceptual role
system (potentially expanded in v2) answers:

- Is this variable a grouping unit?
- Is this variable manipulated/observed *within* the grouping unit?
- Is this variable between-group?
- Is this variable a condition being compared?
- Is this variable time-like?
- Is this variable nested or crossed?

Aspirational v2 surface:

```r
roles(
  subject   = sampled_unit(),
  item      = sampled_unit(),
  time      = repeated_measure(within = subject),
  condition = fixed_condition(within = subject, crossed_with = item)
)
```

In the absence of declared roles, the compiler infers within/between status
from the data and tags the inference on the random term card:

```json
"role_origin": { "declared_by_user": false, "observed_from_data": true }
```

This keeps the package useful without forcing users into a questionnaire, and
makes it explicit when a structural claim ("`condition` does not vary within
subject") was inferred rather than declared.

## 10. Phased roadmap

Each phase ships when `R CMD check` is clean (0 errors / 0 warnings; size-related NOTEs
allowed) **and** the listed tests pass **and** the listed vignette renders.

- **Phase 0 — extendr skeleton + manifest.** `parse_formula()` round-trip,
  `mm_formula_manifest()`, schema negotiation, one typed condition, interrupt bridge.
  Tests: `test-manifest`, `test-schema-versioning`, `test-interrupts`, `test-namespace`.
  Vignette: `intro.Rmd`.

- **Phase 1 — LMM fit, audit-first.** `lmm()` (REML default), `compile_model`,
  `explain_model` (auto-printed by `lmm()`, with the named-form paraphrase,
  scope notes, the mandatory `No random slopes were added.` sentinel for
  intercept-only random terms, the split-block explanation, and "possible
  repairs" wording for refusals — see §9.5), `audit_design` (emits
  `ScopeNote`/`SupportNote`/`SyntaxExpansion`/`CovarianceAssumption`/
  `StructuralRefusal` codes per §9.7), `compare_covariance`, `random_options`
  (the *map*, not a ranking), `roles()` (string-form v1; observed-from-data
  fallback per §9.8), `print/summary`, `audit`, `changes`, `diagnostics`,
  `fit_status`, `parameterization`, `fixef`, `ranef(condVar=FALSE)`, `coef`,
  `VarCorr`, `sigma`, `logLik`, `deviance`, `AIC`, `BIC`, `nobs`, `formula`,
  `model.frame`, basic `predict()` (no SE), `fitted`,
  `residuals(type="response")`, `as_json(spec)`. The Rust side must ship the
  random term card per §9.6 and the extended diagnostic enum per §5.2 item 4.
  Snapshot tests for the five worked examples in `compiler_contract_v0_prd.md`
  lines 907–922 plus the §9.5.7 syntax coverage list (eight formula patterns,
  each with explicit snapshot of the named-form translation, design status,
  and any pedagogical diagnostic codes emitted). Vignettes: `lmm-basics.Rmd`
  and `demystifying-formulas.Rmd` (the latter must walk through the
  `(1 | subject)` punt explanation, the split-block pattern, and the three
  kinds of help).

- **Phase 2 — saveRDS round-trip + lazy extractors.** `revive`, `fit_handle_alive`,
  `getME`, `model.matrix(type=)`, `vcov(type="fixed")`, `random_blocks`,
  `optimizer_certificate`, `inference_table`, `reproducibility`, `is_singular`, lazy cache,
  `ranef(condVar=TRUE)`, and the first `model_report()` / `reporting_table()`
  implementation for sections whose Rust artifacts are available. Reporting
  completeness depends on the upstream field checklist in
  `planning/reporting_artifact_requirements.md`; unavailable sections are
  reported with explicit reasons. Tests include cross-session restart via `callr`. Vignette:
  `saving-and-reviving.Rmd`.

- **Phase 3 — Inference (LMM).** `contrast`, `test_effect`, `df_for_contrast`,
  `estimability`, single-model `anova(type=, method=)`, Satterthwaite,
  `vcov(type="theta")`, `predict(se.fit=TRUE, interval=)`, `confint(method="wald")`,
  honesty tests. Vignette: `inference.Rmd`.

- **Phase 4 — GLMM + simulate + refit + compare.** `glmm()` (Laplace + AGQ), `simulate`,
  `refit`, `compare`, multi-model `anova`, `drop1`, parametric bootstrap. Vignette:
  `glmm.Rmd`.

  > **GLMM estimator default (decision, 2026-06-02).** `glmm()` keeps
  > `method = "pirls_profiled"` (the fast profiled PIRLS path) as its default —
  > this preserves the package's speed positioning and does not silently change
  > existing results. The profiled path is **not** `glmer()`'s estimator and its
  > coefficients differ; the native `method = "joint_laplace"` route is the
  > glmer-equivalent estimator and is certified against `glmer` within tolerance
  > (`test-glmm-joint-laplace-parity.R`: ~5e-4 fixef on cbpp/poisson). To honour
  > "no silent surgery" on the estimator choice, `glmm()` emits an informational
  > notice (class `mm_estimator_notice`) when `method` is left at its default and
  > `verbose >= 0`, pointing users to `joint_laplace` for glmer parity. We chose
  > the notice over flipping the default so the change is opt-in and visible
  > rather than a silent numerical shift for current users.

- **Phase 5 — emmeans + multivariate + profile CIs.** `recover_data.mm_fit`,
  `emm_basis.mm_fit`, `profile()`, `cbind(y1, y2) ~ ...` shared-theta multivariate.
  Vignette: `emmeans-and-multivariate.Rmd`.

## 11. Test strategy day-1

- **Schema snapshots** under `tests/testthat/_snaps/` and validated against
  `inst/schemas/*.json` using `jsonvalidate`.
- **Random-effects guidance surface** (`test-demystify.R`): snapshot tests for
  - `explain_model()` auto-print on each pattern in §9.5.7 (formula →
    named-form → English gloss, including the mandatory
    `No random slopes were added.` sentinel for intercept-only terms);
  - the split-block explanation for `(1 | subject) + (0 + time | subject)`,
    including the explicit "same model, different font" annotation versus
    `(1 + time || subject)`;
  - the three kinds of help (§9.5.4) — structural impossibility (refusal
    text, fires only when a slope variable does not vary within group), low
    information budget (factual reporting, no moralizing language, no
    "should"), unmodeled-but-possible (single-line `Design note`, never
    surfaces as a warning);
  - "possible repairs, not applied automatically" wording for design
    refusals (§9.5.5) — the strings "suggested starting model" and "we
    recommend" must not appear anywhere in package output;
  - singularity rendering (§9.5.6) — describes effective rank, points to
    `changes()` and `random_options()`, never says "Try ... instead." or
    "Drop the random slope.";
  - `compare_covariance()` table layout across `full`/`diag`/`scalar` with no
    "recommended" column;
  - `random_options()` map rendering — rung 0 is first-class, no preference
    ordering, no recommended row, current-model marker `<- this is what you
    wrote` present;
  - pedagogical diagnostic codes (§9.7) round-trip from the artifact JSON
    intact (Rust authors wording, R formats only).
- **Random term card contract** (`test-artifact-cards.R`): assert each of the
  eight §9.5.7 patterns produces a card matching the §9.6 schema, validated
  with `jsonvalidate` against `inst/schemas/random-term-card.schema.json`.
- **Reporting contract** (`test-reporting.R`): snapshot `model_report()` and
  `reporting_table()` sections defined in `planning/reporting_contract.md`.
  Assert fixed-effect rows preserve Rust `method`/`status`/`reason`, random
  terms preserve Rust-authored English and covariance constraints, grouping
  summaries include levels and min/median/max rows where available, boundary
  and reduced-rank states are reported as model state, and unavailable sections
  are explicit rather than omitted. Forbidden recommendation strings from R9
  must not appear in report output.
- **lme4 parity** on `sleepstudy`, `Pastes`, `cake` (`skip_if_not_installed("lme4")`),
  tolerances: `fixef` 1e-4, `theta` 1e-3, `logLik` 1e-3, `sigma` 1e-4. Documented expected
  divergences in `inst/extdata/expected-mismatches.json` (one row per
  `(case_id, field)` with regression-detector bounds; see
  bd-01KQF83XAN4CGAS176JRX8CR7E for the classification ledger contract). The
  legacy `expected-mismatches.yaml` stub is retained as a deprecated pointer
  and is not consulted by tests.
- **`saveRDS` revival without live handle**: serialize, force GC, deserialize, assert
  `print/audit/changes/parameterization/fixef/predict` work.
- **Cross-session restart** via `callr::r()` (`skip_on_cran()`).
- **Interrupt** test (`skip_on_cran()`).
- **No-mask** test: `lme4::lmer` is not displaced after `library(mixeff)`.

## 12. Risk register

| # | Risk | Mitigation |
| --- | --- | --- |
| R1 | `nlopt` blocks CRAN (C lib + CMake at build time) | Upstream PR to feature-gate `nlopt`; pure-Rust default routes to existing `cobyla` and `pattern_search` paths (no new optimizer dep). Decision A resolves to v0.2 CRAN once the upstream PR lands; v0.1 ships R-Universe regardless. Draft at `planning/upstream-nlopt-issue.md`. |
| R2 | Optimizer parity drift vs lme4 once nlopt is gone | Document tolerances; non-goal §3 disclaims bit-exactness. |
| R3 | Tarball size with vendored Rust | Aggressive `.Rbuildignore`, `strip` + `lto`, ≤ 10 MB target. |
| R4 | Windows UCRT toolchain mismatch | Pin `x86_64-pc-windows-gnu`; CI on `windows-latest` with Rtools43. |
| R5 | Schema drift as Rust crate evolves | Nightly `vendor-sync` PR + checksum gate + snapshot tests. |
| R6 | `ExternalPtr` finalizer ordering during `saveRDS` | Drop handle before serialize; serializer hook in `revive.R`. |
| R7 | macOS fork + Rust threads deadlock | Force single-threaded fits when forked; document. |
| R8 | Naming collisions (`compare`, `contrast`) | Namespace-only verbs; opt-in `use_mixeff_verbs()`. |
| R9 | "Advice creep" — the demystification surface drifts toward implicit recommendation (e.g., a "recommended" column, a default-sorted ladder, "should" language in scope notes) | Non-goal in §3 disclaims it. Tests assert the strings "suggested starting model", "we recommend", "you should", "try ... instead", "drop the random slope" never appear in package output. `random_options()` is contract-tested as a map (no ordering, no recommended row). All English wording for blocks/constraints lives in Rust (§9.6) so reviewers see drift in one place. |

## 13. Decisions (resolved 2026-04-27)

All four Phase-0 preconditions are resolved. The upstream `mixedmodels` PR
feature-gating `nlopt` (§5.2 item 1) remains the only open precondition; it is
not a project decision but an upstream dependency tracked in R1.

- **Decision A — CRAN is the ultimate distribution goal.** Resolved: R-Universe
  is the dev channel for v0.1; CRAN submission is targeted for v0.2 once the
  upstream `nlopt` feature-gate PR has merged. The upstream PR is committed as
  a hard precondition for CRAN but does not block v0.1 to R-Universe. R1 stays
  active in the risk register until the upstream PR lands.
- **Decision B — Package name is `mixeff`.** Resolved: matches the existing
  target directory `/Users/bbuchsbaum/code/mixeff`.
- **Decision C — S3-only is the default.** Resolved: S3 unless a strong,
  documented case for S4 emerges (e.g., Matrix/sparse method dispatch).
- **Decision D — No bit-exact lme4 numerical reproduction.** Resolved as an
  explicit non-goal (§3). Pair with the positive commitment: a *close
  interface* to lme4 — formula syntax, S3 generics (`fixef`, `ranef`,
  `VarCorr`, `getME`, etc.), familiar `print`/`summary` surface — so users
  can pick up `mixeff` quickly, then discover the audit-first surface
  (`explain_model`, `audit_design`, `changes`, `random_options`). Statistical
  equivalence within documented tolerances on parity datasets is the
  numerical bar.

## 14. Critical files (path map)

- `/Users/bbuchsbaum/code/mixeff/DESCRIPTION` — `SystemRequirements: Cargo (>= 1.78.0)`,
  `Imports: Matrix, cli, rlang, jsonlite`, `Suggests: lme4, testthat, callr, jsonvalidate,
  rextendr`.
- `/Users/bbuchsbaum/code/mixeff/src/rust/Cargo.toml` — `extendr_api` dep, vendored
  `mixedmodels` path dep, `default-features = false` for upstream's `nlopt` feature,
  `rust-version = "1.78"`.
- `/Users/bbuchsbaum/code/mixeff/src/rust/src/lib.rs` — `extendr_module! { ... }` listing
  the low-level `mm_*` entry points enumerated in `r_layer_proposal.md` lines 649–665.
- `/Users/bbuchsbaum/code/mixeff/R/json.R` — JSON parsing layer; gates the contract.
- `/Users/bbuchsbaum/code/mixeff/R/revive.R` — `saveRDS` round-trip + handle revival.
- `/Users/bbuchsbaum/code/mixeff/tools/vendor-rust.R` — CRAN reproducibility gate.
- `/Users/bbuchsbaum/code/mixeff/.github/workflows/R-CMD-check.yaml` — multi-platform
  matrix with rustup setup.
- `/Users/bbuchsbaum/code/rust/mixedmodels/Cargo.toml` — **upstream change required** to
  feature-gate `nlopt`.

## 15. Verification (end-to-end)

After Phase 1 ships, the following must succeed end-to-end on a clean checkout:

1. `Rscript tools/vendor-rust.R` — re-vendors and regenerates `LICENSE.note` cleanly.
2. `R CMD INSTALL .` — builds Rust + R, no warnings, on each of macOS arm64, Ubuntu LTS,
   Windows UCRT.
3. `R CMD check --as-cran .` — 0 errors / 0 warnings, NOTEs limited to
   first-submission/size.
4. `Rscript -e 'devtools::test()'` — all `testthat` suites green.
5. `Rscript -e 'rmarkdown::render("vignettes/lmm-basics.Rmd")'` — renders.
6. Manual smoke: `library(mixeff); fit <- lmm(Reaction ~ Days + (Days|Subject), sleepstudy);
   summary(fit); audit(fit); changes(fit); saveRDS(fit, tf <- tempfile()); rm(fit); gc();
   fit2 <- readRDS(tf); fixef(fit2); predict(fit2)` — every step works without a live Rust
   handle being assumed.
7. Schema-drift test: rebuild the `mixedmodels` upstream at HEAD, regenerate snapshots,
   and assert the diff is either empty or matches a recorded PRD-update commit.
