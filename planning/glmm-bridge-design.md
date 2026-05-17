# GLMM Bridge — Implementation Design

Status: design only (no code). Tracking bead: `bd-01KRCKCYZ51H7H2WN8C5D7FNGT`
(Stage B.1: `mm_fit_glmm_json` FFI). Produced by the 7-agent deep-dive audit,
remediation bead `bd-01KRV13QPQ6CKKDZ2PJB8A2SF9`, 2026-05-17.

## Problem

The pinned Rust crate fully implements GLMMs, but the R package cannot reach
that capability. Concretely:

- `src/rust/src/lib.rs` `extendr_module!` (lines ~1169–1187) registers
  `mm_compile_model_json`, `mm_fit_lmm_json`, the bootstrap/contrast
  entry points, and audit helpers — **but no GLMM fit primitive**.
- `R/glmm.R` validates `family`/`link`/`nAGQ` and then **unconditionally
  aborts** (`R/glmm.R:71`, `mm_fit_error`, "GLMM fitting is not available in
  this mixeff build because the Rust bridge does not expose a GLMM fit
  primitive yet.").
- The capability manifest honestly reports `fit_glmm = FALSE`
  (`src/rust/src/lib.rs:201`); `tests/testthat/test-manifest.R:79` and
  `tests/testthat/test-phase4.R` lock that in.

This is **not** a pin-staleness problem. The GLMM engine is byte-identical
between the pinned snapshot and `origin/main`; bumping the pin changes
nothing here. The missing piece is entirely in the **mixeff-owned
hand-written extendr glue** (`src/rust/src/lib.rs`), which has no upstream
and must be written in this repository.

### Path corrections (the tracking bead's body is stale)

The bead body references `/Users/bbuchsbaum/code/rust/mixedmodels/...`. That
path/name is obsolete. Correct locations:

- Live peer crate: `/Users/bbuchsbaum/code/rust/mixeff-rs`
- Vendored pinned snapshot (what actually compiles):
  `src/rust/upstream/mixeff-rs/` (pin `PINNED_REV` in `tools/vendor-rust.R`,
  currently `59b2fb3a269356b6ea62bf383adbeba95091fc31`)
- GLMM engine: `src/rust/upstream/mixeff-rs/src/model/generalized.rs`
- Engine contract: `/Users/bbuchsbaum/code/rust/mixeff-rs/docs/glmm_support_contract.md`
- A stray pre-rename `src/rust/upstream/mixedmodels/` directory was found and
  is now wiped on every re-vendor (hardened `tools/vendor-rust.R`). Do not
  reference it.

## Engine surface to wrap (pinned `generalized.rs`)

- `GeneralizedLinearMixedModelBuilder::new(formula, data, family)` (line ~137)
  → `.fit()` (line ~195) returning `GeneralizedLinearMixedModel`.
- Direct constructors: `new`, `new_with_weights`, `new_with_offset`,
  `new_with_weights_and_offset` (lines ~228–280); fitting via
  `fit(&mut self)` (line ~1116) and `fit_with_options(...)` honoring the
  `fast` contract (line ~1167).
- `Family`: Bernoulli, Binomial (with trial weights), Poisson, Gamma.
  `LinkFunction`: Logit, Probit, Cloglog (binomial family); Log, Sqrt
  (Poisson); Log (Gamma). `InverseGaussian` / Normal-as-GLMM exist but are
  **uncertified for 1.0** — keep refused.
- AGQ: `n_agq <= 1` ⇒ Laplace; `n_agq > 1` ⇒ adaptive Gauss–Hermite,
  accepted only for exactly one scalar random-effects term
  (`validate_agq`).
- `fast = true` is the supported mode; `fast = false` must surface an
  explicit unsupported error, never a silent algorithm switch.

## FFI design — `mm_fit_glmm_json`

Mirror `mm_fit_lmm_json` (`src/rust/src/lib.rs:314`) exactly in structure so
the R bridge, error tagging, and artifact handling stay uniform.

```rust
/// @noRd
#[extendr]
fn mm_fit_glmm_json(
    formula: &str,
    family: &str,           // "bernoulli" | "binomial" | "poisson" | "gamma"
    link: &str,             // "logit" | "probit" | "cloglog" | "log" | "sqrt"
    n_agq: i32,             // <=1 Laplace; >1 AGQ (single scalar RE)
    fast: bool,             // true supported; false -> explicit unsupported
    column_order: Strings,
    numeric_columns: List,
    categorical_values: List,
    categorical_levels: List,
    weights: Doubles,       // trial/case weights; empty => none
    offset: Doubles,        // fixed linear-predictor offset; empty => none
    control_json: &str,
) -> std::result::Result<String, String>
```

Behavior:

1. Parse `control_json` (same as LMM path).
2. `parse_formula(formula)` → `mm_formula_error:` on failure.
3. `data::build_dataframe(...)` (reuse existing helper).
4. Map `family`/`link` strings → engine `Family`/`LinkFunction`; unknown or
   uncertified combinations → `mm_fit_error:` (or a dedicated
   `mm_glmm_unsupported:` tag routed to `mm_arg_error`/`mm_fit_error` in R).
5. Build via `GeneralizedLinearMixedModelBuilder::new(...)`, threading
   `weights`/`offset` through the appropriate constructor; set `n_agq` and
   `fast` per contract; reject `fast = false` explicitly.
6. `.fit()`; on success serialize the post-fit `CompiledModelArtifact`
   **plus the same flat extractor-friendly numeric duplicates**
   `mm_fit_lmm_json` emits, under a new result schema
   `mixeff.glmm_fit_result/1` (versioned; negotiate via the existing schema
   machinery — add to the known-schemas table).
7. Register `fn mm_fit_glmm_json;` in `extendr_module!`; regenerate
   `R/extendr-wrappers.R` via `rextendr::document()` (never hand-edit).

Error contract: every refusal crosses as a tagged string the R side routes
to a typed condition (`mm_formula_error:`, `mm_data_error:`,
`mm_fit_error:`). No silent surgery — an unsupported `fast=false`,
family/link, or AGQ request must be an explicit typed refusal.

## R-side changes (`R/glmm.R`)

- Replace the unconditional `mm_abort` stub with a real call to
  `mm_fit_glmm_json(...)`, parsing the artifact exactly as `lmm()` parses
  `mm_fit_lmm_json` (reuse the shared artifact/handle/lazy-cache plumbing;
  return an object with the GLMM analogue of the `mm_lmm` S3 surface).
- Extend `mm_glmm_family_info()` allow-list to include `gamma = c("log")`
  (currently only `binomial = c("logit","probit","cloglog")`,
  `poisson = c("log","sqrt")`); fold Bernoulli through `binomial()`.
- Stop blanket-rejecting `weights`/`offset` (`R/glmm.R:50`): thread them to
  the FFI (binomial trial weights are mandatory per contract §16).
- Keep `nAGQ`/`approximation` validation (already `mm_arg_error` after the
  audit remediation); keep the Inverse-Gaussian/Normal refusal.
- Preserve the typed-refusal contract for any still-unsupported request.

## Manifest, tests, docs (must change in lockstep)

- `src/rust/src/lib.rs:201`: `fit_glmm = TRUE` once the primitive lands.
- `tests/testthat/test-manifest.R:79`, `tests/testthat/test-phase4.R`
  (the `expect_false(cap$fit_glmm)` / "GLMM fitting is not available"
  assertions and the `test-brown-2021` removed-glmer note) flip to assert
  the working capability.
- New tests: a real binomial-logit fit; Poisson-log; Gamma-log; AGQ vs
  Laplace; weights/offset; the `fast=false` explicit-refusal; and
  cross-engine parity vs `lme4::glmer`.
  - **Parity fixtures are not vendored.** `tools/vendor-rust.R` whitelists
    compile-time inputs only and excludes the crate's `tests/`. The R
    package must carry its own fixtures generated/copied from the peer:
    `cbpp_glmm_artifact_v1.json`, `gamma_glmm_engines.json`,
    `glmm_link_lme4.json`, `cbpp_agq5.json`, `glmm_fast_oracles.json`
    (under `/Users/bbuchsbaum/code/rust/mixeff-rs/tests/fixtures/`). Decide a
    fixtures-provenance approach (copy + provenance file, mirroring the
    existing `tests/fixtures/lme4_baseline_cases.json` pattern).
- `vignettes/glmm.Rmd`: replace the `error = TRUE` aspirational chunk with a
  real binomial GLMM fit; keep the prose accurate about certified families.

## Non-goals / related work

- `ranef(condVar = TRUE)` postVar from Rust `cond_var()` is a **separate**
  bead: `bd-01KRCKCZMYRFVEH39RFR72ZW3F` (Stage C.2). Not in this scope.
- `predict(newdata=)` / `predict_new()` exposure is separate
  (`bd-01KRCKCZJ5B5AQS5BV77VMM8ZF`, Stage C.1).
- No pin bump is required for GLMM (engine already in the pin).

## Effort estimate

Rust FFI: ~120–200 lines mirroring `mm_fit_lmm_json` + family/link mapping.
R: ~80–150 lines (rewrite `glmm()` + artifact wiring + allow-list +
weights/offset). Tests + fixtures: the larger share. Plus one `cargo`
recompile cycle. Single focused milestone; coordinate the manifest/test flip
as one atomic change so the suite never advertises a capability it cannot
reach (the audit's "no lying manifest" principle).
