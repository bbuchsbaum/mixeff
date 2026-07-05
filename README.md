# mixeff

> **Audit-first mixed-effects models in R, powered by the `mixeff-rs` Rust crate.**

`mixeff` fits linear and generalized linear mixed-effects models from
lme4-style formulas. It aims to be *functionally equivalent* to
`lme4`: the formula syntax is the same, the extractor surface
(`fixef`, `ranef`, `VarCorr`, `predict`, `simulate`, `anova`,
`summary`, `update`, `broom::tidy`) is the same, and statistical
answers agree within documented tolerances on the parity datasets
shipped with the package. It is not a literal *drop-in* replacement,
by design: you call `lmm()` / `glmm()` (not `lmer()` / `glmer()`),
results are not bit-exact, and the package is audit-first — it
refuses or reports rather than silently transforming a model.

What it adds is *provenance*. Every printed claim — a coefficient,
a standard error, a variance component, a p-value — traces back to
a versioned JSON artifact produced by a named Rust compiler at a
known schema version. It is the package to reach for when you need
to *defend* a mixed-model analysis, not just run one.

Documentation: <https://bbuchsbaum.github.io/mixeff/>

The R surface exposes the engine's formula parser, semantic IR,
design auditor, ThetaMap parameterization, optimizer, and inference
contract as first-class verbs.

## Why a different package?

| Problem in current R practice                          | What `mixeff` does                                                                 |
| ------------------------------------------------------ | ---------------------------------------------------------------------------------- |
| Convergence warnings that scroll off-screen            | `optimizer_certificate(fit)` is structured; status, objective, iterations are fields, not text |
| Singular fits printed without ceremony                 | Singularity is reported model state with effective rank, not a shameful failure    |
| p-values from methods the software never names         | Each inference row carries `method`, `status`, `reliability`, and a stable `reason_code` |
| Refusals (non-identifiable design, unsupported slope) buried in warnings | `audit_design()` raises structured `mm_design_refusal` *before* the fit |
| Reproducibility tied to a live optimizer state         | The fitted object is a serializable record; `saveRDS()` survives session restarts and reloads via `revive()` |

## Two-line install

R-Universe (recommended once 0.1.0 ships):

```r
install.packages("mixeff", repos = c("https://bbuchsbaum.r-universe.dev", getOption("repos")))
```

Development install (needs Rust toolchain ≥ 1.78 and `rextendr`):

```r
remotes::install_github("bbuchsbaum/mixeff")
```

## A six-line tour

```r
library(mixeff)

fit <- lmm(Reaction ~ Days + (Days | Subject), data = lme4::sleepstudy)
summary(fit)        # familiar lme4-style summary
audit(fit)          # Rust-authored audit report
changes(fit)        # what the compiler did to the requested model
fixef(fit); VarCorr(fit); ranef(fit)

saveRDS(fit, tf <- tempfile()); rm(fit); gc()
fit2 <- readRDS(tf)
fixef(fit2)         # works without a live Rust handle — artifact is the source of truth
```

## Audit-first workflow

`mixeff` exposes the contract as first-class verbs:

```r
spec <- compile_model(Reaction ~ Days + (Days | Subject), lme4::sleepstudy)

audit_design(spec)        # structured design audit, before any optimization
explain_model(spec)       # named-form translation of every random term
random_options(spec, Subject)   # map of nearby random-effect spellings (not a ranking)
compare_covariance(spec)  # full / diagonal / scalar comparison per term
```

Once fit:

```r
diagnostics(fit)          # structured diagnostics list
parameterization(fit)     # ThetaMap details
optimizer_certificate(fit)# convergence status, iterations, objective trace
inference_table(fit)      # per-coefficient inference with method + reliability
estimability(fit, L)      # certificate-backed estimability of contrasts
```

Where the engine cannot certify a number, the wrapper returns `NA` with
a stable reason code — never a fabricated value.

## Random-effects guidance, never recommendation

`mixeff` adds a *guidance* layer for random-effects syntax. It explains
what each formula spelling actually estimates, what the data can
support, and which covariances the syntax fixes at zero — but it never
ranks, recommends, or substitutes models.

```r
# What does this formula actually model?
explain_model(compile_model(score ~ week + (1 | clinic), df))
#> Random effects:
#>   clinic:
#>     wrote:        (1 | clinic)
#>     named form:   re(group = clinic, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:        clinics may differ in average outcome.
#> Design notes:
#>   scope_note: week varies within clinic and could be a clinic-level slope.
```

The split-block, double-bar, and nested forms are all explained
explicitly. The `(1 | a/b)` expansion to `(1 | a) + (1 | a:b)` is
labeled as `syntax_expansion`, not silently rewritten.

## Faster than lme4 on the parity benchmark

On the `mixeff` lme4-scaling benchmark (3 reps per cell, harness and
CSV under `benchmarks/lme4-scaling/` in the source tree):

| scenario             | scale            | mixeff median | lme4 median | speedup |
| -------------------- | ---------------- | ------------: | ----------: | ------: |
| <code>(1 &#124; subject)</code> rows         | 5000 rows        | 6 ms          | 13 ms       | 2.2×    |
| <code>(1 &#124; subject)</code> levels       | 200 subjects     | 3 ms          | 9 ms        | 3.0×    |
| <code>(1 + x &#124; subject)</code> slopes   | 200 subjects     | 5 ms          | 17 ms       | 3.4×    |
| <code>(1 &#124; subject) + (1 &#124; item)</code> crossed | 30 each | 5 ms          | 18 ms       | 3.6×    |
| <code>(1 + x &#124; subject) + (1 &#124; item)</code> crossed slope | 30 each | 7 ms | 37 ms | 5.3× |

(Median seconds per fit; full table at
`benchmarks/lme4-scaling/lme4-scaling-summary.csv` and an installed
copy at `inst/extdata/lme4-scaling-summary.csv` for the benchmarking
vignette to plot.)

## Numerical parity with lme4

`mixeff` does **not** target bit-exact reproduction of `lme4`. The
upstream Rust compiler defaults to a pure-Rust optimizer (cobyla /
pattern_search) rather than the C `nlopt` library that lme4 uses for
CRAN compatibility reasons.

Statistical equivalence within documented tolerances on parity datasets
is the bar. Every divergence from `lme4` is classified and bounded in
`inst/extdata/expected-mismatches.json`,
with regression-detector limits enforced by the test suite.

## Status

- **Phase 0** — extendr bridge, vendoring, formula round-trip, schema negotiation: **shipped**.
- **Phase 1** — `lmm()`, audit-first surface (`compile_model`, `audit_design`, `explain_model`, `random_options`, `compare_covariance`, `changes`, `diagnostics`, `parameterization`, `roles`), lme4-style extractors: **shipped**.
- **Phase 2** — `saveRDS` round-trip, `revive()`, lazy extractors, `model.matrix`, `vcov`: **shipped**.
- **Phase 3** — `contrast`, `test_effect`, `inference_table`, `confint`, `anova`, `drop1`, parametric bootstrap, bootstrap LRT: **shipped**.
- **Phase 4** — `glmm()` profiled-PIRLS bridge, `simulate`, `refit`, `compare`: **in progress** (joint Laplace / AGQ still explicit boundaries).
- **Phase 5** — `emmeans` integration, profile-likelihood CIs, multivariate shared-theta: deferred.

The audit-first contract is what the package is for. Read
`planning/PRD.md` and `planning/mission.md` before contributing.

## Acknowledgements

The `mixeff-rs` Rust crate that powers `mixeff` is itself modelled
on Julia's [`MixedModels.jl`](https://juliastats.org/MixedModels.jl/),
whose formula-to-fit pipeline — formula parser, semantic IR,
ThetaMap parameterization, optimizer, and inference contract — is
the basis for the corresponding stages here.

## License

MIT, plus the upstream Rust crate license bundle in `inst/LICENSE.note`.
