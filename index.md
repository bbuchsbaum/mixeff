# mixeff

> Mixed-effects models in R, with lme4-style formulas and a Rust engine.

`mixeff` fits linear and generalized linear mixed-effects models through
the [`mixeff-rs`](https://github.com/bbuchsbaum/mixeff-rs) engine. The R
interface uses familiar formulas and extractors (`fixef`, `ranef`,
`VarCorr`, `predict`, `simulate`, `anova`, `summary`, `update`, and
[`broom::tidy`](https://generics.r-lib.org/reference/tidy.html)), while
keeping the compiled model, optimizer result, and inference metadata
available for inspection.

It is not a drop-in replacement for `lme4`: use
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) or
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) rather
than `lmer()` or `glmer()`, and do not expect bit-identical estimates.
The target is statistical agreement within documented tolerances on the
package’s parity datasets.

Documentation: <https://bbuchsbaum.github.io/mixeff/>

## Why mixeff?

Three parts of the design are useful in practice:

- **Fast repeated fitting.** In the LMM scaling benchmark included with
  the package, `mixeff` was 2.2 to 5.3 times faster than `lme4` at the
  largest tested scale for five common random-effects structures. This
  matters most for bootstrap, simulation, and sensitivity analyses,
  where the same model may be refit hundreds or thousands of times.
- **The model can be inspected before it is fitted.**
  [`compile_model()`](https://bbuchsbaum.github.io/mixeff/reference/compile_model.md)
  turns the formula into an explicit model specification.
  [`audit()`](https://bbuchsbaum.github.io/mixeff/reference/audit.md)
  shows the canonical formula, random-effects structure, covariance
  parameterization, and design diagnostics before optimization starts.
- **Fit diagnostics are data, not only console text.** Convergence
  status, optimizer details, inference method, reliability, and reason
  codes live in structured objects. A fitted model is serializable with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and can be
  inspected after reloading without a live Rust handle.

These features do not remove the need to evaluate a model statistically.
They make the computation and the software’s decisions easier to
examine.

## Installation

From R-Universe:

``` r

install.packages(
  "mixeff",
  repos = c("https://bbuchsbaum.r-universe.dev", getOption("repos"))
)
```

From GitHub (requires Rust 1.78 or newer and `rextendr`):

``` r

remotes::install_github("bbuchsbaum/mixeff")
```

## A short example

Compile and inspect a model before fitting it:

``` r

library(mixeff)

spec <- compile_model(
  Reaction ~ Days + (Days | Subject),
  lme4::sleepstudy
)
audit(spec)
#> Audit Summary:
#>   overall [OK]: clean: no warnings or attention items
#>   attention [OK]: no warnings or unchecked inference-critical items
#>
#> Requested Model:
#>   formula [INFO]: Reaction ~ 1 + Days + (1 + Days | Subject)
#>   model kind [INFO]: linear_mixed_model
#>   distribution/link [INFO]: gaussian/identity
#>   objective [INFO]: exact_gaussian
#>   convergence certificate [INFO]: exact_objective
#>   fixed terms [INFO]: 1, Days
#>   random terms [INFO]: 1
#>   covariance parameter maps [INFO]: 1 map(s)
```

Then fit the model with the Rust engine and use the usual R extractors:

``` r

fit <- lmm(
  Reaction ~ Days + (Days | Subject),
  data = lme4::sleepstudy,
  control = mm_control(verbose = -1)
)

summary(fit)
fixef(fit)
VarCorr(fit)
ranef(fit)
```

The optimizer certificate is programmatic. Here is a compact view of the
sleep-study fit:

``` r

cert <- optimizer_certificate(fit)$table
subset(cert, metric %in% c("status", "optimizer", "iterations"))
#>      metric              value
#>      status converged_interior
#>   optimizer           trust_bq
#>  iterations                457
```

Other inspection functions expose different parts of the same fitted
artifact:

``` r

changes(fit)                 # changes made while compiling the request
diagnostics(fit)             # structured warnings and status information
parameterization(fit)        # random-effects covariance parameterization
inference_table(fit)         # method, status, and reliability by coefficient
reproducibility(fit)         # schema and engine metadata
```

Because the fitted quantities are stored in the R object, ordinary R
serialization works:

``` r

saveRDS(fit, path <- tempfile())
restored <- readRDS(path)
fixef(restored)
optimizer_certificate(restored)
```

`revive(restored)` resets the process-local lazy R-side cache after
serialization; the current bridge intentionally leaves the Rust handle
absent.

## Performance

The committed scaling benchmark fits the same LMM with
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
and [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) three times
per cell. The table below reports the largest cell for each design.

| design | largest tested scale | mixeff median | lme4 median | speedup |
|----|---:|---:|---:|---:|
| random intercept, varying rows | 5,000 rows | 6 ms | 13 ms | 2.2x |
| random intercept, varying groups | 200 subjects | 3 ms | 9 ms | 3.0x |
| correlated random slope | 200 subjects | 5 ms | 17 ms | 3.4x |
| crossed random intercepts | 30 subjects and 30 items | 5 ms | 18 ms | 3.6x |
| crossed design with random slope | 30 subjects and 30 items | 7 ms | 37 ms | 5.3x |

These are small absolute timings from one benchmark run, with only three
replications per cell. They show the behavior of this harness, not a
universal speed guarantee. The scripts and full CSV are included so the
comparison can be rerun on a relevant machine and model:

- `inst/benchmarks/lme4-scaling.R`
- `inst/extdata/lme4-scaling-summary.csv`
- the [benchmarking
  vignette](https://bbuchsbaum.github.io/mixeff/articles/benchmarking.html)

## Inspecting random-effects structure

[`explain_model()`](https://bbuchsbaum.github.io/mixeff/reference/explain_model.md)
describes what a random-effects formula estimates and what the design
can support. It does not rank or replace candidate models.

``` r

explain_model(compile_model(score ~ week + (1 | clinic), df))
#> Random effects:
#>   clinic:
#>     wrote:        (1 | clinic)
#>     named form:   re(group = clinic, intercept = TRUE, slopes = NULL, cov = "scalar")
#>     scope:        clinics may differ in average outcome.
#> Design notes:
#>   scope_note: week varies within clinic and could be a clinic-level slope.
```

Split-block, double-bar, and nested formulas are expanded explicitly.
For example, `(1 | a/b)` is shown as `(1 | a) + (1 | a:b)` and labelled
as a syntax expansion.

## Numerical compatibility with lme4

`mixeff` does not target bit-exact reproduction of `lme4`. Its bundled
Rust engine and optimizer can take a different numerical path while
arriving at a statistically equivalent result. Expected differences on
the parity datasets are classified in
`inst/extdata/expected-mismatches.json`, with tolerances enforced by the
test suite.

The current GLMM surface supports binomial, Poisson, Gamma, and negative
binomial models for documented links. The default profiled PIRLS
estimator is fast but is not the same estimator as `glmer()`; use
`method = "joint_laplace"` when glmer-equivalent Laplace estimates are
needed for a supported family. See the [GLMM
vignette](https://bbuchsbaum.github.io/mixeff/articles/glmm.html) for
the current boundaries.

## Acknowledgements

The `mixeff-rs` engine is modelled on Julia’s
[`MixedModels.jl`](https://juliastats.org/MixedModels.jl/). Its staged
formula-to-fit design—formula parser, semantic representation,
covariance parameterization, optimizer, and inference contract—informed
the corresponding parts of `mixeff`.

## License

MIT, plus the upstream Rust crate license bundle in `inst/LICENSE.note`.
