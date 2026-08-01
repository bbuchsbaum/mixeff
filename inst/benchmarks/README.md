# mixeff benchmarks

The pkgdown article `vignette("benchmarking", package = "mixeff")` explains how
to use these scripts and how to read the CSV outputs.

## mixeff/lme4 scaling benchmark

`lme4-scaling.R` benchmarks `mixeff::lmm()` against `lme4::lmer()` on synthetic
linear mixed-model designs that vary row count, number of grouping levels,
crossed random effects, and random slopes.

Run from the package root:

```sh
Rscript inst/benchmarks/lme4-scaling.R
```

The script writes plot-ready CSV files and ggplot2 PDF curves to
`benchmarks/lme4-scaling/` by default:

- `lme4-scaling-raw.csv`: one row per engine, design, and timing repetition.
- `lme4-scaling-summary.csv`: median elapsed seconds, fits per second, and
  speed relative to `lme4`.
- `lme4-scaling-fits-per-sec.pdf`: speed curves where higher is faster.
- `lme4-scaling-median-sec.pdf`: elapsed-time curves where lower is faster.
- `lme4-scaling-speedup-vs-lme4.pdf`: relative speed curves where values above
  1 mean `mixeff` is faster than `lme4`.

Useful options:

```sh
Rscript inst/benchmarks/lme4-scaling.R \
  --reps=5 \
  --rows=1000,5000,10000,25000 \
  --groups=50,100,250,500 \
  --crossed-levels=20,40,60 \
  --out=benchmarks/lme4-scaling-large
```

Use `--scenarios=rows,slopes,crossed_slope` to run a subset, or `--no-plots`
when only CSV output is needed.

## Bootstrap inference benchmark

`bootstrap-inference.R` times the bootstrap routes that support the
"audit-then-bootstrap" inference story:

- `mixeff::test_effect(method = "bootstrap")`
- `confint(method = "bootstrap")`
- `mixeff::compare(method = "bootstrap")`
- `lme4::bootMer()` fixed-effect bootstrap distribution
- `pbkrtest::PBmodcomp()` when `pbkrtest` is installed

Run from the package root:

```sh
Rscript inst/benchmarks/bootstrap-inference.R --nsim=200 --reps=3
```

The script writes:

- `benchmarks/bootstrap-inference/bootstrap-inference-raw.csv`
- `benchmarks/bootstrap-inference/bootstrap-inference-summary.csv`

Use larger `--nsim` and `--reps` values for publishable timing numbers. The
defaults are intentionally small enough for local regression checks.

## inference-method-simulation.R

Calibration study for the inference route table (Type I error, power at a
small alternative, interval coverage), with Wilson 95% intervals and full
replicate accounting (requested / evaluated / refused / failed / boundary —
nothing silently dropped). Covers four Gaussian LMM regimes (interior,
boundary, reduced-rank, small-group) across the analytic and bootstrap
routes, plus a binomial GLMM fixture with the certified joint-Laplace Wald
route, the opt-in working-Hessian approximation, and the GLMM parametric
bootstrap.

Two modes with different epistemic status:

- `--mode=smoke` (default): 3 replications; proves the pipeline executes.
  NOT calibration evidence, and the vignette refuses to present it as such.
- `--mode=release`: 2,500 replications per analytic route, 500 per bootstrap
  route, 1,000 for the joint-GLMM route. This is the mode that produces the
  shipped artifact.

The shipped artifact is regenerated with exactly:

```sh
Rscript inst/benchmarks/inference-method-simulation.R --mode=release
```

The script writes:

- `inst/extdata/inference-method-simulation-summary.csv`
- `inst/extdata/inference-method-simulation-manifest.json` — provenance
  (package commit, engine pin, lockfile checksum, seed, exact invocation,
  elapsed time); the vignette displays it next to the results, and the
  invocation recorded there reproduces the CSV bit-for-bit.
