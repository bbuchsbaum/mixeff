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
