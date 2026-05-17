# mixeff examples

This folder contains runnable examples for fitting linear mixed models with
`mixeff::lmm()` on standard `lme4` datasets and comparing the same models
against `lme4::lmer()`.

## Standard lme4 dataset parity

Run the main example from the package root:

```sh
Rscript examples/lme4-standard-datasets.R
```

The script fits these baseline LMMs:

- `sleepstudy`: random intercept and random intercept/slope models
- `Dyestuff`: balanced one-way random-effects model
- `Dyestuff2`: boundary/singular random-intercept model
- `Penicillin`: crossed random intercepts
- `Pastes`: two random-intercept grouping factors
- `cake`: fixed interaction with a recipe-by-replicate grouping factor

For each model, the script prints fixed-effect differences, scalar fit
statistics, variance-component tables, and random-effect mode summaries where
row and column labels align.

The script also times each fit with both engines. By default it runs three
timing repetitions per engine and reports median, minimum, and maximum elapsed
seconds plus the median mixeff/lme4 ratio. Override the repeat count with:

```sh
MIXEFF_EXAMPLE_BENCHMARK_REPS=10 Rscript examples/lme4-standard-datasets.R
```

## Scaling benchmark

For plot-ready performance curves against `lme4::lmer()` across row counts,
grouping levels, crossed random effects, and random slopes, run:

```sh
Rscript inst/benchmarks/lme4-scaling.R
```

The benchmark writes raw timings, summaries, and ggplot2 PDF curves under
`benchmarks/lme4-scaling/` by default. See
`inst/benchmarks/README.md` for larger grids and scenario-selection options.

## SDAMR mixed-model companion examples

Chapter 9 of Maarten Speekenbrink's SDAM R companion includes practical
`lme4::lmer()` examples on `sdamr::anchoring` and `sdamr::speeddate`. Run the
mixeff companion script from the package root with:

```sh
Rscript examples/sdamr-lme4-companion.R
```

The script loads `sdamr` when installed and otherwise uses the vendored test
fixtures. It fits the anchoring random-intercept, correlated random-slope,
independent random-slope, and likelihood-ratio examples with `mixeff::lmm()`
next to `lme4::lmer()`. The crossed `speeddate` model remains opt-in because it
is one of the known slow parity fixtures:

```sh
MIXEFF_EXAMPLE_RUN_SLOW=true Rscript examples/sdamr-lme4-companion.R
```
