# Benchmarking mixeff

``` r

library(mixeff)
```

Speed matters most when a workflow repeats the same fit many times.
Parametric bootstrap is the canonical case: a 1,000-replicate bootstrap
means roughly 1,000 refits, plus the original fit and summary work. The
same logic applies to cluster bootstrap, simulation studies, multi-seed
reanalysis, and anything whose inner loop is “fit the same model again
with new data”. A small per-refit speedup compounds into the difference
between an overnight run and a coffee break.

`mixeff` is roughly two to five times faster than `lme4` on the parity
benchmark — enough to make bootstrap a routine option rather than a
special-occasion calculation. This vignette shows the figures from a
single run of that benchmark and then documents the harness used to make
them. The numbers themselves are machine- and model-dependent, so the
scripts also write CSV files that you can re-plot and rerun on the
designs you actually care about.

## What does the benchmark show?

The benchmark sweeps five mixed-model designs — varying rows, varying
subject levels, adding random slopes, crossing subject and item, and
combining crossed grouping with a random slope — across a log-spaced
grid of scale values. At each cell,
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
and [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) fit the
same model three times. The summary table records median elapsed seconds
and speedup over `lme4` per cell.

``` r

benchmark_path <- system.file(
  "extdata", "lme4-scaling-summary.csv",
  package = "mixeff"
)
benchmark <- read.csv(benchmark_path, stringsAsFactors = FALSE)

benchmark$scenario_label <- factor(
  benchmark$scenario,
  levels = c("rows", "groups", "slopes", "crossed", "crossed_slope"),
  labels = c(
    "vary rows (1 | subject)",
    "vary subject levels (1 | subject)",
    "add random slope (1 + x | subject)",
    "crossed (1 | subject) + (1 | item)",
    "crossed + slope (1 + x | subject) + (1 | item)"
  )
)

head(benchmark[, c("scenario", "scale_value", "engine",
                   "median_sec", "speedup_vs_lme4")])
#>   scenario scale_value engine median_sec speedup_vs_lme4
#> 1  crossed          10   lme4      0.008        1.000000
#> 2  crossed          10 mixeff      0.002        4.000000
#> 3  crossed          15   lme4      0.010        1.000000
#> 4  crossed          15 mixeff      0.003        3.333333
#> 5  crossed          20   lme4      0.011        1.000000
#> 6  crossed          20 mixeff      0.004        2.750000
```

The first figure plots the median fit time against design scale, on
log-log axes, with one panel per scenario. `mixeff` runs below `lme4`
everywhere on the grid; the gap widens as design complexity grows.

``` r

library(ggplot2)

ggplot(benchmark,
       aes(scale_value, median_sec * 1000, colour = engine,
           group = engine)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ scenario_label, scales = "free_x", ncol = 2) +
  scale_x_log10() +
  scale_y_log10() +
  scale_colour_manual(values = c(lme4 = "#888888", mixeff = "#3366cc")) +
  labs(
    x = "design scale (rows, subject levels, etc.)",
    y = "median fit time (ms, log scale)",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())
```

![Median fit time vs design scale, log-log, by scenario, for mixeff and
lme4.](benchmarking_files/figure-html/fit-time-plot-1.png)

The second figure folds the same data into one panel: the speedup of
`mixeff` over `lme4` at each cell. The dashed reference line at 1× would
mean the two engines are equally fast.

``` r

mixeff_rows <- subset(benchmark, engine == "mixeff")

ggplot(mixeff_rows,
       aes(scale_value, speedup_vs_lme4,
           colour = scenario_label, group = scenario_label)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  scale_y_continuous(breaks = seq(0, 8, by = 1)) +
  labs(
    x = "design scale",
    y = "speedup over lme4 (×)",
    colour = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        legend.direction = "vertical",
        panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(ncol = 1))
```

![Speedup of mixeff over lme4 vs design scale, by
scenario.](benchmarking_files/figure-html/speedup-plot-1.png)

Read the speedup figure with the caveats it deserves. Three replications
per cell give a stable median but a wide envelope; the absolute timings
are sub-15 ms throughout, so an extra millisecond of system noise moves
the ratio visibly. The qualitative pattern is robust: `mixeff` is faster
everywhere, and the speedup grows with the design’s complexity — exactly
the regime where bootstrap workflows spend most of their time.

## What should you measure?

There are two separate speed questions.

First, measure model fitting. This asks how quickly each engine solves
the same mixed model as the number of rows, groups, slopes, or crossed
random effects changes.

Second, measure bootstrap routes. This asks how long the user-facing
inference verbs take when they repeatedly refit models.

``` r

data.frame(
  bootstrap_replicates = c(50L, 200L, 1000L),
  approximate_refits = c(50L, 200L, 1000L)
)
#>   bootstrap_replicates approximate_refits
#> 1                   50                 50
#> 2                  200                200
#> 3                 1000               1000
```

The arithmetic is simple, but it is the core product point: small
per-refit speed differences compound quickly when `nsim` is large.

## How do you benchmark fitting?

Use the scaling benchmark when you want to compare
[`mixeff::lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)
with [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) over
several synthetic LMM designs.

``` r

system2(
  "Rscript",
  c(
    "inst/benchmarks/lme4-scaling.R",
    "--reps=5",
    "--scenarios=rows,groups,slopes,crossed",
    "--no-plots"
  )
)
```

That script writes one row per engine, design, scale value, and timing
repetition, plus a summary table with median elapsed seconds and speed
relative to `lme4`.

``` r

data.frame(
  file = c(
    "benchmarks/lme4-scaling/lme4-scaling-raw.csv",
    "benchmarks/lme4-scaling/lme4-scaling-summary.csv"
  ),
  contents = c(
    "one timing row per engine/design/repetition",
    "median seconds, fits per second, and relative speed"
  )
)
#>                                               file
#> 1     benchmarks/lme4-scaling/lme4-scaling-raw.csv
#> 2 benchmarks/lme4-scaling/lme4-scaling-summary.csv
#>                                              contents
#> 1         one timing row per engine/design/repetition
#> 2 median seconds, fits per second, and relative speed
```

## How do you benchmark bootstrap?

Use the bootstrap benchmark when you want to time the inference routes
that make the audit-then-bootstrap story practical.

``` r

system2(
  "Rscript",
  c(
    "inst/benchmarks/bootstrap-inference.R",
    "--nsim=200",
    "--reps=3"
  )
)
```

The benchmark includes `test_effect(method = "bootstrap")`,
`confint(method = "bootstrap")`, `compare(method = "bootstrap")`, and
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) as an
R-engine baseline. If `pbkrtest` is installed, it also times
`pbkrtest::PBmodcomp()`.

``` r

data.frame(
  route = c(
    "mixeff_test_effect_bootstrap",
    "mixeff_confint_bootstrap",
    "mixeff_compare_bootstrap_lrt",
    "lme4_bootMer_fixef_distribution",
    "pbkrtest_PBmodcomp"
  ),
  target = c(
    "term-level fixed-effect p-value",
    "fixed-effect confidence interval",
    "model comparison p-value",
    "fixed-effect bootstrap distribution",
    "parametric bootstrap model comparison"
  )
)
#>                             route                                target
#> 1    mixeff_test_effect_bootstrap       term-level fixed-effect p-value
#> 2        mixeff_confint_bootstrap      fixed-effect confidence interval
#> 3    mixeff_compare_bootstrap_lrt              model comparison p-value
#> 4 lme4_bootMer_fixef_distribution   fixed-effect bootstrap distribution
#> 5              pbkrtest_PBmodcomp parametric bootstrap model comparison
```

## What does inference method choice cost?

Speed is only one part of the inference decision. The route table also
needs calibration evidence: how often a method rejects under a null, how
often it detects a small alternative, and whether interval routes cover
that small alternative.

The lightweight simulation scaffold records those three quantities over
four fixtures: an interior fit, a boundary-like fit, a reduced-rank fit,
and a small-group fit. The checked-in summary is produced by the default
fast mode, which excludes bootstrap routes so the script runs in
seconds. Slow mode adds bootstrap and bootstrap-LRT routes for local
evidence gathering.

``` r

simulation_path <- system.file(
  "extdata", "inference-method-simulation-summary.csv",
  package = "mixeff"
)
if (!nzchar(simulation_path)) {
  simulation_file <- file.path("inst", "extdata",
                               "inference-method-simulation-summary.csv")
  simulation_path <- if (file.exists(simulation_file)) {
    simulation_file
  } else {
    file.path("..", simulation_file)
  }
}
simulation <- read.csv(simulation_path, stringsAsFactors = FALSE)

head(simulation)
#>              method  fixture mode n_reps type_I_error power_at_alt
#> 1 asymptotic_wald_z interior fast      2            0            1
#> 2     satterthwaite interior fast      2            0            1
#> 3     kenward_roger interior fast      2           NA           NA
#> 4         bootstrap interior fast      0           NA           NA
#> 5           profile interior fast      2           NA           NA
#> 6     bootstrap_lrt interior fast      0           NA           NA
#>   coverage_at_alt
#> 1               1
#> 2               1
#> 3              NA
#> 4              NA
#> 5              NA
#> 6              NA
```

``` r

system2(
  "Rscript",
  c(
    "inst/benchmarks/inference-method-simulation.R",
    "--mode=slow",
    "--reps=100",
    "--nsim=199"
  )
)
```

The plot below shows whichever metrics are available for each route.
P-value routes contribute Type I error and power; interval routes
contribute coverage. Routes not wired in the current package version
remain in the CSV with `n_reps = 0`, so downstream reports can
distinguish “not run in fast mode” from “not part of the contract”.

``` r

library(ggplot2)

simulation_long <- stats::reshape(
  simulation,
  varying = c("type_I_error", "power_at_alt", "coverage_at_alt"),
  v.names = "value",
  timevar = "metric",
  times = c("Type I error", "Power at small alternative", "Coverage at small alternative"),
  direction = "long"
)
simulation_long <- subset(simulation_long, n_reps > 0 & is.finite(value))
simulation_long$metric <- factor(
  simulation_long$metric,
  levels = c("Type I error", "Power at small alternative",
             "Coverage at small alternative")
)

ggplot(simulation_long, aes(method, value, fill = method)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_grid(metric ~ fixture) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "empirical proportion") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
```

![Inference simulation metrics by method and
fixture.](benchmarking_files/figure-html/inference-simulation-plot-1.png)

## How do you read the result?

For a local speed claim, read the summary CSV and compute ratios from
elapsed seconds. Values above 1 mean the baseline route took longer than
the `mixeff` route.

``` r

summary <- read.csv("benchmarks/bootstrap-inference/bootstrap-inference-summary.csv")
mixeff_ci <- subset(summary, route == "mixeff_confint_bootstrap")
lme4_boot <- subset(summary, route == "lme4_bootMer_fixef_distribution")

data.frame(
  comparison = "lme4 bootMer / mixeff bootstrap CI",
  speedup = lme4_boot$median_sec / mixeff_ci$median_sec
)
```

Treat that ratio as local evidence, not a universal constant. It depends
on the formula, optimizer, random-effects structure, `nsim`, and
machine. The point of the benchmark scripts is to make the claim
reproducible for the model class you are actually fitting.

## What should go in a report?

When you report benchmark results, include the model formula, number of
rows, number of grouping levels, random-effects structure, `nsim`,
repetitions, and the package versions. The raw CSV has the timing rows;
the summary CSV is for tables and figures.

``` r

data.frame(
  field = c("formula", "nobs", "groups", "random_terms", "nsim", "reps", "versions"),
  reason = c(
    "defines the fixed and random structure being refit",
    "controls the amount of data processed per fit",
    "controls the random-effect dimension",
    "controls covariance structure and optimizer work",
    "controls bootstrap cost",
    "shows timing stability",
    "makes the run reproducible"
  )
)
#>          field                                             reason
#> 1      formula defines the fixed and random structure being refit
#> 2         nobs      controls the amount of data processed per fit
#> 3       groups               controls the random-effect dimension
#> 4 random_terms   controls covariance structure and optimizer work
#> 5         nsim                            controls bootstrap cost
#> 6         reps                             shows timing stability
#> 7     versions                         makes the run reproducible
```
