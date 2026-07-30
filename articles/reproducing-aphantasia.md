# Reproducing the Loo Aphantasia GLMMs

``` r

library(mixeff)
```

The fixture bundled with `mixeff` is an anonymized copy of the
trial-level data used by the revision 3 Loo aphantasia manuscript
analysis. It is intended as a real GLMM reproduction target: large
enough to expose optimizer and naming drift, but small enough to ship as
package test data.

``` r

fixture_candidates <- c(
  system.file("extdata", "aphantasia", package = "mixeff"),
  file.path("..", "inst", "extdata", "aphantasia"),
  file.path("inst", "extdata", "aphantasia")
)
fixture_dir <- fixture_candidates[
  dir.exists(fixture_candidates) & nzchar(fixture_candidates)
][1L]
trials <- readRDS(file.path(fixture_dir, "trials.rds"))
metadata <- readRDS(file.path(fixture_dir, "metadata.rds"))
reference <- jsonlite::fromJSON(file.path(fixture_dir, "reference.json"))

c(
  trials = nrow(trials),
  participants = length(unique(trials$participant)),
  metadata_rows = nrow(metadata)
)
#>        trials  participants metadata_rows 
#>         25916            76            76
```

The primary analysis uses the occluded trials, excludes the four
intermediate VVIQ controls, and models accuracy with crossed participant
and item effects.

``` r

excluded <- unlist(reference$excluded_participants, use.names = FALSE)

primary <- subset(
  trials,
  bubbled == "yes" & !is.na(correct) & !participant %in% excluded
)

prepare_model_data <- function(dat, stimtype = FALSE) {
  out <- transform(
    dat,
    participant = factor(participant),
    item = factor(trial_image),
    group = factor(ifelse(aphantasia == "yes", "aphant", "control"),
                   levels = c("control", "aphant")),
    mask = factor(ifelse(back_masked == "yes", "masked", "unmasked"),
                  levels = c("unmasked", "masked")),
    block = factor(block_num),
    soa_log = log(SOA)
  )
  out$soa_s <- as.numeric(scale(out$soa_log))
  if (stimtype) {
    out$stimtype <- factor(
      ifelse(out$bubbled == "yes", "occluded", "intact"),
      levels = c("intact", "occluded")
    )
  }
  out
}

primary_dat <- prepare_model_data(primary)
table(primary_dat$group, primary_dat$mask)
#>          
#>           unmasked masked
#>   control     4920   4920
#>   aphant      3720   3720
```

The live fits are opt-in because the primary GLMM alone takes several
minutes on a laptop. Set `MIXEFF_RUN_APHANTASIA_VIGNETTE=true` to
execute the model chunks while rendering. The test suite uses the same
principle: the core reproduction runs under
`MIXEFF_RUN_APHANTASIA=true`, and the slower S1 random-effects stability
variants run under `MIXEFF_RUN_APHANTASIA_STRESS=true`.

``` r

primary_fit <- glmm(
  correct ~ group * mask * soa_s + block +
    (1 + mask + soa_s || participant) + (1 | item),
  primary_dat,
  family = binomial(),
  control = mm_control(verbose = -1)
)

c(
  logLik = as.numeric(logLik(primary_fit)),
  AIC = AIC(primary_fit)
)
fixef(primary_fit)
```

Without live refitting, the frozen lme4 reference records the target
values the integration test compares against.

``` r

unlist(reference$models$primary[c("nobs", "logLik", "AIC")])
#>      nobs    logLik       AIC 
#> 17280.000 -9966.062 19962.124
unlist(reference$models$primary$fixef)
#>                  (Intercept)                  groupaphant 
#>                   0.45590105                   0.19961426 
#>                   maskmasked                        soa_s 
#>                  -0.30038083                   0.39045247 
#>                       block2       groupaphant:maskmasked 
#>                   0.10639773                  -0.16663973 
#>            groupaphant:soa_s             maskmasked:soa_s 
#>                  -0.05751370                   0.07416383 
#> groupaphant:maskmasked:soa_s 
#>                   0.10641824
```

The same fixture supports the manuscript sensitivity and specificity
fits: the sensitivity model assigns the four intermediate VVIQ
participants to the control group. The estimator matters for two of the
models that follow. The intact high-baseline Bernoulli model reaches
near-exact lme4 fixed-effect and log-likelihood parity only under
`method = "joint_laplace"`, which is how the opt-in reproduction test
fits it (about 40 seconds per fit on a release build). On the default
profiled path its coefficients drift well past the strict tolerances, so
the chunks below do not claim parity for it. The combined model is
fitted by the default profiled path here and in the test suite; like
every profiled GLMM fit, it carries the engine’s support note that a
fast-PIRLS fit is not certified as glmer joint-Laplace parity.

``` r

sensitivity <- subset(trials, bubbled == "yes" & !is.na(correct))
sensitivity$aphantasia[sensitivity$participant %in% excluded] <- "no"
sensitivity_dat <- prepare_model_data(sensitivity)

intact <- subset(
  trials,
  bubbled == "no" & !is.na(correct) & !participant %in% excluded
)
intact_dat <- prepare_model_data(intact)

combined <- subset(trials, !is.na(correct) & !participant %in% excluded)
combined_dat <- prepare_model_data(combined, stimtype = TRUE)

sensitivity_fit <- glmm(
  correct ~ group * mask * soa_s + block +
    (1 + mask + soa_s || participant) + (1 | item),
  sensitivity_dat,
  family = binomial(),
  control = mm_control(verbose = -1)
)

intact_fit <- glmm(
  correct ~ group * mask * soa_s + block +
    (1 + mask + soa_s || participant) + (1 | item),
  intact_dat,
  family = binomial(),
  control = mm_control(verbose = -1)
)

combined_fit <- glmm(
  correct ~ group * mask * soa_s * stimtype + block +
    (1 + mask + soa_s || participant) + (1 | item),
  combined_dat,
  family = binomial(),
  control = mm_control(verbose = -1)
)

rbind(
  sensitivity = c(logLik = as.numeric(logLik(sensitivity_fit)),
                  AIC = AIC(sensitivity_fit)),
  intact = c(logLik = as.numeric(logLik(intact_fit)), AIC = AIC(intact_fit)),
  combined = c(logLik = as.numeric(logLik(combined_fit)),
               AIC = AIC(combined_fit))
)
```

The RT sensitivity is a Gaussian LMM over correct trials with finite
positive reaction times.

``` r

rt_dat <- subset(primary_dat, correct == 1 & is.finite(rt) & rt > 0)
rt_dat$log_rt <- log(rt_dat$rt)

rt_fit <- lmm(
  log_rt ~ group * mask * soa_s + block +
    (1 | participant) + (1 | item),
  rt_dat,
  REML = FALSE,
  control = mm_control(verbose = -1)
)

c(logLik = as.numeric(logLik(rt_fit)), AIC = AIC(rt_fit))
fixef(rt_fit)
```

## What the manuscript reports, computed three ways

Since the engine began serializing the GLMM fixed-effect covariance
artifact, the quantities the manuscript actually reports are available
through three `mixeff` functions. The chunks below run only when live
fitting is enabled (`MIXEFF_RUN_APHANTASIA_VIGNETTE=true`); otherwise
the same calls remain valid against any locally-built fit.

`summary(fit, tests = "coefficients")` returns a Wald-z fixed-effect
table built from the PIRLS/Laplace working-Hessian covariance:

``` r

sm <- summary(primary_fit, tests = "coefficients")
round(sm$coefficients, 3)
sm$vcov_status
```

The status block flags `reliability = "moderate"`: the working-Hessian
covariance is not the same estimator as `lme4::vcov(glmer_fit)`. On this
dataset its Wald standard errors for the difference-in-differences
contrasts run about 11% *smaller* than glmer’s — anti-conservative, so
p-values lean optimistic — which the package’s parity ledger classifies
as an expected mismatch with a 15% bound. The manuscript’s qualitative
conclusions do not flip, but a reader re-using these SEs should know
which way they err.

The manuscript’s primary estimand — the difference-in-differences
contrast at the centered SOA and at the focal 25 ms SOA — is a linear
combination of fixed effects, which
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
computes:

``` r

soa_s_25 <- (log(0.025) - mean(primary_dat$soa_log)) /
              sd(primary_dat$soa_log)

dd_center <- mm_lincomb(
  primary_fit,
  c("groupaphant:maskmasked" = 1)
)
dd_25 <- mm_lincomb(
  primary_fit,
  c("groupaphant:maskmasked"       = 1,
    "groupaphant:maskmasked:soa_s" = soa_s_25)
)
list(centered_soa = dd_center, ms25 = dd_25)
```

Compare directly against the `lme4` reference frozen in the fixture:

``` r

as.data.frame(reference$inference$primary_dd)
#>          where   estimate         SE         z          p      lower
#> 1 centered_soa -0.1666397 0.08079478 -2.062506 0.03915958 -0.3249975
#> 2        25_ms -0.3391448 0.14001919 -2.422131 0.01542979 -0.6135824
#>          upper
#> 1 -0.008281958
#> 2 -0.064707203
```

Both contrasts come out negative (larger masking cost in aphantasia), as
in the manuscript. Be precise about what is pinned. The opt-in test
suite enforces agreement with the frozen reference at ±0.02 on the
estimates and ±0.01 on the SEs; at the worst corner of those bands the
centered-SOA p-value — .039 in the reference — could reach .106, so the
significance verdict there is not guaranteed by the tests. In practice
mixeff lands on the other side: its SEs run smaller than glmer’s (see
above), so its centered-SOA p-value comes out below the reference value,
not above it. The sign and approximate size of both contrasts are
pinned.

`emmeans` works on `mm_glmm` via `emm_basis.mm_glmm`. Population-level
cell means at the centered SOA, on the response (probability) scale,
recover the same group × mask × group pattern the manuscript reports:

``` r

em <- emmeans::emmeans(
  primary_fit, ~ mask | group,
  at   = list(soa_s = 0),
  type = "response"
)
as.data.frame(summary(em))
```

The reverse-pairwise contrast inside each group is the group-conditional
masking cost on the log-odds scale:

``` r

emmeans::contrast(
  emmeans::emmeans(primary_fit, ~ mask | group,
                   at = list(soa_s = 0)),
  method = "revpairwise"
)
```

## Out of scope here

The full sensitivity, intact-stimulus, combined-stimtype, log-RT, S1
random-effects-spec stability, S7 age-covariate, and S9 folder-based
age-matched analyses are reproducible against the same fixture and
reference via `tests/testthat/test-aphantasia-reproduction.R` when
`MIXEFF_RUN_APHANTASIA=true` is set. The S3 leave-one-participant-out
sweep, S4 specification curve, and S5 rstanarm posterior are
intentionally left out of the regular reproduction run: the first two
are heavy opt-in jobs, and `mixeff` is not a Bayesian engine.

## Caveats

- GLMM Wald inference in `mixeff` uses the PIRLS/Laplace working-Hessian
  covariance, graded `reliability = "moderate"`. Its SEs run about 11%
  smaller than `lme4::vcov()` on this dataset (anti-conservative; ledger
  bound 15%). A controlled decomposition attributes roughly 73% of the
  gap to the native `||` random-effect family, 26% to a uniform
  working-Hessian scale factor (tracked upstream), and 1% to optimizer
  drift.
- Population-level GLMM prediction (`re.form = NA`, `type = "link"` or
  `"response"`) is supported and matches `predict(glmer, re.form = NA)`
  on joint-Laplace fits; `emmeans(..., type = "response")` remains the
  route for averaged marginal probabilities.
- Coefficient names are lme4-identical since 0.2.0 (`"groupaphant"`), so
  hand-written lincombs copy straight across from `lme4` code.

## Citation

If you use the bundled fixture or reproduce these analyses in downstream
work, please cite:

> Loo, C., & Buchsbaum, B. R. (2026). Fragile recurrent processing in
> aphantasia: Evidence from visual pattern completion. *Manuscript
> submitted for publication.*
