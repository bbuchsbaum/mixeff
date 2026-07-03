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
participants to the control group. The intact high-baseline Bernoulli
model defaults to the full-budget joint-Laplace route
(`method = "joint_laplace"`) in the opt-in reproduction gate, reaching
near-exact lme4 fixed-effect and log-likelihood parity on a release
build (~40 s per fit; only the AIC parameter-count semantics for the
double-bar factor expansion remain ledgered). The combined model stays
on the profiled ledger path: the engine rejects its joint candidate for
that case and falls back to fast-PIRLS with an explicit
`documented_divergence` diagnostic.

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

## Inferential surfaces

Since `mixeff-rs` started serializing the GLMM fixed-effect covariance
artifact (pin 5e72e0b), the inferential surfaces the manuscript actually
reports are available through three `mixeff` primitives. The chunks
below run only when live fitting is enabled
(`MIXEFF_RUN_APHANTASIA_VIGNETTE=true`); otherwise the same calls remain
valid against any locally-built fit.

`summary(fit, tests = "coefficients")` returns a Wald-z fixed-effect
table built from the PIRLS/Laplace working-Hessian covariance:

``` r

sm <- summary(primary_fit, tests = "coefficients")
round(sm$coefficients, 3)
sm$vcov_status
```

The status block flags `reliability = "moderate"`: the working-Hessian
flavor is close to but not bit-identical with `lme4::vcov(glmer_fit)`.
SE estimates on this dataset drift by ~5-10% in absolute terms, without
flipping any of the manuscript’s qualitative conclusions.

The manuscript’s primary estimand — the difference-in-differences
contrast at the centered SOA and at the focal 25 ms SOA — is a linear
combination of fixed effects, and
[`mm_lincomb()`](https://bbuchsbaum.github.io/mixeff/reference/mm_lincomb.md)
is its front door:

``` r

soa_s_25 <- (log(0.025) - mean(primary_dat$soa_log)) /
              sd(primary_dat$soa_log)

dd_center <- mm_lincomb(
  primary_fit,
  c("group: aphant:mask: masked" = 1)
)
dd_25 <- mm_lincomb(
  primary_fit,
  c("group: aphant:mask: masked"       = 1,
    "group: aphant:mask: masked:soa_s" = soa_s_25)
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

Both contrasts reproduce the manuscript’s sign and significance class:
negative (larger masking cost in aphantasia), CIs excluding zero, and
*p* below the conventional α = .05 at both 25 ms and the centered SOA.

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
intentionally not part of the regular reproduction surface: the first
two are heavy opt-in jobs, and `mixeff` is not a Bayesian engine.

## Caveats

- GLMM Wald inference in `mixeff` is the PIRLS/Laplace working-Hessian
  flavor, advertised as `mm_reliability = "moderate"`. Absolute SEs
  drift 5–10% versus `lme4::vcov()` on this dataset.
- `predict(mm_glmm_fit, ...)` is not yet certified for population-level
  prediction; that lives behind its own bead. Use
  `emmeans(..., type = "response")` for marginal probabilities.
- Coefficient names follow `mixeff`’s `"group: aphant"` scheme, not
  `lme4`’s compact `"groupaphant"`. For copy-pasting hand-written
  lincombs from `lme4` user code, `gsub(": ", "", x, fixed = TRUE)`
  recovers the compact form.

## Citation

If you use the bundled fixture or reproduce these analyses in downstream
work, please cite:

> Loo, C., & Buchsbaum, B. R. (2026). Fragile recurrent processing in
> aphantasia: Evidence from visual pattern completion. *Manuscript
> submitted for publication.*
