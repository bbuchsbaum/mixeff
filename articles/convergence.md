# Convergence, Boundaries, and Singularity

``` r

library(mixeff)
```

[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
reports the convergence status of a model fitted by
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) or
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md). This
vignette explains interior and reduced-rank solutions, the test for
singularity, and the use of
[`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
to check a fitted optimum. It also describes what happens when the
engine substitutes a different GLMM estimator.

## Interior solutions

The first example is a small longitudinal study with twelve clinics, a
week effect, and a treatment contrast.

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
fit_status(fit)
#> [1] "converged_interior"
```

`converged_interior` means the engine accepted a finite interior
solution under its convergence rules, with every variance component away
from zero. Individual derivative checks can still be unavailable; the
certificate below records which checks completed. The status does not
establish global optimality, correct model specification, or that the
data carry enough information.

Three statuses distinguish the main outcomes:

- `converged_interior` — the engine accepted a finite interior solution
  under its convergence rules.
- `converged_reduced_rank` — the engine accepted a solution whose
  random-effect covariance has lower rank than requested.
- `not_optimized` — the run stopped without certifying a stationary
  point, for example in a flat or boundary region. Estimates are
  present, but convergence has not been established.
  [`summary()`](https://rdrr.io/r/base/summary.html) also records this
  status in a `Notes:` line.

Optimization is silent and runs within one native call under bounded
evaluation budgets; the
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
help pages and
[`?mm_control`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
document the budgets and how to cap them.

[`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
reports the optimizer, stopping condition, and checks used to assign the
convergence status:

``` r

optimizer_certificate(fit)
#> Optimizer certificate:
#>                    metric                value
#>                    status   converged_interior
#>                 optimizer       pattern_search
#>                 objective     95.0585066127678
#>                iterations                   23
#>        free_gradient_norm 0.000357840060177139
#>   projected_gradient_norm 0.000357840060177139
#>         hessian_eigen_min     20.7651381551018
#>              hessian_rank                    1
#>          information_rank                    1
#>         optimizer_hessian          approximate
#>  derivative_certification          approximate
#>   optimizer_hessian: active-subspace Hessian from finite differences of the
#>     analytic gradient is positive semidefinite
#>   derivative_certification: analytic-gradient KKT check and
#>     finite-difference-of-gradient Hessian check passed
```

For joint-Laplace GLMMs, optimizer evidence and inference covariance are
separate records. The fit-time optimizer certificate does not include a
full Hessian assessment over fixed effects (beta) and covariance
parameters (theta); it labels that check `not_assessed`, with the stage
in the reason. Optimization may still compute a Hessian for polishing.
The label does not mean that a Hessian computation failed.

The later inference stage attempts a central finite-difference Hessian
of minus twice the log likelihood over beta and interior theta
coordinates. It omits theta coordinates too close to their lower bounds
for central differences. When this Hessian passes its checks,
fixed-effect covariance is twice the beta block of its inverse. If that
route fails, an available RX fallback conditions on theta and carries a
different method and lower reliability. Inspect
`attr(vcov(gfit), "mm_method")`, `"mm_status"`, `"mm_reliability"`, and
`"mm_notes"` for the actual route and its scope. A successful inference
calculation does not rewrite the fit-time certificate.

Raw gradient norms depend on parameter coordinates and are not measures
of coefficient precision. Curvature-scaled stationarity evidence, when
available, is recorded separately in
`optimizer_certificate(gfit)$raw$stationarity_decrement`. Read both with
their stated tolerances; a small scaled quantity does not change the
qualification recorded in an approximate raw-gradient certificate.

## Boundary and reduced-rank fits

The next example produces a singular fit. We simulate ten daily
reaction-time measurements for each of eighteen subjects, with nearly
perfectly correlated subject intercepts and slopes.

[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) reports
`boundary (singular) fit` for this model:

``` r

m <- suppressMessages(lme4::lmer(
  rt ~ days + (1 + days | subj),
  data = sleep_like
))
m
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: rt ~ days + (1 + days | subj)
#>    Data: sleep_like
#> REML criterion at convergence: 1647.697
#> Random effects:
#>  Groups   Name        Std.Dev. Corr 
#>  subj     (Intercept) 23.682        
#>           days         3.221   1.00 
#>  Residual             20.065        
#> Number of obs: 180, groups:  subj, 18
#> Fixed Effects:
#> (Intercept)         days  
#>     239.879        9.251  
#> optimizer (nloptwrap) convergence code: 0 (OK) ; 0 optimizer warnings; 1 lme4 warnings
lme4::isSingular(m)
#> [1] TRUE
```

In `mixeff`,
[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
reports the convergence outcome and
[`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
tests for singularity. For more detail,
[`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
returns codes with severity and stage, while
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
compares the requested model with the effective specification and fitted
result.

``` r

sfit <- lmm(
  rt ~ days + (1 + days | subj),
  sleep_like,
  control = mm_control(verbose = -1)
)

fit_status(sfit)
#> [1] "converged_reduced_rank"
is_singular(sfit)
#> [1] TRUE
diagnostics(sfit)$table[, c("code", "severity", "stage", "message")]
#>                 code severity         stage
#> 1 optimizer_recovery     info certification
#> 2 covariance_reduced     info certification
#>                                                                            message
#> 1            optimizer recovered after covariance KKT-guided restart (2x2 block 0)
#> 2 fitted covariance for (1 + days | subj) has effective rank 1 of requested rank 2
changes(sfit)
#> Model changes:
#>   Fitted covariance for (1 + days | subj): requested rank 2, fitted rank 1 [reduced_rank].
#> Stage-by-stage records available via $table.
```

The requested random-effect covariance for `(1 + days | subj)` has two
dimensions: a subject intercept and a subject slope. The fitted
covariance has effective rank one. This is recorded as the
`covariance_reduced` diagnostic and the `converged_reduced_rank` status,
which also appears when the object is printed.

``` r

print(sfit)
#> Linear mixed model fit by REML
#> Formula: rt ~ days + (1 + days | subj)
#> Fit status: converged_reduced_rank
#> Optimizer: trust_bq; iterations: 30; objective: 1647.7
#> nobs: 180, sigma: 20.065, logLik: -823.849
#> Fixed effects:
#> (Intercept)        days 
#>   239.87900     9.25095 
#> 
#> Fitted covariance state:
#> The fitted covariance matrix is rank-deficient.
#>   r0: requested rank 2; fitted effective rank 1.
#> Use changes(fit) to see which dimension was unsupported.
#> Use random_options(spec, group = subj) to inspect lower-dimensional covariance choices.
#> Audit verbs: audit(), diagnostics(), inference_table(), model_report()
```

A reduced-rank covariance does not automatically invalidate the model,
but it restricts the inference methods available in `mixeff`. These
restrictions are described at the end of this vignette.

## How singularity relates to θ and rank

To read the boundary structure itself, look at the model’s internal
coordinates. The random-effect covariance Σ is not optimized directly.
The engine searches over a vector θ that fills a lower-triangular factor
Λ, with Σ = σ²ΛΛᵀ; θ is by construction exactly the set of free entries
of Λ. This is the parameterization Bates et al. (2015) describe for
`lme4`; `MixedModels.jl` and the engine behind `mixeff` inherit it.
[`parameterization()`](https://bbuchsbaum.github.io/mixeff/reference/parameterization.md)
shows the mapping for the singular fit:

``` r

parameterization(sfit)$table[, c("theta_name", "theta_value", "lambda_value")]
#>                     theta_name  theta_value lambda_value
#> 1 theta[0:intercept,intercept] 1.180242e+00 1.180242e+00
#> 2      theta[0:days,intercept] 1.605414e-01 1.605414e-01
#> 3           theta[0:days,days] 1.678443e-08 1.678443e-08
```

The *diagonal* entry for the `days` slope is pinned at numerical zero —
that is the boundary event behind the `converged_reduced_rank` status.

A zero entry in θ does not always imply singularity. A zero
**off-diagonal** entry of Λ can occur in a full-rank covariance. In the
two-dimensional example below, it sets the covariance to zero. In larger
blocks, covariance entries are sums of products, so one zero factor
entry need not make the corresponding covariance zero.

Singularity depends on the rank of the whole covariance matrix. For this
square triangular factor, it occurs exactly when at least one
**diagonal** entry of Λ is zero. The following example contrasts the two
cases:

``` r

l11 <- 1.2
l21 <- 0.4
l22 <- 0.9

# Zero the OFF-DIAGONAL entry (l21 = 0): zero covariance, full rank.
L_offdiag_zero <- matrix(c(l11, 0, 0, l22), nrow = 2)
Sigma_offdiag_zero <- L_offdiag_zero %*% t(L_offdiag_zero)
Sigma_offdiag_zero
#>      [,1] [,2]
#> [1,] 1.44 0.00
#> [2,] 0.00 0.81
eigen(Sigma_offdiag_zero, symmetric = TRUE)$values
#> [1] 1.44 0.81

# Zero a DIAGONAL entry (l22 = 0): the rank collapses to 1.
L_diag_zero <- matrix(c(l11, l21, 0, 0), nrow = 2)
Sigma_diag_zero <- L_diag_zero %*% t(L_diag_zero)
Sigma_diag_zero
#>      [,1] [,2]
#> [1,] 1.44 0.48
#> [2,] 0.48 0.16
eigen(Sigma_diag_zero, symmetric = TRUE)$values
#> [1] 1.600000e+00 4.241052e-17
```

With the off-diagonal at zero, both eigenvalues are positive: the
intercept and slope are uncorrelated, and each still has its own
variance — a full-rank fit with the covariance structure specified by
`(1 + days || subj)`. With a diagonal entry at zero, one eigenvalue is
zero and the two random effects collapse onto a single dimension: that
is singularity.

[`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
follows lme4’s tolerance rule: it checks whether a zero-bounded diagonal
θ is below `tol` (default `1e-4`). For older fits without a theta map,
it falls back to the recorded engine covariance and fit statuses. A
tolerance-based result can differ from an exact rank verdict.

## Verifying an optimum: `verify_convergence()`

A convergence status describes the run that produced the fit. To check
the stability of that solution,
[`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
re-runs the fit under the engine’s bounded verification workflow: a
restart from the fitted optimum, one or more jittered restarts, and
(LMM-only, opt-in) an alternate-optimizer consensus pass. The engine
judges agreement against explicit tolerances on the objective, θ, and
the fixed effects, and returns the verdict for R to format.

The R bridge reuses a live fitted handle or restores the saved fitted
state before invoking that workflow. Verification compares its
deliberate searches against the recorded solution, including after
reload. If neither a valid handle nor a compatible snapshot is
available, verification refuses rather than optimizing a replacement
baseline.

``` r

vc <- verify_convergence(fit)
vc
#> Convergence verification (status: restart_agrees)
#>   restart from fitted theta agrees with the recorded optimum
#> Runs:
#>                 label      optimizer return_code objective_value
#>  restart_from_optimum pattern_search     SUCCESS         95.0585
#>      jitter_restart_1 pattern_search     SUCCESS         95.0585
#>  objective_delta theta_delta beta_delta agrees
#>      0.00000e+00 0.00000e+00 0.0000e+00   TRUE
#>      1.18234e-10 3.20312e-07 3.9968e-15   TRUE
#> Tolerances: objective 1e-05; theta 0.001; beta 1e-04
```

The verdict summarizes agreement with the reference optimum within the
printed tolerances (`restart_agrees`, `optimizer_consensus`, `fragile`,
`unstable`, or `not_run`). A verdict of `restart_agrees` is evidence of
a *stable* optimum — the bounded re-runs came back to the same place —
not proof of a global one. The `consensus` pass defaults to `FALSE`
because this vendored build compiles without the optional `nlopt`
backend; for some models the engine’s alternate optimizer choice is an
nlopt optimizer, and its absence would read as a spurious `fragile`
verdict that reflects the build, not the fit.

Verification and rank are different questions. Running the verifier on
the singular fit from above also comes back `restart_agrees`: a
reduced-rank optimum can be stable across restarts.

``` r

vs <- verify_convergence(sfit)
vs$status
#> [1] "restart_agrees"
```

The GLMM method uses the same workflow with wider default tolerances —
objective `1e-4` and beta `1e-3`, against the LMM defaults of `1e-5` and
`1e-4` — because GLMM objective evaluations include inner-PIRLS
numerical error and the joint path’s derivative-free beta search bounds
beta reproducibility at its stopping tolerance. The consensus pass is
refused for GLMM fits (the GLMM fitting routines select their own
optimizers).

``` r

gfit <- glmm(
  y ~ x + (1 | g),
  trial,
  family = binomial(),
  control = mm_control(verbose = -1)
)

gv <- verify_convergence(gfit)
gv
#> Convergence verification (status: restart_agrees)
#>   restart from fitted theta agrees with the recorded optimum
#> Runs:
#>                 label optimizer                return_code objective_value
#>  restart_from_optimum  trust_bq JOINT_LAPLACE:FTOL_REACHED         179.457
#>      jitter_restart_1  trust_bq JOINT_LAPLACE:FTOL_REACHED         179.457
#>  objective_delta theta_delta  beta_delta agrees
#>      2.84217e-14 3.70270e-12 1.83527e-11   TRUE
#>      5.68434e-14 4.26725e-12 2.10366e-11   TRUE
#> Tolerances: objective 1e-04; theta 0.001; beta 0.001
```

## GLMM estimator substitution

`glmm(method = "joint_laplace")` requests the joint-Laplace estimator,
which supports Wald inference when the fitted covariance passes its
checks. On some designs, the engine cannot certify a joint-Laplace
optimum and returns a fit from the profiled-PIRLS estimator instead. The
substitution is recorded in several places:

- `optimizer_certificate(fit)$raw$estimator_substitution` contains the
  requested and effective methods, the requested method’s fit status and
  return code, and the reason for substitution. Use this record in
  scripts.
- A condition of class `mm_estimator_substitution_notice` reports the
  change when the fit returns. Suppress it with
  `mm_control(verbose = -1)`.
- [`summary()`](https://rdrr.io/r/base/summary.html) prints an
  `estimator substitution:` note naming both estimators, regardless of
  verbosity settings or the fallback fit’s convergence status.
- [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
  restores and checks the effective estimator that produced the returned
  estimates.

When no substitution occurred, the record is absent:

``` r

optimizer_certificate(gfit)$raw$estimator_substitution
#> NULL
```

For example, in a reproduction of a published `glmer()` analysis of OSF
study data, an error-rate model included an open-practice × year
interaction. With raw, uncentered `Year`, the joint-Laplace fit did not
certify and the engine returned a profiled-PIRLS fallback with Wald
inference withheld. Centering the predictor as `Year - 2015` allowed the
joint fit to certify and agree with `glmer()` on the same data. A
scaling advisory identifies the poorly scaled columns; centering or
rescaling them can improve optimization.

## Inference at a boundary

On singular fits, `mixeff` withholds Satterthwaite and Kenward–Roger
degrees of freedom under its validated-support policy. The results
include reason codes, and parametric bootstrap tests remain available.
See
[`vignette("inference")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md)
for the methods, restrictions, and examples.

## References

- Bates, D., Mächler, M., Bolker, B., and Walker, S. (2015). Fitting
  linear mixed-effects models using lme4. *Journal of Statistical
  Software*, 67(1), 1–48. <https://doi.org/10.18637/jss.v067.i01>
