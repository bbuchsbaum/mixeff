# Convergence, Boundaries, and Singularity

``` r

library(mixeff)
```

Every fit returned by
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) or
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
carries a labeled convergence status, and the label is a claim with a
precise scope. This vignette is the reference for that vocabulary: what
`converged_interior` does and does not establish, how boundary and
reduced-rank outcomes are recorded, how singularity is actually decided,
how to check an optimum with
[`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md),
and what happens when a requested GLMM estimator does not certify. The
other vignettes cross-reference these semantics rather than restating
them.

## What `converged_interior` establishes — and what it does not

Start with a fit that goes well: a small longitudinal study with twelve
clinics, a week effect, and a treatment contrast.

``` r

fit <- lmm(
  score ~ week + treatment + (1 | clinic),
  clinic_visits,
  control = mm_control(verbose = -1)
)
fit_status(fit)
#> [1] "converged_interior"
```

`converged_interior` means the optimizer returned a finite interior
solution that passed mixeff’s declared numerical checks, with every
variance component away from zero. It does not establish global
optimality, correct model specification, or that the data carry enough
information — no numerical status can — but it is the status every other
block of a [`summary()`](https://rdrr.io/r/base/summary.html) assumes.

The status changes when something needs your attention. The vocabulary
has three values you will meet in practice:

- `converged_interior` — a finite interior optimum, as above.
- `converged_reduced_rank` — a genuine optimum at which the
  random-effect covariance has lower rank than the model requested; the
  next section walks through one.
- `not_optimized` — the conservative label for a run that stopped in a
  flat or boundary region *without* certifying a stationary point.
  Estimates are present but the engine does not call them an optimum,
  and [`summary()`](https://rdrr.io/r/base/summary.html) repeats that
  state in a `Notes:` line rather than relying on the header alone.

Optimization is silent and runs within one native call under bounded
evaluation budgets; the
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md)/[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md)
help pages and
[`?mm_control`](https://bbuchsbaum.github.io/mixeff/reference/mm_control.md)
document the budgets and how to cap them.

What the fit ran is recorded on the object.
[`optimizer_certificate()`](https://bbuchsbaum.github.io/mixeff/reference/optimizer_certificate.md)
reports the optimizer, the stopping condition, and the certificate-time
facts the status label summarizes:

``` r

optimizer_certificate(fit)
#> Optimizer certificate:
#>                   metric              value
#>                   status converged_interior
#>                optimizer     pattern_search
#>                objective   95.0585066127678
#>               iterations                 23
#>       free_gradient_norm 0.0003577914994027
#>  projected_gradient_norm 0.0003577914994027
#>        hessian_eigen_min   20.7641716529313
#>             hessian_rank                  1
#>         information_rank                  1
```

## Boundary and reduced-rank fits are labeled facts

Now a fit that does not go well — deliberately. Eighteen subjects, ten
daily reaction-time measurements each, with subject intercepts and
slopes that are nearly perfectly correlated by construction. This design
produces a singular fit in `lme4` and in `mixeff` alike.

[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) fits this model
and reports the situation in one short line: `boundary (singular) fit`.
The fact is correct. What it leaves implicit is *which* variance
component reached the boundary, what the effective rank of the
random-effect covariance is, and which downstream inference methods the
boundary takes off the table.

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

`mixeff` reports the same fact, and then unpacks it.
[`fit_status()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
names the convergence outcome,
[`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
answers the yes/no question,
[`diagnostics()`](https://bbuchsbaum.github.io/mixeff/reference/diagnostics.md)
returns stable codes with severity and stage, and
[`changes()`](https://bbuchsbaum.github.io/mixeff/reference/changes.md)
shows the requested-to-effective-to-fitted arc.

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
#> 1 boundary_parameter     info certification
#> 2 covariance_reduced     info certification
#>                                                                            message
#> 1           standard deviation for days in (1 + days | subj) is on its lower bound
#> 2 fitted covariance for (1 + days | subj) has effective rank 1 of requested rank 2
changes(sfit)
#> Model changes:
#>   Fitted covariance for (1 + days | subj): requested rank 2, fitted rank 1 [reduced_rank].
#> Stage-by-stage records available via $table.
```

The requested random-effect covariance for `(1 + days | subj)` has rank
2: one dimension for the subject intercept and one for the subject
slope. The fitted covariance has effective rank 1 — the model asked for
two random-effect dimensions, but the data support only one effective
dimension. That is recorded as the `covariance_reduced` diagnostic and
the `converged_reduced_rank` status, and the printed object carries the
tag in its header every time it prints — where `lme4` warns once at fit
time and then falls silent.

``` r

print(sfit)
#> Linear mixed model fit by REML
#> Formula: rt ~ days + (1 + days | subj)
#> Fit status: converged_reduced_rank
#> Optimizer: trust_bq; iterations: 314; objective: 1647.7
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

A reduced-rank covariance does not automatically invalidate the model.
It is a labeled fact about the fit, recorded in the object rather than
printed once and lost — and it changes which inference routes the
package will certify (one sentence on that at the end of this vignette).

## How singularity is decided: rank, not θ entries

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
#>                     theta_name theta_value lambda_value
#> 1 theta[0:intercept,intercept]   1.1802493    1.1802493
#> 2      theta[0:days,intercept]   0.1605448    0.1605448
#> 3           theta[0:days,days]   0.0000000    0.0000000
```

The *diagonal* entry for the `days` slope is pinned at numerical zero —
that is the boundary event behind the `converged_reduced_rank` status.

It is tempting to generalize this to “a zero θ entry means the
covariance is reduced-rank”. That reading is wrong. A zero
**off-diagonal** entry of Λ only says the corresponding covariance is
zero; the covariance matrix can keep full rank. Singularity is a
property of the covariance matrix’s rank — equivalently its
eigenstructure — and in this parameterization the rank collapses exactly
when a **diagonal** entry of Λ is zero. A two-by-two example makes the
difference concrete:

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
#> [1] 1.600000e+00 5.551115e-17
```

With the off-diagonal at zero, both eigenvalues are positive: the
intercept and slope are uncorrelated, and each still has its own
variance — an ordinary full-rank fit (the model `(1 + days || subj)`
*asks* for exactly this shape). With a diagonal entry at zero, one
eigenvalue is zero and the two random effects collapse onto a single
dimension: that is singularity.

Consistent with this,
[`is_singular()`](https://bbuchsbaum.github.io/mixeff/reference/is_singular.md)
does not scan θ for zeros. It reads the covariance status the engine
recorded on the fit — the fit status plus the per-block
effective-covariance status — so its answer is the engine’s rank
verdict, not an R-side reconstruction.

## Verifying an optimum: `verify_convergence()`

A convergence status describes the run that produced the fit. If you
want independent evidence,
[`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
re-runs the fit under the engine’s bounded verification workflow: a
restart from the fitted optimum, one or more jittered restarts, and
(LMM-only, opt-in) an alternate-optimizer consensus pass. The engine
judges agreement against explicit tolerances on the objective, θ, and
the fixed effects, and owns the verdict; R only formats it.

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

Each run either agrees with the reference optimum within the printed
tolerances or it does not; the verdict summarizes the runs
(`restart_agrees`, `optimizer_consensus`, `fragile`, `unstable`, or
`not_run`). A verdict of `restart_agrees` is evidence of a *stable*
optimum — the bounded re-runs came back to the same place — not proof of
a global one. The `consensus` pass defaults to `FALSE` because this
vendored build compiles without the optional `nlopt` backend; for some
models the engine’s alternate optimizer choice is an nlopt optimizer,
and its absence would read as a spurious `fragile` verdict that reflects
the build, not the fit.

Verification and rank are different questions. Running the verifier on
the singular fit from above also comes back `restart_agrees`: a
reduced-rank optimum can be a perfectly stable optimum.

``` r

vs <- verify_convergence(sfit)
vs$status
#> [1] "restart_agrees"
```

The GLMM method uses the same workflow with wider default tolerances —
objective `1e-4` and beta `1e-3`, against the LMM defaults of `1e-5` and
`1e-4` — because every GLMM objective evaluation carries inner-PIRLS
noise and the joint path’s derivative-free beta search bounds beta
reproducibility at its stopping tolerance. The consensus pass is refused
for GLMM fits (the GLMM fit drivers own optimizer selection).

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
#>                 label optimizer  return_code objective_value objective_delta
#>  restart_from_optimum    cobyla FTOL_REACHED         179.484     0.00000e+00
#>      jitter_restart_1    cobyla FTOL_REACHED         179.484     8.06031e-09
#>  theta_delta  beta_delta agrees
#>  0.00000e+00 0.00000e+00   TRUE
#>  7.40564e-06 5.80241e-07   TRUE
#> Tolerances: objective 1e-04; theta 0.001; beta 0.001
```

## When the requested estimator does not certify (GLMM)

`glmm(method = "joint_laplace")` requests the joint-Laplace estimator
(the route that certifies Wald inference;
[`vignette("glmm")`](https://bbuchsbaum.github.io/mixeff/articles/glmm.md)
is the operational quickstart). On some designs the joint route does not
certify an optimum. The engine then returns a *labelled* fallback fit
from the fast profiled-PIRLS estimator instead of failing or — worse —
silently presenting the fallback as the thing you asked for. The
substitution crosses the boundary as a typed record, in keeping with the
package’s no-silent-surgery contract:

- **The typed record.**
  `optimizer_certificate(fit)$raw$estimator_substitution` carries the
  requested method, the effective (fallback) method, the requested
  route’s fit status and return code, and the engine’s reason. It is the
  machine-readable contract surface; scripts should gate on it, not on
  message text.
- **The fit-time notice.** A condition of class
  `mm_estimator_substitution_notice` announces the substitution when the
  fit returns (suppressible with `mm_control(verbose = -1)`).
- **The ungated summary note.**
  [`summary()`](https://rdrr.io/r/base/summary.html) prints an
  `estimator substitution:` note naming both estimators regardless of
  verbosity settings, and independent of the fallback’s own
  `converged_*` status — a clean-looking status line must not hide the
  fact that a different estimator than requested produced the numbers.
- **The checking verb.**
  [`verify_convergence()`](https://bbuchsbaum.github.io/mixeff/reference/verify_convergence.md)
  on a substituted fit re-requests the joint route; when the
  verification refit falls back the same way — the expected case — every
  run verifies the profiled objective the fitted numbers actually came
  from.

On an ordinary fit the record is absent, and its absence is meaningful
(no substitution happened):

``` r

optimizer_certificate(gfit)$raw$estimator_substitution
#> NULL
```

A real example, from reproducing a published `glmer()` analysis of OSF
study data (an error-rate model with an open-practice × year
interaction): requesting `method = "joint_laplace"` with the raw,
uncentered `Year` predictor hands the engine an ill-scaled design, the
joint route does not certify, and the result comes back as a labelled
fast-PIRLS fallback — fit-time notice, summary note, and typed record
all present, with Wald inference withheld. Centering the year
(`Year - 2015`) makes the joint route certify and match `glmer()` on the
same data. mixeff flags this situation at fit time: a scaling advisory
names the ill-scaled columns, and centering or rescaling the predictor —
as the year centering does here — is the standard remedy.

## What a boundary refuses

A boundary fit also changes which inference routes the package will
certify: Satterthwaite and Kenward–Roger degrees of freedom, for
example, are refused on singular fits by mixeff’s validated-support
policy, with stable reason codes and a bootstrap alternative offered
under the same labeling rules — the route map, the refusal codes, and
the alternatives are documented in
[`vignette("inference")`](https://bbuchsbaum.github.io/mixeff/articles/inference.md).

## References

- Bates, D., Mächler, M., Bolker, B., and Walker, S. (2015). Fitting
  linear mixed-effects models using lme4. *Journal of Statistical
  Software*, 67(1), 1–48. <https://doi.org/10.18637/jss.v067.i01>
