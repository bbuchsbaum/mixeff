# Inspect the optimizer certificate

Inspect the optimizer certificate

## Usage

``` r
optimizer_certificate(object, ...)

# S3 method for class 'mm_compiled'
optimizer_certificate(object, ...)
```

## Arguments

- object:

  A compiled `mm_spec` or fitted `mm_fit`.

- ...:

  Reserved for future methods.

## Value

An `mm_optimizer_certificate` object containing the raw certificate and
a compact table view.

## Details

The fit status records acceptance under the engine's convergence rules;
it does not imply global optimality or completion of every derivative
check. The table includes the optimizer Hessian evidence and overall
derivative certification quality, with reasons; `raw` preserves the
engine record.

For joint-Laplace GLMMs, optimizer certification and fixed-effect
covariance are separate stages. The fit-time optimizer certificate does
not include a full Hessian assessment over the joint fixed-effect (beta)
and covariance (theta) parameters, even if optimization computed a
Hessian for polishing. The later inference stage attempts a central
finite-difference Hessian of minus twice the log likelihood over beta
and interior theta coordinates. Theta coordinates too near their bounds
for central differences are omitted. The fixed-effect covariance is
twice the beta block of the inverse Hessian. Read
`attr(vcov(fit), "mm_method")`, `"mm_status"`, `"mm_reliability"`, and
`"mm_notes"` for that stage's result and any fallback. A successful
inference Hessian does not retroactively complete the optimizer
certificate.

Gradient norms are expressed in the optimizer's parameter coordinates;
they are not curvature-scaled measures of coefficient precision.
Available Newton-decrement evidence is recorded separately under
`raw$stationarity_decrement`. Neither quantity alone establishes
inferential reliability; inspect the inference method's own status and
reliability.
