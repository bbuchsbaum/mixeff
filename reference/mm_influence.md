# Leverage, Cook's distance and deletion influence for mixeff fits

[`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html) returns
the diagonal of the hat matrix of the fitted mixed model,
[`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
the Cook's distances computed from them, and
[`influence()`](https://rdrr.io/r/stats/lm.influence.html) refits the
model with each case (or each level of a grouping variable) deleted, as
`lme4`'s methods of the same names do.

## Usage

``` r
# S3 method for class 'mm_lmm'
hatvalues(model, fullHatMatrix = FALSE, ...)

# S3 method for class 'mm_glmm'
hatvalues(model, fullHatMatrix = FALSE, ...)

# S3 method for class 'mm_lmm'
cooks.distance(model, ...)

# S3 method for class 'mm_glmm'
cooks.distance(model, ...)

# S3 method for class 'mm_lmm'
influence(model, groups = NULL, data = NULL, do.coef = TRUE, ncores = 1L, ...)

# S3 method for class 'mm_glmm'
influence(model, groups = NULL, data = NULL, do.coef = TRUE, ncores = 1L, ...)

# S3 method for class 'mm_influence'
dfbeta(model, which = c("fixed", "var.cov"), ...)

# S3 method for class 'mm_influence'
dfbetas(model, ...)

# S3 method for class 'mm_influence'
cooks.distance(model, ...)
```

## Arguments

- model:

  A fitted `mm_lmm` or `mm_glmm` (or, for the
  [`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html),
  [`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html) and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  methods, an `mm_influence` result).

- fullHatMatrix:

  If `TRUE`, return the full `n x n` hat matrix instead of its diagonal.

- ...:

  Unused.

- groups:

  For [`influence()`](https://rdrr.io/r/stats/lm.influence.html): `NULL`
  (default) to delete single cases, or the name(s) of grouping column(s)
  whose levels are deleted one at a time. Columns are looked up in
  `data` if given, else in the model frame. Several names delete the
  levels of their interaction.

- data:

  Optional data frame with one row per observation (in model frame
  order) used only to look up `groups`.

- do.coef:

  If `FALSE`, [`influence()`](https://rdrr.io/r/stats/lm.influence.html)
  returns only the hat values without refitting.

- ncores:

  Number of cores for the deletion refits (via
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) on
  non-Windows platforms); default 1.

- which:

  For [`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html):
  `"fixed"` (default) for the fixed effects or `"var.cov"` for the
  variance-covariance components.

## Value

[`hatvalues()`](https://rdrr.io/r/stats/influence.measures.html) and
[`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html): a
named numeric vector (or a matrix for `fullHatMatrix = TRUE`).
[`influence()`](https://rdrr.io/r/stats/lm.influence.html): an
`mm_influence` list.

## Details

The hat matrix is that of the penalized least-squares problem at the
optimum, `H = W^{1/2} C M^{-1} C' W^{1/2}` with `C = [X, Z Lambda]` and
`M = C' W C + diag(0_p, I_q)`, where `Lambda` is the relative covariance
factor and `W` the weights. For an LMM, `W` holds the prior weights and
the result equals `lme4::hatvalues()` on the same fit. For a GLMM, `W`
holds the working (IRLS) weights at convergence, as
[`stats::hatvalues()`](https://rdrr.io/r/stats/influence.measures.html)
does for a `glm`. For GLMMs with a free dispersion parameter `phi`
(Gamma, inverse Gaussian, Gaussian with a non-identity link) the working
weights are divided by `phi`, as lme4 \>= 2.1-0's
`weights(fit, "working")` are, to match the absolute-scale `Lambda`;
lme4 instead combines the working weights inside the factorization with
the prior weights outside it, so the two differ for GLMMs (lme4 warns
that its GLMM hat matrix "may not make sense").

[`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
follows `lme4`: `D_i = (r_i / (1 - h_i))^2 h_i / (phi p)` with Pearson
residuals `r_i` (`(y - mu) sqrt(w) / sqrt(V(mu))`), dispersion `phi`
(`sigma^2` for an LMM, 1 for binomial, Poisson and negative-binomial
GLMMs, `sigma^2` for Gamma, inverse-Gaussian and Gaussian
non-identity-link GLMMs) and `p` the rank of the fixed-effect design.

[`influence()`](https://rdrr.io/r/stats/lm.influence.html) deletes one
case (default) or one level of `groups` at a time and refits with the
original settings (`REML`, family, weights, offset, method, control).
Deletion refits that fail are recorded as `NA` rows with
`converged = FALSE` rather than dropped. The result has class
`mm_influence` and the component names of lme4's `influence.merMod`
result;
[`stats::dfbeta()`](https://rdrr.io/r/stats/influence.measures.html),
[`stats::dfbetas()`](https://rdrr.io/r/stats/influence.measures.html)
and
[`stats::cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
have methods for it (Cook's distance there is the deletion version based
on the change in fixed effects).

## Examples

``` r
set.seed(1)
df <- data.frame(
  y = rnorm(40), x = rnorm(40),
  g = factor(rep(seq_len(8), each = 5))
)
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
head(hatvalues(fit))
#>          1          2          3          4          5          6 
#> 0.02742619 0.02917593 0.03494873 0.03069686 0.04457918 0.04549671 
head(cooks.distance(fit))
#>            1            2            3            4            5            6 
#> 0.0081591221 0.0006112095 0.0272905337 0.0423290513 0.0055096227 0.0174343442 
infl <- influence(fit, groups = "g")
dfbeta(infl)
#>     (Intercept)            x
#> 1 -0.0090546516  0.006611479
#> 2 -0.0009239648 -0.014618308
#> 3  0.0210233106 -0.115020427
#> 4 -0.0582753866  0.076273800
#> 5  0.0126881858 -0.003502144
#> 6  0.0625019587  0.137122729
#> 7  0.0049906926 -0.069092917
#> 8 -0.0258143784  0.016525553
cooks.distance(infl)
#>           1           2           3           4           5           6 
#> 0.002055626 0.003645842 0.217986616 0.142726249 0.003154150 0.430822891 
#>           7           8 
#> 0.078508804 0.015619794 
```
