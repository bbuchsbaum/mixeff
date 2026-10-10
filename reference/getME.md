# Extract low-level model components

`getME()` mirrors
[`lme4::getME()`](https://rdrr.io/pkg/lme4/man/getME.html) for
[`lmm()`](https://bbuchsbaum.github.io/mixeff/reference/lmm.md) and
[`glmm()`](https://bbuchsbaum.github.io/mixeff/reference/glmm.md) fits.
Every component is rebuilt R-side from quantities the fit stores (fixed
effects, theta and its map, conditional modes, fitted values, the model
frame), so it works on fits restored with
[`readRDS()`](https://rdrr.io/r/base/readRDS.html).

## Usage

``` r
getME(object, name, ...)

# S3 method for class 'mm_lmm'
getME(object, name, ...)

# S3 method for class 'mm_glmm'
getME(object, name, ...)
```

## Arguments

- object:

  A fitted `mm_lmm` or `mm_glmm` object.

- name:

  Component name, or a character vector of names, or `"ALL"`.

- ...:

  Reserved for future methods.

## Value

The requested component, or a named list for multiple names.

## Details

Supported components: `"X"`, `"Z"`, `"Zt"`, `"Ztlist"`, `"mmList"`,
`"y"`, `"mu"`, `"u"`, `"b"`, `"Gp"`, `"Tp"`, `"L"`, `"Lambda"`,
`"Lambdat"`, `"Lind"`, `"Tlist"`, `"ST"`, `"A"`, `"RX"`, `"RZX"`,
`"sigma"`, `"flist"`, `"fixef"`, `"beta"`, `"theta"`, `"REML"`,
`"is_REML"`, `"n_rtrms"`, `"n_rfacs"`, `"N"`, `"n"`, `"p"`, `"q"`,
`"p_i"`, `"l_i"`, `"q_i"`, `"k"`, `"m_i"`, `"m"`, `"cnms"`, `"devcomp"`,
`"offset"`, `"weights"`, `"lower"`, `"glmer.nb.theta"`
(negative-binomial GLMMs), and `"ALL"` (a named list of every available
component).

`"L"` is the sparse Cholesky factor (a `Matrix` `CHMfactor`, with the
same fill-reducing permutation lme4 uses) of `Lambda'Z'WZ Lambda + I`,
where `W` holds the prior weights for an LMM and the final PIRLS working
weights for a GLMM. `"RZX"`, `"RX"`, `"u"` and `"devcomp"` (`cmp` and
`dims`) are computed from the same decomposition at the fitted
parameters. `"u"` is the minimum-norm spherical mode with
`b = Lambda u`, which is the penalized optimum also when `Lambda` is
singular.

`"devfun"` is refused: the optimizer runs inside the Rust engine and no
R deviance function exists. Refusals, and requests for unknown
components, raise an `mm_inference_unavailable` condition with
`reason_code = "getme_component_unavailable"` naming the component.

## Examples

``` r
set.seed(1)
df <- data.frame(y = rnorm(60), x = rnorm(60),
                 g = factor(rep(seq_len(10), each = 6)))
fit <- lmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1))
getME(fit, "theta")
#> g.(Intercept) 
#>             0 
getME(fit, "lower")
#> g.(Intercept) 
#>             0 
getME(fit, "devcomp")$cmp
#>        ldL2       ldRX2        wrss        ussq       pwrss       drsum 
#>   0.0000000   7.9926871  43.0619244   0.0000000  43.0619244          NA 
#>        REML         dev     sigmaML   sigmaREML 
#> 155.3169348          NA   0.8471710   0.8616536 
```
