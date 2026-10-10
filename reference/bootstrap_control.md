# Fixed-effect bootstrap control

Fixed-effect bootstrap control

## Usage

``` r
bootstrap_control(
  nsim = 999L,
  seed = NULL,
  failed_refit_policy = c("exclude", "count_extreme", "abort"),
  threads = 1L
)
```

## Arguments

- nsim:

  Requested bootstrap replicate count.

- seed:

  Optional integer seed. `NULL` (the default) draws a seed from R's
  random number generator when the bootstrap runs, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes the bootstrap
  reproducible; the seed used is recorded in the result.

- failed_refit_policy:

  How failed refits are accounted for. Stable Rust wire labels are
  `"exclude"`, `"count_extreme"`, and `"abort"`.

- threads:

  Number of worker threads for the replicate refits (default `1`,
  serial). See the "Threads" section.

## Value

A list used by `contrast(..., method = "bootstrap")`,
[`test_effect()`](https://bbuchsbaum.github.io/mixeff/reference/test_effect.md),
and `confint(..., method = "bootstrap")`.

## Threads

The engine simulates every replicate's response serially, in the serial
random-number order, and refits the replicates on `threads` worker
threads. The result is identical for every `threads` value; only the
wall-clock time changes. Workers never call into R: the R thread alone
polls for interrupts, so Ctrl-C / Esc still stops the run. The default
is `1` because CRAN policy limits packages to at most two threads unless
the user explicitly asks for more; passing `threads` is that request
(for example `threads = parallel::detectCores()`).

## Examples

``` r
bootstrap_control(nsim = 199, seed = 1, threads = 2)
#> $requested_replicates
#> [1] 199
#> 
#> $seed
#> [1] 1
#> 
#> $failed_refit_policy
#> [1] "exclude"
#> 
#> $threads
#> [1] 2
#> 
#> attr(,"class")
#> [1] "mm_bootstrap_control"
```
