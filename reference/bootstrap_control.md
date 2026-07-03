# Fixed-effect bootstrap control

Fixed-effect bootstrap control

## Usage

``` r
bootstrap_control(
  nsim = 999L,
  seed = NULL,
  failed_refit_policy = c("exclude", "count_extreme", "abort")
)
```

## Arguments

- nsim:

  Requested bootstrap replicate count.

- seed:

  Optional integer seed. `NULL` leaves the Rust RNG seed unspecified and
  records that state in row details.

- failed_refit_policy:

  How failed refits are accounted for. Stable Rust wire labels are
  `"exclude"`, `"count_extreme"`, and `"abort"`.

## Value

A list used by `contrast(..., method = "bootstrap")`.
