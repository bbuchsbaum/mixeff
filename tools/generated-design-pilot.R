# Deterministic pilot cases for generated lme4-style parity designs.
#
# This file is intentionally sourceable, not wired into testthat yet. It gives
# future parity tests a small, reviewable generator surface before any larger
# generated-case suite is added.

mm_generated_design_pilot <- function(seed = 20260504L) {
  list(
    schema_version = 1L,
    purpose = paste(
      "Small generated-design pilot for mixeff/lme4 parity planning.",
      "Cases are deterministic and classified before being added to CI."
    ),
    default_tolerances = list(
      fixef = 1e-5,
      sigma = 1e-4,
      logLik = 1e-4,
      fitted = 1e-4,
      residuals = 1e-4,
      varcorr = 1e-3,
      ranef = 1e-4
    ),
    cases = list(
      mm_gen_case_balanced_intercept(seed + 1L),
      mm_gen_case_unbalanced_slope(seed + 2L),
      mm_gen_case_crossed_intercepts(seed + 3L),
      mm_gen_case_boundary_slope(seed + 4L)
    )
  )
}

mm_gen_case_balanced_intercept <- function(seed) {
  set.seed(seed)
  n_group <- 12L
  n_per <- 6L
  subject <- factor(rep(seq_len(n_group), each = n_per))
  x <- rep(seq(-1, 1, length.out = n_per), n_group)
  condition <- factor(rep(rep(c("control", "trained"), each = n_per / 2L),
                          n_group),
                      levels = c("control", "trained"))
  b0 <- stats::rnorm(n_group, sd = 0.45)
  y <- 2 + 0.7 * x - 0.35 * (condition == "trained") +
    b0[as.integer(subject)] + stats::rnorm(length(x), sd = 0.35)

  mm_gen_case(
    id = "gen_balanced_random_intercept",
    formula = "y ~ x + condition + (1 | subject)",
    reml = TRUE,
    expected_status = "match",
    axes = c("balanced_groups", "fixed_factor", "random_intercept"),
    data = data.frame(y, x, condition, subject),
    seed = seed,
    notes = "Balanced random-intercept case with one numeric and one categorical fixed effect."
  )
}

mm_gen_case_unbalanced_slope <- function(seed) {
  set.seed(seed)
  group_sizes <- c(4L, 5L, 6L, 7L, 8L, 4L, 6L, 8L, 5L, 7L)
  subject <- factor(rep(seq_along(group_sizes), group_sizes))
  x <- unlist(lapply(group_sizes, function(n) seq(-1, 1, length.out = n)),
              use.names = FALSE)
  b0 <- stats::rnorm(length(group_sizes), sd = 0.35)
  b1 <- stats::rnorm(length(group_sizes), sd = 0.18)
  y <- 1.4 + 0.55 * x + b0[as.integer(subject)] +
    b1[as.integer(subject)] * x + stats::rnorm(length(x), sd = 0.3)

  mm_gen_case(
    id = "gen_unbalanced_random_slope",
    formula = "y ~ x + (1 + x | subject)",
    reml = TRUE,
    expected_status = "match",
    axes = c("unbalanced_groups", "random_slope", "correlated_random_slope"),
    data = data.frame(y, x, subject),
    seed = seed,
    notes = "Unbalanced random-intercept/slope case with within-group x variation."
  )
}

mm_gen_case_crossed_intercepts <- function(seed) {
  set.seed(seed)
  n_subject <- 10L
  n_item <- 8L
  grid <- expand.grid(
    subject = factor(seq_len(n_subject)),
    item = factor(seq_len(n_item))
  )
  x <- stats::rnorm(nrow(grid))
  subject_shift <- stats::rnorm(n_subject, sd = 0.4)
  item_shift <- stats::rnorm(n_item, sd = 0.25)
  y <- 0.8 + 0.5 * x + subject_shift[as.integer(grid$subject)] +
    item_shift[as.integer(grid$item)] + stats::rnorm(nrow(grid), sd = 0.35)

  mm_gen_case(
    id = "gen_crossed_random_intercepts",
    formula = "y ~ x + (1 | subject) + (1 | item)",
    reml = TRUE,
    expected_status = "match",
    axes = c("crossed_grouping", "random_intercept", "balanced_groups"),
    data = data.frame(y, x, subject = grid$subject, item = grid$item),
    seed = seed,
    notes = "Crossed subject/item random-intercept case."
  )
}

mm_gen_case_boundary_slope <- function(seed) {
  set.seed(seed)
  n_group <- 14L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_group), each = n_per))
  x <- rep(seq(-1, 1, length.out = n_per), n_group)
  b0 <- stats::rnorm(n_group, sd = 0.45)
  b1 <- rep(0, n_group)
  y <- 1.2 + 0.6 * x + b0[as.integer(subject)] +
    b1[as.integer(subject)] * x + stats::rnorm(length(x), sd = 0.3)

  mm_gen_case(
    id = "gen_boundary_random_slope",
    formula = "y ~ x + (1 + x | subject)",
    reml = TRUE,
    expected_status = "known_fragile",
    axes = c("boundary_regime", "random_slope", "singular_fit"),
    data = data.frame(y, x, subject),
    seed = seed,
    notes = "Random-slope variance is generated as zero; parity should classify boundary behavior rather than require interior VarCorr equality."
  )
}

mm_gen_case <- function(id, formula, reml, expected_status, axes, data, seed,
                        notes) {
  stopifnot(
    is.character(id), length(id) == 1L, nzchar(id),
    is.character(formula), length(formula) == 1L, nzchar(formula),
    is.logical(reml), length(reml) == 1L,
    expected_status %in% c("match", "known_fragile", "expected_unavailable"),
    is.data.frame(data),
    is.numeric(seed), length(seed) == 1L
  )
  list(
    id = id,
    formula = formula,
    reml = reml,
    expected_status = expected_status,
    axes = axes,
    seed = as.integer(seed),
    data = data,
    notes = notes
  )
}
