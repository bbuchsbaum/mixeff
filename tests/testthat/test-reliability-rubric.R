# Reliability grade rubric (Audit1 WI-5.5/5.6): every reliability value any
# public surface emits must be a member of the closed four-value scale the
# glossary rubric documents. P2's arbiter fixed one route that wrote a
# method name into the grade column ("asymptotic_wald"); this test keeps
# that class of bug out permanently.

mm_reliability_enum <- c("high", "moderate", "low", "not_available")

mm_collect_reliability <- function(...) {
  vals <- unlist(list(...), use.names = FALSE)
  vals <- as.character(vals)
  vals[!is.na(vals) & nzchar(vals)]
}

test_that("LMM surfaces emit only rubric grades", {
  set.seed(11)
  d <- data.frame(y = rnorm(96), x = rnorm(96),
                  g = factor(rep(1:12, each = 8)))
  d$y <- 1 + 0.4 * d$x + rep(rnorm(12, sd = 0.5), each = 8) + d$y * 0.4
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))

  vals <- mm_collect_reliability(
    inference_table(fit)$table$reliability,
    contrast(fit, c(0, 1))$table$reliability,
    attr(vcov(fit), "mm_reliability")
  )
  expect_true(length(vals) > 0L)
  expect_true(all(vals %in% mm_reliability_enum),
              info = paste("off-scale grades:",
                           paste(setdiff(vals, mm_reliability_enum),
                                 collapse = ", ")))
})

test_that("GLMM surfaces emit only rubric grades across all three routes", {
  set.seed(12)
  g <- factor(rep(1:12, each = 10))
  x <- rnorm(120)
  eta <- -0.2 + 0.7 * x + rep(rnorm(12, sd = 0.5), each = 10)
  d <- data.frame(y = rbinom(120, 1, plogis(eta)), x = x, g = g)

  profiled <- glmm(y ~ x + (1 | g), d, family = binomial(),
                   control = mm_control(verbose = -1))
  opted <- glmm(y ~ x + (1 | g), d, family = binomial(),
                inference = "working_hessian",
                control = mm_control(verbose = -1))
  joint <- glmm(y ~ x + (1 | g), d, family = binomial(),
                method = "joint_laplace", control = mm_control(verbose = -1))

  L <- c(0, 1)
  names(L) <- names(fixef(profiled))
  boot_ci <- confint(profiled, method = "bootstrap", nsim = 25L, seed = 3L)

  vals <- mm_collect_reliability(
    inference_table(profiled)$table$reliability,
    inference_table(opted)$table$reliability,
    inference_table(joint)$table$reliability,
    contrast(profiled, L)$table$reliability,
    contrast(opted, L)$table$reliability,
    contrast(joint, L)$table$reliability,
    attr(vcov(profiled), "mm_reliability"),
    attr(vcov(joint), "mm_reliability"),
    attr(boot_ci, "mm_bootstrap")$reliability$reliability
  )
  expect_true(length(vals) > 0L)
  expect_true(all(vals %in% mm_reliability_enum),
              info = paste("off-scale grades:",
                           paste(setdiff(vals, mm_reliability_enum),
                                 collapse = ", ")))
})

test_that("the bootstrap moderate floor has one shared definition", {
  expect_identical(mixeff:::mm_bootstrap_moderate_min(), 999L)
  # Below the floor: low; at the floor with finite MCSE: moderate.
  low <- mixeff:::mm_bootstrap_reliability(TRUE, 500L, 0.01)
  expect_identical(low$reliability, "low")
  mod <- mixeff:::mm_bootstrap_reliability(TRUE, 999L, 0.01)
  expect_identical(mod$reliability, "moderate")
  # The options route map consumes the SAME constant.
  expect_identical(
    mixeff:::mm_inference_options_bootstrap_lrt_reliability_reason(999L),
    "bootstrap_monte_carlo_replicates"
  )
  expect_identical(
    mixeff:::mm_inference_options_bootstrap_lrt_reliability_reason(998L),
    "bootstrap_insufficient_replicates"
  )
})

test_that("workload cost strings replaced invented seconds", {
  set.seed(13)
  d <- data.frame(y = rnorm(60), x = rnorm(60),
                  g = factor(rep(1:6, each = 10)))
  fit <- lmm(y ~ x + (1 | g), d, REML = FALSE,
             control = mm_control(verbose = -1))
  opts <- inference_options(fit, nsim = 500L)
  costs <- opts$table$approx_cost
  expect_true(any(grepl("model refits", costs, fixed = TRUE)))
  # No fabricated wall-clock estimates anywhere in the table.
  expect_false(any(grepl("~[0-9.]+s", costs)))
  # The retired route-map enum spelling must not resurface.
  expect_false(any(grepl("well_specified",
                         opts$table$expected_reliability_reason, fixed = TRUE)))
})