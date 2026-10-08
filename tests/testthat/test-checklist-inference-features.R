# Regression tests for the 2026-10 pre-CRAN audit items handled on the
# checklist/inference-features branch (R layer inference/refit features).

mm_cif_sleep <- function() {
  skip_if_not_installed("lme4")
  utils::data("sleepstudy", package = "lme4", envir = environment())
  sleepstudy
}

mm_cif_fit <- function(formula = Reaction ~ Days + (1 | Subject), data = NULL,
                       ...) {
  if (is.null(data)) data <- mm_cif_sleep()
  lmm(formula, data, control = mm_control(verbose = -1), ...)
}

test_that("LMM bootstraps with seed = NULL are governed by set.seed()", {
  fit <- mm_cif_fit()
  run <- function() {
    suppressMessages(confint(fit, method = "bootstrap",
                             bootstrap = bootstrap_control(nsim = 20)))
  }
  set.seed(11)
  a <- run()
  set.seed(11)
  b <- run()
  expect_identical(unclass(a)[seq_len(4)], unclass(b)[seq_len(4)])
  set.seed(12)
  c <- run()
  expect_false(identical(unclass(a)[seq_len(4)], unclass(c)[seq_len(4)]))
  # The resolved seed reaches the engine and is recorded.
  seeds <- vapply(attr(a, "bootstrap"), function(p) {
    as.numeric(p$metadata$seed_record$seed %||% NA)
  }, numeric(1))
  expect_true(all(is.finite(seeds)))
})

test_that("bootstrap_control() wire form resolves a NULL seed from R's RNG", {
  set.seed(5)
  w1 <- mixeff:::mm_bootstrap_wire(bootstrap_control(nsim = 10))
  set.seed(5)
  w2 <- mixeff:::mm_bootstrap_wire(bootstrap_control(nsim = 10))
  expect_identical(w1$seed, w2$seed)
  expect_true(is.numeric(w1$seed) && w1$seed >= 0)
  fixed <- mixeff:::mm_bootstrap_wire(bootstrap_control(nsim = 10, seed = 7))
  expect_identical(fixed$seed, 7L)
})

test_that("simulate.mm_lmm seed attribute follows stats::simulate", {
  fit <- mm_cif_fit()
  set.seed(3)
  state <- .Random.seed
  s0 <- simulate(fit, nsim = 2)
  expect_identical(attr(s0, "seed"), state)

  s1 <- simulate(fit, nsim = 2, seed = 42)
  expect_equal(as.numeric(attr(s1, "seed")), 42)
  expect_identical(attr(attr(s1, "seed"), "kind"), as.list(RNGkind()))
  s2 <- simulate(fit, nsim = 2, seed = 42)
  expect_identical(s1[[1]], s2[[1]])

  # Same shape as simulate.lm.
  lm_sim <- simulate(stats::lm(Reaction ~ Days, mm_cif_sleep()), nsim = 1,
                     seed = 42)
  expect_identical(names(attributes(attr(lm_sim, "seed"))),
                   names(attributes(attr(s1, "seed"))))
})

test_that("drop1 and random-term LRT keep a no-intercept model no-intercept", {
  d <- mm_cif_sleep()
  d$g2 <- factor(d$Days %% 3)
  fit <- mm_cif_fit(Reaction ~ 0 + Days + g2 + (1 | Subject), d, REML = FALSE)
  reduced <- mixeff:::mm_drop_fixed_term_formula(fit, "Days")
  expect_identical(deparse1(reduced), "Reaction ~ 0 + g2 + (1 | Subject)")
  fit2 <- mm_cif_fit(Reaction ~ 0 + Days + (1 | Subject) + (0 + Days | Subject),
                     d, REML = FALSE)
  rr <- mixeff:::mm_drop_random_term_formula(fit2, 2L)
  expect_identical(deparse1(rr), "Reaction ~ 0 + Days + (1 | Subject)")
  # Intercept models are unchanged.
  fit3 <- mm_cif_fit(Reaction ~ Days + (1 | Subject), d, REML = FALSE)
  expect_identical(deparse1(mixeff:::mm_drop_fixed_term_formula(fit3, "Days")),
                   "Reaction ~ 1 + (1 | Subject)")

  # drop1's deletion log-likelihood matches lme4's no-intercept reduced fit.
  skip_if_not_installed("lme4")
  dt <- suppressMessages(drop1(fit))
  ref <- lme4::lmer(Reaction ~ 0 + g2 + (1 | Subject), d, REML = FALSE)
  row <- dt$table %||% dt
  ll <- row$logLik[row$dropped == "Days"]
  expect_equal(ll, as.numeric(logLik(ref)), tolerance = 1e-5)
})

test_that("summary() inference recomputes all coefficients in one call", {
  fit <- mm_cif_fit(Reaction ~ Days + (Days | Subject))
  calls <- 0L
  trace_env <- environment(mixeff:::mm_rust_contrast_table)
  orig <- get("mm_rust_contrast_table", envir = trace_env)
  local_mocked_bindings(
    mm_rust_contrast_table = function(...) {
      calls <<- calls + 1L
      orig(...)
    },
    .package = "mixeff"
  )
  tbl <- inference_table(fit, method = "satterthwaite")$table
  expect_identical(calls, 1L)
  expect_identical(tbl$term, names(fixef(fit)))
  skip_if_not_installed("lmerTest")
  ref <- summary(lmerTest::lmer(Reaction ~ Days + (Days | Subject),
                                mm_cif_sleep()))$coefficients
  expect_equal(tbl$df, unname(ref[, "df"]), tolerance = 1e-3)
  expect_equal(tbl$p_value, unname(ref[, "Pr(>|t|)"]), tolerance = 1e-4)
})

test_that("bridge refits warm start from the fitted theta and keep control", {
  fit <- mm_cif_fit(Reaction ~ Days + (Days | Subject))
  ctl <- jsonlite::fromJSON(mixeff:::mm_refit_control_json(fit, warm_start = TRUE))
  expect_equal(ctl$start, unname(fit$theta), tolerance = 1e-12)
  expect_identical(ctl$verbose, -1L)
  cold <- jsonlite::fromJSON(mixeff:::mm_refit_control_json(fit))
  expect_null(cold$start)
  # A user-supplied start is never overridden.
  fit$control$start <- c(1, 0, 1)
  user <- jsonlite::fromJSON(mixeff:::mm_refit_control_json(fit, warm_start = TRUE))
  expect_equal(user$start, c(1, 0, 1))
})
