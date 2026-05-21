## Tests for mm_lincomb(): Wald inference on a named linear combination
## of fixed effects, for both mm_lmm and mm_glmm fits.

mk_lincomb_lmm <- function(seed = 11L) {
  set.seed(seed)
  n_subj <- 8L
  n_per <- 6L
  subject <- factor(rep(seq_len(n_subj), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subj)
  g <- factor(rep(c("a", "b"), length.out = n_subj * n_per))
  b0 <- rnorm(n_subj, sd = 0.5)
  y <- 1 + 0.4 * x - 0.3 * (g == "b") + 0.15 * x * (g == "b") +
    b0[as.integer(subject)] + rnorm(n_subj * n_per, sd = 0.4)
  lmm(
    y ~ x * g + (1 | subject),
    data.frame(y = y, x = x, g = g, subject = subject),
    control = mm_control(verbose = -1)
  )
}

mk_lincomb_glmm <- function(seed = 23L) {
  set.seed(seed)
  n_subj <- 12L
  n_per <- 10L
  subject <- factor(rep(seq_len(n_subj), each = n_per))
  x <- rnorm(n_subj * n_per)
  g <- factor(rep(c("ctrl", "treat"), length.out = n_subj * n_per))
  b0 <- rnorm(n_subj, sd = 0.3)
  eta <- -0.2 + 0.6 * x + 0.4 * (g == "treat") - 0.5 * x * (g == "treat") +
    b0[as.integer(subject)]
  y <- rbinom(length(eta), 1L, plogis(eta))
  glmm(
    y ~ x * g + (1 | subject),
    data.frame(y = y, x = x, g = g, subject = subject),
    family = binomial(),
    method = "pirls_profiled",
    nAGQ = 1L,
    control = mm_control(verbose = -1)
  )
}

hand_wald <- function(beta, V, weights, level = 0.95, df = NA_real_) {
  bnms <- names(beta)
  w <- setNames(rep(0, length(beta)), bnms)
  w[names(weights)] <- as.numeric(weights)
  est <- sum(w * beta)
  se <- sqrt(drop(t(w) %*% V %*% w))
  if (is.finite(df) && df > 0) {
    stat <- est / se
    p <- 2 * pt(abs(stat), df = df, lower.tail = FALSE)
    q <- qt((1 + level) / 2, df = df)
  } else {
    stat <- est / se
    p <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    q <- qnorm((1 + level) / 2)
  }
  list(estimate = est, std_error = se, statistic = stat, p_value = p,
       lower = est - q * se, upper = est + q * se)
}

test_that("mm_lincomb() reproduces hand-rolled Wald z on mm_glmm", {
  fit <- mk_lincomb_glmm()
  beta <- as.numeric(fixef(fit)); names(beta) <- names(fixef(fit))
  V <- as.matrix(unclass(vcov(fit)))
  dimnames(V) <- list(names(beta), names(beta))

  ## Pick a two-coefficient interaction lincomb that depends on the
  ## off-diagonal entries of V (this is the manuscript's DiD pattern).
  inter_name <- grep(":", names(beta), value = TRUE)[1L]
  main_name  <- grep(":", names(beta), value = TRUE, invert = TRUE)
  main_name  <- setdiff(main_name, "(Intercept)")[1L]
  weights <- c(0.7, 1.3); names(weights) <- c(main_name, inter_name)

  ref <- hand_wald(beta, V, weights)
  out <- mm_lincomb(fit, weights)

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_identical(out$statistic_name, "z")
  expect_true(is.na(out$df))
  expect_equal(out$estimate,   ref$estimate,   tolerance = 1e-12)
  expect_equal(out$std_error,  ref$std_error,  tolerance = 1e-12)
  expect_equal(out$statistic,  ref$statistic,  tolerance = 1e-10)
  expect_equal(out$p_value,    ref$p_value,    tolerance = 1e-10)
  expect_equal(out$lower,      ref$lower,      tolerance = 1e-12)
  expect_equal(out$upper,      ref$upper,      tolerance = 1e-12)
})

test_that("mm_lincomb() exposes the underlying vcov status as an attribute", {
  fit <- mk_lincomb_glmm()
  out <- mm_lincomb(fit, c("x" = 1))
  st <- attr(out, "mm_status")
  expect_true(is.list(st))
  expect_true(all(c("status", "method", "reliability", "reason") %in% names(st)))
  expect_identical(st$status, "available")
})

test_that("mm_lincomb() with method='asymptotic' on mm_lmm matches hand-rolled Wald z", {
  fit <- mk_lincomb_lmm()
  beta <- as.numeric(fixef(fit)); names(beta) <- names(fixef(fit))
  V <- as.matrix(unclass(vcov(fit)))
  dimnames(V) <- list(names(beta), names(beta))

  g_name <- grep("^g", names(beta), value = TRUE)[1L]
  weights <- setNames(c(1, -0.5), c("x", g_name))
  ref <- hand_wald(beta, V, weights)
  out <- mm_lincomb(fit, weights, method = "asymptotic")

  expect_identical(out$statistic_name, "z")
  expect_true(is.na(out$df))
  expect_equal(out$estimate,  ref$estimate,  tolerance = 1e-12)
  expect_equal(out$std_error, ref$std_error, tolerance = 1e-12)
  expect_equal(out$p_value,   ref$p_value,   tolerance = 1e-10)
})

test_that("mm_lincomb() default on mm_lmm uses Satterthwaite df via df_for_contrast()", {
  fit <- mk_lincomb_lmm()
  bnms <- names(fixef(fit))
  weights <- setNames(1, "x")
  out <- mm_lincomb(fit, weights) # default method="auto"
  expect_identical(out$statistic_name, "t")
  expect_true(is.finite(out$df))
  expect_true(out$df > 0)
  ## Sanity: matches contrast()'s Satterthwaite output on the same lincomb
  L <- matrix(0, nrow = 1L, ncol = length(bnms))
  L[, which(bnms == "x")] <- 1
  ct <- contrast(fit, L, method = "satterthwaite")
  expect_equal(out$estimate,  ct$table$estimate,  tolerance = 1e-10)
  expect_equal(out$std_error, ct$table$std_error, tolerance = 1e-10)
  expect_equal(out$statistic, ct$table$statistic, tolerance = 1e-10)
  expect_equal(out$df,        ct$table$df,        tolerance = 1e-6)
  expect_equal(out$p_value,   ct$table$p_value,   tolerance = 1e-8)
})

test_that("mm_lincomb() accepts named list and 1-row data.frame", {
  fit <- mk_lincomb_glmm()
  g_name <- grep("^g", names(fixef(fit)), value = TRUE)[1L]
  vec  <- setNames(c(0.5, 1), c("x", g_name))
  out_vec  <- mm_lincomb(fit, vec)
  out_list <- mm_lincomb(fit, as.list(vec))
  out_df   <- mm_lincomb(fit, as.data.frame(rbind(vec), check.names = FALSE))
  expect_equal(out_vec$estimate,  out_list$estimate, tolerance = 1e-12)
  expect_equal(out_vec$estimate,  out_df$estimate,   tolerance = 1e-12)
  expect_equal(out_vec$std_error, out_list$std_error, tolerance = 1e-12)
})

test_that("mm_lincomb() rejects malformed weights", {
  fit <- mk_lincomb_glmm()
  expect_error(mm_lincomb(fit, NULL),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c(1, 2, 3)),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c("x" = NA_real_)),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c("x" = 1, "x" = 2)),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c("not_a_coef" = 1)),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, data.frame(x = c(1, 2), gtreat = c(1, 1))),
               class = "mm_arg_error")
})

test_that("mm_lincomb() rejects unsupported method on mm_glmm", {
  fit <- mk_lincomb_glmm()
  expect_error(mm_lincomb(fit, c("x" = 1), method = "satterthwaite"),
               class = "mm_arg_error")
  expect_error(mm_lincomb(fit, c("x" = 1), method = "kenward_roger"),
               class = "mm_arg_error")
})

test_that("mm_lincomb() default method errors on non-fit input", {
  expect_error(mm_lincomb(list(), c("x" = 1)),
               class = "mm_arg_error")
})

test_that("mm_lincomb() level argument moves the CI as expected", {
  fit <- mk_lincomb_glmm()
  out_95 <- mm_lincomb(fit, c("x" = 1), level = 0.95)
  out_99 <- mm_lincomb(fit, c("x" = 1), level = 0.99)
  expect_lt(out_95$lower, out_99$upper)  ## sanity
  ## 99% interval must be wider than 95%
  width_95 <- out_95$upper - out_95$lower
  width_99 <- out_99$upper - out_99$lower
  expect_gt(width_99, width_95)
  ## Estimate and SE unchanged
  expect_equal(out_95$estimate,  out_99$estimate,  tolerance = 1e-12)
  expect_equal(out_95$std_error, out_99$std_error, tolerance = 1e-12)
})
