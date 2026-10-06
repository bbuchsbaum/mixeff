status_contract_data <- function() {
  set.seed(7)
  d <- data.frame(g = factor(rep(1:12, each = 12)), x = rnorm(144))
  d$y <- rbinom(144, 1, plogis(-.3 + .8 * d$x + rep(rnorm(12, sd=.5), each=12)))
  d
}

test_that("RDS and downstream GLMM output retain actual inference status", {
  skip_if_not_installed("generics")
  for (request in c("auto", "working_hessian")) {
    fit <- glmm(y ~ x + (1 | g), status_contract_data(), family = binomial(),
                inference = request, control = mm_control(verbose = -1))
    path <- tempfile(fileext = ".rds")
    saveRDS(fit, path)
    restored <- revive(readRDS(path))
    unlink(path)
    sm <- summary(restored, tests = "coefficients")$inference$table
    td <- generics::tidy(restored, effects = "fixed", conf.int = TRUE)
    L <- matrix(c(0,1),1,dimnames=list("slope",names(fixef(restored))))
    ct <- contrast(restored, L)$table
    expect_identical(td$method, sm$method)
    expect_identical(td$status, sm$status)
    expect_identical(td$reliability, sm$reliability)
    expect_identical(ct$method, unique(td$method))
    expect_identical(ct$status, unique(td$status))
    if (request == "auto") {
      expect_true(all(is.na(td$std.error) & is.na(td$p.value) & is.na(td$conf.low)))
      expect_true(all(is.na(ct$std_error)))
      expect_error(confint(restored), class = "mm_inference_unavailable")
    } else {
      expect_true(all(is.finite(td$std.error) & is.finite(td$p.value)))
      expect_true(all(td$status == "available_noninferential"))
      expect_true(all(td$method == "wald_z_working_hessian"))
      lc <- mm_lincomb(restored, c(x=1))
      expect_identical(attr(lc, "mm_status")$status, "available_noninferential")
      expect_true(all(is.finite(confint(restored))))
      expect_output(print(summary(restored)), "UNCERTIFIED")
    }
  }
})

test_that("missing inference and covariance evidence does not enable a route", {
  fit <- glmm(y ~ x + (1 | g), status_contract_data(), family = binomial(),
              control = mm_control(verbose = -1))
  fit$artifact$fixed_effect_inference_table <- NULL
  expect_false(mixeff:::mm_glmm_inference_capability(fit)$wald)
  expect_error(confint(fit), class = "mm_inference_unavailable")
  fit$inference_request <- "working_hessian"
  # Malformed/missing covariance fails closed even after an explicit opt-in.
  fit$artifact$fixed_effect_covariance_matrix <- NULL
  fit$fixed_effect_vcov <- NULL
  fit$std_errors[] <- NA_real_
  expect_false(mixeff:::mm_glmm_inference_capability(fit)$wald)
  row <- subset(inference_options(fit)$table, method == "wald_z_working_hessian")
  expect_identical(row$expected_status, "not_assessed")
})
