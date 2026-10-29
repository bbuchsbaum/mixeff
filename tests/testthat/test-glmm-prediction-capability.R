test_that("population GLMM Wald uncertainty respects inference capability", {
  set.seed(11)
  d <- data.frame(x = rnorm(240), g = factor(rep(1:12, each = 20)))
  d$y <- rbinom(240, 1, plogis(.2 + .6 * d$x))
  fit <- glmm(y ~ x + (1 | g), d, family = binomial(),
               method = "pirls_profiled", control = mm_control(verbose = -1))
  expect_false(mixeff:::mm_glmm_inference_capability(fit)$wald)
  for (type in c("link", "response")) {
    ci <- predict(fit, re.form = NA, type = type, interval = "confidence")
    expect_true(all(is.finite(ci[, "fit"])))
    expect_true(all(is.na(ci[, c("lwr", "upr")])))
    expect_identical(attr(ci, "mm_status")$method, "not_computed")
    expect_true(nzchar(attr(ci, "mm_reason")))
    se <- predict(fit, re.form = NA, type = type, se.fit = TRUE)$se.fit
    expect_true(all(is.na(se)))
    expect_identical(attr(se, "mm_status")$method, "not_computed")
  }
  fit$inference_request <- "working_hessian"
  ci <- predict(fit, re.form = NA, interval = "confidence")
  expect_true(all(is.finite(ci)))
  expect_identical(attr(ci, "mm_status")$method, "wald_z_working_hessian")
  expect_identical(attr(ci, "mm_status")$status, "available_noninferential")
  # The provenance also survives newdata missing-row expansion.
  nd <- d[1:3, ]; nd$x[2] <- NA
  out <- predict(fit, nd, re.form = NA, interval = "confidence", se.fit = TRUE)
  expect_true(all(is.na(out$fit[2, ])))
  expect_identical(attr(out$fit, "mm_status")$method, "wald_z_working_hessian")
  expect_identical(attr(out$se.fit, "mm_status")$status, "available_noninferential")
})
