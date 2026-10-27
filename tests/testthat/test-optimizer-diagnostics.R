test_that("certificate table preserves evidence grades and reasons", {
  for (grade in c("not_assessed", "unavailable", "approximate", "failed")) {
    evidence <- setNames(list(list(reason = "stage-specific reason")), grade)
    cert <- list(status = "not_optimized", evidence = list(
      hessian = list(quality = evidence), certification_quality = evidence))
    tbl <- mixeff:::mm_optimizer_certificate_table(cert)
    values <- setNames(tbl$value, tbl$metric)
    expect_identical(values[["optimizer_hessian"]], grade)
    expect_identical(values[["derivative_certification"]], grade)
    expect_identical(values[["optimizer_hessian_reason"]], "stage-specific reason")
  }
  cert <- list(status = "converged_interior", evidence = list(
    hessian = list(quality = "certified"), certification_quality = "certified"))
  tbl <- mixeff:::mm_optimizer_certificate_table(cert)
  expect_identical(tbl$value[tbl$metric == "optimizer_hessian"], "certified")
  expect_identical(tbl$value[tbl$metric == "derivative_certification"], "certified")
  expect_equal(nrow(mixeff:::mm_optimizer_certificate_table(list())), 0L)
})

test_that("joint optimizer and covariance Hessians report distinct stages", {
  set.seed(20261027)
  d <- expand.grid(subject = factor(seq_len(24)), item = factor(seq_len(20)))
  d$x <- runif(nrow(d), -1, 1)
  u <- rnorm(24, 0, 0.6)
  v <- rnorm(20, 0, 0.4)
  eta <- 0.2 + 0.4 * d$x + u[d$subject] + v[d$item]
  # Consume the Gaussian draws to reproduce issue8 exactly.
  d$y <- eta + rnorm(nrow(d), 0, 0.7)
  d$binary <- rbinom(nrow(d), 1, plogis(eta))
  g <- glmm(binary ~ x + (1 | subject) + (1 | item), d,
            family = binomial(), method = "joint_laplace", nAGQ = 1L,
            inference = "asymptotic",
            control = mm_control(max_feval = 5000, ftol_abs = 1e-8,
                                 ftol_rel = 1e-12, xtol_rel = 1e-8, verbose = -1))
  cert <- optimizer_certificate(g)
  expect_identical(cert$raw$status, "converged_interior")
  expect_named(cert$raw$evidence$hessian$quality, "not_assessed")
  expect_named(cert$raw$evidence$certification_quality, "unavailable")
  expect_match(cert$raw$evidence$hessian$quality$not_assessed$reason,
               "optimizer", ignore.case = TRUE)
  expect_null(cert$raw$estimator_substitution)
  out <- paste(capture.output(print(cert)), collapse = "\n")
  expect_match(out, "optimizer_hessian", fixed = TRUE)
  expect_match(out, "not_assessed", fixed = TRUE)
  expect_match(out, "derivative_certification", fixed = TRUE)
  expect_lte(max(nchar(strsplit(out, "\n", fixed = TRUE)[[1L]])),
             max(40L, getOption("width", 80L)))
  V <- vcov(g)
  expect_identical(attr(V, "mm_method"), "joint_laplace_active_hessian")
  expect_identical(attr(V, "mm_status"), "available")
  expect_identical(attr(V, "mm_reliability"), "moderate")
  notes <- paste(attr(V, "mm_notes"), collapse = " ")
  expect_match(notes, "finite", ignore.case = TRUE)
  expect_match(notes, "theta", ignore.case = TRUE)
  expect_match(notes, "inference", ignore.case = TRUE)
  restored <- unserialize(serialize(g, NULL))
  expect_identical(optimizer_certificate(restored), cert)
  expect_equal(vcov(restored), V, tolerance = 1e-12)
})
