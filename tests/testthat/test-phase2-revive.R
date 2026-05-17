mk_phase2_fit <- function(seed = 20L) {
  set.seed(seed)
  n_subjects <- 8L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.6)
  y <- 1.5 + 0.35 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.25)
  lmm(
    y ~ x + (1 | subject),
    data.frame(y = y, x = x, subject = subject),
    control = mm_control(verbose = -1)
  )
}

test_that("Phase 2 fits carry cache and handle metadata", {
  fit <- mk_phase2_fit()

  expect_false(fit_handle_alive(fit))
  expect_true(is.environment(fit$lazy_cache))
  expect_named(fit$schema,
               c("schema_name", "schema_version", "crate_version", "package_version"))
  expect_identical(fit$schema$schema_name, "mixedmodels.compiled_model_artifact")
})

test_that("lazy extractors rebuild fixed and random design components", {
  fit <- mk_phase2_fit()

  X <- stats::model.matrix(fit, type = "fixed")
  Z <- stats::model.matrix(fit, type = "random")
  Lambda <- getME(fit, "Lambda")

  expect_equal(dim(X), c(nobs(fit), length(fixef(fit))))
  expect_s4_class(Z, "dgCMatrix")
  expect_equal(nrow(Z), nobs(fit))
  expect_equal(ncol(Z), nrow(ranef(fit)$subject))
  expect_s4_class(Lambda, "sparseMatrix")
  expect_equal(dim(Lambda), c(ncol(Z), ncol(Z)))
  expect_true(exists("X", envir = fit$lazy_cache, inherits = FALSE))
  expect_true(exists("Z", envir = fit$lazy_cache, inherits = FALSE))

  parts <- getME(fit, c("theta", "beta", "flist", "cnms", "y"))
  expect_equal(parts$theta, fit$theta)
  expect_equal(parts$beta, fit$beta)
  expect_named(parts$flist, "subject")
  expect_named(parts$cnms, "subject")
  expect_equal(length(parts$y), nobs(fit))
})

test_that("Phase 2 report verbs expose certificate, inference, and reproducibility data", {
  fit <- mk_phase2_fit()

  rb <- random_blocks(fit)
  cert <- optimizer_certificate(fit)
  inf <- inference_table(fit)
  repro <- reproducibility(fit)

  expect_s3_class(rb, "mm_random_blocks")
  expect_true(all(c("term_id", "group", "basis", "status") %in% names(rb$table)))
  expect_equal(rb$table$group, "subject")
  expect_s3_class(cert, "mm_optimizer_certificate")
  expect_true("status" %in% cert$table$metric)
  expect_s3_class(inf, "mm_inference_table")
  expect_equal(inf$table$term, names(fixef(fit)))
  expect_true(any(inf$table$status == "available"))
  expect_true(any(is.finite(inf$table$p_value)))
  expect_true(all(inf$table$method %in% c("satterthwaite", "asymptotic_wald_z",
                                          "not_computed")))
  expect_s3_class(repro, "mm_reproducibility")
  expect_true(nrow(repro$thresholds) >= 1L)
  expect_type(is_singular(fit), "logical")
})

test_that("vcov() exposes the honest fixed-effect covariance surface", {
  fit <- mk_phase2_fit()
  V <- stats::vcov(fit)

  expect_equal(unname(diag(V)), unname(fit$std_errors)^2, tolerance = 1e-12)
  expect_equal(dim(V), c(length(fixef(fit)), length(fixef(fit))))
  expect_identical(attr(V, "mm_schema_name"),
                   "mixedmodels.fixed_effect_covariance_matrix")
  expect_identical(attr(V, "mm_schema_version"), "1.0.0")
  expect_identical(attr(V, "mm_method"), "model_based")
  expect_identical(attr(V, "mm_status"), "available")
  expect_identical(attr(V, "mm_reliability"), "high")
  expect_null(attr(V, "mm_unavailable_reason"))
  expect_true(isSymmetric(unname(V)))
  expect_true(any(abs(V[upper.tri(V)]) > 0))
})

test_that("vcov() falls back to stored standard errors when covariance artifact is absent", {
  fit <- mk_phase2_fit()
  fit$fixed_effect_vcov <- NULL
  fit$artifact$fixed_effect_covariance_matrix <- NULL

  V <- stats::vcov(fit)

  expect_equal(unname(diag(V)), unname(fit$std_errors)^2, tolerance = 1e-12)
  expect_true(all(unname(V[row(V) != col(V)]) == 0))
  expect_identical(attr(V, "mm_method"), "diagonal_from_stored_standard_errors")
  expect_identical(attr(V, "mm_status"), "unavailable")
  expect_identical(attr(V, "mm_reliability"), "not_available")
  expect_identical(attr(V, "mm_unavailable_reason"),
                   "fixed_effect_covariance_payload_unavailable")
})

test_that("rank-deficient fixed-effect covariance is explicitly unavailable", {
  set.seed(21)
  subject <- factor(rep(seq_len(8L), each = 4L))
  x <- rep(0:3, 8L)
  x2 <- x
  y <- 1 + 0.5 * x + rnorm(length(x))
  fit <- lmm(
    y ~ x + x2 + (1 | subject),
    data.frame(y = y, x = x, x2 = x2, subject = subject),
    control = mm_control(verbose = -1)
  )

  V <- stats::vcov(fit)
  details <- attr(V, "mm_details")

  expect_equal(dim(V), c(length(fixef(fit)), length(fixef(fit))))
  expect_true(all(is.na(V)))
  expect_identical(attr(V, "mm_schema_name"),
                   "mixedmodels.fixed_effect_covariance_matrix")
  expect_identical(attr(V, "mm_status"), "unavailable")
  expect_match(attr(V, "mm_unavailable_reason"), "rank_deficient")
  expect_equal(details$rank, 2)
  expect_equal(unlist(details$aliased, use.names = FALSE), "x2")
})

test_that("fixed-effect covariance payload parser rejects contract drift", {
  beta <- c("(Intercept)" = 1, x = 2)
  se <- c("(Intercept)" = 0.1, x = 0.2)
  good <- list(
    schema_name = "mixedmodels.fixed_effect_covariance_matrix",
    schema_version = "1.0.0",
    method = "model_based",
    status = "available",
    reliability = "high",
    coef_names = list("x", "(Intercept)"),
    matrix = list(c(4, 1), c(1, 9)),
    reason = NULL,
    details = list(),
    notes = list()
  )

  V <- mm_fixed_effect_vcov_from_payload(good, beta, se)
  expect_equal(dimnames(V), list(names(beta), names(beta)))
  expect_equal(matrix(as.numeric(V), nrow = 2), matrix(c(9, 1, 1, 4), nrow = 2))

  nonfinite <- good
  nonfinite$matrix <- list(c(1, NA_real_), c(NA_real_, 1))
  expect_error(
    mm_fixed_effect_vcov_from_payload(nonfinite, beta, se),
    "must be finite",
    class = "mm_schema_error"
  )

  asymmetric <- good
  asymmetric$matrix <- list(c(1, 0.2), c(0.1, 1))
  expect_error(
    mm_fixed_effect_vcov_from_payload(asymmetric, beta, se),
    "must be symmetric",
    class = "mm_schema_error"
  )

  bad_names <- good
  bad_names$coef_names <- list("(Intercept)", "z")
  expect_error(
    mm_fixed_effect_vcov_from_payload(bad_names, beta, se),
    "coefficient names do not match",
    class = "mm_schema_error"
  )

  unavailable_with_matrix <- good
  unavailable_with_matrix$status <- "unavailable"
  unavailable_with_matrix$reason <- "rank deficient"
  expect_error(
    mm_fixed_effect_vcov_from_payload(unavailable_with_matrix, beta, se),
    "must not contain a matrix",
    class = "mm_schema_error"
  )

  unavailable_without_reason <- good
  unavailable_without_reason$status <- "unavailable"
  unavailable_without_reason$matrix <- NULL
  unavailable_without_reason$reason <- NULL
  expect_error(
    mm_fixed_effect_vcov_from_payload(unavailable_without_reason, beta, se),
    "non-empty reason",
    class = "mm_schema_error"
  )
})

test_that("saveRDS/readRDS revival keeps durable extractors usable", {
  fit <- mk_phase2_fit()
  tf <- tempfile(fileext = ".rds")
  saveRDS(fit, tf)
  restored <- revive(readRDS(tf))

  expect_false(fit_handle_alive(restored))
  expect_true(is.environment(restored$lazy_cache))
  expect_equal(fixef(restored), fixef(fit))
  expect_equal(predict(restored), predict(fit), ignore_attr = TRUE)
  expect_s3_class(audit(restored), "mm_audit")
  expect_true(nrow(changes(restored)$table) >= 1L)
  expect_true(nrow(random_blocks(restored)$table) >= 1L)
})

test_that("serialized fits restart in a fresh R process", {
  testthat::skip_on_cran()

  fit <- mk_phase2_fit()
  rds <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  saveRDS(fit, rds)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(strsplit(args[[2]], .Platform$path.sep, fixed = TRUE)[[1]])",
    "library(mixeff)",
    "fit <- revive(readRDS(args[[1]]))",
    "stopifnot(!fit_handle_alive(fit))",
    "stopifnot(length(fixef(fit)) == 2L)",
    "stopifnot(length(predict(fit)) == nobs(fit))",
    "stopifnot(nrow(random_blocks(fit)$table) >= 1L)",
    "cat('cross-session OK\\n')"
  ), script)

  out <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(script, rds, paste(.libPaths(), collapse = .Platform$path.sep)),
    stdout = TRUE,
    stderr = TRUE
  )
  expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(out, collapse = "\n"))
  expect_true(any(grepl("cross-session OK", out, fixed = TRUE)))
})
