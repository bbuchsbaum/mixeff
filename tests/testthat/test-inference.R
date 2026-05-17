mk_inference_fit <- function(seed = 30L) {
  set.seed(seed)
  n_subjects <- 9L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.6)
  y <- 2 + 0.25 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.35)
  lmm(
    y ~ x + (1 | subject),
    data.frame(y = y, x = x, subject = subject),
    control = mm_control(verbose = -1)
  )
}

test_that("contrast() formats Rust fixed-effect contrast inference rows", {
  fit <- mk_inference_fit()
  ct <- contrast(fit, c(0, 1), method = "satterthwaite")

  expect_s3_class(ct, "mm_contrast")
  expect_equal(ct$table$estimate, unname(fixef(fit)[["x"]]))
  expect_true(all(is.finite(ct$table$std_error)))
  expect_true(all(is.finite(ct$table$statistic)))
  expect_true(all(is.finite(ct$table$p_value)))
  expect_identical(ct$table$method, "satterthwaite")
  expect_identical(ct$table$requested_method, "satterthwaite")
  expect_identical(ct$table$status, "available")
  expect_identical(ct$table$statistic_name, "t")
  expect_false(is.null(ct$table$details[[1L]]$contrast_family))
  expect_identical(ct$table$details[[1L]]$contrast_family$family_id, "c1")
})

test_that("contrast() preserves explicit Rust method outcomes", {
  fit <- mk_inference_fit()

  auto <- contrast(fit, c(0, 1), method = "auto")
  kr <- contrast(fit, c(0, 1), method = "kenward_roger")
  asymptotic <- contrast(fit, c(0, 1), method = "asymptotic")
  bootstrap <- contrast(fit, c(0, 1), method = "bootstrap")
  none <- contrast(fit, c(0, 1), method = "none")

  expect_identical(auto$table$requested_method, "auto")
  expect_identical(auto$table$method, "satterthwaite")
  expect_identical(auto$table$status, "available")

  expect_identical(kr$table$requested_method, "kenward_roger")
  expect_identical(kr$table$method, "kenward_roger")
  expect_identical(kr$table$status, "available")
  expect_false(identical(kr$table$method, "satterthwaite"))
  expect_true(all(is.finite(kr$table$p_value)))

  expect_identical(asymptotic$table$requested_method, "asymptotic")
  expect_identical(asymptotic$table$method, "asymptotic_wald_z")
  expect_identical(asymptotic$table$status, "available")
  expect_identical(asymptotic$table$statistic_name, "z")
  expect_true(all(is.finite(asymptotic$table$p_value)))

  expect_identical(bootstrap$table$requested_method, "bootstrap")
  expect_identical(bootstrap$table$method, "bootstrap")
  expect_identical(bootstrap$table$status, "not_assessed")
  expect_true(is.na(bootstrap$table$p_value))
  expect_match(bootstrap$table$reason, "bootstrap.*payload")

  expect_identical(none$table$requested_method, "none")
  expect_identical(none$table$method, "not_computed")
  expect_identical(none$table$status, "not_assessed")
  expect_true(is.na(none$table$p_value))
  expect_identical(none$table$reason, "inference_not_requested")
})

test_that("contrast() can request Rust fixed-effect-null bootstrap rows", {
  fit <- mk_inference_fit()
  ct <- contrast(
    fit,
    c(0, 1),
    method = "bootstrap",
    bootstrap = bootstrap_control(nsim = 30, seed = 1)
  )

  expect_identical(ct$table$method, "bootstrap")
  expect_identical(ct$table$status, "available")
  expect_true(is.finite(ct$table$p_value))
  expect_equal(ct$table$details[[1L]]$bootstrap$requested_replicates, 30)
  expect_equal(ct$table$details[[1L]]$bootstrap$successful_replicates, 30)
  expect_identical(ct$table$details[[1L]]$bootstrap$failed_refit_policy, "exclude")
  expect_identical(ct$table$details[[1L]]$bootstrap$seed_rng, "StdRng")
  expect_equal(ct$table$details[[1L]]$bootstrap$seed, 1)
  expect_identical(
    ct$table$details[[1L]]$bootstrap$null_target$covariance_policy,
    "reuse_fitted_covariance"
  )
})

test_that("contrast() preserves matrix rows, labels, and right-hand sides", {
  fit <- mk_inference_fit()
  L <- rbind(intercept = c(1, 0), slope_minus_point_one = c(0, 1))
  rhs <- c(0, 0.1)

  ct <- contrast(fit, L, rhs = rhs, method = "asymptotic")

  expect_identical(ct$table$contrast, rownames(L))
  expect_identical(ct$table$rhs, rhs)
  expect_equal(ct$table$estimate, as.numeric(L %*% fixef(fit)) - rhs)
  expect_true(all(is.finite(ct$table$std_error)))
  expect_true(all(is.finite(ct$table$statistic)))
  expect_true(all(is.finite(ct$table$p_value)))
  expect_identical(ct$table$method, rep("asymptotic_wald_z", 2L))
})

test_that("estimability consumes the Rust fixed-contrast assessment", {
  fit <- mk_inference_fit()
  L <- rbind(slope = c(0, 1))

  est <- estimability(fit, L)

  expect_s3_class(est, "mm_estimability")
  expect_identical(est$table$status, "estimable")
  expect_true(isTRUE(est$table$estimable[[1L]]))
  expect_identical(est$table$rank, 1L)
  expect_identical(est$table$requested_rank, 1L)
  expect_true(is.na(est$table$reason[[1L]]))
})

test_that("estimability reports a real engine status, never the unavailable placeholder", {
  fit <- mk_inference_fit()
  # Default L: the fixed-effect coefficient basis. Each row should carry the
  # Rust closed-enum status (estimable / not_estimable / aliased / ...).
  est <- estimability(fit)

  expect_s3_class(est, "mm_estimability")
  expect_true(all(nzchar(est$table$status)))
  expect_false(any(est$table$status == "not_assessed"))
  # Pre-C.3 placeholder reason string must no longer appear on any row.
  expect_false(any(grepl("rust_estimability_certificate_unavailable",
                         as.character(est$table$reason), fixed = TRUE)))
})

test_that("df_for_contrast pipes through the Rust inference-table df values", {
  fit <- mk_inference_fit()
  L <- rbind(slope = c(0, 1))

  df_satt <- df_for_contrast(fit, L, method = "satterthwaite")
  df_kr   <- df_for_contrast(fit, L, method = "kenward_roger")
  df_none <- df_for_contrast(fit, L, method = "none")

  expect_s3_class(df_satt, "mm_df_for_contrast")
  expect_true(all(is.finite(df_satt)))
  expect_identical(attr(df_satt, "method"), "satterthwaite")
  expect_identical(attr(df_satt, "requested_method"), "satterthwaite")

  expect_true(all(is.finite(df_kr)))
  expect_identical(attr(df_kr, "method"), "kenward_roger")

  # method = "none" still returns NA with a reason -- not_requested, not
  # unavailable; the engine never asked for df.
  expect_true(all(is.na(df_none)))
  expect_identical(attr(df_none, "method"), "not_requested")
  expect_false(is.null(attr(df_none, "mm_unavailable_reason")))
})

test_that("test_effect() and single-model anova() consume Rust term rows", {
  fit <- mk_inference_fit()

  te <- test_effect(fit, "x", method = "kenward_roger")
  av <- stats::anova(fit, type = "III", method = "bootstrap")
  kr_av <- stats::anova(fit, type = "III", method = "kenward_roger")

  expect_s3_class(te, "mm_effect_test")
  expect_identical(te$table$term, "x")
  expect_identical(te$table$method, "kenward_roger")
  expect_identical(te$table$requested_method, "kenward_roger")
  expect_identical(te$table$status, "available")
  expect_true(is.finite(te$table$p_value))
  expect_false(is.null(te$table$details[[1L]]$contrast_family))
  printed <- paste(capture.output(print(te)), collapse = "\n")
  expect_match(printed, "Effect tests:", fixed = TRUE)
  expect_match(printed, "Full audit columns available", fixed = TRUE)
  expect_false(grepl("list\\(", printed))

  expect_s3_class(av, "mm_anova")
  expect_identical(av$table$term, "x")
  expect_true(is.na(av$table$p_value))
  expect_identical(av$table$requested_method, "bootstrap")
  expect_identical(av$table$status, "not_assessed")
  expect_s3_class(stats::anova(fit, fit), "mm_model_comparison")

  expect_identical(kr_av$table$term, "x")
  expect_identical(kr_av$table$method, "kenward_roger")
  expect_identical(kr_av$table$status, "available")
  expect_true(is.finite(kr_av$table$p_value))
  expect_false(is.null(kr_av$table$details[[1L]]$kenward_roger))
})

test_that("R inference surfaces preserve Rust detail payloads", {
  fit <- mk_inference_fit()

  sm <- summary(fit, tests = "coefficients", method = "auto")
  inf <- inference_table(fit)
  expect_identical(sm$inference$table$details, inf$table$details)

  boot <- contrast(
    fit,
    c(0, 1),
    method = "bootstrap",
    bootstrap = bootstrap_control(nsim = 12, seed = 42)
  )
  boot_details <- boot$table$details[[1L]]
  expect_identical(boot_details$bootstrap$target_kind, "fixed_effect_null")
  expect_identical(boot_details$bootstrap$requested_replicates, 12L)
  expect_identical(boot_details$bootstrap$seed_rng, "StdRng")
  expect_identical(boot_details$bootstrap$null_target$covariance_policy,
                   "reuse_fitted_covariance")
  expect_identical(boot_details$contrast_family$family_id, "c1")
  expect_identical(boot_details$contrast_family$effective_rank, 1L)

  term <- test_effect(fit, "x", method = "kenward_roger")
  term_details <- term$table$details[[1L]]
  expect_identical(term_details$contrast_family$family_id, "x")
  expect_identical(term_details$contrast_family$restriction_rows, 1L)
  expect_identical(term_details$kenward_roger$restriction_rank, 1L)

  av <- stats::anova(fit, type = "III", method = "kenward_roger")
  expect_identical(av$table$details[[1L]], term_details)
})

test_that("Phase 3 prediction and covariance requests do not fabricate uncertainty", {
  fit <- mk_inference_fit()

  theta_vcov <- stats::vcov(fit, type = "theta")
  expect_true(all(is.na(theta_vcov)))
  expect_identical(attr(theta_vcov, "mm_unavailable_reason"),
                   "theta_covariance_unavailable")

  se <- stats::predict(fit, se.fit = TRUE)
  expect_true(all(is.na(se$se.fit)))
  expect_error(
    stats::predict(fit, interval = "confidence"),
    class = "mm_inference_unavailable"
  )
})

test_that("confint(method = 'wald') is labelled as uncertified asymptotic output", {
  fit <- mk_inference_fit()
  ci <- stats::confint(fit, method = "wald")

  expect_equal(nrow(ci), length(fixef(fit)))
  expect_true(all(is.finite(ci)))
  expect_s3_class(ci, "mm_confint")
  expect_identical(attr(ci, "method"),
                   "wald_asymptotic_from_stored_standard_errors")
  expect_identical(attr(ci, "status"),
                   "not_certified_by_rust_inference_contract")
})

test_that("confint(method = 'bootstrap') consumes full-model bootstrap intervals", {
  fit <- mk_inference_fit()
  ci <- stats::confint(
    fit,
    parm = "x",
    method = "bootstrap",
    bootstrap = bootstrap_control(nsim = 30, seed = 3)
  )

  expect_equal(nrow(ci), 1L)
  expect_identical(rownames(ci), "x")
  expect_true(all(is.finite(ci)))
  expect_s3_class(ci, "mm_confint")
  expect_identical(attr(ci, "method"), "bootstrap_full_model_distribution")
  expect_identical(attr(ci, "interval"), "percentile")
  expect_identical(attr(ci, "status"), "available")

  payload <- attr(ci, "bootstrap")[[1L]]
  expect_identical(payload$metadata$target$kind, "full_model_distribution")
  expect_identical(payload$metadata$requested_replicates, 30L)
  expect_identical(payload$metadata$seed_record$rng, "StdRng")
  expect_equal(length(payload$replicate_statistics),
               payload$metadata$completed_replicates)
  expect_true(any(grepl("do not certify fixed-effect hypothesis-test p-values",
                        payload$metadata$notes, fixed = TRUE)))

  printed <- paste(capture.output(print(ci)), collapse = "\n")
  expect_match(printed, "Bootstrap run:", fixed = TRUE)
  expect_match(printed, "Full bootstrap payload available", fixed = TRUE)
  expect_false(any(grepl("replicates\\$fits|attr\\(,\"bootstrap\"|\\$objective",
                         printed)))
})

test_that("inference_table() consumes Rust artifact rows when available", {
  fit <- mk_inference_fit()
  inf <- inference_table(fit)

  expect_s3_class(inf, "mm_inference_table")
  expect_true("mixedmodels.fixed_effect_inference_table" %in%
                inf$raw$schema_name)
  expect_true(all(inf$table$kind == "coefficient"))
  expect_true(any(inf$table$method %in% c("satterthwaite", "asymptotic_wald_z")))
  available <- inf$table$status == "available"
  expect_true(any(available))
  expect_true(all(is.finite(inf$table$p_value[available])))
  expect_true(all(is.finite(inf$table$statistic[available])))
  expect_true("details" %in% names(inf$table))
})

test_that("summary() renders Rust coefficient inference rows", {
  fit <- mk_inference_fit()
  sm <- summary(fit, tests = "coefficients", method = "auto")

  expect_s3_class(sm, "summary.mm_lmm")
  expect_s3_class(sm$inference, "mm_inference_table")
  expect_identical(sm$requested_method, "auto")
  expect_true(any(sm$coefficients$method %in% c("satterthwaite", "asymptotic_wald_z")))
  p_cols <- grep("^Pr\\(|^p\\.value$", names(sm$coefficients), value = TRUE)
  expect_length(p_cols, 1L)
  expect_true(any(is.finite(sm$coefficients[[p_cols]])))
})

test_that("summary() does not compute p-values missing from Rust rows", {
  fit <- mk_inference_fit()
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$p_value <- NULL
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$status <- "p_value_unavailable"
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$reason <- "unit_test_missing_p_value"

  sm <- summary(fit, tests = "coefficients")
  p_cols <- grep("^Pr\\(|^p\\.value$", names(sm$coefficients), value = TRUE)
  expect_length(p_cols, 1L)
  expect_true(is.na(sm$coefficients[[p_cols]][[1L]]))
  expect_identical(sm$inference$table$reason[[1L]], "unit_test_missing_p_value")
})

test_that("saved fits preserve Rust artifact inference rows", {
  fit <- mk_inference_fit()
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  revived <- readRDS(path)

  expect_identical(inference_table(revived)$table, inference_table(fit)$table)
})

test_that("saved fits preserve stored inference row details", {
  fit <- mk_inference_fit()
  detail <- contrast(fit, c(0, 1), method = "kenward_roger")$table$details[[1L]]
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$details <- detail

  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  revived <- readRDS(path)

  expect_identical(
    inference_table(revived)$table$details[[1L]],
    inference_table(fit)$table$details[[1L]]
  )
  expect_identical(
    inference_table(revived)$table$details[[1L]]$contrast_family$family_id,
    "c1"
  )
})

test_that("legacy fits without artifact inference table use unavailable fallback", {
  fit <- mk_inference_fit()
  fit$artifact$fixed_effect_inference_table <- NULL
  inf <- inference_table(fit)

  expect_true(all(inf$table$method == "not_computed"))
  expect_true(all(inf$table$status == "not_assessed"))
  expect_true(all(inf$table$reason == "fixed_effect_inference_table_unavailable_legacy_object"))
})
