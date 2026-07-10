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

  # API stabilization 2026-07-09: mm_* object with $table + $df, sibling-shaped.
  expect_s3_class(df_satt, "mm_df_for_contrast")
  expect_true(all(c("contrast", "df", "method", "requested_method", "reason")
                  %in% names(df_satt$table)))
  expect_true(all(is.finite(df_satt$df)))
  expect_identical(df_satt$method, "satterthwaite")
  expect_identical(df_satt$requested_method, "satterthwaite")

  expect_true(all(is.finite(df_kr$df)))
  expect_identical(df_kr$method, "kenward_roger")

  # method = "none" still returns NA with a reason -- not_requested, not
  # unavailable; the engine never asked for df.
  expect_true(all(is.na(df_none$df)))
  expect_identical(df_none$method, "not_requested")
  expect_false(anyNA(df_none$table$reason))
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

  # summary(method = "auto") resolves to satterthwaite on this feasible fit,
  # so compare against the same recomputed surface, not the cached table.
  sm <- summary(fit, tests = "coefficients", method = "auto")
  inf <- inference_table(fit, method = "satterthwaite")
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

  # conditional prediction SEs/intervals are now engine-certified (status
  # "available" rows from the prediction-variance payload), not fabricated
  se <- stats::predict(fit, se.fit = TRUE)
  expect_true(all(is.finite(se$se.fit) & se$se.fit > 0))
  ci <- stats::predict(fit, interval = "confidence")
  expect_true(is.matrix(ci) && all(is.finite(ci)))
  expect_true(all(ci[, "lwr"] < ci[, "fit"] & ci[, "fit"] < ci[, "upr"]))
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

test_that("summary(method = 'auto') resolves to satterthwaite on a feasible fit", {
  fit <- mk_inference_fit()
  sm <- summary(fit, tests = "coefficients", method = "auto")
  inf <- sm$inference$table

  expect_true(all(inf$method == "satterthwaite"))
  expect_true(all(inf$status == "available"))
  expect_true(all(is.finite(inf$df)))
  expect_true(all(inf$statistic_name == "t"))
  # The grade always names its closed-enum warrant (engine-authored).
  expect_true(all(inf$reliability_reason ==
                    "satterthwaite_finite_difference_approximation"))

  # df is finite, so the coefficient table keeps its df column.
  expect_true("df" %in% names(sm$coefficients))
  expect_true(all(is.finite(sm$coefficients$df)))

  # The engine's prose warrant is printed beneath the inference block.
  printed <- paste(capture.output(print(sm)), collapse = "\n")
  expect_match(printed, "Notes:", fixed = TRUE)
  expect_match(printed, "Satterthwaite denominator df", fixed = TRUE)
})

test_that("summary auto keeps the labeled asymptotic table when satterthwaite is refused", {
  # Variance pinned at zero: satterthwaite finite-difference df are refused,
  # but is_singular() does not flag the fit, so the refusal-aware path is
  # exercised (not the is_singular gate).
  set.seed(1)
  g <- factor(rep(1:6, each = 5))
  x <- rep(0:4, 6)
  df <- data.frame(y = rnorm(30), x = x, g = g)
  fit <- suppressMessages(
    lmm(y ~ x + (1 + x | g), df, control = mm_control(verbose = -1))
  )
  satt <- inference_table(fit, method = "satterthwaite")$table
  skip_if(any(satt$status == "available"),
          "satterthwaite unexpectedly available on this fit")

  sm <- summary(fit, tests = "coefficients", method = "auto")
  inf <- sm$inference$table
  expect_true(all(inf$method == "asymptotic_wald_z"))
  expect_true(all(inf$status == "available"))
  expect_true(all(inf$reliability_reason == "asymptotic_wald_z_fallback"))

  # df is undefined for z rows; the all-NA column is dropped, not printed.
  expect_false("df" %in% names(sm$coefficients))
  printed <- paste(capture.output(print(sm)), collapse = "\n")
  expect_match(printed, "labeled fallback", fixed = TRUE)
})

test_that("inference rows never fabricate reliability_reason", {
  # Engine row without the field: the parser must leave NA, not invent
  # a "not_available" string the engine never authored.
  row <- mixeff:::mm_fixed_effect_inference_row(list(label = "x"))
  expect_true(is.na(row$reliability_reason))
})

test_that("summary rendering does not compute p-values missing from Rust rows", {
  fit <- mk_inference_fit()
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$p_value <- NULL
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$status <- "p_value_unavailable"
  fit$artifact$fixed_effect_inference_table$rows[[1L]]$reason <- "unit_test_missing_p_value"

  # summary(method = "auto") no longer reads the cached table on a feasible
  # fit, so exercise the contract on the rendering layer directly: given
  # engine rows with a missing p-value, R must render NA, never fill it in.
  inf <- inference_table(fit)
  coef <- mixeff:::mm_summary_coefficients(fit, inf)
  p_cols <- grep("^Pr\\(|^p\\.value$", names(coef), value = TRUE)
  expect_length(p_cols, 1L)
  expect_true(is.na(coef[[p_cols]][[1L]]))
  expect_identical(inf$table$reason[[1L]], "unit_test_missing_p_value")
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

test_that("inference_table rows follow lme4 coefficient order on permuted designs", {
  # Regression for the engine-vs-R interaction column order: with two
  # multi-level factors the engine varies the LAST factor fastest, R the
  # first, so unordered rows silently mispair with fixef() positionally.
  set.seed(11)
  d <- expand.grid(a = factor(c("A1", "A2", "A3")),
                   b = factor(c("B1", "B2", "B3")),
                   g = factor(seq_len(8)))
  d$y <- rnorm(nrow(d))
  fit <- lmm(y ~ a * b + (1 | g), d, control = mm_control(verbose = -1))
  tbl <- inference_table(fit)$table
  coef_rows <- tbl[tbl$kind == "coefficient", , drop = FALSE]
  expect_identical(coef_rows$label, names(fixef(fit)))
  expect_equal(coef_rows$estimate, unname(fixef(fit)), tolerance = 1e-10)
})

test_that("profile.mm_lmm returns a usable mm_profile object", {
  fit <- lmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy,
             REML = FALSE, control = mm_control(verbose = -1))
  prof <- profile(fit)
  expect_s3_class(prof, "mm_profile")
  expect_true(all(c("parameter", "estimate", "lower", "upper") %in%
                    names(prof$table)))
  expect_true(all(c("(Intercept)", "Days", "sigma") %in% prof$table$parameter))
  # confint on the profile reproduces confint(fit, method = "profile")
  direct <- confint(fit, method = "profile")
  via_prof <- confint(prof)
  common <- intersect(rownames(direct), rownames(via_prof))
  expect_true(length(common) >= 2L)
  expect_equal(via_prof[common, ], direct[common, ], tolerance = 1e-10)
  expect_output(print(prof), "Profile-likelihood intervals")
  # which= filters
  p2 <- profile(fit, which = "Days")
  expect_identical(p2$table$parameter, "Days")
})

test_that("lmm refuses multivariate cbind responses with a plain error", {
  d <- data.frame(y1 = rnorm(20), y2 = rnorm(20), x = rnorm(20),
                  g = factor(rep(1:5, 4)))
  err <- tryCatch(
    lmm(cbind(y1, y2) ~ x + (1 | g), d, control = mm_control(verbose = -1)),
    error = function(e) e
  )
  expect_s3_class(err, "mm_inference_unavailable")
  expect_match(conditionMessage(err), "Multivariate responses")
  expect_match(conditionMessage(err), "own model")
})
