mk_lmm_fit_fixture <- function(seed = 10L) {
  set.seed(seed)
  n_subjects <- 10L
  n_per <- 5L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.8)
  y <- 3 + 0.4 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.3)
  data.frame(y = y, x = x, subject = subject)
}

test_that("lmm() fits an LMM and stores flat extractor fields", {
  df <- mk_lmm_fit_fixture()
  expect_silent(
    fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
  )

  expect_s3_class(fit, "mm_lmm")
  expect_s3_class(fit, "mm_fit")
  expect_s3_class(fit, "mm_compiled")
  expect_s3_class(fit$model_frame, "data.frame")
  expect_identical(names(fit$model_frame), c("y", "x", "subject"))
  expect_identical(fit$artifact$schema$schema_name,
                   "mixedmodels.compiled_model_artifact")

  expect_named(fixef(fit), c("(Intercept)", "x"))
  expect_length(fit$theta, 1L)
  expect_true(is.finite(sigma(fit)))
  expect_true(is.finite(fit$logLik))
  expect_equal(nobs(fit), nrow(df))
  expect_equal(df.residual(fit), fit$df_residual)
})

test_that("print.mm_lmm exposes artifact provenance and audit entry points", {
  df <- mk_lmm_fit_fixture()
  fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
  output <- paste(capture.output(print(fit)), collapse = "\n")

  expect_match(output, "Artifact: mixedmodels.compiled_model_artifact v1", fixed = TRUE)
  expect_match(output, "crate:", fixed = TRUE)
  expect_match(output, "Audit verbs: audit(), diagnostics(), inference_table(), model_report()", fixed = TRUE)
})

test_that("lmm() auto-prints explain_model unless verbose is -1", {
  df <- mk_lmm_fit_fixture()
  expect_output(
    fit <- lmm(y ~ x + (1 | subject), df),
    "Random effects",
    fixed = TRUE
  )
  expect_s3_class(fit, "mm_lmm")
  expect_silent(
    fit2 <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))
  )
  expect_s3_class(fit2, "mm_lmm")
})

test_that("standard extractors return stored fit quantities", {
  df <- mk_lmm_fit_fixture()
  fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "nobs"), nrow(df))
  expect_identical(attr(ll, "df"), fit$dof)
  expect_equal(AIC(fit), -2 * as.numeric(ll) + 2 * attr(ll, "df"))
  expect_equal(BIC(fit), fit$BIC)
  expect_equal(deviance(fit), fit$deviance)
  expect_identical(formula(fit), y ~ x + (1 | subject))
  expect_identical(model.frame(fit), fit$model_frame)

  expect_equal(fitted(fit), stats::predict(fit), ignore_attr = TRUE)
  expect_equal(residuals(fit), df$y - fitted(fit), tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_equal(stats::predict(fit, re.form = NA), fit$fixed_fitted,
               ignore_attr = TRUE)

  se <- stats::predict(fit, se.fit = TRUE)
  expect_named(se, c("fit", "se.fit"))
  expect_true(all(is.na(se$se.fit)))
  expect_identical(attr(se, "mm_unavailable_reason"),
                   "conditional_prediction_se_unavailable")
})

test_that("lmm() accepts positive case weights and preserves them for inference", {
  df <- mk_lmm_fit_fixture()
  df$w <- rep(c(1, 2, 0.5, 1.5, 3), length.out = nrow(df))

  fit <- lmm(y ~ x + (1 | subject), df, weights = w,
             control = mm_control(verbose = -1))
  V <- stats::vcov(fit)
  ct <- contrast(fit, c(0, 1), method = "asymptotic")
  ref <- refit(fit, df$y)

  expect_equal(fit$weights, df$w)
  expect_equal(ref$weights, df$w)
  expect_identical(attr(V, "mm_status"), "available")
  expect_equal(unname(diag(V)), unname(fit$std_errors)^2, tolerance = 1e-12)
  expect_equal(ct$table$std_error[[1L]], sqrt(V["x", "x"]), tolerance = 1e-8)
})

test_that("lmm() rejects invalid case weights", {
  df <- mk_lmm_fit_fixture()
  expect_error(
    lmm(y ~ x + (1 | subject), df, weights = rep(1, nrow(df) - 1L),
        control = mm_control(verbose = -1)),
    "one value per row",
    class = "mm_data_error"
  )
  expect_error(
    lmm(y ~ x + (1 | subject), df, weights = c(-1, rep(1, nrow(df) - 1L)),
        control = mm_control(verbose = -1)),
    "finite positive",
    class = "mm_data_error"
  )
})

test_that("random-effect and variance-component extractors are shaped like lme4 basics", {
  df <- mk_lmm_fit_fixture()
  fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))

  re <- ranef(fit)
  expect_s3_class(re, "mm_ranef")
  expect_named(re, "subject")
  expect_s3_class(re$subject, "data.frame")
  expect_identical(names(re$subject), "(Intercept)")
  expect_setequal(rownames(re$subject), levels(df$subject))

  co <- coef(fit)
  expect_s3_class(co, "mm_coef")
  expect_named(co, "subject")
  expect_equal(co$subject[["(Intercept)"]],
               re$subject[["(Intercept)"]] + fixef(fit)[["(Intercept)"]])

  vc <- VarCorr(fit)
  expect_s3_class(vc, "mm_varcorr")
  expect_true(nrow(vc$table) >= 1L)
  expect_equal(vc$residual_sd, sigma(fit))
})

test_that("revived extractor paths return typed values", {
  df <- mk_lmm_fit_fixture()
  fit <- lmm(y ~ x + (1 | subject), df, control = mm_control(verbose = -1))

  re <- ranef(fit, condVar = TRUE)
  expect_s3_class(re, "mm_ranef")
  expect_null(attr(re, "mm_unavailable_reason"))
  pv <- attr(re$subject, "postVar")
  expect_true(is.array(pv))
  expect_identical(dim(pv), c(1L, 1L, nrow(re$subject)))
  expect_true(all(is.finite(pv)))

  pred <- stats::predict(fit, newdata = df)
  expect_equal(unname(pred), unname(fitted(fit)), tolerance = 1e-8)
  expect_error(stats::predict(fit, re.form = ~(1 | subject)),
               class = "mm_inference_unavailable")
})
