mk_audit_fit <- function(seed = 11L) {
  set.seed(seed)
  n_subjects <- 10L
  n_per <- 4L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1, n_subjects)
  b0 <- rnorm(n_subjects, sd = 0.7)
  y <- 2 + 0.25 * x + b0[as.integer(subject)] + rnorm(length(x), sd = 0.35)
  lmm(
    y ~ x + (1 + x | subject),
    data.frame(y = y, x = x, subject = subject),
    control = mm_control(verbose = -1)
  )
}

test_that("audit() is the post-fit audit_design() surface", {
  fit <- mk_audit_fit()
  a <- audit(fit)

  expect_s3_class(a, "mm_audit")
  expect_match(a$text, "Audit Summary", fixed = TRUE)
  expect_identical(audit_design(fit)$text, a$text)
})

test_that("diagnostics() and fit_status() expose artifact status fields", {
  fit <- mk_audit_fit()
  d <- diagnostics(fit)

  expect_s3_class(d, "mm_diagnostics")
  expect_s3_class(d$table, "data.frame")
  expect_named(d$table,
               c("code", "severity", "stage", "message", "affected_terms"))
  expect_identical(fit_status(fit), fit$fit_status)

  warning_diag <- diagnostics(fit, severity = "warning")
  expect_s3_class(warning_diag, "mm_diagnostics")
  expect_true(all(warning_diag$table$severity == "warning"))

  printed <- paste(capture.output(print(d)), collapse = "\n")
  expect_match(printed, "Messages:", fixed = TRUE)
  expect_false(grepl("affected_terms.*message", printed))
})

test_that("changes() reports requested/effective/fitted transitions", {
  fit <- mk_audit_fit()
  ch <- changes(fit)

  expect_s3_class(ch, "mm_change_log")
  expect_s3_class(ch$table, "data.frame")
  expect_named(ch$table,
               c("stage", "term_id", "group", "requested", "effective",
                 "fitted", "status", "detail"))
  expect_true("semantic_ir" %in% ch$table$stage)
  expect_true("certificate_time" %in% ch$table$stage)
  expect_true(any(ch$table$status %in% c("full_rank", "reduced_rank")))
})

test_that("parameterization() exposes theta/Lambda trace rows", {
  fit <- mk_audit_fit()
  p <- parameterization(fit)

  expect_s3_class(p, "mm_theta_map")
  expect_s3_class(p$table, "data.frame")
  expect_equal(nrow(p$table), length(fit$theta))
  expect_equal(p$table$theta_value, fit$theta, tolerance = 1e-8,
               ignore_attr = TRUE)
  expect_true(all(c("theta_index", "theta_name", "lambda_row",
                    "lambda_col", "lambda_value") %in% names(p$table)))

  printed <- paste(capture.output(print(p)), collapse = "\n")
  expect_match(printed, "Full theta/Lambda columns available", fixed = TRUE)
})

test_that("roles() supports declared strings and observed fallback", {
  declared <- roles(subject = "sampled_unit", x = "fixed_condition")
  expect_s3_class(declared, "mm_roles")
  expect_identical(declared$table$origin,
                   c("declared_by_user", "declared_by_user"))
  expect_identical(declared$table$role,
                   c("sampled_unit", "fixed_condition"))

  observed <- roles(mk_audit_fit())
  expect_s3_class(observed, "mm_roles")
  expect_true("observed_from_data" %in% observed$table$origin)
  expect_true("subject" %in% observed$table$variable)
  expect_true("x" %in% observed$table$variable)
})

test_that("as_json() serializes specs and fits with raw artifact JSON", {
  fit <- mk_audit_fit()
  spec <- compile_model(y ~ x + (1 | subject), fit$model_frame)

  spec_json <- as_json(spec)
  fit_json <- as_json(fit)
  spec_payload <- jsonlite::fromJSON(spec_json, simplifyVector = FALSE)
  fit_payload <- jsonlite::fromJSON(fit_json, simplifyVector = FALSE)

  expect_identical(spec_payload$schema$schema_name, "mixeff.r_object")
  expect_identical(spec_payload$object_type, "spec")
  expect_identical(fit_payload$object_type, "fit")
  expect_identical(fit_payload$fit$fit_status, fit_status(fit))
  expect_true(nzchar(fit_payload$artifact_json))
})
