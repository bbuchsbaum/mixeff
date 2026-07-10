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

test_that("audit() works on fits and specs; audit_design() is a deprecated alias", {
  fit <- mk_audit_fit()
  a <- audit(fit)

  expect_s3_class(a, "mm_audit")
  expect_match(a$text, "Audit Summary", fixed = TRUE)
  # The collapsed surface: audit_design() forwards with a deprecation note.
  expect_warning(alias <- audit_design(fit), "deprecated")
  expect_identical(alias$text, a$text)
})

test_that("print.mm_audit defaults to the upstream-rendered compact summary", {
  set.seed(19)
  # 20 levels keeps the information budget above the v0 full-covariance
  # threshold (15), so the only model change is the NotAReduction
  # canonicalization and the audit is clean.
  n_subjects <- 20L
  n_per <- 4L
  subject <- factor(rep(seq_len(n_subjects), each = n_per))
  x <- rep(seq_len(n_per) - 1L, n_subjects)
  df <- data.frame(y = rnorm(length(x)), x = x, subject = subject)
  audit <- audit(compile_model(y ~ x + (x | subject), df))

  compact_lines <- capture.output(print(audit))
  compact <- paste(compact_lines, collapse = "\n")
  full_lines <- strsplit(audit$text, "\n", fixed = TRUE)[[1L]]
  full <- paste(capture.output(print(audit, full = TRUE)), collapse = "\n")

  # The compact default is the upstream render_summary, verbatim.
  expect_identical(sub("\n$", "", audit$summary_text), compact)
  expect_match(compact, "Audit Summary:", fixed = TRUE)
  expect_match(compact, "Requested Model:", fixed = TRUE)
  expect_false(grepl("Model State:", compact, fixed = TRUE))
  expect_false(grepl("Optimizer:", compact, fixed = TRUE))

  # A correctly specified pre-fit model is clean: canonicalization is INFO,
  # and pre-fit optimizer/inference lines are "not applicable", not warnings.
  expect_match(compact, "overall [OK]", fixed = TRUE)
  expect_match(
    compact,
    "attention [OK]: no warnings or unchecked inference-critical items",
    fixed = TRUE
  )
  expect_false(grepl("NOT CHECKED", compact, fixed = TRUE))

  expect_match(full, "Model State:", fixed = TRUE)
  expect_match(full, "changes [INFO]: Diagnostic:NotAReduction", fixed = TRUE)
  expect_match(full, "not applicable before fitting", fixed = TRUE)
  expect_false(grepl("attention \\[NOT CHECKED\\]", full))
  expect_gt(nchar(full), nchar(compact))
  expect_true(all(compact_lines[nzchar(compact_lines)] %in% full_lines))
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

test_that("print(changes()) renders sentences, not the raw stage table", {
  fit <- mk_audit_fit()
  printed <- paste(capture.output(print(changes(fit))), collapse = "\n")
  expect_match(printed, "Model changes:", fixed = TRUE)
  expect_match(printed, "Stage-by-stage records available via $table.",
               fixed = TRUE)
  # the misleading formula-display row must not print
  expect_false(grepl("formula display", printed, fixed = TRUE))
  expect_false(grepl("semantic_ir", printed, fixed = TRUE))
})

test_that("print(changes()) on a converged unchanged fit says so plainly", {
  set.seed(42)
  df <- data.frame(
    t = rep(seq_len(6), 12),
    s = factor(rep(seq_len(12), each = 6))
  )
  df$y <- 2 + 0.5 * df$t + rnorm(12)[as.integer(df$s)] * 2 + rnorm(72)
  fit <- lmm(y ~ t + (1 | s), df, control = mm_control(verbose = -1))
  skip_if(!startsWith(fit$fit_status, "converged_interior"),
          "fixture did not converge interior on this build")
  printed <- paste(capture.output(print(changes(fit))), collapse = "\n")
  expect_match(printed, "none: the model was fitted as requested.",
               fixed = TRUE)
})

test_that("print(changes()) on a stopped-early fit names the optimizer state", {
  set.seed(42)
  df <- data.frame(
    t = rep(seq_len(6), 12),
    s = factor(rep(seq_len(12), each = 6))
  )
  df$y <- 2 + 0.5 * df$t + rnorm(12)[as.integer(df$s)] * 2 + rnorm(72)
  fit <- suppressWarnings(
    lmm(y ~ t + (1 | s), df,
        control = mm_control(verbose = -1, max_feval = 2))
  )
  skip_if(startsWith(fit$fit_status, "converged"),
          "max_feval = 2 unexpectedly converged on this build")
  printed <- paste(capture.output(print(changes(fit))), collapse = "\n")
  expect_match(
    printed,
    "no structural change was made; the optimizer stopped early",
    fixed = TRUE
  )
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
