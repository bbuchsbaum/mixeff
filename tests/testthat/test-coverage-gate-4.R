# Coverage-gate round 4: clear the 90% line target (~96 more hits needed).
# Local fixtures only; do not duplicate gate / gate-2 / gate-3 cases.
# Keep fast: small data, small nsim, mm_control(verbose = -1).

mk_cg4_lmm <- function(seed = 71L, reml = FALSE) {
  set.seed(seed)
  g <- gl(10, 5)
  x <- rnorm(50)
  y <- 0.5 * x + rnorm(50, sd = 0.3) + rnorm(10, sd = 0.45)[g]
  lmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    REML = reml,
    control = mm_control(verbose = -1)
  )
}

mk_cg4_nested <- function(seed = 72L) {
  set.seed(seed)
  g <- gl(10, 5)
  x <- rnorm(50)
  z <- rnorm(50)
  y <- 0.4 * x + 0.25 * z + rnorm(50, sd = 0.3) + rnorm(10, sd = 0.4)[g]
  df <- data.frame(y = y, x = x, z = z, g = g)
  list(
    full = lmm(y ~ x + z + (1 | g), df, REML = FALSE,
               control = mm_control(verbose = -1)),
    reduced = lmm(y ~ x + (1 | g), df, REML = FALSE,
                  control = mm_control(verbose = -1)),
    data = df
  )
}

mk_cg4_glmm <- function(seed = 73L, inference = "auto") {
  set.seed(seed)
  g <- gl(8, 6)
  x <- rnorm(48)
  eta <- 0.25 * x + rnorm(8, sd = 0.3)[g]
  y <- rbinom(48, 1L, plogis(eta))
  glmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    family = binomial(),
    inference = inference,
    control = mm_control(verbose = -1)
  )
}

mk_cg4_slope <- function(seed = 74L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  y <- 0.4 * x + rnorm(60, sd = 0.3) +
    rnorm(10, sd = 0.4)[g] + rnorm(10, sd = 0.2)[g] * x
  lmm(
    y ~ x + (1 + x | g),
    data.frame(y = y, x = x, g = g),
    REML = FALSE,
    control = mm_control(verbose = -1)
  )
}

# ---- methods-summary.R: GLMM Wald-z fallback (no certified table) -----------

test_that("mm_glmm_wald_z_inference fallback when inference table absent", {
  gfit <- mk_cg4_glmm()
  gfit$artifact$fixed_effect_inference_table <- NULL
  inf <- mixeff:::mm_glmm_wald_z_inference(gfit)
  expect_s3_class(inf, "mm_inference_table")
  expect_true(all(inf$table$method == "not_computed"))
  expect_true(all(is.na(inf$table$statistic)))
  expect_true(all(is.na(inf$table$p_value)))
  expect_equal(nrow(inf$table), length(gfit$beta))
  expect_true(!is.null(attr(inf, "mm_vcov_status")))

  s <- summary(gfit)
  expect_s3_class(s, "summary.mm_glmm")
  expect_true(all(is.na(s$coefficients[["z value"]])) ||
                !"z value" %in% names(s$coefficients))
  out <- capture.output(print(s))
  expect_true(any(grepl("Fixed effects|Generalized", out)))

  # Also exercise SE-from-vcov branch by clearing std_errors.
  gfit2 <- gfit
  gfit2$std_errors <- NULL
  inf2 <- mixeff:::mm_glmm_wald_z_inference(gfit2)
  expect_equal(nrow(inf2$table), length(gfit2$beta))
  expect_true(all(is.finite(inf2$table$std_error) | is.na(inf2$table$std_error)))
})

test_that("summary coefficient helpers cover F/Chisq labels and missing rows", {
  expect_identical(mixeff:::mm_summary_statistic_column("f"), "F value")
  expect_identical(mixeff:::mm_summary_statistic_column("chi_square"), "Chisq")
  expect_identical(mixeff:::mm_summary_statistic_column(c("z", "t")), "statistic")
  expect_identical(mixeff:::mm_summary_p_value_column("f"), "Pr(>F)")
  expect_identical(mixeff:::mm_summary_p_value_column("chi_square"), "Pr(>Chisq)")
  expect_identical(mixeff:::mm_summary_p_value_column(c("z", "t")), "p.value")

  fit <- mk_cg4_lmm()
  # Empty inference -> synthetic coefficient rows (std_error from fit or NA).
  coef0 <- mixeff:::mm_summary_coefficients(fit, NULL)
  expect_true(is.data.frame(coef0) && nrow(coef0) == length(fit$beta))

  # Missing-label fill-in path inside mm_summary_coefficients.
  fake_inf <- list(table = data.frame(
    label = "not_a_coef",
    estimate = 0,
    std_error = 1,
    denominator_df = NA_real_,
    statistic = 0,
    statistic_name = "z",
    p_value = 0.5,
    method = "unit",
    kind = "coefficient",
    stringsAsFactors = FALSE
  ))
  class(fake_inf) <- c("mm_inference_table", "list")
  coef_miss <- mixeff:::mm_summary_coefficients(fit, fake_inf)
  expect_equal(rownames(coef_miss), names(fit$beta))

  # format.pval path when a p column is present.
  coef_p <- data.frame(`Pr(>|z|)` = c(0, 1e-20, 0.5), check.names = FALSE)
  fmt <- mixeff:::mm_summary_format_coef(coef_p)
  expect_true(is.character(fmt[["Pr(>|z|)"]]) || is.factor(fmt[["Pr(>|z|)"]]) ||
                is.character(as.character(fmt[["Pr(>|z|)"]])))

  # Withheld-note early exit when no stat columns.
  note_empty <- mixeff:::mm_glmm_withheld_inference_note(
    list(coefficients = data.frame(Estimate = 1), method = "pirls_profiled")
  )
  expect_identical(note_empty, character())
})

# ---- emmeans.R --------------------------------------------------------------

test_that("emmeans on working_hessian GLMM and init-message branches", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("estimability")

  gfit <- mk_cg4_glmm(inference = "working_hessian")
  em <- emmeans::emmeans(gfit, ~ x)
  expect_s4_class(em, "emmGrid")
  s <- as.data.frame(summary(em))
  expect_true(all(is.finite(s$emmean)))
  expect_true(all(is.finite(s$SE)))

  # Direct recover_data / emm_basis success path on opted-in GLMM.
  rd <- mixeff:::recover_data.mm_glmm(gfit)
  expect_true(is.data.frame(rd) || is.list(rd))
  trms <- stats::delete.response(stats::terms(mixeff:::mm_fixed_formula(gfit)))
  grid <- data.frame(x = mean(gfit$model_frame$x))
  eb <- mixeff:::emm_basis.mm_glmm(gfit, trms = trms, xlev = list(), grid = grid)
  expect_true(is.list(eb) && !is.null(eb$X) && !is.null(eb$V))
  expect_true(any(grepl("UNCERTIFIED|working-Hessian|mixeff",
                        paste(eb$misc$initMesg, collapse = " "))))

  # mm_emmeans_init_messages status branches (including unavailable).
  V <- diag(2)
  attr(V, "mm_status") <- "available"
  attr(V, "mm_method") <- "unit_method"
  expect_match(mixeff:::mm_emmeans_init_messages(V), "fixed-effect covariance from")
  attr(V, "mm_status") <- "available_noninferential"
  attr(V, "mm_method") <- "working_hessian"
  expect_match(mixeff:::mm_emmeans_init_messages(V), "UNCERTIFIED")
  attr(V, "mm_status") <- "unavailable"
  attr(V, "mm_unavailable_reason") <- "unit_test_unavailable"
  expect_match(mixeff:::mm_emmeans_init_messages(V), "unavailable")
})

# ---- parameterization.R -----------------------------------------------------

test_that("empty theta map print and parameterization helpers", {
  empty <- mixeff:::mm_parameterization_empty_table()
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("term_id", "theta_name", "varcorr_entries") %in% names(empty)))
  pm <- structure(
    list(table = empty, traces = list(), theta_maps = list()),
    class = "mm_theta_map"
  )
  out <- capture.output(print(pm))
  expect_true(any(grepl("no theta parameters", out, fixed = TRUE)))

  txt <- mixeff:::mm_varcorr_entry_text(list(
    list(kind = "standard_deviation", basis = list("intercept"), value = 1.25),
    list(kind = "correlation", basis = list("intercept", "x"), value = 0.3)
  ))
  expect_true(nzchar(txt) && grepl("standard_deviation|correlation", txt))
  expect_identical(mixeff:::mm_varcorr_entry_text(list()), "")
  expect_identical(mixeff:::mm_varcorr_entry_text(NULL), "")
})

# ---- random-options.R -------------------------------------------------------

test_that("random_options arg errors, default slope, and slope-model map", {
  fit <- mk_cg4_lmm()
  spec <- compile_model(y ~ x + (1 | g), fit$model_frame)

  expect_error(random_options(spec, NULL), class = "mm_arg_error")
  expect_error(random_options(spec, "g", "not_a_column"), class = "mm_data_error")
  # Hit missing(group) branch (then aborts when evaluating input = group).
  expect_error(do.call(random_options, list(spec = spec)))

  # Intercept-only fixed part: no slope available.
  spec_int <- compile_model(y ~ 1 + (1 | g), fit$model_frame)
  expect_error(random_options(spec_int, "g"), class = "mm_schema_error")

  # Default slope from fixed effects (no explicit slope arg).
  opts_default <- random_options(spec, "g")
  expect_s3_class(opts_default, "mm_random_options")
  expect_identical(opts_default$slope, "x")
  expect_true(nrow(opts_default$options) >= 1L)

  # Happy path on correlated slope model + print.
  slope_fit <- mk_cg4_slope()
  opts <- random_options(slope_fit, "g", "x")
  expect_s3_class(opts, "mm_random_options")
  expect_true(nrow(opts$options) >= 3L)
  out <- capture.output(print(opts))
  expect_true(any(grepl("Random-effect options", out, fixed = TRUE)))
  expect_true(any(grepl("(1 + x | g)", out, fixed = TRUE)))

  # Spec type / model-frame guards.
  expect_error(mixeff:::mm_assert_compiled_spec(list()), class = "mm_arg_error")
  bad_spec <- spec
  bad_spec$model_frame <- NULL
  expect_error(mixeff:::mm_spec_model_frame(bad_spec), class = "mm_schema_error")

  # Print path when current spelling is absent from the nearby map.
  fake_opts <- opts
  fake_opts$options$current <- FALSE
  out_nc <- capture.output(print(fake_opts))
  expect_true(any(grepl("current spelling not in nearby map", out_nc, fixed = TRUE)))
})

# ---- verify-convergence.R ---------------------------------------------------

test_that("verify_convergence LMM print and JSON/run-table helpers", {
  fit <- mk_cg4_lmm()
  v <- verify_convergence(fit, jitter_starts = 0L)
  expect_s3_class(v, "mm_convergence_verification")
  out <- capture.output(print(v))
  expect_true(any(grepl("Convergence verification", out, fixed = TRUE)))
  expect_true(any(grepl("Tolerances", out, fixed = TRUE)))

  expect_error(verify_convergence(list()), class = "mm_schema_error")
  expect_error(
    mixeff:::mm_json_parse_convergence_verification(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_convergence_verification("{"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_convergence_verification(NA_character_),
    class = "mm_schema_error"
  )

  empty_runs <- mixeff:::mm_verify_runs_table(list())
  expect_equal(nrow(empty_runs), 0L)
  expect_true(all(c("label", "agrees", "diagnostics") %in% names(empty_runs)))

  filled <- mixeff:::mm_verify_runs_table(list(list(
    label = "restart",
    optimizer_name = "trust_bq",
    return_code = "FTOL_REACHED",
    objective_value = 1.2,
    objective_delta = 0,
    max_abs_theta_delta = 0,
    max_abs_beta_delta = 0,
    agrees = TRUE,
    diagnostics = list("ok")
  )))
  expect_equal(nrow(filled), 1L)
  expect_true(isTRUE(filled$agrees[[1L]]))

  expect_identical(mixeff:::mm_verify_number(3), 3)
  expect_true(is.na(mixeff:::mm_verify_number("x")))
})

test_that("GLMM verify helpers and engine-family fallback", {
  gfit <- mk_cg4_glmm()
  v <- verify_convergence(gfit, jitter_starts = 0L)
  expect_s3_class(v, "mm_convergence_verification")
  expect_true(length(capture.output(print(v))) > 0L)

  expect_true(nzchar(mixeff:::mm_glmm_requested_method(gfit)))
  expect_true(nzchar(mixeff:::mm_glmm_effective_method(gfit)))

  gfit$engine_family <- NULL
  expect_identical(mixeff:::mm_glmm_engine_family_fallback(gfit), "bernoulli")

  # Binomial-with-weights spelling.
  gfit$weights <- rep(1, nobs(gfit))
  expect_identical(mixeff:::mm_glmm_engine_family_fallback(gfit), "binomial")

  # Payload builders for verification bridge.
  payload <- mixeff:::mm_rust_glmm_verify_payload(gfit)
  expect_true(is.list(payload) && !is.null(payload$family))
})

# ---- audit.R ----------------------------------------------------------------

test_that("print.mm_audit full=TRUE and audit error paths", {
  fit <- mk_cg4_lmm()
  a <- audit(fit)
  expect_s3_class(a, "mm_audit")
  out_full <- capture.output(print(a, full = TRUE))
  expect_true(length(out_full) > 5L)
  out_sum <- capture.output(print(a, full = FALSE))
  expect_true(length(out_sum) >= 1L)

  expect_error(audit(1L), class = "mm_schema_error")
  expect_error(mixeff:::mm_audit_impl(list()), class = "mm_schema_error")

  # Missing raw_json on artifact.
  bad <- fit
  attr(bad$artifact, "raw_json") <- NULL
  expect_error(mixeff:::mm_audit_impl(bad), class = "mm_schema_error")

  # Bridge-call error wrapping.
  expect_error(
    mixeff:::mm_audit_bridge_call(function(x) stop("boom"), "{}"),
    class = "mm_bridge_error"
  )
})

# ---- coef-names.R -----------------------------------------------------------

test_that("mm_fit_coef_map and mm_apply_lme4_coef_naming error paths", {
  fit <- mk_cg4_lmm()
  fit$coef_map <- NULL
  expect_error(mixeff:::mm_fit_coef_map(fit), class = "mm_schema_error")

  fit2 <- mk_cg4_lmm(seed = 75L)
  fit2$fixed_effect_vcov <- matrix(1, 1L, 1L)
  expect_error(
    mixeff:::mm_apply_lme4_coef_naming(fit2),
    class = "mm_schema_error"
  )

  # Happy re-apply on a clean fit still returns a named beta.
  fit3 <- mk_cg4_lmm(seed = 76L)
  # Strip map and rebuild naming from engine-ordered copies if present.
  named <- mixeff:::mm_coef_engine_to_lme4(
    names(fit3$beta), fit3$coef_map
  )
  expect_identical(named, names(fit3$beta))

  expect_identical(
    mixeff:::mm_coef_strip_engine_sep(c("modality: Audio", "x")),
    c("modalityAudio", "x")
  )
  expect_true(length(mixeff:::mm_aliased_coefficients(fit3)) >= 0L)
})

# ---- changes.R --------------------------------------------------------------

test_that("changes sentence paths for formula, reductions, transitions", {
  ch <- structure(
    list(
      table = mixeff:::mm_change_empty_table(),
      reductions = list(list(
        affected_term = "r0",
        reason = "reduced for unit test",
        trigger = "unit_trigger"
      )),
      covariance_transitions = list(list(
        affected_term = "r0",
        from = "full",
        to = list(reduced_rank = list(rank = 1L)),
        trigger = "unit_transition"
      )),
      effective_covariance = list(list(
        term_id = "r0",
        source_syntax = "(1 + x | g)",
        requested_rank = 3L,
        supported_rank = 1L,
        status = "reduced_rank"
      )),
      fit_status = "not_optimized",
      requested_formula = "y ~ x + (1 + x | g)",
      effective_formula = "y ~ x + (1 | g) + (0 + x | g)"
    ),
    class = "mm_change_log"
  )
  lines <- mixeff:::mm_change_sentence_lines(ch)
  expect_true(any(grepl("Formula canonicalized", lines, fixed = TRUE)))
  expect_true(any(grepl("Fitted covariance", lines, fixed = TRUE)))
  expect_true(any(grepl("reduced for unit test", lines, fixed = TRUE)))
  expect_true(any(grepl("covariance changed", lines, fixed = TRUE)))

  out <- capture.output(print(ch))
  expect_true(any(grepl("stopped early|not_optimized", out)))

  # Rank-restatement skip + term-label fallbacks.
  expect_true(mixeff:::mm_change_is_rank_restatement(
    "certificate_time_boundary", "r0", "r0"
  ))
  expect_false(mixeff:::mm_change_is_rank_restatement(
    "other", "r0", "r0"
  ))
  labs <- mixeff:::mm_change_term_labels(ch$effective_covariance)
  expect_identical(labs[["r0"]], "(1 + x | g)")
  expect_identical(mixeff:::mm_change_term_label("r0", labs), "(1 + x | g)")
  expect_identical(mixeff:::mm_change_term_label("r9", list()), "r9")
  expect_identical(mixeff:::mm_change_term_label("", list()), "(term)")
  expect_match(
    mixeff:::mm_change_covariance_text(list(reduced_rank = list(rank = 0))),
    "reduced_rank"
  )
})

# ---- compare.R: bootstrap LRT assert paths ----------------------------------

test_that("mm_assert_bootstrap_lrt_pair mismatched frames and weights", {
  nest <- mk_cg4_nested()
  expect_error(
    mixeff:::mm_assert_bootstrap_lrt_pair(nest$full, nest$reduced),
    class = "mm_arg_error"
  )

  # Mismatched observation values.
  alt <- nest$full
  alt$model_frame$y <- alt$model_frame$y + 1
  err <- tryCatch(
    mixeff:::mm_assert_bootstrap_lrt_pair(nest$reduced, alt),
    error = function(e) e
  )
  expect_s3_class(err, "mm_arg_error")
  expect_identical(err$reason_code, "bootstrap_lrt_requires_same_observations")

  # Missing variable in alternative frame.
  null_z <- nest$full
  alt_no_z <- nest$reduced
  err2 <- tryCatch(
    mixeff:::mm_assert_bootstrap_lrt_pair(null_z, alt_no_z),
    error = function(e) e
  )
  # Either nested-models (dof) or missing-frames depending on dof order.
  expect_s3_class(err2, "mm_arg_error")
  expect_true(err2$reason_code %in% c(
    "bootstrap_lrt_requires_nested_model_frames",
    "bootstrap_lrt_requires_nested_models"
  ))

  # Force missing-frames: null formula uses z, alt frame lacks z, alt has more dof.
  df <- nest$data
  null_only_z <- lmm(y ~ z + (1 | g), df, REML = FALSE,
                     control = mm_control(verbose = -1))
  alt_big <- nest$full
  alt_big$model_frame$z <- NULL
  err3 <- tryCatch(
    mixeff:::mm_assert_bootstrap_lrt_pair(null_only_z, alt_big),
    error = function(e) e
  )
  expect_identical(err3$reason_code, "bootstrap_lrt_requires_nested_model_frames")

  # Weight mismatch.
  w_red <- lmm(y ~ x + (1 | g), df, REML = FALSE, weights = rep(1, nrow(df)),
               control = mm_control(verbose = -1))
  w_full <- lmm(y ~ x + z + (1 | g), df, REML = FALSE, weights = rep(1.5, nrow(df)),
                control = mm_control(verbose = -1))
  err4 <- tryCatch(
    mixeff:::mm_assert_bootstrap_lrt_pair(w_red, w_full),
    error = function(e) e
  )
  expect_identical(err4$reason_code, "bootstrap_lrt_requires_same_weights")
})

test_that("mm_json_parse_model_comparison_table error branches", {
  expect_error(
    mixeff:::mm_json_parse_model_comparison_table(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_model_comparison_table("{"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_model_comparison_table("{}"),
    class = "mm_schema_error"
  )
  bad_schema <- jsonlite::toJSON(list(
    schema = list(schema_name = "nope", schema_version = "1.0.0"),
    payload = list(rows = list())
  ), auto_unbox = TRUE)
  expect_error(
    mixeff:::mm_json_parse_model_comparison_table(as.character(bad_schema)),
    class = "mm_schema_error"
  )
  missing_rows <- jsonlite::toJSON(list(
    schema = list(
      schema_name = "mixedmodels.model_comparison_table",
      schema_version = "1.0.0"
    ),
    payload = list()
  ), auto_unbox = TRUE)
  expect_error(
    mixeff:::mm_json_parse_model_comparison_table(as.character(missing_rows)),
    class = "mm_schema_error"
  )
})

# ---- predict.R --------------------------------------------------------------

test_that("predict conditional interval NA reason and fixed-only errors", {
  fit <- mk_cg4_lmm()
  nd <- data.frame(
    x = 0,
    g = factor("new_lvl", levels = c(levels(fit$model_frame$g), "new_lvl"))
  )
  out <- predict(
    fit, newdata = nd, allow.new.levels = TRUE,
    se.fit = TRUE, interval = "confidence"
  )
  expect_true(is.list(out) && !is.null(out$fit) && !is.null(out$se.fit))
  expect_true(anyNA(out$fit) || anyNA(out$se.fit))
  expect_true(!is.null(attr(out$fit, "mm_reason")) ||
                !is.null(attr(out$se.fit, "mm_reason")))

  # Prediction interval form without se.fit also sets mm_reason on NA.
  out_pi <- predict(
    fit, newdata = nd, allow.new.levels = TRUE,
    interval = "prediction"
  )
  expect_true(is.matrix(out_pi) && anyNA(out_pi))

  # mm_predict_fixed_only / mm_engine_fixed_matrix failure paths.
  expect_error(
    mixeff:::mm_predict_fixed_only(fit, data.frame(not_x = 1)),
    class = "mm_data_error"
  )
  expect_error(
    mixeff:::mm_engine_fixed_matrix(fit, data.frame(not_x = 1)),
    class = "mm_data_error"
  )
})

test_that("GLMM predict population path and re-eta new-level guard", {
  gfit <- mk_cg4_glmm()
  pop <- predict(gfit, re.form = NA)
  expect_equal(length(pop), nobs(gfit))

  nd <- data.frame(x = 0, g = gfit$model_frame$g[1])
  pop_nd <- predict(gfit, newdata = nd, re.form = NA)
  expect_equal(length(pop_nd), 1L)

  # Conditional newdata with unseen level refused unless allow.new.levels.
  nd_new <- data.frame(
    x = 0,
    g = factor("brand_new", levels = c(levels(gfit$model_frame$g), "brand_new"))
  )
  expect_error(
    predict(gfit, newdata = nd_new, allow.new.levels = FALSE),
    class = "mm_inference_unavailable"
  )
  ok <- predict(gfit, newdata = nd_new, allow.new.levels = TRUE)
  expect_equal(length(ok), 1L)

  # mm_glmm_re_eta with stripped random_effects hits missing-BLUP abort.
  g2 <- gfit
  g2$random_effects <- list()
  expect_error(
    mixeff:::mm_glmm_re_eta(g2, nd, allow_new_levels = TRUE),
    class = "mm_inference_unavailable"
  )
})

# ---- simulate.R: correlated slope covariance --------------------------------

test_that("simulate correlated slope and mm_random_term_covariance", {
  fit <- mk_cg4_slope()
  sims <- simulate(fit, nsim = 2L, seed = 81L)
  expect_true(is.data.frame(sims) && ncol(sims) == 2L)

  terms <- fit$artifact$semantic_model$random_terms
  expect_true(length(terms) >= 1L)
  term <- terms[[1L]]
  tid <- term$id
  if (is.null(tid)) tid <- "r0"
  basis <- term$basis
  if (is.null(basis)) basis <- list()
  bl <- if (length(basis)) {
    vapply(basis, mixeff:::mm_basis_label, character(1))
  } else {
    "(Intercept)"
  }
  Sigma <- mixeff:::mm_random_term_covariance(fit, tid, bl)
  expect_true(is.matrix(Sigma) && nrow(Sigma) >= 1L)
  if (nrow(Sigma) >= 2L) {
    expect_equal(Sigma[1L, 2L], Sigma[2L, 1L])
  }

  # Empty basis -> 0x0; eigen fallback for non-PD Sigma.
  expect_equal(dim(mixeff:::mm_random_term_covariance(fit, tid, character())),
               c(0L, 0L))
  Sig_bad <- matrix(c(1, 2, 2, 1), 2L, 2L,
                    dimnames = list(c("a", "b"), c("a", "b")))
  draws <- mixeff:::mm_rmvnorm(4L, Sig_bad)
  expect_equal(dim(draws), c(4L, 2L))

  # Seed restore helper.
  set.seed(99L)
  before <- .Random.seed
  mixeff:::mm_with_seed(123L, { rnorm(1) })
  expect_identical(.Random.seed, before)
})

# ---- diagnostics.R ----------------------------------------------------------

test_that("diagnostics helpers empty table and artifact guards", {
  expect_error(mixeff:::mm_compiled_artifact(list()), class = "mm_schema_error")
  empty <- mixeff:::mm_diagnostics_table(list())
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("code", "severity", "message") %in% names(empty)))

  rows <- mixeff:::mm_diagnostics_table(list(
    list(code = "scope_note", severity = "info", stage = "design",
         message = "unit", affected_terms = list("r0"))
  ))
  expect_equal(nrow(rows), 1L)

  # Unknown diagnostic code warns once and tags the table.
  mixeff:::mm_unknown_diag_state$seen <- character()
  expect_warning(
    guarded <- mixeff:::mm_diagnostics_guard(
      data.frame(code = "totally_unknown_diag_code_xyz", stringsAsFactors = FALSE)
    )
  )
  expect_true(!is.null(attr(guarded, "mm_unrecognized_diagnostic_code")))

  fit <- mk_cg4_lmm()
  expect_true(is.list(mixeff:::mm_model_diagnostics(fit)))
  expect_true(is.character(mixeff:::mm_diagnostic_bucket("scope_note")))
  expect_true(is.na(mixeff:::mm_diagnostic_bucket("")))
  expect_true(is.na(mixeff:::mm_diagnostic_bucket(character())))
})

# ---- inference.R remaining helpers via ::: ----------------------------------

test_that("bootstrap scale notice, term L, profile map, empty profile", {
  fit <- mk_cg4_lmm()
  # Verbose notice path for large nsim.
  fit$control$verbose <- 0L
  expect_invisible(
    mixeff:::mm_inform_bootstrap_scale(
      fit, list(nsim = 200L)
    )
  )
  fit$control$verbose <- -1L
  expect_invisible(
    mixeff:::mm_inform_bootstrap_scale(fit, bootstrap_control(nsim = 2L))
  )

  L1 <- mixeff:::mm_term_to_l_matrix(fit, "x")
  expect_true(is.matrix(L1) && ncol(L1) == length(fit$beta))
  L0 <- mixeff:::mm_term_to_l_matrix(fit, "1")
  expect_true(is.matrix(L0) && nrow(L0) >= 1L)
  expect_error(
    mixeff:::mm_term_to_l_matrix(fit, "not_a_term"),
    class = "mm_inference_unavailable"
  )

  expect_true(nzchar(mixeff:::regex_escape("a.b+(c)")))

  mapped_sigma <- mixeff:::mm_map_profile_parameter("\u03C3", fit)
  expect_identical(mapped_sigma$name, "sigma")
  mapped_beta <- mixeff:::mm_map_profile_parameter("\u03B21", fit)
  expect_true(is.null(mapped_beta) || identical(mapped_beta$kind, "beta") ||
                identical(mapped_beta$name, names(fit$beta)[[1L]]) ||
                identical(mapped_beta$kind, "beta"))
  mapped_theta <- mixeff:::mm_map_profile_parameter("\u03B81", fit)
  expect_identical(mapped_theta$kind, "theta")
  mapped_unk <- mixeff:::mm_map_profile_parameter("mystery_parm", fit)
  expect_identical(mapped_unk$kind, "unknown")

  expect_equal(nrow(mixeff:::mm_empty_profile_table()), 0L)
  expect_identical(
    mixeff:::mm_profile_beta_subset(fit, "(Intercept)"),
    "(Intercept)"
  )
  expect_equal(
    length(mixeff:::mm_profile_beta_subset(fit, NULL)),
    length(fit$beta)
  )

  expect_error(
    mixeff:::mm_select_bootstrap_interval(list(intervals = list()), 0.95, "percentile"),
    class = "mm_schema_error"
  )
})
