# Coverage-gate tests: exercise under-hit print methods, arg validation,
# default methods, and a few high-ROI happy paths. Keep fixtures small and
# avoid MIXEFF_RUN_SLOW_PARITY / aphantasia flags.

mk_cov_lmm <- function(seed = 41L, reml = TRUE) {
  set.seed(seed)
  g <- gl(20, 5)
  x <- rnorm(100)
  y <- 0.5 * x + rnorm(100, sd = 0.4) + rnorm(20, sd = 0.6)[g]
  lmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    REML = reml,
    control = mm_control(verbose = -1)
  )
}

mk_cov_lmm_nested <- function(seed = 42L, reml = FALSE) {
  set.seed(seed)
  g <- gl(16, 5)
  x <- rnorm(80)
  z <- rnorm(80)
  y <- 0.4 * x + 0.3 * z + rnorm(80, sd = 0.35) + rnorm(16, sd = 0.5)[g]
  df <- data.frame(y = y, x = x, z = z, g = g)
  list(
    full = lmm(y ~ x + z + (1 | g), df, REML = reml,
               control = mm_control(verbose = -1)),
    reduced = lmm(y ~ x + (1 | g), df, REML = reml,
                  control = mm_control(verbose = -1)),
    data = df
  )
}

mk_cov_lmm_uncorrelated <- function(seed = 43L) {
  set.seed(seed)
  g <- gl(12, 6)
  x <- rnorm(72)
  y <- 0.5 * x + rnorm(72, sd = 0.4) +
    rnorm(12, sd = 0.5)[g] + rnorm(12, sd = 0.3)[g] * x
  lmm(
    y ~ x + (1 | g) + (0 + x | g),
    data.frame(y = y, x = x, g = g),
    control = mm_control(verbose = -1)
  )
}

mk_cov_glmm <- function(seed = 44L) {
  set.seed(seed)
  g <- gl(12, 8)
  x <- rnorm(96)
  eta <- 0.2 * x + rnorm(12, sd = 0.4)[g]
  y <- rbinom(96, 1L, plogis(eta))
  glmm(
    y ~ x + (1 | g),
    data.frame(y = y, x = x, g = g),
    family = binomial(),
    control = mm_control(verbose = -1)
  )
}

mk_cov_between_data <- function(seed = 45L) {
  set.seed(seed)
  g <- gl(10, 5)
  # Predictor constant within each group -> observed_between_group role.
  x <- rnorm(10)[as.integer(g)]
  y <- 0.5 * x + rnorm(50, sd = 0.4) + rnorm(10, sd = 0.5)[as.integer(g)]
  data.frame(y = y, x = x, g = g)
}

# ---- Print methods ----------------------------------------------------------

test_that("print.mm_lmm and extractors produce character output", {
  fit <- mk_cov_lmm()
  out_fit <- capture.output(print(fit))
  expect_true(is.character(out_fit) && length(out_fit) > 0L)
  expect_match(paste(out_fit, collapse = "\n"), "Fixed effects")

  out_vc <- capture.output(print(VarCorr(fit)))
  expect_true(is.character(out_vc) && length(out_vc) > 0L)
  expect_match(paste(out_vc, collapse = "\n"), "Variance components")

  out_re <- capture.output(print(ranef(fit)))
  expect_true(is.character(out_re) && length(out_re) > 0L)
  expect_match(paste(out_re, collapse = "\n"), "Random effects")

  out_coef <- capture.output(print(coef(fit)))
  expect_true(is.character(out_coef) && length(out_coef) > 0L)
  expect_match(paste(out_coef, collapse = "\n"), "Conditional coefficients")
})

test_that("print.mm_glmm produces character output", {
  fit <- mk_cov_glmm()
  out <- capture.output(print(fit))
  expect_true(is.character(out) && length(out) > 0L)
  expect_match(paste(out, collapse = "\n"), "Generalized linear mixed model")
  expect_match(paste(out, collapse = "\n"), "Family/link")
})

test_that("print.mm_anova and print.mm_drop1 produce character output", {
  fit <- mk_cov_lmm()
  av <- anova(fit)
  out_av <- capture.output(print(av))
  expect_true(is.character(out_av) && length(out_av) > 0L)
  expect_match(paste(out_av, collapse = "\n"), "analysis of fixed effects")

  d1 <- drop1(fit, test = "Chisq")
  out_d1 <- capture.output(print(d1))
  expect_true(is.character(out_d1) && length(out_d1) > 0L)
  expect_match(paste(out_d1, collapse = "\n"), "Single-term deletion")
})

test_that("print.mm_parametric_bootstrap and print.mm_model_comparison work", {
  nest <- mk_cov_lmm_nested(reml = FALSE)
  boot <- parametric_bootstrap(nest$reduced, nest$full, nsim = 2L, seed = 1L)
  out_boot <- capture.output(print(boot))
  expect_true(is.character(out_boot) && length(out_boot) > 0L)
  expect_match(paste(out_boot, collapse = "\n"), "Parametric bootstrap")

  cmp <- compare(nest$reduced, nest$full)
  out_cmp <- capture.output(print(cmp))
  expect_true(is.character(out_cmp) && length(out_cmp) > 0L)
  expect_match(paste(out_cmp, collapse = "\n"), "Model comparison")
})

test_that("print.mm_roles covers empty and declared roles", {
  out_empty <- capture.output(print(roles()))
  expect_true(is.character(out_empty) && length(out_empty) > 0L)
  expect_match(paste(out_empty, collapse = "\n"), "none declared")

  out_decl <- capture.output(print(roles(subject = "sampled_unit")))
  expect_true(is.character(out_decl) && length(out_decl) > 0L)
  expect_match(paste(out_decl, collapse = "\n"), "sampled_unit")
})

test_that("print.mm_spec for compile_model produces character output", {
  df <- mk_cov_between_data()
  spec <- compile_model(y ~ x + (1 | g), df)
  out <- capture.output(print(spec))
  expect_true(is.character(out) && length(out) > 0L)
  expect_match(paste(out, collapse = "\n"), "<mm_spec>")
})

test_that("print reproducibility / random_blocks / optimizer_certificate", {
  fit <- mk_cov_lmm()
  out_repro <- capture.output(print(reproducibility(fit)))
  expect_true(is.character(out_repro) && length(out_repro) > 0L)

  out_blocks <- capture.output(print(random_blocks(fit)))
  expect_true(is.character(out_blocks) && length(out_blocks) > 0L)
  expect_match(paste(out_blocks, collapse = "\n"), "Random-effect blocks")

  out_cert <- capture.output(print(optimizer_certificate(fit)))
  expect_true(is.character(out_cert) && length(out_cert) > 0L)
})

test_that("print reporting_table, inference_table, audit, diagnostics, options", {
  fit <- mk_cov_lmm()

  rt <- reporting_table(fit, "overview")
  out_rt <- capture.output(print(rt))
  expect_true(is.character(out_rt) && length(out_rt) > 0L)

  out_inf <- capture.output(print(inference_table(fit)))
  expect_true(is.character(out_inf) && length(out_inf) > 0L)

  out_audit <- capture.output(print(audit(fit)))
  expect_true(is.character(out_audit) && length(out_audit) > 0L)

  out_diag <- capture.output(print(diagnostics(fit)))
  expect_true(is.character(out_diag) && length(out_diag) > 0L)

  out_opt <- capture.output(print(inference_options(fit)))
  expect_true(is.character(out_opt) && length(out_opt) > 0L)
  expect_match(paste(out_opt, collapse = "\n"), "Inference options")
})

test_that("print.mm_contrast and print.mm_confint produce character output", {
  fit <- mk_cov_lmm()
  ct <- contrast(fit, c(0, 1), method = "asymptotic")
  out_ct <- capture.output(print(ct))
  expect_true(is.character(out_ct) && length(out_ct) > 0L)
  expect_match(paste(out_ct, collapse = "\n"), "Fixed-effect contrasts")

  ci <- confint(fit, method = "asymptotic")
  out_ci <- capture.output(print(ci))
  expect_true(is.character(out_ci) && length(out_ci) > 0L)
  expect_match(paste(out_ci, collapse = "\n"), "Confidence intervals")
})

test_that("print.mm_effect_test produces character output", {
  fit <- mk_cov_lmm()
  te <- test_effect(fit, "x", method = "asymptotic")
  out <- capture.output(print(te))
  expect_true(is.character(out) && length(out) > 0L)
  expect_match(paste(out, collapse = "\n"), "Effect tests")
})

# ---- Arg validation ---------------------------------------------------------

test_that("bootstrap_control rejects invalid nsim and seed", {
  expect_error(bootstrap_control(nsim = 0), class = "mm_arg_error")
  expect_error(bootstrap_control(nsim = -3), class = "mm_arg_error")
  expect_error(bootstrap_control(seed = -1), class = "mm_arg_error")
  expect_error(bootstrap_control(seed = c(1, 2)), class = "mm_arg_error")
})

test_that("test_effect rejects bad term arguments", {
  fit <- mk_cov_lmm()
  expect_error(test_effect(fit, character()), class = "mm_arg_error")
  expect_error(test_effect(fit, 1L), class = "mm_arg_error")
  expect_error(test_effect(fit, "not_a_term"), class = "mm_arg_error")
})

test_that("contrast rejects bad L", {
  fit <- mk_cov_lmm()
  expect_error(contrast(fit, NULL), class = "mm_arg_error")
  expect_error(contrast(fit, "bad"), class = "mm_arg_error")
  expect_error(contrast(fit, c(1, 0, 0)), class = "mm_arg_error")
})

test_that("confint rejects bad level and parm", {
  fit <- mk_cov_lmm()
  expect_error(confint(fit, level = 0), class = "mm_arg_error")
  expect_error(confint(fit, level = 1), class = "mm_arg_error")
  expect_error(confint(fit, level = c(0.9, 0.95)), class = "mm_arg_error")
  expect_error(confint(fit, parm = "not_a_coef"), class = "mm_arg_error")
})

test_that("compare refuses incompatible objects", {
  nest <- mk_cov_lmm_nested(reml = FALSE)
  expect_error(compare(nest$full, 1), class = "mm_arg_error")

  other <- mk_cov_lmm(seed = 99L, reml = FALSE)
  # Same response name but different nobs.
  expect_error(compare(nest$full, other), class = "mm_arg_error")
})

test_that("parametric_bootstrap validates arguments", {
  nest <- mk_cov_lmm_nested(reml = FALSE)
  expect_error(
    parametric_bootstrap(1, nest$full, nsim = 2),
    class = "mm_arg_error"
  )
  expect_error(
    parametric_bootstrap(nest$reduced, nest$full, nsim = 0),
    class = "mm_arg_error"
  )

  reml_nest <- mk_cov_lmm_nested(reml = TRUE)
  expect_error(
    parametric_bootstrap(reml_nest$reduced, reml_nest$full, nsim = 2, seed = 1),
    class = "mm_inference_unavailable"
  )
})

test_that("roles() rejects unnamed and empty role strings", {
  expect_error(roles("sampled_unit"), class = "mm_arg_error")
  expect_error(roles(subject = ""), class = "mm_arg_error")
  expect_error(roles(subject = NA_character_), class = "mm_arg_error")
  expect_error(roles(subject = c("a", "b")), class = "mm_arg_error")
})

test_that("glmm validates family and nAGQ", {
  df <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40), g = gl(8, 5))
  expect_error(
    glmm(y ~ x + (1 | g), df, control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = "binomial",
         control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(), nAGQ = 0,
         control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
  expect_error(
    glmm(y ~ x + (1 | g), df, family = binomial(),
         method = "joint_laplace", nAGQ = 5,
         control = mm_control(verbose = -1)),
    class = "mm_arg_error"
  )
})

test_that("summary(method = 'bootstrap') refuses with typed error", {
  fit <- mk_cov_lmm()
  expect_error(
    summary(fit, method = "bootstrap"),
    class = "mm_inference_unavailable"
  )
})

test_that("reporting_table rejects unknown section", {
  fit <- mk_cov_lmm()
  expect_error(
    reporting_table(fit, section = "not_a_section"),
    class = "mm_schema_error"
  )
  expect_error(
    reporting_table(fit, section = 1L),
    class = "mm_schema_error"
  )
})

test_that("default methods refuse non-mixeff objects", {
  expect_error(fixef(1), class = "mm_arg_error")
  expect_error(ranef(1), class = "mm_arg_error")
  expect_error(VarCorr(1), class = "mm_arg_error")
  expect_error(getME(1, "x"), class = "mm_arg_error")
  expect_error(revive(1), class = "mm_arg_error")
  expect_error(refit(1, 1:3), class = "mm_arg_error")
  expect_error(is_singular(1), class = "mm_arg_error")
  expect_error(random_blocks(1), class = "mm_arg_error")
})

test_that("audit.default refuses non-compiled objects", {
  expect_error(audit(1), class = "mm_schema_error")
})

# ---- Happy paths ------------------------------------------------------------

test_that("ranef(glmm, condVar=TRUE) returns unavailable postVar path", {
  fit <- mk_cov_glmm()
  re <- ranef(fit, condVar = TRUE)
  expect_identical(
    attr(re, "mm_unavailable_reason"),
    "random_effect_conditional_variance_unavailable_for_glmm"
  )
  expect_true(is.array(attr(re[[1L]], "postVar")))
  expect_true(all(is.na(attr(re[[1L]], "postVar"))))
  out <- capture.output(print(re))
  expect_true(is.character(out) && length(out) > 0L)
})

test_that("roles on between-group predictor marks observed_between_group", {
  df <- mk_cov_between_data()
  spec <- compile_model(y ~ x + (1 | g), df)
  observed <- roles(spec)
  expect_s3_class(observed, "mm_roles")
  x_row <- observed$table[observed$table$variable == "x", , drop = FALSE]
  expect_equal(nrow(x_row), 1L)
  expect_identical(x_row$role[[1L]], "observed_between_group")
  expect_match(x_row$evidence[[1L]], "constant within")
  out <- capture.output(print(observed))
  expect_true(is.character(out) && length(out) > 0L)
})

test_that("ranef(condVar=TRUE) on uncorrelated slope merges block-diag postVar", {
  fit <- mk_cov_lmm_uncorrelated()
  re <- ranef(fit, condVar = TRUE)
  expect_true("g" %in% names(re))
  pv <- attr(re[["g"]], "postVar")
  expect_true(is.array(pv))
  expect_identical(length(dim(pv)), 3L)
  # Two RE terms on same group -> 2 x 2 conditional variance per level.
  expect_equal(dim(pv)[[1L]], 2L)
  expect_equal(dim(pv)[[2L]], 2L)
  expect_equal(dim(pv)[[3L]], nrow(re[["g"]]))
  expect_true(all(is.finite(pv)))
  # Independent RE terms: off-diagonal of the merged block-diag should be 0.
  expect_equal(as.numeric(pv[1L, 2L, ]), rep(0, dim(pv)[[3L]]))
})

test_that("simulate on correlated/uncorrelated slope model prints", {
  fit <- mk_cov_lmm_uncorrelated()
  sims <- simulate(fit, nsim = 2L, seed = 1L)
  expect_s3_class(sims, "data.frame")
  expect_equal(ncol(sims), 2L)
  out <- capture.output(print(sims))
  expect_true(is.character(out) && length(out) > 0L)
})

test_that("mm_json_parse_* guards empty and invalid JSON", {
  expect_error(
    mixeff:::mm_json_parse_artifact(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_artifact(NULL),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_artifact("{not-json"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_audit_report(""),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_audit_report("[]"),
    class = "mm_schema_error"
  )
  expect_error(
    mixeff:::mm_json_parse_fixed_effect_inference_table(1L),
    class = "mm_schema_error"
  )
  expect_null(mixeff:::mm_json_parse_fixed_effect_inference_table(NULL))
})

test_that("print VarCorr/ranef/coef also work on GLMM", {
  fit <- mk_cov_glmm()
  expect_true(length(capture.output(print(VarCorr(fit)))) > 0L)
  expect_true(length(capture.output(print(ranef(fit)))) > 0L)
  expect_true(length(capture.output(print(coef(fit)))) > 0L)
})

test_that("reporting_table section=all prints all sections", {
  fit <- mk_cov_lmm()
  rt <- reporting_table(fit, section = "all")
  out <- capture.output(print(rt))
  expect_true(is.character(out) && length(out) > 0L)
  expect_match(paste(out, collapse = "\n"), "overview")
})

test_that("refit and simulate arg guards on LMM/GLMM", {
  fit <- mk_cov_lmm()
  expect_error(refit(fit, 1:3), class = "mm_arg_error")
  expect_error(simulate(fit, nsim = 0), class = "mm_arg_error")

  gfit <- mk_cov_glmm()
  expect_error(refit(gfit, rbinom(nobs(gfit), 1, 0.5)),
               class = "mm_inference_unavailable")
  expect_error(simulate(gfit, nsim = 1), class = "mm_inference_unavailable")
})
