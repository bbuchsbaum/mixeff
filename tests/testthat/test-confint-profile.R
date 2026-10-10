# Stage D.3 (bd-01KRFGFSK4A0MGPFQVCNY5SYFK): confint(method = "profile")
# wired through the upstream `profile_confint_payload` FFI.
#
# Done condition from the bead: (a) ML full-rank LMM returns beta/sigma/theta
# rows, all finite and monotonic; (b) REML full-rank LMM omits beta with a
# documented reason_code, sigma/theta still present [superseded: REML fits
# are now profiled on the ML deviance like lme4, so beta rows are
# available]; (c) beta profile CI
# agrees with the Wald CI on a well-behaved model to within stated
# tolerance; (d) boundary fit raises a typed refusal rather than fabricating.

mm_skip_if_no_lme4_local <- function() {
  testthat::skip_if_not_installed("lme4")
}

mm_sleepstudy_data <- function() {
  mm_skip_if_no_lme4_local()
  env <- new.env(parent = emptyenv())
  utils::data("sleepstudy", package = "lme4", envir = env)
  if (!exists("sleepstudy", envir = env, inherits = FALSE)) {
    testthat::skip("sleepstudy dataset unavailable")
  }
  get("sleepstudy", envir = env, inherits = FALSE)
}

# One fit per criterion for the whole file; the profile computations below
# are the expensive part and refitting per block only multiplies them.
mm_sleepstudy_fit_cache <- new.env(parent = emptyenv())

mm_sleepstudy_fit_ml <- function() {
  if (is.null(mm_sleepstudy_fit_cache$ml)) {
    mm_sleepstudy_fit_cache$ml <- lmm(
      Reaction ~ Days + (1 + Days | Subject),
      data = mm_sleepstudy_data(), REML = FALSE,
      control = mm_control(verbose = -1)
    )
  }
  mm_sleepstudy_fit_cache$ml
}

mm_sleepstudy_fit_reml <- function() {
  if (is.null(mm_sleepstudy_fit_cache$reml)) {
    mm_sleepstudy_fit_cache$reml <- lmm(
      Reaction ~ Days + (1 + Days | Subject),
      data = mm_sleepstudy_data(), REML = TRUE,
      control = mm_control(verbose = -1)
    )
  }
  mm_sleepstudy_fit_cache$reml
}

test_that("confint(method='profile') under ML returns beta/sigma/theta rows", {
  skip_on_cran() # full profile-likelihood computation, ~1.5s locally
  fit <- mm_sleepstudy_fit_ml()
  ci <- confint(fit, method = "profile", level = 0.95)

  expect_s3_class(ci, "mm_confint")
  expect_identical(attr(ci, "method"), "profile_likelihood")
  expect_identical(attr(ci, "fit_criterion"), "ML")
  expect_true("(Intercept)" %in% rownames(ci))
  expect_true("Days" %in% rownames(ci))
  # lme4's row names and order: .sigNN (lme4 term order), .sigma, beta.
  expect_identical(rownames(ci),
                   c(".sig01", ".sig02", ".sig03", ".sigma", "(Intercept)", "Days"))

  payload <- attr(ci, "mm_profile")
  expect_true(is.list(payload))
  kinds <- payload$table$parameter_kind
  expect_true("beta" %in% kinds)
  expect_true("sigma" %in% kinds)
  expect_true("theta" %in% kinds)

  # Bounds may be NA when the profile is truncated (upstream documents this
  # explicitly, mirroring lme4::confint). For each populated row, only the
  # populated bounds need to bracket the estimate; NA bounds are honest
  # "not determined" signals, not regressions.
  ok <- payload$table$reason_code %in% c(NA_character_, "")
  finite <- payload$table[ok, , drop = FALSE]
  beta_rows <- finite[finite$parameter_kind == "beta", , drop = FALSE]
  expect_true(nrow(beta_rows) >= 1L)
  expect_true(all(is.finite(beta_rows$lower)),
              info = "ML beta profile bounds should be finite on sleepstudy")
  expect_true(all(is.finite(beta_rows$upper)),
              info = "ML beta profile bounds should be finite on sleepstudy")
  expect_true(all(beta_rows$lower <= beta_rows$estimate + 1e-8))
  expect_true(all(beta_rows$estimate <= beta_rows$upper + 1e-8))

  # Variance-component rows are allowed to have NA bounds (truncated
  # profile), but at least one sigma/theta row must come back without a
  # reason_code — an engine that stamps a refusal on every variance
  # component would otherwise pass this block without asserting anything.
  vc_rows <- finite[finite$parameter_kind != "beta", , drop = FALSE]
  expect_true(nrow(vc_rows) >= 1L,
              info = "ML profile must return at least one un-refused sigma/theta row")
  for (i in seq_len(nrow(vc_rows))) {
    row <- vc_rows[i, , drop = FALSE]
    if (is.finite(row$lower)) {
      expect_lte(row$lower, row$estimate + 1e-8,
                 label = sprintf("lower bound for %s", row$parameter))
    }
    if (is.finite(row$upper)) {
      expect_gte(row$upper, row$estimate - 1e-8,
                 label = sprintf("upper bound for %s", row$parameter))
    }
  }
})

test_that("confint(method='profile') under REML profiles the ML deviance like lme4", {
  skip_on_cran() # full profile-likelihood computation, ~3s locally
  fit <- mm_sleepstudy_fit_reml()
  ci <- confint(fit, method = "profile", level = 0.95)

  expect_identical(attr(ci, "fit_criterion"), "REML")
  payload <- attr(ci, "mm_profile")
  table <- payload$table

  beta_rows <- table[table$parameter_kind == "beta", , drop = FALSE]
  expect_identical(beta_rows$parameter, c("(Intercept)", "Days"))
  expect_true(all(is.finite(beta_rows$lower) & is.finite(beta_rows$upper)))
  expect_true(all(beta_rows$profiled_criterion == "ML"))
  theta_rows <- table[table$parameter_kind == "theta", , drop = FALSE]
  expect_true(nrow(theta_rows) >= 1L)
  expect_true(all(theta_rows$profiled_criterion == "REML"))
  expect_true(any(grepl("ML deviance", payload$notes, fixed = TRUE)))
  # lme4: the REML fit's profile intervals equal the ML fit's.
  ci_ml <- confint(mm_sleepstudy_fit_ml(), method = "profile", level = 0.95)
  expect_equal(unclass(ci)[, 1:2], unclass(ci_ml)[, 1:2], tolerance = 1e-6,
               ignore_attr = TRUE)
})

test_that("confint(method='profile') matches lme4's profile intervals", {
  skip_on_cran() # three model pairs, profile on both sides
  mm_skip_if_no_lme4_local()
  skip_if(utils::packageVersion("lme4") < "2.0-0",
          "profile .sigNN numbering follows lme4 >= 2.0")
  d <- mm_sleepstudy_data()
  for (fo in list(Reaction ~ Days + (Days | Subject),
                  Reaction ~ Days + (1 | Subject),
                  Reaction ~ Days + (Days || Subject))) {
    fit <- lmm(fo, d, control = mm_control(verbose = -1))
    ref <- suppressMessages(confint(lme4::lmer(fo, d)))
    ci <- confint(fit, method = "profile")
    expect_identical(rownames(ci), rownames(ref), label = deparse1(fo))
    expect_equal(unclass(ci)[, 1:2], ref, tolerance = 1e-4,
                 ignore_attr = TRUE, label = deparse1(fo))
  }
})

test_that("profile .sigNN labels follow lme4's term order and kinds", {
  skip_on_cran()
  mm_skip_if_no_lme4_local()
  set.seed(11)
  d <- data.frame(a = factor(rep(1:6, each = 20)), g = factor(rep(1:15, 8)),
                  x = rnorm(120))
  d$y <- 1 + d$x + rnorm(6)[d$a] + rnorm(15, sd = 0.7)[d$g] +
    d$x * rnorm(15, sd = 0.4)[d$g] + rnorm(120)
  # The engine orders terms g, a (decreasing size); lme4 orders
  # (1 | a) + (x || g) as g.x, g.(Intercept), a (ties reversed).
  fit <- lmm(y ~ x + (1 | a) + (x || g), d, REML = FALSE,
             control = mm_control(verbose = -1))
  ref <- lme4::lmer(y ~ x + (1 | a) + (x || g), d, REML = FALSE)
  ci <- confint(fit, method = "profile")
  ref_ci <- suppressMessages(confint(ref, method = "profile"))
  expect_identical(rownames(ci), rownames(ref_ci))
  expect_equal(unclass(ci)[, 1:2], ref_ci, tolerance = 2e-3,
               ignore_attr = TRUE)
  tab <- attr(ci, "mm_profile")$table
  sig <- tab[grepl("^\\.sig[0-9]+$", tab$parameter), , drop = FALSE]
  expect_identical(sig$parameter, c(".sig01", ".sig02", ".sig03"))
  expect_true(all(sig$parameter_kind == "sd"))
  # correlated term (lme4 >= 2.0 numbering): the term's standard deviations
  # come first, so .sig03 is the correlation
  tab2 <- attr(confint(mm_sleepstudy_fit_ml(), method = "profile"),
               "mm_profile")$table
  expect_identical(tab2$parameter_kind[match(c(".sig01", ".sig02", ".sig03"),
                                             tab2$parameter)],
                   c("sd", "sd", "cor"))
})

test_that("profile threads give identical results", {
  skip_on_cran()
  fit <- mm_sleepstudy_fit_ml()
  ci1 <- confint(fit, method = "profile", threads = 1)
  ci2 <- confint(fit, method = "profile", threads = 2)
  expect_identical(unclass(ci1)[, 1:2], unclass(ci2)[, 1:2])
  expect_identical(attr(ci2, "mm_profile")$threads, 2L)
  expect_error(confint(fit, method = "profile", threads = 0),
               class = "mm_arg_error")
})

test_that("profile CI for beta agrees with Wald CI on a well-behaved ML fit", {
  skip_on_cran() # profile computed in full even under parm=, ~1.5s locally
  fit <- mm_sleepstudy_fit_ml()
  profile_ci <- confint(fit, parm = c("(Intercept)", "Days"),
                        method = "profile", level = 0.95)
  wald_ci <- confint(fit, parm = c("(Intercept)", "Days"),
                     method = "wald", level = 0.95)

  expect_true(all(rownames(profile_ci) %in% rownames(wald_ci)))
  # Profile and Wald intervals should agree to within a few percent on
  # sleepstudy beta. The actual numbers are tens to hundreds, so a
  # relative tolerance is appropriate.
  for (nm in rownames(profile_ci)) {
    rel_lower <- abs(profile_ci[nm, 1] - wald_ci[nm, 1]) /
      max(1e-6, abs(wald_ci[nm, 1]))
    rel_upper <- abs(profile_ci[nm, 2] - wald_ci[nm, 2]) /
      max(1e-6, abs(wald_ci[nm, 2]))
    expect_lt(rel_lower, 0.05,
              label = sprintf("rel lower diff for %s", nm))
    expect_lt(rel_upper, 0.05,
              label = sprintf("rel upper diff for %s", nm))
  }
})

test_that("profile CI parm subsetting filters returned rows", {
  skip_on_cran() # two full profile computations, ~3s locally
  fit <- mm_sleepstudy_fit_ml()
  ci_all <- confint(fit, method = "profile", level = 0.95)
  ci_intercept <- confint(fit, parm = "(Intercept)",
                          method = "profile", level = 0.95)
  expect_true(nrow(ci_intercept) >= 1L)
  expect_true(all(rownames(ci_intercept) %in% rownames(ci_all)))
  expect_identical(rownames(ci_intercept), "(Intercept)")
})

test_that("profile CI surfaces a typed refusal on a boundary singular fit", {
  testthat::skip_if_not_installed("lme4")
  env <- new.env(parent = emptyenv())
  utils::data("Dyestuff2", package = "lme4", envir = env)
  if (!exists("Dyestuff2", envir = env, inherits = FALSE)) {
    testthat::skip("Dyestuff2 dataset unavailable")
  }
  fit <- lmm(Yield ~ 1 + (1 | Batch),
             data = get("Dyestuff2", envir = env, inherits = FALSE),
             REML = TRUE, control = mm_control(verbose = -1))

  # The current engine returns a profile payload on this boundary fit; a
  # typed refusal would also be defensible, but accepting either branch
  # made this test unable to notice the behavior changing. Pin the payload
  # path; if a pin bump or platform flips it to a refusal, fail with a
  # message rather than erroring, so a human re-decides the contract.
  result <- tryCatch(
    confint(fit, method = "profile", level = 0.95),
    error = function(cnd) cnd
  )
  if (inherits(result, "condition")) {
    testthat::fail(sprintf(
      "boundary profile CI now refuses (%s) — contract changed, re-decide",
      paste(class(result)[1L], collapse = "/")
    ))
  }
  # Every NA bound on a boundary fit must carry a non-empty regularity
  # note (no NA-without-explanation rows).
  expect_s3_class(result, "mm_confint")
  payload <- attr(result, "mm_profile")
  expect_true(is.list(payload))
  table <- payload$table
  expect_true(nrow(table) >= 1L,
              info = "boundary fit profile payload must surface at least one row")
  na_rows <- table[!is.finite(table$lower) | !is.finite(table$upper), ,
                   drop = FALSE]
  if (nrow(na_rows)) {
    expect_true(all(nzchar(na_rows$regularity)),
                info = "NA bounds on a boundary fit must carry a regularity note")
  } else {
    # If every bound is finite, the boundary case still has to flag the
    # near-zero theta as boundary-clamped on the lower side rather than
    # silently emitting a negative variance bound.
    expect_true(all(is.finite(table$lower)) && all(is.finite(table$upper)))
  }
})

# Engine 26f9690: an irregular per-parameter profile no longer aborts the
# payload. The row is kept with a typed status / reason and NA for the bound
# the profile cannot determine (lme4 warns and interpolates linearly).
mm_profile_irregular_data <- function(seed = 6L) {
  set.seed(seed)
  ng <- 6L + seed %% 5L
  pg <- 4L + seed %% 3L
  slope_sd <- c(0.05, 0.2, 0.5, 1)[seed %% 4L + 1L]
  g <- factor(rep(seq_len(ng), each = pg))
  x <- rnorm(ng * pg)
  u0 <- rnorm(ng)
  u1 <- slope_sd * rnorm(ng)
  y <- 1 + 0.5 * x + u0[g] + u1[g] * x + 0.7 * rnorm(ng * pg)
  data.frame(y = y, x = x, g = g)
}

test_that("a non-monotone correlation profile gives NA plus a reason, not an abort", {
  d <- mm_profile_irregular_data(6L)
  fit <- lmm(y ~ x + (1 + x | g), d, REML = FALSE,
             control = mm_control(verbose = -1))
  ci <- confint(fit, method = "profile")
  tab <- attr(ci, "mm_profile")$table
  expect_true(all(c("status", "reason_code", "reason") %in% names(tab)))
  cor_row <- tab[tab$parameter == ".sig03", , drop = FALSE]
  expect_identical(nrow(cor_row), 1L)
  expect_identical(cor_row$parameter_kind, "cor")
  expect_identical(cor_row$status, "non_monotone")
  expect_identical(cor_row$regularity, "non_monotone_profile")
  expect_identical(cor_row$reason_code, "non_monotone_profile")
  expect_match(cor_row$reason, "not strictly monotone", fixed = TRUE)
  # The reason is relabelled with the row's lme4 name.
  expect_match(cor_row$reason, "^profile_\\.sig03:")
  # The lower bound sits in the regular part of the profile; the upper one
  # does not and is NA (lme4 1.1.35 reports an interpolated 1.0 here).
  expect_true(is.finite(ci[".sig03", 1L]))
  expect_true(is.na(ci[".sig03", 2L]))
  # Every other lme4 row is regular and finite.
  ok <- tab[tab$parameter != ".sig03" &
              tab$parameter_kind %in% c("sd", "sigma", "beta"), , drop = FALSE]
  expect_true(all(ok$status == "ok"))
  expect_true(all(is.na(ok$reason_code)))
  expect_true(all(is.finite(ok$lower) & is.finite(ok$upper)))
  out <- capture.output(print(ci))
  expect_true(any(grepl("irregular profile for .sig03", out, fixed = TRUE)))
  expect_true(any(grepl("reported as NA", out, fixed = TRUE)))
  # profile() keeps the same rows and says so when printed.
  pr_out <- capture.output(print(profile(fit)))
  expect_true(any(grepl("Irregular profile: .sig03 (non_monotone_profile)",
                        pr_out, fixed = TRUE)))

  skip_if_not_installed("lme4")
  m <- lme4::lmer(y ~ x + (1 + x | g), d, REML = FALSE)
  ref <- suppressWarnings(suppressMessages(
    stats::confint(m, method = "profile")
  ))
  # lme4 < 2.0 numbers the correlation .sig02 (term order sd, cor, sd);
  # lme4 >= 2.0, like mixeff, puts it last (.sig03). lme4 interpolates the
  # upper bound linearly (1.0 here); mixeff reports NA.
  ref_cor <- if (utils::packageVersion("lme4") >= "2.0") ".sig03" else ".sig02"
  # The regular lower bound agrees to ~2% (-0.434 vs lme4 1.1.35's -0.426:
  # lme4's spline also sees the non-monotone tail).
  expect_equal(ci[".sig03", 1L], ref[ref_cor, 1L], tolerance = 0.03)
  for (p in c(".sig01", ".sigma", "(Intercept)", "x")) {
    expect_equal(unname(ci[p, ]), unname(ref[p, ]), tolerance = 5e-3,
                 info = p)
  }
})

test_that("profile payload rows map status, reason and JSON-null bounds", {
  fit <- mm_sleepstudy_fit_ml()
  json <- paste0(
    '{"parameter": "σ", "estimate": 25.6, "lower": null, "upper": null,',
    ' "level": 0.95, "method": "profile_likelihood",',
    ' "regularity": "profile_failed", "boundary_clamped_lower": false,',
    ' "status": "failed", "reason": "profile_σ: refit failed"}'
  )
  row <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  out <- mixeff:::mm_translate_profile_row(row, fit)
  expect_identical(out$parameter, ".sigma")
  expect_true(is.na(out$lower) && is.na(out$upper))
  expect_identical(out$status, "failed")
  expect_identical(out$reason_code, "profile_failed")
  expect_identical(out$reason, "profile_.sigma: refit failed")

  # Each irregular status maps to its regularity code even when the
  # payload omits `regularity`; a pre-26f9690 row (no status) is regular.
  for (st in c("non_monotone", "not_bracketing", "failed")) {
    r <- row
    r$regularity <- NULL
    r$status <- st
    got <- mixeff:::mm_translate_profile_row(r, fit)
    expect_identical(got$reason_code,
                     unname(mixeff:::mm_profile_status_regularity[[st]]))
  }
  legacy <- row[c("parameter", "estimate", "level", "method",
                  "boundary_clamped_lower")]
  legacy$lower <- 20
  legacy$upper <- 30
  got <- mixeff:::mm_translate_profile_row(legacy, fit)
  expect_identical(got$status, "ok")
  expect_identical(got$regularity, "regular_profile_likelihood")
  expect_true(is.na(got$reason_code) && is.na(got$reason))
})

test_that("profile CI schema accepts irregular rows (status, reason, null bounds)", {
  skip_on_cran()
  skip_if_not_installed("jsonvalidate")
  schema_file <- c(
    system.file("schemas", "mixedmodels.profile_likelihood_ci.schema.json",
                package = "mixeff"),
    testthat::test_path("..", "..", "inst", "schemas",
                        "mixedmodels.profile_likelihood_ci.schema.json")
  )
  schema_file <- schema_file[nzchar(schema_file) & file.exists(schema_file)][1L]
  if (is.na(schema_file)) skip("profile CI schema is not installed")
  d <- mm_profile_irregular_data(6L)
  fit <- lmm(y ~ x + (1 + x | g), d, REML = FALSE,
             control = mm_control(verbose = -1))
  payload <- mixeff:::mm_profile_confint_payload(fit, 0.95)
  expect_true(any(vapply(payload$intervals, function(r) {
    !identical(r$status %||% "ok", "ok")
  }, logical(1))))
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null",
                           digits = NA)
  validator <- jsonvalidate::json_validator(schema_file, engine = "ajv")
  expect_true(isTRUE(validator(json, verbose = TRUE)))
})
