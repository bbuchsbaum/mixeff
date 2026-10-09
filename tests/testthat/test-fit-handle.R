# Native fitted-model handles (pre-CRAN checklist §5.2(c), §5.9).
#
# lmm()/glmm() keep the fitted engine model behind `fit$rust_handle`;
# follow-on computations reuse it instead of refitting cold. The handle does
# not survive serialization, and the cold path must then give IDENTICAL
# results. `mm_handle_hits()` counts bridge calls served from a live handle,
# so each test also proves which path ran.

mm_handle_hits <- function() mixeff:::mm_handle_hits()

# A serialize/unserialize round trip: the handle comes back as a NULL
# external pointer (dead) and the lazy cache environment is a fresh copy, so
# nothing computed on the live fit leaks into the cold one.
mm_round_trip <- function(fit) unserialize(serialize(fit, NULL))

# Results that embed fit objects (compare(), formulas whose environment holds
# the fits) differ only in the handle itself: drop external pointers and
# formula environments before comparing.
mm_strip_handles <- function(x) {
  if (identical(typeof(x), "externalptr")) return(NULL)
  if (is.environment(x)) return(NULL)
  at <- attributes(x)
  if (is.list(x) && length(x)) {
    x[] <- lapply(x, function(el) {
      v <- mm_strip_handles(el)
      if (is.null(v)) list(NULL) else list(v)
    }) |> unlist(recursive = FALSE)
  }
  for (nm in setdiff(names(at), c("names", "class", "dim", "dimnames",
                                  "row.names"))) {
    attr(x, nm) <- if (identical(nm, ".Environment")) NULL else
      mm_strip_handles(at[[nm]])
  }
  x
}

expect_same_result <- function(live, dead, info = NULL) {
  expect_identical(mm_strip_handles(live), mm_strip_handles(dead), info = info)
}

# Run `expr` and report how many live-handle hits it produced.
mm_count_hits <- function(expr) {
  before <- mm_handle_hits()
  value <- force(expr)
  list(value = value, hits = mm_handle_hits() - before)
}

mm_handle_sleep <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(101)
      n_g <- 12L
      n_per <- 8L
      g <- factor(rep(seq_len(n_g), each = n_per))
      x <- rep(seq(0, 1, length.out = n_per), n_g)
      f <- factor(rep(c("a", "b", "c", "d"), length.out = n_g * n_per))
      b0 <- rnorm(n_g, sd = 1)
      b1 <- rnorm(n_g, sd = 0.5)
      y <- 2 + 1.5 * x + 0.4 * (f == "b") + b0[g] + b1[g] * x +
        rnorm(n_g * n_per, sd = 0.6)
      cache <<- data.frame(y = y, x = x, f = f, g = g)
    }
    cache
  }
})

mm_handle_binom <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      set.seed(7)
      g <- factor(rep(1:12, each = 12))
      x <- rnorm(144)
      eta <- -0.3 + 0.8 * x + rep(rnorm(12, sd = 0.5), each = 12)
      cache <<- data.frame(y = rbinom(144, 1, plogis(eta)), x = x, g = g)
    }
    cache
  }
})

test_that("lmm() and glmm() fits carry a live handle that dies on serialization", {
  d <- mm_handle_sleep()
  fit <- lmm(y ~ x + (1 + x | g), d, control = mm_control(verbose = -1))
  expect_true(fit_handle_alive(fit))
  expect_identical(typeof(fit$rust_handle), "externalptr")
  expect_identical(attr(fit$rust_handle, "mm_handle_key")$kind, "lmm")
  expect_false("handle" %in% names(fit$fit))

  dead <- mm_round_trip(fit)
  expect_identical(typeof(dead$rust_handle), "externalptr")
  expect_false(fit_handle_alive(dead))

  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  restored <- readRDS(path)
  expect_false(fit_handle_alive(restored))
  expect_null(revive(restored)$rust_handle)
  # revive() keeps a live handle (it only clears dead ones).
  expect_true(fit_handle_alive(revive(fit)))

  # Foreign / NULL external pointers are never dereferenced.
  live <- fit
  live$rust_handle <- new("externalptr")
  expect_false(fit_handle_alive(live))

  gfit <- glmm(y ~ x + (1 | g), mm_handle_binom(), family = binomial(),
               control = mm_control(verbose = -1))
  expect_true(fit_handle_alive(gfit))
  expect_identical(attr(gfit$rust_handle, "mm_handle_key")$kind, "glmm")
  expect_false(fit_handle_alive(mm_round_trip(gfit)))
})

test_that("options(mixeff.keep_handle = FALSE) fits without a handle", {
  d <- mm_handle_sleep()
  old <- options(mixeff.keep_handle = FALSE)
  on.exit(options(old), add = TRUE)
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  expect_null(fit$rust_handle)
  expect_true("rust_handle" %in% names(fit))
  expect_false(fit_handle_alive(fit))
  options(old)
  live <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  # Same numbers either way: the handle is a cache, not part of the result.
  for (nm in c("beta", "theta", "sigma", "logLik", "fitted", "residuals",
               "std_errors", "fixed_effect_vcov", "random_effects")) {
    expect_identical(live[[nm]], fit[[nm]], info = nm)
  }
  h <- mm_count_hits(summary(fit))
  expect_identical(h$hits, 0)
  expect_identical(mm_count_hits(summary(live))$hits, 1)
})

test_that("compile_model() output is unchanged and carries no handle", {
  d <- mm_handle_sleep()
  spec <- compile_model(y ~ x + f + (1 + x | g), d)
  expect_identical(names(spec), c("call", "formula", "vars", "model_frame",
                                  "artifact"))
  expect_identical(class(spec), c("mm_spec", "mm_compiled"))
  expect_identical(spec$call, quote(compile_model(formula = y ~ x + f + (1 + x | g),
                                                  data = d)))
  # Same artifact as the low-level one-step compile primitive.
  sd <- mixeff:::mm_translate_data(spec$model_frame)
  raw <- mixeff:::mm_compile_model_json(
    "y ~ x + f + (1 + x | g)", sd$column_order, sd$numeric_columns,
    sd$categorical_values, sd$categorical_levels, sd$categorical_ordered
  )
  expect_identical(spec$artifact, mixeff:::mm_json_parse_artifact(raw))
  # The compile-for-fit route yields the same pre-fit artifact for a plain
  # formula, plus a spec handle the fit consumes.
  for_fit <- mixeff:::mm_compile_model_spec(
    "y ~ x + f + (1 + x | g)", sd$column_order, sd$numeric_columns,
    sd$categorical_values, sd$categorical_levels, sd$categorical_ordered
  )
  expect_identical(for_fit$json, raw)
  expect_identical(typeof(for_fit$handle), "externalptr")
  # A spec the engine refuses (no random effects) falls back to the one-step
  # artifact without a handle, so lmm() keeps raising its usual typed error.
  no_re <- mixeff:::mm_compile_model_spec(
    "y ~ x", sd$column_order, sd$numeric_columns,
    sd$categorical_values, sd$categorical_levels, sd$categorical_ordered
  )
  expect_null(no_re$handle)
})

test_that("lmm() fitted from the compiled spec matches the one-step fit", {
  d <- mm_handle_sleep()
  sd <- mixeff:::mm_translate_data(d)
  ctl <- as.character(jsonlite::toJSON(unclass(mm_control(verbose = -1)),
                                       auto_unbox = TRUE, null = "null",
                                       digits = NA))
  form <- "y ~ x + f + (1 + x | g)"
  one_step <- mixeff:::mm_fit_lmm_json(
    form, TRUE, sd$column_order, sd$numeric_columns, sd$categorical_values,
    sd$categorical_levels, sd$categorical_ordered, numeric(), ctl
  )
  spec <- mixeff:::mm_compile_model_spec(
    form, sd$column_order, sd$numeric_columns, sd$categorical_values,
    sd$categorical_levels, sd$categorical_ordered
  )
  e <- mixeff:::mm_empty_spec_data()
  from_spec <- mixeff:::mm_fit_lmm_json(
    form, TRUE, e$column_order, e$numeric_columns, e$categorical_values,
    e$categorical_levels, e$categorical_ordered, numeric(), ctl,
    spec$handle, TRUE
  )
  strip <- function(json) {
    x <- jsonlite::fromJSON(json, simplifyVector = FALSE)
    a <- jsonlite::fromJSON(x$artifact_json, simplifyVector = FALSE)
    x$artifact_json <- NULL
    list(x, a$optimizer_certificate$status, a$theta_map)
  }
  expect_identical(strip(from_spec$json), strip(one_step$json))
  expect_identical(from_spec$fitted, one_step$fitted)
  expect_identical(from_spec$residuals, one_step$residuals)
  expect_null(one_step$handle)
  expect_identical(typeof(from_spec$handle), "externalptr")
  # The spec is consumed: a second fit from the same handle falls back to the
  # data path (here: no data, so it fails cleanly instead of reusing it).
  expect_error(mixeff:::mm_fit_lmm_json(
    form, TRUE, e$column_order, e$numeric_columns, e$categorical_values,
    e$categorical_levels, e$categorical_ordered, numeric(), ctl,
    spec$handle, TRUE
  ))
})

test_that("LMM follow-ons on the live handle are identical to cold refits", {
  d <- mm_handle_sleep()
  fit <- lmm(y ~ x + f + (1 + x | g), d, control = mm_control(verbose = -1))
  cold <- mm_round_trip(fit)
  expect_false(fit_handle_alive(cold))

  check <- function(expr_live, expr_cold, label, min_hits = 1) {
    live <- mm_count_hits(expr_live)
    dead <- mm_count_hits(expr_cold)
    expect_gte(live$hits, min_hits)
    expect_identical(dead$hits, 0, info = label)
    expect_same_result(live$value, dead$value, info = label)
  }

  check(summary(fit), summary(cold), "summary satterthwaite")
  check(summary(fit, method = "kenward_roger"),
        summary(cold, method = "kenward_roger"), "summary KR")
  check(anova(fit), anova(cold), "anova")
  check(contrast(fit, c(0, 1, 0, 0, 0)),
        contrast(cold, c(0, 1, 0, 0, 0)), "contrast")
  check(test_effect(fit, "f", method = "kenward_roger"),
        test_effect(cold, "f", method = "kenward_roger"), "joint F")
  check(confint(fit, method = "profile"), confint(cold, method = "profile"),
        "profile")
  nd <- d[c(1, 9, 17), ]
  check(predict(fit, newdata = nd), predict(cold, newdata = nd),
        "predict newdata")
  check(predict(fit, newdata = nd, interval = "prediction"),
        predict(cold, newdata = nd, interval = "prediction"),
        "predict interval")
  check(as.data.frame(ranef(fit, condVar = TRUE)),
        as.data.frame(ranef(cold, condVar = TRUE)), "condVar")
  check(verify_convergence(fit, jitter_starts = 1L),
        verify_convergence(cold, jitter_starts = 1L), "verify_convergence")
  boot <- bootstrap_control(nsim = 4L, seed = 3L)
  check(suppressMessages(contrast(fit, c(0, 1, 0, 0, 0), method = "bootstrap",
                                  bootstrap = boot)),
        suppressMessages(contrast(cold, c(0, 1, 0, 0, 0), method = "bootstrap",
                                  bootstrap = boot)),
        "null bootstrap")
  check(suppressMessages(confint(fit, parm = 2L, method = "bootstrap",
                                 bootstrap = boot)),
        suppressMessages(confint(cold, parm = 2L, method = "bootstrap",
                                 bootstrap = boot)),
        "full-model bootstrap")

  # The live model is never mutated: after profiling and verification on the
  # handle (clones), it still reproduces the cold summary.
  expect_same_result(summary(fit), summary(cold))
})

test_that("model comparison and boundary LRT reuse the handles", {
  d <- mm_handle_sleep()
  small <- lmm(y ~ x + (1 | g), d, REML = FALSE,
               control = mm_control(verbose = -1))
  big <- lmm(y ~ x + f + (1 | g), d, REML = FALSE,
             control = mm_control(verbose = -1))
  small_c <- mm_round_trip(small)
  big_c <- mm_round_trip(big)
  live <- mm_count_hits(compare(small, big))
  dead <- mm_count_hits(compare(small_c, big_c))
  expect_gte(live$hits, 2)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)

  rs <- lmm(y ~ x + (1 + x | g), d, REML = FALSE,
            control = mm_control(verbose = -1))
  rs_c <- mm_round_trip(rs)
  live <- mm_count_hits(test_random_effect(rs, "(1 + x | g)"))
  dead <- mm_count_hits(test_random_effect(rs_c, "(1 + x | g)"))
  expect_gte(live$hits, 1)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)

  live <- mm_count_hits(parametric_bootstrap(small, big, nsim = 3L, seed = 5L))
  dead <- mm_count_hits(parametric_bootstrap(small_c, big_c, nsim = 3L,
                                             seed = 5L))
  expect_gte(live$hits, 1)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)
})

test_that("a refit under different settings never uses the handle", {
  d <- mm_handle_sleep()
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  # Profile under ML of a REML fit: key mismatch, so a cold ML refit.
  key <- attr(fit$rust_handle, "mm_handle_key")
  expect_true(key$reml)
  expect_null(mixeff:::mm_lmm_live_handle(fit, reml = FALSE))
  expect_identical(mixeff:::mm_lmm_live_handle(fit), fit$rust_handle)
  # A copy shares the handle (follow-ons only read or clone it).
  copy <- fit
  expect_identical(mixeff:::mm_lmm_live_handle(copy), fit$rust_handle)
  # Changing the stored control changes the refit, so the handle is unused.
  copy$control$max_feval <- 7L
  expect_null(mixeff:::mm_lmm_live_handle(copy))
})

test_that("GLMM follow-ons on the live handle are identical to cold refits", {
  d <- mm_handle_binom()
  fit <- glmm(y ~ x + (1 | g), d, family = binomial(),
              method = "pirls_profiled", control = mm_control(verbose = -1))
  cold <- mm_round_trip(fit)
  live <- mm_count_hits(confint(fit, method = "bootstrap", nsim = 5L,
                                seed = 2L))
  dead <- mm_count_hits(confint(cold, method = "bootstrap", nsim = 5L,
                                seed = 2L))
  expect_gte(live$hits, 1)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)

  live <- mm_count_hits(verify_convergence(fit, jitter_starts = 1L))
  dead <- mm_count_hits(verify_convergence(cold, jitter_starts = 1L))
  expect_gte(live$hits, 1)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)

  joint <- glmm(y ~ x + (1 | g), d, family = binomial(),
                method = "joint_laplace", control = mm_control(verbose = -1))
  joint_c <- mm_round_trip(joint)
  nd <- d[c(1, 13, 25), ]
  live <- mm_count_hits(predict(joint, newdata = nd, se.fit = TRUE))
  dead <- mm_count_hits(predict(joint_c, newdata = nd, se.fit = TRUE))
  expect_gte(live$hits, 1)
  expect_identical(dead$hits, 0)
  expect_same_result(live$value, dead$value)
})

test_that("a follow-on on the live handle stays interruptible", {
  skip_on_cran()
  skip_if_not_installed("processx")
  skip_if_not_installed("pkgload")
  skip_on_os("windows")
  ns_path <- getNamespaceInfo("mixeff", "path")
  load_line <- if (dir.exists(file.path(ns_path, "Meta"))) {
    sprintf(".libPaths(c(%s, .libPaths())); suppressMessages(library(mixeff))",
            deparse(dirname(ns_path)))
  } else {
    skip_if_not(file.exists(file.path(ns_path, "DESCRIPTION")),
                "package source tree is not available")
    sprintf("suppressMessages(pkgload::load_all(%s, compile = FALSE, quiet = TRUE))",
            deparse(ns_path))
  }
  script <- tempfile(fileext = ".R")
  writeLines(c(
    load_line,
    "set.seed(1)",
    "n <- 20000; d <- data.frame(s = factor(sample(400, n, TRUE)),",
    "  it = factor(sample(200, n, TRUE)), x = rnorm(n))",
    "d$y <- 1 + 0.3 * d$x + rnorm(400)[d$s] + rnorm(200)[d$it] + rnorm(n)",
    "fit <- lmm(y ~ x + (1 | s) + (1 | it), d,",
    "           control = mm_control(verbose = -1))",
    "stopifnot(fit_handle_alive(fit))",
    "before <- summary(fit)",
    "hits <- mixeff:::mm_handle_hits()",
    "cat('READY\\n')",
    # A long full-model bootstrap on the live handle (~10 ms per replicate
    # here); the engine polls the interrupt callback the fit captured.
    "r <- tryCatch({",
    "  confint(fit, parm = 2L, method = 'bootstrap',",
    "          bootstrap = bootstrap_control(nsim = 20000, seed = 1))",
    "  'COMPLETED'",
    "}, mm_interrupted = function(e) 'MM_INTERRUPTED', error = function(e) paste('ERROR', conditionMessage(e)))",
    "cat(r, '\\n')",
    "cat('USED_HANDLE', mixeff:::mm_handle_hits() > hits, '\\n')",
    # The cached model is untouched by the interrupted run.
    "cat('UNCHANGED', fit_handle_alive(fit) && identical(summary(fit), before), '\\n')"
  ), script)
  # stderr goes to a file: an undrained stderr pipe can block the child.
  errfile <- tempfile(fileext = ".txt")
  p <- processx::process$new(file.path(R.home("bin"), "Rscript"), script,
                             stdout = "|", stderr = errfile)
  on.exit(p$kill(), add = TRUE)
  out <- ""
  deadline <- Sys.time() + 180
  while (!grepl("READY", out) && p$is_alive() && Sys.time() < deadline) {
    p$poll_io(500)
    out <- paste0(out, p$read_output())
  }
  skip_if_not(grepl("READY", out), "child process did not start")
  Sys.sleep(2)
  p$interrupt()
  p$wait(120000)
  skip_if(p$is_alive(), "child process did not finish")
  out <- paste0(out, p$read_all_output())
  skip_if(grepl("COMPLETED", out), "bootstrap finished before the interrupt arrived")
  expect_match(out, "MM_INTERRUPTED", fixed = TRUE)
  expect_match(out, "USED_HANDLE TRUE", fixed = TRUE)
  expect_match(out, "UNCHANGED TRUE", fixed = TRUE)
})
