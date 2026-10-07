# Consumer acceptance for engine-owned fitted-state persistence.
snapshot_fixture <- function(glmm_method = NULL, offset = FALSE) {
  set.seed(813)
  d <- data.frame(g = factor(rep(seq_len(12), each = 10)), x = rnorm(120))
  eta <- 0.3 + 0.4 * d$x + rep(rnorm(12, sd = 0.7), each = 10)
  d$y <- if (is.null(glmm_method)) eta + rnorm(120) else rpois(120, exp(eta))
  if (is.null(glmm_method)) {
    lmm(y ~ x + (1 | g), d, weights = seq(0.7, 1.3, length.out = 120),
        control = mm_control(verbose = -1))
  } else {
    glmm(y ~ x + (1 | g), d, family = poisson(), method = glmm_method,
         offset = if (offset) rep(0.1, 120) else NULL, control = mm_control(verbose = -1))
  }
}

snapshot_lmm_queries <- function(fit, nd = fit$model_frame[1:5, ]) {
  list(beta = fixef(fit), theta = fit$theta, objective = fit$fit$optimizer$objective,
       certificate = optimizer_certificate(fit),
       contrast = contrast(fit, c(0, 1), method = "satterthwaite")$table,
       terms = anova(fit, method = "satterthwaite"),
       conditional = ranef(fit, condVar = TRUE),
       prediction = predict(fit, newdata = nd),
       uncertainty = predict(fit, newdata = nd, se.fit = TRUE))
}

test_that("native queries restore fitted state without reconstructing a new fit", {
  fit <- snapshot_fixture()
  expect_type(fit$fit$fitted_state, "character")
  expect_null(fit$fit$restoration_error)
  nd <- fit$model_frame[1:5, ]
  expected <- snapshot_lmm_queries(fit, nd)
  # These inputs would make the old reconstruction-and-fit path fail.
  # They are not part of the engine-owned fitted solution queried here.
  restored <- revive(fit)
  restored$model_frame$y <- 0
  restored$control$optimizer <- "not_an_optimizer"
  expect_equal(snapshot_lmm_queries(restored, nd), expected, tolerance = 1e-10)
  expect_identical(restored$fit$fitted_state, fit$fit$fitted_state)
})

test_that("legacy and incompatible objects keep stored extraction and refuse native work", {
  fit <- snapshot_fixture()
  for (kind in c("missing", "version", "corrupt")) {
    old <- revive(fit)
    if (kind == "missing") {
      old$fit$fitted_state <- NULL
    } else if (kind == "version") {
      envelope <- jsonlite::fromJSON(old$fit$fitted_state, simplifyVector = FALSE)
      envelope$crate_version <- "incompatible-engine"
      old$fit$fitted_state <- as.character(jsonlite::toJSON(envelope, auto_unbox = TRUE))
    } else {
      old$fit$fitted_state <- "invalid JSON"
    }
    expect_equal(fixef(old), fixef(fit))
    expect_equal(logLik(old), logLik(fit))
    expect_equal(ranef(old), ranef(fit))
    expect_error(contrast(old, c(0, 1), method = "satterthwaite"),
                 class = "mm_inference_unavailable")
    conditional <- ranef(old, condVar = TRUE)
    expect_true(all(is.na(attr(conditional[[1]], "postVar"))))
    expect_identical(attr(conditional, "mm_unavailable_reason"),
                     "random_effect_conditional_variance_unavailable")
    expect_match(attr(conditional[[1]], "mm_cond_var_error"), "restoration")
    expect_error(predict(old, newdata = fit$model_frame[1:5, ], se.fit = TRUE),
                 class = "mm_inference_unavailable")
  }
})

test_that("verification uses the recorded LMM and GLMM baseline", {
  for (fit in list(snapshot_fixture(), snapshot_fixture("pirls_profiled", offset = TRUE),
                   snapshot_fixture("joint_laplace", offset = TRUE))) {
    original_objective <- fit$fit$optimizer$objective
    fit$control$optimizer <- "not_an_optimizer"
    fit$model_frame$y <- 0
    v <- verify_convergence(fit, restart = FALSE, jitter_starts = 0L)
    expect_equal(v$reference$objective, original_objective, tolerance = 1e-12)
    expect_equal(unname(v$reference$beta), unname(fixef(fit)), tolerance = 1e-10)
    expect_equal(v$reference$theta, fit$theta, tolerance = 1e-12)
    expect_equal(nrow(v$table), 0L)
  }
})

test_that("RDS queries survive a fresh process with cold and warm caches", {
  testthat::skip_on_cran()
  fit <- snapshot_fixture()
  fast <- snapshot_fixture("pirls_profiled")
  joint <- snapshot_fixture("joint_laplace")
  nd <- fit$model_frame[1:5, ]
  query <- function(lmm_fit, fast, joint, nd) {
    list(lmm = snapshot_lmm_queries(lmm_fit, nd),
         fast = predict(fast, newdata = nd, se.fit = TRUE),
         joint = predict(joint, newdata = nd, se.fit = TRUE))
  }
  expected <- query(fit, fast, joint, nd) # populate any query caches before saving
  input <- tempfile(fileext = ".rds")
  output <- tempfile(fileext = ".rds")
  script <- tempfile(fileext = ".R")
  saveRDS(list(fit = fit, fast = fast, joint = joint, nd = nd), input)
  # Serialized formula environments can load a namespace during readRDS().
  # Select the child installation before reading any fitted objects.
  child_package <- find.package("mixeff", lib.loc = .libPaths())
  query_lines <- deparse(snapshot_lmm_queries)
  writeLines(c(
    "args <- commandArgs(TRUE)",
    ".libPaths(args[-seq_len(3L)])",
    "library(mixeff)",
    "stopifnot(identical(normalizePath(getNamespaceInfo('mixeff', 'path')),",
    "                    normalizePath(args[[3]])))",
    "saved <- readRDS(args[[1]])",
    paste0("snapshot_lmm_queries <- ", query_lines[1]), query_lines[-1],
    "results <- lapply(c(FALSE, TRUE), function(cold) {",
    "  f <- saved$fit; fast <- saved$fast; joint <- saved$joint",
    "  if (cold) { f <- revive(f); fast <- revive(fast); joint <- revive(joint) }",
    "  list(lmm = snapshot_lmm_queries(f, saved$nd),",
    "       fast = predict(fast, newdata = saved$nd, se.fit = TRUE),",
    "       joint = predict(joint, newdata = saved$nd, se.fit = TRUE))",
    "})", "saveRDS(results, args[[2]])"
  ), script)
  on.exit(unlink(c(input, output, script)), add = TRUE)
  out <- system2(file.path(R.home("bin"), "Rscript"),
                 c("--vanilla", shQuote(script), shQuote(input), shQuote(output),
                   shQuote(child_package), shQuote(.libPaths())),
                 stdout = TRUE, stderr = TRUE)
  expect_equal(attr(out, "status") %||% 0L, 0L, info = paste(out, collapse = "\n"))
  if (file.exists(output)) {
    actual <- readRDS(output)
    expect_equal(actual[[1]], expected, tolerance = 1e-9)
    expect_equal(actual[[2]], expected, tolerance = 1e-9)
  }
})
