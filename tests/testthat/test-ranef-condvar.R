# Stage C.2 (bd-01KRCKCZMYRFVEH39RFR72ZW3F): ranef(condVar=TRUE) returns
# block-diagonal conditional variance arrays from Rust cond_var().
#
# Done condition from the bead: attr(ranef(fit, condVar=TRUE)[[group]],
# "postVar") is a non-NA c(p, p, n) array whose values agree with lme4
# within tolerance.

mm_skip_if_no_lme4_local <- function() {
  testthat::skip_if_not_installed("lme4")
}

mm_sleepstudy_fit <- function() {
  mm_skip_if_no_lme4_local()
  env <- new.env(parent = emptyenv())
  utils::data("sleepstudy", package = "lme4", envir = env)
  if (!exists("sleepstudy", envir = env, inherits = FALSE)) {
    testthat::skip("sleepstudy dataset unavailable")
  }
  data <- get("sleepstudy", envir = env, inherits = FALSE)
  fit <- lmm(
    Reaction ~ Days + (1 + Days | Subject),
    data = data,
    REML = TRUE,
    control = mm_control(verbose = -1)
  )
  ref <- suppressMessages(suppressWarnings(
    lme4::lmer(Reaction ~ Days + (1 + Days | Subject),
               data = data, REML = TRUE)
  ))
  list(fit = fit, ref = ref, data = data)
}

test_that("ranef(condVar=TRUE) returns finite p x p x n postVar arrays", {
  pair <- mm_sleepstudy_fit()
  re <- ranef(pair$fit, condVar = TRUE)
  expect_true("Subject" %in% names(re))

  pv <- attr(re[["Subject"]], "postVar")
  expect_true(is.array(pv))
  expect_identical(length(dim(pv)), 3L)
  expect_equal(dim(pv)[[1L]], dim(pv)[[2L]])
  expect_equal(dim(pv)[[3L]], nrow(re[["Subject"]]))
  expect_true(all(is.finite(pv)))

  # Symmetry on every level slice and PSD diagonals.
  for (i in seq_len(dim(pv)[[3L]])) {
    slice <- pv[, , i]
    expect_equal(slice, t(slice), tolerance = 1e-9,
                 info = sprintf("postVar slice %d is not symmetric", i))
    expect_true(all(diag(slice) >= 0),
                info = sprintf("postVar slice %d has negative diagonal", i))
  }

  # Unavailable-reason attributes must NOT be set on success.
  expect_null(attr(re[["Subject"]], "mm_unavailable_reason"))
  expect_null(attr(re, "mm_unavailable_reason"))
})

test_that("ranef(condVar=TRUE) postVar agrees with lme4 within 1e-3", {
  pair <- mm_sleepstudy_fit()
  observed_pv <- attr(ranef(pair$fit, condVar = TRUE)[["Subject"]], "postVar")
  expected_pv <- attr(lme4::ranef(pair$ref, condVar = TRUE)[["Subject"]],
                      "postVar")
  expect_identical(dim(observed_pv), dim(expected_pv))
  expect_equal(observed_pv, expected_pv, tolerance = 1e-3,
               ignore_attr = TRUE)
})

test_that("ranef(condVar=TRUE) caches across repeated calls", {
  pair <- mm_sleepstudy_fit()
  first <- attr(ranef(pair$fit, condVar = TRUE)[["Subject"]], "postVar")
  # If the cache is wired, a second call should not refit — and must return
  # the identical array.
  second <- attr(ranef(pair$fit, condVar = TRUE)[["Subject"]], "postVar")
  expect_identical(first, second)
  expect_true(exists("cond_var", envir = pair$fit$lazy_cache,
                     inherits = FALSE))
})

test_that("ranef(condVar=FALSE) is the default and returns no postVar", {
  pair <- mm_sleepstudy_fit()
  re <- ranef(pair$fit)
  expect_null(attr(re[["Subject"]], "postVar"))
  expect_null(attr(re[["Subject"]], "mm_unavailable_reason"))
})
