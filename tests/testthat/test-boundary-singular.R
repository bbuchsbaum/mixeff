# Tests for PRD §9.4 boundary tagging and §9.5.6 singular rendering.
#
# Boundary covariance components are reported model state, not a warning:
# print(VarCorr(fit)) carries a [boundary] flag on near-zero std_devs and
# explains the marker, and print(fit) appends a singular-fit summary that
# names the effective rank and points to changes() / random_options() —
# never to a folk fix.

mk_boundary_fit <- function() {
  testthat::skip_if_not_installed("lme4")
  env <- new.env(parent = emptyenv())
  utils::data("Dyestuff2", package = "lme4", envir = env)
  if (!exists("Dyestuff2", envir = env, inherits = FALSE)) {
    testthat::skip("Dyestuff2 dataset is unavailable")
  }
  lmm(
    Yield ~ 1 + (1 | Batch),
    data = get("Dyestuff2", envir = env),
    REML = TRUE,
    control = mm_control(verbose = -1)
  )
}

test_that("VarCorr table carries a boundary flag on near-zero std_devs", {
  fit <- mk_boundary_fit()
  vc <- VarCorr(fit)
  expect_s3_class(vc, "mm_varcorr")
  expect_true("boundary" %in% names(vc$table))
  expect_true(any(vc$table$boundary),
              info = "Dyestuff2 should produce at least one boundary component")
})

test_that("snapshot: print(VarCorr(fit)) tags boundary components with [boundary]", {
  fit <- mk_boundary_fit()
  expect_snapshot(print(VarCorr(fit)))
})

test_that("snapshot: print(fit) on singular fit names rank and points to audit verbs", {
  fit <- mk_boundary_fit()
  printed <- paste(capture.output(print(fit)), collapse = "\n")

  # Hard-rule checks (machine-readable contract):
  if (isTRUE(is_singular(fit))) {
    expect_match(printed, "covariance matrix is rank-deficient",
                 fixed = TRUE)
    expect_match(printed, "Use changes(fit)", fixed = TRUE)
    expect_match(printed, "Use random_options(spec, group", fixed = TRUE)
    # R9 contract: never advise the user to change their model.
    forbidden_literal <- c(
      "Try (1 | ", "Drop the random slope", "suggested starting model",
      "we recommend", "you should"
    )
    for (pattern in forbidden_literal) {
      expect_false(grepl(pattern, printed, fixed = TRUE),
                   info = sprintf("forbidden phrase matched: %s", pattern))
    }
    expect_false(grepl("try .* instead", printed, perl = TRUE),
                 info = "forbidden phrase matched: try .* instead")
    expect_snapshot(cat(printed))
  } else {
    testthat::skip("Dyestuff2 did not converge to a boundary on this build")
  }
})
