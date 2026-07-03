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

test_that("singular print advertises random_options when a slope candidate exists", {
  testthat::skip_if_not_installed("lme4")
  env <- new.env(parent = emptyenv())
  utils::data("Dyestuff2", package = "lme4", envir = env)
  d <- get("Dyestuff2", envir = env)
  # A non-intercept fixed effect gives random_options() its default slope
  # candidate, so the printed pointer is expected AND must actually run.
  d$x <- rep_len(c(-1, 0, 1), nrow(d))
  fit <- lmm(Yield ~ x + (1 | Batch), data = d, REML = TRUE,
             control = mm_control(verbose = -1))
  if (!isTRUE(is_singular(fit))) {
    testthat::skip("Dyestuff2 + x did not converge to a boundary on this build")
  }
  printed <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(printed, "Use random_options(spec, group = Batch)", fixed = TRUE)
  expect_s3_class(random_options(fit, group = Batch), "mm_random_options")
})

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
    # This fit is intercept-only with no fixed effects, so random_options()
    # has no slope candidate and would refuse -- the printed footer must not
    # advertise a verb that errors on the very fit that printed it.
    expect_false(grepl("random_options", printed, fixed = TRUE),
                 info = "random_options hint printed for a fit it cannot run on")
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
    # Scrub the crate version so engine pin bumps don't churn the snapshot.
    expect_snapshot(
      cat(printed),
      transform = function(lines) {
        sub("crate: [0-9][^ ;]*", "crate: <version>", lines)
      }
    )
  } else {
    testthat::skip("Dyestuff2 did not converge to a boundary on this build")
  }
})

test_that("print(changes(fit)) states the fitted-rank change in one sentence", {
  fit <- mk_boundary_fit()
  if (!isTRUE(is_singular(fit))) {
    testthat::skip("Dyestuff2 did not converge to a boundary on this build")
  }
  printed <- paste(capture.output(print(changes(fit))), collapse = "\n")
  expect_match(
    printed,
    "Fitted covariance for (1 | Batch): requested rank 1, fitted rank 0",
    fixed = TRUE
  )
  # The certificate sentence is the canonical record; its reduction /
  # transition restatements stay in $table but must not repeat here.
  expect_identical(
    lengths(regmatches(printed, gregexpr("rank 0", printed, fixed = TRUE))),
    1L
  )
  expect_false(grepl("formula display", printed, fixed = TRUE))
  expect_match(printed, "Stage-by-stage records available via $table.",
               fixed = TRUE)
})
