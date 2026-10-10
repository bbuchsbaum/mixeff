# PR-4 coverage: hat values / Cook's distance for NB and Gamma GLMMs (vs lme4
# 2.1-0's factors),
# multi-column and failing deletion groups, parallel deletion, and the
# dfbeta()/dfbetas()/print() edge cases.

cp4i_ctrl <- function() mm_control(verbose = -1)

test_that("Gamma hat values and Cook's distances use the dispersion", {
  case <- mm_disp_case("influence_gamma_log")
  d <- case$data()
  ga <- glmm(case$formula, d, family = case$family, control = cp4i_ctrl())
  h <- hatvalues(ga)
  expect_true(all(h > 0 & h < 1))
  expect_gt(sum(h), 2)
  expect_lt(sum(h), 2 + 8)
  cg <- cooks.distance(ga)
  expect_true(all(is.finite(cg) & cg >= 0))
  # Cook's distance scales the squared Pearson residuals by the dispersion.
  parts <- mixeff:::mm_influence_parts(ga)
  expect_equal(parts$dispersion, ga$dispersion^2)
  expect_equal(unname(cg),
               unname((parts$pearson / (1 - h))^2 * h / (parts$dispersion * 2)))
  # The PIRLS working weights carry 1/phi (lme4 >= 2.1-0).
  expect_equal(parts$working, unname(weights(ga, type = "working")))
  # Reference: the working-weight hat matrix (mixeff's documented GLMM
  # convention; lme4's own hatvalues() puts the prior weights outside the
  # factorization) built from lme4 2.1-0's factors at its fit -- live with
  # lme4 >= 2.1-0, else the stored lme4 2.1-0 values.
  ref <- mm_disp_ref("influence_gamma_log")
  expect_equal(unname(h), ref$hat_working, tolerance = 1e-3)
  expect_equal(unname(cg), ref$cooks_working, tolerance = 1e-2)
})

cp4i_data <- function(seed = 482L) {
  set.seed(seed)
  g <- gl(6, 4)
  f <- factor(rep(c("a", "b"), 12))
  x <- rnorm(24)
  y <- 1 + 0.5 * x + rnorm(6, sd = 0.8)[g] + rnorm(24, sd = 0.4)
  data.frame(y = y, x = x, g = g, f = f)
}

test_that("deletion over two grouping columns and failing refits", {
  d <- cp4i_data()
  fit <- lmm(y ~ x + f + (1 | g), d, control = cp4i_ctrl())
  inf <- influence(fit, groups = c("g", "f"))
  expect_identical(inf$groups, "g.f")
  expect_length(inf$deleted, 12L)
  expect_identical(inf$deleted[[1]], "1.a")
  expect_true(all(inf$converged))
  # Deleting a whole level of the fixed factor `f` leaves a one-level
  # factor: those refits fail and are recorded as NA rows.
  bad <- influence(fit, groups = "f")
  expect_identical(unname(bad$fit.status), c("refit_error", "refit_error"))
  expect_false(any(bad$converged))
  expect_true(all(is.na(bad[["fixed.effects[-f]"]])))
  out <- capture.output(print(bad))
  expect_match(out[[1]], "deletion of 2 f level\\(s\\); 2 refit\\(s\\) not converged")
})

test_that("parallel deletion gives the same refits", {
  skip_on_os("windows")
  d <- cp4i_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4i_ctrl())
  serial <- influence(fit, groups = "g")
  par <- influence(fit, groups = "g", ncores = 2L)
  expect_equal(par[["fixed.effects[-g]"]], serial[["fixed.effects[-g]"]])
})

test_that("dfbeta()/dfbetas()/print() edge cases", {
  d <- cp4i_data()
  fit <- lmm(y ~ 1 + (1 | g), d, control = cp4i_ctrl())
  hat_only <- influence(fit, do.coef = FALSE)
  out <- capture.output(print(hat_only))
  expect_match(out[[1]], "hat values only")
  inf <- influence(fit, groups = "g")
  # One fixed effect: dfbetas() keeps a one-column matrix.
  dfs <- dfbetas(inf)
  expect_identical(dim(dfs), c(6L, 1L))
  se <- vapply(inf[["vcov[-g]"]], function(v) sqrt(v[1, 1]), 0)
  expect_equal(unname(dfs[, 1]), unname(dfbeta(inf)[, 1] / se))
  # An empty random-term basis is labelled as the intercept.
  nb <- fit
  nb$artifact$semantic_model$random_terms[[1]]$basis <- list()
  expect_identical(names(mixeff:::mm_influence_vc(nb)),
                   c("sigma^2", "v[(Intercept)]"))
})
