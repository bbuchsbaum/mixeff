# lme4 2.1-0's glmerControl(disp_method, disp_dof_correction, maxPhiIter)
# through mm_control(disp_method, disp_dof_correction, max_phi_iter).
# References: live glmer(..., control = glmerControl(...)) under
# lme4 >= 2.1-0, otherwise the stored lme4 2.1-0 values in
# fixtures/lme4-2.1-dispersion-refs.json (see helper-lme4-dispersion.R).

# Evaluate `expr`, muffling its messages; return the value and whether an
# mm_control_ignored_notice was among them.
with_ignored_notice <- function(expr) {
  seen <- FALSE
  value <- withCallingHandlers(expr, message = function(m) {
    if (inherits(m, "mm_control_ignored_notice")) seen <<- TRUE
    invokeRestart("muffleMessage")
  })
  list(value = value, notice = seen)
}

ctrl_fit <- function(key, verbose = -1L) {
  case <- mm_disp_case(key)
  glmm(case$formula, case$data(), family = case$family,
       control = mm_disp_mm_control(case, verbose = verbose))
}

test_that("mm_control() validates the dispersion controls with typed errors", {
  bad <- list(
    list(disp_method = "buggy"), list(disp_method = c("moment", "old/buggy")),
    list(disp_method = NA_character_), list(disp_method = 1),
    list(disp_dof_correction = NA), list(disp_dof_correction = "yes"),
    list(disp_dof_correction = c(TRUE, FALSE)),
    list(max_phi_iter = 0), list(max_phi_iter = 2.5), list(max_phi_iter = -1),
    list(max_phi_iter = NA_real_), list(max_phi_iter = Inf),
    list(max_phi_iter = "10"), list(max_phi_iter = 1:2)
  )
  for (args in bad) {
    expect_error(do.call(mm_control, args), class = "mm_arg_error",
                 info = deparse(args))
  }
  ctl <- mm_control(disp_method = "old/buggy", disp_dof_correction = FALSE,
                    max_phi_iter = 7)
  expect_identical(ctl$disp_method, "old/buggy")
  expect_identical(ctl$disp_dof_correction, FALSE)
  expect_identical(ctl$max_phi_iter, 7L)
  # Unset controls stay off the wire (engine defaults).
  expect_false(any(c("disp_method", "disp_dof_correction", "max_phi_iter") %in%
                     names(mm_control())))
  # mm_validate_control() (used by lmm()/glmm()) keeps them.
  expect_identical(mixeff:::mm_validate_control(unclass(ctl))[names(ctl)],
                   ctl[names(ctl)])
})

test_that("dispersion controls reproduce lme4 2.1-0 glmerControl() fits", {
  keys <- c("ctrl_gamma_oldbuggy", "ctrl_gamma_nodof", "ctrl_gamma_maxphi3",
            "ctrl_gamma_nodof_maxphi5", "ctrl_ig_oldbuggy", "ctrl_ig_nodof")
  for (key in keys) {
    m <- ctrl_fit(key)
    ref <- mm_disp_ref(key)
    # Looser than the default-control parity: under "old/buggy" the
    # engine's optimum agrees with lme4's to ~2e-4 in the estimates (both
    # sit at the same lme4 deviance to 1e-8 relative) and its reported
    # logLik to ~3e-6 relative.
    expect_equal(unname(fixef(m)), ref$fixef, tolerance = 1e-3, info = key)
    expect_equal(unname(getME(m, "theta")), ref$theta, tolerance = 1e-3,
                 info = key)
    expect_equal(sigma(m), ref$sigma, tolerance = 2e-4, info = key)
    expect_equal(as.numeric(logLik(m)), ref$logLik, tolerance = 1e-5,
                 info = key)
    expect_equal(deviance(m), ref$deviance, tolerance = 1e-4, info = key)
    expect_equal(AIC(m), ref$AIC, tolerance = 1e-5, info = key)
    # lme4 2.1-0 reports theta itself as the SD for every method.
    expect_equal(as.data.frame(VarCorr(m))$sdcor, ref$varcorr_sd,
                 tolerance = 1e-3, info = key)
    dims <- getME(m, "devcomp")$dims[mm_disp_dims]
    expect_equal(unname(as.numeric(dims)), ref$dims, info = key)
  }
})

test_that("the controls change the fit, and explicit defaults change nothing", {
  base <- ctrl_fit("scale_gamma_log")
  case <- mm_disp_case("scale_gamma_log")
  explicit <- glmm(case$formula, case$data(), family = case$family,
                   control = mm_control(verbose = -1, disp_method = "moment",
                                        disp_dof_correction = TRUE,
                                        max_phi_iter = 100))
  expect_identical(fixef(explicit), fixef(base))
  expect_identical(sigma(explicit), sigma(base))
  expect_identical(logLik(explicit), logLik(base))
  for (key in c("ctrl_gamma_oldbuggy", "ctrl_gamma_nodof",
                "ctrl_gamma_maxphi3")) {
    expect_false(isTRUE(all.equal(sigma(ctrl_fit(key)), sigma(base))),
                 info = key)
  }
})

test_that("controls are ignored, with a typed notice, without a free dispersion", {
  set.seed(31)
  d <- data.frame(g = gl(10, 8), x = rnorm(80))
  d$y <- rpois(80, exp(0.5 + 0.3 * d$x + rnorm(10, sd = 0.3)[d$g]))
  plain <- glmm(y ~ x + (1 | g), d, family = poisson(),
                control = mm_control(verbose = -1))
  res <- with_ignored_notice(
    glmm(y ~ x + (1 | g), d, family = poisson(),
         control = mm_control(disp_method = "old/buggy", max_phi_iter = 3))
  )
  expect_true(res$notice)
  expect_identical(fixef(res$value), fixef(plain))
  expect_identical(logLik(res$value), logLik(plain))
  # verbose = -1 silences it, like the other fit notices.
  expect_false(with_ignored_notice(
    glmm(y ~ x + (1 | g), d, family = poisson(),
         control = mm_control(verbose = -1, disp_dof_correction = FALSE))
  )$notice)
  expect_true(with_ignored_notice(
    lmm(x ~ 1 + (1 | g), d, control = mm_control(disp_method = "old/buggy"))
  )$notice)
  # Free-dispersion families take them without a notice.
  case <- mm_disp_case("ctrl_gamma_nodof")
  expect_false(with_ignored_notice(
    glmm(case$formula, case$data(), family = case$family,
         control = mm_disp_mm_control(case, verbose = 0L))
  )$notice)
})

test_that("refit(), update() and simulate()-based refits keep the controls", {
  key <- "ctrl_gamma_oldbuggy"
  case <- mm_disp_case(key)
  d <- case$data()
  m <- ctrl_fit(key)
  expect_identical(m$control$disp_method, "old/buggy")
  # refit() on the original response reproduces the fit.
  r <- refit(m, d$y)
  expect_identical(r$control$disp_method, "old/buggy")
  expect_equal(fixef(r), fixef(m), tolerance = 1e-8)
  expect_equal(sigma(r), sigma(m), tolerance = 1e-8)
  # ... and on a new response equals a fresh fit with the same control.
  ynew <- simulate(m, seed = 2)[[1]]
  r2 <- refit(m, ynew)
  d2 <- transform(d, y = ynew)
  fresh <- glmm(case$formula, d2, family = case$family,
                control = mm_disp_mm_control(case))
  expect_equal(fixef(r2), fixef(fresh), tolerance = 1e-8)
  expect_equal(sigma(r2), sigma(fresh), tolerance = 1e-8)
  # update() carries the control unless a new one is given.
  u <- update(m, data = d2)
  expect_identical(u$control$disp_method, "old/buggy")
  expect_equal(fixef(u), fixef(fresh), tolerance = 1e-8)
  u0 <- update(m, control = mm_control(verbose = -1))
  expect_null(u0$control$disp_method)
  expect_equal(sigma(u0), sigma(ctrl_fit("scale_gamma_log")), tolerance = 1e-8)
  expect_identical(unname(getME(u, "devcomp")$dims[["dispProfile"]]), 0L)
})

test_that("live glmer() agrees for the dispersion controls", {
  skip_if_not(mm_disp_lme4_is_2_1(), "needs lme4 >= 2.1-0")
  case <- mm_disp_case("ctrl_gamma_nodof_maxphi5")
  d <- case$data()
  g <- lme4::glmer(case$formula, d, family = case$family,
                   control = lme4::glmerControl(disp_method = "moment",
                                                disp_dof_correction = FALSE,
                                                maxPhiIter = 5))
  m <- ctrl_fit("ctrl_gamma_nodof_maxphi5")
  expect_equal(sigma(m), sigma(g), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(g)), tolerance = 1e-6)
  expect_equal(unname(weights(m, type = "working")),
               unname(stats::weights(g, type = "working")), tolerance = 1e-4)
  # "old/buggy": working weights without 1/phi, and simulate() draws the
  # random effects from theta itself, as lme4 2.1-0 does.
  case <- mm_disp_case("ctrl_gamma_oldbuggy")
  d <- case$data()
  g <- suppressWarnings(lme4::glmer(case$formula, d, family = case$family,
                   control = lme4::glmerControl(disp_method = "old/buggy")))
  m <- ctrl_fit("ctrl_gamma_oldbuggy")
  expect_equal(unname(weights(m, type = "working")),
               unname(stats::weights(g, type = "working")), tolerance = 1e-3)
  expect_equal(unname(as.matrix(simulate(m, nsim = 2, seed = 3))),
               unname(as.matrix(stats::simulate(g, nsim = 2, seed = 3))),
               tolerance = 1e-3)
})
