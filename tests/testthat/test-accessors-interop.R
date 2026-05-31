# Base/lme4-style accessors used by downstream tooling (broom.mixed, merTools,
# base step()). The contract here is mostly SHAPE parity with lme4; numeric
# values are checked loosely because the REML optimizer drifts sub-tolerance
# from lme4 on correlated-slope models (documented).

mk_acc_fit <- function() {
  suppressMessages(
    lmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy,
        control = mm_control(verbose = -1))
  )
}

test_that("ngrps returns per-grouping-factor level counts", {
  m <- mk_acc_fit()
  g <- ngrps(m)
  expect_type(g, "integer")
  expect_named(g, "Subject")
  expect_equal(unname(g[["Subject"]]), 18L)
})

test_that("ngrps counts each crossed grouping factor", {
  m <- suppressMessages(
    lmm(diameter ~ 1 + (1 | plate) + (1 | sample), lme4::Penicillin,
        control = mm_control(verbose = -1))
  )
  g <- ngrps(m)
  expect_equal(sort(unname(g)), c(6L, 24L))  # 6 samples, 24 plates
})

test_that("extractAIC returns c(edf, AIC) consistent with AIC()", {
  m <- mk_acc_fit()
  ea <- extractAIC(m)
  expect_length(ea, 2L)
  expect_equal(ea[[1L]], m$dof)            # edf
  expect_equal(ea[[2L]], unname(AIC(m)))   # internally consistent
})

test_that("terms() returns the fixed-effect terms", {
  m <- mk_acc_fit()
  tt <- terms(m)
  expect_s3_class(tt, "terms")
  expect_true("Days" %in% all.vars(tt))
})

test_that("weights() returns prior weights or NULL", {
  m <- mk_acc_fit()
  expect_null(weights(m))
  w <- runif(nrow(lme4::sleepstudy), 0.5, 2)
  mw <- suppressMessages(
    lmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy,
        weights = w, control = mm_control(verbose = -1))
  )
  expect_equal(weights(mw), w)
})

test_that("as.data.frame(VarCorr) matches lme4's long shape", {
  testthat::skip_if_not_installed("lme4")
  m  <- mk_acc_fit()
  m4 <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject),
                                    lme4::sleepstudy))
  d  <- as.data.frame(VarCorr(m))
  d4 <- as.data.frame(lme4::VarCorr(m4))

  expect_named(d, c("grp", "var1", "var2", "vcov", "sdcor"))
  expect_equal(nrow(d), nrow(d4))                 # 2 var + 1 cov + Residual
  expect_true("Residual" %in% d$grp)
  expect_true(any(!is.na(d$var2)))                # off-diagonal covariance row
  # loose numeric parity (optimizer drift)
  expect_equal(sort(d$vcov), sort(d4$vcov), tolerance = 1e-2)
  expect_equal(sort(d$sdcor), sort(d4$sdcor), tolerance = 1e-2)
})

test_that("as.data.frame(ranef, condVar=TRUE) matches lme4's long shape", {
  testthat::skip_if_not_installed("lme4")
  m  <- mk_acc_fit()
  m4 <- suppressMessages(lme4::lmer(Reaction ~ Days + (Days | Subject),
                                    lme4::sleepstudy))
  d  <- as.data.frame(ranef(m, condVar = TRUE))
  d4 <- as.data.frame(lme4::ranef(m4))

  expect_named(d, c("grpvar", "term", "grp", "condval", "condsd"))
  expect_equal(nrow(d), nrow(d4))
  expect_true(all(is.finite(d$condsd)))
  expect_true(all(d$condsd > 0))
  expect_equal(sort(d$condval), sort(d4$condval), tolerance = 5e-2)
})

test_that("as.data.frame(ranef) without condVar yields NA condsd, not an error", {
  m <- mk_acc_fit()
  d <- as.data.frame(ranef(m))
  expect_named(d, c("grpvar", "term", "grp", "condval", "condsd"))
  expect_true(all(is.na(d$condsd)))
})
