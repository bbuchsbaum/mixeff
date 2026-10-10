# lme4-style missing-value handling: lmm()/glmm() honour `na.action`
# (default getOption("na.action")), count only the variables the model uses,
# announce dropped rows with a typed `mm_rows_dropped` message (no silent
# surgery), and record them like lme4 (na.action(fit), model-frame attribute).

quiet <- mm_control(verbose = -1)

na_sleep <- function() {
  testthat::skip_if_not_installed("lme4")
  d <- lme4::sleepstudy
  d$Reaction[c(3, 40, 77)] <- NA
  d$Days[5] <- NA
  d$unused <- NA_real_ # NA in a column the model does not use: rows kept
  d
}

na_cbpp <- function() {
  testthat::skip_if_not_installed("lme4")
  d <- lme4::cbpp
  d$period[c(2, 9)] <- NA
  d$size[4] <- NA
  d$unused <- NA
  d
}

test_that("lmm() na.omit (the default) matches lmer()", {
  d <- na_sleep()
  m <- lmm(Reaction ~ Days + (Days | Subject), d, control = quiet)
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  expect_identical(nobs(m), nobs(r))
  expect_identical(nobs(m), 176L)
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-5)
  expect_equal(unname(m$theta), unname(lme4::getME(r, "theta")),
               tolerance = 1e-3)
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(r)), tolerance = 1e-6)
  expect_length(residuals(m), length(residuals(r)))
  expect_length(fitted(m), 176L)
  na <- stats::na.action(m)
  expect_s3_class(na, "omit")
  expect_equal(unname(as.integer(na)), c(3L, 5L, 40L, 77L))
  expect_identical(na, attr(stats::model.frame(r), "na.action"))
  expect_identical(attr(stats::model.frame(m), "na.action"), na)
})

test_that("lmm() na.exclude pads like lmer()", {
  d <- na_sleep()
  m <- lmm(Reaction ~ Days + (Days | Subject), d, na.action = na.exclude,
           control = quiet)
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), d,
                  na.action = na.exclude)
  expect_identical(nobs(m), nobs(r))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(r)), tolerance = 1e-6)
  expect_length(residuals(m), nrow(d))
  expect_length(residuals(m), length(residuals(r)))
  expect_identical(which(is.na(residuals(m))), which(is.na(residuals(r))))
  expect_equal(unname(fitted(m)), unname(fitted(r)), tolerance = 1e-4)
  expect_length(predict(m), nrow(d))
  expect_s3_class(stats::na.action(m), "exclude")
  # Downstream consumers use the fitted rows.
  expect_equal(dim(simulate(m, nsim = 2, seed = 1)), c(176L, 2L))
  expect_length(hatvalues(m), 176L)
  expect_length(cooks.distance(m), 176L)
  expect_identical(nobs(refit(m, simulate(m, seed = 2)[[1L]])), 176L)
  # update() of the stored frame keeps the record (padding persists).
  u <- update(m, . ~ . - Days)
  expect_s3_class(stats::na.action(u), "exclude")
  expect_length(residuals(u), nrow(d))
})

test_that("predict(newdata) is unaffected by the fit's na.action", {
  d <- na_sleep()
  m1 <- lmm(Reaction ~ Days + (Days | Subject), d, na.action = na.exclude,
            control = quiet)
  m2 <- lmm(Reaction ~ Days + (Days | Subject), d, control = quiet)
  nd <- lme4::sleepstudy[1:6, ]
  expect_equal(predict(m1, nd), predict(m2, nd))
  expect_length(predict(m1, nd), 6L)
})

test_that("glmm() na.omit / na.exclude match glmer()", {
  d <- na_cbpp()
  f <- cbind(incidence, size - incidence) ~ period + (1 | herd)
  g <- glmm(f, d, family = binomial, control = quiet)
  r <- lme4::glmer(f, d, family = binomial)
  expect_identical(nobs(g), nobs(r))
  expect_identical(nobs(g), 53L)
  expect_equal(fixef(g), lme4::fixef(r), tolerance = 1e-3)
  expect_equal(unname(g$theta), unname(lme4::getME(r, "theta")),
               tolerance = 1e-3)
  expect_equal(as.numeric(logLik(g)), as.numeric(logLik(r)), tolerance = 1e-5)
  expect_length(residuals(g), length(residuals(r)))
  expect_s3_class(stats::na.action(g), "omit")

  ge <- glmm(f, d, family = binomial, na.action = na.exclude, control = quiet)
  re <- lme4::glmer(f, d, family = binomial, na.action = na.exclude)
  expect_length(residuals(ge), nrow(d))
  expect_identical(which(is.na(residuals(ge))), which(is.na(residuals(re))))
  expect_identical(which(is.na(fitted(ge))), which(is.na(fitted(re))))
  expect_equal(dim(simulate(ge, nsim = 2, seed = 1)), c(53L, 2L))
})

test_that("NA in weights / offset counts as missing, like model.frame()", {
  d <- na_sleep()
  d$w <- 1
  d$w[10] <- NA
  m <- lmm(Reaction ~ Days + (1 | Subject), d, weights = w, control = quiet)
  r <- lme4::lmer(Reaction ~ Days + (1 | Subject), d, weights = w)
  expect_identical(nobs(m), nobs(r))
  expect_identical(nobs(m), 175L)
  cb <- na_cbpp()
  cb$off <- 0
  cb$off[20] <- NA
  g <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cb,
            family = binomial, offset = off, control = quiet)
  expect_identical(nobs(g), 52L)
})

test_that("dropped rows are announced with a typed message", {
  d <- na_sleep()
  ctl <- mm_control(verbose = 0)
  cnd <- tryCatch(
    lmm(Reaction ~ Days + (Days | Subject), d, control = ctl),
    mm_rows_dropped = function(c) c
  )
  expect_s3_class(cnd, "mm_rows_dropped")
  expect_s3_class(cnd, "mm_condition")
  expect_identical(cnd$n_dropped, 4L)
  expect_identical(cnd$n_total, 180L)
  expect_identical(cnd$columns, c("Reaction", "Days"))
  expect_match(conditionMessage(cnd), "Dropped 4 of 180 rows")
  expect_match(conditionMessage(cnd), "`Reaction` (3 NA)", fixed = TRUE)
  expect_match(conditionMessage(cnd), "na.omit", fixed = TRUE)
  expect_no_match(conditionMessage(cnd), "unused")

  cb <- na_cbpp()
  cnd <- tryCatch(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cb,
         family = binomial, na.action = na.exclude, control = ctl),
    mm_rows_dropped = function(c) c
  )
  expect_s3_class(cnd, "mm_rows_dropped")
  expect_match(conditionMessage(cnd), "padded with NA")

  # Silenced by verbose = -1, like the other mixeff notices.
  expect_no_message(lmm(Reaction ~ Days + (Days | Subject), d,
                        control = quiet))
  # No message when nothing is dropped (NA only in an unused column).
  d2 <- lme4::sleepstudy
  d2$unused <- NA
  expect_no_message(
    suppressMessages(
      lmm(Reaction ~ Days + (Days | Subject), d2, control = ctl),
      classes = "mm_explanation_notice"
    ),
    class = "mm_rows_dropped"
  )
})

test_that("na.fail and na.pass are refused with a typed mm_data_error", {
  d <- na_sleep()
  e <- expect_error(
    lmm(Reaction ~ Days + (Days | Subject), d, na.action = na.fail,
        control = quiet),
    class = "mm_data_error"
  )
  expect_identical(e$reason_code, "na_action_failed")
  expect_match(conditionMessage(e), "na.fail", fixed = TRUE)
  e <- expect_error(
    lmm(Reaction ~ Days + (Days | Subject), d, na.action = "na.pass",
        control = quiet),
    class = "mm_data_error"
  )
  expect_identical(e$reason_code, "na_action_keeps_missing")
  expect_error(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         na_cbpp(), family = binomial, na.action = na.fail, control = quiet),
    class = "mm_data_error"
  )
  # na.fail on complete data is fine (lme4 too).
  expect_identical(
    nobs(lmm(Reaction ~ Days + (1 | Subject), lme4::sleepstudy,
             na.action = na.fail, control = quiet)),
    180L
  )
})

test_that("options(na.action = \"na.fail\") is respected", {
  d <- na_sleep()
  op <- options(na.action = "na.fail")
  on.exit(options(op), add = TRUE)
  expect_error(lmm(Reaction ~ Days + (Days | Subject), d, control = quiet),
               class = "mm_data_error")
  expect_error(lme4::lmer(Reaction ~ Days + (Days | Subject), d))
  expect_error(
    glmm(cbind(incidence, size - incidence) ~ period + (1 | herd),
         na_cbpp(), family = binomial, control = quiet),
    class = "mm_data_error"
  )
  # An explicit argument overrides the option.
  expect_identical(
    nobs(lmm(Reaction ~ Days + (Days | Subject), d, na.action = na.omit,
             control = quiet)),
    176L
  )
})

test_that("compare()/anova() refuse models fitted to different rows", {
  d <- na_sleep()
  set.seed(3)
  d$z <- rnorm(nrow(d))
  d$z[20:24] <- NA
  m1 <- lmm(Reaction ~ Days + (1 | Subject), d, REML = FALSE, control = quiet)
  m2 <- lmm(Reaction ~ Days + z + (1 | Subject), d, REML = FALSE,
            control = quiet)
  expect_false(identical(nobs(m1), nobs(m2)))
  e <- expect_error(compare(m1, m2), class = "mm_arg_error")
  expect_identical(e$reason_code, "different_nobs")
  expect_match(conditionMessage(e), "same size of dataset")
  expect_error(anova(m1, m2), class = "mm_arg_error")
  # lme4 refuses the same comparison.
  r1 <- lme4::lmer(Reaction ~ Days + (1 | Subject), d, REML = FALSE)
  r2 <- lme4::lmer(Reaction ~ Days + z + (1 | Subject), d, REML = FALSE)
  expect_error(anova(r1, r2), "same size of dataset")
})

test_that("update() carries the fit's na.action over to new data", {
  d <- na_sleep()
  m <- lmm(Reaction ~ Days + (Days | Subject), d, na.action = na.exclude,
           control = quiet)
  u <- update(m, data = d)
  expect_s3_class(stats::na.action(u), "exclude")
  expect_length(residuals(u), nrow(d))
  o <- update(m, data = d, na.action = na.omit)
  expect_s3_class(stats::na.action(o), "omit")
  expect_length(residuals(o), 176L)
})
