# lme4 / lmerTest parity for extractor shapes and semantics (pre-CRAN audit
# §3 residuals/sigma/deviance and §4 items 3, 4, 8, 9, 10, 12).

shape_sleep <- function() {
  skip_if_not_installed("lme4")
  d <- lme4::sleepstudy
  d$f <- factor(rep(c("a", "b", "c"), length.out = nrow(d)))
  d
}

shape_close <- function(a, b, tol = 1e-3) {
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)), tolerance = tol,
               ignore_attr = TRUE)
}

test_that("LMM residual types and scaling match lme4", {
  d <- shape_sleep()
  d$w <- rep(c(1, 2, 0.5), length.out = nrow(d))
  fit <- lmm(Reaction ~ Days + (Days | Subject), d, weights = w,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d, weights = w)
  for (type in c("response", "pearson", "deviance", "working")) {
    shape_close(residuals(fit, type = type), residuals(m, type = type), 1e-4)
    shape_close(residuals(fit, type = type, scaled = TRUE),
                residuals(m, type = type, scaled = TRUE), 1e-4)
  }
  # Pearson is sqrt(w) * (y - mu): no division by sigma unless scaled.
  expect_equal(unname(residuals(fit, "pearson")),
               unname(sqrt(d$w) * residuals(fit, "response")))
  expect_equal(unname(residuals(fit, "pearson", scaled = TRUE)),
               unname(residuals(fit, "pearson") / sigma(fit)))
})

test_that("GLMM residuals, deviance, sigma, weights match glmer", {
  skip_if_not_installed("lme4")
  cb <- lme4::cbpp
  fit <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cb,
              family = binomial(), method = "joint_laplace",
              control = mm_control(verbose = -1))
  m <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cb, family = binomial())
  expect_equal(unname(residuals(fit)), unname(residuals(m)), tolerance = 5e-3)
  for (type in c("deviance", "pearson", "working", "response")) {
    expect_equal(unname(residuals(fit, type = type)),
                 unname(residuals(m, type = type)), tolerance = 5e-3)
  }
  expect_equal(deviance(fit), deviance(m), tolerance = 1e-3)
  expect_equal(deviance(fit), sum(residuals(fit, "deviance")^2))
  expect_identical(sigma(fit), 1)
  expect_equal(weights(fit), weights(m))
  expect_equal(weights(fit, type = "working"), weights(m, type = "working"),
               tolerance = 5e-3)
  expect_equal(family(fit)$family, "binomial")
  expect_equal(family(fit)$link, "logit")
  # The anova LRT still uses -2 logLik, not the residual deviance.
  fit0 <- glmm(cbind(incidence, size - incidence) ~ 1 + (1 | herd), cb,
               family = binomial(), method = "joint_laplace",
               control = mm_control(verbose = -1))
  m0 <- lme4::glmer(cbind(incidence, size - incidence) ~ 1 + (1 | herd),
                    cb, family = binomial())
  a <- anova(fit0, fit)
  a4 <- anova(m0, m)
  expect_s3_class(a, "anova")
  expect_s3_class(a, "data.frame")
  expect_identical(names(a), names(a4))
  expect_identical(rownames(a), c("fit0", "fit"))
  expect_equal(a$Chisq[[2L]], a4$Chisq[[2L]], tolerance = 1e-2)
  expect_equal(a$deviance, -2 * c(fit0$logLik, fit$logLik))
  expect_true(is.data.frame(a$table))
})

test_that("Poisson GLMM residual types match glmer", {
  skip_if_not_installed("lme4")
  set.seed(11)
  g <- factor(rep(seq_len(15), each = 8))
  x <- rnorm(120)
  d <- data.frame(y = rpois(120, exp(0.5 + 0.4 * x + rnorm(15, 0, .4)[g])),
                  x = x, g = g)
  fit <- glmm(y ~ x + (1 | g), d, family = poisson(),
              method = "joint_laplace", control = mm_control(verbose = -1))
  m <- lme4::glmer(y ~ x + (1 | g), d, family = poisson())
  for (type in c("deviance", "pearson", "working", "response")) {
    expect_equal(unname(residuals(fit, type = type)),
                 unname(residuals(m, type = type)), tolerance = 5e-3)
  }
  expect_equal(deviance(fit), deviance(m), tolerance = 1e-3)
  expect_identical(sigma(fit), 1)
  # GLMM conditional variances (Laplace) now available and match lme4.
  re <- ranef(fit, condVar = TRUE)
  pv <- attr(re$g, "postVar")
  expect_false(anyNA(pv))
  expect_equal(as.numeric(pv), as.numeric(attr(lme4::ranef(m, condVar = TRUE)$g,
                                                "postVar")), tolerance = 5e-3)
  expect_false(is_singular(fit))
  expect_identical(is_singular(fit, tol = 10), TRUE)
  expect_equal(getME(fit, "devcomp")$cmp[["ldL2"]],
               getME(m, "devcomp")$cmp[["ldL2"]], tolerance = 1e-3)
})

test_that("negative-binomial sigma is 1 and theta is a getME component", {
  set.seed(5)
  g <- factor(rep(seq_len(12), each = 10))
  x <- rnorm(120)
  d <- data.frame(y = rnbinom(120, size = 2,
                             mu = exp(1 + 0.3 * x)),
                  x = x, g = g)
  fit <- glmm(y ~ x + (1 | g), d, family = mm_negative_binomial(),
              control = mm_control(verbose = -1))
  expect_identical(sigma(fit), 1)
  expect_equal(getME(fit, "glmer.nb.theta"), fit$family$nb_theta)
  expect_true(grepl("^Negative Binomial", family(fit)$family))
  r <- residuals(fit, "deviance")
  fam <- family(fit)
  expect_equal(unname(r^2),
               unname(fam$dev.resids(d$y, fitted(fit), rep(1, 120))),
               tolerance = 1e-8)
  expect_equal(deviance(fit), sum(r^2))
  expect_error(getME(lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1)),
                     "glmer.nb.theta"),
               class = "mm_inference_unavailable")
})

test_that("coef() repeats fixed-only columns in lme4's order", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + f + (1 | Subject), d,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + f + (1 | Subject), d)
  cf <- coef(fit)$Subject
  cf4 <- coef(m)$Subject
  expect_identical(colnames(cf), colnames(cf4))
  expect_identical(rownames(cf), rownames(cf4))
  shape_close(cf, cf4, 1e-4)
  # random-only column: (0 + Days | Subject) with no fixed Days
  fit2 <- lmm(Reaction ~ 1 + (1 | Subject) + (0 + Days | Subject), d,
              control = mm_control(verbose = -1))
  m2 <- lme4::lmer(Reaction ~ 1 + (1 | Subject) + (0 + Days | Subject), d)
  expect_setequal(colnames(coef(fit2)$Subject), colnames(coef(m2)$Subject))
})

test_that("VarCorr() has lme4's shape, attributes, print and long form", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  vc <- VarCorr(fit)
  vc4 <- lme4::VarCorr(m)
  expect_identical(names(vc), names(vc4))
  expect_true(is.matrix(vc$Subject))
  expect_identical(dimnames(vc$Subject), dimnames(vc4$Subject))
  shape_close(vc$Subject, vc4$Subject, 1e-3)
  expect_equal(attr(vc$Subject, "stddev"), attr(vc4$Subject, "stddev"),
               tolerance = 1e-3)
  shape_close(attr(vc$Subject, "correlation"),
              attr(vc4$Subject, "correlation"), 1e-2)
  expect_equal(attr(vc, "sc"), attr(vc4, "sc"), tolerance = 1e-4)
  expect_identical(attr(vc, "useSc"), TRUE)
  # backwards-compatible mixeff views
  expect_true(is.data.frame(vc$table))
  expect_equal(vc$residual_sd, sigma(fit))
  expect_equal(as.data.frame(vc)[, c("grp", "var1", "var2")],
               as.data.frame(vc4)[, c("grp", "var1", "var2")])
  out <- capture.output(print(vc))
  out4 <- capture.output(print(vc4))
  expect_identical(out[[1L]], out4[[1L]])
  expect_identical(length(out), length(out4))
  expect_match(out[[2L]], "^ Subject  \\(Intercept\\)")
  # GLMM without a free scale: useSc FALSE and no Residual row
  cb <- lme4::cbpp
  g <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cb,
            family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1))
  expect_identical(attr(VarCorr(g), "useSc"), FALSE)
  expect_identical(attr(VarCorr(g), "sc"), 1)
  expect_false("Residual" %in% as.data.frame(VarCorr(g))$grp)
})

test_that("AIC()/BIC() with several models return stats' data frame", {
  d <- shape_sleep()
  f1 <- lmm(Reaction ~ Days + (Days | Subject), d,
            control = mm_control(verbose = -1))
  f0 <- lmm(Reaction ~ 1 + (Days | Subject), d,
            control = mm_control(verbose = -1))
  m1 <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  m0 <- lme4::lmer(Reaction ~ 1 + (Days | Subject), d)
  a <- AIC(f0, f1)
  expect_identical(names(a), c("df", "AIC"))
  expect_identical(rownames(a), c("f0", "f1"))
  expect_equal(unname(as.matrix(a)), unname(as.matrix(AIC(m0, m1))),
               tolerance = 1e-4)
  b <- BIC(f0, f1)
  expect_equal(unname(as.matrix(b)), unname(as.matrix(BIC(m0, m1))),
               tolerance = 1e-4)
  expect_equal(AIC(f1, k = 3), -2 * f1$logLik + 3 * f1$dof)
})

test_that("REMLcrit(), weights(), family(), model.frame() match lme4", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  fml <- lmm(Reaction ~ Days + (Days | Subject), d, REML = FALSE,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  mml <- lme4::lmer(Reaction ~ Days + (Days | Subject), d, REML = FALSE)
  expect_equal(mixeff::REMLcrit(fit), lme4::REMLcrit(m), tolerance = 1e-5)
  expect_equal(mixeff::REMLcrit(fml), lme4::REMLcrit(mml), tolerance = 1e-5)
  expect_equal(mixeff::REMLcrit(m), lme4::REMLcrit(m))
  expect_equal(weights(fit), rep(1, nobs(fit)))
  expect_error(weights(fit, type = "working"), class = "mm_arg_error")
  expect_identical(family(fit)$family, "gaussian")
  fitlog <- lmm(log(Reaction) ~ Days + (1 | Subject), d,
                control = mm_control(verbose = -1))
  mlog <- lme4::lmer(log(Reaction) ~ Days + (1 | Subject), d)
  mf <- model.frame(fitlog)
  expect_identical(names(mf), names(model.frame(mlog)))
  expect_false(is.null(attr(mf, "terms")))
  expect_identical(attr(attr(mf, "terms"), "varnames.fixed"),
                   attr(attr(model.frame(mlog), "terms"), "varnames.fixed"))
  expect_equal(mf[["log(Reaction)"]], log(d$Reaction))
  expect_identical(names(model.frame(fitlog, fixed.only = TRUE)),
                   names(model.frame(mlog, fixed.only = TRUE)))
  expect_equal(unname(getME(fitlog, "y")), log(d$Reaction))
})

test_that("getME() components match lme4 for a correlated-slope LMM", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  for (nm in c("X", "Z", "u", "b", "Lambdat", "Lind", "theta", "lower", "RX",
               "Gp", "Tp", "p_i", "l_i", "q_i", "m_i", "y", "mu", "offset",
               "A", "n_rtrms", "n_rfacs", "N", "p", "q", "k", "m")) {
    shape_close(getME(fit, nm), getME(m, nm), 2e-3)
  }
  expect_identical(names(getME(fit, "theta")), names(getME(m, "theta")))
  expect_identical(names(getME(fit, "lower")), names(getME(m, "lower")))
  shape_close(crossprod(getME(fit, "RZX")), crossprod(getME(m, "RZX")), 2e-3)
  expect_s4_class(getME(fit, "L"), "CHMfactor")
  expect_equal(
    as.numeric(Matrix::determinant(getME(fit, "L"), sqrt = TRUE)$modulus),
    as.numeric(Matrix::determinant(getME(m, "L"), sqrt = TRUE)$modulus),
    tolerance = 1e-4
  )
  dc <- getME(fit, "devcomp")
  dc4 <- getME(m, "devcomp")
  expect_identical(names(dc$cmp), names(dc4$cmp))
  expect_identical(dc$dims, dc4$dims)
  expect_equal(dc$cmp, dc4$cmp, tolerance = 1e-4)
  expect_equal(getME(fit, "Tlist")$Subject, getME(m, "Tlist")[[1L]],
               tolerance = 2e-3, ignore_attr = TRUE)
  expect_equal(getME(fit, "ST")$Subject, getME(m, "ST")$Subject,
               tolerance = 2e-3, ignore_attr = TRUE)
  expect_identical(names(getME(fit, "Ztlist")), names(getME(m, "Ztlist")))
  expect_identical(getME(fit, "cnms"), getME(m, "cnms"))
  expect_identical(getME(fit, "is_REML"), TRUE)
  expect_equal(getME(fit, "sigma"), sigma(m), tolerance = 1e-5)
  all <- getME(fit, "ALL")
  expect_true(all(c("u", "b", "L", "RX", "RZX", "devcomp", "lower") %in%
                    names(all)))
  err <- tryCatch(getME(fit, "devfun"), error = identity)
  expect_s3_class(err, "mm_inference_unavailable")
  expect_identical(err$component, "devfun")
  expect_error(getME(fit, "no_such_thing"), class = "mm_arg_error")
})

test_that("getME() orders crossed terms like lme4 (theta order)", {
  skip_if_not_installed("lme4")
  set.seed(2)
  n <- 240
  d <- data.frame(g1 = factor(sample(1:5, n, TRUE)),
                  g2 = factor(sample(1:30, n, TRUE)), x = rnorm(n))
  d$y <- 1 + d$x + rnorm(5)[d$g1] + rnorm(30)[d$g2] + rnorm(n)
  f <- y ~ x + (1 | g1) + (x | g2)
  fit <- lmm(f, d, control = mm_control(verbose = -1))
  m <- lme4::lmer(f, d)
  shape_close(getME(fit, "Lambdat"), getME(m, "Lambdat"), 5e-3)
  shape_close(getME(fit, "Z"), getME(m, "Z"), 1e-8)
  shape_close(getME(fit, "b"), getME(m, "b"), 5e-3)
  expect_identical(names(getME(fit, "flist")), names(getME(m, "flist")))
  expect_equal(getME(fit, "devcomp")$cmp, getME(m, "devcomp")$cmp,
               tolerance = 1e-4)
})

test_that("is_singular() follows lme4::isSingular() and honours tol", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  for (tol in c(1e-4, 0.1, 0.25, 1)) {
    expect_identical(is_singular(fit, tol = tol), lme4::isSingular(m, tol = tol))
  }
  expect_error(is_singular(fit, tol = -1), class = "mm_arg_error")
})

test_that("summary() coefficients are lmerTest's matrix; ddf= is honoured", {
  skip_if_not_installed("lmerTest")
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days * f + (Days | Subject), d,
             control = mm_control(verbose = -1))
  m <- lmerTest::lmer(Reaction ~ Days * f + (Days | Subject), d)
  s <- summary(fit)
  s4 <- summary(m)
  expect_true(is.matrix(s$coefficients))
  expect_identical(colnames(s$coefficients), colnames(s4$coefficients))
  expect_equal(s$coefficients[, 1:4], s4$coefficients[, 1:4],
               tolerance = 1e-4)
  expect_true(is.data.frame(s$coef_table))
  expect_true("method" %in% names(s$coef_table))
  skr <- summary(fit, ddf = "Kenward-Roger")
  skr4 <- summary(m, ddf = "Kenward-Roger")
  expect_equal(skr$coefficients[, "df"], skr4$coefficients[, "df"],
               tolerance = 1e-4)
  expect_true(all(skr$coef_table$method == "kenward_roger"))
  sl <- summary(fit, ddf = "lme4")
  expect_identical(colnames(sl$coefficients),
                   colnames(coef(summary(m, ddf = "lme4"))))
  expect_equal(sl$coefficients, coef(summary(m, ddf = "lme4")),
               tolerance = 1e-4)
  expect_error(summary(fit, ddf = "bogus"), class = "mm_arg_error")
  expect_error(summary(fit, ddf = "lme4", method = "satterthwaite"),
               class = "mm_arg_error")
  # GLMM summary coefficients are lme4's z matrix
  cb <- lme4::cbpp
  g <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cb,
            family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1))
  g4 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                    cb, family = binomial())
  expect_identical(colnames(summary(g)$coefficients),
                   colnames(coef(summary(g4))))
})

test_that("anova() returns lme4/lmerTest-shaped data frames", {
  skip_if_not_installed("lmerTest")
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days * f + (Days | Subject), d,
             control = mm_control(verbose = -1))
  fit0 <- lmm(Reaction ~ 1 + (Days | Subject), d,
              control = mm_control(verbose = -1))
  m <- lmerTest::lmer(Reaction ~ Days * f + (Days | Subject), d)
  m0 <- lmerTest::lmer(Reaction ~ 1 + (Days | Subject), d)
  a <- anova(fit)
  a4 <- anova(m)
  expect_s3_class(a, "data.frame")
  expect_s3_class(a, "anova")
  expect_identical(names(a), names(a4))
  expect_identical(rownames(a), rownames(a4))
  expect_equal(unname(as.matrix(a)), unname(as.matrix(a4)), tolerance = 1e-3)
  expect_identical(attr(a, "heading"), attr(a4, "heading"))
  expect_true(is.data.frame(a$table))
  expect_identical(a$type, "III")
  akr <- anova(fit, ddf = "Kenward-Roger")
  expect_equal(unname(as.matrix(akr)),
               unname(as.matrix(anova(m, ddf = "Kenward-Roger"))),
               tolerance = 1e-3)
  expect_equal(unname(as.matrix(anova(fit, type = 1))),
               unname(as.matrix(anova(m, type = 1))), tolerance = 1e-3)
  al <- anova(fit, ddf = "lme4")
  al4 <- anova(m, ddf = "lme4")
  expect_identical(names(al), names(al4))
  expect_equal(unname(as.matrix(al)), unname(as.matrix(al4)), tolerance = 1e-3)
  lr <- suppressMessages(anova(fit0, fit))
  lr4 <- suppressMessages(anova(m0, m))
  expect_identical(names(lr), names(lr4))
  expect_identical(rownames(lr), c("fit0", "fit"))
  expect_equal(unname(as.matrix(lr)), unname(as.matrix(lr4)), tolerance = 1e-3)
  expect_true(is.data.frame(lr$table))
  expect_output(print(lr), "Models:")
})

test_that("confint() accepts lme4's method spellings", {
  d <- shape_sleep()
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  expect_equal(unclass(confint(fit, method = "Wald")),
               unclass(confint(fit)))
  m <- lme4::lmer(Reaction ~ Days + (Days | Subject), d)
  w4 <- confint(m, method = "Wald")
  expect_equal(unname(unclass(confint(fit, method = "Wald"))[, 1:2]),
               unname(w4[c("(Intercept)", "Days"), ]), tolerance = 1e-4,
               ignore_attr = TRUE)
  expect_error(confint(fit, method = "nope"))
})
