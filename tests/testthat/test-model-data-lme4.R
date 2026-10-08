# lme4 parity for the shared lmm()/glmm() data preparation (R/model-data.R):
# na.exclude padding, unused-level dropping, glmm() subset/na.action/
# contrasts, newdata robustness in predict(), lme4 grouping labels for
# interaction/nested grouping, logical predictors, ordered grouping factors,
# stateful formula terms, offsets, and partial re.form.

skip_if_not_installed("lme4")

q <- function(expr) suppressMessages(suppressWarnings(expr))
quiet <- mm_control(verbose = -1)

md_data <- function(seed = 11, n = 160) {
  set.seed(seed)
  d <- data.frame(
    x = rnorm(n),
    z = rnorm(n),
    g = factor(rep(seq_len(16), each = n / 16)),
    f = factor(sample(letters[1:3], n, TRUE)),
    a = factor(rep(1:4, n / 4)),
    b = factor(rep(1:2, each = n / 2)),
    l = rep(c(TRUE, FALSE), n / 2),
    t = runif(n, 1, 4)
  )
  u <- rnorm(16, sd = 0.7)
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + 0.4 * (d$f == "b") + d$l + u[d$g] +
    rnorm(n)
  d$cnt <- rpois(n, exp(0.2 + 0.3 * d$x + log(d$t) + u[d$g] / 2))
  d$yb <- rbinom(n, 1, plogis(0.3 * d$x + u[d$g]))
  d
}

lmer_q <- function(...) {
  cl <- match.call()
  cl[[1L]] <- quote(lme4::lmer)
  q(eval(cl, parent.frame()))
}

test_that("na.exclude pads residuals/fitted/predict like lme4", {
  d <- md_data()
  d$y[c(3, 40)] <- NA
  d$x[7] <- NA
  m <- q(lmm(y ~ x + (1 | g), d, na.action = na.exclude, control = quiet))
  r <- lmer_q(y ~ x + (1 | g), d, na.action = na.exclude)
  expect_length(residuals(m), nrow(d))
  expect_length(fitted(m), nrow(d))
  expect_length(predict(m), nrow(d))
  expect_equal(unname(which(is.na(residuals(m)))), c(3L, 7L, 40L))
  expect_equal(residuals(m), residuals(r), tolerance = 1e-4)
  expect_equal(fitted(m), fitted(r), tolerance = 1e-4)
  expect_equal(predict(m), predict(r), tolerance = 1e-4)
  expect_s3_class(m$na.action, "exclude")
  expect_identical(attr(m$model_frame, "na.action"), m$na.action)
  expect_equal(nobs(m), nrow(d) - 3L)
  # na.omit keeps the compact (lme4) shape
  m2 <- q(lmm(y ~ x + (1 | g), d, na.action = na.omit, control = quiet))
  expect_length(residuals(m2), nrow(d) - 3L)
})

test_that("emmeans reads the na.action record of an na.omit fit", {
  skip_if_not_installed("emmeans")
  d <- md_data()
  d$y[c(2, 5)] <- NA
  m <- q(lmm(y ~ f + (1 | g), d, na.action = na.omit, control = quiet))
  r <- lmer_q(y ~ f + (1 | g), d)
  em <- as.data.frame(emmeans::emmeans(m, ~ f))
  er <- as.data.frame(emmeans::emmeans(r, ~ f))
  expect_equal(em$emmean, er$emmean, tolerance = 1e-4)
})

test_that("unused factor levels are dropped after subset / NA removal", {
  d <- md_data()
  m <- q(lmm(y ~ f + (1 | g), d, subset = f != "b", control = quiet))
  r <- lmer_q(y ~ f + (1 | g), d, subset = f != "b")
  expect_identical(names(fixef(m)), names(lme4::fixef(r)))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-4)
  # grouping factor: an emptied level does not inflate ngrps
  m2 <- q(lmm(y ~ x + (1 | g), d, subset = g != "3", control = quiet))
  expect_identical(ngrps(m2), c(g = 15L))
  # ordered factors get contr.poly over the remaining levels
  d$o <- factor(as.integer(d$a), ordered = TRUE)
  m3 <- q(lmm(y ~ o + (1 | g), d, subset = o != "4", control = quiet))
  r3 <- lmer_q(y ~ o + (1 | g), d, subset = o != "4")
  expect_identical(names(fixef(m3)), names(lme4::fixef(r3)))
  expect_equal(fixef(m3), lme4::fixef(r3), tolerance = 1e-4)
})

test_that("glmm() honours subset, na.action and contrasts like lmm()", {
  d <- md_data()
  d$x2 <- d$x
  d$x2[5] <- NA
  # default: NA refused, consistently with lmm()
  expect_error(glmm(cnt ~ x2 + (1 | g), d, family = poisson, control = quiet),
               class = "mm_data_error")
  m <- q(glmm(cnt ~ x2 + (1 | g), d, family = poisson, subset = t > 2,
              na.action = na.exclude, method = "joint_laplace",
              contrasts = list(f = "contr.treatment"), control = quiet))
  r <- q(lme4::glmer(cnt ~ x2 + (1 | g), d, family = poisson,
                     subset = t > 2, na.action = na.exclude))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-3)
  expect_equal(nobs(m), nobs(r))
  expect_length(fitted(m), sum(d$t > 2))
  expect_equal(fitted(m), fitted(r), tolerance = 1e-3)
  expect_error(
    glmm(cnt ~ x + (1 | g), d, family = poisson,
         contrasts = list(f = "contr.sum"), control = quiet),
    class = "mm_arg_error"
  )
})

test_that("predict(newdata) coerces grouping/character columns and handles NA", {
  d <- md_data()
  d$gi <- as.integer(d$g)
  d$fc <- as.character(d$f)
  d$when <- Sys.Date() + seq_len(nrow(d))
  m <- q(lmm(y ~ x + fc + (1 | gi), d, control = quiet))
  r <- lmer_q(y ~ x + fc + (1 | gi), d)
  nd <- d[1:6, ]
  nd$x[2] <- NA
  expect_equal(predict(m, nd), predict(r, nd), tolerance = 1e-4)
  expect_true(is.na(predict(m, nd)[[2L]]))
  expect_equal(predict(m, nd, na.action = na.omit),
               predict(r, nd, na.action = na.omit), tolerance = 1e-4)
  expect_length(predict(m, nd, na.action = na.exclude), 6L)
  expect_error(predict(m, nd, na.action = na.fail))
  # character grouping column in newdata, integer at fit time
  nd2 <- nd[3:4, ]
  nd2$gi <- as.character(nd2$gi)
  expect_equal(predict(m, nd2), predict(r, nd[3:4, ]), tolerance = 1e-4)
  # population prediction ignores grouping columns and unused Date columns
  expect_equal(predict(m, nd[, c("x", "fc", "when")], re.form = NA),
               predict(r, nd, re.form = NA), tolerance = 1e-4)
  # unseen fixed-effect level is a typed refusal
  nd3 <- nd[1, ]
  nd3$fc <- "zzz"
  expect_error(predict(m, nd3), class = "mm_data_error")
})

test_that("interaction and nested grouping use lme4 names and levels", {
  d <- md_data()
  m <- q(lmm(y ~ x + (1 | a:b), d, control = quiet))
  r <- lmer_q(y ~ x + (1 | a:b), d)
  expect_identical(names(ranef(m)), names(lme4::ranef(r)))
  expect_identical(rownames(ranef(m)[["a:b"]]),
                   rownames(lme4::ranef(r)[["a:b"]]))
  expect_equal(ngrps(m), lme4::ngrps(r))
  expect_identical(unique(VarCorr(m)$table$group), "a:b")
  m2 <- q(lmm(y ~ x + (1 | a/b), d, control = quiet))
  r2 <- lmer_q(y ~ x + (1 | a/b), d)
  expect_identical(names(ranef(m2)), names(lme4::ranef(r2)))
  expect_identical(rownames(ranef(m2)[["b:a"]]),
                   rownames(lme4::ranef(r2)[["b:a"]]))
  expect_equal(ranef(m2)[["b:a"]][[1L]], lme4::ranef(r2)[["b:a"]][[1L]],
               tolerance = 1e-3)
  pv <- attr(ranef(m2, condVar = TRUE)[["b:a"]], "postVar")
  expect_identical(dimnames(pv)[[3L]], rownames(ranef(m2)[["b:a"]]))
  # GLMM predict(newdata) looks the interaction group up
  g <- q(glmm(yb ~ x + (1 | a:b), d, family = binomial,
              method = "joint_laplace", control = quiet))
  gr <- q(lme4::glmer(yb ~ x + (1 | a:b), d, family = binomial))
  nd <- d[c(1, 50, 90), ]
  expect_equal(predict(g, nd, type = "link"),
               predict(gr, nd, type = "link"), tolerance = 1e-3)
})

test_that("logical predictors are coded like model.matrix() (no intercept)", {
  d <- md_data()
  m <- q(lmm(y ~ 0 + l + x + (1 | g), d, control = quiet))
  r <- lmer_q(y ~ 0 + l + x + (1 | g), d)
  expect_identical(names(fixef(m)), c("lFALSE", "lTRUE", "x"))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-4)
  m2 <- q(lmm(y ~ x:l + (1 | g), d, control = quiet))
  r2 <- lmer_q(y ~ x:l + (1 | g), d)
  expect_equal(fixef(m2), lme4::fixef(r2), tolerance = 1e-4)
  nd <- data.frame(x = c(1, -1), l = c(TRUE, FALSE), g = c("1", "2"))
  expect_equal(predict(m2, nd), predict(r2, nd), tolerance = 1e-3)
})

test_that("ordered grouping factors are grouped as unordered factors", {
  d <- md_data()
  d$og <- factor(as.integer(d$g), ordered = TRUE)
  m <- q(lmm(y ~ x + (1 | og), d, control = quiet))
  r <- lmer_q(y ~ x + (1 | og), d)
  expect_false(is.ordered(m$model_frame$og))
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(r)), tolerance = 1e-5)
})

test_that("stateful fixed-effect terms match lme4 names, estimates, predictions", {
  d <- md_data()
  forms <- list(
    y ~ poly(x, 2) * f + (1 | g),
    y ~ scale(x) + factor(a) + (1 | g),
    y ~ splines::ns(x, df = 3) + (1 | g),
    y ~ (x + z + f)^2 + (1 | g),
    y ~ relevel(f, "b") + as.numeric(t) + (1 | g),
    y ~ x * (z + f) + (1 | g),
    y ~ I(x > 0) + log(t) + (1 | g)
  )
  nd <- d[c(2, 30, 77, 140), ]
  for (fm in forms) {
    m <- q(lmm(fm, d, control = quiet))
    r <- lmer_q(fm, d)
    expect_identical(names(fixef(m)), names(lme4::fixef(r)), info = deparse1(fm))
    expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-4, info = deparse1(fm))
    expect_equal(as.numeric(logLik(m)), as.numeric(logLik(r)), tolerance = 1e-6)
    expect_equal(predict(m, nd), predict(r, nd), tolerance = 1e-4,
                 info = deparse1(fm))
    expect_equal(predict(m, nd, re.form = NA), predict(r, nd, re.form = NA),
                 tolerance = 1e-4, info = deparse1(fm))
  }
  # cut() has no predvars (as in lm, new data are re-cut), so only the fit
  # is compared
  m <- q(lmm(y ~ cut(z, 3) + (1 | g), d, control = quiet))
  r <- lmer_q(y ~ cut(z, 3) + (1 | g), d)
  expect_identical(names(fixef(m)), names(lme4::fixef(r)))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-4)
  # %in% / nesting in the fixed part
  m <- q(lmm(y ~ f %in% a + (1 | g), d, control = quiet))
  r <- lmer_q(y ~ f %in% a + (1 | g), d)
  expect_equal(as.numeric(logLik(m)), as.numeric(logLik(r)), tolerance = 1e-6)
  # the formula and terms the user wrote are what the fit reports
  m <- q(lmm(y ~ poly(x, 2) + factor(a) + (1 | g), d, REML = FALSE,
             control = quiet))
  r <- lmer_q(y ~ poly(x, 2) + factor(a) + (1 | g), d, REML = FALSE)
  expect_identical(formula(m), y ~ poly(x, 2) + factor(a) + (1 | g))
  expect_identical(attr(terms(m), "term.labels"), c("poly(x, 2)", "factor(a)"))
  expect_identical(colnames(model.matrix(m)), names(lme4::fixef(r)))
  dr <- drop1(m)$table
  expect_identical(dr$dropped, c("poly(x, 2)", "factor(a)"))
  expect_equal(dr$AIC, stats::drop1(r)$AIC[-1L], tolerance = 1e-4)
  # internal refits reuse the training basis exactly
  m_reml <- q(lmm(y ~ splines::ns(x, df = 3) + (1 | g), d, subset = x < 1.5,
                  control = quiet))
  m_ml <- q(lmm(y ~ splines::ns(x, df = 3) + (1 | g), d, subset = x < 1.5,
                REML = FALSE, control = quiet))
  expect_equal(as.numeric(logLik(update(m_reml, REML = FALSE))),
               as.numeric(logLik(m_ml)), tolerance = 1e-6)
})

test_that("stateful terms work for glmm() and emmeans", {
  d <- md_data()
  m <- q(glmm(cnt ~ poly(x, 2) + factor(a) + (1 | g), d, family = poisson,
              method = "joint_laplace", control = quiet))
  r <- q(lme4::glmer(cnt ~ poly(x, 2) + factor(a) + (1 | g), d,
                     family = poisson))
  expect_identical(names(fixef(m)), names(lme4::fixef(r)))
  expect_equal(fixef(m), lme4::fixef(r), tolerance = 1e-3)
  nd <- d[1:3, ]
  expect_equal(predict(m, nd, type = "link"), predict(r, nd, type = "link"),
               tolerance = 1e-3)
  skip_if_not_installed("emmeans")
  ml <- q(lmm(y ~ poly(x, 2) + f + (1 | g), d, control = quiet))
  rl <- lmer_q(y ~ poly(x, 2) + f + (1 | g), d)
  expect_equal(as.data.frame(emmeans::emmeans(ml, ~ f))$emmean,
               as.data.frame(emmeans::emmeans(rl, ~ f))$emmean,
               tolerance = 1e-4)
})

test_that("random-effect slopes on transformed terms are refused clearly", {
  d <- md_data()
  expect_error(lmm(y ~ x + (poly(x, 2) | g), d, control = quiet),
               class = "mm_condition")
})

test_that("LMM offsets (argument and formula) match lmer", {
  d <- md_data()
  m1 <- q(lmm(y ~ x + offset(2 * z) + (1 | g), d, control = quiet))
  r1 <- lmer_q(y ~ x + offset(2 * z) + (1 | g), d)
  m2 <- q(lmm(y ~ x + (1 | g), d, offset = 2 * z, control = quiet))
  r2 <- lmer_q(y ~ x + (1 | g), d, offset = 2 * z)
  for (p in list(list(m1, r1), list(m2, r2))) {
    expect_equal(fixef(p[[1]]), lme4::fixef(p[[2]]), tolerance = 1e-4)
    expect_equal(as.numeric(logLik(p[[1]])), as.numeric(logLik(p[[2]])),
                 tolerance = 1e-6)
    expect_equal(fitted(p[[1]]), fitted(p[[2]]), tolerance = 1e-4)
    expect_equal(residuals(p[[1]]), residuals(p[[2]]), tolerance = 1e-4)
    expect_equal(predict(p[[1]], re.form = NA), predict(p[[2]], re.form = NA),
                 tolerance = 1e-4)
  }
  nd <- d[1:4, ]
  expect_equal(predict(m1, nd), predict(r1, nd), tolerance = 1e-4)
  expect_equal(predict(m1, nd, re.form = NA), predict(r1, nd, re.form = NA),
               tolerance = 1e-4)
  # an offset= argument cannot be reconstructed for newdata: refuse unless
  # it is supplied again (lme4 silently uses 0 here)
  expect_error(predict(m2, nd), class = "mm_inference_unavailable")
  expect_equal(predict(m2, nd, offset = 2 * z),
               predict(r2, nd) + 2 * nd$z, tolerance = 1e-4)
  expect_error(predict(m1, nd, offset = nd$z), class = "mm_arg_error")
  # REML -> ML refit keeps the offset
  ml <- q(lmm(y ~ x + (1 | g), d, offset = 2 * z, REML = FALSE,
              control = quiet))
  expect_equal(as.numeric(logLik(update(m2, REML = FALSE))),
               as.numeric(logLik(ml)), tolerance = 1e-6)
})

test_that("GLMM predict(newdata) uses formula offsets from newdata", {
  d <- md_data()
  m <- q(glmm(cnt ~ x + offset(log(t)) + (1 | g), d, family = poisson,
              method = "joint_laplace", control = quiet))
  r <- q(lme4::glmer(cnt ~ x + offset(log(t)) + (1 | g), d, family = poisson))
  nd <- d[c(1, 20, 33), ]
  expect_equal(predict(m, nd, type = "link"), predict(r, nd, type = "link"),
               tolerance = 1e-3)
  expect_equal(predict(m, nd, re.form = NA, type = "response"),
               predict(r, nd, re.form = NA, type = "response"),
               tolerance = 1e-3)
  se <- predict(m, nd, type = "link", se.fit = TRUE)
  expect_length(se$se.fit, 3L)
  expect_error(predict(m, nd, se.fit = TRUE), class = "mm_inference_unavailable")
  # argument offsets must be resupplied for newdata
  m2 <- q(glmm(cnt ~ x + (1 | g), d, family = poisson, offset = log(t),
               method = "joint_laplace", control = quiet))
  expect_error(predict(m2, nd), class = "mm_inference_unavailable")
  expect_equal(predict(m2, nd, offset = log(t), type = "link"),
               predict(r, nd, type = "link"), tolerance = 1e-3)
})

test_that("partial re.form keeps only the named terms' BLUPs (lme4)", {
  d <- md_data()
  d$h <- factor(rep(1:8, length.out = nrow(d)))
  m <- q(lmm(y ~ x + (1 + x | g) + (1 | h), d, control = quiet))
  r <- lmer_q(y ~ x + (1 + x | g) + (1 | h), d)
  nd <- d[c(1, 25, 99), ]
  for (rf in list(~ (1 | g), ~ (1 | h), ~ (1 + x | g))) {
    expect_equal(predict(m, nd, re.form = rf), predict(r, nd, re.form = rf),
                 tolerance = 1e-3)
    expect_equal(predict(m, re.form = rf), predict(r, re.form = rf),
                 tolerance = 1e-3)
  }
  expect_equal(predict(m, nd, re.form = ~ (1 + x | g) + (1 | h)),
               predict(m, nd), tolerance = 1e-6)
  expect_error(predict(m, nd, re.form = ~ (1 | g), se.fit = TRUE),
               class = "mm_inference_unavailable")
  g <- q(glmm(cnt ~ x + (1 + x | g), d, family = poisson,
              method = "joint_laplace", control = quiet))
  gr <- q(lme4::glmer(cnt ~ x + (1 + x | g), d, family = poisson))
  expect_equal(predict(g, nd, re.form = ~ (1 | g), type = "link"),
               predict(gr, nd, re.form = ~ (1 | g), type = "link"),
               tolerance = 2e-3)
})

test_that("single-model anova tests a poly() term as one multi-df row, like lmerTest", {
  skip_if_not_installed("lmerTest")
  set.seed(93)
  d <- data.frame(x = runif(150, -2, 2), z = rnorm(150),
                  g = factor(rep(1:15, each = 10)))
  d$y <- 1 + d$x - 0.5 * d$x^2 + 0.3 * d$z + rnorm(15)[d$g] + rnorm(150)
  fit <- lmm(y ~ poly(x, 2) + z + (1 | g), d, control = mm_control(verbose = -1))
  ref <- lmerTest::lmer(y ~ poly(x, 2) + z + (1 | g), d)
  a <- anova(fit)
  r <- stats::anova(ref)
  expect_identical(rownames(a), rownames(r))
  expect_equal(a[["NumDF"]], r[["NumDF"]])
  expect_equal(a[["F value"]], r[["F value"]], tolerance = 1e-4)
  expect_equal(a[["DenDF"]], r[["DenDF"]], tolerance = 1e-3)
})

test_that("explain_model() and audit() show the user's terms, not synthetic columns", {
  set.seed(94)
  d <- data.frame(x = runif(60), f = sample(c("a", "b"), 60, TRUE),
                  g = factor(rep(1:6, 10)), y = rnorm(60))
  fit <- lmm(y ~ poly(x, 2) + factor(f) + (1 | g), d,
             control = mm_control(verbose = -1))
  txt <- c(format(explain_model(fit)), audit(fit)$text, audit(fit)$summary_text)
  expect_false(any(grepl(".poly_x_2", txt, fixed = TRUE)))
  expect_false(any(grepl(".factor_f", txt, fixed = TRUE)))
  msgs <- testthat::capture_messages(
    lmm(y ~ poly(x, 2) + (1 | g), d, control = mm_control(verbose = 0))
  )
  expect_false(any(grepl(".poly_x_2", msgs, fixed = TRUE)))
})
