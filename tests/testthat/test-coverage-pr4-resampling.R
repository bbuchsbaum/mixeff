# PR-4 coverage: mm_bootmer() validation/failure paths and printing, rePCA()
# dispatch and GLMM scaling, mm_lmlist() and mm_allfit() edge cases.

cp4r_data <- function(seed = 441L) {
  set.seed(seed)
  g <- gl(8, 6)
  x <- rep(0:5, 8)
  y <- 1 + 0.5 * x + rnorm(8)[g] + rnorm(8, sd = 0.3)[g] * x + rnorm(48)
  data.frame(y = y, x = x, g = g, w = rep(c(1, 2), 24))
}

cp4r_ctrl <- function() mm_control(verbose = -1)

test_that("mm_bootmer() validates its arguments", {
  d <- cp4r_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4r_ctrl())
  expect_error(mm_bootmer(list(), fixef), class = "mm_arg_error")
  expect_error(mm_bootmer(fit, fixef, nsim = 0), class = "mm_arg_error")
  expect_error(mm_bootmer(fit, fixef, nsim = 1.5), class = "mm_arg_error")
  expect_error(mm_bootmer(fit, fixef, use.u = NA), class = "mm_arg_error")
  expect_error(mm_bootmer(fit, function(f) "a"), class = "mm_arg_error")
  expect_error(mm_bootmer(fit, function(f) numeric()), class = "mm_arg_error")
})

test_that("mm_bootmer() records failed replicates and prints", {
  d <- cp4r_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4r_ctrl())
  calls <- 0L
  # Original fit gives two values; every refit gives one -> all fail.
  stat <- function(f) {
    calls <<- calls + 1L
    if (calls == 1L) fixef(f) else fixef(f)[1]
  }
  expect_message(
    b <- mm_bootmer(fit, stat, nsim = 2L, seed = 1, verbose = TRUE),
    class = "mm_bootstrap_failure_notice"
  ) |> expect_output("\\.\\.")
  expect_identical(attr(b, "bootFail"), 2L)
  expect_true(all(is.na(b$t)))
  expect_match(names(attr(b, "boot.fail.msgs")), "different length")
  expect_output(print(b), "failed replicates: 2")

  ok <- mm_bootmer(fit, fixef, nsim = 3L, seed = 2, use.u = TRUE)
  expect_identical(dim(ok$t), c(3L, 2L))
  expect_true(all(is.finite(ok$t)))
  out <- capture.output(print(ok))
  expect_match(out[[1]], "3 replicates, conditional on u")
  expect_true(any(grepl("bias", out)))
  expect_false(any(grepl("failed replicates", out)))
})

test_that("rePCA() dispatch and GLMM scaling", {
  skip_if_not_installed("lme4")
  expect_error(rePCA(1), class = "mm_arg_error")
  d <- cp4r_data()
  ref <- lme4::lmer(y ~ x + (x | g), d)
  fwd <- rePCA(ref)
  expect_equal(fwd$g$sdev, lme4::rePCA(ref)$g$sdev)

  set.seed(442)
  gd <- data.frame(g = gl(12, 8), x = rnorm(96))
  gd$y <- rpois(96, exp(0.3 + 0.3 * gd$x + rnorm(12, sd = 0.5)[gd$g]))
  gf <- glmm(y ~ x + (1 | g), gd, family = poisson(), control = cp4r_ctrl())
  pc <- rePCA(gf)
  expect_s3_class(pc, "prcomplist")
  # Poisson has no residual scale: rePCA is on the theta scale itself.
  expect_equal(pc$g$sdev, unname(gf$theta), tolerance = 1e-6)
  gref <- lme4::glmer(y ~ x + (1 | g), gd, family = poisson())
  expect_equal(pc$g$sdev, lme4::rePCA(gref)$g$sdev, tolerance = 1e-3)
  expect_type(summary(pc), "list")

  mu <- exp(1 + 0.2 * gd$x + rnorm(12, sd = 0.3)[gd$g])
  gd$yg <- rgamma(96, shape = 5, rate = 5 / mu)
  gg <- glmm(yg ~ x + (1 | g), gd, family = Gamma(link = "log"),
             control = cp4r_ctrl())
  blocks <- mixeff:::mm_relative_covariance_blocks(gg)
  sd_g <- attr(VarCorr(gg)$g, "stddev")
  expect_equal(sqrt(blocks$g[1, 1]), unname(sd_g / sigma(gg)),
               tolerance = 1e-6)
  # A non-finite scale falls back to the unscaled covariance.
  bad <- gg
  bad$sigma <- NA_real_
  expect_equal(mixeff:::mm_relative_covariance_blocks(bad)$g[1, 1],
               unname(sd_g)^2, tolerance = 1e-6)
  # An empty basis is the intercept.
  nb <- gg
  nb$artifact$semantic_model$random_terms[[1]]$basis <- list()
  expect_identical(rownames(mixeff:::mm_relative_covariance_blocks(nb)$g),
                   "(Intercept)")
})

test_that("mm_lmlist() matches lme4::lmList and handles options", {
  d <- cp4r_data()
  expect_error(mm_lmlist(y ~ x, d), class = "mm_formula_error")
  expect_error(mm_lmlist(y ~ x | g, as.list(d)), class = "mm_data_error")
  ll <- mm_lmlist(y ~ x | g, d)
  out <- capture.output(print(ll))
  expect_match(out[[1]], "Per-group fits")
  expect_true(any(grepl("pooled", out)))
  ci <- confint(ll, parm = "x")
  expect_identical(dim(ci), c(8L, 2L, 1L))
  ci2 <- confint(ll, parm = 2)
  expect_equal(unclass(ci), unclass(ci2))
  expect_output(print(ci), "x")
  skip_if_not_installed("lme4")
  ref <- lme4::lmList(y ~ x | g, d)
  expect_equal(as.matrix(coef(ll)), as.matrix(coef(ref)), ignore_attr = TRUE)
  pooled <- sqrt(sum(vapply(ref, function(f) sum(resid(f)^2), 0)) /
                   sum(vapply(ref, df.residual, 0)))
  expect_equal(sigma(ll), pooled)
  # Pooled intervals use the pooled residual SD and degrees of freedom.
  f1 <- ref[[1]]
  se <- sqrt(diag(vcov(f1))) / sigma(f1) * pooled
  q <- qt(0.975, sum(vapply(ref, df.residual, 0)))
  expect_equal(unname(unclass(confint(ll))[1, , "x"]),
               unname(coef(f1)[["x"]] + c(-1, 1) * q * se[["x"]]))
  # Unpooled intervals are each group's own lm() intervals.
  expect_equal(unname(unclass(confint(ll, pool = FALSE))[1, , "x"]),
               unname(confint(f1)["x", ]))

  # weights, subset, family-as-string, unpooled confint
  lw <- mm_lmlist(y ~ x | g, d, weights = w, subset = x > 0)
  f1 <- lm(y ~ x, d[d$g == "1" & d$x > 0, ], weights = w)
  expect_equal(unname(unlist(coef(lw)[1, ])), unname(coef(f1)))
  set.seed(443)
  d$cnt <- rpois(48, exp(0.5 + 0.1 * d$x))
  lp <- mm_lmlist(cnt ~ x | g, d, family = "poisson")
  expect_false(attr(lp, "pool"))
  expect_identical(attr(lp, "family")$family, "poisson")
  cip <- confint(lp)
  fp <- glm(cnt ~ x, poisson, d[d$g == "2", ])
  expect_equal(unclass(cip)["2", , "x"],
               unname(confint.default(fp)["x", ]), ignore_attr = TRUE)
  expect_output(print(lp), "Coefficients")

  # A group whose fit fails is reported and kept as NULL.
  d$z <- d$x
  d$z[d$g == "3"] <- NA
  expect_message(lf <- mm_lmlist(y ~ z | g, d, na.action = na.fail),
                 class = "mm_lmlist_failure_notice")
  expect_null(lf[["3"]])
  expect_true(all(is.na(unlist(coef(lf)["3", ]))))
  expect_true(all(is.na(unclass(confint(lf))["3", , ])))
  expect_length(attr(lf, "failed"), 1L)
})

test_that("mm_allfit() defaults, validation and an all-failed summary", {
  expect_identical(mixeff:::mm_allfit_optimizers()[[1]], "auto")
  expect_true(all(c("bobyqa", "pattern_search", "prima_newuoa") %in%
                    mixeff:::mm_allfit_optimizers()))
  expect_error(mm_allfit(list()), class = "mm_arg_error")
  d <- cp4r_data()
  fit <- lmm(y ~ x + (1 | g), d, control = cp4r_ctrl())
  expect_error(mm_allfit(fit, optimizers = character()), class = "mm_arg_error")
  expect_output(af <- mm_allfit(fit, optimizers = c("auto", "pattern_search"),
                                verbose = TRUE), "pattern_search")
  s <- summary(af)
  expect_true(all(s$which.OK))
  expect_equal(unname(s$llik[[1]]), unname(s$llik[[2]]), tolerance = 1e-5)
  dead <- structure(list(a = simpleError("boom")), class = "mm_allfit",
                    times = c(a = 0))
  sd <- summary(dead)
  expect_null(sd$fixef)
  expect_null(sd$theta)
  expect_identical(sd$msgs$a, "boom")
  out <- capture.output(print(dead))
  expect_match(out[[1]], "0 of 1")
  expect_true(any(grepl("a: boom", out)))
})
