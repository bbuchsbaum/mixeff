# hatvalues(), cooks.distance(), influence() and diagnostic plots
# (pre-CRAN audit 4.5).

test_that("hatvalues and cooks.distance match lme4 on LMMs", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  m <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  h <- hatvalues(m)
  expect_named(h, rownames(sleepstudy))
  expect_equal(unname(h), unname(stats::hatvalues(r)), tolerance = 1e-5)
  expect_equal(unname(cooks.distance(m)), unname(stats::cooks.distance(r)),
               tolerance = 1e-4)
  H <- hatvalues(m, fullHatMatrix = TRUE)
  expect_equal(dim(H), c(180L, 180L))
  expect_equal(diag(H), h, ignore_attr = TRUE)
  # H maps y to the conditional fitted values.
  expect_equal(as.numeric(H %*% sleepstudy$Reaction), unname(fitted(m)),
               tolerance = 1e-4)

  # Crossed terms that the engine reorders, and prior weights.
  data("Penicillin", package = "lme4", envir = environment())
  m <- lmm(diameter ~ 1 + (1 | sample) + (1 | plate), Penicillin,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(diameter ~ 1 + (1 | sample) + (1 | plate), Penicillin)
  expect_equal(unname(hatvalues(m)), unname(stats::hatvalues(r)),
               tolerance = 1e-5)
  expect_equal(unname(cooks.distance(m)), unname(stats::cooks.distance(r)),
               tolerance = 1e-4)

  set.seed(2)
  sleepstudy$w <- runif(180, 0.5, 2)
  m <- lmm(Reaction ~ Days + (Days || Subject), sleepstudy, weights = w,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(Reaction ~ Days + (Days || Subject), sleepstudy, weights = w)
  expect_equal(unname(hatvalues(m)), unname(stats::hatvalues(r)),
               tolerance = 1e-5)
  expect_equal(unname(cooks.distance(m)), unname(stats::cooks.distance(r)),
               tolerance = 1e-4)
})

test_that("GLMM hat values use the working weights and stay in [0, 1]", {
  skip_if_not_installed("lme4")
  data("cbpp", package = "lme4", envir = environment())
  f <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp,
            family = binomial(), control = mm_control(verbose = -1))
  h <- hatvalues(f)
  expect_length(h, nrow(cbpp))
  expect_true(all(h > 0 & h < 1))
  expect_lt(sum(h), 4 + 15)
  expect_gt(sum(h), 4)
  cd <- cooks.distance(f)
  expect_true(all(is.finite(cd) & cd >= 0))
})

test_that("influence() deletion refits match lme4", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  m <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  im <- influence(m, groups = "Subject")
  ir <- stats::influence(r, groups = "Subject")
  expect_s3_class(im, "mm_influence")
  expect_identical(names(im)[1:9], names(ir)[1:9])
  expect_identical(colnames(im[["var.cov.comps[-Subject]"]]),
                   colnames(ir[["var.cov.comps[-Subject]"]]))
  expect_true(all(im$converged))
  expect_equal(unname(dfbeta(im)), unname(stats::dfbeta(ir)), tolerance = 1e-4)
  expect_equal(unname(dfbetas(im)), unname(stats::dfbetas(ir)),
               tolerance = 1e-3)
  expect_equal(unname(cooks.distance(im)), unname(stats::cooks.distance(ir)),
               tolerance = 1e-3)
  expect_equal(unname(dfbeta(im, "var.cov")),
               unname(stats::dfbeta(ir, "var.cov")), tolerance = 1e-3)
  expect_output(print(im), "deletion of 18 Subject")
  expect_s3_class(influence(m, do.coef = FALSE), "mm_influence")
})

test_that("influence() handles case deletion, GLMMs and bad input", {
  set.seed(4)
  d <- data.frame(y = rpois(30, 3), x = rnorm(30), g = gl(6, 5))
  f <- glmm(y ~ x + (1 | g), d, family = poisson(),
            control = mm_control(verbose = -1))
  inf <- influence(f)
  expect_identical(inf$groups, "case")
  expect_equal(nrow(dfbeta(inf)), 30L)
  expect_identical(names(inf$converged), rownames(d))
  expect_error(influence(f, groups = "nope"), class = "mm_arg_error")
  expect_error(influence(f, groups = "g", data = d[1:3, ]),
               class = "mm_arg_error")
  expect_error(influence(f, maxfun = 10), class = "mm_arg_error")
  inf_g <- influence(f, groups = "g", data = d)
  expect_identical(inf_g$deleted, levels(d$g))
})

test_that("diagnostic plots draw for LMM, GLMM and ranef objects", {
  set.seed(5)
  d <- data.frame(y = rnorm(40), yc = rpois(40, 2), x = rnorm(40),
                  g = gl(8, 5))
  m <- lmm(y ~ x + (x | g), d, control = mm_control(verbose = -1))
  f <- glmm(yc ~ x + (1 | g), d, family = poisson(),
            control = mm_control(verbose = -1))
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_invisible(plot(m))
  expect_identical(plot(m, which = 1:3, ask = FALSE), m)
  expect_identical(plot(f), f)
  q <- qqnorm(m)
  expect_length(q$x, 40L)
  expect_error(plot(m, which = 7), class = "mm_arg_error")
  re <- ranef(m, condVar = TRUE)
  expect_identical(plot(re), re)
  expect_identical(qqnorm(re), re)
  expect_identical(plot(ranef(f, condVar = TRUE)), ranef(f, condVar = TRUE))
  expect_error(plot(re, level = 2), class = "mm_arg_error")

  skip_if_not_installed("lattice")
  dp <- lattice::dotplot(re)
  expect_named(dp, "g")
  expect_s3_class(dp$g, "trellis")
  qp <- lattice::qqmath(re)
  expect_s3_class(qp$g, "trellis")
  expect_s3_class(lattice::dotplot(ranef(f))$g, "trellis")
})
