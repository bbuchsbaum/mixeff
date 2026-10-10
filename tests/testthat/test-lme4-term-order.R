# Random-effect term order and `||` layout follow lme4's mkReTrms (the
# engine orders terms by decreasing size and keeps `||` as one diagonal
# block; mixeff-rs 2873312 expands `||` with factors like lme4 and codes
# factor bases with model.matrix()).

mm_term_order_data <- function() {
  set.seed(1)
  d <- data.frame(a = factor(rep(1:5, each = 24)), g = factor(rep(1:20, 6)),
                  h = factor(rep(1:5, 24)), x = rnorm(120),
                  f = factor(sample(c("p", "q", "r"), 120, TRUE)))
  d$y <- rnorm(120) + rnorm(20)[d$g] + rnorm(5)[d$a] + d$x * rnorm(20)[d$g]
  d
}

test_that("ranef/VarCorr/getME follow lme4's term order and || layout", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  d <- mm_term_order_data()
  forms <- list(
    y ~ x + (1 | a) + (x || g) + (1 | h),
    y ~ x + (0 + x | g) + (1 | g),
    y ~ x + (1 | h) + (1 | g) + (1 | a),
    y ~ x + (1 + x + f || g),
    y ~ x + (0 + f | g),
    y ~ x + (0 + f + x | g),
    y ~ x + (1 | a / h)
  )
  for (fo in forms) {
    lab <- deparse1(fo)
    fit <- lmm(fo, d, control = mm_control(verbose = -1))
    ref <- suppressMessages(suppressWarnings(lme4::lmer(fo, d)))
    expect_identical(getME(fit, "cnms"), lme4::getME(ref, "cnms"), label = lab)
    expect_identical(names(getME(fit, "theta")), names(lme4::getME(ref, "theta")),
                     label = lab)
    expect_equal(unname(getME(fit, "theta")), unname(lme4::getME(ref, "theta")),
                 tolerance = 2e-2, label = lab)
    vc <- as.data.frame(VarCorr(fit))
    vc_ref <- as.data.frame(lme4::VarCorr(ref))
    expect_identical(vc[, 1:3], vc_ref[, 1:3], label = lab)
    expect_equal(vc$sdcor, vc_ref$sdcor, tolerance = 2e-2, label = lab)
    expect_identical(names(VarCorr(fit)), names(lme4::VarCorr(ref)), label = lab)
    expect_identical(lapply(ranef(fit), names), lapply(lme4::ranef(ref), names),
                     label = lab)
    expect_equal(ngrps(fit), lme4::ngrps(ref), label = lab)
    expect_identical(getME(fit, "Lind"), lme4::getME(ref, "Lind"), label = lab)
    expect_equal(as.matrix(getME(fit, "Zt")), as.matrix(lme4::getME(ref, "Zt")),
                 ignore_attr = TRUE, label = lab)
    expect_identical(attr(logLik(fit), "df"), attr(logLik(ref), "df"), label = lab)
    expect_equal(as.numeric(logLik(fit)), as.numeric(logLik(ref)),
                 tolerance = 1e-5, label = lab)
  }
})

test_that("|| prints lme4's separate VarCorr entries", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  fit <- lmm(Reaction ~ Days + (Days || Subject), sleepstudy,
             control = mm_control(verbose = -1))
  vc <- VarCorr(fit)
  expect_identical(names(vc), c("Subject", "Subject.1"))
  expect_identical(as.data.frame(vc)$grp, c("Subject", "Subject.1", "Residual"))
  out <- paste(utils::capture.output(print(vc)), collapse = "\n")
  expect_false(grepl("Corr", out, fixed = TRUE))
})

test_that("unused grouping levels are dropped consistently", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4", envir = environment())
  d <- sleepstudy
  d$Subject <- factor(d$Subject, levels = c(levels(d$Subject), "999"))
  fit <- lmm(Reaction ~ Days + (Days | Subject), d,
             control = mm_control(verbose = -1))
  expect_identical(unname(ngrps(fit)), 18L)
  expect_identical(nrow(ranef(fit)$Subject), 18L)
  expect_identical(unname(getME(fit, "l_i")), 18L)
  expect_false("999" %in% rownames(ranef(fit)$Subject))
})
