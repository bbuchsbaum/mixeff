# simulate()/refit() for GLMMs and lme4-consistent re.form semantics
# (pre-CRAN audit 4.4).

skip_if_no_lme4_sim <- function() {
  testthat::skip_if_not_installed("lme4")
}

test_that("simulate.mm_lmm matches lme4 draws for re.form = NA and NULL", {
  skip_if_no_lme4_sim()
  data("sleepstudy", package = "lme4", envir = environment())
  m <- lmm(Reaction ~ Days + (Days | Subject), sleepstudy,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  # default re.form = NA: new random effects, same RNG stream as lme4
  a <- simulate(m, nsim = 3, seed = 11)
  b <- simulate(r, nsim = 3, seed = 11)
  expect_identical(names(a), names(b))
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)), tolerance = 1e-4)
  expect_identical(attr(a, "mm_re_form"), "unconditional")
  # re.form = NULL / use.u = TRUE: conditional on the BLUPs
  a <- simulate(m, nsim = 2, seed = 11, re.form = NULL)
  b <- simulate(r, nsim = 2, seed = 11, re.form = NULL)
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)), tolerance = 1e-5)
  expect_equal(simulate(m, nsim = 2, seed = 11, use.u = TRUE), a,
               ignore_attr = TRUE)
  expect_error(simulate(m, use.u = TRUE, re.form = NULL),
               class = "mm_arg_error")
  # ~0 is the same as NA
  expect_equal(simulate(m, seed = 4, re.form = ~0),
               simulate(m, seed = 4), ignore_attr = TRUE)
})

test_that("simulate.mm_lmm follows lme4's term order for crossed effects", {
  skip_if_no_lme4_sim()
  data("Penicillin", package = "lme4", envir = environment())
  m <- lmm(diameter ~ 1 + (1 | sample) + (1 | plate), Penicillin,
           control = mm_control(verbose = -1))
  r <- lme4::lmer(diameter ~ 1 + (1 | sample) + (1 | plate), Penicillin)
  expect_equal(unname(as.matrix(simulate(m, 2, seed = 3))),
               unname(as.matrix(simulate(r, 2, seed = 3))),
               tolerance = 1e-4)
})

test_that("binomial cbind() GLMM simulate returns lme4's matrix shape and draws", {
  skip_if_no_lme4_sim()
  data("cbpp", package = "lme4", envir = environment())
  f <- glmm(cbind(incidence, size - incidence) ~ period + (1 | herd), cbpp,
            family = binomial(), method = "joint_laplace",
            control = mm_control(verbose = -1))
  g <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                   cbpp, family = binomial())
  a <- simulate(f, nsim = 2, seed = 5)
  b <- simulate(g, nsim = 2, seed = 5)
  expect_s3_class(a, "data.frame")
  expect_true(is.matrix(a[[1]]) && ncol(a[[1]]) == 2L)
  expect_identical(colnames(a[[1]]), colnames(b[[1]]))
  expect_identical(unname(a[[1]]), unname(b[[1]]))
  expect_identical(unname(a[[2]]), unname(b[[2]]))
  expect_equal(rowSums(a[[1]]), cbpp$size)
  # refit() accepts the matrix and a one-column data frame
  f2 <- refit(f, a[1])
  g2 <- lme4::refit(g, b[[1]])
  expect_s3_class(f2, "mm_glmm")
  expect_equal(unname(fixef(f2)), unname(lme4::fixef(g2)), tolerance = 1e-4)
})

test_that("weighted proportion binomial and Poisson simulate match lme4", {
  skip_if_no_lme4_sim()
  data("cbpp", package = "lme4", envir = environment())
  cbpp$prop <- cbpp$incidence / cbpp$size
  f <- glmm(prop ~ period + (1 | herd), cbpp, family = binomial(),
            weights = size, method = "joint_laplace",
            control = mm_control(verbose = -1))
  g <- lme4::glmer(prop ~ period + (1 | herd), cbpp, family = binomial(),
                   weights = size)
  a <- simulate(f, nsim = 2, seed = 5)
  b <- simulate(g, nsim = 2, seed = 5)
  expect_equal(unname(as.matrix(a)), unname(as.matrix(b)))
  f2 <- refit(f, a[[1]])
  expect_equal(unname(fixef(f2)), unname(lme4::fixef(lme4::refit(g, b[[1]]))),
               tolerance = 1e-4)

  data("grouseticks", package = "lme4", envir = environment())
  f <- glmm(TICKS ~ YEAR + (1 | BROOD), grouseticks, family = poisson(),
            method = "joint_laplace", control = mm_control(verbose = -1))
  g <- lme4::glmer(TICKS ~ YEAR + (1 | BROOD), grouseticks, family = poisson())
  expect_equal(unname(as.matrix(simulate(f, 2, seed = 5))),
               unname(as.matrix(simulate(g, 2, seed = 5))))
  a <- simulate(f, 2, seed = 5, use.u = TRUE)
  expect_equal(unname(as.matrix(a)),
               unname(as.matrix(simulate(g, 2, seed = 5, use.u = TRUE))))
  expect_identical(attr(a, "mm_re_form"), "conditional")
  expect_identical(attr(a, "mm_method"), "r_side_glmm_parametric")
  expect_true(all(a[[1]] >= 0 & a[[1]] %% 1 == 0))
  f3 <- refit(f, a[[1]])
  expect_s3_class(f3, "mm_glmm")
  expect_identical(f3$method, f$method)
})

test_that("Bernoulli, Gamma and negative-binomial simulate are family-aware", {
  set.seed(9)
  g <- factor(rep(1:15, each = 8))
  x <- rnorm(120)
  u <- rep(rnorm(15, sd = 0.3), each = 8)
  mu <- exp(0.5 + 0.4 * x + u)
  d <- data.frame(
    yb = rbinom(120, 1, plogis(0.2 + x + u)),
    yg = rgamma(120, shape = 5, rate = 5 / mu),
    yn = rnbinom(120, mu = mu, size = 2),
    x = x, g = g
  )
  fb <- glmm(yb ~ x + (1 | g), d, family = binomial(),
             control = mm_control(verbose = -1))
  sb <- simulate(fb, nsim = 3, seed = 1)
  expect_equal(dim(sb), c(120L, 3L))
  expect_true(all(unlist(sb) %in% c(0, 1)))
  expect_identical(simulate(fb, nsim = 3, seed = 1), sb)
  expect_s3_class(refit(fb, sb[[1]]), "mm_glmm")
  expect_s3_class(refit(fb, factor(sb[[1]], labels = c("no", "yes"))),
                  "mm_glmm")

  fg <- glmm(yg ~ x + (1 | g), d, family = Gamma(link = "log"),
             control = mm_control(verbose = -1))
  sg <- simulate(fg, nsim = 400, seed = 2, use.u = TRUE)
  expect_true(all(unlist(sg) > 0))
  # Var(y) = phi * mu^2 with phi = sigma^2
  mu_hat <- fitted(fg)
  emp <- apply(as.matrix(sg), 1, var) / mu_hat^2
  expect_equal(mean(emp), sigma(fg)^2, tolerance = 0.1)

  fn <- glmm(yn ~ x + (1 | g), d, family = mm_negative_binomial(),
             control = mm_control(verbose = -1))
  sn <- simulate(fn, nsim = 400, seed = 3, use.u = TRUE)
  expect_true(all(unlist(sn) %% 1 == 0))
  emp <- (apply(as.matrix(sn), 1, var) - fitted(fn)) / fitted(fn)^2
  expect_equal(mean(emp), 1 / fn$family$nb_theta, tolerance = 0.25)
  rn <- refit(fn, sn[[1]])
  expect_s3_class(rn, "mm_glmm")
  expect_true(is.finite(rn$family$nb_theta))
})

test_that("GLMM simulate/refit refuse what they cannot honour", {
  set.seed(3)
  d <- data.frame(y = rpois(40, 2), x = rnorm(40), g = gl(8, 5))
  f <- glmm(y ~ x + (1 | g), d, family = poisson(),
            control = mm_control(verbose = -1))
  expect_error(simulate(f, newdata = d), class = "mm_arg_error")
  # A formula naming every term is conditional simulation (lme4).
  expect_identical(attr(simulate(f, re.form = ~ (1 | g)), "mm_re_form"),
                   "conditional")
  expect_error(refit(f), class = "mm_arg_error")
  expect_error(refit(f, d[, 1:2]), class = "mm_arg_error")
  expect_error(refit(f, c(1, NA, rep(1, 38))), class = "mm_arg_error")

  d$p <- runif(40)
  fw <- glmm(p ~ x + (1 | g), d, family = binomial(), weights = rep(2.5, 40),
             control = mm_control(verbose = -1))
  expect_error(simulate(fw), class = "mm_inference_unavailable")
})

test_that("refit() reproduces fits with formula or argument offsets (no double counting)", {
  set.seed(91)
  n <- 200
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, each = 10)),
                  expo = runif(n, 0.5, 2))
  d$count <- rpois(n, exp(0.2 + 0.4 * d$x + rnorm(20, sd = 0.3)[d$g] + log(d$expo)))
  ctl <- mm_control(verbose = -1)
  f_formula <- glmm(count ~ x + offset(log(expo)) + (1 | g), d,
                    family = poisson, control = ctl)
  expect_equal(fixef(refit(f_formula, d$count)), fixef(f_formula), tolerance = 1e-6)
  f_arg <- glmm(count ~ x + (1 | g), d, family = poisson,
                offset = log(d$expo), control = ctl)
  expect_equal(fixef(refit(f_arg, d$count)), fixef(f_arg), tolerance = 1e-6)

  d$y <- 1 + 0.5 * d$x + rnorm(20)[d$g] + rnorm(n) + log(d$expo)
  l_formula <- lmm(y ~ x + offset(log(expo)) + (1 | g), d, control = ctl)
  expect_equal(fixef(refit(l_formula, d$y)), fixef(l_formula), tolerance = 1e-6)
  l_arg <- lmm(y ~ x + (1 | g), d, offset = log(d$expo), control = ctl)
  expect_equal(fixef(refit(l_arg, d$y)), fixef(l_arg), tolerance = 1e-6)
})

test_that("simulate() with a partial re.form conditions on the named terms only (lme4 semantics)", {
  set.seed(92)
  d <- expand.grid(subj = factor(1:20), item = factor(1:10), rep = 1:2)
  d$y <- rnorm(20, sd = 2)[d$subj] + rnorm(10, sd = 1)[d$item] + rnorm(nrow(d))
  fit <- lmm(y ~ 1 + (1 | subj) + (1 | item), d,
             control = mm_control(verbose = -1))
  sims <- as.matrix(simulate(fit, nsim = 400, seed = 3, re.form = ~ (1 | subj)))
  cond <- predict(fit, re.form = ~ (1 | subj))
  # Averaging over the redrawn item effects and noise recovers the
  # subject-conditional mean.
  expect_lt(max(abs(rowMeans(sims) - cond)), 0.35)
  # Naming every term is ordinary conditional simulation.
  full <- simulate(fit, nsim = 1, seed = 3, re.form = ~ (1 | subj) + (1 | item))
  expect_identical(attr(full, "mm_re_form"), "conditional")
  # GLMMs take the same route.
  d$count <- rpois(nrow(d), exp(0.5 + rnorm(20, sd = 0.4)[d$subj]))
  g <- glmm(count ~ 1 + (1 | subj) + (1 | item), d, family = poisson,
            control = mm_control(verbose = -1))
  gs <- simulate(g, nsim = 2, seed = 1, re.form = ~ (1 | subj))
  expect_identical(attr(gs, "mm_re_form"), "partial")
  expect_equal(dim(gs), c(nrow(d), 2L))
  expect_error(simulate(fit, re.form = "subj"), class = "mm_arg_error")
})
