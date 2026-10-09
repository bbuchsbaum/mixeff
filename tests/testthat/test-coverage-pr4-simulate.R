# PR-4 coverage: GLMM simulation helpers (target resolution, PSD factor,
# term layout fallbacks, cbind response naming) and refit() newresp checks.

cp4sim_ctrl <- function() mm_control(verbose = -1)

cp4sim_binom <- function(seed = 471L) {
  set.seed(seed)
  g <- gl(10, 6)
  x <- rnorm(60)
  n <- rep(c(5L, 8L), 30)
  k <- rbinom(60, n, plogis(0.2 + 0.4 * x + rnorm(10, sd = 0.5)[g]))
  data.frame(k = k, n = n, x = x, g = g)
}

test_that("simulation target arguments are validated", {
  d <- cp4sim_binom()
  fit <- glmm(cbind(k, n - k) ~ x + (1 | g), d, family = binomial(),
              control = cp4sim_ctrl())
  expect_error(simulate(fit, use.u = NA), class = "mm_arg_error")
  expect_error(simulate(fit, use.u = c(TRUE, FALSE)), class = "mm_arg_error")
  rp <- mixeff:::mm_simulate_resolve_partial
  expect_identical(rp(fit, "partial", NA), list(target = "unconditional"))
  expect_identical(rp(fit, "partial", NULL), list(target = "conditional"))
  expect_identical(rp(fit, "conditional", NULL), list(target = "conditional"))
})

test_that("in-sample simulation drops na.exclude padding", {
  set.seed(472)
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5))
  d$x[c(3, 9)] <- NA
  fit <- lmm(y ~ x + (1 | g), d, na.action = na.exclude,
             control = cp4sim_ctrl())
  padded <- fitted(fit)
  expect_length(padded, 40L)
  ins <- mixeff:::mm_simulate_insample(fit, padded)
  expect_length(ins, 38L)
  expect_false(anyNA(ins))
  expect_identical(mixeff:::mm_simulate_insample(fit, 1:38), as.numeric(1:38))
})

test_that("the PSD lower factor handles empty and singular blocks", {
  ch <- mixeff:::mm_chol_lower_psd
  expect_identical(dim(ch(matrix(0, 0, 0))), c(0L, 0L))
  S <- matrix(c(4, 2, 2, 1), 2)  # rank one
  L <- ch(S)
  expect_equal(L %*% t(L), S)
  expect_identical(L[2, 2], 0)
  full <- matrix(c(2, 0.5, 0.5, 1), 2)
  expect_equal(ch(full), t(chol(full)), ignore_attr = TRUE)
})

test_that("term layout and new random effects fall back sensibly", {
  d <- cp4sim_binom()
  fit <- glmm(cbind(k, n - k) ~ x + (1 | g), d, family = binomial(),
              control = cp4sim_ctrl())
  z <- mixeff:::mm_simulate_new_re(fit, 3L, exclude = "g")
  expect_identical(dim(z), c(60L, 3L))
  expect_true(all(z == 0))
  nb <- fit
  nb$artifact$semantic_model$random_terms[[1]]$basis <- list()
  lay <- mixeff:::mm_simulate_term_layout(nb)
  expect_length(lay[[1]]$basis, 1L)
  expect_equal(lay[[1]]$basis[[1]], rep(1, 60))
  expect_equal(lay[[1]]$L[1, 1], unname(attr(VarCorr(fit)$g, "stddev")),
               tolerance = 1e-8)
})

test_that("cbind response names follow the call, even via a formula object", {
  d <- cp4sim_binom()
  fo <- cbind(k, n - k) ~ x + (1 | g)
  fit <- glmm(fo, d, family = binomial(), control = cp4sim_ctrl())
  expect_identical(mixeff:::mm_glmm_cbind_colnames(fit), c("k", ""))
  sims <- simulate(fit, nsim = 1, seed = 3)
  expect_identical(colnames(sims[[1]]), c("k", ""))
  expect_true(all(rowSums(sims[[1]]) == d$n))
  lost <- fit
  lost$call$formula <- quote(no_such_formula_object)
  expect_identical(mixeff:::mm_glmm_cbind_colnames(lost), c("", ""))
  skip_if_not_installed("lme4")
  ref <- lme4::glmer(fo, d, family = binomial())
  expect_identical(colnames(simulate(ref, seed = 3)[[1]]), c("k", ""))
})

test_that("refit() checks binomial newresp shapes", {
  d <- cp4sim_binom()
  fit <- glmm(cbind(k, n - k) ~ x + (1 | g), d, family = binomial(),
              control = cp4sim_ctrl())
  bad <- cbind(d$k, d$n - d$k)
  bad[1, ] <- 0
  expect_error(refit(fit, bad), class = "mm_arg_error")
  expect_error(refit(fit, d$k / d$n * 2), class = "mm_arg_error")
  ok <- refit(fit, cbind(d$k, d$n - d$k))
  expect_equal(fixef(ok), fixef(fit), tolerance = 1e-5)

  d$b <- as.integer(d$k > d$n / 2)
  fb <- glmm(b ~ x + (1 | g), d, family = binomial(), control = cp4sim_ctrl())
  rl <- refit(fb, d$b == 1L)
  expect_equal(fixef(rl), fixef(fb), tolerance = 1e-5)
})

test_that("refit() forwards foreign fits and refuses other objects", {
  expect_error(refit(1), class = "mm_arg_error")
  expect_error(refit(1, 2), class = "mm_arg_error")
  skip_if_not_installed("lme4")
  set.seed(473)
  d <- data.frame(y = rnorm(40), x = rnorm(40), g = gl(8, 5))
  d$y <- d$y + rnorm(8)[d$g]
  ref <- lme4::lmer(y ~ x + (1 | g), d)
  same <- refit(ref)
  expect_equal(lme4::fixef(same), lme4::fixef(ref))
  y2 <- d$y + 1
  shifted <- refit(ref, y2)
  expect_equal(lme4::fixef(shifted)[["(Intercept)"]],
               lme4::fixef(ref)[["(Intercept)"]] + 1, tolerance = 1e-5)
})

test_that("the simulate seed attribute initialises an unset RNG", {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  if (had) rm(".Random.seed", envir = globalenv())
  st <- mixeff:::mm_rng_state(NULL)
  expect_true(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_identical(st, get(".Random.seed", envir = globalenv()))
  seeded <- mixeff:::mm_rng_state(5L)
  expect_identical(as.integer(seeded), 5L)
  expect_identical(attr(seeded, "kind"), as.list(RNGkind()))
})
