# Engine threading for bootstrap and profile (mixeff-rs 2873312): responses
# are simulated serially in the serial RNG order and refits run on worker
# threads, so every result must be identical for threads = 1 and 2.

mm_threads_lmm_data <- function() {
  set.seed(42)
  d <- data.frame(g = factor(rep(1:12, each = 6)), x = rnorm(72))
  d$y <- 1 + 0.5 * d$x + rnorm(12, sd = 0.8)[d$g] + rnorm(72)
  d
}

test_that("bootstrap_control() validates and carries threads", {
  expect_identical(bootstrap_control()$threads, 1L)
  expect_identical(bootstrap_control(threads = 2)$threads, 2L)
  for (bad in list(0, -1, 1.5, NA, "2", c(1, 2))) {
    expect_error(bootstrap_control(threads = bad), class = "mm_arg_error")
  }
})

test_that("LMM bootstrap confint is identical for threads = 1 and 2", {
  skip_on_cran()
  fit <- lmm(y ~ x + (1 | g), mm_threads_lmm_data(),
             control = mm_control(verbose = -1))
  ci1 <- confint(fit, method = "bootstrap",
                 bootstrap = bootstrap_control(nsim = 40, seed = 3))
  ci2 <- confint(fit, method = "bootstrap",
                 bootstrap = bootstrap_control(nsim = 40, seed = 3, threads = 2))
  expect_identical(unclass(ci1)[, 1:2], unclass(ci2)[, 1:2])
})

test_that("fixed-effect-null bootstrap contrast is identical across threads", {
  skip_on_cran()
  fit <- lmm(y ~ x + (1 | g), mm_threads_lmm_data(),
             control = mm_control(verbose = -1))
  L <- matrix(c(0, 1), nrow = 1, dimnames = list("x", names(fixef(fit))))
  c1 <- contrast(fit, L, method = "bootstrap",
                 bootstrap = bootstrap_control(nsim = 30, seed = 5))
  c2 <- contrast(fit, L, method = "bootstrap",
                 bootstrap = bootstrap_control(nsim = 30, seed = 5, threads = 2))
  expect_identical(c1$table$p_value, c2$table$p_value)
  expect_identical(c1$table$statistic, c2$table$statistic)
})

test_that("bootstrap LRT is identical across threads", {
  skip_on_cran()
  d <- mm_threads_lmm_data()
  m0 <- lmm(y ~ 1 + (1 | g), d, REML = FALSE, control = mm_control(verbose = -1))
  m1 <- lmm(y ~ x + (1 | g), d, REML = FALSE, control = mm_control(verbose = -1))
  b1 <- parametric_bootstrap(m0, m1, nsim = 25, seed = 9)
  b2 <- parametric_bootstrap(m0, m1, nsim = 25, seed = 9, threads = 2)
  expect_identical(b1$simulated, b2$simulated)
  expect_identical(b1$p_value, b2$p_value)
  expect_identical(b2$threads, 2L)
})

test_that("GLMM parametric bootstrap is identical across threads", {
  skip_on_cran()
  set.seed(7)
  d <- data.frame(g = factor(rep(1:10, each = 8)), x = rnorm(80))
  d$y <- rpois(80, exp(0.3 + 0.4 * d$x + rnorm(10, sd = 0.5)[d$g]))
  fit <- glmm(y ~ x + (1 | g), d, family = poisson(),
              method = "pirls_profiled", control = mm_control(verbose = -1))
  ci1 <- confint(fit, method = "bootstrap", nsim = 20, seed = 2)
  ci2 <- confint(fit, method = "bootstrap", nsim = 20, seed = 2, threads = 2)
  expect_identical(unclass(ci1)[, 1:2], unclass(ci2)[, 1:2])
  expect_identical(attr(ci2, "mm_bootstrap")$threads, 2L)
})
