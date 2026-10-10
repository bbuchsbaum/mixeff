# PR-4 coverage: diagnostic plots drawn to a null device (ask = TRUE paging,
# empty and SE-less random-effect panels, lattice panel functions).

cp4pl_fit <- function() {
  set.seed(491)
  d <- data.frame(x = rnorm(48), g = gl(8, 6))
  d$y <- 1 + 0.4 * d$x + rnorm(8)[d$g] + rnorm(48, sd = 0.5)
  lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
}

test_that("plot() with ask = TRUE restores the page prompt", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  fit <- cp4pl_fit()
  before <- grDevices::devAskNewPage()
  expect_identical(plot(fit, which = 1:3, ask = TRUE), fit) |>
    expect_invisible()
  expect_identical(grDevices::devAskNewPage(), before)
})

test_that("random-effect plots handle no panels and missing SEs", {
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  fit <- cp4pl_fit()
  empty <- structure(list(), class = class(ranef(fit)))
  expect_invisible(plot(empty))
  re <- ranef(fit)  # no condVar: SEs are NA
  panels <- mixeff:::mm_ranef_panels(re)
  expect_length(panels, 1L)
  expect_true(all(is.na(panels[[1]]$se)))
  expect_invisible(plot(re))
})

test_that("lattice dotplot()/qqmath() panels draw with and without SEs", {
  skip_if_not_installed("lattice")
  grDevices::pdf(NULL)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev), add = TRUE)
  fit <- cp4pl_fit()
  rc <- ranef(fit, condVar = TRUE)
  dp <- lattice::dotplot(rc)
  expect_s3_class(dp$g, "trellis")
  expect_no_error(print(dp$g))
  qq <- lattice::qqmath(rc)
  expect_s3_class(qq$g, "trellis")
  expect_no_error(print(qq$g))
  plain <- ranef(fit)
  expect_no_error(print(lattice::dotplot(plain)$g))
  expect_no_error(print(lattice::qqmath(plain)$g))
})
