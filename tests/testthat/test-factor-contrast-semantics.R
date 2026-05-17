factor_contrast_design <- function(seed = 8L) {
  set.seed(seed)
  group_count <- 18L
  reps <- 4L
  dat <- expand.grid(
    g = factor(seq_len(group_count)),
    f = factor(c("a", "b")),
    rep = seq_len(reps)
  )
  x <- ifelse(dat$f == "b", 0.5, -0.5)
  b_a <- rnorm(group_count, sd = 0.35)
  b_b <- rnorm(group_count, sd = 0.30)
  dat$y <- 1.8 + 0.7 * x +
    ifelse(dat$f == "a", b_a[as.integer(dat$g)], b_b[as.integer(dat$g)]) +
    rnorm(nrow(dat), sd = 0.2)
  contrasts(dat$f) <- matrix(
    c(-0.5, 0.5),
    ncol = 1L,
    dimnames = list(c("a", "b"), "half")
  )
  dat
}

test_that("no-intercept factor random slopes use cell-means coding", {
  dat <- factor_contrast_design()
  fit <- lmm(y ~ f + (0 + f | g), dat, control = mm_control(verbose = -1))

  vc_names <- VarCorr(fit)$table$name
  re_names <- names(ranef(fit)$g)
  expect_true(all(c("f: a", "f: b") %in% vc_names))
  expect_true(all(c("f: a", "f: b") %in% re_names))
  expect_false(any(grepl("half", vc_names, fixed = TRUE)))
  expect_false(any(grepl("half", re_names, fixed = TRUE)))

  if (requireNamespace("lme4", quietly = TRUE)) {
    ref <- suppressMessages(suppressWarnings(
      lme4::lmer(y ~ f + (0 + f | g), dat)
    ))
    ref_vc <- as.data.frame(lme4::VarCorr(ref))
    ref_names <- ref_vc$var1[ref_vc$grp == "g" & is.na(ref_vc$var2)]
    expect_setequal(ref_names, c("fa", "fb"))
  }
})
