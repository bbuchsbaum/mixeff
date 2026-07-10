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
  # VarCorr display keeps the engine's pretty labels; ranef() columns carry
  # lme4-concatenated names (mm_apply_lme4_coef_naming strips the separator).
  expect_true(all(c("f: a", "f: b") %in% vc_names))
  expect_true(all(c("fa", "fb") %in% re_names))
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

test_that("as.data.frame.mm_varcorr var1/var2 use lme4-compatible concatenated names", {
  dat <- factor_contrast_design()
  fit <- lmm(y ~ f + (0 + f | g), dat, control = mm_control(verbose = -1))

  # Display path keeps human-readable "f: a" format
  expect_true(all(c("f: a", "f: b") %in% VarCorr(fit)$table$name))

  # Serialisation path (as.data.frame) must match lme4's "fa"/"fb" convention
  df <- as.data.frame(VarCorr(fit))
  diag_names <- df$var1[df$grp == "g" & is.na(df$var2)]
  expect_false(any(grepl(": ", diag_names, fixed = TRUE)),
               info = "var1 must not contain ': ' — use lme4 concatenated form")
  expect_setequal(diag_names, c("fa", "fb"))
})

# Ordered-factor (contr.poly) semantics. `o` is a 3-level ordered factor used
# both as a fixed effect and inside a random-slope term.
ordered_factor_design <- function(seed = 11L) {
  set.seed(seed)
  group_count <- 20L
  reps <- 6L
  dat <- expand.grid(
    g = factor(seq_len(group_count)),
    o = factor(c("lo", "mid", "hi"), levels = c("lo", "mid", "hi"),
               ordered = TRUE),
    rep = seq_len(reps)
  )
  bg <- rnorm(group_count, sd = 0.4)[as.integer(dat$g)]
  dat$y <- 1 + 0.5 * as.integer(dat$o) + bg + rnorm(nrow(dat), sd = 0.5)
  dat
}

test_that("ordered fixed factor is coded with contr.poly, matching lme4", {
  mm_skip_if_no_lme4()
  dat <- ordered_factor_design()
  fit <- lmm(y ~ o + (1 | g), dat, REML = TRUE, control = mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(lme4::lmer(y ~ o + (1 | g), dat, REML = TRUE)))

  # Polynomial trend labels in lme4 form, not treatment level labels.
  expect_true(any(grepl("o.L", names(fixef(fit)), fixed = TRUE)))
  expect_false(any(grepl("omid", names(fixef(fit)), fixed = TRUE)))

  ref_beta <- lme4::fixef(ref)
  expect_identical(names(fixef(fit)), names(ref_beta))
  expect_equal(fixef(fit), ref_beta, tolerance = 1e-5)
})

test_that("ordered factor in a random slope uses contr.poly Z coding, matching lme4", {
  mm_skip_if_no_lme4()
  dat <- ordered_factor_design()
  fit <- lmm(y ~ o + (1 + o | g), dat, REML = TRUE, control = mm_control(verbose = -1))
  ref <- suppressMessages(suppressWarnings(lme4::lmer(y ~ o + (1 + o | g), dat, REML = TRUE)))

  # Random-effect design uses polynomial trend columns, not cell-means levels.
  vc <- VarCorr(fit)$table
  expect_true(all(c("(Intercept)", "o: .L", "o: .Q") %in% vc$name))
  expect_false(any(grepl("o: lo", vc$name, fixed = TRUE)))

  # Variance components are basis-dependent: matching lme4's poly-coded values
  # (not merely the likelihood) confirms the Z columns share lme4's contr.poly
  # basis. Treatment/cell-means Z would give wildly different components at the
  # same logLik.
  mm_var <- stats::setNames(vc$variance, gsub(": ", "", vc$name, fixed = TRUE))
  ref_vc <- as.data.frame(lme4::VarCorr(ref))
  ref_diag <- ref_vc[ref_vc$grp == "g" & is.na(ref_vc$var2), ]
  ref_var <- stats::setNames(ref_diag$vcov, ref_diag$var1)
  expect_setequal(names(mm_var), names(ref_var))
  expect_equal(mm_var[names(ref_var)], ref_var, tolerance = 2e-2)

  # Likelihood parity too (necessary but not sufficient on its own).
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-4)
})

test_that("ordered factor fits (poly-coded) under the standard UNNAMED contrasts option", {
  mm_skip_if_no_lme4()
  dat <- ordered_factor_design()
  # The default global option form is unnamed; the ordered-contrast guard must
  # resolve it positionally rather than throwing "subscript out of bounds".
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  fit <- tryCatch(
    lmm(y ~ o + (1 | g), dat, REML = TRUE, control = mm_control(verbose = -1)),
    error = function(e) e
  )
  options(old)
  expect_false(inherits(fit, "error"))
  expect_true(any(grepl("o.L", names(fixef(fit)), fixed = TRUE)))
  ref <- suppressMessages(suppressWarnings(lme4::lmer(y ~ o + (1 | g), dat, REML = TRUE)))
  expect_equal(as.numeric(stats::logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-4)
})
