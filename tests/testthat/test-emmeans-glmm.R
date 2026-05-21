## emmeans support for mm_glmm: recover_data.mm_glmm + emm_basis.mm_glmm,
## with linkinv applied automatically on type = "response".

skip_if_not_installed("emmeans")
skip_if_not_installed("estimability")

mk_emm_glmm <- function(seed = 71L) {
  set.seed(seed)
  n_subj <- 16L
  n_per <- 12L
  subject <- factor(rep(seq_len(n_subj), each = n_per))
  x <- rnorm(n_subj * n_per)
  g <- factor(rep(c("ctrl", "treat"), length.out = n_subj * n_per))
  b0 <- rnorm(n_subj, sd = 0.4)
  eta <- -0.1 + 0.5 * x + 0.3 * (g == "treat") - 0.4 * x * (g == "treat") +
    b0[as.integer(subject)]
  y <- rbinom(length(eta), 1L, plogis(eta))
  glmm(
    y ~ x * g + (1 | subject),
    data.frame(y = y, x = x, g = g, subject = subject),
    family = binomial(),
    method = "pirls_profiled",
    nAGQ = 1L,
    control = mm_control(verbose = -1)
  )
}

test_that("emmeans(mm_glmm) returns a reference grid for a categorical factor", {
  fit <- mk_emm_glmm()
  em <- emmeans::emmeans(fit, ~ g)
  expect_s4_class(em, "emmGrid")
  s <- as.data.frame(summary(em))
  expect_true(all(c("g", "emmean", "SE") %in% colnames(s)))
  expect_equal(nrow(s), 2L)
  ## Link scale by default
  expect_true(all(is.finite(s$emmean)))
})

test_that("emmeans(mm_glmm, type='response') applies linkinv automatically", {
  fit <- mk_emm_glmm()
  em <- emmeans::emmeans(fit, ~ g, type = "response")
  s <- as.data.frame(summary(em))
  ## Response-scale column name from emmeans for binomial-logit is "prob"
  prob_col <- intersect(c("prob", "response", "rate"), colnames(s))
  expect_true(length(prob_col) >= 1L,
              info = "response-scale column not found in emmeans output")
  prob <- s[[prob_col[1L]]]
  expect_true(all(prob > 0 & prob < 1),
              info = "response-scale predictions should be probabilities in (0,1)")
})

test_that("emmeans response-scale matches manual linkinv on link-scale estimate", {
  fit <- mk_emm_glmm()
  em_link <- emmeans::emmeans(fit, ~ g, type = "link")
  em_resp <- emmeans::emmeans(fit, ~ g, type = "response")
  s_link <- as.data.frame(summary(em_link))
  s_resp <- as.data.frame(summary(em_resp))

  expected <- plogis(s_link$emmean)
  prob_col <- intersect(c("prob", "response", "rate"), colnames(s_resp))[1L]
  expect_equal(s_resp[[prob_col]], expected, tolerance = 1e-10)
})

test_that("emmeans pairwise contrast on mm_glmm uses asymptotic z (df = Inf)", {
  fit <- mk_emm_glmm()
  em <- emmeans::emmeans(fit, ~ g)
  pw <- as.data.frame(emmeans::contrast(em, method = "pairwise"))
  expect_true("z.ratio" %in% colnames(pw))
  expect_true(any(c("p.value") %in% colnames(pw)))
  expect_true(all(is.finite(pw$z.ratio)))
})

test_that("emm_basis.mm_glmm carries fixed-effect bhat and full V from mm_glmm", {
  fit <- mk_emm_glmm()
  rd <- recover_data.mm_glmm(fit)
  trms <- attr(rd, "terms")
  xlev <- attr(rd, "xlev") %||% list()
  grid <- expand.grid(g = factor(c("ctrl", "treat"), levels = c("ctrl", "treat")),
                       x = 0)
  b <- emm_basis.mm_glmm(fit, trms = trms, xlev = xlev, grid = grid)
  expect_named(b, c("X", "bhat", "nbasis", "V", "dffun", "dfargs", "misc"),
               ignore.order = TRUE)
  expect_equal(length(b$bhat), length(fixef(fit)))
  expect_equal(nrow(b$V),  length(fixef(fit)))
  expect_equal(ncol(b$V),  length(fixef(fit)))
  ## dffun returns Inf (asymptotic)
  expect_identical(b$dffun(rep(0, length(b$bhat)), b$dfargs), Inf)
})

test_that("emm_basis.mm_glmm misc carries family link info for type='response'", {
  fit <- mk_emm_glmm()
  rd <- recover_data.mm_glmm(fit)
  trms <- attr(rd, "terms")
  xlev <- attr(rd, "xlev") %||% list()
  grid <- expand.grid(g = factor(c("ctrl", "treat"), levels = c("ctrl", "treat")),
                       x = 0)
  b <- emm_basis.mm_glmm(fit, trms = trms, xlev = xlev, grid = grid)
  ## emmeans::.std.link.labels populates either misc$tran or misc$inv.lbl
  expect_true(!is.null(b$misc$tran) || !is.null(b$misc$inv.lbl),
              info = "emm_basis.mm_glmm should annotate misc with link info")
})
