# Pre-fit structural guards (PRD §8.1). These refuse degenerate inputs with a
# typed condition before crossing the Rust bridge, instead of letting the fit
# driver panic (empty data) or silently return a misleading fit
# (one-observation-per-group, n = 1).

ctl <- mm_control(verbose = -1)

test_that("empty data is refused with a typed data error (no Rust panic)", {
  s <- lme4::sleepstudy[0, ]
  expect_error(
    lmm(Reaction ~ Days + (Days | Subject), s, control = ctl),
    class = "mm_data_error"
  )
})

test_that("a single observation is refused as not identifiable", {
  s <- lme4::sleepstudy[1, , drop = FALSE]
  expect_error(
    lmm(Reaction ~ Days + (Days | Subject), s, control = ctl),
    class = "mm_not_identifiable"
  )
})

test_that("one observation per grouping level is refused (lme4's rule)", {
  # 18 rows, Subject relabelled so every row is its own level: n_levels == nobs.
  d <- lme4::sleepstudy[1:18, ]
  d$Subject <- factor(seq_len(18))
  expect_error(
    lmm(Reaction ~ Days + (1 | Subject), d, control = ctl),
    class = "mm_not_identifiable"
  )
})

test_that("empty data is refused for glmm() too", {
  expect_error(
    glmm(r2 ~ 1 + (1 | id), lme4::VerbAgg[0, ], family = binomial, control = ctl),
    class = "mm_data_error"
  )
})

test_that("an observation-level random effect is NOT refused for a GLMM", {
  # OLRE (one grouping level per observation) is a valid overdispersion device
  # for GLMMs and is accepted by glmer; the levels < observations rule is
  # LMM-only. Regression for an over-aggressive guard.
  set.seed(1L)
  d <- data.frame(
    y = rpois(40L, 3),
    x = rnorm(40L),
    g = factor(rep(seq_len(8L), each = 5L)),
    obs = factor(seq_len(40L))
  )
  expect_s3_class(
    glmm(y ~ x + (1 | g) + (1 | obs), d, family = poisson, control = ctl),
    "mm_glmm"
  )
})

test_that("valid models across topologies still fit (no false positives)", {
  expect_s3_class(
    lmm(Reaction ~ Days + (Days | Subject), lme4::sleepstudy, control = ctl),
    "mm_lmm"
  )
  expect_s3_class(
    lmm(angle ~ recipe * temperature + (1 | recipe:replicate), lme4::cake,
        control = ctl),
    "mm_lmm"
  )
  expect_s3_class(
    lmm(strength ~ 1 + (1 | batch/cask), lme4::Pastes, control = ctl),
    "mm_lmm"
  )
  expect_s3_class(
    lmm(diameter ~ 1 + (1 | plate) + (1 | sample), lme4::Penicillin,
        control = ctl),
    "mm_lmm"
  )
  # A singular-but-identifiable fit (n_levels < nobs) must NOT be refused.
  expect_s3_class(
    lmm(Yield ~ 1 + (1 | Batch), lme4::Dyestuff2, control = ctl),
    "mm_lmm"
  )
})
