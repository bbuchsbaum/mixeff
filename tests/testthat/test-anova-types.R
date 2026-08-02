# WI-9.1 (Audit1): numerical references for anova(type = "I"/"II"/"III").
#
# The pre-existing typed-anova test (test-lmerTest-parity.R, "anova(type=...)
# routes to Rust-owned typed term hypotheses") only checks that the requested
# type is echoed in the row notes. This file pins the actual statistics on an
# unbalanced two-factor design where the three types must differ:
#
#   A (2 levels) x B (3 levels), deliberately unequal cell sizes
#   (8, 14, 11, 5, 13, 9), one continuous covariate x, (1 | g) random
#   intercept, REML = FALSE.
#
# Contrast coding assumed throughout: the package defaults enforced by
# mm_translate_data() -- contr.treatment for unordered factors, contr.poly
# for ordered factors (see test-data-translate.R). A and B are unordered, so
# every factor column below is treatment-coded.
#
# Reference computations, all in the F family that anova(method =
# "satterthwaite") produces (single-df t rows are folded to F = t^2 by
# anova.mm_lmm, so every row here is an F row):
#   * Type I  -- own construction: sequential Doolittle contrast basis of
#     X'X (the construction named in the upstream contract,
#     mixeff-rs/docs/fixed_effect_p_values_plan.md), turned into Wald F via
#     the fit's own fixef()/vcov(). Cross-checked externally against
#     lmerTest::anova(type = 1, ddf = "Satterthwaite").
#   * Type II -- own construction: for each term, its model-matrix columns
#     residualized (OLS) against the columns of every term that does not
#     contain it (marginality: main effects adjusted for other main effects
#     and the covariate, never for the interaction), then Wald F on the
#     resulting hypothesis row space. Cross-checked against lmerTest type 2.
#   * Type III -- lmerTest::anova(type = 3, ddf = "Satterthwaite"), all
#     rows (marginal Type III since engine 1f3f689), plus a
#     reference-level invariance test and coverage of the honestly named
#     `type = "block"` (the pre-1f3f689 coefficient-block hypothesis).
#
# Tolerances: the own constructions are the same estimator on the same fit,
# so they must agree to numerical noise (measured ~1e-15 relative; asserted
# at 1e-8). lmerTest refits with its own optimizer, so those comparisons use
# the suite's parity tolerances (1e-3 on F, den_df within 1e-2 relative).

mm_anova_types_data <- function() {
  set.seed(42)
  cells <- expand.grid(
    A = c("a1", "a2"),
    B = c("b1", "b2", "b3"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  # Unequal cell counts: this is what makes Types I/II/III disagree.
  cells$n <- c(8L, 14L, 11L, 5L, 13L, 9L)
  df <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    data.frame(
      A = cells$A[i], B = cells$B[i],
      stringsAsFactors = FALSE
    )[rep(1, cells$n[i]), , drop = FALSE]
  }))
  df$A <- factor(df$A)
  df$B <- factor(df$B)
  n <- nrow(df)
  df$g <- factor(rep_len(sprintf("g%02d", 1:10), n))
  df$x <- round(rnorm(n), 3)
  u <- rnorm(10, sd = 0.8)
  df$y <- round(
    2 + 0.9 * (df$A == "a2") + 0.6 * (df$B == "b2") - 0.4 * (df$B == "b3") +
      1.2 * (df$A == "a2" & df$B == "b3") + 0.5 * df$x +
      u[as.integer(df$g)] + rnorm(n, sd = 0.7),
    4
  )
  rownames(df) <- NULL
  df
}

# Wald F for the hypothesis L beta = 0 computed from the fit's own
# fixef()/vcov(): F = (L b)' (L C L')^{-1} (L b) / nrow(L). The F statistic
# does not depend on the (Satterthwaite) denominator df, so this pins the
# hypothesis row space and numerator exactly.
mm_own_wald_f <- function(L, b, C) {
  L <- matrix(L, ncol = length(b))
  est <- L %*% b
  as.numeric(t(est) %*% solve(L %*% C %*% t(L), est)) / nrow(L)
}

# Unit-lower-triangular Doolittle factor of a symmetric positive-definite
# matrix (M = L U with diag(L) = 1). Rows of t(L) are the sequential
# contrast basis: row j tests column j adjusted for columns 1..j-1 only.
mm_doolittle_l <- function(M) {
  p <- ncol(M)
  L <- diag(p)
  U <- M
  for (i in seq_len(p - 1L)) {
    for (j in (i + 1L):p) {
      L[j, i] <- U[j, i] / U[i, i]
      U[j, ] <- U[j, ] - L[j, i] * U[i, ]
    }
  }
  L
}

# mixeff anova table with rows keyed by term, plus basic shape checks.
mm_anova_rows <- function(fit, type) {
  tab <- stats::anova(fit, type = type, method = "satterthwaite")$table
  expect_true(all(tab$status == "available"))
  expect_true(all(tab$statistic_name == "f"))
  rownames(tab) <- tab$term
  tab
}

mm_lmerTest_type_table <- function(ref, type) {
  out <- as.data.frame(stats::anova(ref, type = type, ddf = "Satterthwaite"))
  out$term <- rownames(out)
  rownames(out) <- out$term
  out
}

test_that("Type I matches its own sequential Doolittle reference and lmerTest type 1", {
  mm_skip_if_no_lmerTest()
  df <- mm_anova_types_data()
  # Covariate first: the literal formula-expansion order (x, A, B, A:B) then
  # coincides with R's terms() order, so the sequential decomposition is
  # unambiguous and lmerTest's Type I is a valid external reference.
  fit <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  ref <- mm_lme4_reference_or_skip(
    lmerTest::lmer(y ~ x + A * B + (1 | g), df, REML = FALSE),
    "anova_types_unbalanced", "lmerTest::lmer()"
  )
  # Both engines must sit on the same optimum before statistics are compared.
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-6)

  obs <- mm_anova_rows(fit, "I")
  terms <- c("x", "A", "B", "A:B")
  expect_setequal(obs$term, terms)

  # Own construction: sequential (type I) contrasts are the rows of
  # t(doolittle(X'X)), partitioned by term; identical estimator, so 1e-8.
  X <- stats::model.matrix(~ x + A * B, df)
  asgn <- attr(X, "assign")
  Lseq <- t(mm_doolittle_l(crossprod(X)))
  b <- fixef(fit)
  expect_identical(names(b), colnames(X)) # coef-map contract guards L layout
  C <- as.matrix(vcov(fit))
  own_f <- vapply(seq_along(terms), function(j) {
    mm_own_wald_f(Lseq[asgn == j, , drop = FALSE], b, C)
  }, numeric(1))
  expect_equal(obs[terms, "statistic"], own_f, tolerance = 1e-8)

  # External reference: lmerTest Type I with Satterthwaite ddf.
  exp_tab <- mm_lmerTest_type_table(ref, 1)
  expect_equal(obs[terms, "statistic"], exp_tab[terms, "F value"],
               tolerance = 1e-3)
  expect_equal(obs[terms, "num_df"], exp_tab[terms, "NumDF"],
               tolerance = 1e-8)
  expect_equal(obs[terms, "den_df"], exp_tab[terms, "DenDF"],
               tolerance = 1e-2)
  expect_equal(obs[terms, "p_value"], exp_tab[terms, "Pr(>F)"],
               tolerance = 1e-3)
})

test_that("Type I sequences terms in terms() order (interaction degree)", {
  # Engine 1f3f689 (upstream bd-01KZ15535G) fixed Type I sequencing: terms
  # are ordered by interaction degree, stable within a degree, matching R's
  # terms() convention -- so `y ~ A * B + x` sequences A, B, x, A:B, and a
  # covariate written after an interaction is no longer adjusted for it.
  # Before the fix the engine used literal expansion order (A, B, A:B, x),
  # giving x an F of 56.0998 instead of 63.3110 on this design.
  df <- mm_anova_types_data()
  fit <- lmm(y ~ A * B + x + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  obs <- mm_anova_rows(fit, "I")
  expect_identical(obs$term, c("A", "B", "x", "A:B"))

  # Own construction in that order: permute the R model-matrix columns into
  # terms() order, take the sequential Doolittle basis there, and map each
  # contrast row back to the beta layout used by fixef()/vcov(). R's assign
  # ids for ~ A * B + x are A = 1, B = 2, x = 3, A:B = 4, which IS
  # terms() order, so no permutation of the assign sequence is needed.
  X <- stats::model.matrix(~ A * B + x, df)
  asgn <- attr(X, "assign")
  seq_terms <- c(1L, 2L, 3L, 4L)
  ord <- order(match(asgn, c(0L, seq_terms)))
  Lseq <- t(mm_doolittle_l(crossprod(X[, ord])))
  b <- fixef(fit)
  expect_identical(names(b), colnames(X))
  C <- as.matrix(vcov(fit))
  own_f <- vapply(seq_terms, function(j) {
    rows <- Lseq[asgn[ord] == j, , drop = FALSE]
    L <- matrix(0, nrow(rows), ncol(X))
    L[, ord] <- rows
    mm_own_wald_f(L, b, C)
  }, numeric(1))
  expect_equal(obs[c("A", "B", "x", "A:B"), "statistic"], own_f,
               tolerance = 1e-8)

  # Writing the covariate before the interaction must now give the same
  # sequence (both spellings canonicalize to terms() order).
  fit2 <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
              control = mm_control(verbose = -1))
  obs2 <- mm_anova_rows(fit2, "I")
  expect_identical(obs2$term, c("x", "A", "B", "A:B"))
})

test_that("Type II matches the marginality-respecting own construction and lmerTest type 2", {
  mm_skip_if_no_lmerTest()
  df <- mm_anova_types_data()
  fit <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  ref <- mm_lme4_reference_or_skip(
    lmerTest::lmer(y ~ x + A * B + (1 | g), df, REML = FALSE),
    "anova_types_unbalanced", "lmerTest::lmer()"
  )
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-6)

  obs <- mm_anova_rows(fit, "II")
  terms <- c("x", "A", "B", "A:B")

  # Own construction of the Type II hypothesis for each term: residualize
  # the term's model-matrix columns against the columns of every term that
  # does NOT contain it (marginality: A is adjusted for B and x but not for
  # A:B; the interaction, contained in no other term, is adjusted for
  # everything). The hypothesis row space in beta coordinates is then
  # rowspace(Xj_perp' X), and the Wald F is invariant to the row basis.
  X <- stats::model.matrix(~ x + A * B, df)
  asgn <- attr(X, "assign")
  trm <- stats::terms(y ~ x + A * B)
  fac <- attr(trm, "factors")[-1L, , drop = FALSE] # drop response row
  labs <- attr(trm, "term.labels")
  b <- fixef(fit)
  expect_identical(names(b), colnames(X))
  C <- as.matrix(vcov(fit))
  own_f <- vapply(seq_along(labs), function(j) {
    contains_j <- vapply(seq_along(labs), function(k) {
      k != j && all(fac[, j] <= fac[, k])
    }, logical(1))
    keep <- asgn %in% c(0L, which(!contains_j & seq_along(labs) != j))
    Xn <- X[, keep, drop = FALSE]
    Xj <- X[, asgn == j, drop = FALSE]
    Xj_perp <- Xj - Xn %*% qr.coef(qr(Xn), Xj)
    mm_own_wald_f(t(Xj_perp) %*% X, b, C)
  }, numeric(1))
  names(own_f) <- labs
  expect_equal(obs[labs, "statistic"], unname(own_f), tolerance = 1e-8)

  # External reference: lmerTest Type II with Satterthwaite ddf.
  exp_tab <- mm_lmerTest_type_table(ref, 2)
  expect_equal(obs[terms, "statistic"], exp_tab[terms, "F value"],
               tolerance = 1e-3)
  expect_equal(obs[terms, "num_df"], exp_tab[terms, "NumDF"],
               tolerance = 1e-8)
  expect_equal(obs[terms, "den_df"], exp_tab[terms, "DenDF"],
               tolerance = 1e-2)
  expect_equal(obs[terms, "p_value"], exp_tab[terms, "Pr(>F)"],
               tolerance = 1e-3)
})

test_that("Type III matches lmerTest type 3 where the hypotheses coincide", {
  mm_skip_if_no_lmerTest()
  df <- mm_anova_types_data()
  fit <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  ref <- mm_lme4_reference_or_skip(
    lmerTest::lmer(y ~ x + A * B + (1 | g), df, REML = FALSE),
    "anova_types_unbalanced", "lmerTest::lmer()"
  )
  expect_equal(as.numeric(logLik(fit)), as.numeric(stats::logLik(ref)),
               tolerance = 1e-6)

  obs <- mm_anova_rows(fit, "III")
  exp_tab <- mm_lmerTest_type_table(ref, 3)

  # The covariate and the highest-order interaction have the same Type III
  # hypothesis under every construction; these rows must agree.
  for (term in c("x", "A:B")) {
    expect_equal(obs[term, "statistic"], exp_tab[term, "F value"],
                 tolerance = 1e-3)
    expect_equal(obs[term, "num_df"], exp_tab[term, "NumDF"],
                 tolerance = 1e-8)
    expect_equal(obs[term, "den_df"], exp_tab[term, "DenDF"],
                 tolerance = 1e-2)
    expect_equal(obs[term, "p_value"], exp_tab[term, "Pr(>F)"],
                 tolerance = 1e-3)
  }

  # Engine 1f3f689 (upstream bd-01KZ1551AE) made Type III MARGINAL: a term's
  # own contrast columns plus the equally-weighted average over the levels
  # of every containing term. Main-effect rows now agree with lmerTest on
  # this unbalanced treatment-coded design, where the previous
  # coefficient-block hypothesis did not (A: 0.2914 block vs 49.9498
  # marginal -- opposite conclusions from the same fit).
  expect_equal(obs["A", "statistic"], exp_tab["A", "F value"],
               tolerance = 1e-3)
  expect_equal(obs["B", "statistic"], exp_tab["B", "F value"],
               tolerance = 1e-3)
})

test_that("Type III is invariant to the reference level; block is not", {
  # The defining property of a marginal Type III hypothesis: relabeling
  # which level of a factor is the reference must not move the test. The
  # coefficient-block hypothesis is reference-dependent by construction
  # (it IS the simple effect at the reference), so it must move.
  df <- mm_anova_types_data()
  df_relevel <- df
  df_relevel$B <- stats::relevel(df_relevel$B, ref = "b3")

  f1 <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
            control = mm_control(verbose = -1))
  f2 <- lmm(y ~ x + A * B + (1 | g), df_relevel, REML = FALSE,
            control = mm_control(verbose = -1))

  a1 <- mm_anova_rows(f1, "III")
  a2 <- mm_anova_rows(f2, "III")
  expect_equal(a1["A", "statistic"], a2["A", "statistic"], tolerance = 1e-6)
  expect_equal(a1["A:B", "statistic"], a2["A:B", "statistic"],
               tolerance = 1e-6)

  b1 <- mm_anova_rows(f1, "block")
  b2 <- mm_anova_rows(f2, "block")
  expect_gt(abs(b1["A", "statistic"] - b2["A", "statistic"]) /
              max(abs(b1["A", "statistic"]), abs(b2["A", "statistic"])), 0.05)
})

test_that("type = 'block' reproduces the coefficient-block hypothesis", {
  # The block rows must equal a Wald F on the identity block over each
  # term's own columns, computed from the fit's own fixef()/vcov().
  df <- mm_anova_types_data()
  fit <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  obs <- mm_anova_rows(fit, "block")
  X <- stats::model.matrix(~ x + A * B, df)
  asgn <- attr(X, "assign")
  b <- fixef(fit)
  C <- as.matrix(vcov(fit))
  labels <- attr(stats::terms(y ~ x + A * B), "term.labels")
  for (j in seq_along(labels)) {
    cols <- which(asgn == j)
    L <- matrix(0, length(cols), length(b))
    for (i in seq_along(cols)) L[i, cols[i]] <- 1
    expect_equal(obs[labels[[j]], "statistic"], mm_own_wald_f(L, b, C),
                 tolerance = 1e-8,
                 info = sprintf("block row %s", labels[[j]]))
  }
})
test_that("Types I, II, and III actually differ on the unbalanced design", {
  # Guards the unbalanced-ness of the fixture: if the cell counts were ever
  # rebalanced, the three types would collapse and the tests above would
  # stop distinguishing them. Rows are chosen so the assertions hold under
  # both the current block-form Type III and a marginal Type III:
  #   * x row separates I from II and from III (sequential-at-position vs
  #     fully adjusted; III's x row is fully adjusted in either semantics).
  #   * A row separates II from III (62.50 vs 11.94 today; vs 56.51 if the
  #     FINDING above is fixed -- still > 5% apart).
  df <- mm_anova_types_data()
  fit <- lmm(y ~ x + A * B + (1 | g), df, REML = FALSE,
             control = mm_control(verbose = -1))
  f1 <- mm_anova_rows(fit, "I")
  f2 <- mm_anova_rows(fit, "II")
  f3 <- mm_anova_rows(fit, "III")

  rel_gap <- function(a, b) abs(a - b) / pmax(abs(a), abs(b))
  expect_gt(rel_gap(f1["x", "statistic"], f2["x", "statistic"]), 0.05)
  expect_gt(rel_gap(f1["x", "statistic"], f3["x", "statistic"]), 0.05)
  expect_gt(rel_gap(f2["A", "statistic"], f3["A", "statistic"]), 0.05)
  expect_gt(rel_gap(f1["A", "statistic"], f2["A", "statistic"]), 0.05)
})
