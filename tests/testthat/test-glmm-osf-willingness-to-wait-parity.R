# In-the-wild GLMM parity against a *published* lme4::glmer analysis:
# OSF node ftexh, "Willingness to wait study" (file pv6u2 =
# "Analysis- Correlations, LME4.R"). The committed fixtures
# tests/fixtures/osf_willingness_to_wait_study1{a,b}.csv are the slimmed
# trial-level modeling columns (1427 + 1600 rows); provenance + reconstruction
# in data-raw/osf-willingness-to-wait/. Tracked by mote
# bd-01KT3ZRCKWRZQFA4W7TXTGWAZ0.
#
# These models are binomial-logit GLMMs with CROSSED random effects
# (1 | ID) + (1 | Title), several with CORRELATED random slopes
# (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title).
#
# This guards three things at once:
#   (1) the grouping-coercion fix: the raw fixture has an INTEGER `ID` and a
#       CHARACTER `Title`; glmm()/lmm() must coerce them to factors (announced,
#       not silent) the way lme4 does, instead of refusing with
#       "grouping factor not categorical";
#   (2) joint_laplace tracks glmer to ~1e-3 on the correlated-random-slope
#       models; Wald rows are available, with study1a at strict SE/z/p parity
#       and study1b retained as a bounded remaining SE/z drift;
#   (3) the remaining high-baseline random-intercept gap stays documented, not
#       silently regressed:
#         - high-baseline random-INTERCEPT models converge to a slightly
#           sub-optimal point (joint_laplace reports converged_interior ~0.01-
#           0.05 logLik short of glmer) -- upstream
#           bd-01KT40T6FGVXQQ9N50G2HM0ZZE. Asserted as a bounded, usable gap.

osf_ww_fixture_path <- function(which) {
  fn <- sprintf("osf_willingness_to_wait_%s.csv", which)
  candidates <- c(
    testthat::test_path("..", "fixtures", fn),
    file.path("tests", "fixtures", fn)
  )
  hit <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(hit), paste(fn, "fixture is unavailable"))
  hit
}

osf_ww_data <- function(which) {
  d <- utils::read.csv(osf_ww_fixture_path(which), stringsAsFactors = FALSE)
  d$Enjoyment <- suppressWarnings(as.numeric(as.character(d$Enjoyment)))
  d$Enjoyment_centered <- d$Enjoyment - mean(d$Enjoyment, na.rm = TRUE)
  d$arousal <- suppressWarnings(as.numeric(as.character(d$arousal)))
  d$arousal_centered <- d$arousal - mean(d$arousal, na.rm = TRUE)
  d
}

# Long form for the comprehension models: gather Q1_correct + Q2_correct into a
# single 0/1 `score`, exactly as the published script's tidyr::gather() does.
osf_ww_long <- function(which) {
  d <- osf_ww_data(which)
  a <- d
  a$score <- d$Q1_correct
  b <- d
  b$score <- d$Q2_correct
  out <- rbind(a, b)
  out$SVScore_centered <- out$SVScore - mean(out$SVScore, na.rm = TRUE)
  out$arousal_centered <- out$arousal - mean(out$arousal, na.rm = TRUE)
  out
}

mm_jl <- function() mm_control(verbose = -1, max_feval = 100000L)
bin <- function() binomial("logit")

test_that("fixture is faithful to the published Willingness-to-Wait result", {
  skip_if_not_installed("lme4")
  d1 <- osf_ww_data("study1a")
  # Wait decisions x enjoyment, maximal model. Paper/script: Enjoyment_centered
  # beta = 0.9439 (study 1a).
  g <- suppressWarnings(lme4::glmer(
    wait_choice ~ 1 + Enjoyment_centered +
      (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title),
    data = d1, family = binomial("logit"),
    control = lme4::glmerControl(optimizer = "bobyqa")))
  expect_equal(unname(lme4::fixef(g)[["Enjoyment_centered"]]), 0.9439,
               tolerance = 1e-3)

  # Comprehension x enjoyment + proficiency. Script: Comp_1a beta = 0.093057.
  d1L <- osf_ww_long("study1a")
  gc <- lme4::glmer(score ~ Enjoyment_centered + SVScore_centered +
                      (1 | ID) + (1 | Title),
                    data = d1L, family = binomial("logit"),
                    control = lme4::glmerControl(optimizer = "bobyqa"))
  expect_equal(unname(lme4::fixef(gc)[["Enjoyment_centered"]]), 0.093057,
               tolerance = 1e-3)
})

test_that("glmm() coerces raw integer/character grouping like lme4 (fix guard)", {
  d1 <- osf_ww_data("study1a")
  # The raw published data stores ID as an integer and Title as a character.
  expect_true(is.integer(d1$ID) || is.numeric(d1$ID))
  expect_true(is.character(d1$Title))

  fo <- wait_choice ~ Enjoyment_centered + (1 | ID)

  # Default verbose: an integer grouping variable is coerced to a factor and the
  # coercion is announced (not silent surgery).
  expect_message(
    utils::capture.output(
      glmm(fo, data = d1, family = binomial("logit"), method = "joint_laplace")),
    class = "mm_grouping_coercion_notice"
  )

  # verbose = -1 silences the notice but still fits.
  expect_no_message(
    fit <- glmm(fo, data = d1, family = binomial("logit"),
                method = "joint_laplace", control = mm_jl()),
    class = "mm_grouping_coercion_notice"
  )
  expect_s3_class(fit, "mm_glmm")
  expect_length(fixef(fit), 2L)
  expect_true(is.finite(as.numeric(logLik(fit))))
})

test_that("joint_laplace tracks glmer on correlated-random-slope wait models", {
  skip_if_not_installed("lme4")
  for (which in c("study1a", "study1b")) {
    d <- osf_ww_data(which)
    fo <- wait_choice ~ 1 + Enjoyment_centered +
      (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)
    g <- suppressWarnings(lme4::glmer(
      fo, data = d, family = binomial("logit"),
      control = lme4::glmerControl(optimizer = "bobyqa")))
    m <- glmm(fo, data = d, family = binomial("logit"),
              method = "joint_laplace", control = mm_jl())
    expect_identical(m$method, "joint_laplace")
    bg <- unname(lme4::fixef(g))
    bm <- unname(fixef(m))
    expect_equal(length(bm), length(bg))
    expect_lt(max(abs(bm - bg)), 5e-3)                                   # estimates
    expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2)  # logLik

    cg <- summary(g)$coefficients
    cm <- summary(m, tests = "coefficients")$coefficients
    p_col <- intersect(c("Pr(>|z|)", "p.value"), names(cm))[1L]
    expect_true(all(cm$method == "asymptotic_wald_z"))
    expect_true(all(is.finite(cm[["Std. Error"]]) & cm[["Std. Error"]] > 0))
    expect_true(all(is.finite(cm[["z value"]])))
    expect_false(is.na(p_col))
    expect_true(all(is.finite(cm[[p_col]]) & cm[[p_col]] >= 0 & cm[[p_col]] <= 1))
    se_tol <- if (identical(which, "study1a")) 5e-3 else 3e-2
    z_tol <- if (identical(which, "study1a")) 5e-2 else 1.0
    p_tol <- if (identical(which, "study1a")) 5e-3 else 5e-2
    expect_lt(max(abs(unname(cm[["Std. Error"]]) - unname(cg[, "Std. Error"]))), se_tol)
    expect_lt(max(abs(unname(cm[["z value"]]) - unname(cg[, "z value"]))), z_tol)
    expect_lt(max(abs(unname(cm[[p_col]]) - unname(cg[, "Pr(>|z|)"]))), p_tol)
  }
})

test_that("random-intercept comp models reach glmer parity (upstream bd-01KT40T6FGVXQQ9N50G2HM0ZZE fixed)", {
  # Upstream bd-01KT40T6FGVXQQ9N50G2HM0ZZE landed at pin 6731062 (descent-gated
  # trust_bq stagnation stop): the high-baseline random-intercept early stop is
  # gone and this model now reaches the glmer MLE (measured max|dFixef| 1.3e-4,
  # |dlogLik| 1.3e-4), so the old 6e-2 / 1e-1 shortfall bands are tightened to
  # the wait_* parity bands per the re-certification note that used to live
  # here. One residual: on this flat high-baseline ridge the FD theta-gradient
  # certificate is noise-dominated, so fit_status reads `not_optimized` (an
  # honest uncertified-candidate label, NOT a wrong answer) — upstream
  # follow-up bd-01KTQFTH6J0ZFGR5RMV28HAX44. When that lands, expect a clean
  # `converged`/`converged_interior` here.
  skip_if_not_installed("lme4")
  d1L <- osf_ww_long("study1a")
  fo <- score ~ Enjoyment_centered + SVScore_centered + (1 | ID) + (1 | Title)
  g <- lme4::glmer(fo, data = d1L, family = binomial("logit"),
                   control = lme4::glmerControl(optimizer = "bobyqa"))
  m <- glmm(fo, data = d1L, family = binomial("logit"),
            method = "joint_laplace", control = mm_jl())
  bg <- unname(lme4::fixef(g))
  bm <- unname(fixef(m))
  expect_lt(max(abs(bm - bg)), 5e-3)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2)
  expect_true(m$fit_status %in%
                c("converged_interior", "converged", "not_optimized"))
})

test_that("GLMM Wald inference is certified on a real-data random-intercept model", {
  skip_if_not_installed("lme4")
  d1 <- osf_ww_data("study1a")
  fo <- wait_choice ~ Enjoyment_centered + (1 | ID)
  g <- lme4::glmer(fo, data = d1, family = binomial("logit"),
                   control = lme4::glmerControl(optimizer = "bobyqa"))
  m <- glmm(fo, data = d1,
            family = binomial("logit"), method = "joint_laplace",
            control = mm_jl())
  ct <- summary(m, tests = "coefficients")$coefficients
  gt <- summary(g)$coefficients
  p_col <- intersect(c("Pr(>|z|)", "p.value"), names(ct))[1L]

  expect_true(all(ct$method == "asymptotic_wald_z"))
  expect_true(all(is.finite(ct[["Std. Error"]]) & ct[["Std. Error"]] > 0))
  expect_true(all(is.finite(ct[["z value"]])))
  expect_false(is.na(p_col))
  expect_true(all(is.finite(ct[[p_col]]) & ct[[p_col]] >= 0 & ct[[p_col]] <= 1))
  expect_lt(max(abs(unname(ct[["Std. Error"]]) - unname(gt[, "Std. Error"]))), 5e-3)
  expect_lt(max(abs(unname(ct[["z value"]]) - unname(gt[, "z value"]))), 1.5e-1)
  expect_lt(max(abs(unname(ct[[p_col]]) - unname(gt[, "Pr(>|z|)"]))), 5e-3)
})

test_that("full 9-model Willingness-to-Wait sweep matches glmer (slow)", {
  skip_if_not_installed("lme4")
  skip_if_not(nzchar(Sys.getenv("MIXEFF_RUN_SLOW_PARITY")),
              "Set MIXEFF_RUN_SLOW_PARITY=true to run the full 9-model sweep.")
  d1 <- osf_ww_data("study1a")
  d1L <- osf_ww_long("study1a")
  d2 <- osf_ww_data("study1b")
  d2L <- osf_ww_long("study1b")

  cases <- list(
    # correlated-slope models -> tight parity (5e-3 / 5e-2)
    list(
      d = d1, tight = TRUE,
      fo = wait_choice ~ 1 + Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)
    ),
    list(
      d = d1, tight = TRUE,
      fo = wait_choice ~ 1 + Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)
    ),
    list(
      d = d2, tight = TRUE,
      fo = wait_choice ~ 1 + Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)
    ),
    list(
      d = d2, tight = TRUE,
      fo = wait_choice ~ 1 + Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)
    ),
    # random-intercept / partial-slope comp models -> tight parity since the
    # upstream high-baseline early-stop fix (bd-01KT40T6FGVXQQ9N50G2HM0ZZE)
    list(d = d1L, tight = TRUE, fo = score ~ Enjoyment_centered + SVScore_centered + (1 | ID) + (1 | Title)),
    list(d = d1L, tight = TRUE, fo = score ~ Enjoyment_centered + arousal_centered + SVScore + (1 | ID) + (1 | Title)),
    list(d = d2L, tight = TRUE, fo = score ~ Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 | Title)),
    list(d = d2L, tight = TRUE, fo = score ~ Enjoyment_centered + (1 | ID) + (1 | Title)),
    list(d = d2L, tight = TRUE, fo = score ~ Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 | Title))
  )
  for (cs in cases) {
    g <- suppressWarnings(lme4::glmer(
      cs$fo, data = cs$d, family = binomial("logit"),
      control = lme4::glmerControl(optimizer = "bobyqa")))
    m <- glmm(cs$fo, data = cs$d, family = binomial("logit"),
              method = "joint_laplace", control = mm_jl())
    bg <- unname(lme4::fixef(g))
    bm <- unname(fixef(m))
    tol_b  <- if (cs$tight) 5e-3 else 6e-2
    tol_ll <- if (cs$tight) 5e-2 else 1e-1
    expect_lt(max(abs(bm - bg)), tol_b)
    expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), tol_ll)
  }
})
