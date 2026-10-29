# In-the-wild GLMM parity against Brannick & LaCroix (2025),
# Scientific Reports 15:27262. OSF project gf3km / analysis node 7ryjh.
# Committed fixtures are the Cue-Onset slices of 8-25_Database.csv
# (Tables 3–4 nobs). Table 2 / 5 / 6 are not in the corpus: listening
# files are missing and the public event-level file does not recover
# those nobs. Tracked by mote bd-01M0842MM295RFXXN4C1ZDE9CV.
#
# Fast tests check fixture shape and paper nobs. Alerting refits are
# skip_on_cran. Orienting is behind MIXEFF_RUN_SLOW_PARITY.

osf_bl_fixture_path <- function(which) {
  fn <- sprintf("osf_bl_%s.csv.gz", which)
  candidates <- c(
    testthat::test_path("..", "fixtures", fn),
    file.path("tests", "fixtures", fn)
  )
  hit <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(hit), paste(fn, "fixture is unavailable; run",
                            "Rscript data-raw/osf-brannick-lacroix-2025/reconstruct.R"))
  hit
}

osf_bl_read <- function(which) {
  utils::read.csv(gzfile(osf_bl_fixture_path(which)), stringsAsFactors = FALSE)
}

mm_jl <- function() mm_control(verbose = -1, max_feval = 100000L)

test_that("Brannick-LaCroix Cue-Onset fixtures match published nobs and coding", {
  a <- osf_bl_prepare_alerting(osf_bl_read("alerting"))
  o <- osf_bl_prepare_orienting(osf_bl_read("orienting"))

  expect_identical(nrow(a), osf_bl_paper_nobs$alerting)
  expect_identical(nrow(o), osf_bl_paper_nobs$orienting)

  expect_identical(levels(a$Music), c("HAPPY", "SAD", "SILENCE"))
  expect_identical(levels(a$ANT), c("ANT_1", "ANT_2"))
  expect_identical(levels(a$Cue), c("Double", "No"))
  expect_identical(levels(o$Cue), c("Center", "Valid"))
  expect_true(all(a$Blink %in% c(0L, 1L)))
  # Blink rate is ~2% in this 500 ms window; some participants have no
  # blinks here. The paper's no-separation check was on the full ANT series.
  expect_gt(mean(a$Blink), 0)
  expect_lt(mean(a$Blink), 0.1)
  expect_gte(length(unique(a$P)), 50L)
})

test_that("alerting fixture refits the locked glmer Table-3 design", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  a <- osf_bl_prepare_alerting(osf_bl_read("alerting"))
  expect_identical(nrow(a), osf_bl_paper_nobs$alerting)
  g <- suppressWarnings(lme4::glmer(
    osf_bl_forms$alerting, data = a, family = binomial("logit"),
    control = lme4::glmerControl(optimizer = "bobyqa")
  ))
  b <- lme4::fixef(g)
  # Locked 2026-08-17 refit on 8-25_Database.csv (nobs exact vs Table 3).
  # The printed Table 3 vector is NOT recovered (Cue-No 0.618 vs 0.781;
  # several Music terms sign-flipped). See data-raw README.
  expect_equal(unname(b[["(Intercept)"]]), -4.498, tolerance = 5e-3)
  expect_equal(unname(b[["CueNo"]]), 0.618, tolerance = 5e-3)
  expect_gt(unname(b[["CueNo"]]), 0)
})

test_that("joint_laplace tracks glmer on the published alerting GLMM", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  a <- osf_bl_prepare_alerting(osf_bl_read("alerting"))
  fo <- osf_bl_forms$alerting
  g <- suppressWarnings(lme4::glmer(
    fo, data = a, family = binomial("logit"),
    control = lme4::glmerControl(optimizer = "bobyqa")
  ))
  m <- glmm(fo, data = a, family = binomial("logit"),
            method = "joint_laplace", control = mm_jl())
  expect_identical(m$method, "joint_laplace")
  bg <- unname(lme4::fixef(g))
  bm <- unname(fixef(m))
  expect_equal(length(bm), length(bg))
  expect_lt(max(abs(bm - bg)), 5e-3)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2)
})

test_that("joint_laplace tracks glmer on the published orienting GLMM (slow)", {
  skip_on_cran()
  skip_if_not_installed("lme4")
  skip_if_not(mm_run_slow_parity(),
              "Set MIXEFF_RUN_SLOW_PARITY=true to run the orienting Brannick-LaCroix case.")

  o <- osf_bl_prepare_orienting(osf_bl_read("orienting"))
  g <- suppressWarnings(lme4::glmer(
    osf_bl_forms$orienting, data = o, family = binomial("logit"),
    control = lme4::glmerControl(optimizer = "bobyqa")
  ))
  m <- glmm(osf_bl_forms$orienting, data = o, family = binomial("logit"),
            method = "joint_laplace", control = mm_jl())
  expect_lt(max(abs(unname(fixef(m)) - unname(lme4::fixef(g)))), 5e-3)
  expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2)
})
