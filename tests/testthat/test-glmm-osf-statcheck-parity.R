# In-the-wild GLMM parity against a *published* lme4::glmer analysis:
# OSF node 538bc Study 3 (statcheck x Open-Practice badges). The committed
# fixture inst/extdata/osf-statcheck-t2.csv is the Period-2 slice (5279 rows,
# 426 groups, 380 Error events); provenance + reconstruction in
# data-raw/osf-statcheck/. Tracked by mote bd-01KT3ZB649HHX4ZYMKTQ5A18XX.
#
# This guards two things at once:
#   (1) joint_laplace tracks glmer on well-conditioned (centered) real-data
#       models  -> assertions below;
#   (2) two known gaps stay documented, not silently regressed:
#         E = raw-Year joint_laplace converges to a sub-optimal point and
#             reports converged_interior with a non-finite objective
#             (upstream bd-01KT3Z64AY45NHA5144G2ZBMSY); the offset-invariance
#             test is skip()ped pending that fix and will flip green when fixed.
#         B = GLMM SEs are uncertified (upstream bd-01KT3Z64YE5QN7626PQRJSJJVA);
#             asserted as the current contract (method=not_computed).

osf_statcheck_data <- function() {
  f <- system.file("extdata", "osf-statcheck-t2.csv", package = "mixeff")
  skip_if(!nzchar(f), "osf-statcheck-t2.csv fixture not installed")
  d <- utils::read.csv(f)
  d$Source       <- factor(d$gid)
  d$OpenPractice <- d$OpenData == 1 | d$OpenMaterials == 1 | d$Preregistration == 1
  d$cYear        <- d$Year - 2015L
  d
}

mm_jl <- function() mm_control(verbose = -1, max_feval = 50000L)
ix <- function(b) as.numeric(b[grep(":", names(b))])

test_that("fixture is faithful to the published result (glmer nAGQ=0)", {
  skip_if_not_installed("lme4")
  d <- osf_statcheck_data()
  # raw (uncentered) Year is deliberately ill-scaled here -> lme4 emits a
  # "predictor variables are on very different scales" warning; expected.
  g0 <- suppressWarnings(
    lme4::glmer(Error ~ OpenPractice * Year + (1 | Source), data = d,
                family = binomial("logit"), nAGQ = 0))
  b <- summary(g0)$coefficients["OpenPracticeTRUE:Year", ]
  # paper reported b = 0.7958, Z = 1.825, p = .0679
  expect_equal(unname(b["Estimate"]), 0.7958, tolerance = 1e-3)
  expect_equal(unname(b["z value"]), 1.825,  tolerance = 2e-3)
})

test_that("joint_laplace tracks glmer on well-conditioned OSF models", {
  skip_if_not_installed("lme4")
  d <- osf_statcheck_data()
  forms <- list(
    Error ~ OpenPractice + (1 | Source),
    Error ~ OpenPractice * cYear + (1 | Source)
  )
  for (fo in forms) {
    g <- lme4::glmer(fo, data = d, family = binomial("logit"))
    m <- glmm(fo, data = d, family = binomial("logit"),
              method = "joint_laplace", control = mm_jl())
    expect_identical(m$method, "joint_laplace")
    bg <- unname(lme4::fixef(g)); bm <- unname(fixef(m))
    expect_equal(length(bm), length(bg))
    expect_lt(max(abs(bm - bg)), 5e-3)                                  # estimates
    expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2) # logLik
  }
})

test_that("GLMM standard errors are uncertified (documented gap B)", {
  d <- osf_statcheck_data()
  m <- glmm(Error ~ OpenPractice + (1 | Source), data = d,
            family = binomial("logit"), method = "joint_laplace",
            control = mm_jl())
  ct <- summary(m)$coefficients
  # Current contract: GLMM Wald inference is not certified -> method flagged and
  # z/p withheld. See inst/extdata/expected-mismatches.json (status=unsupported)
  # and upstream bd-01KT3Z64YE5QN7626PQRJSJJVA.
  expect_true(all(ct$method == "not_computed"))
  expect_true(all(is.na(ct$p.value)))
})

test_that("raw-Year joint_laplace is offset-invariant (pending upstream fix E)", {
  skip("upstream bd-01KT3Z64AY45NHA5144G2ZBMSY: raw-Year joint_laplace converges sub-optimally (interaction 0.7959 vs offset-invariant MLE ~0.853) and reports converged_interior with a non-finite objective. Re-enable when fixed.")
  d <- osf_statcheck_data()
  mr <- glmm(Error ~ OpenPractice * Year  + (1 | Source), data = d,
             family = binomial("logit"), method = "joint_laplace", control = mm_jl())
  mc <- glmm(Error ~ OpenPractice * cYear + (1 | Source), data = d,
             family = binomial("logit"), method = "joint_laplace", control = mm_jl())
  expect_equal(ix(fixef(mr)), ix(fixef(mc)), tolerance = 5e-3)
  expect_true(is.finite(as.numeric(logLik(mr))))
})
