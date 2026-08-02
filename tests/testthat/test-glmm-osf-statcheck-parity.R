# In-the-wild GLMM parity against a *published* lme4::glmer analysis:
# OSF node 538bc Study 3 (statcheck x Open-Practice badges). The committed
# fixture inst/extdata/osf-statcheck-t2.csv is the Period-2 slice (5279 rows,
# 426 groups, 380 Error events); provenance + reconstruction in
# data-raw/osf-statcheck/. Tracked by mote bd-01KT3ZB649HHX4ZYMKTQ5A18XX.
#
# This guards three things at once:
#   (1) joint_laplace tracks glmer on well-conditioned (centered) real-data
#       models  -> assertions below;
#   (2) certified joint-laplace active-Hessian Wald rows match glmer SE/z/p on
#       the same well-conditioned models;
#   (3) the remaining known raw-Year gap stays documented, not silently
#       regressed:
#         E = on the ill-scaled raw-Year design the joint route does not
#             certify and returns a labelled fast-PIRLS fallback whose typed
#             estimator_substitution record is asserted below. The historical
#             non-finite-objective-as-converged bug (bd-01KT3Z64AY45NHA5144G2ZBMSY)
#             is FIXED at the pinned engine and guarded by a live finiteness
#             assertion. Offset invariance is gated on the substitution record
#             (evidence-gated skip, not unconditional) and flips on when the
#             joint frontier work lands (bd-01KX66CNV2QJ5SKPXMTX1BW8SQ).

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
  # paper reported b = 0.7958, Z = 1.825, p = .0679. Loose tolerance by
  # design: nAGQ = 0 on deliberately ill-scaled data is the most
  # optimizer-sensitive glmer path there is, and this line only anchors
  # that we reproduce the published model, not mixeff behavior.
  expect_equal(unname(b["Estimate"]), 0.7958, tolerance = 1e-2)
  expect_equal(unname(b["z value"]), 1.825,  tolerance = 1e-2)
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
    bg <- unname(lme4::fixef(g))
    bm <- unname(fixef(m))
    expect_equal(length(bm), length(bg))
    expect_lt(max(abs(bm - bg)), 5e-3)                                  # estimates
    expect_lt(abs(as.numeric(logLik(m)) - as.numeric(logLik(g))), 5e-2) # logLik
  }
})

test_that("GLMM Wald standard errors match glmer on OSF statcheck", {
  skip_if_not_installed("lme4")
  d <- osf_statcheck_data()
  fo <- Error ~ OpenPractice + (1 | Source)
  g <- lme4::glmer(fo, data = d, family = binomial("logit"))
  m <- glmm(fo, data = d,
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
  expect_equal(unname(ct[["Std. Error"]]), unname(gt[, "Std. Error"]), tolerance = 5e-3)
  expect_equal(unname(ct[["z value"]]), unname(gt[, "z value"]), tolerance = 5e-2)
  expect_lt(max(abs(unname(ct[[p_col]]) - unname(gt[, "Pr(>|z|)"]))), 5e-3)
})

# One raw-Year joint fit shared by the two tests below (it costs ~a minute;
# the second test must not pay for it twice).
osf_raw_joint_fit <- local({
  cache <- NULL
  function(d) {
    if (is.null(cache)) {
      cache <<- glmm(Error ~ OpenPractice * Year + (1 | Source), data = d,
                     family = binomial("logit"), method = "joint_laplace",
                     control = mm_jl())
    }
    cache
  }
})

test_that("raw-Year joint_laplace substitution is typed, finite, and visible", {
  skip_on_cran() # minute-class real-data fit
  d <- osf_statcheck_data()
  mr <- osf_raw_joint_fit(d)

  # The historical bug (upstream bd-01KT3Z64AY45NHA5144G2ZBMSY, fixed at
  # engine 9b3b4ad, vendored since pin 4a2abb3): this fit reported
  # converged_interior with a NON-FINITE objective. The finiteness half is
  # asserted live; a regression here is a stop-ship convergence-honesty bug.
  expect_true(is.finite(as.numeric(logLik(mr))))

  # What remains on the current pin: the joint route does not certify on
  # this ill-scaled design and returns a labelled fast-PIRLS fallback
  # (upstream bd-01KYW1VTQ3YR16SQP29J092QKG, fixed at engine f82c646 by the
  # typed record asserted here). The substitution must be machine-readable
  # on the public certificate surface and visible in summary() despite the
  # fallback's converged_* status -- no silent surgery.
  sub <- optimizer_certificate(mr)$raw$estimator_substitution
  expect_false(is.null(sub))
  expect_identical(as.character(sub$requested_method), "joint_laplace")
  expect_identical(as.character(sub$effective_method), "fast_pirls_profiled")
  expect_true(nzchar(as.character(sub$requested_fit_status %||% "")))

  out <- capture.output(print(summary(mr)))
  expect_true(any(grepl("estimator substitution", out, fixed = TRUE)))
})

test_that("raw-Year joint_laplace is offset-invariant once the joint route certifies", {
  skip_on_cran() # minute-class real-data fits
  d <- osf_statcheck_data()
  mr <- osf_raw_joint_fit(d)
  # Typed, evidence-gated skip (NOT unconditional): while the engine still
  # substitutes fast-PIRLS on this design, the fallback optimum is known to
  # sit off the offset-invariant MLE (interaction 0.7959 vs ~0.8538; ledger
  # case osf_statcheck_t2_error_openpractice_year_raw). The moment the joint
  # route certifies this fit (upstream frontier bead
  # bd-01KX66CNV2QJ5SKPXMTX1BW8SQ: autoscaling/multi-start), the skip
  # predicate flips off and the invariance assertions run unedited.
  skip_if(
    !is.null(optimizer_certificate(mr)$raw$estimator_substitution),
    paste0(
      "joint route substitutes fast-PIRLS on the raw-Year design; offset ",
      "invariance is asserted only for a certified joint fit (upstream ",
      "frontier bead bd-01KX66CNV2QJ5SKPXMTX1BW8SQ)."
    )
  )
  mc <- glmm(Error ~ OpenPractice * cYear + (1 | Source), data = d,
             family = binomial("logit"), method = "joint_laplace", control = mm_jl())
  expect_equal(ix(fixef(mr)), ix(fixef(mc)), tolerance = 5e-3)
})
