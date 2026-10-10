# Free-dispersion GLMMs (Gamma, inverse Gaussian, Gaussian with a non-identity
# link) follow lme4 2.1-0: sigma() = sqrt(phi) with
# phi = deviance / (n - rank([X, Z])), PIRLS working weights divided by phi,
# VarCorr() with useSc = FALSE (absolute SDs, no Residual row) and sc =
# sigma(). The lme4-free checks below always run; the live comparisons need
# lme4 >= 2.1-0 (the stored-reference parity lives in
# test-glmm-family-links.R and friends).

disp_fit <- function(key = "scale_gamma_log", ...) {
  case <- mm_disp_case(key)
  glmm(case$formula, case$data(), family = case$family,
       control = mm_control(verbose = -1), ...)
}

test_that("free-dispersion families are recognised", {
  fd <- function(family, link) {
    mixeff:::mm_glmm_free_dispersion(list(family = list(family = family,
                                                         link = link)))
  }
  expect_true(fd("gamma", "log"))
  expect_true(fd("inverse_gaussian", "inverse"))
  expect_true(fd("gaussian", "log"))
  expect_false(fd("gaussian", "identity"))
  expect_false(fd("poisson", "log"))
  expect_false(fd("binomial", "logit"))
  expect_false(fd("negative_binomial", "log"))
})

test_that("sigma() is sqrt(deviance / (n - rank([X, Z])))", {
  for (key in c("scale_gamma_log", "fl_gaussian_log")) {
    fit <- disp_fit(key)
    dc <- getME(fit, "devcomp")
    q_eff <- dc$dims[["qEff"]]
    X <- model.matrix(fit, type = "fixed")
    Z <- getME(fit, "Z")
    expect_identical(q_eff, as.integer(Matrix::rankMatrix(cbind(X, as.matrix(Z)))))
    phi <- sum(residuals(fit, type = "deviance")^2) / (nobs(fit) - q_eff)
    expect_equal(sigma(fit), sqrt(phi), tolerance = 1e-4, label = key)
    expect_equal(dc$cmp[["sigmaML"]], sigma(fit))
    expect_identical(dc$dims[["dispProfile"]], 1L)
    expect_identical(dc$dims[["maxPhiIter"]], 100L)
    expect_identical(dc$dims[["useSc"]], 1L)
  }
  # Families without a free dispersion: sigma() is 1 and qEff is NA.
  pd <- mm_fl_data(function(x, re) rpois(length(x), exp(0.5 + 0.3 * x + re)))
  pf <- glmm(y ~ x + (1 | g), pd, family = poisson(),
             control = mm_control(verbose = -1))
  expect_identical(sigma(pf), 1)
  expect_true(is.na(getME(pf, "devcomp")$dims[["qEff"]]))
})

test_that("working weights carry 1/phi and VarCorr follows lme4 2.1-0", {
  fit <- disp_fit("scale_gamma_log")
  fam <- family(fit)
  mu <- fitted(fit)
  base <- fam$mu.eta(fam$linkfun(mu))^2 / fam$variance(mu)
  expect_equal(unname(weights(fit, type = "working")),
               unname(base / sigma(fit)^2))
  vc <- VarCorr(fit)
  expect_false(attr(vc, "useSc"))
  expect_equal(attr(vc, "sc"), sigma(fit))
  # Absolute SDs: never multiplied by sigma.
  expect_equal(unname(attr(vc$g, "stddev")), unname(fit$theta),
               tolerance = 1e-6)
  df <- as.data.frame(vc)
  expect_false("Residual" %in% df$grp)
  expect_false(any(grepl("Residual", capture.output(print(vc)))))
  # Conditional variances: the 1/phi weights already scale the factor.
  pv <- attr(ranef(fit, condVar = TRUE)$g, "postVar")
  pls <- mixeff:::mm_pls(fit)
  Lt <- getME(fit, "Lambdat")
  M <- Matrix::solve(pls$L, Matrix::solve(pls$L, Lt, system = "P"),
                     system = "L")
  expect_equal(as.numeric(pv[1, 1, ]), as.numeric(Matrix::colSums(M^2)),
               tolerance = 1e-8)
})

test_that("nAGQ = 0 vcov() is the unscaled RX covariance", {
  fit <- disp_fit("scale_gamma_log", nAGQ = 0)
  RX <- getME(fit, "RX")
  expect_equal(unname(as.matrix(vcov(fit))), unname(chol2inv(RX)),
               tolerance = 1e-3, ignore_attr = TRUE)
})

test_that("lme4-shaped pieces match live lme4 >= 2.1-0", {
  skip_on_cran()
  skip_if_not(mm_disp_lme4_is_2_1(), "needs lme4 >= 2.1-0")
  for (key in c("scale_gamma_log", "scale_ig_log", "fl_gaussian_log")) {
    case <- mm_disp_case(key)
    d <- case$data()
    fit <- glmm(case$formula, d, family = case$family,
                control = mm_control(verbose = -1))
    ref <- suppressWarnings(lme4::glmer(case$formula, d, family = case$family))
    expect_equal(unname(weights(fit, type = "working")),
                 unname(weights(ref, type = "working")), tolerance = 1e-3,
                 label = key)
    dc <- getME(fit, "devcomp")
    dr <- getME(ref, "devcomp")
    expect_identical(dc$dims, dr$dims)
    expect_equal(dc$cmp, dr$cmp[names(dc$cmp)], tolerance = 1e-3, label = key)
    expect_equal(crossprod(getME(fit, "RX")), crossprod(getME(ref, "RX")),
                 tolerance = 1e-3, ignore_attr = TRUE)
    vc <- VarCorr(fit)
    vr <- lme4::VarCorr(ref)
    expect_identical(attr(vc, "useSc"), attr(vr, "useSc"))
    expect_equal(attr(vc, "sc"), attr(vr, "sc"), tolerance = 1e-3)
    expect_equal(as.data.frame(vc)$sdcor, as.data.frame(vr)$sdcor,
                 tolerance = 1e-3)
    expect_equal(rePCA(fit)$g$sdev, lme4::rePCA(ref)$g$sdev, tolerance = 1e-3)
  }
})
