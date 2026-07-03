test_that("inference_options enumerates the six method routes", {
  set.seed(3)
  n_subj <- 8L
  days <- as.numeric(0:4)
  b0 <- rnorm(n_subj, sd = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(subj = factor(i), days = days,
               y = 1 + 0.5 * days + b0[i] + rnorm(length(days), sd = 0.3))
  }))
  fit <- lmm(y ~ days + (1 | subj), d, control = mm_control(verbose = -1))

  opt <- inference_options(fit)
  expect_s3_class(opt, "mm_inference_options")
  expect_true(all(c("asymptotic_wald_z", "satterthwaite", "kenward_roger",
                    "bootstrap", "bootstrap_lrt", "cluster_bootstrap") %in%
                  opt$table$method))
  expect_true("expected_status" %in% names(opt$table))
  expect_true("expected_reliability_reason" %in% names(opt$table))
  expect_true("current" %in% names(opt$table))
  expect_false("recommended" %in% tolower(names(opt$table)))
  expect_equal(sum(opt$table$current), 1L)
})

test_that("inference_options marks satterthwaite/kenward_roger as not_assessed at boundary", {
  set.seed(3)
  n_subj <- 18L
  days <- as.numeric(0:9)
  b0 <- rnorm(n_subj, sd = 30)
  b1 <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(subj = factor(i), days = days,
               rt = 250 + b0[i] + (10 + b1[i]) * days +
                    rnorm(length(days), sd = 20))
  }))
  fit <- lmm(rt ~ days + (1 + days | subj), d, control = mm_control(verbose = -1))
  expect_true(is_singular(fit))

  opt <- inference_options(fit)
  satt <- opt$table[opt$table$method == "satterthwaite", ]
  kr   <- opt$table[opt$table$method == "kenward_roger", ]
  expect_identical(satt$expected_status, "not_assessed")
  expect_identical(kr$expected_status, "not_assessed")
  expect_identical(satt$expected_reliability_reason,
                   "satterthwaite_unavailable_at_boundary")

  wald <- opt$table[opt$table$method == "asymptotic_wald_z", ]
  expect_identical(wald$expected_status, "available")
  # On a singular fit auto resolves to wald, so the wald row reflects the
  # actually-observed reliability_reason from the upstream JSON.
  expect_true(wald$current)
})

test_that("summary() and inference_table() honor an explicit method on a singular fit", {
  set.seed(3)
  n_subj <- 18L
  days <- as.numeric(0:9)
  b0 <- rnorm(n_subj, sd = 30)
  b1 <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(subj = factor(i), days = days,
               rt = 250 + b0[i] + (10 + b1[i]) * days +
                    rnorm(length(days), sd = 20))
  }))
  fit <- lmm(rt ~ days + (1 + days | subj), d, control = mm_control(verbose = -1))
  expect_true(is_singular(fit))

  # auto resolves to asymptotic_wald_z on this singular fit
  auto_inf <- summary(fit, method = "auto")$inference$table
  expect_true(all(auto_inf$method == "asymptotic_wald_z"))
  expect_true(all(auto_inf$status == "available"))

  # explicit satterthwaite request must surface the satterthwaite refusal,
  # not the auto-resolved cached row (regression guard for the latent bug)
  satt_inf <- summary(fit, method = "satterthwaite")$inference$table
  expect_true(all(satt_inf$method == "satterthwaite"))
  expect_true(all(satt_inf$status == "not_assessed"))
  expect_true(all(nzchar(satt_inf$reason)))

  # inference_table() honors the same method argument
  it_satt <- inference_table(fit, method = "satterthwaite")$table
  expect_true(all(it_satt$method == "satterthwaite"))
  expect_true(all(it_satt$status == "not_assessed"))
})

test_that("test_effect(method = 'bootstrap') works on a single-df term (singular fit)", {
  set.seed(3)
  n_subj <- 18L
  days <- as.numeric(0:9)
  b0 <- rnorm(n_subj, sd = 30)
  b1 <- 0.5 * b0 / 30 * 10 + rnorm(n_subj, sd = 0.5)
  d <- do.call(rbind, lapply(seq_len(n_subj), function(i) {
    data.frame(subj = factor(i), days = days,
               rt = 250 + b0[i] + (10 + b1[i]) * days +
                    rnorm(length(days), sd = 20))
  }))
  fit <- lmm(rt ~ days + (1 + days | subj), d, control = mm_control(verbose = -1))

  te <- test_effect(fit, "days", method = "bootstrap",
                    bootstrap = bootstrap_control(nsim = 50, seed = 1))
  expect_identical(te$table$method, "bootstrap")
  expect_identical(te$table$status, "available")
  # The warrant is the engine-authored closed enum, passed through verbatim.
  expect_identical(te$table$reliability_reason,
                   "parametric_bootstrap_monte_carlo")
  expect_identical(te$table$statistic_name, "t")
  expect_true(is.finite(te$table$p_value))
})

test_that("test_effect(method = 'bootstrap') produces a joint F test on a multi-df factor", {
  set.seed(11)
  n_subj <- 12L
  n_per <- 6L
  subj <- factor(rep(seq_len(n_subj), each = n_per))
  cond <- factor(rep(c("A", "B", "C"), times = length(subj) / 3),
                 levels = c("A", "B", "C"))
  b0 <- rnorm(n_subj, sd = 0.4)
  y <- 2 + 0.6 * (cond == "B") + 0.3 * (cond == "C") +
       b0[as.integer(subj)] + rnorm(length(subj), sd = 0.3)
  d <- data.frame(y, cond, subj)
  fit <- lmm(y ~ cond + (1 | subj), d, control = mm_control(verbose = -1))

  te <- test_effect(fit, "cond", method = "bootstrap",
                    bootstrap = bootstrap_control(nsim = 50, seed = 1))
  expect_identical(te$table$method, "bootstrap")
  expect_identical(te$table$status, "available")
  expect_identical(te$table$statistic_name, "f")
  expect_identical(te$table$num_df, 2)  # effective restriction rank for 3-level factor
  expect_true(is.finite(te$table$p_value))
})

test_that("test_effect(method = 'bootstrap_lrt') refuses REML with a stable reason", {
  set.seed(11)
  n_subj <- 12L
  n_per <- 6L
  subj <- factor(rep(seq_len(n_subj), each = n_per))
  cond <- factor(rep(c("A", "B", "C"), times = length(subj) / 3),
                 levels = c("A", "B", "C"))
  b0 <- rnorm(n_subj, sd = 0.4)
  y <- 2 + 0.6 * (cond == "B") + 0.3 * (cond == "C") +
       b0[as.integer(subj)] + rnorm(length(subj), sd = 0.3)
  d <- data.frame(y, cond, subj)
  fit_reml <- lmm(y ~ cond + (1 | subj), d, control = mm_control(verbose = -1))

  te <- test_effect(fit_reml, "cond", method = "bootstrap_lrt",
                    bootstrap = bootstrap_control(nsim = 20, seed = 1))
  expect_identical(te$table$method, "bootstrap_lrt")
  expect_identical(te$table$status, "not_assessed")
  expect_identical(te$table$reason_code, "bootstrap_lrt_requires_ml")
})

test_that("test_effect(method = 'bootstrap_lrt') runs on ML fit and returns a chi-square row", {
  set.seed(11)
  n_subj <- 12L
  n_per <- 6L
  subj <- factor(rep(seq_len(n_subj), each = n_per))
  cond <- factor(rep(c("A", "B", "C"), times = length(subj) / 3),
                 levels = c("A", "B", "C"))
  b0 <- rnorm(n_subj, sd = 0.4)
  y <- 2 + 0.6 * (cond == "B") + 0.3 * (cond == "C") +
       b0[as.integer(subj)] + rnorm(length(subj), sd = 0.3)
  d <- data.frame(y, cond, subj)
  fit_ml <- lmm(y ~ cond + (1 | subj), d, REML = FALSE,
                control = mm_control(verbose = -1))

  te <- test_effect(fit_ml, "cond", method = "bootstrap_lrt",
                    bootstrap = bootstrap_control(nsim = 50, seed = 1))
  expect_identical(te$table$method, "bootstrap_lrt")
  expect_identical(te$table$status, "available")
  expect_identical(te$table$statistic_name, "chi_square")
  # 50 replicates is below the >= 999 moderate threshold: the row must
  # honestly grade itself "low", not the previously hardcoded "moderate".
  expect_identical(te$table$reliability, "low")
  expect_identical(te$table$reliability_reason, "bootstrap_insufficient_replicates")
  expect_true(is.finite(te$table$statistic))
  expect_true(is.finite(te$table$p_value))
  boot <- te$table$details[[1L]]$bootstrap
  expect_equal(boot$successful_replicates, 50)
  expect_true(is.finite(boot$mcse))
  expect_length(boot$replicate_statistics, 50)
})

test_that("inference_options mirrors bootstrap_lrt reliability threshold", {
  set.seed(11)
  n_subj <- 12L
  n_per <- 6L
  subj <- factor(rep(seq_len(n_subj), each = n_per))
  cond <- factor(rep(c("A", "B", "C"), times = length(subj) / 3),
                 levels = c("A", "B", "C"))
  b0 <- rnorm(n_subj, sd = 0.4)
  y <- 2 + 0.6 * (cond == "B") + 0.3 * (cond == "C") +
       b0[as.integer(subj)] + rnorm(length(subj), sd = 0.3)
  d <- data.frame(y, cond, subj)
  fit_ml <- lmm(y ~ cond + (1 | subj), d, REML = FALSE,
                control = mm_control(verbose = -1))

  opt_low <- inference_options(fit_ml, nsim = 50)$table
  lrt_low <- opt_low[opt_low$method == "bootstrap_lrt", ]
  expect_identical(lrt_low$expected_status, "available")
  expect_identical(lrt_low$expected_reliability_reason,
                   "bootstrap_insufficient_replicates")

  opt_mod <- inference_options(fit_ml, nsim = 999)$table
  lrt_mod <- opt_mod[opt_mod$method == "bootstrap_lrt", ]
  expect_identical(lrt_mod$expected_reliability_reason,
                   "bootstrap_monte_carlo_replicates")
})

test_that("test_effect(method = 'cluster_bootstrap') refuses p-values with a stable reason", {
  set.seed(11)
  n_subj <- 12L
  n_per <- 4L
  subj <- factor(rep(seq_len(n_subj), each = n_per))
  x <- rep(seq_len(n_per), n_subj)
  y <- 2 + 0.3 * x + rnorm(n_subj)[as.integer(subj)] + rnorm(length(x), sd = 0.3)
  fit <- lmm(y ~ x + (1 | subj), data.frame(y, x, subj),
             control = mm_control(verbose = -1))

  te <- test_effect(fit, "x", method = "cluster_bootstrap", group = "subj")
  expect_identical(te$table$method, "cluster_bootstrap")
  expect_identical(te$table$status, "not_assessed")
  expect_true(is.na(te$table$p_value))
  expect_identical(te$table$reason_code,
                   "bootstrap_cluster_resample_p_value_unavailable")
  expect_identical(te$table$details[[1L]]$bootstrap$target_kind,
                   "cluster_resample")
  expect_false(te$table$details[[1L]]$bootstrap$p_value_certified)
})

test_that("test_effect(method = 'cluster_bootstrap') requires group for crossed models", {
  set.seed(12)
  d <- expand.grid(subj = factor(seq_len(5)), item = factor(seq_len(4)))
  d$x <- rnorm(nrow(d))
  d$y <- 1 + 0.2 * d$x + rnorm(nrow(d))
  fit <- lmm(y ~ x + (1 | subj) + (1 | item), d,
             control = mm_control(verbose = -1))

  te <- test_effect(fit, "x", method = "cluster_bootstrap")
  expect_identical(te$table$status, "not_assessed")
  expect_identical(te$table$reason_code,
                   "cluster_bootstrap_multifactor_ambiguous")
  expect_error(
    test_effect(fit, "x", method = "cluster_bootstrap", group = "missing"),
    class = "mm_arg_error"
  )
})

test_that("inference_options rejects unknown terms", {
  set.seed(3)
  d <- data.frame(y = rnorm(20), x = rnorm(20),
                  g = factor(rep(1:5, each = 4)))
  fit <- lmm(y ~ x + (1 | g), d, control = mm_control(verbose = -1))
  expect_error(inference_options(fit, term = "definitely_not_a_term"),
               class = "mm_arg_error")
})
