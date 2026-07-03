# Smoke fixture covering Phase 1.A: compile_model() builds an mm_spec and
# audit_design() returns the upstream-rendered audit text for a few worked
# shapes drawn from compiler_contract_v0_prd.md lines 907-922 plus the
# §9.5.7 syntax coverage list. The full snapshot suite (with locked text)
# lands in 1.C alongside explain_model().

mk_sleepstudy_like <- function(n_subjects = 8L, n_per = 4L, seed = 1L) {
  set.seed(seed)
  data.frame(
    Reaction = rnorm(n_subjects * n_per, mean = 250, sd = 30),
    Days     = rep(seq_len(n_per) - 1, n_subjects),
    Subject  = factor(rep(seq_len(n_subjects), each = n_per))
  )
}

mk_crossed <- function(n_subjects = 4L, n_items = 3L, seed = 2L) {
  set.seed(seed)
  grid <- expand.grid(
    Subject = factor(seq_len(n_subjects)),
    Item    = factor(seq_len(n_items))
  )
  grid$y <- rnorm(nrow(grid))
  grid
}

diagnostic_codes <- function(x) {
  vapply(x$diagnostics %||% list(), function(d) d$code, character(1))
}

test_that("compile_model() returns an mm_spec for a sleepstudy-shape design", {
  df <- mk_sleepstudy_like()
  spec <- compile_model(Reaction ~ Days + (1 + Days | Subject), df)

  expect_s3_class(spec, "mm_spec")
  expect_true(is.list(spec$artifact))
  expect_identical(spec$vars, c("Reaction", "Days", "Subject"))
  expect_s3_class(spec$model_frame, "data.frame")
  expect_identical(names(spec$model_frame), spec$vars)

  schema <- spec$artifact$schema
  expect_identical(as.character(schema$schema_name),
                   "mixedmodels.compiled_model_artifact")
  expect_identical(as.character(schema$schema_version), "1")

  # The semantic model recorded one random term grouped by Subject.
  rt <- spec$artifact$semantic_model$random_terms
  expect_length(rt, 1L)
})

test_that("compile_model() handles crossed grouping factors", {
  df <- mk_crossed()
  spec <- compile_model(y ~ 1 + (1 | Subject) + (1 | Item), df)

  expect_s3_class(spec, "mm_spec")
  expect_length(spec$artifact$semantic_model$random_terms, 2L)
})

test_that("audit_design() round-trips upstream text and exposes design_audit", {
  df <- mk_sleepstudy_like()
  spec <- compile_model(Reaction ~ Days + (1 + Days | Subject), df)
  audit <- audit_design(spec)

  expect_s3_class(audit, "mm_audit")
  expect_type(audit$text, "character")
  expect_length(audit$text, 1L)
  expect_match(audit$text, "Audit Summary", fixed = TRUE)
  expect_match(audit$text, "Random Effects",   fixed = TRUE)
  expect_match(audit$text, "Random-Effect Information Budget", fixed = TRUE)

  expect_true(is.list(audit$design_audit))
  expect_true(is.list(audit$report))
  expect_identical(audit$report$schema_name, "mixedmodels.model_audit_report")
  expect_identical(as.character(audit$report$schema_version), "2")
  expect_length(audit$random_term_cards, 1L)
  expect_identical(audit$random_term_cards[[1L]]$schema_name,
                   "mixedmodels.random_term_card")
  expect_identical(audit$random_term_cards[[1L]]$term_id, "r0")
})

test_that("compile/audit surface the upstream pedagogical DiagnosticCode variants", {
  sleep <- mk_sleepstudy_like()
  double_bar <- compile_model(Reaction ~ Days + (1 + Days || Subject), sleep)
  expect_true("covariance_assumption" %in% diagnostic_codes(double_bar$artifact))

  nested <- compile_model(Reaction ~ Days + (1 | Subject / Days), sleep)
  expect_true("syntax_expansion" %in% diagnostic_codes(nested$artifact))

  sparse <- data.frame(
    y       = rnorm(8),
    x       = rep(0:1, 4),
    between = rep(seq_len(4), each = 2),
    g       = factor(rep(seq_len(4), each = 2))
  )

  punt <- compile_model(y ~ x + (1 | g), sparse)
  punt_codes <- diagnostic_codes(audit_design(punt))
  expect_true("support_note" %in% punt_codes)
  expect_true("scope_note" %in% punt_codes)

  refused <- compile_model(y ~ between + (1 + between | g), sparse)
  expect_true("structural_refusal" %in% diagnostic_codes(refused$artifact))
})

# Negative-case tests for scope_note. The upstream generator (audit.rs::
# scope_note_diagnostics) has dedicated Rust tests for the positive trigger
# and split-block suppression; these R tests pin the three remaining filters
# (numeric-only, no-interaction, varies-within-group) so a future Rust change
# cannot regress them silently from the R surface.

test_that("scope_note does not fire for categorical (factor) fixed effects", {
  set.seed(101)
  n_subj <- 8L
  n_per  <- 4L
  # `treat` varies within Subject but is a factor; the upstream generator
  # restricts scope_note candidates to numeric fixed effects.
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    treat   = factor(rep(c("a", "b"), length.out = n_subj * n_per)),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ treat + (1 | Subject), df)
  codes <- diagnostic_codes(audit_design(spec))
  expect_false("scope_note" %in% codes)
})

test_that("scope_note fires for numeric main effects but not for their interaction", {
  set.seed(102)
  n_subj <- 8L
  n_per  <- 4L
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    x       = rnorm(n_subj * n_per),
    z       = rnorm(n_subj * n_per),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ x * z + (1 | Subject), df)
  diags <- audit_design(spec)$diagnostics %||% list()
  scope_terms <- vapply(diags, function(d) {
    if (identical(d$code, "scope_note")) d$payload$fixed_effect else NA_character_
  }, character(1))
  scope_terms <- scope_terms[!is.na(scope_terms)]
  expect_true("x" %in% scope_terms)
  expect_true("z" %in% scope_terms)
  expect_false("x:z" %in% scope_terms)
})

test_that("scope_note does not fire for fixed effects constant within group", {
  set.seed(103)
  n_subj <- 8L
  n_per  <- 4L
  # `dose` is a between-subject covariate: one value per Subject, constant
  # across the within-subject rows. The basis.supported check upstream
  # rejects within-group-constant fixed effects as scope_note candidates.
  dose_per_subj <- rnorm(n_subj)
  df <- data.frame(
    y       = rnorm(n_subj * n_per),
    dose    = rep(dose_per_subj, each = n_per),
    Subject = factor(rep(seq_len(n_subj), each = n_per))
  )
  spec  <- compile_model(y ~ dose + (1 | Subject), df)
  codes <- diagnostic_codes(audit_design(spec))
  expect_false("scope_note" %in% codes)
})

test_that("compile_model() refuses unknown formula variables", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(Reaction ~ Missing + (1 | Subject), df),
    class = "mm_data_error",
    regexp = "not found in `data`.*Missing"
  )
})

test_that("compile_model() refuses NA in design variables", {
  df <- mk_sleepstudy_like()
  df$Days[3L] <- NA_real_
  expect_error(
    compile_model(Reaction ~ Days + (1 | Subject), df),
    class = "mm_data_error",
    regexp = "Missing values"
  )
})

test_that("compile_model() refuses one-sided formulas", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(~ Days + (1 | Subject), df),
    class = "mm_formula_error"
  )
})

test_that("compile_model() refuses non-data.frame data", {
  df <- mk_sleepstudy_like()
  expect_error(
    compile_model(Reaction ~ Days + (1 | Subject),
                  list(Reaction = df$Reaction, Days = df$Days,
                       Subject = df$Subject)),
    class = "mm_data_error"
  )
})

test_that("audit_design() refuses inputs that are not mm_spec / mm_fit", {
  expect_error(audit_design(list(foo = 1)), class = "mm_schema_error")
})

test_that("compile_model() preserves input column order through the FFI", {
  # Re-order the data.frame so Subject comes first; the artifact's
  # requested_formula and semantic model should still describe the
  # original formula. The translator preserves data.frame column order.
  df <- mk_sleepstudy_like()
  df <- df[, c("Subject", "Days", "Reaction")]
  spec <- compile_model(Reaction ~ Days + (1 | Subject), df)
  expect_identical(spec$vars, c("Reaction", "Days", "Subject"))
  # `requested_formula` is the upstream canonical Display form (which makes
  # the implicit intercept explicit), not the original input string.
  expect_identical(spec$artifact$requested_formula,
                   mm_parse_formula(Reaction ~ Days + (1 | Subject)))
})

test_that("a factor inside || emits the double_bar_factor_term diagnostic", {
  # Contract for the || factor-term semantics decision (upstream
  # bd-01KTRQRZKB): mixeff's || fully decorrelates the block, INCLUDING a
  # factor's level contrasts (each treatment-coded contrast gets an
  # independent variance; no within-factor covariances). lme4's || instead
  # leaves factor terms intact with a full within-factor covariance block,
  # so the same formula fits a smaller model family here. The compatibility
  # bridge is this Info diagnostic naming the divergence and the
  # correlated-block rewrite -- a pin bump that drops it should fail here.
  set.seed(3)
  df <- data.frame(
    y = rnorm(120), x = rnorm(120),
    f = factor(rep(c("a", "b"), 60)),
    g = factor(rep(seq_len(10), each = 12))
  )
  spec <- compile_model(y ~ x + f + (1 + f + x || g), df)
  hit <- Filter(function(d) {
    identical(d$code, "covariance_assumption") &&
      identical(d$payload$reason, "double_bar_factor_term")
  }, spec$artifact$diagnostics %||% list())
  expect_length(hit, 1L)
  d <- hit[[1L]]
  expect_identical(d$severity, "info")
  expect_match(d$message, "fully decorrelates factor 'f'")
  expect_identical(d$payload$factor, "f")
  expect_identical(d$payload$group, "g")
  expect_identical(d$payload$correlated_block_equivalent, "(0 + f | g)")

  # numeric-only || terms do not emit it
  spec_num <- compile_model(y ~ x + (1 + x || g), df)
  expect_length(
    Filter(function(d) identical(d$payload$reason, "double_bar_factor_term"),
           spec_num$artifact$diagnostics %||% list()),
    0L
  )

  # the explicit correlated-block expansion (lme4's family) does not either
  spec_exp <- compile_model(y ~ x + f + (1 | g) + (0 + f | g) + (0 + x | g), df)
  expect_length(
    Filter(function(d) identical(d$payload$reason, "double_bar_factor_term"),
           spec_exp$artifact$diagnostics %||% list()),
    0L
  )

  # No error-severity diagnostics on any of these healthy compiles. Guards
  # the upstream split theta-map regression (mixeff-rs bd-01KTVEQFPD, fixed
  # at e7f19a8): pre-fix, every ||-with-factor compile left a spurious
  # severity=error "split random-term missing from optimizer basis" behind.
  error_diags <- function(spec) {
    Filter(function(d) identical(d$severity, "error"),
           spec$artifact$diagnostics %||% list())
  }
  expect_length(error_diags(spec), 0L)
  expect_length(error_diags(spec_num), 0L)
  expect_length(error_diags(spec_exp), 0L)

  # 3+-level factors exercise the multi-contrast Diagonal split map
  df$f3 <- factor(rep(c("a", "b", "c"), 40))
  spec3 <- compile_model(y ~ x + f3 + (1 + f3 || g), df)
  hit3 <- Filter(function(d) identical(d$payload$reason, "double_bar_factor_term"),
                 spec3$artifact$diagnostics %||% list())
  expect_length(hit3, 1L)
  expect_identical(hit3[[1L]]$payload$n_levels, 3L)
  expect_length(error_diags(spec3), 0L)
})

test_that("factor || and explicit correlated expansion fit different covariance families", {
  set.seed(4)
  n_group <- 24L
  df <- expand.grid(
    g = factor(seq_len(n_group)),
    f = factor(c("a", "b", "c")),
    rep = seq_len(3L)
  )
  df$x <- rnorm(nrow(df))
  group_shift <- rnorm(n_group, sd = 0.3)
  f_shift <- matrix(rnorm(n_group * 3L, sd = 0.2), nrow = n_group)
  df$y <- 1 + 0.25 * df$x + group_shift[as.integer(df$g)] +
    f_shift[cbind(as.integer(df$g), as.integer(df$f))] +
    rnorm(nrow(df), sd = 0.5)

  native <- lmm(
    y ~ x + f + (1 + f + x || g),
    df,
    control = mm_control(verbose = -1)
  )
  expanded <- lmm(
    y ~ x + f + (1 | g) + (0 + f | g) + (0 + x | g),
    df,
    control = mm_control(verbose = -1)
  )

  native_theta <- parameterization(native)$table
  expanded_theta <- parameterization(expanded)$table

  native_f <- native_theta[native_theta$source_syntax == "(0 + f | g)", ]
  expanded_f <- expanded_theta[expanded_theta$source_syntax == "(0 + f | g)", ]

  expect_identical(length(native$theta), 4L)
  expect_identical(length(expanded$theta), 8L)

  expect_equal(nrow(native_f), 2L)
  expect_true(all(native_f$covariance_family == "diagonal"))
  expect_true(all(native_f$lambda_row == native_f$lambda_col))
  expect_setequal(native_f$lambda_row_basis, c("f: b", "f: c"))
  expect_false(any(grepl("correlation[", native_f$varcorr_entries, fixed = TRUE)))

  expect_equal(nrow(expanded_f), 6L)
  expect_true(all(expanded_f$covariance_family == "full_cholesky"))
  expect_setequal(expanded_f$lambda_row_basis, c("f: a", "f: b", "f: c"))
  expect_true(any(expanded_f$lambda_row != expanded_f$lambda_col))
  expect_true(any(grepl("correlation[", expanded_f$varcorr_entries, fixed = TRUE)))
})
