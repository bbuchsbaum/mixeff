# Snapshot tests for the random-effects guidance surface (PRD §9.5.7, §9.5.4,
# §9.5.5, §11). These complement test-explain-model.R, which uses targeted
# `expect_match` against individual phrases. The snapshots here are the
# contract: drift in tone, wording, or layout will surface as a snapshot diff
# in one place. R authors no per-block English; Rust does (PRD §9.6 / R9).
#
# Update protocol:
#   1. If the wording change is *intentional* and the new output is correct,
#      run `testthat::snapshot_review("explain-model-snapshots")` and accept.
#   2. If the wording change is *accidental*, fix the source (usually the Rust
#      `RandomTermCard.english` field or a Rust diagnostic message) and re-run.

mk_snapshot_design <- function() {
  set.seed(3L)
  df <- expand.grid(
    s = factor(seq_len(6)),
    i = factor(seq_len(4)),
    b = factor(seq_len(2))
  )
  df$a <- factor(rep(seq_len(3), length.out = nrow(df)))
  df$t <- rep(seq_len(4), length.out = nrow(df))
  df$y <- rnorm(nrow(df))
  df
}

mk_refusal_design <- function() {
  # `between` does not vary within `g` → requesting (1 + between | g) must
  # surface a structural_refusal diagnostic with "possible repairs" wording.
  set.seed(11L)
  data.frame(
    y       = rnorm(8),
    between = rep(seq_len(4), each = 2),
    g       = factor(rep(seq_len(4), each = 2))
  )
}

explain_text <- function(formula, data = mk_snapshot_design()) {
  explain_model(compile_model(formula, data))$text
}

#-- §9.5.7 — eight syntax patterns --------------------------------------------

test_that("snapshot: explain_model() for the eight §9.5.7 syntax patterns", {
  expect_snapshot(cat(explain_text(y ~ t + (1 | s))))                    # punt
  expect_snapshot(cat(explain_text(y ~ t + (0 + t | s))))                # slope-only
  expect_snapshot(cat(explain_text(y ~ t + (1 + t | s))))                # correlated full
  expect_snapshot(cat(explain_text(y ~ t + (1 + t || s))))               # double-bar
  expect_snapshot(cat(explain_text(y ~ t + (1 | s) + (0 + t | s))))      # split-block
  expect_snapshot(cat(explain_text(y ~ t + (1 | a / b))))                # nested
  expect_snapshot(cat(explain_text(y ~ t + (1 | s:i))))                  # interaction
  expect_snapshot(cat(explain_text(y ~ t + (1 | s) + (1 | i))))          # crossed
})

#-- §9.5.4 — three kinds of help ----------------------------------------------

test_that("snapshot: three kinds of help register cleanly", {
  # Unmodeled-but-possible (Design note): scope_note on the intercept-only
  # term for `s` (since `t` varies within `s`).
  expect_snapshot(cat(explain_text(y ~ t + (1 | s))))

  # Structural impossibility: between-group `between` requested as a random
  # slope. Must produce 'Possible repairs, not applied automatically:' and
  # name the unsupported variable + group.
  expect_snapshot(cat(explain_text(y ~ between + (1 + between | g),
                                   mk_refusal_design())))
})

#-- §9.5.5 — refusals use 'Possible repairs', never 'suggested model' --------

test_that("snapshot: structural_refusal renders 'Possible repairs' wording", {
  refusal_audit <- audit_design(compile_model(
    y ~ between + (1 + between | g),
    mk_refusal_design()
  ))
  expect_snapshot(cat(refusal_audit$text))
})

#-- §9.5.3 — random_options() map (rung 0 first-class, no ranking) ----------

test_that("snapshot: random_options() prints rung 0 first-class, no ranking", {
  spec <- compile_model(y ~ t + (1 | s), mk_snapshot_design())
  expect_snapshot(print(random_options(spec, "s", "t")))
})

#-- compare_covariance() layout — no recommended column ---------------------

test_that("snapshot: compare_covariance() layout", {
  spec <- compile_model(y ~ t + (1 + t | s), mk_snapshot_design())
  expect_snapshot(print(compare_covariance(spec)))
})

#-- DiagnosticCode pedagogical variants ---------------------------------------

test_that("snapshot: pedagogical DiagnosticCode variants round-trip from Rust", {
  # `(1 | a/b)` should emit syntax_expansion; (1 + t || s) emits
  # covariance_assumption; intercept-only with within-group fixed effect emits
  # scope_note.
  scope <- diagnostics(compile_model(y ~ t + (1 | s), mk_snapshot_design()))
  cov   <- diagnostics(compile_model(y ~ t + (1 + t || s), mk_snapshot_design()))
  expand <- diagnostics(compile_model(y ~ t + (1 | a / b), mk_snapshot_design()))

  expect_snapshot(print(scope))
  expect_snapshot(print(cov))
  expect_snapshot(print(expand))
})
