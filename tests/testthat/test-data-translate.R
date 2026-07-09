test_that("mm_translate_data accepts numeric / integer / logical / factor / character", {
  df <- data.frame(
    y    = c(1.1, 2.2, 3.3, 4.4),
    n    = 1:4,
    flag = c(TRUE, FALSE, TRUE, FALSE),
    f    = factor(c("a", "b", "a", "b"), levels = c("b", "a")),
    s    = c("x", "y", "x", "z"),
    stringsAsFactors = FALSE
  )

  out <- mixeff:::mm_translate_data(df)

  expect_identical(out$column_order, c("y", "n", "flag", "f", "s"))

  expect_named(out$numeric_columns, c("y", "n", "flag"))
  expect_identical(out$numeric_columns$y, c(1.1, 2.2, 3.3, 4.4))
  expect_identical(out$numeric_columns$n, c(1, 2, 3, 4))
  expect_identical(out$numeric_columns$flag, c(1, 0, 1, 0))

  expect_named(out$categorical_values, c("f", "s"))
  expect_named(out$categorical_levels, c("f", "s"))
  expect_identical(out$categorical_values$f, c("a", "b", "a", "b"))
  # Factor levels are preserved canonically (not first-appearance).
  expect_identical(out$categorical_levels$f, c("b", "a"))
  # Character levels are first-appearance.
  expect_identical(out$categorical_levels$s, c("x", "y", "z"))
  # No ordered factors present, so contr.poly coding applies to nothing.
  expect_identical(out$categorical_ordered, character(0))
})

test_that("mm_translate_data flags ordered factors for contr.poly coding", {
  df <- data.frame(
    y = 1:4,
    o = factor(c("lo", "hi", "mid", "hi"),
               levels = c("lo", "mid", "hi"), ordered = TRUE),
    u = factor(c("a", "b", "a", "b"))
  )

  out <- mixeff:::mm_translate_data(df)

  # The ordered factor is named in `categorical_ordered`; the unordered one is
  # not. Both still appear in the value/level lists (encoding is chosen in Rust).
  expect_identical(out$categorical_ordered, "o")
  expect_true(all(c("o", "u") %in% names(out$categorical_values)))
  expect_identical(out$categorical_levels$o, c("lo", "mid", "hi"))
})

test_that("mm_translate_data refuses ordered factors under a non-poly ordered option", {
  df <- data.frame(
    y = 1:4,
    o = factor(c("lo", "hi", "mid", "hi"),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  old <- options(contrasts = c(unordered = "contr.treatment",
                               ordered = "contr.treatment"))
  err <- tryCatch(mixeff:::mm_translate_data(df), error = function(e) e)
  options(old)
  expect_s3_class(err, "mm_arg_error")
  expect_match(conditionMessage(err), "contr.poly", fixed = TRUE)
})

test_that("mm_translate_data honours the standard UNNAMED contr.poly option", {
  df <- data.frame(
    y = 1:4,
    o = factor(c("lo", "hi", "mid", "hi"),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  # R's default option form is unnamed; ordered coding is resolved positionally
  # (element 2). Reading it by name would throw and abort a correct fit.
  old <- options(contrasts = c("contr.treatment", "contr.poly"))
  out <- tryCatch(mixeff:::mm_translate_data(df), error = function(e) e)
  options(old)
  expect_false(inherits(out, "error"))
  expect_identical(out$categorical_ordered, "o")
})

test_that("mm_translate_data refuses ordered factors under an unnamed non-poly option", {
  df <- data.frame(
    y = 1:4,
    o = factor(c("lo", "hi", "mid", "hi"),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  old <- options(contrasts = c("contr.treatment", "contr.treatment"))
  err <- tryCatch(mixeff:::mm_translate_data(df), error = function(e) e)
  options(old)
  expect_s3_class(err, "mm_arg_error")
  expect_match(conditionMessage(err), "contr.poly", fixed = TRUE)
})

test_that("mm_translate_data refuses ordered factors with a non-poly contrasts attribute", {
  df <- data.frame(
    y = 1:4,
    o = factor(c("lo", "hi", "mid", "hi"),
               levels = c("lo", "mid", "hi"), ordered = TRUE)
  )
  stats::contrasts(df$o) <- stats::contr.helmert(3)
  expect_error(
    mixeff:::mm_translate_data(df),
    class = "mm_arg_error",
    regexp = "contrasts"
  )
})

test_that("mm_translate_data rejects unsupported column types", {
  df <- data.frame(
    y = 1:3,
    d = as.Date(c("2026-01-01", "2026-01-02", "2026-01-03"))
  )
  expect_error(
    mixeff:::mm_translate_data(df),
    class = "mm_data_error",
    regexp = "Unsupported column type for `d`"
  )
})

test_that("mm_translate_data rejects non-data.frame and empty data.frame inputs", {
  expect_error(
    mixeff:::mm_translate_data(list(y = 1:3)),
    class = "mm_data_error",
    regexp = "must be a data.frame"
  )
  expect_error(
    mixeff:::mm_translate_data(data.frame()),
    class = "mm_data_error",
    regexp = "at least one column"
  )
})

test_that("mm_check_no_na refuses NA in design variables", {
  df <- data.frame(y = c(1, NA, 3), x = c(1, 2, 3))
  expect_error(
    mixeff:::mm_check_no_na(df, c("y", "x")),
    class = "mm_data_error",
    regexp = "Missing values in design variable"
  )
  # No-op on complete cases.
  expect_invisible(mixeff:::mm_check_no_na(df[c(1L, 3L), ], c("y", "x")))
})
