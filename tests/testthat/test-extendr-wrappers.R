# Committed extendr wrappers must cover the #[extendr] surface.
#
# Installation no longer regenerates R/extendr-wrappers.R (doing so rewrote
# the package's own sources at install time and, with a --target set, had
# to execute a freshly built binary -- which fails when host and target
# triples differ). The wrappers are therefore committed artifacts, and this
# test is what keeps them honest: every #[extendr]-exported Rust function
# must have a wrapper, and every wrapper must correspond to a real export.
#
# Source-tree-only: src/rust is not on the installed tests path.

test_that("every #[extendr] function has a committed wrapper", {
  root <- test_path("../..")
  lib_rs <- file.path(root, "src", "rust", "src", "lib.rs")
  wrappers <- file.path(root, "R", "extendr-wrappers.R")
  skip_if_not(file.exists(lib_rs) && file.exists(wrappers),
              "package source is not available")

  rust <- readLines(lib_rs, warn = FALSE)
  # `#[extendr]` (optionally with args) immediately precedes `fn name(`,
  # sometimes with doc comments or other attributes between.
  exported <- character()
  for (i in grep("^\\s*#\\[extendr", rust)) {
    for (j in seq(i + 1L, min(i + 12L, length(rust)))) {
      m <- regmatches(rust[[j]],
                      regexec("^\\s*(pub\\s+)?fn\\s+([A-Za-z0-9_]+)", rust[[j]]))[[1L]]
      if (length(m) >= 3L) {
        exported <- c(exported, m[[3L]])
        break
      }
    }
  }
  exported <- unique(exported)
  expect_true(length(exported) > 20L,
              info = "parsed suspiciously few #[extendr] functions")

  wrapper_text <- readLines(wrappers, warn = FALSE)
  missing <- exported[!vapply(exported, function(fn) {
    any(grepl(sprintf("^%s <- function", fn), wrapper_text)) ||
      any(grepl(sprintf("wrap__%s\\b", fn), wrapper_text))
  }, logical(1))]

  expect_identical(
    missing, character(0),
    info = paste0(
      "these #[extendr] functions have no wrapper -- run ",
      "`Rscript tools/document-rust.R` and commit the result: ",
      paste(missing, collapse = ", ")
    )
  )
})

test_that("no wrapper references a Rust function that no longer exists", {
  root <- test_path("../..")
  lib_rs <- file.path(root, "src", "rust", "src", "lib.rs")
  wrappers <- file.path(root, "R", "extendr-wrappers.R")
  skip_if_not(file.exists(lib_rs) && file.exists(wrappers),
              "package source is not available")

  rust <- paste(readLines(lib_rs, warn = FALSE), collapse = "\n")
  wrapper_text <- readLines(wrappers, warn = FALSE)
  referenced <- unique(unlist(regmatches(
    wrapper_text, gregexpr("wrap__[A-Za-z0-9_]+", wrapper_text)
  )))
  referenced <- sub("^wrap__", "", referenced)
  stale <- referenced[!vapply(referenced, function(fn) {
    grepl(sprintf("fn\\s+%s\\s*\\(", fn), rust)
  }, logical(1))]

  expect_identical(
    stale, character(0),
    info = paste0("wrappers reference removed Rust functions: ",
                  paste(stale, collapse = ", "))
  )
})
