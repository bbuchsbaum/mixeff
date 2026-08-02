# Documentation claim guards (Audit1 WI-8.11a).
#
# The audit's "claims to remove" list, enforced mechanically: phrases that
# once carried overclaims must not reappear anywhere in the documentation
# surface. Source-tree-only (vignettes are not on the installed test path
# during R CMD check). NEWS.md is exempt everywhere: it is annotated
# history, and its corrections deliberately quote the retired claims.

mm_doc_files <- function() {
  root <- test_path("../..")
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")),
              "package source is not available")
  c(
    list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$",
               recursive = TRUE, full.names = TRUE),
    file.path(root, "README.md")
  )
}

mm_doc_hits <- function(files, pattern, fixed = TRUE, unless = NULL) {
  # `unless`: a pattern that excuses a line (used for retirement notes that
  # must quote the retired phrase in order to document its retirement).
  hits <- character()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    idx <- grep(pattern, lines, fixed = fixed, ignore.case = TRUE)
    if (!is.null(unless) && length(idx)) {
      idx <- idx[!grepl(unless, lines[idx], ignore.case = TRUE)]
    }
    if (length(idx)) {
      hits <- c(hits, sprintf("%s:%d", basename(f), idx))
    }
  }
  hits
}

test_that("retired overclaims stay retired across vignettes and README", {
  files <- mm_doc_files()
  forbidden <- c(
    "at face value",                        # converged_interior overclaim
    "functionally equivalent",              # blanket-equivalence positioning
    "full replacement",                     # ditto
    "fits always terminate",                # bounded budget != termination proof
    "certified 1.0 surface"                 # "certified" reserved for its contract
  )
  for (phrase in forbidden) {
    hits <- mm_doc_hits(files, phrase)
    expect_true(length(hits) == 0L,
                info = sprintf("forbidden phrase '%s' found at: %s",
                               phrase, paste(hits, collapse = ", ")))
  }
  # The retired enum spelling may appear ONLY on a line that documents its
  # retirement (the glossary's migration note for downstream tooling).
  hits <- mm_doc_hits(files, "interior_converged_well_specified",
                      unless = "retired")
  expect_true(length(hits) == 0L,
              info = sprintf("retired enum used outside a retirement note: %s",
                             paste(hits, collapse = ", ")))
})

test_that("'average treatment effect' appears only in the inference estimand discussion", {
  files <- mm_doc_files()
  offending <- files[basename(files) != "inference.Rmd"]
  hits <- mm_doc_hits(offending, "average treatment effect")
  expect_true(length(hits) == 0L,
              info = paste("causal ATE language outside vignettes/inference.Rmd:",
                           paste(hits, collapse = ", ")))
  # And inside inference.Rmd it appears exactly once: inside the four-way
  # estimand distinction (item 4), never as an identification of a
  # model-based average.
  inf <- files[basename(files) == "inference.Rmd"]
  if (length(inf)) {
    n <- sum(vapply(inf, function(f) {
      length(grep("average treatment effect", readLines(f, warn = FALSE),
                  fixed = TRUE, ignore.case = TRUE))
    }, integer(1)))
    expect_identical(n, 1L)
  }
})

test_that("'drop-in' appears only in explicit negations", {
  files <- mm_doc_files()
  for (f in files) {
    lines <- readLines(f, warn = FALSE)
    idx <- grep("drop-in", lines, fixed = TRUE, ignore.case = TRUE)
    for (i in idx) {
      expect_match(lines[[i]], "not", ignore.case = TRUE,
                   info = sprintf("positive drop-in claim at %s:%d: %s",
                                  basename(f), i, lines[[i]]))
    }
  }
})

test_that("Ctrl-C interruptibility is not claimed in current docs", {
  files <- mm_doc_files()
  hits <- mm_doc_hits(files, "Ctrl-C")
  expect_true(length(hits) == 0L,
              info = paste("Ctrl-C claim outside annotated NEWS history:",
                           paste(hits, collapse = ", ")))
})

test_that("the six-vignette surface is exactly the declared set", {
  root <- test_path("../..")
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")),
              "package source is not available")
  cran_vignettes <- sort(basename(list.files(
    file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = FALSE,
    recursive = FALSE
  )))
  expect_identical(
    cran_vignettes,
    sort(c("mixeff.Rmd", "lmm.Rmd", "glmm.Rmd", "convergence.Rmd",
           "inference.Rmd", "lme4-migration.Rmd"))
  )
  articles <- sort(basename(list.files(
    file.path(root, "vignettes", "articles"), pattern = "\\.Rmd$"
  )))
  expect_identical(
    articles,
    sort(c("benchmarking.Rmd", "inference-method-glossary.Rmd",
           "reproducing-aphantasia.Rmd", "reporting-lmms.Rmd"))
  )
})