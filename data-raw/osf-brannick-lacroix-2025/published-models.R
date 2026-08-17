# Thin wrapper so reconstruct.R / parity-harness.R keep a stable source path.
# Canonical definitions live in tests/testthat/helper-osf-brannick-lacroix.R
# (data-raw is Rbuildignored; the helper ships with the test suite).

helper <- file.path("tests", "testthat", "helper-osf-brannick-lacroix.R")
if (!file.exists(helper)) {
  helper <- file.path("..", "..", "tests", "testthat",
                      "helper-osf-brannick-lacroix.R")
}
if (!file.exists(helper)) {
  stop("cannot find tests/testthat/helper-osf-brannick-lacroix.R")
}
source(helper, local = FALSE)
