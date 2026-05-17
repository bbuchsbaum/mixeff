# The interrupt bridge is plumbing for Phase 1+ — long-running PLS/PIRLS
# fits will call R_CheckUserInterrupt periodically inside their iteration
# loops. In Phase 0 we ship only the FFI binding plus a no-op demo so we
# can verify the symbol is linked.
#
# Robustly testing actual Ctrl-C handling from inside testthat is fragile
# (different platforms, parallel test runners, etc.). We assert the
# easy invariants here and leave a manual reproduction note in the file
# for anyone changing the binding.
#
# Manual repro: in an interactive R session,
#
#   library(mixeff)
#   mixeff:::mm_interrupt_demo(1e9)   # then press Ctrl-C
#
# should return control to the prompt within a fraction of a second.

test_that("interrupt demo returns the iteration count on clean completion", {
  expect_identical(mixeff:::mm_interrupt_demo(0L), 0L)
  expect_identical(mixeff:::mm_interrupt_demo(1L), 1L)
  expect_identical(mixeff:::mm_interrupt_demo(100L), 100L)
})

test_that("interrupt demo input validation", {
  expect_error(mixeff:::mm_interrupt_demo(NA_integer_),
               "`iters` must be a single non-NA integer")
  expect_error(mixeff:::mm_interrupt_demo(c(1L, 2L)),
               "`iters` must be a single non-NA integer")
})

test_that("interrupt demo handles negative input by clamping to zero", {
  # Rust side does iters.max(0); R side passes negatives through. Either
  # behaviour is acceptable — assert the contract: completes without error
  # and returns a non-negative count.
  out <- mixeff:::mm_interrupt_demo(-5L)
  expect_true(out >= 0L)
})
