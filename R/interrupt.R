# Internal helper for tests/testthat/test-interrupts.R only.
#
# Wraps `wrap__mm_interrupt_demo`, the Phase 0 smoke test that the
# `R_CheckUserInterrupt` FFI binding is linked correctly. Long-running
# Rust loops (PLS / PIRLS in Phase 1+) will use the same hook from inside
# their own iteration step. Not exported because it has no end-user purpose.
#
# Returns `iters` on clean completion. Pressing Ctrl-C while the loop runs
# unwinds via R's error handler.

mm_interrupt_demo <- function(iters) {
  iters <- as.integer(iters)
  if (length(iters) != 1L || is.na(iters)) {
    stop("`iters` must be a single non-NA integer.")
  }
  .Call(wrap__mm_interrupt_demo, iters)
}
