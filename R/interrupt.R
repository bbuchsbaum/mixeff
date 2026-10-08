# Internal helper for tests/testthat/test-interrupts.R only.
#
# Wraps `wrap__mm_interrupt_demo`, a smoke test of the interrupt bridge.
# The bridge checks for a pending interrupt with R_CheckUserInterrupt run
# inside R_ToplevelExec, so the longjmp never crosses Rust frames; the fit
# loops use the same check through the engine's progress callback. Not
# exported because it has no end-user purpose.
#
# Returns `iters` on clean completion. Pressing Ctrl-C while the loop runs
# returns an `mm_interrupted` error.

mm_interrupt_demo <- function(iters) {
  iters <- as.integer(iters)
  if (length(iters) != 1L || is.na(iters)) {
    stop("`iters` must be a single non-NA integer.")
  }
  out <- tryCatch(.Call(wrap__mm_interrupt_demo, iters),
                  error = function(cnd) cnd)
  if (inherits(out, "condition")) mm_abort_from_bridge(out)
  out
}
