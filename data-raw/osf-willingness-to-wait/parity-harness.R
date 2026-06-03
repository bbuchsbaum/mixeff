#!/usr/bin/env Rscript
# Offline glmer-vs-mixeff parity harness on the committed Willingness-to-Wait
# fixtures. Reproduces the in-the-wild findings recorded on mote
# bd-01KT3ZRCKWRZQFA4W7TXTGWAZ0 (and upstream bd-01KT40T6FGVXQQ9N50G2HM0ZZE).
#   Rscript data-raw/osf-willingness-to-wait/parity-harness.R

suppressMessages({
  library(lme4)
  library(mixeff)
})

fx <- function(which) {
  f <- file.path("tests", "fixtures",
                 sprintf("osf_willingness_to_wait_%s.csv", which))
  if (!file.exists(f)) {
    stop("fixture not found: ", f, " (run reconstruct.R first)")
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  d$Enjoyment <- suppressWarnings(as.numeric(as.character(d$Enjoyment)))
  d$Enjoyment_centered <- d$Enjoyment - mean(d$Enjoyment, na.rm = TRUE)
  d$arousal <- suppressWarnings(as.numeric(as.character(d$arousal)))
  d$arousal_centered <- d$arousal - mean(d$arousal, na.rm = TRUE)
  d
}
gather_qc <- function(d) {
  a <- d; a$score <- d$Q1_correct
  b <- d; b$score <- d$Q2_correct
  out <- rbind(a, b)
  out$SVScore_centered <- out$SVScore - mean(out$SVScore, na.rm = TRUE)
  out$arousal_centered <- out$arousal - mean(out$arousal, na.rm = TRUE)
  out
}

d1  <- fx("study1a"); d1L <- gather_qc(d1)
d2  <- fx("study1b"); d2L <- gather_qc(d2)
ctl <- mm_control(verbose = -1, max_feval = 100000L)
bin <- binomial("logit")

models <- list(
  list(id = "wait_enj_1a",        d = d1,  fo = wait_choice ~ 1 + Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)),
  list(id = "comp_1a",            d = d1L, fo = score ~ Enjoyment_centered + SVScore_centered + (1 | ID) + (1 | Title)),
  list(id = "wait_enj_arous_1a",  d = d1,  fo = wait_choice ~ 1 + Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)),
  list(id = "comp_arous_1a",      d = d1L, fo = score ~ Enjoyment_centered + arousal_centered + SVScore + (1 | ID) + (1 | Title)),
  list(id = "wait_enj_1b",        d = d2,  fo = wait_choice ~ 1 + Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)),
  list(id = "comp_1b",            d = d2L, fo = score ~ Enjoyment_centered + (1 + Enjoyment_centered | ID) + (1 | Title)),
  list(id = "comp_1b_2",          d = d2L, fo = score ~ Enjoyment_centered + (1 | ID) + (1 | Title)),
  list(id = "wait_enj_arous_1b",  d = d2,  fo = wait_choice ~ 1 + Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 + Enjoyment_centered | Title)),
  list(id = "comp_arous_1b",      d = d2L, fo = score ~ Enjoyment_centered + arousal_centered + (1 + Enjoyment_centered | ID) + (1 | Title))
)

cat(sprintf("%-20s %8s %10s %10s %6s\n",
            "model", "dlogLik", "max|dFix|", "status", "iters"))
for (m in models) {
  g <- suppressWarnings(glmer(m$fo, data = m$d, family = bin,
                              control = glmerControl(optimizer = "bobyqa")))
  fit <- glmm(m$fo, data = m$d, family = bin, method = "joint_laplace",
              control = ctl)
  bg <- unname(lme4::fixef(g)); bm <- unname(fixef(fit))
  cert <- fit$artifact$optimizer_certificate
  cat(sprintf("%-20s %+8.1e %10.1e %10s %6s\n",
              m$id,
              as.numeric(logLik(fit)) - as.numeric(logLik(g)),
              max(abs(bm - bg)),
              fit$fit_status,
              as.character(cert$iterations %||% NA)))
}
cat("\nwait_* (correlated slopes) match glmer to ~1e-3; comp_* (random",
    "intercept, high baseline) land ~0.01-0.05 logLik short -- upstream",
    "bd-01KT40T6FGVXQQ9N50G2HM0ZZE.\n")
