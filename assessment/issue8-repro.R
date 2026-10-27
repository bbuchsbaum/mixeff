# Synthetic reproducer from https://github.com/bbuchsbaum/mixeff/issues/8
# Run in a fresh R session against the package build being investigated.
# Fully synthetic data: no source-study observations, identifiers, or fitted values.
library(mixeff)
set.seed(20261027)
d <- expand.grid(subject = factor(seq_len(24)), item = factor(seq_len(20)))
d$x <- runif(nrow(d), -1, 1)
u <- rnorm(24, 0, 0.6)
v <- rnorm(20, 0, 0.4)
eta <- 0.2 + 0.4 * d$x + u[d$subject] + v[d$item]
d$y <- eta + rnorm(nrow(d), 0, 0.7)
d$binary <- rbinom(nrow(d), 1, plogis(eta))
f <- y ~ x + (1 | subject) + (1 | item) # two theta parameters

cat('mixeff:', as.character(packageVersion('mixeff')), '\n')
cat('RemoteSha:', packageDescription('mixeff')$RemoteSha, '\n')
cat('platform:', R.version$platform, '\n')
cat('\nScalar absolute tolerance, documented as allowed:\n')
scalar_result <- tryCatch({
  lmm(f, d, control = mm_control(verbose = -1,
      optCtrl = list(xtol_abs = 1e-8)))
  'success'
}, error = function(e) conditionMessage(e))
cat(scalar_result, '\n')

cat('\nExplicit length-two tolerance vector:\n')
m <- lmm(f, d, control = mm_control(verbose = -1,
    optCtrl = list(xtol_abs = rep(1e-8, 2))))
print(fixef(m))
cat('optimizer status:', optimizer_certificate(m)$raw$status, '\n')

cat('\nJoint-Laplace diagnostics on a small synthetic binary model:\n')
g <- glmm(binary ~ x + (1 | subject) + (1 | item), d,
    family = binomial(), method = 'joint_laplace', nAGQ = 1L,
    inference = 'asymptotic',
    control = mm_control(max_feval = 5000, ftol_abs = 1e-8,
        ftol_rel = 1e-12, xtol_rel = 1e-8, verbose = -1))
cg <- optimizer_certificate(g)$raw
print(cg$status)
print(cg$evidence$hessian)
print(cg$evidence$certification_quality)
vg <- vcov(g)
print(attributes(vg)[c('mm_method', 'mm_status', 'mm_reliability')])
# This compact example does not reproduce the large-model runtime.

print(optimizer_certificate(g))
print(attr(vcov(g), 'mm_notes'))
sessionInfo()
