# Pre-CRAN audit — October 2026

A review of mixeff 0.2.0 and of the mixeff-rs engine before the first CRAN
submission. It covers:

- functionality gaps relative to lme4 and its ecosystem;
- correctness of the R layer and of the engine;
- CRAN mechanics;
- performance.

Most findings come from reading the code. Items marked **[repro]** were
reproduced with a probe program or a test.

Paths are relative to each repository. In mixeff, `upstream/` means
`src/rust/upstream/mixeff-rs/`.

Status key: `[ ]` open, `[~]` partial, `[x]` fixed or decided (with the commit,
PR or decision). Engine fixes are in bbuchsbaum/mixeff-rs#7; R fixes in
bbuchsbaum/mixeff#4.

## 1. Release blockers

- [x] **The vendored engine is stale.** _Fixed: re-pinned to mixeff-rs `accd4b1` (rc.5 + fixes) in mixeff `f934dff`._
  - mixeff vendors mixeff-rs `1f3f689` (rc.1, 2026-08-02), which is 42
    commits behind rc.5.
  - Fixes that R users currently do not get:
    - Bootstrapping a `joint_laplace` GLMM refits every replicate with the
      fast PIRLS estimator (695ddc2). On OSF the intercept comes out −2.84
      instead of −3.20.
    - Joint GLMM fits without NLopt stop early, and Ctrl-C during the inner
      PIRLS is swallowed (0fbfd83). The R package builds without NLopt, so
      this is the live path.
    - The interaction grouping key is joined with `_`, so `("s_1","x")` and
      `("s","1_x")` collide in `(1|a:b)` and `(1|a/b)` (909e42f).
    - Wald p-values underflow to 0 beyond |z|≈8.3, and the negative
      binomial deviance loses precision at large θ (695ddc2, 977ef11).
    - Bootstrap refits drop the optimizer controls (897dd25).
    - The performance work: kb07 maximal 2.7 s → 0.39 s, grouseticks about
      7×, kb07 bootstrap 6×.
  - Re-pinning needs:
    - inference-table schema 1.2.0 in the bridge (`KNOWN_SCHEMAS` and the
      parsers);
    - the rc.5 nested-REML fix (977ef11), because the gradient-oracle
      TrustBQ from f636282 had a wrong-optimum regression on PW2;
    - a full re-run of the parity suite.
- [x] **The declared Rust minimum is wrong.** _Fixed: 1.85 everywhere, plus a test against the bundled manifests._
  - `DESCRIPTION` and `src/rust/Cargo.toml` say 1.78.
  - The engine declares `rust-version = "1.85"`.
  - The vendored `indexmap 2.14` and `hashbrown 0.17` need 1.85, and
    `cobyla 1.0.0` and `hashbrown` use `edition = "2024"`.
  - Update `DESCRIPTION`, `src/rust/Cargo.toml`, `cran-comments.md` and
    `tests/testthat/test-msrv.R`.
- [x] **`y ~ 0 + f` silently drops a factor level** **[repro]** _Fixed: mixeff-rs `b0037a4`._
  - The engine always treatment-coded factors, even without an intercept.
  - R full-codes the first factor in that case (`model.c` `modelmatrix`).
  - The R name map keeps the matching subset of columns, so nothing
    warned.
- [x] **`predict(newdata=)` is wrong for non-marginal fixed formulas** _Fixed: mixeff-rs `b0037a4`._
  **[repro]**.
  - `linear/predict.rs` built the newdata design with a separate
    all-treatment builder.
  - It then looked coefficients up by name and treated missing columns as
    zero.
  - `y ~ f + f:h` (that is, `f/h`) was off by up to 0.88; `y ~ 0 + f:h` was
    off by 7.6.
- [x] **A character fixed-effect predictor aborts the fit.** _Fixed: converted to sorted factors in `compile_model()`._
  - `R/data-translate.R:88-91` builds levels in first-appearance order.
  - `R/coef-names.R` `mm_coef_name_map()` only treats `is.factor` columns as
    factors.
  - Result: `mm_schema_error`. This hits every `read.csv` user.
- [x] **A factor whose name prefixes another predictor aborts the fit.** _Fixed: suffix must be a level or contrast name._
  - Examples: `group` with `group_size`, `a` with `age`.
  - `mm_engine_encode_names()` uses `startsWith` without checking that the
    remainder is a level, which produces `group: _size`.
- [x] **Requested contrasts are silently ignored.** _Fixed: all three are refused with `mm_arg_error`._
  - Unordered factors carrying `contrasts<-` matrices and
    `options(contrasts = c("contr.sum", …))` are dropped.
  - The model.frame "contrasts dropped" warning is muffled in
    `R/coef-names.R`.
  - `contr.SAS` is accepted, but the engine uses the first level as the
    reference.
  - All three should be refused, as ordered factors already are.
- [x] **The vignettes overclaim parity.** _Fixed: verb map rewritten; bead ID removed._
  - `vignettes/lme4-migration.Rmd:85-91` says `coef`, `VarCorr`,
    `residuals`, `anova`, `getME`, `isSingular` and `emmeans` are
    "identical" to lme4.
  - Line 102 contains an internal bead ID.
  - `mixeff.Rmd` says `simulate()` does "what you expect", but it is
    LMM-only and ignores weights.
  - The verb map at line 84 says optimizer knobs are "engine-chosen", which
    contradicts `mm_control()`.
- [x] **The default GLMM estimator differs from glmer's (decision needed).**
  _Decided: `glmm()` now defaults to `method = "joint_laplace"`; NB, `nAGQ > 1` and `working_hessian` fall back to the profiled path with a notice._
  - `glmm()` defaults to profiled or fast PIRLS.
  - On OSF that is 5.7 log-likelihood units and about 0.35 on the logit
    scale from lme4. You need `method = "joint_laplace"` to match.
  - Either change the default or say so prominently.

## 2. Correctness: engine (mixeff-rs)

- [x] **High: LMM model comparison never checks random-effect nesting.**
  - `random_effect_terms()` is not implemented for `LinearMixedModel`
    (`model/traits.rs:277`).
  - So `random_effect_comparison` (`stats/lrt.rs:1347`) always returns
    `Same`, and non-nested random structures get a "valid" LRT.
- [x] **High: the boundary LRT applies the 50:50 χ² mixture to correlation
  tests** **[repro]** (`stats/lrt.rs:576-663`).
  - `||` vs `|` halves the p-value.
  - A correlation of 0 is interior; only an added variance lies on the
    boundary.
- [x] **High: binomial integer counts with trials fit to `Ok(NaN)`**
  **[repro]** (`generalized/pirls.rs:642`). _Fixed in mixeff-rs `accd4b1`:
  responses above 1 are refused._
- [x] **High: LMM prior weights are not validated** **[repro]**
  (`linear/mod.rs:1364`). _Fixed in mixeff-rs `accd4b1`: length, finiteness
  and positivity are checked up front; zero weights are refused._
  - A wrong length panics on an assert in `FeMat::reweight`.
  - A zero weight gives a misleading "misspecified" error, and an unfitted
    model's logLik is `-inf`.
  - Decide lme4's zero-weight semantics.
- [x] **Medium: weighted-LMM simulation and parametric bootstrap draw
  ε ~ N(0, σ²)** instead of N(0, σ²/w) (`linear/mod.rs:3753`). _Fixed in
  mixeff-rs `accd4b1`._
- [x] **Medium: `||` with a factor creates a different model from lme4**
  (`linear/mod.rs:4937`).
  - It is not documented in `docs/guide/05_what_is_supported.md`.
  - Separately, `(0+f+h|g)` codes every factor as cell means, so Z is
    rank-deficient.
- [x] **Medium: Gamma and inverse-Gaussian scale conventions are
  inconsistent.**
  - The objective, the residual SD and the vcov rescale use different φ.
  - Verify against lme4, or mark these families experimental.
- [x] **Medium: GLMM `stderror()` and the bootstrap replicate SEs fall back to
  the unscaled working-LMM SEs** (`generalized/mod.rs:1775`, `:465`).
- [x] **Medium: profile intervals are on θ, not the SD or correlation
  scale.** Check how the R side labels them.
- [x] **Low: the bootstrap interval helpers can panic.**
  - `quantile_sorted` and `shortest_interval` panic on an empty finite set.
  - `shortest_cov_int` asserts on the level.
- [x] **Low: a collinear X drops the earlier column (R drops the later one),
  and dropped coefficients report 0.0 rather than NA.**
- [x] **Low: the REML log-determinant skips non-positive diagonals
  silently** (`linear/mod.rs:2204`).
- [x] **Low: term removal is order-sensitive** (`a*b - b:a`).
- [x] **Low: unused grouping levels are kept** (this affects `ngrps`,
  `ranef` and `condVar`).
- [x] **Low: `cluster_resample` merges implicitly nested ids across duplicate
  draws.**
- [x] **Low: the design audit's interaction expansion ignores marginality**,
  so its rank and columns can differ from the fitted design.

## 3. Correctness: R layer (mixeff)

- [x] **`na.exclude` is accepted, but `residuals`/`fitted`/`predict` are not
  padded.**
  - `R/fit-lmm.R:449-461` discards the `na.action` attribute, and
    `R/emmeans.R` reads an attribute that is never set.
- [x] **Unused factor levels are kept after `subset` or NA removal**
  (`R/data-translate.R:83`).
  - Unordered factors get NA coefficients.
  - Ordered factors get a different `contr.poly`.
  - Grouping factors inflate `ngrps`.
- [x] **`predict(newdata)` is fragile** (`R/predict.R:540,688,745`).
  - An integer or character grouping column fails, because it is not
    coerced the way it is at fit time.
  - Any unused Date column fails.
  - An NA in any column fails.
- [x] **`(1|a:b)` labels differ from lme4.**
  - Groups are named `"a & b"` with levels `x_y`.
  - GLMM `predict(newdata)` then aborts (`R/predict.R:427`).
- [x] **`simulate.mm_lmm` on singular crossed fits** draws a zero-variance
  term with another term's variance (`R/simulate.R:224-233`). _Fixed in
  `f934dff`._
- [x] **`simulate.mm_lmm` ignores weights.** _Fixed in `f934dff`._
- [x] **`simulate.mm_lmm`'s `attr(, "seed")` does not follow the
  `stats::simulate` convention.**
- [x] **LMM bootstraps ignore `set.seed()`.**
  - `bootstrap_control(seed = NULL)` leads to `StdRng::from_entropy()`.
  - Draw the seed from R's RNG, as the GLMM path already does.
- [x] **`update(fit, . ~ . + z)` fails, because it reuses the narrowed model
  frame.** Unknown arguments are silently dropped (`R/update.R`). _Fixed: it
  re-evaluates `data` and refuses unknown arguments._
- [x] **`glmm()` argument handling differs from `lmm()`.**
  - [x] `weights = col` and `offset = col` are not evaluated in `data`.
    _Fixed._
  - [x] The `na.action = na.omit` default is never applied.
  _Decided: the default is now `NULL` (NA refused with a typed error naming the columns), as for `lmm()`; `na.omit`/`na.exclude` are opt-in._
  - [x] `subset` is forced before it is refused. _Fixed: it is no longer
    forced._
- [x] **Residuals:**
  - LMM Pearson residuals divide by σ and ignore weights, and
    `scaled = TRUE` divides by σ twice (`R/predict.R:175-181`).
  - GLMM residuals are response-only and default to response residuals;
    lme4 defaults to deviance residuals.
- [x] **Low: `sigma()` returns θ for negative binomial fits** (lme4 returns
  1).
- [x] **Low: GLMM `deviance()` is −2·logLik**, not the sum of squared
  deviance residuals.
- [x] **Low: `drop1` and random-term LRTs re-add the intercept on
  no-intercept models** (`R/compare.R:952`, `R/inference.R:461`).
- [x] **Low: `refit.mm_lmm` with a transformed response** such as `log(y)`
  likely returns the original fit.
- [x] **Low: ordered grouping factors get a dense k×k `contr.poly`**
  (overflow or O(k³) for large k).
- [x] **Low: internal refits pass `weights = fit$weights` through the data
  mask**, so a data column named `fit`, `full` or `object` hijacks it.
- [x] **Low: REML→ML refits in `compare`/`anova` drop the user's
  `mm_control()`.**
- [x] **Low: the exported `fixef`/`ranef`/`VarCorr` mask nlme's**, and
  `fixef.default` only forwards `merMod`.
- [x] **Low: emmeans registration covers `mm_lmm` only** (`R/zzz.R:64`).
- [x] **Low: fits cannot be interrupted.** `R_CheckUserInterrupt` longjmps
  through Rust frames in the demo.
- [x] **Low: in the no-intercept case, a logical predictor is a factor in
  R's `model.matrix` (`xFALSE`, `xTRUE`) but numeric in the engine.**
  Check the name map.

## 4. lme4 functionality gaps (ranked by user impact)

1. [x] **High: stateful formula terms are refused** (`formula/parser.rs:381-446`).
   - `factor()`, `scale()`, `poly()`, `ns()`/`bs()`, `offset()`, `^`,
     `%in%`.
   - Parenthesised fixed groups such as `(a+b)^2` and `a*(b+c)`.
2. [x] **High: missing-data and argument handling.**
   - `lmm()` refuses NA by default.
   - `glmm()` cannot drop NA, and refuses `subset=` and `contrasts=`.
   - LMMs have no `offset`.
3. [x] **High: accessor shapes differ from lme4.**
   - `coef()` omits fixed effects that have no random term.
   - `VarCorr(m)$Subject` is NULL, and there are no `sc`/`stddev`
     attributes.
   - `AIC(m1, m2)` is refused.
   - `anova()` and `summary()` are not lme4/lmerTest-shaped.
   - lmerTest's `ddf=` is silently ignored by `summary()`.
4. [x] **High: GLMM methods are thin.**
   - `residuals()` is response-only, which breaks DHARMa and performance.
   - `simulate`, `refit`, `getME`, `is_singular` and `profile` are
     LMM-only.
   - `ranef(condVar = TRUE)` returns NA for GLMMs.
5. [x] **Medium: no diagnostic plots or influence measures.**
   - No `plot()` residual method.
   - No ranef `dotplot`/`qqmath`.
   - No `hatvalues`, `cooks.distance` or `influence`.
6. [x] **Medium: missing resampling tools.**
   - No `bootMer(FUN=)`.
   - No `rePCA`.
   - No GLMM parametric bootstrap through `simulate`.
7. [x] **Medium: families and links.**
   - Gamma's default inverse link is refused.
   - inverse.gaussian and binomial cauchit/log/identity are missing.
   - `nAGQ > 1` is profiled-only and needs a single scalar random effect,
     which is not stated in the docs.
   - `nAGQ = 0` is refused.
8. [x] **Medium: no `performance::icc`/`r2` or insight support.**
   - `model.frame()` lacks `terms` and transformed columns.
   - `weights()` returns NULL when the fit is unweighted.
9. [x] **Medium: `getME` is missing many components.**
   - Missing: u, b, sigma, Gp, L, RX, RZX, devcomp, lower, offset, weights,
     Ztlist, glmer.nb.theta, …
   - There is no GLMM method.
10. [x] **Medium: `confint` differs from lme4.**
  _Profile intervals are now on lme4's `.sig01`/`.sigma` scale with lme4 names (REML fits profiled on the ML deviance, as lme4); `"Wald"`/`"boot"` spellings accepted. Wald stays the default (documented) because profiling cost grows with model size._
    - It defaults to Wald where lme4 defaults to profile.
    - It rejects the `"Wald"` and `"boot"` spellings.
    - There is no GLMM profile.
11. [x] **Low: `predict` and `simulate` limits.**
    - Partial `re.form` formulas such as `~(1|g)` are refused.
    - GLMM `predict(newdata)` is refused when the fit has an offset.
12. [~] **Low: `nlmer`, `lmList`, `allFit`, `REMLcrit`, the modular
  _Done: `mm_lmlist()`, `mm_allfit()`, `REMLcrit()`, `mm_control(optCtrl =)`. Declined: `nlmer` (non-goal), `devFunOnly` (the engine owns the objective), `check.conv.*` (convergence is reported as a typed certificate)._
    `devFunOnly` API, `optCtrl` and `check.conv.*` are missing.**
13. [x] **Low: pbkrtest's KRmodcomp model-comparison form and lmerTest's
    `step()` are missing.**

## 5. Performance

1. [x] **Re-vendor the engine.** This is the largest single win; see §1.
2. [x] **`summary(fit)` refits once per coefficient.**
  _(a) done: one bridge call for all coefficients. (b) warm starts used for bootstrap refits only (same-model refits stay cold so summaries reproduce the fit exactly). (c) done in 4e317d3 (engine mixeff-rs#10, bridge-perf): `fit$rust_handle` keeps the fitted engine model; every follow-on entry point (contrasts, summary/anova term tables, bootstraps, predict and prediction variance, condVar, compare, boundary LRT, profile, verify_convergence) runs on the live model (a clone for the mutating ones), and falls back to the cold refit, with identical results, when the handle is dead after saveRDS()/readRDS() or does not match the requested refit. 200k-row crossed LMM: contrast 0.39 s → 0.002 s, predict(interval) 0.70 s → 0.04 s, summary 5.4 s → 4.8 s (the rest is the Satterthwaite computation itself)._
   - The default resolves to Satterthwaite and then runs
     `mm_inference_table_recompute` (`R/revive.R:998`), which makes one
     `contrast()` call per coefficient.
   - Each call is a cold bridge rebuild and refit.
   - Fixes:
     - (a) one call with `L = diag(p)`;
     - (b) pass `start = fit$theta` on every bridge refit;
     - (c) implement the `rust_handle` cache, which is always NULL today.
   - About 20 `lib.rs` entry points refit cold.
3. [x] **The dense Cholesky is unblocked and cache-hostile**
   (`model/linear/blocks.rs:2816`).
   - A crossed model with 2000 second-factor levels takes 7.8 s per
     evaluation; InstEval takes about 28 s.
   - A blocked version measured 236 ms → 34 ms at n=1000 and
     7.5 s → 0.19 s at n=2000, and agrees to about 1e-13.
   - The same stride problem affects `rdiv_lower_transpose`;
     `rank_k_downdate` should be a symmetric rank-k update.
4. [x] **The fit forces the deferred certificate** (`lib.rs:485`, `:629`).
  _Declined: mixeff fits are self-contained R objects (saved and revived
  without a live engine handle), so the bridge must complete the
  certificate's derivative evidence at fit time; deferring it would ship
  fits whose `optimizer_certificate()` lacks gradient/Hessian evidence.
  Revisit together with a `rust_handle` cache._ The cache landed in
  4e317d3 but fits must still survive serialization without it, so the
  certificate stays complete at fit time (as the engine report recommends).
   That adds 16% on InstEval and about 20% on GLMMs.
5. [x] **GLMM A-block rebuild per PIRLS iteration.**
   - `compute_wtxy_cross_product` and `recompute_*_a_blocks` are
     latency-bound and allocate on every iteration.
   - Estimated −25–35% on verbagg and contraception.
6. [x] **Kenward-Roger materializes dense n×n matrices**
   (`linear/mod.rs:2683-2990`). A Woodbury formulation through L would
   avoid them.
7. [x] **The Satterthwaite Jacobian is finite-difference** with 2(d+1)
   factorizations, and clones L each time.
8. [x] **Bootstrap and profile could run in parallel**: 2 threads under
   CRAN, more on user opt-in.
9. [x] **Bridge overhead.**
  _Per-observation vectors now cross as R doubles and are stored once (200k rows: 1.19 s to 0.61 s, 50 MB to 11.7 MB). The duplicate compile/audit pass and second data translation are gone in 4e317d3: lmm()/glmm() compile once through the engine's `CompiledModelSpec` (mixeff-rs#10) and fit from it (fit time unchanged within noise: 200k rows 0.44 s → 0.43 s, sleepstudy 9 ms → 8 ms)._
   - Data is translated, and the audit run, twice per `lmm()`.
   - Factors cross as character vectors.
   - JSON numeric arrays come back as lists.
   - Each observation is stored about three times in the fit object.
10. [x] **Small GLMM allocation items.** The L-block clones, the u clone, b
    reallocation, and the AGQ Xβ recompute.

## 6. CRAN cosmetics

- [x] `\dontrun` on fast examples (`audit`, `compile_model`,
  `explain_model`).
- [x] The `mm_json_negotiate` error example should use `try()`.
- [x] The `mm_lincomb` example is stale: it uses engine-style names and a
  nonexistent `fit$data`.
- [x] `src/Makevars.win.in` does not add `~/.cargo/bin` to `PATH`.

Already fine:

- the offline vendored build, `-j 2`, the cleanup, and toolchain-version
  reporting;
- no threads;
- `inst/AUTHORS` lists all 69 vendored crates;
- tests write only to `tempdir()`;
- `options()` changes are restored;
- Suggests packages are guarded.
