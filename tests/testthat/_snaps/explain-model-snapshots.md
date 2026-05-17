# snapshot: explain_model() for the eight §9.5.7 syntax patterns

    Code
      cat(explain_text(y ~ t + (1 | s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | s)
      
      Random effects:
        r0:
          wrote:      (1 | s)
          canonical:  (1 | s)
          named form: re(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
      
      Design notes:
        scope_note: `t` varies within `s`, so a `s`-level slope is structurally possible

---

    Code
      cat(explain_text(y ~ t + (0 + t | s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (0 + t | s)
      
      Random effects:
        r0:
          wrote:      (0 + t | s)
          canonical:  (0 + t | s)
          named form: re(group = s, intercept = FALSE, slopes = t, cov = "scalar")
          scope:      `s` units may differ in their `t` slope.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  t=present

---

    Code
      cat(explain_text(y ~ t + (1 + t | s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 + t | s)
      
      Random effects:
        r0:
          wrote:      (1 + t | s)
          canonical:  (1 + t | s)
          named form: re(group = s, intercept = TRUE, slopes = t, cov = "full")
          scope:      `s` units differ in baseline and `t` slope; the model estimates whether these are associated.
          covariance: full; theta parameters: 3
          support:    too_rich; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed; t=present

---

    Code
      cat(explain_text(y ~ t + (1 + t || s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 + t || s)
      
      Random effects:
        s has 2 separate random-effect blocks.
        r0:
          wrote:      (1 + t || s)
          canonical:  (1 | s)
          named form: re(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
        r1:
          wrote:      (1 + t || s)
          canonical:  (0 + t | s)
          named form: re(group = s, intercept = FALSE, slopes = t, cov = "scalar")
          scope:      `s` units may differ in their `t` slope.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  t=present
      
      Relationship between blocks:
        r0 <-> r1 (Intercept <-> t): double-bar syntax fixes the covariance between `Intercept` and `t` to zero.
      
      Design notes:
        covariance_assumption: the covariance between 'Intercept' and 't' is fixed at zero by || syntax

---

    Code
      cat(explain_text(y ~ t + (1 | s) + (0 + t | s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | s) + (0 + t | s)
      
      Random effects:
        s has 2 separate random-effect blocks.
        r0:
          wrote:      (1 | s)
          canonical:  (1 | s)
          named form: re(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
        r1:
          wrote:      (0 + t | s)
          canonical:  (0 + t | s)
          named form: re(group = s, intercept = FALSE, slopes = t, cov = "scalar")
          scope:      `s` units may differ in their `t` slope.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  t=present
      
      Relationship between blocks:
        r0 <-> r1 (Intercept <-> t): separate random-effect blocks fix the covariance between `Intercept` and `t` to zero.
      
      Design notes:
        covariance_assumption: the covariance between 'Intercept' and 't' is fixed at zero by separate random-effect blocks

---

    Code
      cat(explain_text(y ~ t + (1 | a / b)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | a) + (1 | a:b)
      
      Random effects:
        r0:
          wrote:      (1 | a/b)
          canonical:  (1 | a)
          named form: re(group = a, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `a` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    weakly_supported; group levels: 3; min rows/group: 16; median rows/group: 16
          variation:  intercept=not_assessed
        r1:
          wrote:      (1 | a/b)
          canonical:  (1 | a:b)
          named form: re(group = a:b, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `a:b` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
      
      Design notes:
        syntax_expansion: random-effect shorthand expands to (1 | a) + (1 | a:b)
        support_note: the requested covariance structure is information-hungry relative to the observed grouping levels
        scope_note: `t` varies within `a`, so a `a`-level slope is structurally possible
        scope_note: `t` varies within `a:b`, so a `a:b`-level slope is structurally possible

---

    Code
      cat(explain_text(y ~ t + (1 | s:i)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | s:i)
      
      Random effects:
        r0:
          wrote:      (1 | s:i)
          canonical:  (1 | s:i)
          named form: re(group = s:i, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s:i` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 24; min rows/group: 2; median rows/group: 2
          variation:  intercept=not_assessed

---

    Code
      cat(explain_text(y ~ t + (1 | s) + (1 | i)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | s) + (1 | i)
      
      Random effects:
        r0:
          wrote:      (1 | s)
          canonical:  (1 | s)
          named form: re(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
        r1:
          wrote:      (1 | i)
          canonical:  (1 | i)
          named form: re(group = i, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `i` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    weakly_supported; group levels: 4; min rows/group: 12; median rows/group: 12
          variation:  intercept=not_assessed
      
      Design notes:
        support_note: the requested covariance structure is information-hungry relative to the observed grouping levels
        scope_note: `t` varies within `s`, so a `s`-level slope is structurally possible
        scope_note: `t` varies within `i`, so a `i`-level slope is structurally possible
        covariance_assumption: no correlation parameter is estimated between random-intercept groups 's' and 'i'; separate scalar random-effect terms define independent covariance blocks

# snapshot: three kinds of help register cleanly

    Code
      cat(explain_text(y ~ t + (1 | s)))
    Output
      Random effects explanation:
        formula: y ~ 1 + t + (1 | s)
      
      Random effects:
        r0:
          wrote:      (1 | s)
          canonical:  (1 | s)
          named form: re(group = s, intercept = TRUE, slopes = NULL, cov = "scalar")
          scope:      `s` units may differ in average outcome.
          covariance: scalar; theta parameters: 1
          support:    sufficient; group levels: 6; min rows/group: 8; median rows/group: 8
          variation:  intercept=not_assessed
      
      Design notes:
        scope_note: `t` varies within `s`, so a `s`-level slope is structurally possible

---

    Code
      cat(explain_text(y ~ between + (1 + between | g), mk_refusal_design()))
    Output
      Random effects explanation:
        formula: y ~ 1 + between + (1 + between | g)
      
      Random effects:
        r0:
          wrote:      (1 + between | g)
          canonical:  (1 + between | g)
          named form: re(group = g, intercept = TRUE, slopes = between, cov = "full")
          scope:      `g` units differ in baseline and `between` slope; the model estimates whether these are associated.
          covariance: full; theta parameters: 3
          support:    too_rich; group levels: 4; min rows/group: 2; median rows/group: 2
          variation:  between=absent; intercept=not_assessed
      
      Possible repairs, not applied automatically:
        1. structural_refusal: `between` does not vary within `g`, so a `g`-level `between` slope cannot be estimated from this design.

# snapshot: structural_refusal renders 'Possible repairs' wording

    Code
      cat(refusal_audit$text)
    Output
      Audit Summary:
        overall [WARNING]: 8 warning(s), 4 not checked item(s); review attention lines before treating inference as routine
        attention [WARNING]: Model State / supported: status=refused; formula=y ~ 1 + between + (1 + between | g); random_terms=1; reason=design audit found at least one unsupported random-effect distribution
        attention [WARNING]: Model State / changes: Recommended:DesignTime:r0 -> basis direction(s) unsupported by within-group variation: between | Recommended:DesignTime:r0 -> number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable
        attention [WARNING]: Random Effects / r0: group=g, rows=8, levels=4, obs_per_level=2..2, basis=2, covariance=full, params=3, budget=too_rich; reason=number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable
        attention [WARNING]: Random-Effect Information Budget / r0: levels=4, rows=8, obs_per_level=2..2, basis=2, cov_params=3, levels/basis=2.00, levels/param=1.33, rows/param=2.67; total rows can be misleading for covariance support; risk=maximal covariance structure is too rich for the grouping-level budget; recommendation=random-effect coefficients saturate the rows for this term; drop unsupported random slopes, split/simplify the random-effect structure, treat the grouping factor as fixed when appropriate, or collect more observations per grouping level; explanation=8 rows are clustered into 4 grouping levels; covariance support is limited by grouping levels, not by total rows
        attention [WARNING]: Random Term Cards / r0: original=(1 + between | g), canonical=(1 + between | g), group=g, blocks=basis=[intercept, between], covariance=full, params=3
        attention [WARNING]: Policy Recommendations / r0: drop_unsupported_basis: basis direction(s) unsupported by within-group variation: between; recommended covariance=scalar_or_diagonal_on_supported_basis; inference=fixed-effect inference would be conditional on the supported random-effect basis
        attention [WARNING]: Policy Recommendations / r0: refuse_random_term_distribution: number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable; inference=confirmatory fixed-effect p-values should be withheld or recomputed after a declared design-level change
        attention [WARNING]: Diagnostics / covariance_too_rich: number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable; affected=(1 + between | g)
        attention [NOT CHECKED]: Effective Covariance / effective covariance rank: not assessed
        attention [NOT CHECKED]: Optimizer / certificate: model has not been fitted
        attention [NOT CHECKED]: Inference / finite-sample inference: finite-sample inference is not implemented in compiler v0
        attention [NOT CHECKED]: Inference / covariance derivatives: compiler v0 does not expose covariance derivative certificates
      
      Requested Model:
        formula [INFO]: y ~ 1 + between + (1 + between | g)
        model kind [INFO]: linear_mixed_model
        distribution/link [INFO]: gaussian/identity
        objective [INFO]: exact_gaussian
        certificate scope [INFO]: exact_objective
        fixed terms [INFO]: 1, between
        random terms [INFO]: 1
        theta maps [INFO]: 1 map(s)
      
      Model State:
        requested [OK]: status=requested; formula=y ~ 1 + between + (1 + between | g); random_terms=1; reason=formula as requested by the caller
        semantic [OK]: status=canonical; formula=y ~ 1 + between + (1 + between | g); random_terms=1; reason=formula compiled into semantic IR
        supported [WARNING]: status=refused; formula=y ~ 1 + between + (1 + between | g); random_terms=1; reason=design audit found at least one unsupported random-effect distribution
        fitted [NOT CHECKED]: status=not_assessed; formula=y ~ 1 + between + (1 + between | g); random_terms=1; reason=model has not been fitted
        changes [WARNING]: Recommended:DesignTime:r0 -> basis direction(s) unsupported by within-group variation: between | Recommended:DesignTime:r0 -> number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable
      
      Fixed Effects:
        rank [OK]: 2 of 2
        aliased columns [OK]: none
        empty cells [OK]: 0
      
      Random Effects:
        r0 [WARNING]: group=g, rows=8, levels=4, obs_per_level=2..2, basis=2, covariance=full, params=3, budget=too_rich; reason=number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable
      
      Random-Effect Information Budget:
        r0 [WARNING]: levels=4, rows=8, obs_per_level=2..2, basis=2, cov_params=3, levels/basis=2.00, levels/param=1.33, rows/param=2.67; total rows can be misleading for covariance support; risk=maximal covariance structure is too rich for the grouping-level budget; recommendation=random-effect coefficients saturate the rows for this term; drop unsupported random slopes, split/simplify the random-effect structure, treat the grouping factor as fixed when appropriate, or collect more observations per grouping level; explanation=8 rows are clustered into 4 grouping levels; covariance support is limited by grouping levels, not by total rows
      
      Random Term Cards:
        r0 [WARNING]: original=(1 + between | g), canonical=(1 + between | g), group=g, blocks=basis=[intercept, between], covariance=full, params=3
      
      Cross-Card Constraints:
        constraints [OK]: none
      
      Dependence Paths:
        kernels [INFO]: r0=marginal(g, intercept=true, covariance=full, basis=intercept, between)
        repeated units [OK]: g=marginal(g, levels=4, obs_per_level=2..2, covered_by=r0)
        missing paths [OK]: none
      
      Parameterization Trace:
        r0 [INFO]: source=(1 + between | g); group=g; family=FullCholesky; user_basis=intercept, between; optimizer_basis=intercept, between; theta_slots=theta[0], theta[1], theta[2]; lambda_slots=(0, 0), (1, 0), (1, 1); parmap_aligned=0/3; varcorr_entries=sd(intercept), sd(between), corr(between,intercept)
      
      Effective Covariance:
        effective covariance rank [NOT CHECKED]: not assessed
      
      Policy Recommendations:
        r0 [WARNING]: drop_unsupported_basis: basis direction(s) unsupported by within-group variation: between; recommended covariance=scalar_or_diagonal_on_supported_basis; inference=fixed-effect inference would be conditional on the supported random-effect basis
        r0 [WARNING]: refuse_random_term_distribution: number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable; inference=confirmatory fixed-effect p-values should be withheld or recomputed after a declared design-level change
      
      Optimizer:
        certificate [NOT CHECKED]: model has not been fitted
      
      Inference:
        finite-sample inference [NOT CHECKED]: finite-sample inference is not implemented in compiler v0
        covariance derivatives [NOT CHECKED]: compiler v0 does not expose covariance derivative certificates
      
      Diagnostics:
        structural_refusal [INFO]: `between` does not vary within `g`, so a `g`-level `between` slope cannot be estimated from this design; affected=(1 + between | g); suggested=`between` does not vary within `g`, so a `g`-level `between` slope cannot be estimated from this design.
        covariance_too_rich [WARNING]: number of observations (8) is <= random coefficients (8) for grouping factor 'g' with basis dimension 2; random-effect variances and the residual scale are probably not separately identifiable; affected=(1 + between | g)

# snapshot: random_options() prints rung 0 first-class, no ranking

    Code
      print(random_options(spec, "s", "t"))
    Output
      Random-effect options for group: s
      Current model:
        (1 | s) <- this is what you wrote
        `s` units may differ in average outcome.
      Nearby options:
        (1 | s) <- this is what you wrote
          varying coefficients: intercept
          covariance family:    scalar
          theta parameters:     1
          design status:        sufficient
          plain meaning:        `s` units may differ in average outcome.
        (0 + t | s)
          varying coefficients: t
          covariance family:    scalar
          theta parameters:     1
          design status:        sufficient
          plain meaning:        `s` units may differ in their `t` slope.
        (1 | s) + (0 + t | s)
          varying coefficients: intercept, t
          covariance family:    diagonal via separate blocks
          theta parameters:     2
          design status:        sufficient
          plain meaning:        `s` units may differ in average outcome. `s` units may differ in their `t` slope. separate random-effect blocks fix the covariance between `Intercept` and `t` to zero.
        (1 + t || s)
          varying coefficients: intercept, t
          covariance family:    diagonal via separate blocks
          theta parameters:     2
          design status:        sufficient
          plain meaning:        `s` units may differ in average outcome. `s` units may differ in their `t` slope. double-bar syntax fixes the covariance between `Intercept` and `t` to zero.
        (1 + t | s)
          varying coefficients: intercept, t
          covariance family:    full
          theta parameters:     3
          design status:        too_rich
          plain meaning:        `s` units differ in baseline and `t` slope; the model estimates whether these are associated.

# snapshot: compare_covariance() layout

    Code
      print(compare_covariance(spec))
    Output
      Covariance comparison:
        r0 / s / full <- current
          basis:            intercept, t
          theta parameters: 3
          assumes zero:     none
          design status:    too_rich
        r0 / s / diagonal
          basis:            intercept, t
          theta parameters: 2
          assumes zero:     off-diagonal covariances
          design status:    too_rich
        r0 / s / scalar
          basis:            intercept, t
          theta parameters: 1
          assumes zero:     off-diagonal covariances
          design status:    too_rich

# snapshot: pedagogical DiagnosticCode variants round-trip from Rust

    Code
      print(scope)
    Output
      Diagnostics:
             code severity        stage affected_terms
       scope_note     info design_audit             r0
      
      Messages:
        scope_note: `t` varies within `s`, so a `s`-level slope is structurally possible

---

    Code
      print(cov)
    Output
      Diagnostics:
                        code severity       stage affected_terms
       formula_canonicalized     info semantic_ir   (1 + t || s)
       covariance_assumption     info semantic_ir   (1 + t || s)
      
      Messages:
        formula_canonicalized: random-effect term was canonicalized as (1 | s) + (0 + t | s)
        covariance_assumption: the covariance between 'Intercept' and 't' is fixed at zero by || syntax

---

    Code
      print(expand)
    Output
      Diagnostics:
                           code severity        stage affected_terms
          formula_canonicalized     info  semantic_ir      (1 | a/b)
               syntax_expansion     info  semantic_ir      (1 | a/b)
       random_effect_few_levels  warning design_audit        (1 | a)
                   support_note     info design_audit             r0
                     scope_note     info design_audit             r0
                     scope_note     info design_audit             r1
        repeated_unit_unmodeled  warning design_audit              b
      
      Messages:
        formula_canonicalized: random-effect term was canonicalized as (1 | a) + (1 | a:b)
        syntax_expansion: random-effect shorthand expands to (1 | a) + (1 | a:b)
        random_effect_few_levels: 3 levels are fit-eligible for a scalar random intercept but below the v0
              reliability threshold 5
        support_note: the requested covariance structure is information-hungry relative to the
              observed grouping levels
        scope_note: `t` varies within `a`, so a `a`-level slope is structurally possible
        scope_note: `t` varies within `a:b`, so a `a:b`-level slope is structurally possible
        repeated_unit_unmodeled: repeated marginal unit 'b' is not covered by a random-intercept dependence
              path

