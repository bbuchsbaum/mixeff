# snapshot: print(VarCorr(fit)) tags boundary components with [boundary]

    Code
      print(VarCorr(fit))
    Output
      Variance components:
       group        name variance std_dev correlation       note
       Batch (Intercept)        0       0             [boundary]
      [boundary]: variance component is at the boundary of the parameter space.
      Residual std. dev.: 3.71568

# snapshot: print(fit) on singular fit names rank and points to audit verbs

    Code
      cat(printed)
    Output
      Linear mixed model fit by REML
      Formula: Yield ~ 1 + (1 | Batch)
      Fit status: converged_reduced_rank
      Optimizer: pattern_search; iterations: 23; objective: 161.828
      Artifact: mixedmodels.compiled_model_artifact v1; crate: <version>
      nobs: 30, sigma: 3.71568, logLik: -80.9141
      Fixed effects:
      (Intercept) 
           5.6656 
      
      Fitted covariance state:
      The fitted covariance matrix is rank-deficient.
        r0: requested rank 1; fitted effective rank 0.
      Use changes(fit) to see which dimension was unsupported.
      Audit verbs: audit(), diagnostics(), inference_table(), model_report()

