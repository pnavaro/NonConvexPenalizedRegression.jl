# NonConvexPenalizedRegression.jl

Regularization Paths for SCAD and MCP Penalized Regression Models

This is a partial Julia translation of the R package [ncvreg](http://pbreheny.github.io/ncvreg/).
The keyword-based `ncvreg` interface supports Gaussian, binomial, and Poisson
families, with SCAD, MCP, or Lasso penalties. Binomial responses use 0/1
coding; Poisson responses must be nonnegative integer counts.

Start with the [Getting started](@ref) tutorial for a complete analysis of the
Prostate dataset, including coefficient paths, predictions, and diagnostics.

```julia
fit = ncvreg(X, y;
    family=BinomialFamily(),
    penalty=:SCAD,
    lambda=nothing, # generate a decreasing path automatically
)

probabilities = predict(fit, X)
coefficients = coef(fit)
```

For Poisson regression, set `family=PoissonFamily()`; predictions are expected
counts. The existing positional `SCAD(X, y, lambda)`, `MCP(X, y, lambda)`,
and `Lasso(X, y, lambda)` constructors remain available for Gaussian fits.
The model also provides `fitted`, `predict`, `residuals`, `deviance`,
`loglikelihood`, and `summary` methods. RecipesBase recipes enable coefficient
paths, deviance paths, and residual-versus-fitted plots when used with Plots.jl.

## Lambda-path setup and fitting options

```julia
fit = ncvreg(X, y;
    family=GaussianFamily(),
    penalty=:MCP,
    alpha=0.8,
    penalty_factor=ones(size(X, 2)),
    nlambda=100,
    lambda_min_ratio=0,
    dfmax=size(X, 2),
    max_iter=10000,
    return_X=true,
)

summary(fit)
fitted(fit)
residuals(fit; type=:deviance)
```

| Keyword | Behavior |
|---|---|
| `lambda` | A finite, nonnegative scalar or nonincreasing vector. Omit it to generate a log-spaced path from the largest standardized null-model score. |
| `nlambda` | Number of automatically generated values; default 100. One value requests only the maximum. |
| `lambda_min_ratio` | Minimum lambda divided by the maximum. Defaults to 0.001 when `n > p`, and 0.05 otherwise. Zero appends an unpenalized endpoint, using 0.001 for the preceding positive endpoint; at least two values are required. |
| `gamma` | Penalty shape: greater than 2 for SCAD (default 3.7), greater than 1 for MCP (default 3). |
| `alpha` | Mixes SCAD/MCP/Lasso with ridge; default 1. Pure ridge (`alpha=0`) requires explicit lambdas because it has no finite all-zero threshold. |
| `penalty_factor` | Nonnegative multiplier per predictor, applied to both penalty components without normalization. Zero means no shrinkage, not exclusion. |
| `eps` | Positive convergence tolerance for coefficient changes and first-order optimality; default `1e-6`. |
| `max_iter` | Total coordinate sweeps across the whole path, including unpenalized initialization; default 10000. |
| `dfmax` | Stop after retaining the first model exceeding this many nonzero slopes, including unpenalized slopes. The intercept is excluded. Defaults to the number of predictors. |
| `warn` | Whether to warn on early termination; default `true`. Status fields remain available when warnings are suppressed. |
| `convex` | Compute positive definiteness of the active-set penalized Hessian at each fitted lambda; default `true`. |
| `return_X` | Retain copies of the original design and response for methods without explicit data; default `false`. |

Predictors are always standardized, coefficients are returned on the original
scale, and the intercept is never penalized. SCAD/MCP use penalties on the
standardized coefficients, without R ncvreg's GLM adaptive rescaling.
Before generating a path, unpenalized predictors are fitted to establish the
correct null-model scores. Initialization failure raises an error rather than
returning a misleading path. If there are no penalized predictors or all scores
are zero, automatic setup returns an all-zero lambda sequence.

Each fit stores the actual `gamma`, `alpha`, `penalty_factor`, `eps`, `max_iter`,
and `dfmax`. Inspect `iter`, `setup_iter`, `converged`, `status`, and
`stop_reason` to distinguish completed paths from iteration-budget exhaustion,
`dfmax` stopping, numerical stalling, and GLM saturation (deviance below 1% of
the null deviance). Paths are truncated to attempted lambdas; the terminal fit
is retained even if it did not converge. Warnings can be suppressed with
`warn=false`, but errors and status information are never suppressed.

`locally_convex` records active-set positive definiteness, and `convex_min` is
the last index in its initial contiguous convex region (zero if the first
model fails the check). Both are `nothing` when `convex=false`. This is a
local diagnostic, not a guarantee of a unique global optimum.
Existing `SCAD(X, y, lambda; ...)`, `MCP(X, y, lambda; ...)`, and
`Lasso(X, y, lambda; ...)` constructors use the same fitting implementation.

Algorithm is described in **Breheny P and Huang J (2011)** "Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection". *Annals of Applied Statistics*, 5: 232–253

# Index

```@autodocs
Modules = [NonConvexPenalizedRegression]
```

```@index
```
