# NonConvexPenalizedRegression.jl

Regularization Paths for SCAD and MCP Penalized Regression Models

This is a translation in Julia of the R package [ncvreg](http://pbreheny.github.io/ncvreg/).  
Only `gaussian` family is translated.

Algorithm is described in **Breheny P and Huang J (2011)** "Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection". *Annals of Applied Statistics*, 5: 232–253

```@example scad
using LinearAlgebra 
using Random
using NonConvexPenalizedRegression
using RCall
rng = MersenneTwister(1234);
n, p = 50, 5

X = randn(rng, n, p)              # feature matrix
```

```@setup scad
R"install.packages('ncvreg')"
```

```@example scad
a0 = collect(1:p)                # ground truths
```

```@example scad
y = X * a0 + 0.1 * randn(n) # generate response
```

```@example scad
XX = hcat(X, randn(rng, n, p))

@rput XX
@rput y

R"library(ncvreg)"

R"scad <- coef(ncvreg(XX, y, lambda=0.2, penalty='SCAD', eps=.0001))"

@rget scad
```

```@example scad
λ = [0.2]

scad = NonConvexPenalizedRegression.coef(SCAD(XX, y, λ))
```

# Index

```@autodocs
Modules = [NonConvexPenalizedRegression]
```

```@index
```
