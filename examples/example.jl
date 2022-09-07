using LinearAlgebra 
using NonConvexPenalizedRegression 
using Random
using RCall
using Test

rng = MersenneTwister(1234);

n, p = 50, 5

X = randn(rng, n, p)        # feature matrix

a0 = collect(1:p)           # ground truths

y = X * a0 + 0.1 * randn(n) # generate response

XX = hcat(X, randn(rng, n, p))

@rput XX

@rput y

R"library(ncvreg)"

R"rscad <- coef(ncvreg(XX, y, lambda=0.2, penalty='SCAD', eps=.0001))"

@rget rscad

λ = [0.2]

jscad = NonConvexPenalizedRegression.coef(SCAD(XX, y, λ))

@show rscad .- jscad

@test rscad ≈ jscad atol=1e-6
