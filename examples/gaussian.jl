# -*- coding: utf-8 -*-
# +
using LinearAlgebra 
using NonConvexPenalizedRegression 
using Random
using RCall
using Test
using RData
    
data = load(joinpath(@__DIR__, "..", "test", "data", "Prostate.rda"))
X  = data["Prostate"]["X"]
y  = data["Prostate"]["y"]
@rput X
@rput y

# +
R"""
library(ncvreg)
fit <- ncvreg(X, y)
rscad <- coef(fit, lambda=0.05)
"""

@rget rscad
# -

λ = [0.05]
fit = MCP(X, y, λ)
jscad = NonConvexPenalizedRegression.coef(fit)
