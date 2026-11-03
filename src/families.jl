abstract type AbstractFamily end

struct GaussianFamily <: AbstractFamily end
struct BinomialFamily <: AbstractFamily end
struct PoissonFamily <: AbstractFamily end

"""
    RegressionPath

Regularization path with an intercept in the first row of `beta`. Only attempted
lambda values are retained. `iter` counts coordinate sweeps, `status` records
each fit's outcome, and `stop_reason` explains an early path termination.
`setup_iter` counts sweeps used to initialize unpenalized predictors.
"""
struct RegressionPath{F<:AbstractFamily}
    beta::Matrix{Float64}
    lambda::Vector{Float64}
    family::F
    penalty::Symbol
    deviance::Vector{Float64}
    converged::BitVector
    iter::Vector{Int}
    status::Vector{Symbol}
    stop_reason::Symbol
    setup_iter::Int
    gamma::Float64
    alpha::Float64
    penalty_factor::Vector{Float64}
    eps::Float64
    max_iter::Int
    dfmax::Int
    nobs::Int
    locally_convex::Union{Nothing,BitVector}
    convex_min::Union{Nothing,Int}
    X::Union{Nothing,Matrix{Float64}}
    y::Union{Nothing,Vector{Float64}}
end

"""
    coef(path::RegressionPath; index=nothing)

Return the coefficient matrix for the full path, or a coefficient vector for
one path index.
"""
function coef(path::RegressionPath; index=nothing)
    isnothing(index) && return path.beta
    return path.beta[:, index]
end

_response(::GaussianFamily, eta) = eta
_response(::BinomialFamily, eta) = _logistic.(eta)
_response(::PoissonFamily, eta) = exp.(clamp.(eta, -700, 700))

export AbstractFamily, GaussianFamily, BinomialFamily, PoissonFamily
export RegressionPath, coef, predict
