module NonConvexPenalizedRegression

using Statistics

ncol(A::AbstractMatrix) = size(A)[2]

function standardize(X)

    # Declarations
    n, p = size(X)

    XX = similar(X)
    c = zeros(Float64, p)
    s = zeros(Float64, p)

    for j = 1:p

        # Center
        c[j] = 0
        for i = 1:n
            c[j] += X[i, j]
        end
        c[j] = c[j] / n
        for i = 1:n
            XX[i, j] = X[i, j] - c[j]
        end

        # Scale
        s[j] = 0
        for i = 1:n
            s[j] += XX[i, j]^2
        end
        s[j] = sqrt(s[j] / n)
        for i = 1:n
            XX[i, j] = XX[i, j] / s[j]
        end

    end

    any(s .< 1e-6) && throw("Sigular matrix")

    XX, c, s

end


"""
    gaussian_loss( r )

Gaussian loss 
"""
function gaussian_loss(r::Vector{Float64})::Float64
    l::Float64 = 0
    for i in eachindex(r)
        l = l + r[i]^2
    end
    l
end

function scad(z::Float64, l1::Float64, l2::Float64, γ::Float64)::Float64

    s::Float64 = sign(z)

    if abs(z) <= l1
        return 0
    elseif abs(z) <= (l1 * (1 + l2) + l1)
        return s * (abs(z) - l1) / (γ * (1 + l2))
    elseif abs(z) <= γ * l1 * (1 + l2)
        return s * (abs(z) - γ * l1 / (γ - 1)) / (1 - 1 / (γ - 1) + l2)
    else
        return z / (1 + l2)
    end

end

function mcp(z::Float64, l1::Float64, l2::Float64, γ::Float64)::Float64

    s::Float64 = sign(z)

    if abs(z) <= l1
        return 0
    elseif abs(z) <= γ * l1 * (1 + l2)
        return s * (abs(z) - l1) / (1 + l2 - 1 / γ)
    else
        return z / (1 + l2)
    end

end

function lasso(z::Float64, l1::Float64, l2::Float64)::Float64

    s::Float64 = sign(z)
    if abs(z) <= l1
        return 0
    else
        return s * (abs(z) - l1) / (1 + l2)
    end

end

"""
   crossprod( X, y, n, j)

Cross product of y with jth column of X
"""
function crossprod(X, y, j)

    val::Float64 = 0.0

    for i in eachindex(y)
        val += X[i, j] * y[i]
    end

    return val

end

include("gaussian.jl")
include("binomial.jl")
include("poisson.jl")


abstract type AbstractModel end

export SCAD
struct SCAD <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    SCAD(X, y, λ) = new(ncvreg(X, y, λ, :SCAD, 3.7), λ)

end

export MCP
struct MCP <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    MCP(X, y, λ) = new(ncvreg(X, y, λ, :MCP, 3.0), λ)

end

export Lasso
struct Lasso <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    Lasso(X, y, λ) = new(ncvreg(X, y, λ, :Lasso, 3.0), λ)

end

function ncvreg(X, y, λ, penalty, γ)

    α = 1
    eps = 1e-4
    max_iter = 10000
    n, p = size(X)
    penalty_factor = ones(Float64, p)

    ## Set up XX, yy, λ
    XX, center, scale = standardize(X)
    dfmax = p + 1

    yy = y .- mean(y)

    @assert n == length(yy)
    @assert penalty in [:SCAD, :MCP, :Lasso]

    λ_min = ifelse(n > p, 0.001, 0.05)

    nλ = length(λ)

    β, loss, iter =
        cdfit_gaussian(XX, yy, penalty, λ, eps, max_iter, γ, penalty_factor, α, dfmax)

    ## Unstandardize
    b = β ./ scale
    a = [mean(y) for i = 1:nλ] .- vec(center' * b)

    beta = zeros(Float64, (ncol(X) + 1, p))
    beta .= transpose(collect(hcat(a, b')))

    beta

end

function coef(self::AbstractModel)

    which = 1:length(self.λ)

    return self.beta[:, which]

end

end # module
