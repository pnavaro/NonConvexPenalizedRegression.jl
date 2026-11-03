
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

abstract type AbstractModel end

export SCAD

struct SCAD <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    function SCAD(X, y, λ; gamma=3.7, kwargs...)
        fit = ncvreg(X, y; lambda=λ, penalty=:SCAD, gamma, kwargs...)
        return new(fit.beta, fit.lambda)
    end

end
