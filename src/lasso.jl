function lasso(z::Float64, l1::Float64, l2::Float64)::Float64

    s::Float64 = sign(z)
    if abs(z) <= l1
        return 0.0
    else
        return s * (abs(z) - l1) / (1 + l2)
    end

end

export Lasso
struct Lasso <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    Lasso(X, y, λ) = new(ncvreg(X, y, λ, :Lasso, 3.0), λ)

end
