
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

export MCP
struct MCP <: AbstractModel

    beta::Array{Float64,2}
    λ::Vector{Float64}

    MCP(X, y, λ) = new(ncvreg(X, y, λ, :MCP, 3.0), λ)

end

