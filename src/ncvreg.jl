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
