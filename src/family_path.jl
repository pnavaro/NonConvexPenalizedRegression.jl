const _VALID_PENALTIES = (:SCAD, :MCP, :Lasso)

function _penalty_symbol(penalty)
    penalty isa Union{Symbol,AbstractString} ||
        throw(ArgumentError("penalty must be :SCAD, :MCP, or :Lasso"))
    name = Symbol(penalty)
    name == :lasso && (name = :Lasso)
    name in _VALID_PENALTIES ||
        throw(ArgumentError("penalty must be :SCAD, :MCP, or :Lasso"))
    return name
end

function _linear_predictor(X, beta, intercept)
    return intercept .+ X * beta
end

function _logistic(eta::Real)
    if eta >= 0
        return 1 / (1 + exp(-eta))
    end
    e = exp(eta)
    return e / (1 + e)
end

function _family_mean(::GaussianFamily, eta)
    return eta
end

function _family_mean(::BinomialFamily, eta)
    return _logistic.(eta)
end

function _family_mean(::PoissonFamily, eta)
    return exp.(eta)
end

function _family_weights(::GaussianFamily, mu)
    return ones(length(mu))
end

function _family_weights(::BinomialFamily, mu)
    return max.(mu .* (1 .- mu), 1e-8)
end

function _family_weights(::PoissonFamily, mu)
    return max.(mu, 1e-8)
end

function _family_deviance(::GaussianFamily, y, mu)
    return sum(abs2, y .- mu)
end

function _family_deviance(::BinomialFamily, y, mu)
    dev = 0.0
    floor = eps(Float64)
    for i in eachindex(y)
        probability = clamp(mu[i], floor, 1 - floor)
        if y[i] > 0
            dev += 2 * y[i] * log(y[i] / probability)
        end
        if y[i] < 1
            dev += 2 * (1 - y[i]) * log((1 - y[i]) / (1 - probability))
        end
    end
    return dev
end

function _family_deviance(::PoissonFamily, y, mu)
    dev = 0.0
    for i in eachindex(y)
        mean_i = max(mu[i], floatmin(Float64))
        if y[i] > 0
            dev += 2 * (y[i] * log(y[i] / mean_i) - (y[i] - mean_i))
        else
            dev += 2 * mean_i
        end
    end
    return dev
end

function _penalty_value(beta, lambda1, lambda2, penalty, gamma)
    value = abs(beta)
    if penalty == :Lasso
        nonconvex_penalty = lambda1 * value
    elseif penalty == :MCP
        nonconvex_penalty = value <= gamma * lambda1 ?
            lambda1 * value - value^2 / (2 * gamma) : gamma * lambda1^2 / 2
    elseif value <= lambda1
        nonconvex_penalty = lambda1 * value
    elseif value <= gamma * lambda1
        nonconvex_penalty =
            (-value^2 + 2 * gamma * lambda1 * value - lambda1^2) /
            (2 * (gamma - 1))
    else
        nonconvex_penalty = (gamma + 1) * lambda1^2 / 2
    end
    return nonconvex_penalty + lambda2 * beta^2 / 2
end

function _family_loss(::GaussianFamily, y, eta)
    return sum(abs2, y .- eta) / (2 * length(y))
end

function _family_loss(::BinomialFamily, y, eta)
    return sum(max.(eta, 0) .- y .* eta .+ log1p.(exp.(-abs.(eta)))) /
           length(y)
end

function _family_loss(::PoissonFamily, y, eta)
    return sum(exp.(eta) .- y .* eta) / length(y)
end

function _objective(family, y, eta, beta, factors, lambda, alpha, penalty, gamma)
    value = _family_loss(family, y, eta)
    for j in eachindex(beta)
        value += _penalty_value(
            beta[j],
            lambda * factors[j] * alpha,
            lambda * factors[j] * (1 - alpha),
            penalty,
            gamma,
        )
    end
    return value
end

function _penalty_coordinate(z, h, lambda1, lambda2, penalty, gamma)
    ridge = h + lambda2
    if lambda1 == 0
        return z / ridge
    end

    value = abs(z)
    sign_z = sign(z)

    if penalty == :Lasso
        return sign_z * max(value - lambda1, 0) / ridge
    end

    # Evaluate all piecewise minima, including boundaries when the quadratic
    # approximation has less curvature than the concave part of the penalty.
    boundary = gamma * lambda1
    candidates = [0.0, boundary, max(value / ridge, boundary)]
    if penalty == :MCP
        curvature = ridge - 1 / gamma
        if curvature > 0
            push!(candidates, clamp((value - lambda1) / curvature, 0, boundary))
        end
    else
        push!(candidates, lambda1, clamp((value - lambda1) / ridge, 0, lambda1))
        curvature = ridge - 1 / (gamma - 1)
        if curvature > 0
            push!(candidates, clamp(
                (value - gamma * lambda1 / (gamma - 1)) / curvature,
                lambda1, boundary,
            ))
        end
    end
    best = 0.0
    minimum_value = 0.0
    for candidate in candidates
        objective = h * candidate^2 / 2 - value * candidate +
                    _penalty_value(candidate, lambda1, lambda2, penalty, gamma)
        if objective < minimum_value
            best, minimum_value = candidate, objective
        end
    end
    return sign_z * best
end

function _standardize_predictors(X::AbstractMatrix{<:Real})
    n, p = size(X)
    n > 0 && p > 0 || throw(ArgumentError("X must have at least one row and column"))
    all(isfinite, X) || throw(ArgumentError("X must contain only finite values"))

    means = vec(mean(X; dims=1))
    scales = vec(sqrt.(mean((X .- means') .^ 2; dims=1)))
    constant = scales .< 1e-12
    scales[constant] .= 1.0
    Xs = (Matrix{Float64}(X) .- means') ./ scales'
    all(isfinite, Xs) && all(isfinite, means) && all(isfinite, scales) ||
        throw(ArgumentError("X cannot be standardized to finite Float64 values"))
    Xs[:, constant] .= 0.0
    return Xs, means, scales
end

function _validate_response(::GaussianFamily, y)
    all(isfinite, y) || throw(ArgumentError("y must contain only finite values"))
end

function _validate_response(::BinomialFamily, y)
    all(v -> isfinite(v) && (v == 0 || v == 1), y) ||
        throw(ArgumentError("binomial responses must be 0/1 values"))
    any(y .== 0) && any(y .== 1) ||
        throw(ArgumentError("binomial response must contain both classes"))
end

function _validate_response(::PoissonFamily, y)
    all(v -> isfinite(v) && v >= 0 && isinteger(v), y) ||
        throw(ArgumentError("Poisson responses must be nonnegative integer counts"))
    any(y .> 0) || throw(ArgumentError("Poisson response must include a positive count"))
end

function _penalty_derivative(value, lambda1, penalty, gamma)
    lambda1 == 0 && return 0.0
    penalty == :Lasso && return lambda1
    penalty == :MCP && return max(lambda1 - value / gamma, 0)
    value <= lambda1 && return lambda1
    return max((gamma * lambda1 - value) / (gamma - 1), 0)
end

function _stationarity(X, y, family, eta, beta, factors, lam, alpha, penalty, gamma, columns)
    residual = y .- _family_mean(family, eta)
    violation = abs(mean(residual))
    for j in columns
        gradient = dot(view(X, :, j), residual) / length(y)
        lambda1 = lam * factors[j] * alpha
        lambda2 = lam * factors[j] * (1 - alpha)
        if beta[j] == 0
            violation = max(violation, abs(gradient) - lambda1)
        else
            derivative = _penalty_derivative(abs(beta[j]), lambda1, penalty, gamma)
            violation = max(violation, abs(
                gradient - sign(beta[j]) * derivative - lambda2 * beta[j],
            ))
        end
    end
    return violation
end

function _fit_lambda!(
    X, y, family, beta, intercept, factors, lam, alpha, penalty, gamma,
    tolerance, budget, columns; null_deviance=0.0,
)
    n = length(y)
    eta = _linear_predictor(X, beta, intercept)
    isfinite(_family_loss(family, y, eta)) ||
        error("The initial model has a nonfinite loss; rescale the data")
    for iteration in 1:budget
        old_intercept = intercept
        old_beta = copy(beta)
        mu = _family_mean(family, eta)
        weights = _family_weights(family, mu)
        current_loss = _family_loss(family, y, eta)
        intercept_step = sum(y .- mu) / sum(weights)
        fraction = 1.0
        accepted = false
        for backtrack in 1:40
            candidate_eta = eta .+ fraction * intercept_step
            candidate_loss = _family_loss(family, y, candidate_eta)
            rounding_tolerance = 8 * Base.eps(Float64) * max(1.0, abs(current_loss))
            if isfinite(candidate_loss) && candidate_loss <= current_loss + rounding_tolerance
                intercept += fraction * intercept_step
                eta = candidate_eta
                current_loss = candidate_loss
                accepted = true
                break
            end
            fraction /= 2
        end
        blocked = !accepted

        for j in columns
            x = view(X, :, j)
            mu = _family_mean(family, eta)
            weights = _family_weights(family, mu)
            h = dot(weights, abs2.(x)) / n
            h > 0 || continue
            gradient = dot(x, y .- mu) / n
            lambda1 = lam * factors[j] * alpha
            lambda2 = lam * factors[j] * (1 - alpha)
            # z is the unscaled partial score, not the Newton coefficient.
            z = h * beta[j] + gradient
            candidate = _penalty_coordinate(z, h, lambda1, lambda2, penalty, gamma)
            shift = candidate - beta[j]
            shift == 0 && continue
            old_penalty = _penalty_value(beta[j], lambda1, lambda2, penalty, gamma)
            accepted = false
            fraction = 1.0
            for backtrack in 1:40
                candidate_beta = beta[j] + fraction * shift
                candidate_eta = eta .+ (fraction * shift) .* x
                candidate_loss = _family_loss(family, y, candidate_eta)
                penalty_change =
                    _penalty_value(candidate_beta, lambda1, lambda2, penalty, gamma) -
                    old_penalty
                rounding_tolerance =
                    8 * Base.eps(Float64) * max(1.0, abs(current_loss), abs(old_penalty))
                if isfinite(candidate_loss) &&
                   candidate_loss - current_loss + penalty_change <= rounding_tolerance
                    beta[j] = candidate_beta
                    eta = candidate_eta
                    current_loss = candidate_loss
                    accepted = true
                    break
                end
                fraction /= 2
            end
            blocked |= !accepted
        end

        violation = _stationarity(
            X, y, family, eta, beta, factors, lam, alpha, penalty, gamma, columns,
        )
        change = max(maximum(abs.(beta .- old_beta)), abs(intercept - old_intercept))
        if violation <= tolerance &&
           change <= tolerance * max(1.0, maximum(abs.(beta)), abs(intercept))
            return intercept, iteration, :converged
        elseif !(family isa GaussianFamily) && null_deviance > 0 &&
               _family_deviance(family, y, _family_mean(family, eta)) < 0.01 * null_deviance
            return intercept, iteration, :saturated
        elseif blocked && change == 0
            return intercept, iteration, :stalled
        end
    end
    return intercept, budget, :max_iter
end

function _locally_convex(X, family, eta, beta, factors, lam, alpha, penalty, gamma)
    active = findall(j -> beta[j] != 0 || factors[j] == 0, eachindex(beta))
    design = hcat(ones(size(X, 1)), X[:, active])
    weights = _family_weights(family, _family_mean(family, eta))
    hessian = design' * (weights .* design) / size(X, 1)
    for (k, j) in enumerate(active)
        lambda1 = lam * factors[j] * alpha
        curvature = lam * factors[j] * (1 - alpha)
        if lambda1 > 0
            value = abs(beta[j])
            if penalty == :MCP && value < gamma * lambda1
                curvature -= 1 / gamma
            elseif penalty == :SCAD && lambda1 < value < gamma * lambda1
                curvature -= 1 / (gamma - 1)
            end
        end
        hessian[k + 1, k + 1] += curvature
    end
    return eigmin(Symmetric(hessian)) > 0
end

"""
    ncvreg(X, y; family=GaussianFamily(), penalty=:SCAD, lambda=nothing,
           gamma=nothing, alpha=1, penalty_factor=ones(size(X, 2)),
           nlambda=100, lambda_min_ratio=nothing, eps=1e-6, max_iter=10000,
           dfmax=size(X, 2), warn=true, convex=true, return_X=false)

Fit a warm-started regularization path for Gaussian, binomial, or Poisson
regression. Predictors are centered and standardized; coefficients are returned
on the original scale and the intercept is never penalized.

`lambda` accepts a nonnegative scalar or nonincreasing vector. Otherwise,
`nlambda` log-spaced values start at the null-model score divided by `alpha`,
accounting for `penalty_factor`. `lambda_min_ratio` defaults to 0.001 when
`n > p`, and 0.05 otherwise. Setting it to zero appends an unpenalized endpoint.
An all-zero score produces an all-zero path. Automatic paths require `alpha > 0`;
pure ridge (`alpha=0`) is available with explicit lambda values.

`penalty_factor[j]=0` fits predictor `j` without shrinkage. Nonzero factors
multiply both penalty components without normalization. `gamma` defaults to
3.7 for SCAD and 3 for MCP; `alpha` mixes the selected penalty with ridge.
SCAD/MCP are applied to standardized coefficients (without GLM adaptive rescaling).

`eps` controls coefficient-change and score convergence. `max_iter` bounds
coordinate sweeps over the entire path, including initialization of unpenalized
predictors. `dfmax` stops after retaining the first model exceeding this number
of nonzero slopes, including unpenalized ones. A nonconverged or saturated fit
also stops the path. Only attempted lambdas are returned, with `iter`, `status`,
`converged`, and `stop_reason` diagnostics. `warn=false` suppresses warnings,
not diagnostics or errors. `convex=false` skips active-set Hessian diagnostics.
`return_X=true` retains copies of original `X` and `y` for in-sample methods.
"""
function ncvreg(
    X::AbstractMatrix{<:Real},
    y::AbstractVector{<:Real};
    family::AbstractFamily=GaussianFamily(),
    penalty=:SCAD,
    lambda=nothing,
    gamma=nothing,
    alpha::Real=1.0,
    penalty_factor=ones(size(X, 2)),
    nlambda::Integer=100,
    lambda_min_ratio=nothing,
    eps::Real=1e-6,
    max_iter::Integer=10000,
    dfmax::Integer=size(X, 2),
    warn::Bool=true,
    convex::Bool=true,
    return_X::Bool=false,
)
    n, p = size(X)
    n > 0 && p > 0 || throw(ArgumentError("X must have at least one row and column"))
    length(y) == n || throw(DimensionMismatch("length(y) must equal size(X, 1)"))
    _validate_response(family, y)
    mixing = Float64(alpha)
    tolerance = Float64(eps)
    isfinite(mixing) && 0 <= mixing <= 1 ||
        throw(ArgumentError("alpha must be in [0, 1]"))
    isfinite(tolerance) && tolerance > 0 ||
        throw(ArgumentError("eps must be finite and positive"))
    max_iter > 0 || throw(ArgumentError("max_iter must be positive"))
    dfmax >= 0 || throw(ArgumentError("dfmax must be nonnegative"))
    pen = _penalty_symbol(penalty)
    g = isnothing(gamma) ? (pen == :SCAD ? 3.7 : 3.0) : Float64(gamma)
    isfinite(g) || throw(ArgumentError("gamma must be finite"))
    pen == :SCAD && g <= 2 && throw(ArgumentError("SCAD gamma must be greater than 2"))
    pen == :MCP && g <= 1 && throw(ArgumentError("MCP gamma must be greater than 1"))
    penalty_factor isa AbstractVector{<:Real} ||
        throw(ArgumentError("penalty_factor must be a numeric vector"))
    factors = Float64.(penalty_factor)
    length(factors) == p ||
        throw(DimensionMismatch("penalty_factor must have one value per column of X"))
    all(v -> isfinite(v) && v >= 0, factors) ||
        throw(ArgumentError("penalty_factor values must be finite and nonnegative"))

    ratio = isnothing(lambda_min_ratio) ? (n > p ? 0.001 : 0.05) :
            Float64(lambda_min_ratio)
    supplied_lambda = isnothing(lambda) ? nothing : _validate_lambda(lambda)
    if isnothing(supplied_lambda)
        mixing > 0 ||
            throw(ArgumentError("automatic lambda paths require alpha > 0; supply lambda for ridge"))
        _validate_path_options(nlambda, ratio)
    end

    Xs, xmean, xscale = _standardize_predictors(X)
    yy = Float64.(y)
    _validate_response(family, yy)
    ymean = mean(yy)
    intercept = family isa GaussianFamily ? ymean :
                family isa BinomialFamily ? log(ymean / (1 - ymean)) : log(ymean)
    beta = zeros(Float64, p)
    unpenalized = findall(iszero, factors)
    setup_iter = 0
    if isnothing(supplied_lambda) && !isempty(unpenalized)
        intercept, setup_iter, setup_status = _fit_lambda!(
            Xs, yy, family, beta, intercept, factors, 0.0, mixing, pen, g,
            tolerance, max_iter, unpenalized,
        )
        setup_status == :converged || error(
            "Unpenalized initialization did not converge ($setup_status); " *
            "increase max_iter or check for separation/ill-conditioning",
        )
        setup_iter < max_iter || error(
            "max_iter was exhausted during unpenalized initialization; increase max_iter",
        )
    end
    null_eta = _linear_predictor(Xs, beta, intercept)
    null_mu = _family_mean(family, null_eta)
    lambdas = if isnothing(supplied_lambda)
        setup_lambda(Xs, yy, null_mu, factors, mixing, nlambda, ratio)
    else
        supplied_lambda
    end
    null_deviance = _family_deviance(family, yy, null_mu)
    coefficients = zeros(Float64, p + 1, length(lambdas))
    deviances = zeros(Float64, length(lambdas))
    converged = falses(length(lambdas))
    iterations = zeros(Int, length(lambdas))
    statuses = fill(:unattempted, length(lambdas))
    local_convexity = convex ? falses(length(lambdas)) : nothing
    total_iter = setup_iter
    stop_reason = :completed
    last = 0
    for (l, lam) in enumerate(lambdas)
        intercept, iterations[l], status = _fit_lambda!(
            Xs, yy, family, beta, intercept, factors, lam, mixing, pen, g,
            tolerance, max_iter - total_iter, 1:p; null_deviance,
        )
        total_iter += iterations[l]
        eta = _linear_predictor(Xs, beta, intercept)
        deviances[l] = _family_deviance(family, yy, _family_mean(family, eta))
        coefficients[1, l] = intercept - dot(xmean ./ xscale, beta)
        coefficients[2:end, l] = beta ./ xscale
        converged[l] = status == :converged
        if status == :converged && count(!iszero, beta) > dfmax
            status = :dfmax
        elseif status == :converged && !(family isa GaussianFamily) &&
               null_deviance > 0 && deviances[l] < 0.01 * null_deviance
            status = :saturated
        end
        statuses[l] = status
        if convex
            local_convexity[l] = _locally_convex(
                Xs, family, eta, beta, factors, lam, mixing, pen, g,
            )
        end
        last = l
        if status != :converged
            stop_reason = status
            break
        elseif total_iter == max_iter && l < length(lambdas)
            stop_reason = :max_iter
            break
        end
    end
    if warn && stop_reason != :completed
        @warn "Regularization path stopped early" reason=stop_reason index=last lambda=lambdas[last]
    end
    local_convexity = convex ? local_convexity[1:last] : nothing
    first_nonconvex = convex ? findfirst(!, local_convexity) : nothing
    convex_min = convex ? (isnothing(first_nonconvex) ? last : first_nonconvex - 1) : nothing
    return RegressionPath(
        coefficients[:, 1:last], lambdas[1:last], family, pen, deviances[1:last],
        converged[1:last], iterations[1:last], statuses[1:last], stop_reason,
        setup_iter, g, mixing, factors, tolerance, Int(max_iter), Int(dfmax), n,
        local_convexity, convex_min, return_X ? Matrix{Float64}(X) : nothing,
        return_X ? copy(yy) : nothing,
    )
end

export ncvreg
