function _validate_lambda(lambda)
    lambda isa Real || lambda isa AbstractVector{<:Real} ||
        throw(ArgumentError("lambda must be a real scalar or numeric vector"))
    values = lambda isa Real ? [Float64(lambda)] : Float64.(lambda)
    !isempty(values) || throw(ArgumentError("lambda must not be empty"))
    all(v -> isfinite(v) && v >= 0, values) ||
        throw(ArgumentError("lambda values must be finite and nonnegative"))
    all(diff(values) .<= 0) ||
        throw(ArgumentError("lambda values must be in nonincreasing order"))
    return values
end

function _validate_path_options(nlambda, ratio)
    nlambda >= 1 || throw(ArgumentError("nlambda must be positive"))
    isfinite(ratio) && 0 <= ratio <= 1 ||
        throw(ArgumentError("lambda_min_ratio must be in [0, 1]"))
    ratio == 0 && nlambda < 2 && throw(ArgumentError(
        "nlambda must be at least 2 when lambda_min_ratio=0",
    ))
end

function setup_lambda(X, y, null_mu, factors, alpha, nlambda, ratio)
    penalized = findall(factors .> 0)
    score = X' * (y .- null_mu) / length(y)
    lambda_max = isempty(penalized) ? 0.0 :
                 maximum(abs.(score[penalized]) ./ factors[penalized]) / alpha
    isfinite(lambda_max) ||
        throw(ArgumentError("lambda_max is not finite; rescale data or penalty_factor"))
    lambda_max == 0 && return zeros(Float64, nlambda)
    nlambda == 1 && return [lambda_max]
    positive_count = ratio == 0 ? nlambda - 1 : nlambda
    endpoint = ratio == 0 ? 0.001 : ratio
    lambdas = positive_count == 1 ? [lambda_max] : exp.(range(
        log(lambda_max), log(lambda_max) + log(endpoint); length=positive_count,
    ))
    # Keep the null-model boundary exact instead of round-tripping it through log.
    lambdas[1] = lambda_max
    ratio == 0 && push!(lambdas, 0.0)
    return lambdas
end
