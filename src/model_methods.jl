"""
    fitted(path::RegressionPath, X; index=length(path.lambda), type=:response)

Return fitted responses (`type=:response`) or linear predictors (`type=:link`)
for a lambda-path index.
"""
function fitted(
    path::RegressionPath,
    X::AbstractMatrix;
    index=length(path.lambda),
    type=:response,
)
    _check_prediction_inputs(path, X, index)
    eta = _linear_predictor(X, path.beta[2:end, index], path.beta[1, index])
    type == :link && return eta
    type == :response && return _response(path.family, eta)
    throw(ArgumentError("type must be :response or :link"))
end

function _check_prediction_inputs(path::RegressionPath, X, index)
    size(X, 2) + 1 == size(path.beta, 1) ||
        throw(DimensionMismatch("X must have $(size(path.beta, 1) - 1) columns"))
    1 <= index <= length(path.lambda) || throw(BoundsError(path.lambda, index))
    all(isfinite, X) || throw(ArgumentError("X must contain only finite values"))
    return nothing
end

"""
    predict(path::RegressionPath, X; index=length(path.lambda), type=:response)

Predict responses (`type=:response`) or linear predictors (`type=:link`).
Binomial response predictions are probabilities; Poisson response predictions
are expected counts.
"""
function predict(
    path::RegressionPath,
    X::AbstractMatrix;
    index=length(path.lambda),
    type=:response,
)
    _check_prediction_inputs(path, X, index)
    eta = fitted(path, X; index, type=:link)
    type == :link && return eta
    type == :response && return _response(path.family, eta)
    throw(ArgumentError("type must be :response or :link"))
end

"""
    residuals(path::RegressionPath, X, y; index=length(path.lambda), type=:response)

Compute response, Pearson, or signed deviance residuals at one path index.
"""
function residuals(
    path::RegressionPath,
    X::AbstractMatrix,
    y::AbstractVector;
    index=length(path.lambda),
    type=:response,
)
    length(y) == size(X, 1) ||
        throw(DimensionMismatch("length(y) must equal size(X, 1)"))
    _validate_response(path.family, y)
    mu = predict(path, X; index)

    type == :response && return Float64.(y) .- mu

    if type == :pearson
        variance = if path.family isa GaussianFamily
            ones(length(y))
        elseif path.family isa BinomialFamily
            mu .* (1 .- mu)
        else
            mu
        end
        return (Float64.(y) .- mu) ./ sqrt.(max.(variance, eps(Float64)))
    elseif type == :deviance
        return _deviance_residuals(path.family, Float64.(y), mu)
    end
    throw(ArgumentError("type must be :response, :pearson, or :deviance"))
end

function _deviance_residuals(::GaussianFamily, y, mu)
    return y .- mu
end

function _deviance_residuals(::BinomialFamily, y, mu)
    residual = similar(mu)
    floor = eps(Float64)
    for i in eachindex(y)
        probability = clamp(mu[i], floor, 1 - floor)
        contribution = 0.0
        y[i] > 0 && (contribution += y[i] * log(y[i] / probability))
        y[i] < 1 && (contribution += (1 - y[i]) * log((1 - y[i]) / (1 - probability)))
        residual[i] = sign(y[i] - probability) * sqrt(2 * contribution)
    end
    return residual
end

function _deviance_residuals(::PoissonFamily, y, mu)
    residual = similar(mu)
    for i in eachindex(y)
        mean_i = max(mu[i], floatmin(Float64))
        contribution = y[i] > 0 ?
            y[i] * log(y[i] / mean_i) - (y[i] - mean_i) : mean_i
        residual[i] = sign(y[i] - mean_i) * sqrt(2 * contribution)
    end
    return residual
end

"""
    deviance(path::RegressionPath)

Return the model deviance at each lambda in the fitted path.
"""
deviance(path::RegressionPath) = copy(path.deviance)

"""
    deviance(path::RegressionPath, X, y; index=length(path.lambda))

Compute deviance for new observations at one path index.
"""
function deviance(path::RegressionPath, X::AbstractMatrix, y::AbstractVector; index=length(path.lambda))
    length(y) == size(X, 1) ||
        throw(DimensionMismatch("length(y) must equal size(X, 1)"))
    _validate_response(path.family, y)
    mu = predict(path, X; index)
    return _family_deviance(path.family, Float64.(y), mu)
end

"""
    loglikelihood(path::RegressionPath, X, y; index=length(path.lambda))

Evaluate the response log-likelihood at one path index.
"""
function loglikelihood(
    path::RegressionPath,
    X::AbstractMatrix,
    y::AbstractVector;
    index=length(path.lambda),
)
    length(y) == size(X, 1) ||
        throw(DimensionMismatch("length(y) must equal size(X, 1)"))
    _validate_response(path.family, y)
    mu = predict(path, X; index)
    yy = Float64.(y)
    if path.family isa GaussianFamily
        rss = sum(abs2, yy .- mu)
        variance = max(rss / length(yy), floatmin(Float64))
        return -length(yy) / 2 * (log(2pi * variance) + 1)
    elseif path.family isa BinomialFamily
        probability = clamp.(mu, eps(Float64), 1 - eps(Float64))
        return sum(yy .* log.(probability) .+ (1 .- yy) .* log1p.(-probability))
    end

    mean_value = max.(mu, floatmin(Float64))
    return sum(yy .* log.(mean_value) .- mean_value) -
           sum(_logfactorial, yy)
end

function _logfactorial(value::Real)
    n = Int(value)
    return sum(log, 2:n; init=0.0)
end

"""
    summary(path::RegressionPath; index=length(path.lambda))

Return model-family, penalty, degrees-of-freedom, deviance, and convergence
information for one path index.
"""
function summary(path::RegressionPath; index=length(path.lambda))
    1 <= index <= length(path.lambda) || throw(BoundsError(path.lambda, index))
    coefficients = path.beta[2:end, index]
    return (
        family=string(nameof(typeof(path.family))),
        penalty=path.penalty,
        lambda=path.lambda[index],
        degrees_of_freedom=count(!iszero, coefficients),
        deviance=path.deviance[index],
        converged=path.converged[index],
        iterations=path.iter[index],
        status=path.status[index],
        stop_reason=path.stop_reason,
        gamma=path.gamma,
        alpha=path.alpha,
        locally_convex=isnothing(path.locally_convex) ? nothing : path.locally_convex[index],
        convex_min=path.convex_min,
    )
end

function Base.show(io::IO, path::RegressionPath)
    print(
        io,
        "RegressionPath(",
        nameof(typeof(path.family)),
        ", penalty=:",
        path.penalty,
        ", ",
        length(path.lambda),
        " lambda values, ",
        count(path.converged),
        " converged, stop_reason=:",
        path.stop_reason,
        ")",
    )
end

"""
    devianceplot(path)

Create a RecipesBase-compatible diagnostic plot of deviance against lambda.
"""
struct DeviancePlot{P<:RegressionPath}
    path::P
end

devianceplot(path::RegressionPath) = DeviancePlot(path)

"""
    residualplot(path, X, y; index=length(path.lambda), type=:deviance)

Create a RecipesBase-compatible residual-versus-fitted diagnostic plot.
"""
struct ResidualPlot{P<:RegressionPath}
    path::P
    fitted::Vector{Float64}
    residual::Vector{Float64}
    type::Symbol
end

function residualplot(
    path::RegressionPath,
    X::AbstractMatrix,
    y::AbstractVector;
    index=length(path.lambda),
    type=:deviance,
)
    fitted_values = predict(path, X; index)
    residual = residuals(path, X, y; index, type)
    return ResidualPlot(path, Float64.(fitted_values), Float64.(residual), type)
end

export fitted, residuals, deviance, loglikelihood, summary
export DeviancePlot, ResidualPlot, devianceplot, residualplot

function _training_data(path::RegressionPath)
    isnothing(path.X) && throw(ArgumentError(
        "training data were not retained; pass X and y explicitly or fit with return_X=true",
    ))
    return path.X, path.y
end

"""
    fitted(path; kwargs...)
    predict(path; kwargs...)
    residuals(path; kwargs...)
    loglikelihood(path; kwargs...)
    residualplot(path; kwargs...)

Use the training data retained by `ncvreg(...; return_X=true)`.
"""
function fitted(path::RegressionPath; kwargs...)
    X, _ = _training_data(path)
    return fitted(path, X; kwargs...)
end

function predict(path::RegressionPath; kwargs...)
    X, _ = _training_data(path)
    return predict(path, X; kwargs...)
end

function residuals(path::RegressionPath; kwargs...)
    X, y = _training_data(path)
    return residuals(path, X, y; kwargs...)
end

function loglikelihood(path::RegressionPath; kwargs...)
    X, y = _training_data(path)
    return loglikelihood(path, X, y; kwargs...)
end

function residualplot(path::RegressionPath; kwargs...)
    X, y = _training_data(path)
    return residualplot(path, X, y; kwargs...)
end
