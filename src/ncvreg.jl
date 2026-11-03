function ncvreg(X, y, lambda, penalty, gamma; kwargs...)
    return ncvreg(X, y; lambda, penalty, gamma, kwargs...).beta
end
