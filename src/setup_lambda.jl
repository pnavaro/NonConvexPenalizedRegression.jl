using GLM
using InvertedIndices

function setup_lambda(X, y, family, alpha, lambda_min, nlambda, penalty_factor)

  n, p = size(X)

  ## Determine lambda_max

  ind = findall( penalty_factor .!= 0 )

  if length(ind) != p 
    fit = lm( X[:, Not(ind)], y, family=family)
  else 
    fit <- lm(y~1, family=family)
  end

  if (family=="gaussian") {
    zmax <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
  } else {
    zmax <- .Call("maxprod", X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  lambda.max <- zmax/alpha

  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
  }

  if (length(ind)!=p) lambda[1] <- lambda[1] * 1.000001
  lambda

end
