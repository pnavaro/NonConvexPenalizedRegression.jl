using GLM
using InvertedIndices

function maxprod( X, y, v, m)

    n, p  = size(X)

    for j in eachindex(v)
      zz = crossprod(X, y, v[j]-1) / m[v[j]-1]
      if abs(zz) > z 
          z = abs(zz)
      end
    end

    return z

end

function setup_lambda(X, y, family, alpha, lambda_min, nlambda, penalty_factor)

    f = Dict("gaussian"=>Normal(), "binomial"=>Binomial(), "poisson"=>Poisson())

    n, p = size(X)

    ## Determine lambda_max

    ind = findall( penalty_factor .!= 0 )

    if length(ind) != p 
      fit = lm( X[:, Not(ind)], y, family=f[family])
    else 
      fit = lm(ones(n,1), y, family=f[family])
    end

    if family == "gaussian"
        zmax = maxprod( X, fit.residuals, ind, penalty_factor) / n
    else
        zmax = maxprod(X, residuals(fit, "working") * fit.weights, ind, penalty_factor) / n
    end

    lambda_max = zmax/alpha

    if lambda.min==0
      lambda = [exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)); 0]
    else
      lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), len=nlambda))
    end

    if (length(ind)!=p) lambda[1] = lambda[1] * 1.000001
  
    return lambda

end
