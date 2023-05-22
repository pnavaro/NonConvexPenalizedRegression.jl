"""
   crossprod( X, y, n, j)

Cross product of y with jth column of X
"""
function crossprod(X, y, j)

    val::Float64 = 0.0

    for i in eachindex(y)
        val += X[i, j] * y[i]
    end

    return val

end
