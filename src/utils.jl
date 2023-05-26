"""
$(SIGNATURES)

Cross product of y with jth column of X
"""
function crossprod(X, y, j)

    val = 0.0

    for i in eachindex(y)
        val += X[i, j] * y[i]
    end

    return val

end

"""
$(SIGNATURES)

Weighted cross product of y with jth column of x
"""
function wcrossprod(X, y, w, j)

    val = 0.0

    for i in eachindex(y)
        val += X[i, j] * y[i] * w[i]
    end

    return val

end

"""
$(SIGNATURES)
Weighted sum of squares of jth column of X
"""
function wsqsum(X, w, j)
   val=0
   for i in eachindex(w) 
      val += w[i] * X[i, j]^2
   end
   val
end

"""
$(SIGNATURES)
Sum of squares of jth column of X
"""
sqsum(X, j) = sum( view(X, :, j).^2 )

"""
$(SIGNATURES)
Gaussian loss
"""
g_loss(r) = sum( r.^2 )
