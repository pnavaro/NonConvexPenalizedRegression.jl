module NonConvexPenalizedRegression

using Statistics
using DocStringExtensions

ncol(A::AbstractMatrix) = size(A)[2]

include("utils.jl")
include("standardize.jl")
include("gaussian_loss.jl")
include("scad.jl")
include("mcp.jl")
include("lasso.jl")
include("gaussian.jl")
include("binomial.jl")
include("poisson.jl")
include("ncvreg.jl")

end
