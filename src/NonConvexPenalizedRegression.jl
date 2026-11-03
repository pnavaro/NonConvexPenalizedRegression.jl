module NonConvexPenalizedRegression

    using LinearAlgebra
    using Statistics
    using DocStringExtensions
    import RecipesBase
    import Base: summary
    
    ncol(A::AbstractMatrix) = size(A)[2]
    
    include("utils.jl")
    include("standardize.jl")
    include("gaussian_loss.jl")
    include("scad.jl")
    include("mcp.jl")
    include("lasso.jl")
    include("gaussian.jl")
    include("families.jl")
    include("setup_lambda.jl")
    include("family_path.jl")
    include("model_methods.jl")
    include("recipes.jl")
    include("ncvreg.jl")

end
