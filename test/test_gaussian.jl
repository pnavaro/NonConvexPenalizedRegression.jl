using Test

@testset "Gaussian (linear regression)" begin

    using NonConvexPenalizedRegression
    using RData
    
    data = load(joinpath(@__DIR__, "data/Prostate.rda"))
    X  = data["Prostate"]["X"]
    y  = data["Prostate"]["y"]
  
    λ = [0.0]
    mcp = MCP(X, y, λ)
    
    λ = [0.0]
    scad = SCAD(X, y, λ)
    
    # @test maximum(abs.(beta .- scad.beta)) < 2e-5
    # 
    # yp = XX * scad.beta 
    # 
    # rmse = sqrt(mean(abs2.(y .- yp)))
    # println("SCAD rmse = $rmse")

    # 
    # @test maximum(abs.(beta .- mcp.beta)) < 2e-5
    # 
    # yp = XX * mcp.beta 
    # 
    # rmse = sqrt(mean(abs2.(y .- yp)))
    # println("MCP rmse = $rmse")

    # lasso = Lasso(X, y, λ)
    # 
    # @test maximum(abs.(beta .- lasso.beta)) < 2e-5
    # 
    # yp = XX * lasso.beta 
    # 
    # rmse = sqrt(mean(abs2.(y .- yp)))
    # println("Lasso rmse = $rmse")

end
