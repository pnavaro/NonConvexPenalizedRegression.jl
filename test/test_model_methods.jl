using NonConvexPenalizedRegression
using Random
using RecipesBase

@testset "Model methods, diagnostics, and plot recipes" begin
    rng = MersenneTwister(23)
    X = randn(rng, 120, 2)

    y = 0.4 .+ X * [0.7, -0.5] .+ 0.2 .* randn(rng, 120)
    gaussian = ncvreg(
        X, y; family=GaussianFamily(), penalty=:Lasso, lambda=[0.1, 0.03]
    )
    @test fitted(gaussian, X) ≈ predict(gaussian, X)
    @test fitted(gaussian, X; type=:link) ≈ predict(gaussian, X; type=:link)
    @test predict(gaussian, X) ≈ fitted(gaussian, X)
    @test residuals(gaussian, X, y) ≈ y .- predict(gaussian, X)
    @test all(isfinite, residuals(gaussian, X, y; type=:pearson))
    @test all(isfinite, residuals(gaussian, X, y; type=:deviance))
    @test deviance(gaussian) ≈ gaussian.deviance
    @test deviance(gaussian, X, y) ≈ gaussian.deviance[end]
    @test isfinite(loglikelihood(gaussian, X, y))
    @test summary(gaussian).degrees_of_freedom >= 0
    @test occursin("RegressionPath", sprint(show, gaussian))
    @test length(RecipesBase.apply_recipe(Dict{Symbol,Any}(), gaussian)) == 3
    @test length(RecipesBase.apply_recipe(Dict{Symbol,Any}(), devianceplot(gaussian))) == 1
    @test length(
        RecipesBase.apply_recipe(
            Dict{Symbol,Any}(), residualplot(gaussian, X, y)
        )
    ) == 2
    residual_series = RecipesBase.apply_recipe(
        Dict{Symbol,Any}(), residualplot(gaussian, X, y),
    )
    @test residual_series[1].args[1] == collect(extrema(predict(gaussian, X)))
    @test residual_series[1].args[2] == [0.0, 0.0]

    y_binomial = Float64.(X[:, 1] .+ 0.3 .* randn(rng, 120) .> 0)
    binomial = ncvreg(
        X, y_binomial; family=BinomialFamily(), penalty=:Lasso, lambda=[0.1]
    )
    @test all(isfinite, residuals(binomial, X, y_binomial; type=:deviance))
    @test isfinite(loglikelihood(binomial, X, y_binomial))

    y_poisson = Float64.(round.(Int, exp.(0.1 .+ 0.2 .* X[:, 1])))
    poisson = ncvreg(
        X, y_poisson; family=PoissonFamily(), penalty=:MCP, lambda=[0.1]
    )
    @test all(isfinite, residuals(poisson, X, y_poisson; type=:pearson))
    @test all(isfinite, residuals(poisson, X, y_poisson; type=:deviance))
    @test isfinite(loglikelihood(poisson, X, y_poisson))

    @test_throws ArgumentError predict(gaussian, X; type=:invalid)
    @test_throws ArgumentError residuals(gaussian, X, y; type=:invalid)
    @test_throws DimensionMismatch predict(gaussian, randn(5, 3))
end
