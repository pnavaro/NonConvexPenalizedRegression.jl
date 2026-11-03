using NonConvexPenalizedRegression
using Random

@testset "Binomial and Poisson regularization paths" begin
    rng = MersenneTwister(11)
    X = randn(rng, 300, 3)

    eta_binomial = -0.3 .+ X * [0.8, -0.6, 0.0]
    probability = 1.0 ./ (1.0 .+ exp.(-eta_binomial))
    y_binomial = Float64.(rand(rng, 300) .< probability)

    binomial_fit = ncvreg(
        X,
        y_binomial;
        family=BinomialFamily(),
        penalty=:Lasso,
        lambda=[0.1, 0.05, 0.02],
    )
    @test size(coef(binomial_fit)) == (4, 3)
    @test all(binomial_fit.converged)
    probabilities = predict(binomial_fit, X)
    @test all((probabilities .>= 0) .& (probabilities .<= 1))
    @test binomial_fit.deviance[end] < binomial_fit.deviance[1]

    eta_poisson = 0.2 .+ X * [0.25, -0.15, 0.0]
    y_poisson = Float64.(round.(Int, exp.(eta_poisson)))

    poisson_fit = ncvreg(
        X,
        y_poisson;
        family=PoissonFamily(),
        penalty=:MCP,
        lambda=[0.1, 0.05, 0.02],
    )
    @test size(coef(poisson_fit)) == (4, 3)
    @test all(poisson_fit.converged)
    @test all(predict(poisson_fit, X) .> 0)
    @test poisson_fit.deviance[end] < poisson_fit.deviance[1]

    automatic_path = ncvreg(
        X,
        y_binomial;
        family=BinomialFamily(),
        penalty=:Lasso,
        nlambda=5,
    )
    @test length(automatic_path.lambda) == 5
    @test all(diff(automatic_path.lambda) .<= 0)

    @test_throws ArgumentError ncvreg(
        X, fill(0.5, size(X, 1)); family=BinomialFamily()
    )
    @test_throws ArgumentError ncvreg(
        X, fill(-1.0, size(X, 1)); family=PoissonFamily()
end
