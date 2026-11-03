using LinearAlgebra
using NonConvexPenalizedRegression
using Random
using Statistics

@testset "Lambda setup and fitting options" begin
    X = [-1.0 -1.0; -1.0 1.0; 1.0 -1.0; 1.0 1.0]
    y = 0.3 .+ X * [2.0, 0.4]

    @testset "Automatic and explicit paths" begin
        for penalty in (:Lasso, :MCP, :SCAD)
            fit = ncvreg(X, y; penalty, nlambda=7)
            @test length(fit.lambda) == 7
            @test fit.lambda[1] ≈ 2
            @test fit.lambda[end] ≈ 0.002
            @test all(iszero, fit.beta[2:end, 1])
            @test all(fit.converged)
            @test all(diff(log.(fit.lambda)) .≈ log(0.001) / 6)
            @test fit.stop_reason == :completed
        end
        high_dimensional = ncvreg(hcat(X, X), y; nlambda=3, penalty=:Lasso)
        @test high_dimensional.lambda[end] ≈ 0.05 * high_dimensional.lambda[1]

        one = ncvreg(X, y; nlambda=1, penalty=:Lasso)
        @test one.lambda == [2.0]
        zero_end = ncvreg(X, y; nlambda=5, lambda_min_ratio=0, penalty=:Lasso)
        @test length(zero_end.lambda) == 5
        @test zero_end.lambda[end] == 0
        @test zero_end.lambda[end - 1] ≈ 0.002
        @test zero_end.beta[:, end] ≈ [0.3, 2.0, 0.4]
        two = ncvreg(X, y; nlambda=2, lambda_min_ratio=0, penalty=:Lasso)
        @test two.lambda == [2.0, 0.0]
        equal = ncvreg(X, y; nlambda=3, lambda_min_ratio=1, penalty=:Lasso)
        @test equal.lambda ≈ fill(2.0, 3)
        constant = ncvreg(X, fill(2.0, 4); nlambda=4)
        @test constant.lambda == zeros(4)
        @test all(constant.converged)
        scalar = ncvreg(X, y; lambda=0.2, penalty="lasso")
        @test scalar.lambda == [0.2]
        @test scalar.beta[:, 1] ≈ [0.3, 1.8, 0.2]
        supplied = ncvreg(X, y; lambda=[0.8, 0.8, 0.0], penalty=:Lasso)
        @test supplied.lambda == [0.8, 0.8, 0.0]
        @test supplied.beta[:, 1] ≈ supplied.beta[:, 2]
    end

    @testset "Penalty mixing, weights, and shape" begin
        factors = [2.0, 0.5]
        lambda = [0.6, 0.2, 0.0]
        fit = ncvreg(X, y; penalty=:Lasso, alpha=0.4, penalty_factor=factors, lambda)
        expected = hcat([
            vcat(0.3, max.(abs.([2.0, 0.4]) .- lam .* factors .* 0.4, 0) ./
                 (1 .+ lam .* factors .* 0.6))
            for lam in lambda
        ]...)
        @test fit.beta ≈ expected atol=1e-8
        @test fit.alpha == 0.4
        @test fit.penalty_factor == factors
        factors[1] = 100
        @test fit.penalty_factor[1] == 2
        weighted = ncvreg(
            X, y; penalty=:Lasso, alpha=0.5, penalty_factor=[2.0, 1.0], nlambda=4,
        )
        @test weighted.lambda[1] ≈ 2
        @test all(iszero, weighted.beta[2:end, 1])

        for penalty in (:Lasso, :MCP, :SCAD)
            ridge = ncvreg(X, y; penalty, alpha=0, lambda=[2.0, 0.5, 0.0])
            @test ridge.beta[2:end, :] ≈ [2.0, 0.4] ./ (1 .+ [2.0 0.5 0.0])
            @test ridge.beta[1, :] ≈ fill(0.3, 3)
        end
        mcp = ncvreg(X, y; penalty=:MCP, gamma=3.0, lambda=[1.0])
        @test mcp.beta[2, 1] ≈ 1.5
        scad = ncvreg(X, y; penalty=:SCAD, gamma=3.7, lambda=[0.8])
        @test scad.beta[2, 1] ≈ (2.0 - 3.7 * 0.8 / 2.7) / (1 - 1 / 2.7)
        @test scad.gamma == 3.7
        for h in (0.1, 0.5, 1.0, 3.0), penalty in (:MCP, :SCAD)
            for z in (0.0, 0.3, 1.5, 3.0)
                update = NonConvexPenalizedRegression._penalty_coordinate(
                    z, h, 0.7, 0.2, penalty, 3.7,
                )
                objective(b) = h * b^2 / 2 - z * b +
                    NonConvexPenalizedRegression._penalty_value(b, 0.7, 0.2, penalty, 3.7)
                @test objective(update) <= minimum(objective.(0:0.005:15)) + 1e-10
            end
        end
    end

    @testset "Unpenalized initialization for every family" begin
        rng = MersenneTwister(100)
        x = randn(rng, 200)
        design = hcat(x, 0.8 .* x .+ 0.6 .* randn(rng, 200))
        binomial_y = Float64.(rand(rng, 200) .< 1 ./ (1 .+ exp.(-0.2 .- x)))
        poisson_y = Float64.(round.(Int, exp.(0.2 .+ 0.3 .* x)))
        gaussian_y = 0.3 .+ 2 .* x .+ 0.1 .* randn(rng, 200)
        for (family, response) in (
            (GaussianFamily(), gaussian_y),
            (BinomialFamily(), binomial_y),
            (PoissonFamily(), poisson_y),
        )
            baseline = ncvreg(
                design[:, 1:1], response; family, penalty=:Lasso, lambda=0, eps=1e-9,
            )
            fit = ncvreg(
                design, response; family, penalty=:Lasso, penalty_factor=[0, 2],
                alpha=0.7, nlambda=4, eps=1e-9,
            )
            standardized = (design .- mean(design; dims=1)) ./
                           sqrt.(mean((design .- mean(design; dims=1)) .^ 2; dims=1))
            baseline_residual = response .- predict(baseline, design[:, 1:1])
            expected_max = abs(dot(standardized[:, 2], baseline_residual)) /
                           (length(response) * 2 * 0.7)
            @test fit.lambda[1] ≈ expected_max atol=1e-8
            @test fit.beta[1:2, 1] ≈ baseline.beta[:, 1] atol=1e-7
            @test abs(fit.beta[3, 1]) < 1e-7
            @test fit.setup_iter > 0
            @test all(fit.converged)
            @test fit.stop_reason == :completed
            @test sum(fit.iter) + fit.setup_iter <= fit.max_iter

            full = ncvreg(design, response; family, penalty=:Lasso, lambda=0)
            free = ncvreg(
                design, response; family, penalty_factor=[0, 0], lambda=[5.0, 1.0, 0.0],
            )
            @test free.beta ≈ repeat(full.beta, 1, 3) atol=1e-5
            auto_free = ncvreg(
                design, response; family, penalty_factor=[0, 0], nlambda=3,
            )
            @test auto_free.lambda == zeros(3)
            @test auto_free.beta ≈ repeat(full.beta, 1, 3) atol=1e-5
        end
    end

    @testset "GLM convergence and zero-lambda endpoints" begin
        design = reshape(repeat([-1.0, 1.0]; inner=50), :, 1)
        responses = (
            vcat(ones(20), zeros(30), ones(35), zeros(15)),
            vcat(repeat([0.0, 1.0], 25), repeat([1.0, 2.0], 25)),
        )
        for (family, response, left, right) in (
            (BinomialFamily(), responses[1], log(0.4 / 0.6), log(0.7 / 0.3)),
            (PoissonFamily(), responses[2], log(0.5), log(1.5)),
        )
            for penalty in (:Lasso, :MCP, :SCAD)
                fit = ncvreg(design, response; family, penalty, lambda=0, eps=1e-9)
                @test all(fit.converged)
                @test fit.beta[:, 1] ≈ [(left + right) / 2, (right - left) / 2] atol=1e-7
            end
            for alpha in (1.0, 0.6, 0.0)
                fit = ncvreg(
                    design, response; family, penalty=:Lasso, alpha, lambda=0.03, eps=1e-9,
                )
                residual = response .- predict(fit, design)
                @test all(fit.converged)
                @test abs(mean(residual)) < 1e-8
                @test dot(vec(design), residual) / length(response) ≈
                      0.03 * alpha + 0.03 * (1 - alpha) * fit.beta[2, 1] atol=1e-8
            end
        end
    end

    @testset "Stopping, retention, and legacy entry points" begin
        limited = @test_logs (:warn, r"Regularization path stopped early") ncvreg(
            X, y; penalty=:Lasso, lambda=[3.0, 1.0, 0.0], max_iter=1,
        )
        @test limited.lambda == [3.0]
        @test limited.stop_reason == :max_iter
        @test sum(limited.iter) == 1
        nonconverged = @test_logs ncvreg(
            X, y; penalty=:Lasso, lambda=[0.1, 0.0], max_iter=1, warn=false,
        )
        @test nonconverged.status == [:max_iter]
        @test !any(nonconverged.converged)
        @test size(nonconverged.beta, 2) == length(nonconverged.lambda) == 1

        sparse_path = ncvreg(
            X, y; penalty=:Lasso, lambda=[3.0, 1.0, 0.1, 0.0], dfmax=1, warn=false,
        )
        @test sparse_path.stop_reason == :dfmax
        @test sparse_path.status[end] == :dfmax
        @test sparse_path.lambda == [3.0, 1.0, 0.1]
        @test summary(sparse_path).status == :dfmax
        @test summary(sparse_path).iterations == sparse_path.iter[end]
        @test sparse_path.convex_min == 3
        @test all(sparse_path.locally_convex)
        @test NonConvexPenalizedRegression._locally_convex(
            X, GaussianFamily(), fill(0.3, 4), [0.5, 0.5], ones(2),
            0.5, 1.0, :MCP, 3.0,
        )
        @test !NonConvexPenalizedRegression._locally_convex(
            hcat(X[:, 1], 0.9 .* X[:, 1] .+ 0.1 .* X[:, 2]),
            GaussianFamily(), fill(0.3, 4), [0.5, 0.5], ones(2),
            0.5, 1.0, :MCP, 3.0,
        )

        separated = ncvreg(
            X, Float64.(X[:, 1] .> 0); family=BinomialFamily(),
            penalty=:Lasso, lambda=[0.0], warn=false,
        )
        @test separated.stop_reason == :saturated
        @test separated.status == [:saturated]
        @test all(isfinite, separated.beta)
        @test sum(separated.iter) < separated.max_iter

        exhausted_setup = () -> ncvreg(
            X, y; penalty_factor=[0, 1], max_iter=1, warn=false,
        )
        @test_throws ErrorException exhausted_setup()

        retained = ncvreg(X, y; lambda=[0.5, 0.0], return_X=true, convex=false)
        @test retained.X == X
        @test retained.X !== X
        @test retained.y == y
        @test retained.y !== y
        @test isnothing(retained.convex_min)
        @test isnothing(retained.locally_convex)
        @test fitted(retained) ≈ fitted(retained, X)
        @test predict(retained) ≈ predict(retained, X)
        @test residuals(retained) ≈ residuals(retained, X, y)
        @test loglikelihood(retained) ≈ loglikelihood(retained, X, y)
        @test residualplot(retained).residual ≈ residualplot(retained, X, y).residual
        @test_throws ArgumentError fitted(limited)
        @test_throws ArgumentError residuals(limited)

        for constructor in (SCAD, MCP, Lasso)
            fit = constructor(X, y, [0.5, 0.0]; alpha=0.5)
            @test size(fit.beta) == (3, 2)
            @test fit.λ == [0.5, 0.0]
        end
        positional = ncvreg(X, y, [0.5, 0.0], :SCAD, 3.7; alpha=0.5)
        @test positional ≈ SCAD(X, y, [0.5, 0.0]; alpha=0.5).beta
    end

    @testset "Input validation" begin
        for options in (
            (; nlambda=0), (; lambda_min_ratio=-0.1), (; lambda_min_ratio=1.1),
            (; lambda_min_ratio=NaN), (; nlambda=1, lambda_min_ratio=0),
            (; alpha=0), (; alpha=-0.1), (; alpha=1.1), (; alpha=NaN),
            (; gamma=NaN), (; gamma=Inf), (; gamma=2.0, penalty=:SCAD),
            (; gamma=1.0, penalty=:MCP), (; eps=0), (; eps=Inf), (; max_iter=0),
            (; dfmax=-1), (; penalty=:invalid), (; penalty=3),
            (; penalty_factor=[1.0, -1.0]), (; penalty_factor=[1.0, Inf]),
            (; lambda=[]), (; lambda=[0.1, 0.2]), (; lambda=-0.1),
            (; lambda=NaN), (; lambda=[1.0, Inf]),
        )
            @test_throws ArgumentError ncvreg(X, y; options...)
        end
        @test_throws DimensionMismatch ncvreg(X, y; penalty_factor=[1.0])
        @test_throws DimensionMismatch ncvreg(X, y[1:3])
        @test_throws ArgumentError ncvreg(zeros(0, 2), Float64[])
    end
end
