RecipesBase.@recipe function f(path::RegressionPath)
    seriestype := :path
    xscale --> (all(value -> value > 0, path.lambda) ? :log10 : :identity)
    xlabel --> "Lambda"
    ylabel --> "Coefficient"
    title --> "$(nameof(typeof(path.family))) $(path.penalty) coefficient paths"
    legend --> :best

    for j in axes(path.beta, 1)
        RecipesBase.@series begin
            label --> (j == 1 ? "Intercept" : "β$(j - 1)")
            path.lambda, path.beta[j, :]
        end
    end
end

RecipesBase.@recipe function f(diagnostic::DeviancePlot)
    path = diagnostic.path
    seriestype := :line
    xscale --> (all(value -> value > 0, path.lambda) ? :log10 : :identity)
    xlabel --> "Lambda"
    ylabel --> "Deviance"
    title --> "$(nameof(typeof(path.family))) model deviance"
    path.lambda, path.deviance
end

RecipesBase.@recipe function f(diagnostic::ResidualPlot)
    seriestype := :scatter
    xlabel --> "Fitted response"
    ylabel --> "$(diagnostic.type) residual"
    title --> "$(nameof(typeof(diagnostic.path.family))) residuals"
    RecipesBase.@series begin
        seriestype := :line
        label --> false
        color --> :black
        linestyle --> :dash
        x = collect(extrema(diagnostic.fitted))
        x, [0.0, 0.0]
    end
    diagnostic.fitted, diagnostic.residual
end
