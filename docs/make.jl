using Documenter, NonConvexPenalizedRegression

makedocs(;
    modules=[NonConvexPenalizedRegression],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pnavaro/NonConvexPenalizedRegression.jl/blob/{commit}{path}#L{line}",
    sitename="NonConvexPenalizedRegression.jl",
    authors="Pierre Navaro",
)

deploydocs(;
    repo="github.com/pnavaro/NonConvexPenalizedRegression.jl",
    branch = "gh-pages",
    devbranch = "master",
)
