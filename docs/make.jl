using Documenter, NCVREG

makedocs(;
    modules=[NCVREG],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/pnavaro/NCVREG.jl/blob/{commit}{path}#L{line}",
    sitename="NCVREG.jl",
    authors="Pierre Navaro, Institut de Recherche Mathématique de Rennes",
    assets=String[],
)

deploydocs(;
    repo="github.com/pnavaro/NCVREG.jl",
)
