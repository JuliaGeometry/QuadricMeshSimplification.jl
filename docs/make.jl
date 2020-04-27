using Documenter, QuadricMeshSimplification

makedocs(;
    modules=[QuadricMeshSimplification],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/sjkelly/QuadricMeshSimplification.jl/blob/{commit}{path}#L{line}",
    sitename="QuadricMeshSimplification.jl",
    authors="Steve Kelly",
    assets=String[],
)

deploydocs(;
    repo="github.com/sjkelly/QuadricMeshSimplification.jl",
)
