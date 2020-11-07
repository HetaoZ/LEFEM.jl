using LEFEM
using Documenter

makedocs(;
    modules=[LEFEM],
    authors="Hetao Z.",
    repo="https://github.com/HetaoZ/LEFEM.jl/blob/{commit}{path}#L{line}",
    sitename="LEFEM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://HetaoZ.github.io/LEFEM.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/HetaoZ/LEFEM.jl",
)
