using ConvexHulls2d
using Documenter

DocMeta.setdocmeta!(ConvexHulls2d, :DocTestSetup, :(using ConvexHulls2d); recursive=true)

makedocs(;
    modules=[ConvexHulls2d],
    authors="Jan Weidner <jw3126@gmail.com> and contributors",
    repo="https://github.com/jw3126/ConvexHulls2d.jl/blob/{commit}{path}#{line}",
    sitename="ConvexHulls2d.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jw3126.github.io/ConvexHulls2d.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jw3126/ConvexHulls2d.jl",
    devbranch="main",
)
