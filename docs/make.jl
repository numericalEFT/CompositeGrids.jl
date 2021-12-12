using CompositeGrids
using Documenter

DocMeta.setdocmeta!(CompositeGrids, :DocTestSetup, :(using CompositeGrids); recursive=true)

makedocs(;
    modules=[CompositeGrids],
    authors="Kun Chen, Tao Wang, Xiansheng Cai",
    repo="https://github.com/iintSjds/CompositeGrids.jl/blob/{commit}{path}#{line}",
    sitename="CompositeGrids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iintSjds.github.io/CompositeGrids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
        ],
        "Library" => Any[
                map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
                # "Internals" => map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
        ]
    ],
)

deploydocs(;
    repo="github.com/iintSjds/CompositeGrids.jl.git",
    devbranch="dev"
)
