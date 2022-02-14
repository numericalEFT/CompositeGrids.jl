using CompositeGrids
using Documenter

DocMeta.setdocmeta!(CompositeGrids, :DocTestSetup, :(using CompositeGrids); recursive=true)

makedocs(;
    modules=[CompositeGrids],
    authors="Kun Chen, Tao Wang, Xiansheng Cai",
    repo="https://github.com/numericalEFT/CompositeGrids.jl/blob/{commit}{path}#{line}",
    sitename="CompositeGrids.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://numericaleft.github.io/CompositeGrids.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
        ],
        "Library" => Any[
            "lib/simple.md",
            "lib/composite.md",
            "lib/interpolate.md",
                # map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
                # "Internals" => map(s -> "lib/$(s)", sort(readdir(joinpath(@__DIR__, "src/lib"))))
        ]
    ],
)

deploydocs(;
    repo="github.com/numericalEFT/CompositeGrids.jl",
    devbranch="dev"
)
