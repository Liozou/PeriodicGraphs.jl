using Documenter
using PeriodicGraphs, Graphs

DocMeta.setdocmeta!(PeriodicGraphs, :DocTestSetup, quote
    using PeriodicGraphs, Graphs
end; recursive=true)

makedocs(
    sitename = "PeriodicGraphs.jl",
    format = Documenter.HTML(),
    modules = [PeriodicGraphs],
    pages = [
        "Home" => "index.md",
        "Basics" => [
            "Types"        => "types.md",
            "Neighborhood" => "neighborhood.md",
        ],
        "Symmetries" => "symmetries.md",
        "Algorithms" => [
            "Dimensionality" => "dimension.md",
            "Rings"      => "rings.md",
        ],
        "Utilities"  => "utilities.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/Liozou/PeriodicGraphs.jl.git"
)
