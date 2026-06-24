using EcoEvoSim
using Documenter

DocMeta.setdocmeta!(EcoEvoSim, :DocTestSetup, :(using EcoEvoSim); recursive = true)

makedocs(;
    modules = [EcoEvoSim],
    authors = "György Barabás <dysordys@protonmail.com>",
    sitename = "EcoEvoSim.jl",
    format = Documenter.HTML(;
        canonical = "https://dysordys.github.io/EcoEvoSim",
        edit_link = "main",
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    doctest = true,
    # Only the core module's exported docstrings are collected here. Re-exported
    # solver algorithms (DynamicSS, SSRootfind, …) live in other packages, so we
    # don't enforce that every exported name is documented locally.
    checkdocs = :none,
)

deploydocs(;
    repo = "github.com/dysordys/EcoEvoSim",
    devbranch = "main",
    push_preview = true,
)
