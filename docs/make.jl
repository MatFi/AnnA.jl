using Documenter
using AnnABase

makedocs(
    sitename = "AnnABase Documentation",
    modules = [AnnABase],
    doctest = true,
    strict = true,
    checkdocs = :exports,
    repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}",
    pages = Any["Home"=>"index.md", "API"=>"api.md"],
)
