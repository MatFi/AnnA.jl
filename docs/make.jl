using Documenter
using AnnA

makedocs(
    sitename = "AnnA.jl",
    modules = [AnnA],
    doctest = true,
    strict = true,
    checkdocs = :exports,
#    repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}",
    pages = Any["Home"=>"index.md", "API"=>"api.md"],
)


deploydocs(
   repo = "https://github.com/MatFi/AnnA.jl"
)