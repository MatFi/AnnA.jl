#push!(LOAD_PATH, "../src/")

using Documenter
using AnnA

#DocMeta.setdocmeta!(AnnA, :DocTestSetup, :(using AnnA); recursive=true)

makedocs(
    sitename = "AnnA.jl",
    modules = [AnnA],
    doctest = :fix,
  #  strict = true,
    checkdocs = :exports,
#    repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}",
    pages = Any["Home"=>"index.md", "API"=>"api.md"],
    
)


deploydocs(
   repo = "https://github.com/MatFi/AnnA.jl"
)