#push!(LOAD_PATH, "../src/")

using Documenter
using AnnA

#DocMeta.setdocmeta!(AnnA, :DocTestSetup, :(using AnnA); recursive=true)

makedocs(
    sitename = "AnnA.jl",
    modules = [AnnA],
    doctest = true,
  #  strict = true,
    checkdocs = :exports,
#    repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}",
    pages = Any[
        "Home"=>"index.md",
        "Usage"=>Any[
            "usage.md"
            ],
        "Examples"=>Any["examples/iv_sim.md"],
        "API"=>"api.md"
    ],
    
)


deploydocs(
   repo = "https://github.com/MatFi/AnnA.jl"
)