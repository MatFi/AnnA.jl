#push!(LOAD_PATH, "../src/")

using Documenter
using AnnA

DocMeta.setdocmeta!(AnnA, :DocTestSetup, :(using AnnA, Unitful); recursive=true)
makedocs(
    sitename = "AnnA.jl",
    modules = [AnnA],
    doctest = true,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
  #  strict = true,
    checkdocs = :exports,
#    repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}",
    pages = Any[
        "Home"=>"index.md",
        "Usage"=>Any[
            "Interface" => "interface.md"
            ],
        "Examples"=>Any[
            "IV-Problems" =>"examples/iv_sim.md"
            ],
        "API"=>"api.md"
    ],
    
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "https://github.com/MatFi/AnnA.jl"
    )
end