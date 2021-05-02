#push!(LOAD_PATH, "../src/")

using Documenter
using AnnA


doctest = "fix" in ARGS ? :fix : true

DocMeta.setdocmeta!(AnnA, :DocTestSetup, :(using AnnA, Unitful, UnitfulRecipes); recursive=true)
makedocs(
    sitename = "AnnA.jl",
    modules = [AnnA],
    doctest = doctest,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
  #  strict = true,
    checkdocs = :exports,
    pages = Any[
        "Home"=>"index.md",
        "Working principle" => "working_principle.md",
        "Interface"=>Any[
            "Simulation parameters" => "interface/parameters.md",
            "Problem definition" => "interface/problems.md",
            ],
        "Examples"=>Any[
            "Predefined Protocols" => "examples/protos_sim.md",
            ],
        "API"=>"api.md"
    ],
    
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "https://github.com/MatFi/AnnA.jl"
    )
end