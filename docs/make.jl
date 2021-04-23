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
        "Interface"=>Any[
            "Simulation parameters" => "interface/parameters.md",
            "Problem definition" => "interface/problems.md",
            ],
        "Examples"=>Any[
            "IV" => "examples/iv_sim.md",
            "JscVoc" => "examples/jscvoc_sim.md",
            "OCVD"  =>   "examples/ocvd_sim.md",
            "TPV"  =>   "examples/tpv_sim.md",
            ],
        "API"=>"api.md"
    ],
    
)

if get(ENV, "CI", nothing) == "true"
    deploydocs(
        repo = "https://github.com/MatFi/AnnA.jl"
    )
end