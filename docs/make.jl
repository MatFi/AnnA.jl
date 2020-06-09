#using Pkg
#Pkg.build()
#Pkg.instantiate()
#Pkg.develop(PackageSpec(path=pwd()))

using Documenter
using AnnABase

makedocs(
         sitename = "AnnABase Documentation",
         modules = [AnnABase],
         doctest = true,
         strict = true,
         checkdocs = :exports,
         pages = Any[
             "Home" => "index.md",
             "API"  => "api.md"
             ],
         repo = "https://gitlab.com/MatFi/AnnABase.jl/blob/{commit}{path}#{line}"
         )
