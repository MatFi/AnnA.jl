using Test, AnnA

using DiffEqBase
using ForwardDiff
using Setfield
using LinearAlgebra
using BenchmarkTools
using Logging
using Unitful
using Documenter
import AnnA: lcache

include("jacobian_test.jl")
include("operators_test.jl")
include("cell_test.jl")
include("rhs_test.jl")
include("solve_test.jl")
include("problems_test.jl")
include("load_parameter_test.jl")

DocMeta.setdocmeta!(AnnA, :DocTestSetup, :(using AnnA, Unitful); recursive=true)
@testset "Doctests" begin
    doctest(AnnA)
end