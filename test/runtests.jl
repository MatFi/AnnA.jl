using Test, AnnABase

using DiffEqBase
using ForwardDiff
using Setfield
using LinearAlgebra
using BenchmarkTools
using Logging
using Unitful

@testset "Test all" begin
    include("jacobian_test.jl")
    include("operators_test.jl")
    include("cell_test.jl")
    include("rhs_test.jl")
    include("solve_test.jl")
    include("iv_problems_test.jl")
end
