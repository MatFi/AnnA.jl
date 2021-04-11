module AnnA
    using DataFrames:DataFrame
    using DataStructures:OrderedDict
    using DelimitedFiles:readdlm
    using DiffEqBase: get_tmp, ODEFunction, ODEProblem,ODESolution, solve
    using DiffEqCallbacks: TerminateSteadyState, CallbackSet
    import ForwardDiff
    using Interpolations: interpolate, SteffenMonotonicInterpolation
    using LinearAlgebra
    using NumericalIntegration:integrat
    using OrdinaryDiffEq: Rodas4P, Rodas5
    import DiffEqBase, OrdinaryDiffEq
    using RecursiveArrayTools:copyat_or_push!
    using Roots
    using Setfield:setproperties
    using SparseArrays
    using SparseDiffTools: matrix_colors, forwarddiff_color_jacobian!
    using Unitful
    using Unitful.DefaultSymbols
    
  	import Base:length, getproperty, setproperty!, propertynames

    abstract type AbstractProblem end
    abstract type AbstractSolution end
    abstract type AbstractProblemSolution <: AbstractSolution end

    const CurrentDensity = Union{Quantity{T,Unitful.𝐈 * Unitful.𝐋^-2,U},Level{L,S,Quantity{T,Unitful.𝐈 * Unitful.𝐋^-2,U}} where S where L} where U where T

    include("parameters.jl")
    include("./routines/helpers.jl")
    include("./helpers/dataloader.jl")
    include("./routines/nondimensionalise.jl")
    include("./routines/grid.jl")
    include("./routines/operators.jl")
    include("./routines/mass_matrix.jl")
    include("./routines/make_cell.jl")
    include("./routines/rhs.jl")
    include("./routines/initial_conditions.jl")
    include("./routines/jacobian.jl")
    include("./routines/solve.jl")
    include("./routines/redimensionalize.jl")
    include("./routines/calculate_externals.jl")

    include("./problems/IVProblem.jl") 
    include("./problems/JscVocProblem.jl")
    include("./problems/OCVDProblem.jl")
    include("./problems/TPVProblem.jl")
    include("./problems/ProblemSolution.jl")

    export solve, IVProblem, JscVocProblem, OCVDProblem, TPVProblem, Parameters

end # module
