module AnnA
    using DataStructures: OrderedDict
    using OrdinaryDiffEq: Rodas4P, Rodas5
    using DiffEqBase: ODEFunction, ODEProblem,ODESolution, ODEProblem, solve
    import DiffEqBase, OrdinaryDiffEq
    using SparseDiffTools: matrix_colors, forwarddiff_color_jacobian!
    using DiffEqBase: get_tmp
    #using DiffEqCallbacks: AutoAbstol
   # using SteadyStateDiffEq
    import ForwardDiff
    using LinearAlgebra
    using NumericalIntegration: integrate
  #  using NLsolve: nlsolve, OnceDifferentiable
    using Roots
    using SparseArrays
    using Unitful
    using Unitful.DefaultSymbols
    using DataFrames: DataFrame
    using Waveforms: trianglewave
    using RecursiveArrayTools: copyat_or_push!
    using Setfield: setproperties
    using Interpolations: interpolate, SteffenMonotonicInterpolation
    using DelimitedFiles: readdlm

    import DiffEqCallbacks: TerminateSteadyState, CallbackSet
    import DiffEqBase: solve, dualcache # dos not work if using
  #  import DiffEqBase: solve
    import Base:length, getproperty, setproperty!, propertynames

    abstract type AbstractProblem end
    abstract type AbstractSolution end
    abstract type AbstractProblemSolution <: AbstractSolution end

    const CurrentDensity = Union{Quantity{T,Unitful.𝐈*Unitful.𝐋^-2,U}, Level{L,S,Quantity{T,Unitful.𝐈*Unitful.𝐋^-2,U}} where S where L} where U where T

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
