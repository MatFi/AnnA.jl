module AnnABase

    using DataStructures: OrderedDict
    using OrdinaryDiffEq
    using DiffEqBase: ODEFunction, ODEProblem, solve
    using SparseDiffTools: matrix_colors, forwarddiff_color_jacobian!
    using DiffEqBase: get_tmp
    #using DiffEqCallbacks: AutoAbstol
    using SteadyStateDiffEq
    using ForwardDiff
    using LinearAlgebra
    using NumericalIntegration: integrate
    using NLsolve: nlsolve, OnceDifferentiable
    using Roots
    using SparseArrays
    using Unitful
    using Unitful.DefaultSymbols

    using RecursiveArrayTools: copyat_or_push!
    using Setfield: setproperties
    using Interpolations
    using DelimitedFiles: readdlm

    import DiffEqCallbacks: TerminateSteadyState
    import DiffEqBase.dualcache # dos not work if using
    import DiffEqBase.solve
    import Base:length, getproperty, setproperty!, propertynames

    abstract type AbstractProblem end
    abstract type AbstractSolution end
    abstract type AbstractProblemSolution <: AbstractSolution end

    const CurrentDensity = Union{Quantity{T,Unitful.𝐈*Unitful.𝐋^-2,U}, Level{L,S,Quantity{T,Unitful.𝐈*Unitful.𝐋^-2,U}} where S where L} where U where T

    include("./Common/helpers.jl")
    include("./helpers/dataloader.jl")
    include("parameters.jl")
    include("./Common/nondimensionalise.jl")
    include("./Common/grid.jl")
    include("./Common/operators.jl")
    include("./Common/mass_matrix.jl")
    include("./Common/make_cell.jl")
    include("./Common/rhs.jl")
    include("./Common/initial_conditions.jl")
    include("./Common/jacobian.jl")
    include("./Common/solve.jl")
#    include("./Common/saving_callback.jl")
    include("./Common/redimensionalize.jl")
    include("./Common/calculate_externals.jl")

    include("./Problems/IVProblem.jl")
    include("./Problems/JscVocProblem.jl")
    include("./Problems/OCVDProblem.jl")


end # module
