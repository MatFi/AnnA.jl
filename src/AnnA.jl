module AnnA
using DataFrames: DataFrame
using DataStructures: OrderedDict
using DelimitedFiles: readdlm
import DiffEqBase: get_tmp, ODEFunction, ODEProblem, ODESolution, solve
using DiffEqCallbacks: TerminateSteadyState, CallbackSet
import ForwardDiff
using Interpolations: interpolate, SteffenMonotonicInterpolation
using LinearAlgebra
using NumericalIntegration: integrate
using OrdinaryDiffEq: Rodas4P, Rodas4P2, Rodas5, RadauIIA5
using DiffEqBase
using DiffEqBase: get_tmp
	using OrdinaryDiffEq
using Roots
using Setfield:setproperties
using SparseArrays
using SparseDiffTools: matrix_colors, forwarddiff_color_jacobian!
using Unitful
using Unitful.DefaultSymbols
using LoopVectorization
import Base:length, getproperty, setproperty!, propertynames
abstract type AbstractProblem end
abstract type AbstractSolution end
abstract type AbstractProblemSolution <: AbstractSolution end

const CurrentDensity = Union{Quantity{T,Unitful.𝐈 * Unitful.𝐋^-2,U},Level{L,S,Quantity{T,Unitful.𝐈 * Unitful.𝐋^-2,U}} where S where L} where U where T

#
include("parameters.jl")

"""
	load_parameters(
		p::AbstractParameters=Parameters(),
		file::String=joinpath(pwd(),"Parameters.jl");
		kwargs...)

Creates a new `Parameters` object. A list of modifications to the defaults can be passed as keyword arguments. Defaults can be supplied in different ways. By 
	1. placing a `Parameters.jl` file in the current working directory. A template can be found in the root folder of source files or created [`write_template`](@ref). 
	2. passing a already created `Parameters` object `p` as first argument.
	3. specifing a path to a valid file.
	4. doing all together. This way `p` is copied first and modified with the contents of `file` before the `kwargs` are applied. 

!!! note "Parameters.jl file" 

	The `Parameters.jl` file does not have do contain the full specifications of all parameters., it will fallback to the `AbstractParameters` struct when specified or to the internal defaults if not.  

"""	
function load_parameters(
	p::AbstractParameters=Parameters(),
	file::String=joinpath(pwd(),"Parameters.jl");
	kwargs...)

	if isfile(file)
		pd = include(file)
		pn = setproperties(p;pd...)
		return setproperties(pn;kwargs...)
	else
		throw(LoadError(file,0,"'$(file)' not found. Create that file, specify a correct path, or use the internal defaults 'Parameters(;kwargs...)' instead"))
	end
end
load_parameters(file::String; kwargs...) = load_parameters(Parameters(), file; kwargs...)


"""
	write_template(;force=false)

Creates a new `Parameters.jl` template in the current working directory. An existing file can be overwritten using the `force=true` keyword argument.
"""	
function write_template(;force=false)
	cp(joinpath(@__DIR__,"../Parameters.jl"),"./Parameters.jl";force=force)
end

include("./routines/helpers.jl")
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

export solve, IVProblem, JscVocProblem, OCVDProblem, TPVProblem, 
		Parameters, pulse, load_parameters
export @u_str
end # module