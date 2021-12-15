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
using OrdinaryDiffEq: Rodas4P, Rodas4P2, Rodas5, RadauIIA5, QNDF, QBDF
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

using SparseMatricesCSR

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
	write_template([parm::AbstracParameters=Parameters()], filename::String="./Parameters.jl"; force=false)

Creates a new template from an optinally provided `Parameters` struct at the path specified in `filename`, defaults to `Parameters.jl` in the current working directory. An existing file can be overwritten using the `force=true` keyword argument.
"""	
function write_template(parm::AbstractParameters, filename::String;force=false)
	template = joinpath(@__DIR__,"../Parameters.jl")
	cp(joinpath(@__DIR__,"../Parameters.jl"),filename;force=force)
	p= unpac_struct(parm,privates=false)
	lns =readlines(filename,keep=false)
	open(filename,"w") do f
		for line in lns
			#grep parametername
			pname = match(r"\s*(\S+)\s*=",line)
		
			if pname !== nothing
				pval = p[Symbol(pname[1])]

				regx = Regex("\\s+($(pname[1]))\\s*=\\s*([-\\+0-9\\.]+\\s*\\**\\s*[\\S]+)(\\s*#.*)")
	
				if pval isa Unitful.Quantity	
					ustring= string(unit(pval))
					ustring = replace(ustring,r"(\s)" => s"*" )	
					subs = SubstitutionString("    \\1 = $(string(ustrip(pval))*"*u\""*ustring*"\",")\\3")
				elseif pval isa Number
					subs =  SubstitutionString("    \\1 = $(pval),\\3")
				else
					subs =  SubstitutionString("    \\1 = t->0,\\3")
				end
				println(f, replace(line,regx => subs) )
			else
				println(f,line)
			end
		end
	end
end

write_template(filename::String="./Parameters.jl";kwargs...) = write_template(Parameters(), filename; kwargs...)
write_template(parm::AbstractParameters;kwargs...) = write_template(parm, "./Parameters.jl"; kwargs...)



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

include("./problems/ProblemSolution.jl")
include("./problems/IVProblem.jl") 
include("./problems/JscVocProblem.jl")
include("./problems/OCVDProblem.jl")
include("./problems/TPVProblem.jl")


export solve, IVProblem, JscVocProblem, OCVDProblem, TPVProblem, 
		Parameters, pulse, load_parameters, write_template,
		AlgControl
export @u_str
end # module