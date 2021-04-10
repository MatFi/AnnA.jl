var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"DocTestSetup = quote\n    using AnnA\nend","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [AnnA]\nOrder   = [:function, :type]","category":"page"},{"location":"api/#AnnA.NonDimensionalize-Tuple{AnnA.AbstractParameters}","page":"API","title":"AnnA.NonDimensionalize","text":"NonDimensionalize(p::AbstractParameters)\n\nPerforms the nondimensionalisation of the input parameters and returns a NodimParameters struct\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA._resample-Union{Tuple{T}, Tuple{Matrix{T}, Int64, Val{:xlog}}} where T<:Number","page":"API","title":"AnnA._resample","text":"_resample(data::Array{T,2},N::Int,::Val{:xlog}) where T <: Number\n\nresamples the data to a log scaled x-axis by appling peacwise quadratic fit arround the 'N' samle points ( this apporach is comparable to loesssmotthing with non-constant stepp size)\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.calc_ideality-Tuple{AnnA.IVSolution}","page":"API","title":"AnnA.calc_ideality","text":"calc_ideality(s::IVSolution;at_V=:all)\n\nreturns the ideality of an IV solution at the unitful voltage at_V. If no number is provided a array is returned.\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.init_guess-Tuple{AnnA.Grid, AnnA.NodimParameters, Any}","page":"API","title":"AnnA.init_guess","text":"init_guess(g::Grid,ndim::NodimParameters)\n\nReturns an initial guess for NLsolve root finding. The NodimParameter value is needed to guarantee consitancy of the guess.\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.initial_conditions!-Tuple{AnnA.Cell}","page":"API","title":"AnnA.initial_conditions!","text":"initial_conditions!(c::Cell)\n\nInplace mutation variant of initial_conditions(c::Cell)`\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.initial_conditions-Tuple{AnnA.Cell}","page":"API","title":"AnnA.initial_conditions","text":"initial_conditions(c::Cell)\n\nCalculates the initial condition of a given Cell using NLsolve. Returns the steady state solution vector.\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.pulse-Tuple{}","page":"API","title":"AnnA.pulse","text":"pulse(; Δh=1.0, h₀=0.0, w=1.0 - 2.0e-12, Δt=1.0e-12, tₑ=1.0)\n\ncreates a pulse. Monotonic interpolation is used, so the pulse shape is differentiable, allowing for autodiff.\n\nArguments:\n\nΔh: Pulse Amplitude\nh₀: Baseline\nw: Pulse width\nΔt: rise and fall time\ntₑ: end time (end of fall)\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.rdim_sol-Tuple{Any, AnnA.AbstractParameters}","page":"API","title":"AnnA.rdim_sol","text":"rdim_sol(sol,p::AbstractParameters)\n\na single time slice sol of the ODESolution is redimensionalized with p\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.rdim_sol-Tuple{SciMLBase.ODESolution, AbstractArray}","page":"API","title":"AnnA.rdim_sol","text":"rdim_sol(sol::DiffEqBase.ODESolution,t)\n\nreturns a redimensionalized solution vector at time t, the time is treated as seconds if no unit is given.\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.rdim_sol-Tuple{SciMLBase.ODESolution}","page":"API","title":"AnnA.rdim_sol","text":"rdim_sol(sol::DiffEqBase.ODESolution)\n\nRedimensionalizes the complete ODESolution at each time step.\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.unpac_struct-Tuple","page":"API","title":"AnnA.unpac_struct","text":"unpac_struct(s...)\n\nReturns a OrderedDict containing all field => value pairs. Allows for multiple structs as input. The result is then a single Dict. The letter arguments will overwrite previous one is existing already\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.length-Tuple{AnnA.Grid}","page":"API","title":"Base.length","text":"length(g::Grid)\n\ngives the length of the input-array which is 4*g.N+4+2*g.Nₑ+2*g.Nₕ\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.Gen_function","page":"API","title":"AnnA.Gen_function","text":"struct Gen_function{T, F, S} <: Function\n\nafter constructing this type using g! = Gen_function(αb, dir, light, tion) givs access to the generation function g!(R, n, p)\n\n#Arguments:\n\nαb: Nondim absorption coefficient\ndir: Light direction +1/-1\nlight: scalar tempoal function of generation\ntion: ionic timescale, needed for redimensionalize the time\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.Gen_function-Tuple{Any, Any, Any}","page":"API","title":"AnnA.Gen_function","text":"(g!::Gen_function)(G, x, t)\n\nInplace generation function\n\n#Arguments:\n\nG: Buffer Array\nx: Nodimensionalized positions Array\nt: Nondimensionalized time\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.IVProblem-Tuple{AnnA.AbstractParameters, Union{Tuple, AbstractArray}, Unitful.AbstractQuantity}","page":"API","title":"AnnA.IVProblem","text":"IVProblem(\n    parm::Parameters,\n    range::AbstractArray,\n    rate::Unitful.AbstractQuantity;\n    double_sweep = true,\n    alg_control = AlgControl(dtmin = 1e-20,\n        dt = 1e-6,\n        reltol = 1e-4,\n        abstol = 1e-12,\n        tend = (maximum(range) - minimum(range)) / abs(rate)\n    ),\n)\n\nCreates an IVProblem. The voltage parameter V defined in parm gets overwritten by the linear voltage sweep defined as V = t-> first(range) + rate*t. AlgControl.tend is forced to be (maximum(range) - minimum(range)) / abs(rate), an overwrite can be done on the finalized object using Setfield:     p = Setfield.setproperties(p::IVProblem, alg_control=AlgControl(...))\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.JscVocProblem-Tuple{AnnA.AbstractParameters, Number, Number}","page":"API","title":"AnnA.JscVocProblem","text":"\n\n\n\n","category":"method"},{"location":"api/#AnnA.LinOperator","page":"API","title":"AnnA.LinOperator","text":"struct LinOperator{T, S} <: AbstractOperators\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.NodimParameters","page":"API","title":"AnnA.NodimParameters","text":"struct NodimParameters{T, F <: Function, Fl <: Function, Fr <: Function, Fg <: Function, Fv <: Function}\n\nholds all nondimensionalized parameters\n\nArguments:\n\nλ: DESCRIPTION\nλ²: DESCRIPTION\nδ: DESCRIPTION\nχ: DESCRIPTION\nσ: DESCRIPTION\nκₚ: DESCRIPTION\nκₙ: DESCRIPTION\nαb: DESCRIPTION\ndpt: DESCRIPTION\ndpf: DESCRIPTION\nVbi: Vbi optained from bolzmann statistics\nwₑ: DESCRIPTION\nwₕ: DESCRIPTION\nκₑ: DESCRIPTION\nκₕ: DESCRIPTION\nrₑ: DESCRIPTION\nrₕ: DESCRIPTION\nλₑ: DESCRIPTION\nλₑ²: DESCRIPTION\nλₕ²: DESCRIPTION\nλₕ: DESCRIPTION\nkₑ: proportionality factor between ETL / PVSK electron densities\nkₕ: proportionality factor between PVSK / HTL hole densities\nR: Bulk recombination function\nRₗ: ETL/PVSK surface recombnation fnction\nRᵣ: PVSK/HTL surface recombination function\nG: Carrier generation function\nV: Voltage/kt\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.Operators-Tuple{AnnA.Grid}","page":"API","title":"AnnA.Operators","text":"Operators(g::Grid;dense::Bool=true)\n\nCompute the operators based on a given Grid. Use the dense keyword to switch between sparse ord dens representation\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.Pot_function","page":"API","title":"AnnA.Pot_function","text":"struct Pot_function{T, F} <: Function\n\nafter constructing this type using v =  Pot_function(VT, V, tion) givs  access to the generation function v(t)\n\n#Arguments:\n\nVT: thermal voltage\nV: temporal potential function\ntion: ionic timescale, needed for redimensionalize the time\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.Pot_function-Tuple{Any}","page":"API","title":"AnnA.Pot_function","text":"(v::Pot_function)(t)\n\nReturns the nondimensionalized potential at time t\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.Pulse","page":"API","title":"AnnA.Pulse","text":"struct Pulse{A, AA, T, TT, I} <: Function\n\nHolds all information about the pulse formation.\n\nArguments:\n\nΔh: Pulse Amplitude\nh₀: Baseline\nw: Pulse width\nΔt: rise and fall time\ntₑ: end time (end of fall)\nts: array of timepoints\nhs: array of amplitude points\nitp: interploating function\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.Pulse-Tuple{Any}","page":"API","title":"AnnA.Pulse","text":"(p::Pulse)(t)\n\ncalls the interpolating function from Pulse\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.Rec_function","page":"API","title":"AnnA.Rec_function","text":"struct Rec_function{T} <: Function\n\nafter constructing this type using r = Rec_function(nᵢ², kk, γ, τᵣ, rtrap) givs  access to the recombination function r(R, n, p)\n\n#Arguments:\n\nnᵢ²: nondimensionalized intrinsic carrier square\nkk: nondimensionalized bimolecular recombination rate\nγ: nondimensionalized rate constant for hole-dominated recombination\nτᵣ: ratio of SRH carrier lifetime\nrtrap: nondimensionalized k₂ (deep trap) parameter\n\n\n\n\n\n","category":"type"},{"location":"api/#AnnA.Rec_function-Tuple{AbstractArray, AbstractArray, AbstractArray}","page":"API","title":"AnnA.Rec_function","text":"(r!::Rec_function)(R::AbstractArray, n::AbstractArray, p::AbstractArray)\n\nCall the inplace recombination function. This is where the bulk recombination is calculated\n\n#Arguments:\n\nR: Buffer Array\nn: Electron Array\np: Hole Array\n\n\n\n\n\n","category":"method"},{"location":"api/#AnnA.Rec_function-Tuple{Number, Number}","page":"API","title":"AnnA.Rec_function","text":"(r::Rec_function)(n::T, p::T) where T\n\nOut of place definition of the recombination function. Calculates the recombination for a given n and p. This is intended to calculate the sruface crecombinations\n\n\n\n\n\n","category":"method"},{"location":"api/#LinearAlgebra.Tridiagonal-Union{Tuple{T}, Tuple{L}, Tuple{Integer, T, T, T}} where {L<:Number, T<:Union{AbstractVector{L}, Number}}","page":"API","title":"LinearAlgebra.Tridiagonal","text":"(LinearAlgebra.Tridiagonal(N::Integer, a::T, b::T, c::T) where T <: Union{Number, AbstractArray{L, 1}}) where L <: Number\n\nA simple wrapper for creating tridiagonal matrices a bit more convineantly\n\nArguments:\n\nN: Size of the NxN tridiagonal\na: row a\nb: row b\nc: row c\n\nExample\n\njulia> AnnA.Tridiagonal(4,1,-2,-1)\n4×4 LinearAlgebra.Tridiagonal{Int64, Vector{Int64}}:\n -2  -1   ⋅   ⋅\n  1  -2  -1   ⋅\n  ⋅   1  -2  -1\n  ⋅   ⋅   1  -2\n\n\n\n\n\n","category":"method"},{"location":"usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"#AnnA.jl","page":"Home","title":"AnnA.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AnnA.jl","category":"page"}]
}
