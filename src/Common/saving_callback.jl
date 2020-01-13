abstract type AbstractDDVariable end
abstract type AbstractDDSolution end

"""
    struct DDSpacialVariable{vType<:Number} <: AbstractDDVariable

The variables are stored as a `v::Vector{vType}`. So it becoms possible to
dispatch on them.
"""
struct DDSpacialVariable{vType<:Number} <: AbstractDDVariable
    v::Vector{vType}
end


"""
    struct DDTransientVariable{xType<:Number, vType<:DDSpacialVariable}  <: AbstractDDVariable

Stores a transient variable `v` as `Vector` of `DDSpacialVariable` together with its grid information `x`. Thereby the layers can be called seperately via indexing, the first instance of index selects and returns the infomration of the layer which can be indexed as well.
"""
struct DDTransientVariable{xType<:Number, vType<:DDSpacialVariable}  <: AbstractDDVariable
    x::Vector{Vector{xType}}
    v::Vector{Vector{vType}}
end

Base.getindex(var::DDTransientVariable, i::Int) = i > length(var.v) ? throw(BoundsError(var, i)) : var.v[i]

Base.iterate(var::DDTransientVariable, state=1) = state > length(var.v) ? nothing : (var.v[state], state+1)


"""
    DDTransientSolution{tType<:Number, vvType<:DDSpacialVariable, C<:Cell}

A struct used to save values of the time in `t::Vector{tType}` and
solution variables `I::Vector{vvType}, ϕ::Vector{vvType}, n::Vector{vvType}, p::Vector{vvType}`.
"""
struct DDTransientSolution{ C<:Cell,tType<:Number, jType<:Number,VType<:Number,D1,D2,D3,D4}
    cell::C
    t::Vector{tType}
    P::D1 #DDTransientVariable
    ϕ::D2 #DDTransientVariable
    n::D3 #DDTransientVariable
    p::D4 #DDTransientVariable
    J::Vector{jType}
    V::Vector{VType}
    #C
end

"""
    function DDTransientSolution(c::Cell)

Return `DDTransientSolution`
with empty storage vectors.
"""
function DDTransientSolution(c::Cell)
    t = 1.0*c.parameters.τᵢ
    u = c.g(c.u0)
    J = calculate_currents(c, u, 0.1, u) * c.parameters.jay
    V = u[2][3][end]*c.ndim.Vbi*c.parameters.VT
    u = rdim_sol(c,u)
    P = DDSpacialVariable.(u[1])
    ϕ = DDSpacialVariable.(u[2])
    n = DDSpacialVariable.(u[3])
    p = DDSpacialVariable.(u[4])

    P = DDTransientVariable([c.g.x], Vector{typeof(P)}())
    ϕ = DDTransientVariable([c.g.xₑ,c.g.x,c.g.xₕ], Vector{typeof(ϕ)}())
    n = DDTransientVariable([c.g.xₑ,c.g.x], Vector{typeof(n)}())
    p = DDTransientVariable([c.g.x,c.g.xₕ], Vector{typeof(p)}())

    DDTransientSolution(
        c, Vector{typeof(t)}(), P, ϕ, n, p, Vector{typeof(J)}(), Vector{typeof(V)}()
    )
end

mutable struct SavingAffect{Cell}
    saved_values::DDTransientSolution{Cell}
    saveiter::Int
end

function (affect!::SavingAffect)(integrator)

    u = affect!.saved_values.cell.g(integrator.u) #dcompose the soultion vecto
    u_prev = affect!.saved_values.cell.g(integrator.uprev)
    if !iszero(affect!.saveiter)
        J= calculate_currents(integrator.p, u, integrator.t - integrator.tprev, u_prev) * integrator.p.parameters.jay
    else
        J= calculate_currents(integrator.p, u, 0, u) * integrator.p.parameters.jay
    end

    affect!.saveiter += 1
    u = rdim_sol(integrator.p,u)
    V  =u[2][3][end]+integrator.p.ndim.Vbi*integrator.p.parameters.VT
    @debug abs(V)
#=
    if (integrator.p.mode == :oc && integrator.t> 1 && abs(integrator.p.g(integrator.u)[2][3][end]-integrator.p.g(integrator.p.u0)[2][3][end]) < 1e-20)
        terminate!(integrator)
    end
=#
    copyat_or_push!(affect!.saved_values.t, affect!.saveiter, integrator.t*integrator.p.parameters.τᵢ)
    copyat_or_push!(affect!.saved_values.J, affect!.saveiter, J)
    copyat_or_push!(affect!.saved_values.V, affect!.saveiter, V)
    copyat_or_push!(affect!.saved_values.P.v, affect!.saveiter, DDSpacialVariable.(u[1]), Val{false})
    copyat_or_push!(affect!.saved_values.ϕ.v, affect!.saveiter, DDSpacialVariable.(u[2]), Val{false})
    copyat_or_push!(affect!.saved_values.n.v, affect!.saveiter, DDSpacialVariable.(u[3]), Val{false})
    copyat_or_push!(affect!.saved_values.p.v, affect!.saveiter, DDSpacialVariable.(u[4]), Val{false})

    u_modified!(integrator, false)
end

function saving_initialize(cb, u, t, integrator)
    cb.affect!(integrator)
end


"""
    SavingCallback(c::Cell)

returns a tuble (callback, sol)
"""
function SavingCallback(c::Cell)
    tv =DDSpacialVariable(Vector{Float64}())
    sav = DDTransientSolution(c)
    affect! = SavingAffect(sav, 0)
    condtion = (u, t, integrator) -> true
    cb = DiscreteCallback(condtion, affect!;
                     initialize = saving_initialize,
                     save_positions=(false,false))
    return (cb,sav)
end


#=
    @test typeof(sol) <: DDTransientSolution

    @test typeof(sol(0)) <: DDSpacialSolution

    @test typeof(sol.n) <: Vector{DDSpacialVariable}



    @test typeof(sol.n(0)) <: DDTransientLocalVariable

    @test typeof(sol(0).n(0)) <: Number

    @test sol[1] <: DDSpacialSolution

    @test sol.n[1] <: DDTransientLocalVariable

    @test sol[1].n[1] <: Number
=#
